import KratosMultiphysics
import sympy
import KratosMultiphysics.sympy_fe_utilities as KratosSympy

from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src import generate_convective_flux
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src import generate_diffusive_flux
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src import generate_source_term
from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src import generate_stabilization_matrix

from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src \
    .symbolic_parameters import FormulationParameters, ShockCapturingParameters, ShockCapturingNodalParameters

from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src \
    .symbolic_geometry import GeometryDataFactory

from KratosMultiphysics.FluidDynamicsApplication.symbolic_generation.compressible_navier_stokes.src \
    .defines import CompressibleNavierStokesDefines as defs


class CompressibleNavierStokesSymbolicGenerator:
    """
    This class is in charge of generating a element to solve the compressible
    Navier-Stokes using the conservative magnitudes (density, momentum, energy)
    as solution.

    Formulation
    -----------
    Bayona (2017). Adaptive Mesh Simulations of Compressible Flows using Stabilized
    Formulations. Chapter 3.


    Testing
    -------
    If you modify this file without intending to modify its output, run test
    >> compressible_navier_stokes_symbolic_generator_test.py
    in order to ensure that the output stays the same.

    If you intend to modify the output, then run the aforementioned test and update
    the test checksum.

    The test is quite slow so it is only run in nightly CI, not raising any
    alarms in your pull request, so make sure to check manually.
    """

    def __init__(self, settings):
        settings.RecursivelyValidateAndAssignDefaults(self.GetDefaultParameters())

        self.write_mode = settings["mode"].GetString()
        self.is_explicit = settings["explicit"].GetBool()
        self.simplify = settings["do_simplifications"].GetBool()
        self.shock_capturing = settings["shock_capturing"].GetBool()
        self.echo_level = settings["echo_level"].GetInt()

        self.subscales_types = [s for (s, enabled) in settings["subscales"].items() if enabled]

        self.geometry = GeometryDataFactory(settings["geometry"].GetString())

        self.template_filename = settings["template_filename"].GetString()
        self.output_filename = settings["output_filename"].GetString()
        self._GenerateFiles()

        self._print(2, settings)

    def _CollectAndReplace(self, target_substring, expression, name):
        # If integrated during run-time, the assignment has to be an accumulation:
        assignment = " = " if self.geometry.symbolic_integration else " += "

        out = KratosSympy.OutputVector_CollectingFactors(expression, name, self.write_mode, replace_indices=False, assignment_op=assignment)

        self.outstring = self.outstring.replace(target_substring, out)

    @classmethod
    def GetDefaultParameters(cls):
        return KratosMultiphysics.Parameters("""
        {
            "mode": "c",
            "explicit": true,
            "do_simplifications": false,
            "geometry": "triangle",
            "shock_capturing": true,
            "subscales": {
                "ASGS" : true,
                "OSS" : true
            },
            "template_filename" : "PLEASE PROVIDE A template_filename",
            "output_filename"   : "symbolic_generator_name_not_provided.cpp",
            "echo_level" : 1
        }""")

    def _print(self, lvl, *text):
        if self.echo_level >= lvl:
            print(*text)

    def _GenerateFiles(self):
        with open(self.template_filename, "r") as template:
            self.outstring = template.read()

        # Touching the outfile just to make sure it is available
        with open(self.output_filename, "w") as tmp:
            tmp.write("This file is currently being processed by" + __file__)

    def _ComputeNonLinearOperator(self, A, H, S, Ug):
        L = defs.ZeroVector(self.geometry.blocksize)
        for j in range(self.geometry.ndims):
            # Convective operator product (A x grad(U))
            A_j = A[j]
            H_j = H.col(j)
            L += A_j * H_j
            # Diffusive flux
            # Note that the diffusive flux is not added as it will involve 2nd
            # order derivatives that vanish when introducing the linear FE
            # discretization

        # Source term addition
        L -= S * Ug
        return L

    def _ComputeNonLinearAdjointOperator(self, A, H, Q, S, Ug, V):
        L_adj = defs.ZeroVector(self.geometry.blocksize)
        for j in range(self.geometry.ndims):
            Q_j = Q.col(j)
            H_j = H.col(j)
            # Convective operator product
            A_j_trans = A[j].transpose()
            L_adj += A_j_trans * Q_j
            aux_conv = defs.ZeroMatrix(self.geometry.blocksize, self.geometry.blocksize)
            for m in range(self.geometry.blocksize):
                for n in range(self.geometry.blocksize):
                    A_j_trans_mn = A_j_trans[m, n]
                    for p in range(self.geometry.blocksize):
                        aux_conv[m, n] += sympy.diff(A_j_trans_mn, Ug[p]) * H_j[p]
            L_adj += aux_conv * V
            # Diffusive operator product
            # Note that the adjoint diffusive flux is not added as it will
            # involve 2nd order derivatives that vanish when introducing
            # the linear FE discretization
        # Source term addition
        L_adj += S.transpose() * V

        return L_adj

    def _ComputeVariationalFormulation(self, A, acc, G, H, L_adj, Q, S, Ug, V):
        # Mass (inertial) term - FE scale (only computed in the implicit case)
        n1 = defs.ZeroMatrix(1, 1) if self.is_explicit else -V.T*acc

        # Convective term - FE scale
        conv_flux = defs.ZeroVector(self.geometry.blocksize)
        for j in range(self.geometry.ndims):
            conv_flux += A[j] * H.col(j)
        n2 = - V.transpose() * conv_flux

        # Diffusive term - FE scale
        n3 = defs.ZeroMatrix(1, 1)
        for j in range(self.geometry.ndims):
            for k in range(self.geometry.blocksize):
                n3[0, 0] += Q[k, j] * G[k, j]

        # Source term - FE scale
        n4 = V.transpose() * (S * Ug)

        # VMS_adjoint - Subscales
        subscales = defs.Vector('subscales', self.geometry.blocksize)
        n5 = L_adj.transpose() * subscales

        # Variational formulation (Galerkin functional)
        self._print(1, " - Compute variational formulation")
        rv = n1 + n2 + n3 + n4 + n5  # Implicit case with inertial term n1

        return (rv, subscales)

    def _ComputeOSSProjectionsAtGaussPoint(self, acc, bdf, dUdt, f, forcing_terms, H, i_gauss, mg, projections, res, rg, U, Ug, Un, Unn):
        # Get Gauss point geometry data
        N = self.geometry.N_gauss(i_gauss)

        # Data interpolation at the gauss point
        U_gauss = U.transpose() * N
        f_gauss = forcing_terms["force"].transpose() * N
        r_gauss = (forcing_terms["thermal"].transpose()*N)[0]
        mass_gauss = (forcing_terms["mass"].transpose()*N)[0]

        if self.is_explicit:
            # In the explicit case, the acceleration is linearised taking the
            # previous step one. Note that in the explicit case this
            # acceleration is only used in the calculation of the stabilization
            # terms
            acc_gauss = dUdt.transpose()*N
        else:
            # In the implicit case, calculate the time derivatives with the
            # BDF2 formula
            acc_gauss = (bdf[0] * U + bdf[1] * Un + bdf[2] * Unn).transpose()*N

        # Gradients computation
        grad_U = KratosSympy.DfjDxi(self.geometry.DN(), U).transpose()

        # Substitute the symbols in the residual
        res_gauss = res.copy()
        KratosSympy.SubstituteMatrixValue(res_gauss, Ug, U_gauss)
        KratosSympy.SubstituteMatrixValue(res_gauss, acc, acc_gauss)
        KratosSympy.SubstituteMatrixValue(res_gauss, H, grad_U)
        KratosSympy.SubstituteMatrixValue(res_gauss, f, f_gauss)
        KratosSympy.SubstituteScalarValue(res_gauss, rg, r_gauss)
        KratosSympy.SubstituteScalarValue(res_gauss, mg, mass_gauss)

        # Add the projection contributions
        for i_node in range(self.geometry.nnodes):
            # Note that the weights will be added later on in the cpp
            projections["rho"][i_node] += N[i_node] * res_gauss[0]
            for d in range(self.geometry.ndims):
                i_mom = i_node * self.geometry.ndims + d
                projections["momentum"][i_mom] += N[i_node] * res_gauss[1 + d]
            projections["energy"][i_node] += N[i_node] * res_gauss[self.geometry.ndims + 1]

    def _OutputProjections(self, res_rho_proj, res_mom_proj, res_tot_ener_proj):
        dim = self.geometry.ndims
        self._CollectAndReplace("//substitute_rho_proj_{}D".format(dim), res_rho_proj, "rho_proj")
        self._CollectAndReplace("//substitute_mom_proj_{}D".format(dim), res_mom_proj, "mom_proj")
        self._CollectAndReplace("//substitute_tot_ener_proj_{}D".format(dim), res_tot_ener_proj, "tot_ener_proj")

    def _SubstituteSubscales(self, res, res_proj, rv, subscales, subscales_type, Tau):
        rv_gauss = rv.copy()
        if subscales_type == "ASGS":
            subs = Tau * res
        elif subscales_type == "OSS":
            subs = Tau * (res - res_proj)
        else:
            raise ValueError("Unrecognized subscales type: {}".format(subscales_type))

        KratosSympy.SubstituteMatrixValue(rv_gauss, subscales, subs)
        return rv_gauss

    def _ComputeResidualAtGaussPoint(self, acc, bdf, dUdt, f, forcing_terms, H, i_gauss, mg, params, Q, res_proj, ResProj, rg, rv_gauss, sc_nodes, sc_params, subscales_type, Tau, U, Ug, Un, Unn, V, w):
        self._print(1, "    Gauss point: " + str(i_gauss))

        # Get Gauss point geometry data
        Ng = self.geometry.N_gauss(i_gauss)

        # Data interpolation at the gauss point
        U_gauss = U.transpose() * Ng
        w_gauss = w.transpose() * Ng

        f_gauss = forcing_terms["force"].transpose() * Ng
        r_gauss = (forcing_terms["thermal"].transpose()*Ng)[0]
        mass_gauss = (forcing_terms["mass"].transpose()*Ng)[0]

        alpha_sc_gauss = (sc_nodes.alpha.transpose()*Ng)[0]
        mu_sc_gauss    = (sc_nodes.mu.transpose()*Ng)[0]
        beta_sc_gauss  = (sc_nodes.beta.transpose()*Ng)[0]
        lamb_sc_gauss  = (sc_nodes.lamb.transpose()*Ng)[0]

        if not self.is_explicit:
            # In the implicit case, calculate the time derivatives with the
            # BDF2 formula
            acc_gauss = (bdf[0] * U + bdf[1] * Un + bdf[2] * Unn).transpose()*Ng
        else:
            # In the explicit case, the acceleration is linearised taking the
            # previous step one. Note that in the explicit case this
            # acceleration is only used in the calculation of the stabilization
            # terms
            acc_gauss = dUdt.transpose()*Ng

        # Gauss pt. stabilization matrix calculation
        self._print(1, "    - Compute stabilization matrix on Gauss pt.")
        if self.shock_capturing:
            tau_gauss = generate_stabilization_matrix.ComputeStabilizationMatrixOnGaussPoint(params, U_gauss, f_gauss, r_gauss, mu_sc_gauss, lamb_sc_gauss)
        else:
            tau_gauss = generate_stabilization_matrix.ComputeStabilizationMatrixOnGaussPoint(params, U_gauss, f_gauss, r_gauss)

        # If OSS, residual projections interpolation
        res_proj_gauss = ResProj.T * Ng if subscales_type == "OSS" else None

        # Gradients computation
        grad_U = KratosSympy.DfjDxi(self.geometry.DN(), U).transpose()
        grad_w = KratosSympy.DfjDxi(self.geometry.DN(), w).transpose()

        self._print(1, "    - Substitution in the variational formulation")
        KratosSympy.SubstituteMatrixValue(rv_gauss, Ug, U_gauss)
        KratosSympy.SubstituteMatrixValue(rv_gauss, acc, acc_gauss)
        KratosSympy.SubstituteMatrixValue(rv_gauss, H, grad_U)
        KratosSympy.SubstituteMatrixValue(rv_gauss, V, w_gauss)
        KratosSympy.SubstituteMatrixValue(rv_gauss, Q, grad_w)
        KratosSympy.SubstituteMatrixValue(rv_gauss, Tau, tau_gauss)
        KratosSympy.SubstituteMatrixValue(rv_gauss, f, f_gauss)
        KratosSympy.SubstituteScalarValue(rv_gauss, rg, r_gauss)
        KratosSympy.SubstituteScalarValue(rv_gauss, mg, mass_gauss)
        KratosSympy.SubstituteScalarValue(rv_gauss, sc_params.alpha, alpha_sc_gauss)
        KratosSympy.SubstituteScalarValue(rv_gauss, sc_params.mu, mu_sc_gauss)
        KratosSympy.SubstituteScalarValue(rv_gauss, sc_params.beta, beta_sc_gauss)
        KratosSympy.SubstituteScalarValue(rv_gauss, sc_params.lamb, lamb_sc_gauss)
        if subscales_type == "OSS":
            KratosSympy.SubstituteMatrixValue(rv_gauss, res_proj, res_proj_gauss)

        # Accumulate in the total value
        return rv_gauss

    def _ComputeLHSandRHS(self, rv_tot, U, w):
        # Set the DOFs and test function matrices to do the differentiation
        dofs = defs.ZeroVector(self.geometry.ndofs)
        testfunc = defs.ZeroVector(self.geometry.ndofs)
        for i in range(0, self.geometry.nnodes):
            for j in range(0, self.geometry.blocksize):
                dofs[i*self.geometry.blocksize + j] = U[i, j]
                testfunc[i*self.geometry.blocksize + j] = w[i, j]

        # Compute LHS and RHS
        self._print(1, "    Compute RHS")
        rhs = KratosSympy.Compute_RHS(rv_tot.copy(), testfunc, self.simplify)

        if not self.is_explicit:
            self._print(1, "    Compute LHS")
            lhs = KratosSympy.Compute_LHS(rhs, testfunc, dofs, self.simplify)
        else:
            lhs = None

        return(lhs, rhs)

    def _OutputLHSandRHS(self, lhs, rhs, subscales_type):
        # Reading and filling the template file
        self._print(1, "    Substituting outstring from {}".format(self.template_filename))

        target = "//substitute_rhs_{}D_{}".format(self.geometry.ndims, subscales_type)
        self._CollectAndReplace(target, rhs, "rRightHandSideBoundedVector")

        if self.is_explicit:
            return

        target = "//substitute_lhs_{}D_{}".format(self.geometry.ndims, subscales_type)
        self._CollectAndReplace(target, lhs, "rLeftHandSideBoundedVector")

    def _ReplaceWarningMessage(self):
        message = "\n".join([
            "/**",
            " *",
            " *                          WARNING! THIS FILE IS READ-ONLY",
            " *",
            " * This file has been auto-generated by compressible_navier_stokes_symbolic_generator.py",
            " * located in the symbolic_generation directories of the FluidDynamicsApplication",
            " *",
            " * Any modifications to this file will be overwritten if and when that script is run again.",
            " *",
            " * In order to do any lasting changes, modify the template used by the script:",
            " * " + self.template_filename,
            " * located in the symbolic_generation directories of the FluidDynamicsApplication.",
            " *",
            " * In order to change the formulation you will have to modify the script itself.",
            " */"
        ])

        self.outstring = self.outstring.replace("//automated_message", message)

    def Generate(self):
        self._ReplaceWarningMessage()

        self._print(1, "Computing geometry: {}".format(self.geometry))

        dim = self.geometry.ndims
        n_nodes = self.geometry.nnodes
        block_size = self.geometry.blocksize

        params = FormulationParameters(self.geometry)

        # Unknowns
        U = defs.Matrix('data.U', n_nodes, block_size, real=True)

        # Residuals projection
        ResProj = defs.Matrix('data.ResProj', n_nodes, block_size, real=True)

        if self.is_explicit:
            # Unknowns time derivatives
            dUdt = defs.Matrix('data.dUdt', n_nodes, block_size, real=True)
            # Unknowns in previous steps
            Un = None
            Unn = None
        else:
            dUdt = None
            Un = defs.Matrix('data.Un', n_nodes, block_size, real=True)
            Unn = defs.Matrix('data.Unn', n_nodes, block_size, real=True)

        # Test functions defintion
        w = defs.Matrix('w', n_nodes, block_size, real=True)

        # External terms definition
        forcing_terms = {
            "mass":    defs.Vector('data.m_ext', n_nodes, real=True),
            "thermal": defs.Vector('data.r_ext', n_nodes, real=True),
            "force":   defs.Matrix('data.f_ext', n_nodes, dim, real=True)
        }

        # Nodal artificial magnitudes
        sc_nodes = ShockCapturingNodalParameters(self.geometry)

        # Backward differantiation coefficients
        if self.is_explicit:
            bdf = None
        else:
            bdf = [sympy.Symbol('bdf{}'.format(i)) for i in range(3)]

        # Construction of the variational equation
        Ug  = defs.Vector('Ug', block_size)             # Dofs vector
        H   = defs.Matrix('H', block_size, dim)         # Gradient of , real=TrueU
        mg  = sympy.Symbol('mg')                         # Mass source term
        f   = defs.Vector('f',  dim)                    # Body force vector
        rg  = sympy.Symbol('rg')                         # Thermal source/sink term
        V   = defs.Vector('V', block_size)              # Test function
        Q   = defs.Matrix('Q', block_size, dim)         # Gradient of , real=TrueV
        acc = defs.Vector('acc', block_size)            # Derivative of Dofs/Time
        G   = defs.Matrix('G', block_size, dim)         # Diffusive Flux matri, real=Truex
        res_proj = defs.Vector('res_proj', block_size)  # Residuals projection for the OSS

        # Calculate the Gauss point residual
        # Matrix Computation
        self._print(1, " - Compute Source Matrix")
        S = generate_source_term.ComputeSourceMatrix(Ug, mg, f, rg, params)

        self._print(1, " - Compute Euler Jacobian matrix")
        A = generate_convective_flux.ComputeEulerJacobianMatrix(Ug, params)

        if self.shock_capturing:
            sc_params = ShockCapturingParameters()
            self._print(1, " - Compute diffusive flux (shock capturing ON)")
            G = generate_diffusive_flux.ComputeDiffusiveFluxWithShockCapturing(Ug, H, params, sc_params)
        else:
            sc_params = None
            self._print(1, " - Compute diffusive flux (shock capturing OFF)")
            G = generate_diffusive_flux.ComputeDiffusiveFlux(Ug, H, params)

        self._print(1, " - Compute stabilization matrix")
        Tau = generate_stabilization_matrix.ComputeStabilizationMatrix(params)

        # Non-linear operator definition
        self._print(1, " - Compute non-linear operator")
        L = self._ComputeNonLinearOperator(A, H, S, Ug)

        # FE residuals definition
        # Note that we include the DOF time derivatives in both the implicit
        # and the explicit cases. It is required to include it in both cases
        # to calculate the subscale inertial component. In the implicit case
        # it is computed with the BDF formulas. In the explicit case it is
        # linearised by using the values already stored in the database.
        res = - acc - L

        self._print(1, " - Compute non-linear adjoint operator")
        L_adj = self._ComputeNonLinearAdjointOperator(A, H, Q, S, Ug, V)

        (rv, subscales) = self._ComputeVariationalFormulation(A, acc, G, H, L_adj, Q, S, Ug, V)

        # OSS Residual projections calculation
        # Calculate the residuals projection
        self._print(1, " - Calculate the projections of the residuals")
        projections = {
            "rho"      : defs.ZeroVector(n_nodes),
            "momentum" : defs.ZeroVector(n_nodes*dim),
            "energy"   : defs.ZeroVector(n_nodes)
        }

        for i_gauss in self.geometry.SymbolicIntegrationPoints():
            self._print(1, "   - Gauss point: " + str(i_gauss))
            self._ComputeOSSProjectionsAtGaussPoint(acc, bdf, dUdt, f, forcing_terms, H, i_gauss, mg, projections, res, rg, U, Ug, Un, Unn)

        # Output the projections
        self._OutputProjections(*projections.values())

        # Algebraic form calculation #
        for subscales_type in self.subscales_types:
            # Substitution of the discretized values at the gauss points
            # Loop and accumulate the residual in each Gauss point
            rv_tot = defs.ZeroMatrix(1, 1)

            self._print(1, " - Subscales type: " + subscales_type)
            self._print(1, " - Substitution of the discretized values at the gauss points")
            for i_gauss in self.geometry.SymbolicIntegrationPoints():
                # Substitute the subscales model
                rv_gauss = self._SubstituteSubscales(res, res_proj, rv, subscales, subscales_type, Tau)
                rv_tot += self._ComputeResidualAtGaussPoint(acc, bdf, dUdt, f, forcing_terms, H, i_gauss, mg, params, Q, res_proj, ResProj, rg, rv_gauss, sc_nodes, sc_params, subscales_type, Tau, U, Ug, Un, Unn, V, w)

            (lhs, rhs) = self._ComputeLHSandRHS(rv_tot, U, w)
            self._OutputLHSandRHS(lhs, rhs, subscales_type)

    def Write(self):
        self._print(1, " - Writing {}".format(self.output_filename))

        with open(self.output_filename, "w") as outputfile:
            outputfile.write(self.outstring)
            outputfile.close()

        self._print(1, "Geometry {} done".format(self.geometry))
