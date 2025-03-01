//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Joaquin Gonzalez-Usua
//
//

//These formulas are derived from "Derivatives.py" script

// System includes

// External includes

// Project includes
#include "containers/model.h"
#include "includes/checks.h"
#include "utilities/variable_utils.h"
#include "utilities/math_utils.h"

// Application includes
#include "swimming_DEM_application.h"
#include "porosity_solution_and_body_force_process.h"
#include "swimming_dem_application_variables.h"

// Other applications includes
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{

/* Public functions *******************************************************/
PorositySolutionAndBodyForceProcess::PorositySolutionAndBodyForceProcess(
    ModelPart& rModelPart)
    : Process(),
      mrModelPart(rModelPart)
{}

PorositySolutionAndBodyForceProcess::PorositySolutionAndBodyForceProcess(
    ModelPart& rModelPart,
    Parameters& rParameters)
    : Process(),
      mrModelPart(rModelPart)
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);
}

PorositySolutionAndBodyForceProcess::PorositySolutionAndBodyForceProcess(
    Model &rModel,
    Parameters &rParameters)
    : Process(),
      mrModelPart(rModel.GetModelPart(rParameters["model_part_name"].GetString()))
{
    // Check default settings
    this->CheckDefaultsAndProcessSettings(rParameters);

}


void PorositySolutionAndBodyForceProcess::CheckDefaultsAndProcessSettings(Parameters &rParameters)
{
    const Parameters default_parameters = GetDefaultParameters();

    rParameters.ValidateAndAssignDefaults(default_parameters);

    mDensity     = rParameters["benchmark_parameters"]["density"].GetDouble();
    mUchar       = rParameters["benchmark_parameters"]["u_char"].GetDouble();
    mDeltaAlpha  = rParameters["benchmark_parameters"]["delta_alpha"].GetDouble();
    mLength      = rParameters["benchmark_parameters"]["length"].GetDouble();
    mOmega       = rParameters["benchmark_parameters"]["omega"].GetDouble();
    mX1Origin    = rParameters["benchmark_parameters"]["x1_origin"].GetDouble();
    mX2Origin    = rParameters["benchmark_parameters"]["x2_origin"].GetDouble();
    mSqueezeAmplitude = rParameters["benchmark_parameters"]["squeeze_amplitude"].GetDouble();
    mNSafety     = rParameters["benchmark_parameters"]["n_safety"].GetDouble();
    mReynoldsNumber = rParameters["benchmark_parameters"]["n_reynolds"].GetDouble();
    mDamKohlerNumber = rParameters["benchmark_parameters"]["n_dam"].GetDouble();
    mInitialConditions = rParameters["benchmark_parameters"]["use_initial_conditions"].GetBool();
    mAlternativeFormulation = rParameters["benchmark_parameters"]["use_alternative_formulation"].GetBool();

    this->CalculateKinematicViscosity();

    double dynamic_viscosity = mViscosity * mDensity;

    this->CalculatePermeability(dynamic_viscosity);

}

void PorositySolutionAndBodyForceProcess::CalculateKinematicViscosity()
{
    mViscosity = mUchar * mLength / mReynoldsNumber;
}

void PorositySolutionAndBodyForceProcess::CalculatePermeability(double &dynamic_viscosity)
{
    mPermeability = dynamic_viscosity * mUchar / (mDamKohlerNumber * (2 * mViscosity * (mUchar/std::pow(mLength,2))));
}

const Parameters PorositySolutionAndBodyForceProcess::GetDefaultParameters() const
{
    const Parameters default_parameters( R"(
    {
                "model_part_name"          : "please_specify_model_part_name",
                "variable_name"            : "BODY_FORCE",
                "benchmark_name"           : "custom_body_force.vortex",
                "benchmark_parameters"     : {
                                                "velocity"    : 1.0,
                                                "length"      : 1.0,
                                                "density"     : 1.0,
                                                "frequency"   : 1.0,
                                                "alpha0"      : 0.7,
                                                "alpha_min"   : 0.5,
                                                "period"      : 0.1,
                                                "delta_alpha" : 0.25,
                                                "squeeze_amplitude" : 0.425,
                                                "x1_origin"   : 0.5,
                                                "x2_origin"   : 0.5,
                                                "omega"       : 5.0,
                                                "n_safety"    : 1.0,
                                                "n_reynolds"  : 1000.0,
                                                "n_dam"       : 0.0001,
                                                "u_char"      : 100.0,
                                                "use_alternative_formulation" : false
                },
                "compute_nodal_error"      : true,
                "print_convergence_output" : false,
                "output_parameters"        : {}
    }  )" );

    return default_parameters;
}

void PorositySolutionAndBodyForceProcess::Execute()
{
    this->ExecuteInitialize();
    this->ExecuteInitializeSolutionStep();
}

void PorositySolutionAndBodyForceProcess::ExecuteInitialize()
{}

void PorositySolutionAndBodyForceProcess::ExecuteBeforeSolutionLoop()
{
    this->SetFluidProperties();
    if (mInitialConditions == true)
    {
        this->SetInitialBodyForceAndPorosityField();
    }
}

void PorositySolutionAndBodyForceProcess::ExecuteInitializeSolutionStep()
{}

void PorositySolutionAndBodyForceProcess::ExecuteFinalizeSolutionStep() {}

/* Protected functions ****************************************************/

void PorositySolutionAndBodyForceProcess::SetInitialBodyForceAndPorosityField()
{
    const double dim = mrModelPart.GetProcessInfo()[DOMAIN_SIZE];
    const double delta_alpha = mDeltaAlpha;
    const double omega = mOmega;
    const double L = mLength;
    const double rho = mDensity;
    const double nu = mViscosity;
    const double squeeze_amplitude = mSqueezeAmplitude;
    const double n_safety = mNSafety;
    const double x10 = mX1Origin;
    const double x20 = mX2Origin;
    const double c_min = (1 - squeeze_amplitude);
    const double R = (c_min * L/2) / n_safety;
    Matrix inv_permeability = ZeroMatrix(dim,dim);

    const double c = (1 + squeeze_amplitude * std::sin(omega));

    double du1dt, du2dt, du11, du12, du111, du112, du121, du122, du21, du22, du211, du212, du221, du222;

    // Computation of the BodyForce and Porosity fields
    for (auto it_node = mrModelPart.NodesBegin(); it_node != mrModelPart.NodesEnd(); it_node++){

        const double x1 = it_node->X();
        const double x2 = it_node->Y();

        double& r_mass_source = it_node->FastGetSolutionStepValue(MASS_SOURCE);

        double& r_alpha = it_node->FastGetSolutionStepValue(FLUID_FRACTION);
        double& r_dalphat = it_node->FastGetSolutionStepValue(FLUID_FRACTION_RATE);

        double& r_alpha1 = it_node->FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT_X);
        double& r_alpha2 = it_node->FastGetSolutionStepValue(FLUID_FRACTION_GRADIENT_Y);

        double& r_body_force1 = it_node->FastGetSolutionStepValue(BODY_FORCE_X);
        double& r_body_force2 = it_node->FastGetSolutionStepValue(BODY_FORCE_Y);

        double& r_u1 = it_node->FastGetSolutionStepValue(EXACT_VELOCITY_X);
        double& r_u2 = it_node->FastGetSolutionStepValue(EXACT_VELOCITY_Y);

        Matrix& permeability = it_node->FastGetSolutionStepValue(PERMEABILITY);

        double& r_pressure = it_node->FastGetSolutionStepValue(EXACT_PRESSURE);

        if (this->IsInsideEllipticalSupport(x1, x2, c, R)){

            r_u1 = std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1);

            r_u2 = std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1);

            du1dt = 0.0;

            du2dt = 0.0;

            r_alpha = -delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1;

            r_dalphat = 0.0;

            r_alpha1 = delta_alpha*(2*x1 - 2*x10)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))/(std::pow(R,2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2));

            r_alpha2 =  delta_alpha*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2));

            du11 = Globals::Pi*std::sin(Globals::Pi*x2 + 0.2)*std::cos(Globals::Pi*x1 - 0.7)/(-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1) - delta_alpha*(2*x1 - 2*x10)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2));

            du12 = Globals::Pi*std::sin(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1) - delta_alpha*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2));

            du111 = -std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1) - 2*Globals::Pi*delta_alpha*(2*x1 - 2*x10)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x2 + 0.2)*std::cos(Globals::Pi*x1 - 0.7)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2)) - 2*delta_alpha*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2)) + 2*std::pow(delta_alpha,2)*std::pow((2*x1 - 2*x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),4)*std::exp(2 - 2/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),3)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),4)) - 2*delta_alpha*std::pow((2*x1 - 2*x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),4)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),3)) + delta_alpha*std::pow((2*x1 - 2*x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),4)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),4));

            du112 = std::pow(Globals::Pi,2)*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1) - Globals::Pi*delta_alpha*(2*x1 - 2*x10)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2)) - Globals::Pi*delta_alpha*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x2 + 0.2)*std::cos(Globals::Pi*x1 - 0.7)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2)) + 2*std::pow(delta_alpha,2)*(2*x1 - 2*x10)*(2*x2 - 2*x20)*std::exp(2 - 2/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),3)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),4)) - 2*delta_alpha*(2*x1 - 2*x10)*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),3)) + delta_alpha*(2*x1 - 2*x10)*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),4));

            du121 = std::pow(Globals::Pi,2)*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1) - Globals::Pi*delta_alpha*(2*x1 - 2*x10)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2)) - Globals::Pi*delta_alpha*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x2 + 0.2)*std::cos(Globals::Pi*x1 - 0.7)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2)) + 2*std::pow(delta_alpha,2)*(2*x1 - 2*x10)*(2*x2 - 2*x20)*std::exp(2 - 2/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),3)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),4)) - 2*delta_alpha*(2*x1 - 2*x10)*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),3)) + delta_alpha*(2*x1 - 2*x10)*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),4));

            du122 = -std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1) - 2*Globals::Pi*delta_alpha*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2)) - 2*delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2)) + 2*std::pow(delta_alpha,2)*std::pow((2*x2 - 2*x20),2)*std::exp(2 - 2/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),3)*std::pow((squeeze_amplitude*std::sin(omega) + 1),4)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),4)) - 2*delta_alpha*std::pow((2*x2 - 2*x20),2)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),4)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),3)) + delta_alpha*std::pow((2*x2 - 2*x20),2)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),4)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),4));

            du21 = -Globals::Pi*std::sin(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1) - delta_alpha*(2*x1 - 2*x10)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2));

            du22 = -Globals::Pi*std::sin(Globals::Pi*x2 + 0.2)*std::cos(Globals::Pi*x1 - 0.7)/(-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1) - delta_alpha*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2));

            du211 = -std::pow(Globals::Pi,2)*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1) + 2*Globals::Pi*delta_alpha*(2*x1 - 2*x10)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2)) - 2*delta_alpha*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2)) + 2*std::pow(delta_alpha,2)*std::pow((2*x1 - 2*x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),4)*std::exp(2 - 2/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),3)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),4)) - 2*delta_alpha*std::pow((2*x1 - 2*x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),4)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),3)) + delta_alpha*std::pow((2*x1 - 2*x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),4)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),4));

            du212 = std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1) + Globals::Pi*delta_alpha*(2*x1 - 2*x10)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x2 + 0.2)*std::cos(Globals::Pi*x1 - 0.7)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2)) + Globals::Pi*delta_alpha*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2)) + 2*std::pow(delta_alpha,2)*(2*x1 - 2*x10)*(2*x2 - 2*x20)*std::exp(2 - 2/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),3)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),4)) - 2*delta_alpha*(2*x1 - 2*x10)*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),3)) + delta_alpha*(2*x1 - 2*x10)*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),4));

            du221 = std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2)/(-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1) + Globals::Pi*delta_alpha*(2*x1 - 2*x10)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x2 + 0.2)*std::cos(Globals::Pi*x1 - 0.7)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2)) + Globals::Pi*delta_alpha*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2)) + 2*std::pow(delta_alpha,2)*(2*x1 - 2*x10)*(2*x2 - 2*x20)*std::exp(2 - 2/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),3)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),4)) - 2*delta_alpha*(2*x1 - 2*x10)*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),3)) + delta_alpha*(2*x1 - 2*x10)*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),4));

            du222 =  -std::pow(Globals::Pi,2)*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1) + 2*Globals::Pi*delta_alpha*(2*x2 - 2*x20)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::sin(Globals::Pi*x2 + 0.2)*std::cos(Globals::Pi*x1 - 0.7)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2)) - 2*delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,2)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),2)) + 2*std::pow(delta_alpha,2)*std::pow((2*x2 - 2*x20),2)*std::exp(2 - 2/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),3)*std::pow((squeeze_amplitude*std::sin(omega) + 1),4)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),4)) - 2*delta_alpha*std::pow((2*x2 - 2*x20),2)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),4)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),3)) + delta_alpha*std::pow((2*x2 - 2*x20),2)*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))))*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2)/(std::pow(R,4)*std::pow((-delta_alpha*std::exp(1 - 1/(1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)))) + 1),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),4)*std::pow((1 - std::pow((x1 - x10),2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2)/std::pow(R,2) - std::pow((x2 - x20),2)/(std::pow(R,2)*std::pow((squeeze_amplitude*std::sin(omega) + 1),2))),4));

            r_pressure = std::sin(Globals::Pi*x1)*std::cos(Globals::Pi*x2) + (-1.0 + std::cos(1.0))*std::sin(1.0);

            for (unsigned int d = 0; d < dim; ++d){
                permeability(d,d) = mPermeability;
            }

            permeability(dim-1,dim-1) = 1.0e+30;

        }else{

            for (unsigned int d = 0; d < dim; ++d){
                permeability(d,d) = 1.0e+30;
            }

            r_u1 = std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2);

            r_u2 = std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2);

            du1dt = 0.0;

            du2dt = 0.0;

            r_alpha = 1.0;

            r_dalphat = 0.0;

            r_alpha1 = 0.0;

            r_alpha2 = 0.0;

            du11 = Globals::Pi*std::sin(Globals::Pi*x2 + 0.2)*std::cos(Globals::Pi*x1 - 0.7);

            du12 = Globals::Pi*std::sin(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2);

            du111 = -std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2);

            du112 = std::pow(Globals::Pi,2)*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2);

            du121 = std::pow(Globals::Pi,2)*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2);

            du122 = -std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2);

            du21 = -Globals::Pi*std::sin(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2);

            du22 = -Globals::Pi*std::sin(Globals::Pi*x2 + 0.2)*std::cos(Globals::Pi*x1 - 0.7);

            du211 = -std::pow(Globals::Pi,2)*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2);

            du212 = std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2);

            du221 = std::pow(Globals::Pi,2)*std::sin(Globals::Pi*x1 - 0.7)*std::sin(Globals::Pi*x2 + 0.2);

            du222 = -std::pow(Globals::Pi,2)*std::cos(Globals::Pi*x1 - 0.7)*std::cos(Globals::Pi*x2 + 0.2);

            r_pressure = std::sin(Globals::Pi*x1)*std::cos(Globals::Pi*x2) + (-1 + std::cos(1.0))*std::sin(1.0);

        }

        double det_permeability = MathUtils<double>::Det(permeability);
        MathUtils<double>::InvertMatrix(permeability, inv_permeability, det_permeability, -1.0);

        Matrix sigma = nu * rho * inv_permeability;

        const double convective1 = r_u1 * du11 + r_u2 * du12;
        const double convective2 = r_u1 * du21 + r_u2 * du22;

        const double div_of_sym_grad1 = (1.0/2.0) * (2 * du111 + du212 + du122);
        const double div_of_sym_grad2 = (1.0/2.0) * (du121 + du211 + 2 * du222);

        const double grad_of_div1 = du111 + du221;
        const double grad_of_div2 = du112 + du222;

        const double press_grad1 = Globals::Pi*std::cos(Globals::Pi*x1)*std::cos(Globals::Pi*x2);;
        const double press_grad2 = -Globals::Pi*std::sin(Globals::Pi*x1)*std::sin(Globals::Pi*x2);

        if (mAlternativeFormulation){

            const double grad_alpha_sym_grad1 = (1.0/2.0) * (2 * r_alpha1 * du11 + r_alpha2 * (du21 + du12));
            const double grad_alpha_sym_grad2 = (1.0/2.0) * (r_alpha1 * (du12 + du21) + 2 * r_alpha2 * du22);

            const double grad_alpha_div1 = r_alpha1 * (du11 + du22);
            const double grad_alpha_div2 = r_alpha2 * (du11 + du22);

            r_body_force1 = r_alpha * du1dt + r_alpha * convective1 + r_alpha / rho * press_grad1 - 2 * nu * (r_alpha * div_of_sym_grad1 + grad_alpha_sym_grad1) + (2.0/3.0) * nu * (r_alpha * grad_of_div1 + grad_alpha_div1) + sigma(0,0) * r_u1 + sigma(1,0) * r_u1;

            r_body_force2 = r_alpha * du2dt + r_alpha * convective2 + r_alpha / rho * press_grad2 - 2 * nu * (r_alpha * div_of_sym_grad2 + grad_alpha_sym_grad2) + (2.0/3.0) * nu * (r_alpha * grad_of_div2 + grad_alpha_div2) + sigma(0,1) * r_u2 + sigma(1,1) * r_u2;

        }else{
            r_body_force1 = du1dt + convective1 + 1.0/rho * press_grad1 - 2 * nu * div_of_sym_grad1 + (2.0/3.0) * nu * grad_of_div1 + sigma(0,0) * r_u1 + sigma(1,0) * r_u1;

            r_body_force2 = du2dt + convective2 + 1.0/rho * press_grad2 - 2 * nu * div_of_sym_grad2 + (2.0/3.0) * nu * grad_of_div2 + sigma(0,1) * r_u2 + sigma(1,1) * r_u2;
        }

        r_mass_source = (r_dalphat + r_u1 * r_alpha1 + r_u2 * r_alpha2 + r_alpha * (du11 + du22));

        it_node->FastGetSolutionStepValue(VELOCITY_X) = r_u1;
        it_node->FastGetSolutionStepValue(VELOCITY_Y) = r_u2;
        it_node->FastGetSolutionStepValue(PRESSURE) = r_pressure;
    }

}

void PorositySolutionAndBodyForceProcess::SetFluidProperties()
{
    (mrModelPart.pGetProperties(1))->SetValue(DENSITY, mDensity);
    (mrModelPart.pGetProperties(1))->SetValue(DYNAMIC_VISCOSITY, mViscosity * mDensity);
    (mrModelPart.pGetProperties(1))->SetValue(VISCOSITY, mViscosity);

    block_for_each(mrModelPart.Elements(), [&](Element& rElement){
        rElement.SetProperties(mrModelPart.pGetProperties(1));
    });

    block_for_each(mrModelPart.Nodes(), [&](Node<3>& rNode){
        rNode.FastGetSolutionStepValue(VISCOSITY) = mViscosity;
        rNode.FastGetSolutionStepValue(DENSITY) = mDensity;
        rNode.FastGetSolutionStepValue(DYNAMIC_VISCOSITY) = mViscosity * mDensity;
    });
}

bool PorositySolutionAndBodyForceProcess::IsInsideEllipticalSupport(
    const double x1,
    const double x2,
    const double c,
    const double R)
{
    if (std::pow(c*(x1 - mX1Origin), 2) + std::pow(((x2 - mX2Origin) / c), 2) < std::pow(R, 2)){
        return true;
    }
    else{
        return false;
    }
}

/* Private functions ****************************************************/

};  // namespace Kratos.
