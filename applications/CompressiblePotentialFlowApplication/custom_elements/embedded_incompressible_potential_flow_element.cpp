//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  Kratos default license: kratos/license.txt
//
//  Main authors:    Marc Núñez, based on Iñigo Lopez and Riccardo Rossi work
//
#include "embedded_incompressible_potential_flow_element.h"
#include "compressible_potential_flow_application_variables.h"
#include "custom_utilities/potential_flow_utilities.h"
#include "fluid_dynamics_application_variables.h"

namespace Kratos
{
///////////////////////////////////////////////////////////////////////////////////////////////////
// Public Operations

template <int Dim, int NumNodes>
Element::Pointer EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, NodesArrayType const& ThisNodes, typename PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EmbeddedIncompressiblePotentialFlowElement>(
        NewId, this->GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Create(
    IndexType NewId, typename GeometryType::Pointer pGeom, typename PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EmbeddedIncompressiblePotentialFlowElement>(
        NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
Element::Pointer EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Clone(
    IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<EmbeddedIncompressiblePotentialFlowElement>(
        NewId, this->GetGeometry().Create(ThisNodes), this->pGetProperties());
    KRATOS_CATCH("");
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    const EmbeddedIncompressiblePotentialFlowElement& r_this = *this;
    const int wake = r_this.GetValue(WAKE);

    BoundedVector<double,NumNodes> distances;
    for(unsigned int i_node = 0; i_node<NumNodes; i_node++){
        distances[i_node] = this->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);
    }
    const bool is_embedded = PotentialFlowUtilities::CheckIfElementIsCutByDistance<Dim,NumNodes>(distances);

    if (is_embedded && wake == 0) {
        CalculateEmbeddedLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
        if (std::abs(rCurrentProcessInfo[STABILIZATION_FACTOR]) > std::numeric_limits<double>::epsilon()) {
            PotentialFlowUtilities::AddPotentialGradientStabilizationTerm<Dim, NumNodes>(*this,rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
        }
    }
    else {
        if (this->Is(STRUCTURE)) {
            CalculateKuttaWakeLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
        } else {
            BaseType::CalculateLocalSystem(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo);
        }
    }

    if (std::abs(rCurrentProcessInfo[PENALTY_COEFFICIENT]) > std::numeric_limits<double>::epsilon()) {
        PotentialFlowUtilities::AddKuttaConditionPenaltyTerm<Dim, NumNodes>(r_this,rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
    }
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateEmbeddedLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != NumNodes || rLeftHandSideMatrix.size2() != NumNodes)
        rLeftHandSideMatrix.resize(NumNodes, NumNodes, false);
    if (rRightHandSideVector.size() != NumNodes)
        rRightHandSideVector.resize(NumNodes, false);
    rLeftHandSideMatrix.clear();

    array_1d<double, NumNodes> potential;
    Vector distances(NumNodes);
    for(unsigned int i_node = 0; i_node<NumNodes; i_node++)
        distances(i_node) = this->GetGeometry()[i_node].GetSolutionStepValue(GEOMETRY_DISTANCE);

    potential = PotentialFlowUtilities::GetPotentialOnNormalElement<Dim, NumNodes>(*this);

    ModifiedShapeFunctions::Pointer pModifiedShFunc = this->pGetModifiedShapeFunctions(distances);
    Matrix positive_side_sh_func;
    ModifiedShapeFunctions::ShapeFunctionsGradientsType positive_side_sh_func_gradients;
    Vector positive_side_weights;
    pModifiedShFunc -> ComputePositiveSideShapeFunctionsAndGradientsValues(
        positive_side_sh_func,
        positive_side_sh_func_gradients,
        positive_side_weights,
        GeometryData::IntegrationMethod::GI_GAUSS_1);

    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];

    BoundedMatrix<double,NumNodes,Dim> DN_DX;
    for (unsigned int i_gauss=0;i_gauss<positive_side_sh_func_gradients.size();i_gauss++){
        DN_DX=positive_side_sh_func_gradients(i_gauss);
        noalias(rLeftHandSideMatrix) += free_stream_density*prod(DN_DX,trans(DN_DX))*positive_side_weights(i_gauss);
    }

    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, potential);
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::CalculateKuttaWakeLocalSystem(
    MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    // Note that the lhs and rhs have double the size
    if (rLeftHandSideMatrix.size1() != 2 * NumNodes ||
        rLeftHandSideMatrix.size2() != 2 * NumNodes)
        rLeftHandSideMatrix.resize(2 * NumNodes, 2 * NumNodes, false);
    if (rRightHandSideVector.size() != 2 * NumNodes)
        rRightHandSideVector.resize(2 * NumNodes, false);
    rLeftHandSideMatrix.clear();
    rRightHandSideVector.clear();

    PotentialFlowUtilities::ElementalData<NumNodes,Dim> data;

    // Calculate shape functions
    GeometryUtils::CalculateGeometryData(this->GetGeometry(), data.DN_DX, data.N, data.vol);

    const double free_stream_density = rCurrentProcessInfo[FREE_STREAM_DENSITY];

    data.distances = PotentialFlowUtilities::GetWakeDistances<Dim, NumNodes>(*this);

    BoundedMatrix<double, NumNodes, NumNodes> lhs_total = data.vol*free_stream_density*prod(data.DN_DX, trans(data.DN_DX));

    for (unsigned int i = 0; i < NumNodes; ++i)
    {
        for (unsigned int j = 0; j < NumNodes; ++j)
        {
            rLeftHandSideMatrix(i, j) = lhs_total(i, j);
            rLeftHandSideMatrix(i + NumNodes, j + NumNodes) = lhs_total(i, j);
        }
    }

    BoundedVector<double, 2*NumNodes> split_element_values;
    split_element_values = PotentialFlowUtilities::GetPotentialOnWakeElement<Dim, NumNodes>(*this, data.distances);
    noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix, split_element_values);
}

template <>
ModifiedShapeFunctions::Pointer EmbeddedIncompressiblePotentialFlowElement<2,3>::pGetModifiedShapeFunctions(Vector& rDistances) {
    return Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(this->pGetGeometry(), rDistances);
}

template <>
ModifiedShapeFunctions::Pointer EmbeddedIncompressiblePotentialFlowElement<3,4>::pGetModifiedShapeFunctions(Vector& rDistances) {
    return Kratos::make_shared<Tetrahedra3D4ModifiedShapeFunctions>(this->pGetGeometry(), rDistances);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Inquiry

template <int Dim, int NumNodes>
int EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

    // Generic geometry check
    int out = BaseType::Check(rCurrentProcessInfo);
    if (out != 0)
    {
        return out;
    }

    for (unsigned int i = 0; i < this->GetGeometry().size(); i++)
    {
        KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(GEOMETRY_DISTANCE,this->GetGeometry()[i]);
    }

    return out;

    KRATOS_CATCH("");
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Input and output

template <int Dim, int NumNodes>
std::string EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::Info() const
{
    std::stringstream buffer;
    buffer << "EmbeddedIncompressiblePotentialFlowElement #" << this->Id();
    return buffer.str();
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::PrintInfo(std::ostream& rOStream) const
{
    rOStream << "EmbeddedIncompressiblePotentialFlowElement #" << this->Id();
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::PrintData(std::ostream& rOStream) const
{
    this->pGetGeometry()->PrintData(rOStream);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Private functions
///////////////////////////////////////////////////////////////////////////////////////////////////

// serializer

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element);
}

template <int Dim, int NumNodes>
void EmbeddedIncompressiblePotentialFlowElement<Dim, NumNodes>::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Template class instantiation

template class EmbeddedIncompressiblePotentialFlowElement<2, 3>;
template class EmbeddedIncompressiblePotentialFlowElement<3, 4>;


} // namespace Kratos
