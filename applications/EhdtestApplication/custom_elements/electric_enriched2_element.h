//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    , KratosAppGenerator
//

#if !defined(KRATOS_ELECTRIC_ENRICHED2_ELEMENT_H_INCLUDED )
#define KRATOS_ELECTRIC_ENRICHED2_ELEMENT_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/element.h"
#include "ehdtest_application_variables.h"
#include "modified_shape_functions/triangle_2d_3_modified_shape_functions.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{


class ElectricEnriched2Element : public Element
{
public:

    

    ///@name Type Definitions
    ///@{

    typedef Element BaseType;

    ///@}
    ///@name Pointer Definitions
    /// Pointer definition of ElectricEnrichedElement
    KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(ElectricEnriched2Element);

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * Constructor.
     */
    ElectricEnriched2Element(IndexType NewId = 0);

    /**
     * Constructor using an array of nodes
     */
    ElectricEnriched2Element(IndexType NewId, const NodesArrayType& ThisNodes);

    /**
     * Constructor using Geometry
     */
    ElectricEnriched2Element(IndexType NewId, GeometryType::Pointer pGeometry);

    /**
     * Constructor using Properties
     */
    ElectricEnriched2Element(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

    /**
     * Copy Constructor
     */
    ElectricEnriched2Element(ElectricEnriched2Element const& rOther);

    /**
     * Destructor
     */
    ~ElectricEnriched2Element() override;

    ///@}
    ///@name Operators
    ///@{

    /// Assignment operator.
    ElectricEnriched2Element & operator=(ElectricEnriched2Element const& rOther);

    ///@}
    ///@name Operations
    ///@{

    /**
     * ELEMENTS inherited from this class have to implement next
     * Create and Clone methods: MANDATORY
     */

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

    /**
     * creates a new element pointer
     * @param NewId: the ID of the new element
     * @param pGeom: the geometry to be employed
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Create(IndexType NewId, GeometryType::Pointer pGeom, PropertiesType::Pointer pProperties) const override;

    /**
     * creates a new element pointer and clones the previous element data
     * @param NewId: the ID of the new element
     * @param ThisNodes: the nodes of the new element
     * @param pProperties: the properties assigned to the new element
     * @return a Pointer to the new element
     */
    Element::Pointer Clone(IndexType NewId, NodesArrayType const& ThisNodes) const override;

    /**
     * this determines the elemental equation ID vector for all elemental
     * DOFs
     * @param rResult: the elemental equation ID vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const override;

    /**
     * determines the elemental list of DOFs
     * @param ElementalDofList: the list of DOFs
     * @param rCurrentProcessInfo: the current process info instance
     */
    void GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo) const override;

    /**
     * ELEMENTS inherited from this class have to implement next
     * CalculateLocalSystem, CalculateLeftHandSide and CalculateRightHandSide methods
     * they can be managed internally with a private method to do the same calculations
     * only once: MANDATORY
     */

    /**
     * this is called during the assembling process in order
     * to calculate all elemental contributions to the global system
     * matrix and the right hand side
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateLocalSystem(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix only
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
     
     
    
    void CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, 
    std::vector<array_1d<double,3> >& Output, const ProcessInfo& rCurrentProcessInfo);

    void GetValueOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, 
    std::vector<array_1d<double,3> >& rValues, const ProcessInfo& rCurrentProcessInfo);
    

    void Calculate(const Variable<array_1d<double, 3 > >& rVariable,
    array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo);
   
     
    void CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector only
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the first derivatives contributions for the LHS and RHS
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateFirstDerivativesContributions(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix for the first derivatives constributions
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector for the first derivatives constributions
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * ELEMENTS inherited from this class must implement this methods
     * if they need to add dynamic element contributions
     * note: second derivatives means the accelerations if the displacements are the dof of the analysis
     * note: time integration parameters must be set in the rCurrentProcessInfo before calling these methods
     * CalculateSecondDerivativesContributions,
     * CalculateSecondDerivativesLHS, CalculateSecondDerivativesRHS methods are : OPTIONAL
     */


    /**
     * this is called during the assembling process in order
     * to calculate the second derivative contributions for the LHS and RHS
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rRightHandSideVector: the elemental right hand side
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateSecondDerivativesContributions(
        MatrixType& rLeftHandSideMatrix,
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental left hand side matrix for the second derivatives constributions
     * @param rLeftHandSideMatrix: the elemental left hand side matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateSecondDerivativesLHS(
        MatrixType& rLeftHandSideMatrix,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental right hand side vector for the second derivatives constributions
     * @param rRightHandSideVector: the elemental right hand side vector
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateSecondDerivativesRHS(
        VectorType& rRightHandSideVector,
        const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental mass matrix
     * @param rMassMatrix: the elemental mass matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * this is called during the assembling process in order
     * to calculate the elemental damping matrix
     * @param rDampingMatrix: the elemental damping matrix
     * @param rCurrentProcessInfo: the current process info instance
     */
    void CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo) override;

    /**
     * This method provides the place to perform checks on the completeness of the input
     * and the compatibility with the problem options as well as the contitutive laws selected
     * It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input
     * or that no common error is found.
     * @param rCurrentProcessInfo
     * this method is: MANDATORY
     */
    int Check(const ProcessInfo& rCurrentProcessInfo) const override;

    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override;

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    void PrintData(std::ostream& rOStream) const override;

    ///@}
    ///@name Friends
    ///@{

    ///@}

protected:

    ///@name Protected static Member Variables
    ///@{

    ///@}
    ///@name Protected member Variables
    ///@{

    ///@}
    ///@name Protected Operators
    ///@{

    ///@}
    ///@name Protected Operations
    ///@{

    ///@}
    ///@name Protected  Access
    ///@{

    ///@}
    ///@name Protected Inquiry
    ///@{

    ///@}
    ///@name Protected LifeCycle
    ///@{

    ///@}

private:

    ///@name Static Member Variables
    ///@{

    

    ///@}
    ///@name Member Variables
    ///@{

    

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    ///@}
    ///@name Serialization
    ///@{

    friend class Serializer;

    void save(Serializer& rSerializer) const override;
    void load(Serializer& rSerializer) override;

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    ///@}

}; // Class ElectricEnrichedElement

///@}

///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

///@}

} // namespace Kratos.

#endif // KRATOS_ENRICHED2_ELEMENT_H_INCLUDED  defined
