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

// System includes


// External includes


// Include Base h
#include "includes/checks.h"
#include "includes/variables.h"
#include "utilities/geometry_utilities.h"
#include "custom_elements/electric_element.h"


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

/**
 * Constructor.
 */
ElectricElement::ElectricElement(IndexType NewId)
    : Element(NewId) 
{
}

/**
 * Constructor using an array of nodes
 */
ElectricElement::ElectricElement(IndexType NewId, const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes) 
{
}

/**
 * Constructor using Geometry
 */
ElectricElement::ElectricElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) 
{
}

/**
 * Constructor using Properties
 */
ElectricElement::ElectricElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) 
{
}

/**
 * Copy Constructor
 */
ElectricElement::ElectricElement(ElectricElement const& rOther)
    : Element(rOther) 
{
}

/**
 * Destructor
 */
ElectricElement::~ElectricElement()
{
}

///@}
///@name Operators
///@{

/// Assignment operator.
ElectricElement & ElectricElement::operator=(ElectricElement const& rOther)
{
    BaseType::operator=(rOther);
    Flags::operator =(rOther);
    // mpProperties = rOther.mpProperties;
    return *this;
}

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
Element::Pointer ElectricElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<ElectricElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer ElectricElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<ElectricElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer and clones the previous element data
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer ElectricElement::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<ElectricElement>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

/**
 * this determines the elemental equation ID vector for all elemental
 * DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if (rResult.size() != number_of_nodes)
        rResult.resize(number_of_nodes, false);

      for (unsigned int i = 0; i < number_of_nodes; i++)
    rResult[i] = GetGeometry()[i].GetDof(EPOTENCIAL).EquationId();


}

/**
 * determines the elemental list of DOFs
 * @param ElementalDofList: the list of DOFs
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricElement::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo) const
{
    unsigned int number_of_nodes = GetGeometry().PointsNumber();
    if (rElementalDofList.size() != number_of_nodes)
        rElementalDofList.resize(number_of_nodes);

      for (unsigned int i = 0; i < number_of_nodes; i++)
    rElementalDofList[i] = GetGeometry()[i].pGetDof(EPOTENCIAL);


}

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
void ElectricElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
KRATOS_TRY


    bounded_matrix<double,3,2> DN_DX= ZeroMatrix(3,2);  // Gradients matrix 
    bounded_matrix<double,2,2> D= ZeroMatrix(2,2);  
     array_1d<double,3> N; //dimension = number of nodes . Position of the gauss point 
    array_1d<double,3> temp; //dimension = number of nodes . . since we are using a residualbased approach 
    array_1d<double,3> surface_sources;
    array_1d<double,3> Output;

		const unsigned int number_of_points = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    //resizing as needed the LHS
		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false);
			
    //resizing as needed the RHS		
			

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points,false);
   
		//getting data for the given geometry


		double area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area);

		//reading properties and conditions
		double permittivity = GetProperties()[PERMITTIVITYNEG];
		D(0,0)=permittivity;
		D(1,1)=permittivity;
	    surface_sources[0] = (this)->GetValue(SCHARGE);
		surface_sources[1] = (this)->GetValue(SCHARGE);
		surface_sources[2] = (this)->GetValue(SCHARGE);

		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(this->GetIntegrationMethod());
        KRATOS_WATCH(DN_DX)
  	    noalias(rLeftHandSideMatrix) = prod(DN_DX,Matrix(prod(D,trans(DN_DX))));
		rLeftHandSideMatrix *= area;
		noalias(rRightHandSideVector) = surface_sources*area/3.0; //to get the continuty the voltage 

	
		for(unsigned int iii = 0; iii<number_of_points; iii++)
			temp[iii] = GetGeometry()[iii].FastGetSolutionStepValue(EPOTENCIAL);
  
			noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp);
		        KRATOS_CATCH("");
}

void ElectricElement::CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, std::vector<array_1d<double,3> >& Output, const ProcessInfo& rCurrentProcessInfo)
	{
		bounded_matrix<double,3,2> DN_DX;  // Gradients matrix 
    		bounded_matrix<double,2,2> D;      // Conductivity matrix 
    		D = ZeroMatrix(2,2); //initializing the matrix as zero
		array_1d<double,3> N; //dimension = number of nodes . Position of the gauss point 		

		
		IntegrationMethod mThisIntegrationMethod;

        mThisIntegrationMethod= GetGeometry().GetDefaultIntegrationMethod();
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints(mThisIntegrationMethod);
		
        bounded_matrix<double,3,3> I= ZeroMatrix(3,3);

        Vector vect_tmp(3);
        Vector efield_tmp(3);
		vect_tmp[0] = GetGeometry()[0].FastGetSolutionStepValue(EPOTENCIAL);
		vect_tmp[1] = GetGeometry()[1].FastGetSolutionStepValue(EPOTENCIAL);
		vect_tmp[2] = GetGeometry()[2].FastGetSolutionStepValue(EPOTENCIAL);

		//reading properties and conditions
		double permittivity = GetProperties()[PERMITTIVITYNEG];

		double area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, area);

        efield_tmp = -prod(trans(DN_DX),vect_tmp);
        //auto& efield_tmp = this->SetValue(EDISPLACEMENT)

        I(0,0)=efield_tmp[0]*efield_tmp[0];
        I(1,1)=efield_tmp[1]*efield_tmp[1];

		if(Output.size() != GetGeometry().IntegrationPoints(mThisIntegrationMethod).size())
            Output.resize(GetGeometry().IntegrationPoints(mThisIntegrationMethod).size());

		
        for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
		{
			if (rVariable == EFIELD){
            noalias(Output[PointNumber])= efield_tmp;    
			//noalias(Output[PointNumber])=-prod(trans(DN_DX),vect_tmp);
			Output[PointNumber][2]=0.0;
            }
        /*    else if (rVariable == EDISPLACEMENT){
            noalias(Output[PointNumber])= -0.5 * permittivity * area * prod(I,trans(N));     
			//noalias(Output[PointNumber])=-prod(trans(DN_DX),vect_tmp);
			Output[PointNumber][2]=0.0;
            }*/
		}
 			
	}

void ElectricElement::GetValueOnIntegrationPoints(const Variable<array_1d<double,3> >& rVariable, 
		 std::vector<array_1d<double,3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
    {
		//KRATOS_WATCH("GiD Post Electrostatic - GetValueOnIntegrationPoints");

        if(rVariable==EFIELD)
        {
//			KRATOS_WATCH("GiD Post Electrostatic - get - elec-field");
            CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        }
        else if(rVariable==EDISPLACEMENT)
        {
//			KRATOS_WATCH("GiD Post Electrostatic - get - elec-field");
            CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
        }
    }





















/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix only
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector only
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the first derivatives contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricElement::CalculateFirstDerivativesContributions(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix for the first derivatives constributions
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricElement::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector for the first derivatives constributions
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricElement::CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

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
void ElectricElement::CalculateSecondDerivativesContributions(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix for the second derivatives constributions
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricElement::CalculateSecondDerivativesLHS(
    MatrixType& rLeftHandSideMatrix,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rLeftHandSideMatrix.size1() != 0)
        rLeftHandSideMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector for the second derivatives constributions
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricElement::CalculateSecondDerivativesRHS(
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    if (rRightHandSideVector.size() != 0)
        rRightHandSideVector.resize(0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental mass matrix
 * @param rMassMatrix: the elemental mass matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if (rMassMatrix.size1() != 0)
        rMassMatrix.resize(0, 0, false);
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental damping matrix
 * @param rDampingMatrix: the elemental damping matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
{
    if (rDampingMatrix.size1() != 0)
        rDampingMatrix.resize(0, 0, false);
}

/**
 * This method provides the place to perform checks on the completeness of the input
 * and the compatibility with the problem options as well as the contitutive laws selected
 * It is designed to be called only once (or anyway, not often) typically at the beginning
 * of the calculations, so to verify that nothing is missing from the input
 * or that no common error is found.
 * @param rCurrentProcessInfo
 * this method is: MANDATORY
 */
int ElectricElement::Check(const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_TRY

       // Base class checks for positive Jacobian and Id > 0
      int ierr = Element::Check(rCurrentProcessInfo);
      if(ierr != 0) return ierr;
  
      // Check that all required variables have been registered
      KRATOS_CHECK_VARIABLE_KEY(EPOTENCIAL)
 
      unsigned const int number_of_points = GetGeometry().size();  //added cornejo
      // Check that the element's nodes contain all required SolutionStepData and Degrees of freedom
      for ( unsigned int i = 0; i < number_of_points; i++ )
      {
          auto &rnode = this->GetGeometry()[i];
          KRATOS_CHECK_VARIABLE_IN_NODAL_DATA(EPOTENCIAL,rnode)
          KRATOS_CHECK_DOF_IN_NODE(EPOTENCIAL,rnode)
      }
  
      return ierr;
  
      KRATOS_CATCH("");
}

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

std::string ElectricElement::Info() const {
    std::stringstream buffer;
    buffer << "ElectricElement #" << Id();
    return buffer.str();
}

/// Print information about this object.

void ElectricElement::PrintInfo(std::ostream& rOStream) const {
    rOStream << "ElectricElement #" << Id();
}

/// Print object's data.

void ElectricElement::PrintData(std::ostream& rOStream) const {
    pGetGeometry()->PrintData(rOStream);
}

///@}
///@name Friends
///@{

///@}

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

void ElectricElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );

    // List
    // To be completed with the class member list
}

void ElectricElement::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Element );

    // List
    // To be completed with the class member list
}

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
///@name Type Definitions
///@{

///@}
///@name Input and output
///@{

/// input stream function
inline std::istream & operator >> (std::istream& rIStream, ElectricElement& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const ElectricElement& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.
