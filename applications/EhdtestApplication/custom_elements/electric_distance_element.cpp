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
#include "custom_elements/electric_distance_element.h"


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
ElectricDistanceElement::ElectricDistanceElement(IndexType NewId)
    : Element(NewId) 
{
}

/**
 * Constructor using an array of nodes
 */
ElectricDistanceElement::ElectricDistanceElement(IndexType NewId, const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes) 
{
}

/**
 * Constructor using Geometry
 */
ElectricDistanceElement::ElectricDistanceElement(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) 
{
}

/**
 * Constructor using Properties
 */
ElectricDistanceElement::ElectricDistanceElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) 
{
}

/**
 * Copy Constructor
 */
ElectricDistanceElement::ElectricDistanceElement(ElectricDistanceElement const& rOther)
    : Element(rOther) 
{
}

/**
 * Destructor
 */
ElectricDistanceElement::~ElectricDistanceElement()
{
}

///@}
///@name Operators
///@{

/// Assignment operator.
ElectricDistanceElement & ElectricDistanceElement::operator=(ElectricDistanceElement const& rOther)
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
Element::Pointer ElectricDistanceElement::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<ElectricDistanceElement>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer ElectricDistanceElement::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<ElectricDistanceElement>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer and clones the previous element data
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer ElectricDistanceElement::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<ElectricDistanceElement>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

/**
 * this determines the elemental equation ID vector for all elemental
 * DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricDistanceElement::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
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
void ElectricDistanceElement::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo) const
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
void ElectricDistanceElement::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY

    const int num_dim  = 2;
    const int num_nodes  = num_dim + 1;
    const int num_dof = num_nodes; 

    GeometryData::ShapeFunctionsGradientsType DN_DX;
    Matrix N;
    Vector DetJ;
    Vector weights;

    Vector distance(num_dof);
    Vector values(num_dof);

    BoundedMatrix<double,num_dof,num_dof> lhs = ZeroMatrix(num_dof,num_dof);
    Vector rhs = ZeroVector(num_dof);
    array_1d<double,3> surface_sources;

    array_1d<double,3> efieldpos; //dimension = number of nodes . Position of the gauss point 
    array_1d<double,3> efieldneg; //dimension = number of nodes . Position of the gauss point

    /////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////////
    
    const double permittivity_pos = GetProperties()[PERMITTIVITYPOS];
    const double permittivity_neg = GetProperties()[PERMITTIVITYNEG];

    surface_sources[0] = (this)->GetValue(SCHARGE);
	surface_sources[1] = (this)->GetValue(SCHARGE);
	surface_sources[2] = (this)->GetValue(SCHARGE);
    ///////////////////////////////////////////////////////////////////////////////////
    const GeometryData::IntegrationMethod integration_method = GeometryData::GI_GAUSS_1;
    GeometryType::Pointer p_geometry = this->pGetGeometry();

    unsigned int nneg=0, npos=0;

    for (unsigned int i_node=0; i_node < num_nodes; ++i_node){

        const double dist = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE);
        distance[i_node] = dist;
        values[i_node] = (*p_geometry)[i_node].FastGetSolutionStepValue(EPOTENCIAL);

        if (dist == 0.0){
            (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE) = 1.0e-12;
        }

        if (dist > 0.0) npos += 1;
        else nneg += 1;
    }
      
    if(rLeftHandSideMatrix.size1() != num_dof)
        rLeftHandSideMatrix.resize(num_dof,num_dof,false); //resizing the system in case it does not have the right size 

    if(rRightHandSideVector.size() != num_dof)
        rRightHandSideVector.resize(num_dof,false);
     
    if (nneg == num_nodes || npos == num_nodes){
        
        const unsigned int number_of_gauss_points = p_geometry->IntegrationPointsNumber(integration_method);
        
        // Getting data for the given geometry
       
        p_geometry->ShapeFunctionsIntegrationPointsGradients(DN_DX, DetJ, integration_method);
         
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = p_geometry->IntegrationPoints(integration_method);
        if (weights.size() != number_of_gauss_points) {
            weights.resize(number_of_gauss_points,false);
        }
        //calculating the area of the elementes
        for (unsigned int gp = 0; gp < number_of_gauss_points; gp++){
            weights[gp] = DetJ[gp] * IntegrationPoints[gp].Weight();
        }
            

        if (npos == 3){
           lhs+= permittivity_pos*weights(0)*prod(DN_DX[0],trans(DN_DX[0]));
        }else {
           lhs+= permittivity_neg*weights(0)*prod(DN_DX[0],trans(DN_DX[0]));
        }
        rhs+= surface_sources*weights(0)/3.0;
        //KRATOS_WATCH("FINISHED COMPUTING NON-CUT ELEMENTS...........")
    } else { //(npos != 0 && nneg != 0)//if (npos > 0 )
        ModifiedShapeFunctions::Pointer p_modified_sh_func =
            Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geometry, distance);
        Matrix pos_N;
        GeometryType::ShapeFunctionsGradientsType pos_DN_DX;
        Vector pos_weights;

        p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
        pos_N,        // N
        pos_DN_DX,    // DN_DX
        pos_weights,  // weight * detJ
        integration_method); 
        const std::size_t number_of_pos_gauss_points = pos_weights.size();
        for (unsigned int pos_gp = 0; pos_gp < number_of_pos_gauss_points; pos_gp++){
            
            //lhs+= permittivity_pos*pos_weights(pos_gp)*prod(pos_DN_DX[0],trans(pos_DN_DX[0]));
            lhs+= permittivity_pos*pos_weights(pos_gp)*prod(pos_DN_DX[0],trans(pos_DN_DX[0]));         
        }     

             
        Matrix neg_N;    
        GeometryType::ShapeFunctionsGradientsType neg_DN_DX;
        Vector neg_weights;
      
        p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
            neg_N,        // N
            neg_DN_DX,    // DN_DX
            neg_weights,  // weight * detJ
            integration_method);
        const std::size_t number_of_neg_gauss_points = neg_weights.size();

            for (unsigned int neg_gp = 0; neg_gp < number_of_neg_gauss_points; neg_gp++){
            
            //lhs+= permittivity_neg*neg_weights(neg_gp)*prod(neg_DN_DX[0],trans(neg_DN_DX[0]));
            lhs+= permittivity_neg*neg_weights(neg_gp)*prod(neg_DN_DX[0],trans(neg_DN_DX[0]));
            }    
        
        
        //KRATOS_WATCH("SHOULD NOT ENETER HERE IF NOT CUT ELEMENTS ARE PRESENT..........")

    }
    
    noalias(rLeftHandSideMatrix) = lhs;
    noalias(rRightHandSideVector) = rhs - prod(rLeftHandSideMatrix,values); 
    
    if (npos == num_nodes || nneg == num_nodes){
    }else    {
        KRATOS_WATCH(rLeftHandSideMatrix)
        KRATOS_WATCH(rRightHandSideVector)
    }
    KRATOS_CATCH("");
    
}

void ElectricDistanceElement::GetValueOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, 
std::vector<array_1d<double,3> >& rValues, const ProcessInfo& rCurrentProcessInfo)
{
 
            if(rVariable==EFIELD)
            {
            CalculateOnIntegrationPoints (rVariable, rValues, rCurrentProcessInfo);    
            }  
            else if(rVariable==EFORCEPOS)
            {
            CalculateOnIntegrationPoints (rVariable, rValues, rCurrentProcessInfo);    
            }  
            else if(rVariable==EFORCENEG)
            {
            CalculateOnIntegrationPoints (rVariable, rValues, rCurrentProcessInfo);    
            } 

}


///////************************* ELECTRIC FORCE CALCULATION ************//////////////////////////////

void ElectricDistanceElement::CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, 
std::vector<array_1d<double,3> >& Output, const ProcessInfo& rCurrentProcessInfo)
{ 
 
    const int num_dim  = 2;
    const int num_nodes  = num_dim + 1;
    const int num_dof = num_nodes; 

    Vector Fepos(num_dof);
    Vector Feneg(num_dof);

    Vector Ne(num_dof);
    Ne[0] = 1;
    Ne[1] = 1;
    Ne[2] = 1;  
   

    GeometryData::ShapeFunctionsGradientsType DN_DX;
    Matrix N;
    Vector DetJ;
    Vector weights;

    Vector distance(num_dof);
    Vector values(num_dof);

    //Vector efield_pos;
    //Vector efield_neg;

    bounded_matrix<double,3,3> I= ZeroMatrix(3,3);
    bounded_matrix<double,3,3> Ipos= ZeroMatrix(3,3);
    bounded_matrix<double,3,3> Ineg= ZeroMatrix(3,3);
    

    array_1d<double,3> efield_pos = ZeroVector(3); //dimension = number of nodes . Position of the gauss point 
    array_1d<double,3> efield_neg = ZeroVector(3); //dimension = number of nodes . Position of the gauss point 
    array_1d<double,3> efield = ZeroVector(3); //dimension = number of nodes . Position of the gauss point 
    /////////////////////////////////////////////////////////////////////////////////
    const double permittivity_pos = GetProperties()[PERMITTIVITYPOS];
    const double permittivity_neg = GetProperties()[PERMITTIVITYNEG];
    //////////////////////////////////////////////////////////////////////////////////

    const GeometryData::IntegrationMethod integration_method = GeometryData::GI_GAUSS_1;
    GeometryType::Pointer p_geometry = this->pGetGeometry();

    unsigned int nneg=0, npos=0;
      
    for (unsigned int i_node=0; i_node < num_nodes; ++i_node){

        const double dist = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE);
        distance[i_node] = dist;
        values[i_node] = (*p_geometry)[i_node].FastGetSolutionStepValue(EPOTENCIAL);

        if (dist == 0.0){
            (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE) = 1.0e-12;
        }

        if (dist > 0.0) npos += 1;
        else nneg += 1;
    }
   
    if (nneg == num_nodes || npos == num_nodes){
            
            //Vector efield; //dimension = number of nodes . Position of the gauss point 
            const unsigned int number_of_gauss_points = p_geometry->IntegrationPointsNumber(integration_method);
        
            p_geometry->ShapeFunctionsIntegrationPointsGradients(DN_DX, DetJ, integration_method);
            const GeometryType::IntegrationPointsArrayType& IntegrationPoints = p_geometry->IntegrationPoints(integration_method);
            if (weights.size() != number_of_gauss_points) {
                weights.resize(number_of_gauss_points,false);
            }
        //calculating the area of the elementes
            for (unsigned int gp = 0; gp < number_of_gauss_points; gp++){
                weights[gp] = DetJ[gp] * IntegrationPoints[gp].Weight();
            }
            ;                  
            efield= -prod(trans(DN_DX[0]),values);
          
            I(0,0)=efield[0]*efield[0];
            I(1,1)=efield[1]*efield[1];

            //if (npos == 3){
            //    Fepos =  -0.1667 * permittivity_pos * weights(0)*prod(I,trans(Ne)); 
            //} else {
            //    Feneg =  -0.1667 * permittivity_neg * weights(0)*prod(I,trans(Ne)); 
            // }

    } else  { 
            ModifiedShapeFunctions::Pointer p_modified_sh_func =
                Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geometry, distance);
        
            Matrix pos_N;
            GeometryType::ShapeFunctionsGradientsType pos_DN_DX;
            Vector pos_weights;
            //Vector efield_pos;
            
            p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
            pos_N,        // N
            pos_DN_DX,    // DN_DX
            pos_weights,  // weight * detJ
            integration_method);              
            
            //CALCULATION THE ELECTRIC FIELD ON THE SPITED ELEMENTS /////

            efield_pos = -prod(trans(pos_DN_DX[0]),values); //the values at the interface have two DN_DX

            Ipos(0,0)=efield_pos[0]*efield_pos[0];
            Ipos(1,1)=efield_pos[1]*efield_pos[1]; 
            //compute the electric force 
            //const std::size_t number_of_pos_gauss_points = pos_weights.size();
            //for (unsigned int pos_gp = 0; pos_gp < number_of_pos_gauss_points; pos_gp++){
            //Fepos =  -0.1667 * permittivity_pos * pos_weights(pos_gp)*prod(Ipos,trans(Ne));
            //}    
            Matrix neg_N;
            GeometryType::ShapeFunctionsGradientsType neg_DN_DX;
            Vector neg_weights;
            //Vector efield_neg;
            
            p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
            neg_N,        // N
            neg_DN_DX,    // DN_DX
            neg_weights,  // weight * detJ
            integration_method);

            efield_neg = -prod(trans(neg_DN_DX[0]),values);

            //KRATOS_WATCH(this->Id())
            //Ineg(0,0)=efield_neg[0]*efield_neg[0];
            //Ineg(1,1)=efield_neg[1]*efield_neg[1]; 

            //const std::size_t number_of_neg_gauss_points = neg_weights.size();
            //for (unsigned int neg_gp = 0; neg_gp < number_of_neg_gauss_points; neg_gp++){
            //Feneg =  -0.1667 * permittivity_neg * neg_weights(neg_gp)*prod(Ineg,trans(Ne));
            //}    

        }
    
    const unsigned int number_of_gauss_points = p_geometry->IntegrationPointsNumber(integration_method);
    for (unsigned int i = 0; i < number_of_gauss_points; i++ ) {
    
        if(rVariable == EFORCEPOS)//electric force form the negative side 
        {
                noalias(Output[i]) = Fepos;
                Output[i][2]=0.0;
        }
        else if(rVariable == EFORCENEG)//electric force form the negative side 
        {
                noalias(Output[i]) = Feneg;
                Output[i][2]=0.0;
        }
        else if(rVariable == EFIELD)//electric force form the negative side 
        {
                noalias(Output[i]) = efield + efield_pos;
                Output[i][2]=0.0;
        }
        
    }     

}

void ElectricDistanceElement::Calculate(const Variable<array_1d<double, 3 > >& rVariable,
array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo)
{
    
   
    const int num_dim  = 2;
    const int num_nodes  = num_dim + 1;
    const int num_dof = num_nodes; 

    Vector Fepos(num_dof);
    Vector Feneg(num_dof);

    Vector Ne(num_dof);
    Ne[0] = 1;
    Ne[1] = 1;
    Ne[2] = 1;  
   

    GeometryData::ShapeFunctionsGradientsType DN_DX;
    Matrix N;
    Vector DetJ;
    Vector weights;

    Vector distance(num_dof);
    Vector values(num_dof);

    //Vector efield_pos;
    //Vector efield_neg;


    bounded_matrix<double,3,3> Ipos= ZeroMatrix(3,3);
    bounded_matrix<double,3,3> Ineg= ZeroMatrix(3,3);

    //array_1d<double,3> efield_pos; //dimension = number of nodes . Position of the gauss point 
    //array_1d<double,3> efield_neg; //dimension = number of nodes . Position of the gauss point 
    /////////////////////////////////////////////////////////////////////////////////
    const double permittivity_pos = GetProperties()[PERMITTIVITYPOS];
    const double permittivity_neg = GetProperties()[PERMITTIVITYNEG];
    //////////////////////////////////////////////////////////////////////////////////

    const GeometryData::IntegrationMethod integration_method = GeometryData::GI_GAUSS_2;
    GeometryType::Pointer p_geometry = this->pGetGeometry();

    unsigned int nneg=0, npos=0;
      
    for (unsigned int i_node=0; i_node < num_nodes; ++i_node){

        const double dist = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE);
        distance[i_node] = dist;
        values[i_node] = (*p_geometry)[i_node].FastGetSolutionStepValue(EPOTENCIAL);

        if (dist == 0.0){
            (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE) = 1.0e-12;
        }

        if (dist > 0.0) npos += 1;
        else nneg += 1;
    }
   
    if (nneg == num_nodes || npos == num_nodes){
            Vector efield_pos; //dimension = number of nodes . Position of the gauss point 
            const unsigned int number_of_gauss_points = p_geometry->IntegrationPointsNumber(integration_method);
        
            p_geometry->ShapeFunctionsIntegrationPointsGradients(DN_DX, DetJ, integration_method);
            const GeometryType::IntegrationPointsArrayType& IntegrationPoints = p_geometry->IntegrationPoints(integration_method);
            if (weights.size() != number_of_gauss_points) {
                weights.resize(number_of_gauss_points,false);
            }
        //calculating the area of the elementes
            for (unsigned int gp = 0; gp < number_of_gauss_points; gp++){
                weights[gp] = DetJ[gp] * IntegrationPoints[gp].Weight();
            }
                   
            efield_pos = -prod(trans(DN_DX[0]),values);
            
            Ipos(0,0)=efield_pos[0]*efield_pos[0];
            Ipos(1,1)=efield_pos[1]*efield_pos[1];
            Fepos =  -0.1667 * permittivity_pos * weights(0)*prod(Ipos,trans(Ne)); 
    } else  { 
            ModifiedShapeFunctions::Pointer p_modified_sh_func =
                Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geometry, distance);
        
            Matrix pos_N;
            GeometryType::ShapeFunctionsGradientsType pos_DN_DX;
            Vector pos_weights;
            Vector efield_pos;
            
            p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
            pos_N,        // N
            pos_DN_DX,    // DN_DX
            pos_weights,  // weight * detJ
            integration_method);   
 
            efield_pos = -prod(trans(pos_DN_DX[0]),values); //the values at the interface have two DN_DX
            Ipos(0,0)=efield_pos[0]*efield_pos[0];
            Ipos(1,1)=efield_pos[1]*efield_pos[1]; 
            //compute the electric force 

            const std::size_t number_of_pos_gauss_points = pos_weights.size();
            for (unsigned int pos_gp = 0; pos_gp < number_of_pos_gauss_points; pos_gp++){
            Fepos =  -0.1667 * permittivity_pos * pos_weights(pos_gp)*prod(Ipos,trans(Ne));
            }    
            Matrix neg_N;
            GeometryType::ShapeFunctionsGradientsType neg_DN_DX;
            Vector neg_weights;
            Vector efield_neg;
      
            p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
            neg_N,        // N
            neg_DN_DX,    // DN_DX
            neg_weights,  // weight * detJ
            integration_method);
            /////// *****  NEW ******
            //auto& neighbour_elems = this->GetValue(NEIGHBOUR_ELEMENTS);
            

            ////// ****** END NEW ******
            efield_neg = -prod(trans(neg_DN_DX[0]),values);
            Ineg(0,0)=efield_neg[0]*efield_neg[0];
            Ineg(1,1)=efield_neg[1]*efield_neg[1]; 

            const std::size_t number_of_neg_gauss_points = neg_weights.size();
            for (unsigned int neg_gp = 0; neg_gp < number_of_neg_gauss_points; neg_gp++){
            Feneg =  -0.1667 * permittivity_neg * neg_weights(neg_gp)*prod(Ineg,trans(Ne));
            }    
            
        }
    
    if(rVariable == EFORCEPOS) //electric force form the positive side 
    {    
        GetGeometry().GetValue(EFORCEPOS) = Fepos;
        //Output[i][2]=0.0;
    }
    else if(rVariable == EFORCENEG)//electric force form the negative side 
    {
        GetGeometry().GetValue(EFORCENEG) = Feneg;
        //Output[i][2]=0.0;
    }
}


/**
 * this is called during the assembling process in order
 * to calculate the elemental left hand side matrix only
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricDistanceElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector only
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricDistanceElement::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the first derivatives contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricDistanceElement::CalculateFirstDerivativesContributions(
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
void ElectricDistanceElement::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
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
void ElectricDistanceElement::CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
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
void ElectricDistanceElement::CalculateSecondDerivativesContributions(
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
void ElectricDistanceElement::CalculateSecondDerivativesLHS(
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
void ElectricDistanceElement::CalculateSecondDerivativesRHS(
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
void ElectricDistanceElement::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
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
void ElectricDistanceElement::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
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
int ElectricDistanceElement::Check(const ProcessInfo& rCurrentProcessInfo) const
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

std::string ElectricDistanceElement::Info() const {
    std::stringstream buffer;
    buffer << "ElectricDistanceElement #" << Id();
    return buffer.str();
}

/// Print information about this object.

void ElectricDistanceElement::PrintInfo(std::ostream& rOStream) const {
    rOStream << "ElectricDistanceElement #" << Id();
}

/// Print object's data.

void ElectricDistanceElement::PrintData(std::ostream& rOStream) const {
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

void ElectricDistanceElement::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );

    // List
    // To be completed with the class member list
}

void ElectricDistanceElement::load(Serializer& rSerializer)
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
inline std::istream & operator >> (std::istream& rIStream, ElectricDistanceElement& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const ElectricDistanceElement& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.
