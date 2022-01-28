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
#include "custom_elements/electric_enriched2_element.h"
#include "utilities/math_utils.h"
#include "utilities/enrichment_utilities.h"
#include "utilities/enrich_2d_2dofs.h"


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
ElectricEnriched2Element::ElectricEnriched2Element(IndexType NewId)
    : Element(NewId) 
{
}

/**
 * Constructor using an array of nodes
 */
ElectricEnriched2Element::ElectricEnriched2Element(IndexType NewId, const NodesArrayType& ThisNodes)
    : Element(NewId, ThisNodes) 
{
}

/**
 * Constructor using Geometry
 */
ElectricEnriched2Element::ElectricEnriched2Element(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) 
{
}

/**
 * Constructor using Properties
 */
ElectricEnriched2Element::ElectricEnriched2Element(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) 
{
}

/**
 * Copy Constructor
 */
ElectricEnriched2Element::ElectricEnriched2Element(ElectricEnriched2Element const& rOther)
    : Element(rOther) 
{
}

/**
 * Destructor
 */
ElectricEnriched2Element::~ElectricEnriched2Element()
{
}

///@}
///@name Operators
///@{

/// Assignment operator.
ElectricEnriched2Element & ElectricEnriched2Element::operator=(ElectricEnriched2Element const& rOther)
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
Element::Pointer ElectricEnriched2Element::Create(
    IndexType NewId,
    NodesArrayType const& ThisNodes,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<ElectricEnriched2Element>(NewId, GetGeometry().Create(ThisNodes), pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer
 * @param NewId: the ID of the new element
 * @param pGeom: the geometry to be employed
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer ElectricEnriched2Element::Create(
    IndexType NewId,
    GeometryType::Pointer pGeom,
    PropertiesType::Pointer pProperties) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<ElectricEnriched2Element>(NewId, pGeom, pProperties);
    KRATOS_CATCH("");
}

/**
 * creates a new element pointer and clones the previous element data
 * @param NewId: the ID of the new element
 * @param ThisNodes: the nodes of the new element
 * @param pProperties: the properties assigned to the new element
 * @return a Pointer to the new element
 */
Element::Pointer ElectricEnriched2Element::Clone(IndexType NewId, NodesArrayType const& ThisNodes) const
{
    KRATOS_TRY
    return Kratos::make_intrusive<ElectricEnriched2Element>(NewId, GetGeometry().Create(ThisNodes), pGetProperties());
    KRATOS_CATCH("");
}

/**
 * this determines the elemental equation ID vector for all elemental
 * DOFs
 * @param rResult: the elemental equation ID vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricEnriched2Element::EquationIdVector(EquationIdVectorType& rResult, const ProcessInfo& CurrentProcessInfo) const
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
void ElectricEnriched2Element::GetDofList(DofsVectorType& rElementalDofList, const ProcessInfo& CurrentProcessInfo) const
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
void ElectricEnriched2Element::CalculateLocalSystem(
    MatrixType& rLeftHandSideMatrix,
    VectorType& rRightHandSideVector,
    const ProcessInfo& rCurrentProcessInfo)
{
    KRATOS_TRY
    
    /// general sets  /////
    const int num_dim  = 2;
    const int num_nodes  = num_dim + 1;
    const int num_dof = num_nodes; 
    const int num_enrch = num_nodes+1; 

    GeometryData::ShapeFunctionsGradientsType DN_DX;
    Matrix N;
    Vector DetJ;
    Vector weights;
    
    //////////////////////////////////////////////////////////////////////////
    ///// matrix /////////////////////////////
    Vector distance(num_dof);
    Vector values(num_dof);
    BoundedMatrix<double,num_dof,num_dof> lhs = ZeroMatrix(num_dof,num_dof);
    BoundedMatrix<double,num_dof,num_dof> lhspos = ZeroMatrix(num_dof,num_dof);
    BoundedMatrix<double,num_dof,num_dof> lhsneg = ZeroMatrix(num_dof,num_dof);
    Vector rhs = ZeroVector(num_dof);
    array_1d<double,3> surface_sources;
    /////////////////////////////////////////////////////////////////////////////////
    const double permittivity_pos = GetProperties()[PERMITTIVITYPOS];
    const double permittivity_neg = GetProperties()[PERMITTIVITYNEG];
    surface_sources[0] = (this)->GetValue(SCHARGE);
	surface_sources[1] = (this)->GetValue(SCHARGE);
	surface_sources[2] = (this)->GetValue(SCHARGE);
    ///////////////////////////////////////////////////////////////////////////////////
    const GeometryData::IntegrationMethod integration_method = GeometryData::GI_GAUSS_1;
    GeometryType::Pointer p_geometry = this->pGetGeometry();
    const unsigned int number_of_gauss_points = p_geometry->IntegrationPointsNumber(integration_method);
    p_geometry->ShapeFunctionsIntegrationPointsGradients(DN_DX, DetJ, integration_method);
        //READ ALSO THE SHAPE FUNCTIONS THEMSELVES

    if(rLeftHandSideMatrix.size1() != num_dof)
        rLeftHandSideMatrix.resize(num_dof,num_dof,false); //resizing the system in case it does not have the right size 

    if(rRightHandSideVector.size() != num_dof)
        rRightHandSideVector.resize(num_dof,false);
    
    unsigned int nneg=0, npos=0;
    std::vector<unsigned int> pos_indices;
    std::vector<unsigned int> neg_indices;

    for (unsigned int i_node=0; i_node < num_nodes; ++i_node){

        const double dist = (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE);
        distance[i_node] = dist;
        values[i_node] = (*p_geometry)[i_node].FastGetSolutionStepValue(EPOTENCIAL);
        
        if (dist == 0.0){
            (*p_geometry)[i_node].FastGetSolutionStepValue(DISTANCE) = 1.0e-12;
        }

        if (dist > 0.0){
             npos += 1;
             pos_indices.push_back(i_node);
        }
        else{
            nneg += 1;
            neg_indices.push_back(i_node);
        }
    }
    
    
    if (nneg == num_nodes || npos == num_nodes){
          
        const GeometryType::IntegrationPointsArrayType& IntegrationPoints = p_geometry->IntegrationPoints(integration_method);
        if (weights.size() != number_of_gauss_points) {
            weights.resize(number_of_gauss_points,false);
        }
        //calculating the area of the elementes
        for (unsigned int gp = 0; gp < number_of_gauss_points; gp++){
            weights[gp] = DetJ[gp] * IntegrationPoints[gp].Weight();
           
           
            if (npos == 3){ //same results as enriched_n_element code
                lhs = -permittivity_pos*weights[gp]*prod(DN_DX[gp],trans(DN_DX[gp])); 
            }else {
                lhs = -permittivity_neg*weights[gp]*prod(DN_DX[gp],trans(DN_DX[gp]));
            }
            rhs += surface_sources*weights[gp]/3.0; 
        }
    //KRATOS_WATCH("FINISHED COMPUTING NON-CUT ELEMENTS...........")
    } else { 
            
            ///////////////////////// SPLITTED SHAPE FUNCTIONS  //////////////////////////////////////  
            ModifiedShapeFunctions::Pointer p_modified_sh_func =
                Kratos::make_shared<Triangle2D3ModifiedShapeFunctions>(p_geometry, distance);

            Matrix neg_N, pos_N;
            GeometryType::ShapeFunctionsGradientsType neg_DN_DX, pos_DN_DX;
            Vector neg_weights, pos_weights;

            Matrix int_N;
            GeometryType::ShapeFunctionsGradientsType int_DN_DX;
            Vector int_weights;

            ///////////////////////// SPLITTED SHAPE FUNCTIONS  //////////////////////////////////////  
            p_modified_sh_func->ComputePositiveSideShapeFunctionsAndGradientsValues(
                pos_N,        // N
                pos_DN_DX,    // DN_DX
                pos_weights,  // weight * detJ
                integration_method);

            p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
                neg_N,        // N
                neg_DN_DX,    // DN_DX
                neg_weights,  // weight * detJ
                integration_method);
                //KRATOS_WATCH(neg_weights)       
            /////////////////////////////////////////////////////////////////////////////////////////////////
            p_modified_sh_func->ComputeInterfaceNegativeSideShapeFunctionsAndGradientsValues(
                    int_N,
                    int_DN_DX,
                    int_weights,
                    integration_method);  
                    
            //////////////////////////////////////////////////////////////////////////////////////////////// 
            ////// ENRICHMENT FUNCTIONS//////
            const double unsigned_distance0=fabs(distance(0));
		    const double unsigned_distance1=fabs(distance(1));
		    const double unsigned_distance2=fabs(distance(2));
		    //finding the largest distance:

		    double longest_distance=fabs(unsigned_distance0);
		    if (unsigned_distance1>longest_distance)
			    longest_distance=unsigned_distance1;
		    if (unsigned_distance2>longest_distance)
		    	longest_distance=unsigned_distance2;

            const double tolerable_distance =longest_distance*0.00001;
                      
            if (unsigned_distance0<tolerable_distance)
	    		distance[0]=tolerable_distance*(distance[0]/fabs(distance[0]));
    	    if (unsigned_distance1<tolerable_distance)        
			    distance[1]=tolerable_distance*(distance[1]/fabs(distance[1]));
		    if (unsigned_distance2<tolerable_distance)
			    distance[2]=tolerable_distance*(distance[2]/fabs(distance[2]));
                      
            std::vector< Matrix > enrich_gradient(3);
            for (unsigned int i = 0; i < 3; i++)
                enrich_gradient[i].resize(2, 2, false);
            enrich_gradient[0]=ZeroMatrix(2,2);
            enrich_gradient[1]=ZeroMatrix(2,2);
            enrich_gradient[2]=ZeroMatrix(2,2); 

            ////////////////////////////////////
            double K1 = DN_DX[0](0,0)*fabs(distance[0]) + DN_DX[0](1,0)*fabs(distance[1]) + DN_DX[0](2,0)*fabs(distance[2]) - fabs(DN_DX[0](0,0)*(distance[0]) + DN_DX[0](1,0)*(distance[1]) + DN_DX[0](2,0)*(distance[2]));
            double K2 = DN_DX[0](1,0)*fabs(distance[0]) + DN_DX[0](1,1)*fabs(distance[1]) + DN_DX[0](2,1)*fabs(distance[2]) - fabs(DN_DX[0](1,0)*(distance[0]) + DN_DX[0](1,1)*(distance[1]) + DN_DX[0](2,1)*(distance[2]));
            
            KRATOS_WATCH(K1)
            KRATOS_WATCH(K2)

            double adm_Nenr_i = 0;
                           
                if (npos ==1){
                    adm_Nenr_i = 1 / (int_N(pos_indices[0]));
                } else {
                    adm_Nenr_i = 1 / (int_N(neg_indices[0]));
                }  
                   
            double adm_Nenr_j = adm_Nenr_i*K1/(1-K1);
            double adm_Nenr_k = adm_Nenr_i*K2/(1-K2);
            //for the jump, we will create a shape function that holds a constant difference of 2 along the interfase: 
            double adm_Nenr_j_b = (2.0 -K1*adm_Nenr_i)/(1-K1);
            double adm_Nenr_k_b = (2.0 -K2*adm_Nenr_i)/(1-K2);
            /*
            enrich_gradient[0](0,0) = adm_Nenr_j*DN_DX[0](order_partitions[1],0) + adm_Nenr_k*DN_DX[0](order_partitions[2],0);
            enrich_gradient[0](0,1) = adm_Nenr_j*DN_DX[0](order_partitions[1],1) + adm_Nenr_k*DN_DX[0](order_partitions[2],1); 
               
            enrich_gradient[0](1,0) = adm_Nenr_j_b*DN_DX[0](order_partitions[1],0) + adm_Nenr_k_b*DN_DX[0](order_partitions[2],0);
            enrich_gradient[0](1,1) = adm_Nenr_j_b*DN_DX[0](order_partitions[1],1) + adm_Nenr_k_b*DN_DX[0](order_partitions[2],1);

            enrich_gradient[1](0,0) = adm_Nenr_i*DN_DX[0](order_partitions[0],0);
            enrich_gradient[1](0,1) = adm_Nenr_i*DN_DX[0](order_partitions[0],1);
               
            enrich_gradient[1](1,0) = - enrich_gradient[1](0,0);
            enrich_gradient[1](1,1) = - enrich_gradient[1](0,1);

            enrich_gradient[2](0,0) = adm_Nenr_i*DN_DX[0](order_partitions[0],0);
            enrich_gradient[2](0,1) = adm_Nenr_i*DN_DX[0](order_partitions[0],1);
                   
            enrich_gradient[2](1,0) = - enrich_gradient[2](0,0);
            enrich_gradient[2](1,1) = - enrich_gradient[2](0,1);
            */

            BoundedMatrix<double, 1, 1 > K_enrich;
			noalias(K_enrich) = ZeroMatrix(1,1);
			BoundedMatrix<double, 1, 3 > B_node_enrich;
			noalias(B_node_enrich) = ZeroMatrix(1,3);
            
            /// STANDARD CONTRIBUTION
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //// POSITIVE SIDE NORMAL
            
            const std::size_t number_splitted_pos_elements = pos_weights.size();
                for (unsigned int pos_e = 0; pos_e < number_splitted_pos_elements; pos_e++){   
                    lhs -= permittivity_pos*pos_weights(pos_e)*prod(pos_DN_DX[pos_e],trans(pos_DN_DX[pos_e]));       
                }
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //// NEGATIVE SIDE NORMAL
            const std::size_t number_splitted_neg_elements = neg_weights.size();
                for (unsigned int neg_e = 0; neg_e < number_splitted_neg_elements; neg_e++){           
                    lhs -= permittivity_neg*neg_weights(neg_e)*prod(neg_DN_DX[neg_e],trans(neg_DN_DX[neg_e]));       
                }
           
            /*if (npos==1){
                lhspos -= prod(DN_DX[0],trans(DN_DX[0]))*(permittivity_pos*pos_weights[0]+permittivity_neg*neg_weights[0]+permittivity_neg*neg_weights[1]);  
            } else {
                lhsneg -= prod(DN_DX[0],trans(DN_DX[0]))*(permittivity_pos*pos_weights[0]+permittivity_pos*pos_weights[1]+permittivity_neg*neg_weights[0]); 

            }*/

            /// ENRICH CONTRIBUTION
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //// POSITIVE SIDE NORMAL
            array_1d<double,3> areas; 
            array_1d<double,3> permittivities; 
            if (npos ==1){
                areas[0] = pos_weights[0];
                areas[1] = neg_weights[0];
                areas[2] = neg_weights[1];
                permittivities[0] = permittivity_pos;
                permittivities[1] = permittivity_neg;
                permittivities[2] = permittivity_neg;
            }else {
                areas[0] = neg_weights[0];
                areas[1] = pos_weights[0];
                areas[2] = pos_weights[1];
                permittivities[0] = permittivity_neg;
                permittivities[1] = permittivity_pos;
                permittivities[2] = permittivity_pos;
            }

            for (unsigned int i = 0; i < 3; i++){   
                    K_enrich(0,0) -= (pow(enrich_gradient[i](0,0),2)+pow(enrich_gradient[i](0,1),2))*permittivities(i)*areas(i);      
                    for (unsigned int j = 0; j < 3; j++){
                    B_node_enrich(0,j) -= DN_DX[0](j,0)*(enrich_gradient[i](0,0)*permittivities(i)*areas(i)) + DN_DX[0](j,1)*(enrich_gradient[i](0,1)*permittivities(i)*areas(i));
                    
                    }
            }
           
            const double inv_K_enrich_weighted = 1/K_enrich(0,0);
            lhs -= prod(trans(B_node_enrich),B_node_enrich)*inv_K_enrich_weighted;
            ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
    noalias(rLeftHandSideMatrix) = lhs;
    noalias(rRightHandSideVector) = rhs - prod(rLeftHandSideMatrix,values);
    KRATOS_CATCH(""); 
}

void ElectricEnriched2Element::GetValueOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, 
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

void ElectricEnriched2Element::CalculateOnIntegrationPoints(const Variable<array_1d<double,3>>& rVariable, 
std::vector<array_1d<double,3> >& Output, const ProcessInfo& rCurrentProcessInfo)
{ 

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



   /////////////////////////////////////////////////////////////////////////////////

    array_1d<double,3> efield_pos = ZeroVector(3); //dimension = number of nodes . Position of the gauss point 
    array_1d<double,3> efield_neg = ZeroVector(3); //dimension = number of nodes . Position of the gauss point 
    array_1d<double,3> efield = ZeroVector(3); //dimension = number of nodes . Position of the gauss point 
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
            
        efield = -prod(trans(DN_DX[0]),values);

    //KRATOS_WATCH("FINISHED COMPUTING NON-CUT ELEMENTS...........")
    } else { //(npos != 0 && nneg != 0)//if (npos > 0 )
    
        ///////////////////////// SPLITTED SHAPE FUNCTIONS  //////////////////////////////////////    
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

        Matrix neg_N;
        GeometryType::ShapeFunctionsGradientsType neg_DN_DX;
        Vector neg_weights;
      
        p_modified_sh_func->ComputeNegativeSideShapeFunctionsAndGradientsValues(
            neg_N,        // N
            neg_DN_DX,    // DN_DX
            neg_weights,  // weight * detJ
            integration_method);
        /////////////////////////////////////////////////////////////////////////////////////////////////
        /// POSITIVE SIDE 
            efield_pos = -prod(trans(pos_DN_DX[0]),values);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //// NEGATIVE SIDE 
            efield_neg = -prod(trans(neg_DN_DX[0]),values);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
    //KRATOS_WATCH("SHOULD NOT ENETER HERE IF NOT CUT ELEMENTS ARE PRESENT..........")
    }
    
    const unsigned int number_of_gauss_points = p_geometry->IntegrationPointsNumber(integration_method);
    for (unsigned int i = 0; i < number_of_gauss_points; i++ ) {
    
       if(rVariable == EFORCEPOS)//electric force form the negative side 
        {
                noalias(Output[i]) = efield_pos;
                Output[i][2]=0.0;
        }
        else if(rVariable == EFORCENEG)//electric force form the negative side 
        {
                noalias(Output[i]) = efield_neg;
                Output[i][2]=0.0;
        }
        else
        if(rVariable == EFIELD)//electric force form the negative side 
        {
                noalias(Output[i]) = efield + efield_pos + efield_neg;
                Output[i][2]=0.0;
        }
        
    }     

}

void ElectricEnriched2Element::Calculate(const Variable<array_1d<double, 3 > >& rVariable,
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
void ElectricEnriched2Element::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the elemental right hand side vector only
 * @param rRightHandSideVector: the elemental right hand side vector
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricEnriched2Element::CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
{
}

/**
 * this is called during the assembling process in order
 * to calculate the first derivatives contributions for the LHS and RHS
 * @param rLeftHandSideMatrix: the elemental left hand side matrix
 * @param rRightHandSideVector: the elemental right hand side
 * @param rCurrentProcessInfo: the current process info instance
 */
void ElectricEnriched2Element::CalculateFirstDerivativesContributions(
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
void ElectricEnriched2Element::CalculateFirstDerivativesLHS(MatrixType& rLeftHandSideMatrix, const ProcessInfo& rCurrentProcessInfo)
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
void ElectricEnriched2Element::CalculateFirstDerivativesRHS(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo)
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
void ElectricEnriched2Element::CalculateSecondDerivativesContributions(
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
void ElectricEnriched2Element::CalculateSecondDerivativesLHS(
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
void ElectricEnriched2Element::CalculateSecondDerivativesRHS(
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
void ElectricEnriched2Element::CalculateMassMatrix(MatrixType& rMassMatrix, const ProcessInfo& rCurrentProcessInfo)
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
void ElectricEnriched2Element::CalculateDampingMatrix(MatrixType& rDampingMatrix, const ProcessInfo& rCurrentProcessInfo)
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
int ElectricEnriched2Element::Check(const ProcessInfo& rCurrentProcessInfo) const
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

std::string ElectricEnriched2Element::Info() const {
    std::stringstream buffer;
    buffer << "ElectricEnrichedElement #" << Id();
    return buffer.str();
}

/// Print information about this object.

void ElectricEnriched2Element::PrintInfo(std::ostream& rOStream) const {
    rOStream << "ElectricEnrichedElement #" << Id();
}

/// Print object's data.

void ElectricEnriched2Element::PrintData(std::ostream& rOStream) const {
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

void ElectricEnriched2Element::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Element );

    // List
    // To be completed with the class member list
}

void ElectricEnriched2Element::load(Serializer& rSerializer)
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
inline std::istream & operator >> (std::istream& rIStream, ElectricEnriched2Element& rThis);

/// output stream function
inline std::ostream & operator << (std::ostream& rOStream, const ElectricEnriched2Element& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << " : " << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace Kratos.
