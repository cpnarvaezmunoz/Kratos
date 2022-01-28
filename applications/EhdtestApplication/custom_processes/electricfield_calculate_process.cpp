// System includes

// External includes

// Project includes
#include "custom_processes/electricfield_calculate_process.h"
#include "ehdtest_application_variables.h"
#include "utilities/geometry_utilities.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
namespace Kratos
{
double ElectricfieldCalculateProcess::CalculateElectricfield(Element& rElement, const std::size_t DomainSize)
{
    KRATOS_TRY

    // We get the element geometry
    auto& r_this_geometry = rElement.GetGeometry();
    const std::size_t local_space_dimension = r_this_geometry.LocalSpaceDimension();
    const std::size_t number_of_nodes = r_this_geometry.size();
    std::vector<array_1d<double, 3>> gem(number_of_nodes);
    
   
   KRATOS_CATCH("")
}
void ElectricfieldCalculateProcess::Execute()
{
KRATOS_TRY

        ProcessInfo& proc_info = mrThisModelPart.GetProcessInfo();
        bounded_matrix<double,3,2> DN_DX= ZeroMatrix(3,2);  // Gradients matrix 
        array_1d<double,3> N; //dimension = number of nodes . Position of the gauss point 
        array_1d<double,3> efieldx; //dimension = number of nodes . Position of the gauss point 
        array_1d<double,3> efieldy; //dimension = number of nodes . Position of the gauss point 
        array_1d<double,3> temp; //dimension = number of nodes . . since we are using a residualbased approach 
        auto& nodes_array = mrThisModelPart.GetCommunicator().LocalMesh().Nodes();

        const int number_of_nodes = static_cast<int>(nodes_array.size());        
        double nodal_area=0;
        double area;
        for(ModelPart::ElementsContainerType::iterator im = mrThisModelPart.ElementsBegin() ;
                 im != mrThisModelPart.ElementsEnd() ; ++im)
        {
        				//get the geometry
		Geometry< Node<3> >& geometry = im->GetGeometry();
        	GeometryUtils::CalculateGeometryData(geometry, DN_DX, N, area);

        for(unsigned int iii = 0; iii<number_of_nodes; iii++)
        {
		temp[iii] = geometry[iii].FastGetSolutionStepValue(EPOTENCIAL);
        //geometry[iii].FastGetSolutionStepValue(EFIELD) += -prod(trans(DN_DX),temp); //with this line i got values even along the Z axis 
        }
       
        KRATOS_WATCH(DN_DX)
        //KRATOS_WATCH(temp)
 	 	//// I used these lines 
        geometry[0].FastGetSolutionStepValue(EFIELD_X) = DN_DX(0,0)*temp[0] + DN_DX(1,0)*temp[1]  + DN_DX(2,0)*temp[2];
		geometry[0].FastGetSolutionStepValue(EFIELD_Y) = DN_DX(0,1)*temp[0] + DN_DX(1,1)*temp[1]  + DN_DX(2,1)*temp[2];

        geometry[1].FastGetSolutionStepValue(EFIELD_X) = DN_DX(0,0)*temp[0] + DN_DX(1,0)*temp[1]  + DN_DX(2,0)*temp[2];
		geometry[1].FastGetSolutionStepValue(EFIELD_Y) = DN_DX(0,1)*temp[0] + DN_DX(1,1)*temp[1]  + DN_DX(2,1)*temp[2];

		geometry[2].FastGetSolutionStepValue(EFIELD_X) = DN_DX(0,0)*temp[0] + DN_DX(1,0)*temp[1]  + DN_DX(2,0)*temp[2];
		geometry[2].FastGetSolutionStepValue(EFIELD_Y) = DN_DX(0,1)*temp[0] + DN_DX(1,1)*temp[1]  + DN_DX(2,1)*temp[2];


        //efield[3] = 0.0;
        
        //geometry[0].FastGetSolutionStepValue(EFIELD) = efield;
              
        }

        KRATOS_CATCH("")
}
} // namespace Kratos
 