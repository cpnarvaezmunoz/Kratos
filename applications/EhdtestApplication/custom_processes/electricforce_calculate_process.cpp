// System includes

// External includes

// Project includes
#include "custom_processes/electricforce_calculate_process.h"
#include "ehdtest_application_variables.h"
#include "utilities/geometry_utilities.h"
#include "includes/node.h"
#include "includes/element.h"
#include "includes/model_part.h"
#include "utilities/variable_utils.h"
#include "ehdtest_application.h"
namespace Kratos
{
double ElectricforceCalculateProcess::CalculateElectricforce(Element& rElement, const std::size_t DomainSize)
{
    KRATOS_TRY

    // We get the element geometry
    auto& r_this_geometry = rElement.GetGeometry();
    const std::size_t local_space_dimension = r_this_geometry.LocalSpaceDimension();
    const std::size_t number_of_nodes = r_this_geometry.size();
    std::vector<array_1d<double, 3>> gem(number_of_nodes);
    
   
   KRATOS_CATCH("")
}
void ElectricforceCalculateProcess::Execute()
{
    KRATOS_TRY

    ProcessInfo& proc_info = mrThisModelPart.GetProcessInfo();
    array_1d<double,3> dummy;
   
        /*for(ModelPart::NodesContainerType::iterator in = mrThisModelPart.NodesBegin() ;
            in != mrThisModelPart.NodesEnd() ; ++in)
        {
        	in->GetValue(EFORCEPOS)= 0.0;		
        }*/
 
    for(ModelPart::ElementsContainerType::iterator im = mrThisModelPart.ElementsBegin() ;
            im != mrThisModelPart.ElementsEnd() ; ++im)
    {
      	im->Calculate(EFORCEPOS,dummy,proc_info);
        im->Calculate(EFORCENEG,dummy,proc_info);  		
    }

    KRATOS_CATCH("")
}
} // namespace Kratos
 
