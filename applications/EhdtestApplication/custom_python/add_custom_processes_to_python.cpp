// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes

// External includes


// Project includes
#include "custom_python/add_custom_processes_to_python.h"
#include "ehdtest_application_variables.h"


//Processes
#include "custom_processes/electricfield_calculate_process.h"
#include "custom_processes/electricforce_calculate_process.h"

namespace Kratos {
namespace Python {

void  AddCustomProcessesToPython(pybind11::module& m)
{
    namespace py = pybind11;

    /// Processes
       py::class_<ElectricfieldCalculateProcess, ElectricfieldCalculateProcess::Pointer, Process>(m,"ElectricfieldCalculateProcess")
        .def(py::init<ModelPart&>())
        ;
       py::class_<ElectricforceCalculateProcess, ElectricforceCalculateProcess::Pointer, Process>(m,"ElectricforceCalculateProcess")
        .def(py::init<ModelPart&>())
        ;
}

}  // namespace Python.
} // Namespace Kratos

