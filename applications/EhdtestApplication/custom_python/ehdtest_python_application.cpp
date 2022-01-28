//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


// System includes

#if defined(KRATOS_PYTHON)
// External includes
#include <pybind11/pybind11.h>


// Project includes
#include "includes/define_python.h"
#include "ehdtest_application.h"
#include "ehdtest_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"
#include "custom_python/add_custom_processes_to_python.h"

namespace Kratos {
namespace Python {

PYBIND11_MODULE(KratosEhdtestApplication,m)
{
    namespace py = pybind11;

    py::class_<KratosEhdtestApplication,
        KratosEhdtestApplication::Pointer,
        KratosApplication>(m, "KratosEhdtestApplication")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);
    AddCustomProcessesToPython(m);

    //registering variables in python
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, EPOTENCIAL )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PERMITTIVITYPOS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PERMITTIVITYNEG )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONDUCTIVITYPOS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, CONDUCTIVITYNEG )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, PCHARGE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, SCHARGE )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, VCHARGE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EFIELD )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EDISPLACEMENT )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EFORCE )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EFORCENEG )
    KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, EFORCEPOS )
    KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ICOEFFICIENT )

    //	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);

}

} // namespace Python.
} // namespace Kratos.

#endif // KRATOS_PYTHON defined
