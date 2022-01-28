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

#include "ehdtest_application_variables.h"

namespace Kratos
{
KRATOS_CREATE_VARIABLE( double, EPOTENCIAL )
KRATOS_CREATE_VARIABLE( double, PERMITTIVITYPOS )
KRATOS_CREATE_VARIABLE( double, PERMITTIVITYNEG )
KRATOS_CREATE_VARIABLE( double, CONDUCTIVITYPOS )
KRATOS_CREATE_VARIABLE( double, CONDUCTIVITYNEG )
KRATOS_CREATE_VARIABLE( double, PCHARGE )
KRATOS_CREATE_VARIABLE( double, SCHARGE )
KRATOS_CREATE_VARIABLE( double, VCHARGE )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EFIELD )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EDISPLACEMENT )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EFORCE )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EFORCEPOS )
KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( EFORCENEG )
KRATOS_CREATE_VARIABLE( double, ICOEFFICIENT )

}
