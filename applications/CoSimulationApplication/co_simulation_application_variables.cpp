#include "co_simulation_application_variables.h"

namespace Kratos
{

KRATOS_CREATE_VARIABLE(double, SCALAR_DISPLACEMENT);
KRATOS_CREATE_VARIABLE(double, SCALAR_ROOT_POINT_DISPLACEMENT);
KRATOS_CREATE_VARIABLE(double, SCALAR_REACTION);
KRATOS_CREATE_VARIABLE(double, SCALAR_FORCE);
KRATOS_CREATE_VARIABLE(double, SCALAR_VOLUME_ACCELERATION);

KRATOS_CREATE_VARIABLE(int, COUPLING_ITERATION_NUMBER);

KRATOS_CREATE_VARIABLE(int, INTERFACE_EQUATION_ID);
KRATOS_CREATE_VARIABLE(int, EXPLICIT_EQUATION_ID);

KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS(MIDDLE_VELOCITY)
}