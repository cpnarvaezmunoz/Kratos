# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosEhdtestApplication import *
application = KratosEhdtestApplication()
application_name = "KratosEhdtestApplication"

_ImportApplication(application, application_name)
