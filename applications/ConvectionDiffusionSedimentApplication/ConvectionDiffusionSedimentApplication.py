
# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosConvectionDiffusionSedimentApplication import *
application = KratosConvectionDiffusionSedimentApplication()
application_name = "KratosConvectionDiffusionSedimentApplication"

_ImportApplication(application, application_name)
