# Application dependent names and paths
from KratosMultiphysics import _ImportApplication, python_registry_utilities
from KratosMyExample2Application import *

application = KratosMyExample2Application()
application_name = "KratosMyExample2Application"

_ImportApplication(application, application_name)

from . import python_registry_lists
python_registry_utilities.RegisterModelersList("KratosMultiphysics.MyExample2Application", python_registry_lists)
python_registry_utilities.RegisterOperationsList("KratosMultiphysics.MyExample2Application", python_registry_lists)
python_registry_utilities.RegisterProcessesList("KratosMultiphysics.MyExample2Application", python_registry_lists)