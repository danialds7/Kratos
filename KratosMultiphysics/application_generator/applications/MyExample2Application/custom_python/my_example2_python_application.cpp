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

// External includes

// Project includes
#include "includes/define_python.h"
#include "my_example2_application.h"
#include "my_example2_application_variables.h"
#include "custom_python/add_custom_strategies_to_python.h"
#include "custom_python/add_custom_utilities_to_python.h"


namespace Kratos::Python {

PYBIND11_MODULE(KratosMyExample2Application,m)
{
    namespace py = pybind11;

    py::class_<KratosMyExample2Application,
        KratosMyExample2Application::Pointer,
        KratosApplication>(m, "KratosMyExample2Application")
        .def(py::init<>())
        ;

    AddCustomStrategiesToPython(m);
    AddCustomUtilitiesToPython(m);

    //registering variables in python
      KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DOF_1 )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, DOF_2 )
  KRATOS_REGISTER_IN_PYTHON_VARIABLE(m, ScalarVariable )
  KRATOS_REGISTER_IN_PYTHON_3D_VARIABLE_WITH_COMPONENTS(m, VectorVariable )

    //	KRATOS_REGISTER_IN_PYTHON_VARIABLE(NODAL_AREA);

}

} // namespace Kratos::Python
