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
#include "my_example2_application.h"
#include "my_example2_application_variables.h"


namespace Kratos {

KratosMyExample2Application::KratosMyExample2Application():
    KratosApplication("MyExample2Application")
    {}

void KratosMyExample2Application::Register()
{
     KRATOS_INFO("") << "Initializing KratosMyExample2Application..." << std::endl;

      KRATOS_REGISTER_VARIABLE( DOF_1 )
  KRATOS_REGISTER_VARIABLE( DOF_2 )
  KRATOS_REGISTER_VARIABLE( ScalarVariable )
  KRATOS_REGISTER_3D_VARIABLE_WITH_COMPONENTS( VectorVariable )

}

}  // namespace Kratos.
