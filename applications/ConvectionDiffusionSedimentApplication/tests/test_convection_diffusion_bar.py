from KratosMultiphysics import *
import KratosMultiphysics.ConvectionDiffusionSedimentApplication
from KratosMultiphysics.ConvectionDiffusionSedimentApplication.convection_diffusion_sediment_analysis import ConvectionDiffusionSedimentAnalysis
import KratosMultiphysics.KratosUnittest as KratosUnittest

from math import *

class TestConvectionDiffusionSedimentBar(KratosUnittest.TestCase):

    def _SetUpBarSimulation(self, ProjectParametersFileName):
        with open(ProjectParametersFileName,'r') as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())
        model = KratosMultiphysics.Model()
        bar_simulation = ConvectionDiffusionSedimentAnalysis(model, parameters)

        return bar_simulation

    def _CalculateErrorNorm(self, BarSimulation):
        # Check L2 error via midpoint rule
        error = 0
        KratosMultiphysics.CalculateNonHistoricalNodalAreaProcess(BarSimulation._GetSolver().main_model_part,2).Execute()
        for node in BarSimulation._GetSolver().main_model_part.Nodes:
            # L2 norm
            x = node.X
            u_analytical = x - sin(BarSimulation.time)
            u_numerical = node.GetSolutionStepValue(KratosMultiphysics.TEMPERATURE)
            nodal_area = node.GetValue(KratosMultiphysics.NODAL_AREA)
            error += nodal_area * ((u_analytical - u_numerical)**2)
        error = sqrt(error)

        return error

    def testConvectionDiffusionSedimentBarExplicitElementUnsteadyDOSS(self):
        project_parameters_file_name = "test_convection_diffusion_sediment_bar/project_parameters_bar_DOSS.json"
        bar_simulation = self._SetUpBarSimulation(project_parameters_file_name)
        bar_simulation.Run()

        error = self._CalculateErrorNorm(bar_simulation)
        self.assertAlmostEqual(error,0.0017678016487118946,delta=1e-12)

    def testConvectionDiffusionSedimentBarSemiImplicit(self):
        project_parameters_file_name = "test_convection_diffusion_sediment_bar/project_parameters_bar_semi_implicit.json"
        bar_simulation = self._SetUpBarSimulation(project_parameters_file_name)
        bar_simulation.Run()

        error = self._CalculateErrorNorm(bar_simulation)
        self.assertAlmostEqual(error, 0.05151520448578559, delta=1e-12)

    def testConvectionDiffusionSedimentBarExplicitElementUnsteadyQOSS(self):
        project_parameters_file_name = "test_convection_diffusion_sediment_bar/project_parameters_bar_QOSS.json"
        bar_simulation = self._SetUpBarSimulation(project_parameters_file_name)
        bar_simulation.Run()

        error = self._CalculateErrorNorm(bar_simulation)
        self.assertAlmostEqual(error,0.0017678016487242358,delta=1e-12)

    def testConvectionDiffusionSedimentBarExplicitElementUnsteadyDASGS(self):
        project_parameters_file_name = "test_convection_diffusion_sediment_bar/project_parameters_bar_DASGS.json"
        bar_simulation = self._SetUpBarSimulation(project_parameters_file_name)
        bar_simulation.Run()

        error = self._CalculateErrorNorm(bar_simulation)
        self.assertAlmostEqual(error,0.04953302968545053,delta=1e-12)

    def testConvectionDiffusionSedimentBarExplicitElementUnsteadyQASGS(self):
        project_parameters_file_name = "test_convection_diffusion_sediment_bar/project_parameters_bar_QASGS.json"
        bar_simulation = self._SetUpBarSimulation(project_parameters_file_name)
        bar_simulation.Run()

        error = self._CalculateErrorNorm(bar_simulation)
        self.assertAlmostEqual(error,0.045232415452161175,delta=1e-12)

if __name__ == '__main__':
    KratosUnittest.main()

