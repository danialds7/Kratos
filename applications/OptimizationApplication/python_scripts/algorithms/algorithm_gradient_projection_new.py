import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.opt_convergence import CreateConvergenceCriteria
from KratosMultiphysics.OptimizationApplication.utilities.opt_line_search import CreateLineSearch
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_constraint import StandardizedConstraint
from KratosMultiphysics.LinearSolversApplication.dense_linear_solver_factory import ConstructSolver
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmSteepestDescent(model, parameters, optimization_problem)

class AlgorithmSteepestDescent(Algorithm):
    """
        A classical steepest descent algorithm to solve unconstrainted optimization problems.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "model_part_name"   : "OptimizationModelPart",
            "objective"         : {},
            "constraints"       : [],
            "controls"          : [],
            "echo_level"        : 0,
            "settings"          : {
                "echo_level"      : 0,
                "gradient_scaling": "inf_norm",
                "line_search"     : {},
                "conv_settings"   : {}
            }
        }""")
    
    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        self.model = model
        self.parameters = parameters
        self._optimization_problem = optimization_problem
        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.master_control = MasterControl() # Need to fill it with controls

        for control_name in parameters["controls"].GetStringArray():
            control = optimization_problem.GetControl(control_name)
            self.master_control.AddControl(control)


        settings = parameters["settings"]
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["settings"])

        self.gradient_scaling = settings["gradient_scaling"].GetString()
        self.echo_level = settings["echo_level"].GetInt()

        ComponentDataView("algorithm", self._optimization_problem).SetDataBuffer(self.GetMinimumBufferSize())

        self.__convergence_criteria = CreateConvergenceCriteria(settings["conv_settings"], self._optimization_problem)
        self.__line_search_method = CreateLineSearch(settings["line_search"], self._optimization_problem)

        self.__objective = StandardizedObjective(parameters["objective"], self.master_control, self._optimization_problem)
        for constraint_param in parameters["constraints"]:
            constraint = StandardizedConstraint(constraint_param, self.master_control, optimization_problem)
            self.__constraints_list.append(constraint)
        self.__control_field = None
        self.__obj_val = None

        default_linear_solver_settings = Kratos.Parameters("""{
            "solver_type": "LinearSolversApplication.dense_col_piv_householder_qr"
        }""")
        settings["linear_solver_settings"].ValidateAndAssignDefaults(default_linear_solver_settings)
        self.linear_solver = ConstructSolver(settings["linear_solver_settings"])

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self):
        pass
        
    def Initialize(self):
        self.converged = False
        self.__obj_val = None
        self.__objective.GetReponse().Initialize()
        self.__objective.Initialize()
        CallOnAll(self.optimization_problem.__constraints_list(), StandardizedConstraint.Initialize())
        self.__objective.Check()
        CallOnAll(self.optimization_problem.__constraints_list(), StandardizedConstraint.Check())
        self.__control_field = self.master_control.GetControlField() # GetInitialControlFields() later

    def Finalize(self):
        pass
    
    def ComputeSearchDirection(self, obj_grad) -> KratosOA.ContainerExpression.CollectiveExpressions:
        return obj_grad * -1.0
    
    def GetCurrentObjValue(self) -> float:
        return self.__obj_val
    
    def GetCurrentControlField(self):
        return self.__control_field

    def SolveOptimizationProblem(self) -> bool:
        self.Initialize()
        with TimeLogger("Solve Optimization problem", "Start", "End"):
            self.Solve()
        return self.converged
    
    def Solve(self):
        algorithm_data = ComponentDataView("algorithm", self._optimization_problem)
        while not self.converged:
            print("")
            with TimeLogger("Optimization", f" Start Iteration {self._optimization_problem.GetStep()}", f"End Iteration {self._optimization_problem.GetStep()}"):

                with TimeLogger("Calculate objective value", "Start", "End"):
                    self.__obj_val = self.__objective.CalculateStandardizedValue(self.__control_field) 
                    algorithm_data.GetBufferedData()["std_obj_value"] = self.__obj_val
                    print(self.__objective.GetInfo())

                with TimeLogger("Calculate gradient", "Start", "End"):
                    obj_grad = self.__objective.CalculateStandardizedGradient()
                
                with TimeLogger("Calculate design update", "Start", "End"):
                    search_direction = self.ComputeSearchDirection(obj_grad)
                    algorithm_data.GetBufferedData()["search_direction"] = search_direction

                    alpha = self.__line_search_method.ComputeStep()

                    update = search_direction * alpha
                    self.__control_field += update

                algorithm_data.GetBufferedData()["parameter_update"] = update
                algorithm_data.GetBufferedData()["control_field"] = self.__control_field

                self.converged = self.__convergence_criteria.IsConverged()

                self.CallOnAllProcesses(["output_processes"], Kratos.OutputProcess.PrintOutput)

                self._optimization_problem.AdvanceStep()

            self.Finalize()
    
    def GetOptimizedObjectiveValue(self) -> float:
        if self.converged:
            return self.__obj_val
        else:
            raise RuntimeError("Optimization problem hasn't been solved.")