import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.python_solver import PythonSolver
from KratosMultiphysics.OptimizationApplication.optimization_info import OptimizationInfo
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.controls.control_wrapper import ControlWrapper
from KratosMultiphysics.OptimizationApplication.responses.response_function_wrapper import ResponseFunctionBaseWrapper
from KratosMultiphysics.OptimizationApplication.responses.response_function_wrapper import ObjectiveResponseFunctionWrapper
from KratosMultiphysics.OptimizationApplication.responses.response_function_wrapper import ConstraintResponseFunctionWrapper
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import GetSensitivityContainer
from KratosMultiphysics.OptimizationApplication.utilities.helper_utils import ContainerEnum

class GradientProjectionSolver(PythonSolver):
    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "max_correction_share": 0.75,
            "relative_tolerance"  : 1e-3,
            "step_size"           : 1e-3,
            "control_names_list"  : []
        }""")

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters, optimization_info: OptimizationInfo):
        super().__init__(model, settings)
        self.optimization_info = optimization_info
        if self.optimization_info.GetBufferSize() < 2:
            self.optimization_info.SetBufferSize(2)

        self.max_correction_share = settings["max_correction_share"].GetDouble()
        self.relative_tolerance = settings["relative_tolerance"].GetDouble()
        self.step_size = settings["step_size"].GetDouble()

        # objectives and constraints
        self.objectives_list = []
        self.constraints_list = []
        for response in optimization_info.GetRoutines(ResponseFunctionBaseWrapper):
            if isinstance(response, ObjectiveResponseFunctionWrapper):
                self.objectives_list.append(response)
            elif isinstance(response, ConstraintResponseFunctionWrapper):
                self.constraints_list.append(response)
            else:
                raise RuntimeError("Unsupproted response type found.")

        if len(self.objectives_list) == 0:
            raise RuntimeError("Atleast one objective should be defined for optimization.")

        if len(self.objectives_list) > 1:
            raise RuntimeError("Only one objective is allowed in the optimization.")

        self.objective_response_function_wrapper: ResponseFunctionBaseWrapper = self.objectives_list[0]
        self.control_wrappers_list = [control_wrapper for control_wrapper in self.optimization_info.GetRoutines(ControlWrapper) if control_wrapper.GetName() in self.settings["control_names_list"].GetStringArray()]

    def SolveSolutionStep(self) -> bool:
        # calculate objective value
        self.optimization_info["objective"]["value"] = self.objective_response_function_wrapper.GetValue()
        self.optimization_info["objective"]["standardized_value"] = self.objective_response_function_wrapper.GetStandardizedValue()

        # calculate constraint values
        constraint_values = []
        for constraint_wrapper in self.constraints_list:
            constraint_wrapper: ResponseFunctionBaseWrapper = constraint_wrapper
            constraint_data = {
                "value": constraint_wrapper.GetValue(),
                "standardized_value": constraint_wrapper.GetStandardizedValue(),
                "is_active": constraint_wrapper.IsActive()
            }
            constraint_value = constraint_wrapper.GetStandardizedValue()
            constraint_values.append(constraint_data)
        self.optimization_info["constraints"] = constraint_values

        # calculate objective gradients
        controls_data = []
        for control_wrapper in self.control_wrappers_list:
            control_wrapper: ControlWrapper = control_wrapper
            control: Control = ControlWrapper.GetControl()

            control_model_part = control.GetModelPart()
            control_container_type = control.GetContainerType()
            control_sensitivty_variable = control.GetControlSensitivityVariable()

            modified_objective_sensitivities = GradientProjectionSolver.__GetModifiedSensitivities(
                control_wrapper, self.objective_response_function_wrapper, control_sensitivty_variable, control_model_part, control_container_type)

            self.optimization_info["objective"]["modified_sensitivity"] = modified_objective_sensitivities

            active_constraint_values = []
            modified_active_constraints_sensitivities = []
            for i, constraint_wrapper in enumerate(self.constraints_list):
                constraint_wrapper: ResponseFunctionBaseWrapper = constraint_wrapper

                if self.optimization_info["constraints"][i]["is_active"]:
                    active_constraint_values.append(constraint_value)
                    modified_active_constraint_sensitivities = GradientProjectionSolver.__GetModifiedSensitivities(
                        control_wrapper, constraint_wrapper, control_sensitivty_variable, control_model_part, control_container_type)
                    modified_active_constraints_sensitivities.append(modified_active_constraint_sensitivities)
                    self.optimization_info["constraints"][i]["modified_sensitivity"] = modified_active_constraint_sensitivities

            active_constraint_values_vector = Kratos.Vector(len(active_constraint_values))
            for i, v in enumerate(active_constraint_values):
                active_constraint_values_vector[i] = v

            control_container = GetSensitivityContainer(control_model_part)
            control_domain_size = control_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]

            search_direction_variable, search_correction_variable = GradientProjectionSolver.__GetSearchVariables(control_sensitivty_variable)

            KratosOA.GradientProjectionSolverUtils.CalculateProjectedSearchDirectionAndCorrection(
                control_container,
                control_domain_size,
                search_direction_variable,
                search_correction_variable,
                active_constraint_values_vector,
                modified_objective_sensitivities,
                modified_active_constraints_sensitivities)

            KratosOA.GradientProjectionSolverUtils.CalculateControlChange(
                control_container,
                control_model_part.GetCommunicator().GetDataCommunicator(),
                search_direction_variable,
                search_correction_variable,
                control.GetControlUpdateVariable(),
                self.step_size,
                self.max_correction_share)

            control_update_vector = Kratos.Vector()
            KratosOA.OptimizationUtils.GetContainerVariableToVector(
                control_container,
                control.GetControlUpdateVariable(),
                control_domain_size,
                control_update_vector)

            control.SetControlUpdatesVector(control_wrapper.ModifyControlUpdates(control_update_vector))
            controls_data.append({
                "modified_control_update": control.GetControlUpdatesVector()
            })



        self.optimization_info["controls"] = controls_data
        return self._CheckConvergence()

    def _CheckConvergence(self) -> bool:
        if self.optimization_info["step"] > 1:
            # check for objective convergence
            is_converged = abs(self.optimization_info["objective"]["value"] / self.optimization_info["objective"]["value", 1] - 1) < self.relative_tolerance

            # check for constraint convergence
            for constraint_data in self.optimization_info["constraints"]:
                is_converged = is_converged and not constraint_data["is_active"]

            return is_converged
        else:
            return False

    @staticmethod
    def __GetModifiedSensitivities(control_wrapper: ControlWrapper, response_function_wrapper: ResponseFunctionBaseWrapper, sensitivity_variable, sensitivity_model_part: Kratos.ModelPart, sensitivity_container_type: ContainerEnum) -> Kratos.Vector:
        sensitivities = response_function_wrapper.GetStandardizedSensitivity(
                                        sensitivity_variable, sensitivity_model_part, sensitivity_container_type)
        return control_wrapper.ModifySensitivities(sensitivities)

    @staticmethod
    def __GetSearchVariables(sensitivty_variable):
        sensitivity_variable_type = Kratos.KratosGlobals.GetVariableType(sensitivty_variable.Name())
        match sensitivity_variable_type:
            case "Double":
                search_direction_variable = KratosOA.SCALAR_SEARCH_DIRECTION
                search_correction_variable = KratosOA.SCALAR_SEARCH_CORRECTION
            case "Array":
                search_direction_variable = KratosOA.VECTOR_SEARCH_DIRECTION
                search_correction_variable = KratosOA.VECTOR_SEARCH_CORRECTION
            case _:
                raise RuntimeError(f"Unsupported {sensitivty_variable.Name()} of type {sensitivity_variable_type} is used. This only supports double and array variable types.")

        return search_direction_variable, search_correction_variable
