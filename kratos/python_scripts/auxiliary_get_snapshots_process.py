import KratosMultiphysics
import numpy as np

def Factory(settings, model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return AuxiliaryGetSnapshotsProcess(model, settings["Parameters"])

## All the processes python should be derived from "Process"
class AuxiliaryGetSnapshotsProcess(KratosMultiphysics.Process):
    """This class is a dummy-process that shows how the functions that can be implemented
    in order to customize the behavior

    Public member variables:
    model -- the container of the different model parts.
    settings -- Kratos parameters containing process settings.
    """

    def __init__(self, model, settings):
        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing process settings.
        """
        KratosMultiphysics.Process.__init__(self) # calling the baseclass constructor

        default_settings = KratosMultiphysics.Parameters("""{
            "computing_model_part_name" : "ComputingModelPartName",
            "list_of_variables" : ["VARIABLE_1","VARIABLE_2"],
            "output_name" : "output_name",
            "output_frequency" : 0.0
        }""")

        # Add missing settings that the user did not provide but that
        # are necessary for this process
        settings.ValidateAndAssignDefaults(default_settings)
        
        self.model = model
        self.settings = settings

    def ExecuteInitialize(self):
        self.computing_model_part_name = self.settings["computing_model_part_name"].GetString()
        self.computing_model_part = self.model.GetModelPart(self.computing_model_part_name)
        self.list_of_variables = self.settings["list_of_variables"].GetStringArray()
        self.output_name = self.settings["output_name"].GetString()
        self.results = [[] for i in range(len(self.list_of_variables))]
        self.output_frequency = self.settings["output_frequency"].GetDouble()
        self.output_frequency_counter = 0
        self.output_frequency_counter_global = 0
        

    def Check(self):
        """ This method verifies that the input is correct

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteBeforeSolutionLoop(self):
        """ This method is executed just before the solution-loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteInitializeSolutionStep(self):
        """ This method is executed in order to initialize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteBeforeOutputStep(self):
        """ This method is executed before writing the output (if output
        is being written in this step)

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteAfterOutputStep(self):
        """ This method is executed after writing the output (if output
        is being written in this step)

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        pass

    def ExecuteFinalizeSolutionStep(self):
        """ This method is executed in order to finalize the current step

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        self.output_frequency_counter += 1
        self.output_frequency_counter_global += 1
        var_counter = 0
        for var_name in self.list_of_variables:
            var = KratosMultiphysics.KratosGlobals.GetVariable(var_name)
            current_snapshot = []
            for node in self.computing_model_part.Nodes:
                current_snapshot.append(node.GetSolutionStepValue(var))
            self.results[var_counter].append(np.array(current_snapshot))
            var_counter += 1
        if self.output_frequency_counter>=self.output_frequency:
            var_counter = 0
            for var_name in self.list_of_variables:
                file_name = self.output_name + "_" + var_name + "_" + str(self.output_frequency_counter_global) + ".npy"
                current_numpy = np.array(self.results[var_counter]).T
                np.save(file_name, current_numpy)
                var_counter += 1
            self.results = [[] for i in range(len(self.list_of_variables))]
            self.output_frequency_counter = 0
        debug = True

    def ExecuteFinalize(self):
        """ This method is executed after the computations, at the end of the solution-loop

        Keyword arguments:
        self -- It signifies an instance of a class.
        """
        var_counter = 0
        for var_name in self.list_of_variables:
            file_name = self.output_name + "_" + var_name + ".npy"
            current_numpy = np.array(self.results[var_counter]).T
            np.save(file_name, current_numpy)
            var_counter += 1