# Importing the Kratos Library
import KratosMultiphysics as KM

# Importing the MPI
if KM.IsDistributedRun():
    import KratosMultiphysics.mpi as KratosMPI

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as KratosUnittest

# Importing the Kratos utilities
import KratosMultiphysics.kratos_utilities as kratos_utils

# Importing the read model part utility from the testing utilities
from KratosMultiphysics.testing.utilities import ReadModelPart

# Importing the OS module
import os

def GetFilePath(fileName):
    """
    Get the absolute path of the given file name.
    """
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def RemoveFiles(mdpa_name):
    """
    Remove the time-related files associated with a given .mdpa file.
    """
    kratos_utils.DeleteFileIfExisting(mdpa_name + ".time")

class TestSearchWrapper(KratosUnittest.TestCase):
    """
    A test class for the SearchWrapperGeometricalObjectBins class.

    This class is used to test the functionality of the GeometricalObjectBins search wrapper
    for searching and querying geometric objects within a specified radius.

    Attributes:
        current_model (KM.Model): The Kratos Model used for testing.
        model_part (KM.ModelPart): The main model part for testing.
        sub_model_part (KM.ModelPart): The sub-model part for testing.
        mdpa_name (str): The name of the MDPA file used for testing.
        data_comm (KM.DataCommunicator): The data communicator for parallel execution.
        node_id (int): The ID of the test node.
        node (KM.Node): The test node used for searching.
        search (KM.SearchWrapperGeometricalObjectBins): The search wrapper instance for testing.

    Methods:
        setUpClass(cls):
            Set up the model part and its related entities for the test.

        tearDownClass(cls):
            Clean up after all tests are run.

        setUp(self):
            Set up for each individual test.

        _create_search(self, search_type):
            Create an appropriate search wrapper based on the given search type.

        test_SearchWrapperGeometricalObjectBins_SearchInRadius(self):
            Test for the 'SearchInRadius' method of the GeometricalObjectBins search wrapper.

        test_SearchWrapperGeometricalObjectBins_SearchNearestInRadius(self):
            Test for the 'SearchNearestInRadius' method of the GeometricalObjectBins search wrapper.

        test_SearchWrapperGeometricalObjectBins_SearchNearest(self):
            Test for the 'SearchNearest' method of the GeometricalObjectBins search wrapper.

        test_SearchWrapperGeometricalObjectBins_SearchIsInside(self):
            Test for the 'SearchIsInside' method of the GeometricalObjectBins search wrapper.
    """
    @classmethod
    def setUpClass(cls):
        """
        Setting up the model part and its related entities for the test.
        """
        # Create model and model parts
        cls.current_model = KM.Model()
        cls.model_part = cls.current_model.CreateModelPart("Main")
        cls.sub_model_part = cls.model_part.CreateSubModelPart("SubModelPart")

        # Define model part variables
        cls.model_part.ProcessInfo[KM.DOMAIN_SIZE] = 3
        cls.model_part.AddNodalSolutionStepVariable(KM.BULK_MODULUS)
        cls.model_part.AddNodalSolutionStepVariable(KM.NODAL_VAUX)
        cls.model_part.AddNodalSolutionStepVariable(KM.EXTERNAL_FORCES_VECTOR)
        cls.model_part.AddNodalSolutionStepVariable(KM.LOCAL_AXES_MATRIX)

        # If distributed run, add PARTITION_INDEX to model part
        if KM.IsDistributedRun():
            cls.model_part.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)
        
        # Define mdpa file name and read it into the model part
        cls.mdpa_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere_with_conditions")
        ReadModelPart(cls.mdpa_name, cls.model_part)

    @classmethod
    def tearDownClass(cls):
        """
        Cleanup after all tests are run.
        """
        RemoveFiles(cls.mdpa_name)

    def setUp(self):
        """
        Setup for each individual test.
        """
        # Create search
        self.data_comm = self.model_part.GetCommunicator().GetDataCommunicator()

        # Create node for search
        self.node_id = 100
        if KM.IsDistributedRun():
            # Only added to first rank to actually check it works in all ranks
            if self.data_comm.Rank() == 0:
                self.node = self.sub_model_part.CreateNewNode(self.node_id, 0.0, 0.0, 0.15)
                self.node.SetSolutionStepValue(KM.PARTITION_INDEX, 0)
            ParallelFillCommunicator = KratosMPI.ParallelFillCommunicator(self.model_part, self.data_comm)
            ParallelFillCommunicator.Execute()
        else:
            self.node = self.sub_model_part.CreateNewNode(self.node_id, 0.0, 0.0, 0.15)

    def _create_search(self, search_type):
        """
        Create an appropriate search wrapper based on the given search type.
        """
        # if search_type == "Octree":
        #     self.search = KM.SearchWrapperOctree(self.model_part.Conditions, self.data_comm)
        # elif search_type == "KDTree":
        #     self.search = KM.SearchWrapperKDTree(self.model_part.Conditions, self.data_comm)
        # elif search_type == "BinBased":
        #     self.search = KM.SearchWrapperBinBased(self.model_part.Conditions, self.data_comm)
        # el
        if search_type == "GeometricalObjectBins":
            self.search = KM.SearchWrapperGeometricalObjectBins(self.model_part.Conditions, self.data_comm)
        else:
            raise Exception("Invalid search type: " + search_type)

    def test_SearchWrapperGeometricalObjectBins_SearchInRadius(self):
        """
        Test for the 'SearchInRadius' method of the GeometricalObjectBins search wrapper.
        """
        # Create search
        self._create_search("GeometricalObjectBins")

        # Define radius
        radius = 0.35

        # Reference solution
        cond_id_ref = [125,78,117,18,68,1,41,119]

        # Nodes array search
        results = self.search.SearchInRadius(self.sub_model_part.Nodes, radius)
        self.assertEqual(results.NumberOfSearchResults(), 1)
        node_results = results[self.node_id]
        self.assertEqual(node_results.NumberOfGlobalResults(), 8)
        ids = node_results.GetResultIndices()
        self.assertEqual(len(ids), 8)
        for id in ids:
            self.assertTrue(id in cond_id_ref)

    def test_SearchWrapperGeometricalObjectBins_SearchNearestInRadius(self):
        """
        Test for the 'SearchNearestInRadius' method of the GeometricalObjectBins search wrapper.
        """
        # Create search
        self._create_search("GeometricalObjectBins")

        # Define radius
        radius = 0.35

        # Nodes array search
        results = self.search.SearchNearestInRadius(self.sub_model_part.Nodes, radius)   
        self.assertEqual(results.NumberOfSearchResults(), 1)
        node_results = results[self.node_id]
        self.assertEqual(node_results.NumberOfGlobalResults(), 1)
        # Local result
        if node_results.NumberOfLocalResults() == 1:
            self.assertTrue(node_results[0].IsObjectFound())
        # Global result
        ids = node_results.GetResultIndices()
        self.assertEqual(len(ids), 1)
        self.assertTrue(1 in ids)

    def test_SearchWrapperGeometricalObjectBins_SearchNearest(self):
        """
        Test for the 'SearchNearest' method of the GeometricalObjectBins search wrapper.
        """
        # Create search
        self._create_search("GeometricalObjectBins")

        # Nodes array search
        results = self.search.SearchNearest(self.sub_model_part.Nodes) 
        self.assertEqual(results.NumberOfSearchResults(), 1)
        node_results = results[self.node_id] 
        self.assertEqual(node_results.NumberOfGlobalResults(), 1)
        # Local result
        if node_results.NumberOfLocalResults() == 1:
            self.assertTrue(node_results[0].IsObjectFound())
        # Global result
        ids = node_results.GetResultIndices()
        self.assertEqual(len(ids), 1)
        self.assertTrue(1 in ids)

    def test_SearchWrapperGeometricalObjectBins_SearchIsInside(self):
        """
        Test for the 'SearchIsInside' method of the GeometricalObjectBins search wrapper.
        """       
        # Create search
        self._create_search("GeometricalObjectBins")

        # Nodes array search
        results = self.search.SearchIsInside(self.sub_model_part.Nodes) 
        self.assertEqual(results.NumberOfSearchResults(), 1)
        node_results = results[self.node_id] 
        self.assertFalse(node_results.IsObjectFound())

if __name__ == '__main__':
    # Set logging severity and start the unittest
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()