//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// System includes
#include <limits>

// External includes

// Project includes
#include "testing/testing.h"
#include "containers/model.h"
#include "spatial_containers/search_wrapper.h"
#include "utilities/cpp_tests_utilities.h"
#include "mpi/utilities/parallel_fill_communicator.h"

namespace Kratos::Testing 
{

// Definition of the geometrical object bins search wrapper
using SearchWrapperGeometricalObjectsBins = SearchWrapper<GeometricalObjectsBins, GeometricalObject>;

ModelPart& CreateCubeSkinModelPart(
    Model& rCurrentModel,
    const double HalfX = 0.6,
    const double HalfY = 0.9,
    const double HalfZ = 0.3
    )
{
    // Generate the cube skin
    const DataCommunicator& r_data_communicator = Testing::GetDefaultDataCommunicator();
    ModelPart& r_skin_part = CppTestsUtilities::CreateCubeSkinModelPart(rCurrentModel, HalfX, HalfY, HalfZ, r_data_communicator);

    // Compute communication plan and fill communicator meshes correctly
    ParallelFillCommunicator(r_skin_part, r_data_communicator).Execute();

    // Return the skin model part
    return r_skin_part;
}

ModelPart& CreateCubeModelPart(Model& rCurrentModel)
{
    // Generate the cube
    const DataCommunicator& r_data_communicator = Testing::GetDefaultDataCommunicator();
    ModelPart& r_model_part = CppTestsUtilities::CreateCubeModelPart(rCurrentModel, r_data_communicator);

    // Compute communication plan and fill communicator meshes correctly
    ParallelFillCommunicator(r_model_part, r_data_communicator).Execute();

    // Return the model part
    return r_model_part;
}

/** Checks search_wrapper_bins search in radius
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperGeometricalObjectsBinsSearchInRadius, KratosMPICoreFastSuite)
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    SearchWrapperGeometricalObjectsBins search_wrapper_bins(r_skin_part.Elements(), r_data_comm);

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    r_point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    const std::size_t point_id = 1;
    if (rank == 0) {
        auto p_node = r_point_model_part.CreateNewNode(point_id, 0.0, 0.0, 0.0);
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    SearchWrapperGeometricalObjectsBins::ResultContainerVectorType results;

    // 0.29 radius
    search_wrapper_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.29, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[point_id].IsObjectFound());
    KRATOS_EXPECT_EQ(results[point_id].NumberOfGlobalResults(), 0);

    // 0.3 radius
    search_wrapper_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.3, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[point_id].IsObjectFound());
    KRATOS_EXPECT_EQ(results[point_id].NumberOfGlobalResults(), 4);

    // 0.4 radius
    search_wrapper_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.4, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[point_id].IsObjectFound());
    KRATOS_EXPECT_EQ(results[point_id].NumberOfGlobalResults(), 4);

    // 0.6 radius
    search_wrapper_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.6, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[point_id].IsObjectFound());
    KRATOS_EXPECT_EQ(results[point_id].NumberOfGlobalResults(), 8);

    // 0.7 radius
    search_wrapper_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.7, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[point_id].IsObjectFound());
    KRATOS_EXPECT_EQ(results[point_id].NumberOfGlobalResults(), 8);

    // 0.9 radius
    search_wrapper_bins.SearchInRadius(r_array_nodes.begin(), r_array_nodes.end(), 0.9, results);
    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[point_id].IsObjectFound());
    KRATOS_EXPECT_EQ(results[point_id].NumberOfGlobalResults(), 12);
}

/** Checks search_wrapper_bins search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperGeometricalObjectsBinsSearchNearestInRadius, KratosMPICoreFastSuite) 
{
    constexpr double tolerance = 1e-12;

    Model current_model;

    // Cube coordinates
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    SearchWrapperGeometricalObjectsBins search_wrapper_bins(r_skin_part.Elements(), r_data_comm);

    double epsilon = 1.0e-6;
    const std::size_t near_point_id = 1;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    r_point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    if (rank == 0) {
        auto p_node = r_point_model_part.CreateNewNode(near_point_id, epsilon,epsilon,epsilon);
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    SearchWrapperGeometricalObjectsBins::ResultContainerVectorType results;
    search_wrapper_bins.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z - 1.e-4, results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[near_point_id].IsObjectFound());

    search_wrapper_bins.SearchNearestInRadius(r_array_nodes.begin(), r_array_nodes.end(), cube_z + 1.e-4, results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[near_point_id].IsObjectFound());
    KRATOS_EXPECT_EQ(results[near_point_id].NumberOfGlobalResults(), 1);

    // Distances are just local
    const auto distances = results[near_point_id].GetDistances();
    KRATOS_EXPECT_NEAR(distances[0], (cube_z - epsilon), tolerance);

    // Compute indices
    auto indices = results[near_point_id].GetResultIndices();
    const std::size_t id = indices[0];
    KRATOS_EXPECT_EQ(id, 3);
}

/** Checks search_wrapper_bins search nearest
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperGeometricalObjectsBinsSearchNearest, KratosMPICoreFastSuite) 
{
    constexpr double tolerance = 1e-12;

    Model current_model;
    
    // Cube coordinates
    const double cube_z = 0.3;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeSkinModelPart(current_model, 0.6, 0.9, cube_z);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    SearchWrapperGeometricalObjectsBins search_wrapper_bins(r_skin_part.Elements(), r_data_comm);

    double epsilon = 1.0e-6;
    const std::size_t near_point_id = 1;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    r_point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    if (rank == 0) {
        auto p_node = r_point_model_part.CreateNewNode(near_point_id, epsilon,epsilon,epsilon);
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    SearchWrapperGeometricalObjectsBins::ResultContainerVectorType results;
    search_wrapper_bins.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[near_point_id].IsObjectFound());
    KRATOS_EXPECT_EQ(results[near_point_id].NumberOfGlobalResults(), 1);

    // Distances are just local
    const auto distances = results[near_point_id].GetDistances();
    KRATOS_EXPECT_NEAR(distances[0], (cube_z - epsilon), tolerance);

    // Compute indices
    auto indices = results[near_point_id].GetResultIndices();
    const std::size_t id = indices[0];
    KRATOS_EXPECT_EQ(id, 3);
}

/** Checks search_wrapper_bins empty search nearest 
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperGeometricalObjectsBinsEmptySearchNearest, KratosMPICoreFastSuite) 
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = current_model.CreateModelPart("Skin");
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    SearchWrapperGeometricalObjectsBins search_wrapper_bins(r_skin_part.Elements(), r_data_comm);

    const std::size_t point_id = 1;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    r_point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    if (rank == 0) {
        auto p_node = r_point_model_part.CreateNewNode(point_id, 0.0,0.0,0.0);
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    SearchWrapperGeometricalObjectsBins::ResultContainerVectorType results;
    search_wrapper_bins.SearchNearest(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[point_id].IsObjectFound());
}

/** Checks search_wrapper_bins search is inside 
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperGeometricalObjectsBinsSearchIsInside, KratosMPICoreFastSuite) 
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeModelPart(current_model);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    SearchWrapperGeometricalObjectsBins search_wrapper_bins(r_skin_part.Elements(), r_data_comm);

    const std::size_t inside_point_id = 1;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    r_point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    if (rank == 0) {
        auto p_node = r_point_model_part.CreateNewNode(inside_point_id, 0.5,0.5,0.5);
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    SearchWrapperGeometricalObjectsBins::ResultContainerVectorType results;
    search_wrapper_bins.SearchIsInside(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_TRUE(results[inside_point_id].IsObjectFound());
    KRATOS_EXPECT_EQ(results[inside_point_id].NumberOfGlobalResults(), 1);
}

/** Checks search_wrapper_bins search is inside = not found
*/
KRATOS_DISTRIBUTED_TEST_CASE_IN_SUITE(MPISearchWrapperGeometricalObjectsBinsSearchIsNotInside, KratosMPICoreFastSuite) 
{
    Model current_model;

    // Generate the cube skin
    ModelPart& r_skin_part = CreateCubeModelPart(current_model);
    const DataCommunicator& r_data_comm = Testing::GetDefaultDataCommunicator();

    SearchWrapperGeometricalObjectsBins search_wrapper_bins(r_skin_part.Elements(), r_data_comm);

    const std::size_t outside_point_id = 1;

    // Generate new model part
    ModelPart& r_point_model_part = current_model.CreateModelPart("PointModelPart");
    r_point_model_part.AddNodalSolutionStepVariable(PARTITION_INDEX);

    // We generate only in first rank
    const int rank = r_data_comm.Rank();
    if (rank == 0) {
        auto p_node = r_point_model_part.CreateNewNode(outside_point_id, 100.0,100.0,100.0);
        p_node->FastGetSolutionStepValue(PARTITION_INDEX) = 0;
    }
    auto& r_array_nodes = r_point_model_part.Nodes();

    SearchWrapperGeometricalObjectsBins::ResultContainerVectorType results;
    search_wrapper_bins.SearchIsInside(r_array_nodes.begin(), r_array_nodes.end(), results);

    KRATOS_EXPECT_EQ(results.NumberOfSearchResults(), 1);
    KRATOS_EXPECT_FALSE(results[outside_point_id].IsObjectFound());
}

} // namespace Kratos::Testing.