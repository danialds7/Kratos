// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Franziska Wahl
//

// System includes
#include <unordered_set>

// External includes

// Project includes
#include "processes/find_nodal_neighbours_process.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "embedded_local_constraint_process.h"

namespace Kratos
{

using MLSShapeFunctionsFunctionType = std::function<void(const Matrix&, const array_1d<double,3>&, const double, Vector&)>;

/* Public functions *******************************************************/

EmbeddedLocalConstraintProcess::EmbeddedLocalConstraintProcess(
    Model& rModel,
    Parameters ThisParameters)
    : Process()
{
    // Validate input settings with defaults
    ThisParameters.ValidateAndAssignDefaults(GetDefaultParameters());

    // Retrieve the required model parts
    const std::string model_part_name = ThisParameters["model_part_name"].GetString();
    mpModelPart = &rModel.GetModelPart(model_part_name);

    // Set how constraints are applied
    mApplyToAllNegativeCutNodes = ThisParameters["apply_to_all_negative_cut_nodes"].GetBool();
    mUseMLSShapeFunctions = ThisParameters["use_mls_shape_functions"].GetBool();
    mIncludeIntersectionPoints = ThisParameters["include_intersection_points"].GetBool();

    // Set whether nodal distances will be modified to avoid levelset zeros
    mAvoidZeroDistances = ThisParameters["avoid_zero_distances"].GetBool();

    // Set which elements will not be active
    mDeactivateNegativeElements = ThisParameters["deactivate_negative_elements"].GetBool();
    mDeactivateIntersectedElements = ThisParameters["deactivate_intersected_elements"].GetBool();
}

void EmbeddedLocalConstraintProcess::Execute()
{
    // Declare extension operator map for negative nodes of split elements
    NodesCloudMapType clouds_map;
    NodesOffsetMapType offsets_map;

    if (mAvoidZeroDistances) { ModifyDistances(); }

    // Deactivate intersected and negative elements as well as their nodes depending on the process settings
    DeactivateElementsAndNodes();

    // Calculate the negative nodes' clouds
    CalculateNodeClouds(clouds_map, offsets_map);

    // Apply a constraint to negative nodes of split elements
    ApplyConstraints(clouds_map, offsets_map);
}

/* Protected functions ****************************************************/

/* Private functions ****************************************************/

void EmbeddedLocalConstraintProcess::CalculateNodeClouds(NodesCloudMapType& rCloudsMap, NodesOffsetMapType& rOffsetsMap)
{
    //Pointer for which function to call
    void(Kratos::EmbeddedLocalConstraintProcess::*fptrNodeClouds)(NodesCloudMapType&, NodesOffsetMapType&, std::vector<NodeType::Pointer>, std::vector<NodeType::Pointer>);
    if (mUseMLSShapeFunctions) {
        if (mIncludeIntersectionPoints) {
            fptrNodeClouds = &Kratos::EmbeddedLocalConstraintProcess::AddMLSNodeCloudsIncludingBC;
        } else {
            fptrNodeClouds = &Kratos::EmbeddedLocalConstraintProcess::AddMLSNodeClouds;
        }
    } else {
        if (mIncludeIntersectionPoints) {
            fptrNodeClouds = &Kratos::EmbeddedLocalConstraintProcess::AddAveragedNodeCloudsIncludingBC;
        } else {
            fptrNodeClouds = &Kratos::EmbeddedLocalConstraintProcess::AddAveragedNodeClouds;
        }
    }

    // Containers for negative and positive side nodes of split element
    std::vector<std::vector<NodeType::Pointer>> neg_nodes = {};
    std::vector<std::vector<NodeType::Pointer>> pos_nodes = {};

    // Loop through all elements to get negative and positive nodes of split elements
    for (auto& rElement : mpModelPart->Elements()) {
    auto& r_geom = rElement.GetGeometry();
        if (IsSplit(r_geom)) {
            if (mApplyToAllNegativeCutNodes || IsSmallCut(r_geom)) {
                // Containers for negative and positive side nodes of split element
                std::vector<NodeType::Pointer> neg_nodes_element = {};
                std::vector<NodeType::Pointer> pos_nodes_element = {};
                for (auto& rNode : r_geom) {
                    if (rNode.FastGetSolutionStepValue(DISTANCE) > 0.0) {
                        pos_nodes_element.push_back(&rNode);
                    } else {
                        neg_nodes_element.push_back(&rNode);
                    }
                }
                pos_nodes.push_back(pos_nodes_element);
                neg_nodes.push_back(neg_nodes_element);
            }
        }
    }

    //TODO delete "used" entries for faster looping?
    // Loop through cut elements to add single negative nodes first
    for(std::size_t it = 0; it < neg_nodes.size(); ++it) {
        if (neg_nodes[it].size() == 1) {
            (this->*fptrNodeClouds)(rCloudsMap, rOffsetsMap, neg_nodes[it], pos_nodes[it]);
        }
    }
    // Loop for two negative nodes of cut elements
    for(std::size_t it = 0; it < neg_nodes.size(); ++it) {
        if (neg_nodes[it].size() == 2) {
            (this->*fptrNodeClouds)(rCloudsMap, rOffsetsMap, neg_nodes[it], pos_nodes[it]);
        }
    }
    // Loop for three negative nodes of cut elements if 3D
    if (mpModelPart->GetProcessInfo()[DOMAIN_SIZE] == 3) {
        for(std::size_t it = 0; it < neg_nodes.size(); ++it) {
            if (neg_nodes[it].size() == 3) {
                (this->*fptrNodeClouds)(rCloudsMap, rOffsetsMap, neg_nodes[it], pos_nodes[it]);
            }
        }
    }
}

void EmbeddedLocalConstraintProcess::AddAveragedNodeClouds(NodesCloudMapType& rCloudsMap, NodesOffsetMapType& rOffsetsMap, std::vector<NodeType::Pointer> neg_nodes_element, std::vector<NodeType::Pointer> pos_nodes_element)
{

    for (auto slave_node : neg_nodes_element) {
        // Check whether node was already added because of a previous element
        if (!rCloudsMap.count(slave_node)) {
            // Get number of positive element nodes
            const std::size_t n_cloud_nodes = pos_nodes_element.size();
            // Calculate weight as averaging of all master nodes
            const double cloud_node_weight = 1.0 / n_cloud_nodes;

            // Save positive element nodes (cloud nodes) and weights (averaged)
            CloudDataVectorType cloud_data_vector(n_cloud_nodes);
            for (std::size_t i_pos_node = 0; i_pos_node < n_cloud_nodes; ++i_pos_node) {
                auto i_data = std::make_pair(pos_nodes_element[i_pos_node], cloud_node_weight);
                cloud_data_vector(i_pos_node) = i_data;
            }

            // Pair slave node and constraint offset as well as slave node and cloud vector and add to respective map
            rOffsetsMap.insert(std::make_pair(slave_node, 0.0));
            rCloudsMap.insert(std::make_pair(slave_node, cloud_data_vector));
            KRATOS_WATCH("Added local averaged cloud for a negative node of a cut element");
        }
    }
}

void EmbeddedLocalConstraintProcess::AddAveragedNodeCloudsIncludingBC(NodesCloudMapType& rCloudsMap, NodesOffsetMapType& rOffsetsMap, std::vector<NodeType::Pointer> neg_nodes_element, std::vector<NodeType::Pointer> pos_nodes_element)
{
    for (auto slave_node : neg_nodes_element) {
        // Check whether node was already added because of a previous element
        if (!rCloudsMap.count(slave_node)) {
            // Get number of positive element nodes
            const std::size_t n_cloud_nodes = pos_nodes_element.size();

            // Get nodal distances
            const double dist_slave = slave_node->FastGetSolutionStepValue(DISTANCE);
            double sum_dist_cloud_nodes = 0.0;
            for (auto pos_node : pos_nodes_element) {
                sum_dist_cloud_nodes += pos_node->FastGetSolutionStepValue(DISTANCE);
            }

            // Calculate averaged coordinates of intersection points
            array_1d<double,3> avg_int_pt_coord;
            avg_int_pt_coord[0] = 0.0;
            avg_int_pt_coord[1] = 0.0;
            avg_int_pt_coord[2] = 0.0;

            std::size_t n_intersections = 0;
            for (auto neg_node : neg_nodes_element) {
                double dist_neg_node = neg_node->FastGetSolutionStepValue(DISTANCE);
                for (auto pos_node : pos_nodes_element) {
                    double dist_pos_node = pos_node->FastGetSolutionStepValue(DISTANCE);
                    // 2D: Calculate intersection point
                    const double int_pt_rel_location = std::abs(dist_neg_node/(dist_pos_node-dist_neg_node));
                    avg_int_pt_coord += *neg_node + int_pt_rel_location * (*pos_node - *neg_node);
                    n_intersections++;
                }
            }
            avg_int_pt_coord /= n_intersections;

            // Get boundary condition value for averaged intersection point  // TODO  different boundary condition?!?
            //const double temp_bc = rElement.GetValue(EMBEDDED_SCALAR);  // TODO =0.0 if value_bc = 0.0 ???
            const double value_bc = std::pow(avg_int_pt_coord[0],2) + std::pow(avg_int_pt_coord[1],2);

            //Calculate offset for master-slave constraint
            double offset = value_bc - dist_slave * value_bc * n_cloud_nodes / sum_dist_cloud_nodes;

            // Calculate weight
            double cloud_node_weight = dist_slave / sum_dist_cloud_nodes;
            KRATOS_WATCH(cloud_node_weight);

            //TODO: case: small cut for positive side - change!!!  // TODO: scale slave_distance using edge length/elements size?
            if (sum_dist_cloud_nodes < 1e-10) {
                offset = value_bc;
                cloud_node_weight = 0.0;
            }

            // Save positive element nodes (cloud nodes) and weights (linear interpolation of averaged positive nodes and averaged intersection points)
            CloudDataVectorType cloud_data_vector(n_cloud_nodes);
            for (std::size_t i_pos_node = 0; i_pos_node < n_cloud_nodes; ++i_pos_node) {
                auto i_data = std::make_pair(pos_nodes_element[i_pos_node], cloud_node_weight);
                cloud_data_vector(i_pos_node) = i_data;
            }

            // Pair slave node and constraint offset as well as slave node and cloud vector and add to respective map
            rOffsetsMap.insert(std::make_pair(slave_node, offset));
            rCloudsMap.insert(std::make_pair(slave_node, cloud_data_vector));
            KRATOS_WATCH("Added local averaged cloud including intersection points for a negative node of a cut element");
        }
    }
}

void EmbeddedLocalConstraintProcess::AddMLSNodeClouds(NodesCloudMapType& rCloudsMap, NodesOffsetMapType& rOffsetsMap, std::vector<NodeType::Pointer> neg_nodes_element, std::vector<NodeType::Pointer> pos_nodes_element)
{
    for (auto slave_node : neg_nodes_element) {
        // Check whether node was already added because of a previous element
        if (!rCloudsMap.count(slave_node)) {
            // Get the MLS shape functions function
            auto p_mls_sh_func = GetMLSShapeFunctionsFunction();

            // Get number of positive element nodes
            const std::size_t n_cloud_nodes = pos_nodes_element.size();

            // Get the current negative node support cloud coordinates (positive nodes of cut element)
            Matrix cloud_nodes_coordinates;
            cloud_nodes_coordinates.resize(n_cloud_nodes, 3);
            IndexPartition<std::size_t>(n_cloud_nodes).for_each(array_1d<double,3>(), [&pos_nodes_element, &cloud_nodes_coordinates](std::size_t i_pos_node, array_1d<double,3>& rAuxCoordTLS){
                noalias(rAuxCoordTLS) = pos_nodes_element[i_pos_node]->Coordinates();
                cloud_nodes_coordinates(i_pos_node, 0) = rAuxCoordTLS[0];
                cloud_nodes_coordinates(i_pos_node, 1) = rAuxCoordTLS[1];
                cloud_nodes_coordinates(i_pos_node, 2) = rAuxCoordTLS[2];
            });

            // Calculate the MLS basis in the current negative node
            Vector N_container;
            const array_1d<double,3> r_coords = slave_node->Coordinates();
            const double mls_kernel_rad = CalculateKernelRadius(cloud_nodes_coordinates, r_coords);
            p_mls_sh_func(cloud_nodes_coordinates, r_coords, 1.01 * mls_kernel_rad, N_container);

            /// Save positive element nodes (cloud nodes) and weights (MLS shape function values)  // =0.0 if value_bc = 0.0 ???
            CloudDataVectorType cloud_data_vector(n_cloud_nodes);
            for (std::size_t i_pos_node = 0; i_pos_node < n_cloud_nodes; ++i_pos_node) {
                auto i_data = std::make_pair(pos_nodes_element[i_pos_node], N_container[i_pos_node]);
                cloud_data_vector(i_pos_node) = i_data;
            }

            // Pair slave node and constraint offset as well as slave node and cloud vector and add to respective map
            rOffsetsMap.insert(std::make_pair(slave_node, 0.0));
            rCloudsMap.insert(std::make_pair(slave_node, cloud_data_vector));
            KRATOS_WATCH("Added local MLS cloud for a negative node of a cut element");
        }
    }
}

void EmbeddedLocalConstraintProcess::AddMLSNodeCloudsIncludingBC(NodesCloudMapType& rCloudsMap, NodesOffsetMapType& rOffsetsMap, std::vector<NodeType::Pointer> neg_nodes_element, std::vector<NodeType::Pointer> pos_nodes_element)
{
    for (auto slave_node : neg_nodes_element) {
        // Check whether node was already added because of a previous element
        if (!rCloudsMap.count(slave_node)) {
            // Get the MLS shape functions function
            auto p_mls_sh_func = GetMLSShapeFunctionsFunction();

            // Get number of positive element nodes
            const std::size_t n_cloud_nodes = pos_nodes_element.size();

            // Calculate intersection points
            std::vector<array_1d<double,3>> int_pt_coord;
            for (auto neg_node : neg_nodes_element) {
                double dist_neg_node = neg_node->FastGetSolutionStepValue(DISTANCE);
                for (auto pos_node : pos_nodes_element) {
                    double dist_pos_node = pos_node->FastGetSolutionStepValue(DISTANCE);
                    const double int_pt_rel_location = std::abs(dist_neg_node/(dist_pos_node-dist_neg_node));
                    array_1d<double,3> int_pt = *neg_node + int_pt_rel_location * (*pos_node - *neg_node);
                    int_pt_coord.push_back(int_pt);
                }
            }
            const std::size_t n_intersections = int_pt_coord.size();

            // Get the current negative node support cloud coordinates (positive nodes of cut element)
            Matrix cloud_nodes_coordinates;
            cloud_nodes_coordinates.resize(n_cloud_nodes+n_intersections, 3);
            IndexPartition<std::size_t>(n_cloud_nodes).for_each(array_1d<double,3>(), [&pos_nodes_element, &cloud_nodes_coordinates](std::size_t i_pos_node, array_1d<double,3>& rAuxCoordTLS){
                noalias(rAuxCoordTLS) = pos_nodes_element[i_pos_node]->Coordinates();
                cloud_nodes_coordinates(i_pos_node, 0) = rAuxCoordTLS[0];
                cloud_nodes_coordinates(i_pos_node, 1) = rAuxCoordTLS[1];
                cloud_nodes_coordinates(i_pos_node, 2) = rAuxCoordTLS[2];
            });

            //Add intersection points coordinates to cloud
            for (std::size_t i_int_pt = n_cloud_nodes; i_int_pt < n_cloud_nodes+n_intersections; ++i_int_pt) {
                cloud_nodes_coordinates(i_int_pt, 0) = int_pt_coord[i_int_pt][0];
                cloud_nodes_coordinates(i_int_pt, 1) = int_pt_coord[i_int_pt][1];
                cloud_nodes_coordinates(i_int_pt, 2) = int_pt_coord[i_int_pt][2];
            }

            // Calculate the MLS basis in the current negative node
            Vector N_container;
            const array_1d<double,3> r_coords = slave_node->Coordinates();
            const double mls_kernel_rad = CalculateKernelRadius(cloud_nodes_coordinates, r_coords);
            p_mls_sh_func(cloud_nodes_coordinates, r_coords, 1.01 * mls_kernel_rad, N_container);

            /// Save positive element nodes (cloud nodes) and weights (MLS shape function values)
            CloudDataVectorType cloud_data_vector(n_cloud_nodes);
            for (std::size_t i_pos_node = 0; i_pos_node < n_cloud_nodes; ++i_pos_node) {
                auto i_data = std::make_pair(pos_nodes_element[i_pos_node], N_container[i_pos_node]);
                cloud_data_vector(i_pos_node) = i_data;
            }

            // Calculate offset for master-slave constraint
            double offset = 0.0;
            for (std::size_t i_int_pt = n_cloud_nodes; i_int_pt < n_cloud_nodes+n_intersections; ++i_int_pt) {
                // Get boundary condition value for intersection point  // TODO  different boundary condition?!?
                //const double temp_bc = rElement.GetValue(EMBEDDED_SCALAR);
                const double value_bc = std::pow(int_pt_coord[i_int_pt][0],2) + std::pow(int_pt_coord[i_int_pt][1],2);
                offset += N_container[i_int_pt] * value_bc;
            }

            // Pair slave node and constraint offset as well as slave node and cloud vector and add to respective map
            rOffsetsMap.insert(std::make_pair(slave_node, offset));
            rCloudsMap.insert(std::make_pair(slave_node, cloud_data_vector));
            KRATOS_WATCH("Added local MLS cloud including intersection points for a negative node of a cut element");
        }
    }
}

double EmbeddedLocalConstraintProcess::CalculateKernelRadius(
    const Matrix& rCloudCoordinates,
    const array_1d<double,3>& rOrigin)
{
    const std::size_t n_nodes = rCloudCoordinates.size1();
    const double squared_rad = IndexPartition<std::size_t>(n_nodes).for_each<MaxReduction<double>>([&](std::size_t I){
        return std::pow(rCloudCoordinates(I,0) - rOrigin(0), 2) + std::pow(rCloudCoordinates(I,1) - rOrigin(1), 2) + std::pow(rCloudCoordinates(I,2) - rOrigin(2), 2);
    });
    return std::sqrt(squared_rad);
}

void EmbeddedLocalConstraintProcess::ApplyConstraints(NodesCloudMapType& rCloudsMap, NodesOffsetMapType& rOffsetsMap)
{
    // Initialize counter of master slave constraints
    ModelPart::IndexType id = mpModelPart->NumberOfMasterSlaveConstraints()+1;
    // Get variable to constrain  //TODO: different variable??
    const auto& r_var = KratosComponents<Variable<double>>::Get("TEMPERATURE");

    // Loop through all negative nodes of split elements (slave nodes)
    for (auto it_slave = rCloudsMap.begin(); it_slave != rCloudsMap.end(); ++it_slave) {
        auto p_slave_node = std::get<0>(*it_slave);
        auto& r_cloud_data = rCloudsMap[p_slave_node];
        double offset = rOffsetsMap[p_slave_node];

        // Add one master slave constraint for every node of the support cloud (master) of the negative node (slave)
        // The contributions of each master will be summed up in the BuilderAndSolver to give an equation for the slave dof
        for (auto it_data = r_cloud_data.begin(); it_data != r_cloud_data.end(); ++it_data) {
            auto& r_node_data = *it_data;
            auto p_support_node = std::get<0>(r_node_data);
            const double support_node_N = std::get<1>(r_node_data);

            // Add master slave constraint, the support node MLS shape function value N serves as weight of the constraint
            mpModelPart->CreateNewMasterSlaveConstraint("LinearMasterSlaveConstraint", id++,
            *p_support_node, r_var, *p_slave_node, r_var,
            support_node_N, offset);
            offset = 0.0;
        }
    }
}

void EmbeddedLocalConstraintProcess::DeactivateElementsAndNodes()
{
    // Initialize flags to true //TODO: necessary??
    block_for_each(mpModelPart->Nodes(), [](NodeType& rNode){
        rNode.Set(ACTIVE, true);  // Nodes that belong to the elements to be assembled
    });
    block_for_each(mpModelPart->Elements(), [](Element& rElement){
        rElement.Set(ACTIVE, true);  // Elements in the positive distance region (the ones to be assembled)
    });

    // Deactivate intersected elements and their nodes  //TODO: move to where the constraints are applied??
    if ( mDeactivateIntersectedElements ) {
        for (auto& rElement : mpModelPart->Elements()) {
            // Check if the element is split
            const auto& r_geom = rElement.GetGeometry();
            if (IsSplit(r_geom)) {
                rElement.Set(ACTIVE, false);
                for (auto& rNode : r_geom) {
                    rNode.Set(ACTIVE, false);
                }
            }
        }
    }
    // Deactivate negative elements and their nodes
    if ( mDeactivateNegativeElements ) {
        for (auto& rElement : mpModelPart->Elements()) {
            // Check if the element is negative
            const auto& r_geom = rElement.GetGeometry();
            if (IsNegative(r_geom)) {
                rElement.Set(ACTIVE, false);
                for (auto& rNode : r_geom) {
                    rNode.Set(ACTIVE, false);
                }
            }
        }
    }
}

void EmbeddedLocalConstraintProcess::ModifyDistances()
{
    auto& r_nodes = mpModelPart->Nodes();
    const double tol_d = 1.0e-12;

    #pragma omp parallel for
    for (int i = 0; i < static_cast<int>(r_nodes.size()); ++i) {
        auto it_node = r_nodes.begin() + i;
        double& d = it_node->FastGetSolutionStepValue(DISTANCE);

        // Check if the distance values are close to zero, if so set the tolerance as distance value
        if (std::abs(d) < tol_d) {
            d = (d > 0.0) ? tol_d : -tol_d;
        }
    }
}

bool EmbeddedLocalConstraintProcess::IsSplit(const GeometryType& rGeometry)
{
    std::size_t n_neg = 0;
    std::size_t n_pos = 0;
    for (const auto& r_node : rGeometry) {
        if (r_node.FastGetSolutionStepValue(DISTANCE) < 0.0) {
            n_neg++;
        } else {
            n_pos++;
        }
    }
    return (n_pos != 0 && n_neg != 0);
}

bool EmbeddedLocalConstraintProcess::IsSmallCut(const GeometryType& rGeometry)
{
    const double tol_d = 0.001;
    for (const auto& r_node : rGeometry) {
        if (abs(r_node.FastGetSolutionStepValue(DISTANCE)) < tol_d) {
            return true;
        }
    }
    return false;
}

bool EmbeddedLocalConstraintProcess::IsNegative(const GeometryType& rGeometry)
{
    std::size_t n_neg = 0;
    for (const auto& r_node : rGeometry) {
        if (r_node.FastGetSolutionStepValue(DISTANCE) < 0.0) {
            n_neg++;
        }
    }
    return (n_neg == rGeometry.PointsNumber());
}

MLSShapeFunctionsFunctionType EmbeddedLocalConstraintProcess::GetMLSShapeFunctionsFunction()
{
    switch (mpModelPart->GetProcessInfo()[DOMAIN_SIZE]) {
        case 2:
            return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN){
                MLSShapeFunctionsUtility::CalculateShapeFunctions<2,1>(rPoints, rX, h, rN);};
        case 3:
            return [&](const Matrix& rPoints, const array_1d<double,3>& rX, const double h, Vector& rN){
                MLSShapeFunctionsUtility::CalculateShapeFunctions<3,1>(rPoints, rX, h, rN);};
        default:
            KRATOS_ERROR << "Wrong domain size. MLS shape functions utility cannot be set.";
    }
}

}; // namespace Kratos.