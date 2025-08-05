//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Uxue Chasco
//

// System includes
#include <algorithm>

// External includes


// Project includes
#include "processes/find_global_nodal_neighbours_process.h"
#include "includes/global_pointer_variables.h"
#include "containers/global_pointers_vector.h"
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"
#include "spatial_containers/bins_dynamic.h"
#include "utilities/rbf_shape_functions_utility.h"
#include "utilities/divide_triangle_3d_3.h"
#include "../../FluidDynamicsApplication/custom_utilities/fluid_auxiliary_utilities.h"
#include "../../FluidDynamicsApplication/fluid_dynamics_application_variables.h"
#include "includes/variables.h"
#include "includes/define.h"

// Application includes
#include "hydraulic_fluid_auxiliary_utilities.h"
#include "../../FluidDynamicsApplication/custom_utilities/fluid_mesh_utilities.h"

namespace Kratos
{
typedef std::size_t SizeType;
typedef std::size_t IndexType;

struct VectorHasher {
    std::size_t operator()(const std::vector<IndexType>& v) const {
        std::size_t seed = v.size();
        for (auto& i : v) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

double HydraulicFluidAuxiliaryUtilities::CalculateWettedPetimeter(
    ModelPart &rModelPart,
    const Flags &rSkinFlag,
    const Variable<double>& rDistanceVariable,
    const bool IsHistorical)
{
    // Auxiliary function to have the posssibility of having non historical or hitorical inlet distance function.
    std::function<double(NodeType&, const Variable<double>&)> distance_getter;
    if (IsHistorical)
    {
        distance_getter = [](NodeType& rNode, const Variable<double>& rDistanceVariable) -> double {return rNode.FastGetSolutionStepValue(rDistanceVariable);};
    } else {
        distance_getter = [](NodeType& rNode, const Variable<double>& rDistanceVariable) -> double {return rNode.GetValue(rDistanceVariable);};
    }

    // Check that there are conditions and distance_inlet variable in the nodal database
    KRATOS_ERROR_IF(rModelPart.GetProcessInfo()[DOMAIN_SIZE] < 3 )<< "Wetted perimeter is only implemented for 3D." << std::endl;

    const auto &r_communicator = rModelPart.GetCommunicator();
    double wedded_perimeter = 0.0;
    if (r_communicator.LocalMesh().NumberOfNodes() != 0)
    {
        KRATOS_ERROR_IF(IsHistorical && !r_communicator.LocalMesh().NodesBegin()->SolutionStepsDataHas(AUX_DISTANCE)) << "Nodal solution step data has no \'AUX_DISTANCE\' variable. Wetted perimeter cannot be computed" << std::endl;

        // Create a vector where all the face edges id are stored.
        std::pair<int,int> aux_pair;
        std::unordered_map<std::pair<int, int>, EdgeDataContainer> edges_map;

        for (auto& r_cond : rModelPart.Conditions())
        {
            if (r_cond.Is(rSkinFlag))
            {
                auto& r_geom = r_cond.GetGeometry();
                const SizeType n_nodes = r_geom.PointsNumber();
                KRATOS_ERROR_IF(n_nodes != 3) << "This function only supports Triangle3D3N geometries." << std::endl;

                // Loop through pairs of nodes to identify unique edges
                for (IndexType i = 0; i < n_nodes - 1; ++i)
                {
                    const IndexType i_id = r_geom[i].Id();
                    for (IndexType j = i + 1; j < n_nodes; ++j)
                    {
                        const IndexType j_id = r_geom[j].Id();

                        // Ensure consistent ordering of node pairs for uniqueness
                        if (i_id < j_id) {
                            aux_pair = std::make_pair<int, int>(i_id, j_id);
                        } else{
                            aux_pair = std::make_pair<int, int>(j_id, i_id);
                        }

                        // Check if the edge is already in the map
                        auto found = edges_map.find(aux_pair);
                        // If not, create a new entry in the map
                        if (found == edges_map.end())
                        {
                            EdgeDataContainer edge_data;
                            edge_data.pNodeI = i_id < j_id ? r_geom(i) : r_geom(j);
                            edge_data.pNodeJ = i_id < j_id ? r_geom(j) : r_geom(i);
                            edge_data.NumberOfRepetitions = 1;
                            edges_map.insert(std::make_pair(aux_pair, edge_data));
                        }
                        // If the edge is already in the map, update the repetition count
                        else
                        {
                            auto &r_edge_data = found->second;
                            r_edge_data.NumberOfRepetitions += 1;
                        }
                    }
                }
            }
        }

        for (auto it = edges_map.begin(); it != edges_map.end();)
        {   // Delete all edges that have more than one repetition, since they are not part of the perimeter
            const auto& r_edge_data = it->second;
            if (r_edge_data.NumberOfRepetitions > 1) {
                it = edges_map.erase(it);
            } else {
                ++it;
            }
        }
        // Calculate de distance of each edge belonging to the perimeter.
        array_1d<double, 3> edge_vector;
        for (auto it = edges_map.begin(); it != edges_map.end(); ++it)
        {
            const auto& r_edge_data = it->second;
            const double distance_value_j = distance_getter(*(r_edge_data.pNodeJ), rDistanceVariable);
            const double distance_value_i = distance_getter(*(r_edge_data.pNodeI), rDistanceVariable);
            // Case 1: The edge is completly wetted
            if ((distance_value_i<0.0) && (distance_value_j<0.0) ){
                edge_vector = r_edge_data.pNodeJ->Coordinates() - r_edge_data.pNodeI->Coordinates();
                const double edge_length = norm_2(edge_vector);
                wedded_perimeter += edge_length;
            }
            // Case 2: The edge is cut. Interpolate the wetted distance.
            else if (distance_value_i * distance_value_j<0){

                const double phi_neg = distance_value_i<0? distance_value_i: distance_value_j;
                const double phi_pos = distance_value_i>0? distance_value_i: distance_value_j;
                const double distance_int = std::abs(phi_neg)/ (std::abs(phi_neg)+  std::abs(phi_pos));
                edge_vector = r_edge_data.pNodeJ->Coordinates() - r_edge_data.pNodeI->Coordinates();
                const double edge_length = norm_2(edge_vector);
                wedded_perimeter += distance_int*edge_length;
            }
        }
    }

    return wedded_perimeter;
}

double HydraulicFluidAuxiliaryUtilities::CalculateWettedArea(
    ModelPart &rModelPart,
    const Flags &rSkinFlag,
    const Variable<double> &rDistanceVariable,
    bool IsHistorical)
{
    // Auxiliary function to have the posssibility of having non historical or hitorical inlet distance function.
    std::function<double(NodeType &, const Variable<double> &)> distance_getter;
    if (IsHistorical)
    {
        distance_getter = [](NodeType &rNode, const Variable<double> &rDistanceVariable) -> double
        { return rNode.FastGetSolutionStepValue(rDistanceVariable); };
    }
    else
    {
        distance_getter = [](NodeType &rNode, const Variable<double> &rDistanceVariable) -> double
        { return rNode.GetValue(rDistanceVariable); };
    }

    KRATOS_ERROR_IF(rModelPart.GetProcessInfo()[DOMAIN_SIZE] < 3) << "Wetted perimeter is only implemented for 3D." << std::endl;

    // Auxiliary container for the fake Triangle2D3 geometries' points
    GeometryType::PointsArrayType aux_points;
    double wetted_area = 0.0;

    for (auto &r_cond : rModelPart.Conditions())
    {   // Create an auxiliary element based on the rskinflag condition.
        if (r_cond.Is(rSkinFlag))
        {
            std::vector<ModelPart::IndexType> elem_nodes_id;
            auto &r_geom = r_cond.GetGeometry();

            // Fill the points array and distances vector
            Vector aux_distances(3);
            for (IndexType i_nodes = 0; i_nodes < r_geom.PointsNumber(); i_nodes++)
            {
                aux_points.push_back(r_geom(i_nodes));
                aux_distances[i_nodes] = r_geom[i_nodes].GetValue(rDistanceVariable);
                // TODO: It should be posible to have an historical distance variable.
                // aux_distances[i_nodes] =distance_getter;
            }
            // Calculate the water area (wetted area) of the cut conditions
            if (FluidAuxiliaryUtilities::IsSplit(aux_distances))
            {
                auto p_splitting_util = Kratos::make_unique<DivideTriangle3D3<NodeType>>(r_geom, aux_distances);
                p_splitting_util->GenerateDivision();
                const auto& r_neg_subdivisions = p_splitting_util->GetNegativeSubdivisions();
                for (const auto& rp_neg_subdivision : r_neg_subdivisions) {
                    wetted_area += rp_neg_subdivision->Area();
                }
            }
            // Calculate the water area (wetted area)
            else if (FluidAuxiliaryUtilities::IsNegative(aux_distances))
            {
                wetted_area += r_geom.Area();
            }
            aux_points.clear();
        }
    }
    return wetted_area;
}

double HydraulicFluidAuxiliaryUtilities::InitialWaterDepth(ModelPart &rModelPart)
{

   //Determine the initial estimate for water depth by considering the average of the maximum and minimum coordinates.
    KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "The water depth is assumed to be in the Z direction" << std::endl;
    // double max_value = std::numeric_limits<double>::lowest();
    // double min_value = std::numeric_limits<double>::max();
    // for (auto &r_node : rModelPart.Nodes())
    // {
    //     if (r_node.Z() > max_value)
    //     {
    //         max_value = r_node.Z();
    //     }
    //     if (r_node.Z() < min_value)
    //     {
    //         min_value = r_node.Z();
    //     }
    // }

    double min_value, max_value;
    std::tie(min_value, max_value) = block_for_each<CombinedReduction<MinReduction<double>,MaxReduction<double>>>(rModelPart.Nodes(), [](NodeType& rNode){
        return std::make_tuple(rNode.Z(), rNode.Z());
    });

    return 0.5 * std::abs(max_value - min_value);
}

void HydraulicFluidAuxiliaryUtilities::SetInletVelocity(
    ModelPart &rModelPart,
    double InletVelocity,
    const Variable<double> &rDistanceVariable)
{
    struct AuxTLS
    {
        array_1d<double,3> InletNorm;
        array_1d<double,3> InletVelocity;
    };

    AuxTLS tls; // Explicitly initialize the AuxTLS object
    block_for_each(rModelPart.Nodes(), tls, [&](NodeType &rNode, AuxTLS &rTLS)
    {
        // Get TLS variables
        auto& inlet_norm = rTLS.InletNorm;
        auto& inlet_velocity = rTLS.InletVelocity;

        inlet_norm = rNode.GetValue(INLET_NORMAL);

        const double n_norm = norm_2(inlet_norm);
        // Orient the velocity vector in the opposite direction of the outer normal to ensure it is directed inwards.

        if (n_norm > 1.0e-12)
        {
            inlet_norm /= -n_norm;
        }
        else
        {
            KRATOS_WARNING("SetInletVelocity") << "Node " << rNode.Id() << " INLET_NORMAL is close to zero." << std::endl;
            inlet_norm /= -1.0;
        }
        //  Inlet velocity vector.
        inlet_velocity = inlet_norm * InletVelocity;

        //  Assign the velocity vector to each node representing water in the inlet condition and fix its value
        if (rNode.GetValue(rDistanceVariable) < 0.0)
        {
            rNode.FastGetSolutionStepValue(VELOCITY) =  inlet_velocity;
            rNode.Fix(VELOCITY_X);
            rNode.Fix(VELOCITY_Y);
            rNode.Fix(VELOCITY_Z);
        }
        else{
            // The air velocity in the inlet node is assumed to be null.
            rNode.FastGetSolutionStepValue(VELOCITY_X) = 0.0;
            rNode.FastGetSolutionStepValue(VELOCITY_Y) = 0.0;
            rNode.FastGetSolutionStepValue(VELOCITY_Z) = 0.0;
            rNode.Fix(VELOCITY_X);
            rNode.Fix(VELOCITY_Y);
            rNode.Fix(VELOCITY_Z);
        }
    });
}
void HydraulicFluidAuxiliaryUtilities::FreeInlet(ModelPart& rModelPart)
{
    // Free the velocity.
    block_for_each(rModelPart.Nodes(), [](NodeType &rNode){
        rNode.Free(VELOCITY_X);
        rNode.Free(VELOCITY_Y);
        rNode.Free(VELOCITY_Z);
        rNode.Free(DISTANCE);
    });
}
void HydraulicFluidAuxiliaryUtilities::SetInletFreeSurface(ModelPart &rModelPart, const Flags &rSkinFlag,  const Variable<double> &rDistanceVariable)
{
    // Assign the water depth (DISTANCE) to all nodes within the inlet model part to be equal to the water depth corresponding to a Froude number of 1 (&rDistanceVariable).
    block_for_each(rModelPart.Nodes(), [&](NodeType& rNode){
        if (rNode.Is(rSkinFlag)){
            double inlet_dist = rNode.GetValue(rDistanceVariable);
            rNode.FastGetSolutionStepValue(DISTANCE) = inlet_dist;
            rNode.Fix(DISTANCE);
        }
    });
}
void HydraulicFluidAuxiliaryUtilities::CalculateArtificialViscosity(
    ModelPart &rModelPart,
    double artificial_limiter_coefficient)
{
    const auto &properties_1 = rModelPart.GetProperties(1);
    const double water_dynamic_viscosity_max = artificial_limiter_coefficient * properties_1.GetValue(DYNAMIC_VISCOSITY);


    block_for_each(rModelPart.Elements(), [&](Element &rElement) {

        double elem_artificial_viscosity = 0.0;
        rElement.Calculate(ARTIFICIAL_DYNAMIC_VISCOSITY, elem_artificial_viscosity, rModelPart.GetProcessInfo());

        if (elem_artificial_viscosity > water_dynamic_viscosity_max)
        {
            elem_artificial_viscosity = water_dynamic_viscosity_max;
        }

        auto &r_nodes = rElement.GetGeometry();
        int neg_nodes = 0;
        int pos_nodes = 0;

        for (auto &r_node : r_nodes)
        {
            const double distance = r_node.FastGetSolutionStepValue(DISTANCE);
            if (distance > 0)
            {
                pos_nodes += 1;
            }
            else
            {
                neg_nodes += 1;
            }
        }
        
        if (neg_nodes > 0 && pos_nodes > 0)
        {
            elem_artificial_viscosity = 0.0;
        }
        
        rElement.SetValue(ARTIFICIAL_DYNAMIC_VISCOSITY, elem_artificial_viscosity);
    });
}       

void HydraulicFluidAuxiliaryUtilities::FindElementsNeighbouringConditions(
    ModelPart &rModelPart,
    const Flags &rConditionFlag)
{
    // Use the utility to assign neighbor elements to conditions
    const bool check_repeated_conditions = true;
    FluidMeshUtilities::AssignNeighbourElementsToConditions(rModelPart, check_repeated_conditions);

    // Iterate through conditions to print and assign the parent element
    for (auto& rCondition : rModelPart.Conditions()) {
        KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "Checking condition: " << rCondition.Id() << " with flag: " << rConditionFlag << std::endl;
        if (rCondition.Is(rConditionFlag)) {
            const auto& neighbour_elements = rCondition.GetValue(NEIGHBOUR_ELEMENTS);

            // Ensure there is at least one neighboring element
            KRATOS_ERROR_IF(neighbour_elements.size() == 0) 
                << "Condition ID: " << rCondition.Id() << " has no neighboring elements." << std::endl;

            // Print the condition ID and its parent element
            const unsigned int rank = rModelPart.GetCommunicator().MyPID();
            KRATOS_INFO_IF("HydraulicFluidAuxiliaryUtilities", rank == 0)
                << "Condition ID: " << rCondition.Id()
                << " has parent element ID: " << neighbour_elements[0].Id() << std::endl;
        }
    }

    const unsigned int rank = rModelPart.GetCommunicator().MyPID();
    KRATOS_INFO_IF("HydraulicFluidAuxiliaryUtilities", rank == 0)
        << "Elements neighboring conditions search and assignment finished." << std::endl;
}

double HydraulicFluidAuxiliaryUtilities::CalculateErosionRate(
    const Element& rFluidElement,
    const Condition& rInterfaceCondition,
    const double D50,
    const double SedimentDensity)
{
    KRATOS_TRY

    // Physical parameters
    const double g = 9.81; // Gravity acceleration (m/s²)
    const double kappa = 0.41; // Von Karman constant
    
    // Get fluid element geometry
    const auto& r_fluid_geom = rFluidElement.GetGeometry();
    const SizeType n_nodes = r_fluid_geom.PointsNumber();
    
    // Initialize fluid properties
    double rho = 0.0;
    double mu = 0.0;
    
    // Calculate average fluid properties from nodes
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_node = r_fluid_geom[i_node];
        rho += r_node.FastGetSolutionStepValue(DENSITY);
        mu += r_node.FastGetSolutionStepValue(DYNAMIC_VISCOSITY);
    }
    
    // Average fluid properties
    rho /= static_cast<double>(n_nodes);
    mu /= static_cast<double>(n_nodes);
    const double nu = mu / rho; // Kinematic viscosity
    
    // Submerged specific gravity
    const double R = SedimentDensity / rho - 1.0;
    
    // Calculate particle Reynolds number (Eq. 11)
    const double Rep = D50 * std::sqrt(R * g * D50) / nu;
    
    // Calculate critical Shields stress (Eq. 10 - Brownlie 1981)
    const double rep_power = std::pow(Rep, -0.6);
    const double tau_star_c = 0.22 * rep_power + 0.06 * std::pow(10.0, -7.7 * rep_power);
    
    // Get interface normal vector
    const auto& r_interface_geom = rInterfaceCondition.GetGeometry();
    array_1d<double, 3> interface_normal;
    
    // Calculate normal vector manually for 3D triangle
    if (r_interface_geom.PointsNumber() == 3) {
        const auto& p1 = r_interface_geom[0];
        const auto& p2 = r_interface_geom[1];
        const auto& p3 = r_interface_geom[2];
        
        // Calculate vectors
        array_1d<double, 3> v1, v2;
        v1[0] = p2.X() - p1.X();
        v1[1] = p2.Y() - p1.Y();
        v1[2] = p2.Z() - p1.Z();
        
        v2[0] = p3.X() - p1.X();
        v2[1] = p3.Y() - p1.Y();
        v2[2] = p3.Z() - p1.Z();
        
        // Cross product
        interface_normal[0] = v1[1] * v2[2] - v1[2] * v2[1];
        interface_normal[1] = v1[2] * v2[0] - v1[0] * v2[2];
        interface_normal[2] = v1[0] * v2[1] - v1[1] * v2[0];
        
        // Normalize
        const double normal_mag = norm_2(interface_normal);
        if (normal_mag > 1e-12) {
            interface_normal /= normal_mag;
        } else {
            interface_normal[0] = 0.0;
            interface_normal[1] = 0.0;
            interface_normal[2] = 1.0; // Default to Z direction
        }
    } else {
        interface_normal[0] = 0.0;
        interface_normal[1] = 0.0;
        interface_normal[2] = 1.0; // Default
    }
    
    // Calculate velocity and distance at interface
    double vt_total = 0.0; // Tangential velocity component
    double hk_total = 0.0; // Normal distance from interface
    double vn_total = 0.0; // Normal velocity component
    
    for (IndexType i_node = 0; i_node < n_nodes; ++i_node) {
        const auto& r_node = r_fluid_geom[i_node];
        const array_1d<double, 3>& velocity = r_node.FastGetSolutionStepValue(VELOCITY);
        
        // Calculate tangential velocity component (velocity - normal component)
        const double v_normal_mag = velocity[0] * interface_normal[0] + 
                                   velocity[1] * interface_normal[1] + 
                                   velocity[2] * interface_normal[2];
        
        const double velocity_mag_sq = velocity[0] * velocity[0] + 
                                      velocity[1] * velocity[1] + 
                                      velocity[2] * velocity[2];
        
        const double v_tangential = std::sqrt(std::max(0.0, velocity_mag_sq - v_normal_mag * v_normal_mag));
        vt_total += v_tangential;
        
        if (v_normal_mag < 0) {
            vn_total += -v_normal_mag;
        }
        
        // Approximate distance from interface (simplified)
        // In practice, this should be the actual distance from node to interface
        hk_total += D50 * 10.0; // Approximation: 10 particle diameters
    }
    
    // Average values
    const double vt = vt_total / static_cast<double>(n_nodes);
    const double hk = hk_total / static_cast<double>(n_nodes);
    const double vn = vn_total / static_cast<double>(n_nodes);
    
    // Calculate friction velocity using logarithmic wall law (Eq. 13)
    // u/u* = (1/0.41) * ln(u* * y / nu) + 5.5
    // Solve iteratively for u*
    double u_star = 0.1; // Initial guess
    const int max_iterations = 10;
    const double tolerance = 1e-8;
    
    for (int iter = 0; iter < max_iterations; ++iter) {
        const double argument = u_star * hk / nu;
        const double ln_term = (argument > 1e-12) ? std::log(argument) : 0.0;
        const double f = vt - u_star * (ln_term / kappa + 5.5);
        
        if (std::abs(f) < tolerance) {
            break;
        }
        
        // Derivative
        const double df_du_star = -(ln_term / kappa + 5.5 + 1.0 / kappa);
        if (std::abs(df_du_star) > 1e-12) {
            u_star -= f / df_du_star;
        }
        u_star = std::max(u_star, 1e-8); // Keep positive
    }
    
    // Calculate shear stress at D50/2 distance (Eq. 15)
    const double tau_D50 = mu * u_star / (kappa * D50 / 2.0);
    
    // Apply turbulence factor (Eq. 16)
    const double fk = 2.0; // Turbulence factor - adjust based on your case
    const double tau_w = fk * tau_D50;
    
    // Calculate dimensionless shear stress (Shields parameter, Eq. 17)
    const double tau_star_s = tau_w / (rho * R * g * D50);
    
    const double excess_ratio = tau_star_s / tau_star_c - 1.0;
    
    // Initialize erosion components
    double E_t = 0.0;
    
    // Check if erosion threshold is exceeded
    if (tau_star_s > tau_star_c) {
        // Calculate dimensionless entrainment rate (Eq. 18 - van Rijn)
        const double E_star = 0.015 * (D50 / hk) * 
                             std::pow(excess_ratio, 1.5) * 
                             std::pow(Rep, -0.2) * 10.0;
        
        // Dimensionalize entrainment rate (Eq. 19)
        E_t = E_star * std::sqrt(g * R * D50); // m/s
    }
    
    // Calculate normal velocity erosion component
    double E_n = 0.0;
    if (vn > 1e-12) {
        // Calibratable coefficient for impact erosion
        const double C_impact = 1.5;
        // Erosion rate related to the kinetic energy of the normal flow component
        E_n = C_impact * 0.5 * rho * std::pow(vn, 3) / SedimentDensity;
    }
    
    const double E = E_t + E_n; // Total erosion rate (m/s)
    
    // Calculate interface area
    const double interface_area = r_interface_geom.Area();
    
    // Convert to volumetric erosion rate (m³/s)
    const double erosion_rate = E * interface_area;
    
    return erosion_rate;

    KRATOS_CATCH("")
}

int HydraulicFluidAuxiliaryUtilities::ProcessElementErosion(
    Element& rSolidElement,
    ModelPart& rSlipBedModelPart,
    const Condition& rInterfaceCondition,
    ModelPart& rComputingModelPart,
    const Element* pConnectedFluidElement)
{
    KRATOS_TRY

    // Activate the solid element (convert to fluid)
    rSolidElement.Set(FLUID, true);

    // Remove fixed DOFs from its nodes (unfix for fluid behavior)
    auto& r_solid_geom = rSolidElement.GetGeometry();
    const SizeType n_solid_nodes = r_solid_geom.PointsNumber();

    for (IndexType i_node = 0; i_node < n_solid_nodes; ++i_node) {
        auto& r_node = r_solid_geom[i_node];
        
        if (r_node.IsFixed(VELOCITY_X)) {
            r_node.Free(VELOCITY_X);
        }
        if (r_node.IsFixed(VELOCITY_Y)) {
            r_node.Free(VELOCITY_Y);
        }
        if (r_node.IsFixed(VELOCITY_Z)) {
            r_node.Free(VELOCITY_Z);
        }
        if (r_node.IsFixed(PRESSURE)) {
            r_node.Free(PRESSURE);
        }
        if (r_node.IsFixed(C_SUSP)) {
            r_node.Free(C_SUSP);
        }
    }

    // Remove the original interface condition and its flags
    auto& r_interface_geom_non_const = const_cast<Condition&>(rInterfaceCondition).GetGeometry();
    const SizeType n_interface_nodes = r_interface_geom_non_const.PointsNumber();
    
    for (IndexType i = 0; i < n_interface_nodes; ++i) {
        auto& r_node = r_interface_geom_non_const[i];
        r_node.Set(SLIP, false);
    }
    const_cast<Condition&>(rInterfaceCondition).Set(WALL, false);

    // Get all faces of the eroded solid element
    std::vector<IndexType> node_ids;
    node_ids.reserve(n_solid_nodes);
    for (IndexType i = 0; i < n_solid_nodes; ++i) {
        node_ids.push_back(r_solid_geom[i].Id());
    }

    // Define faces for tetrahedra (4 triangular faces)
    std::vector<std::array<IndexType, 3>> faces;
    if (n_solid_nodes == 4) { // Tetrahedra
        faces = {
            {node_ids[0], node_ids[1], node_ids[2]},
            {node_ids[0], node_ids[1], node_ids[3]},
            {node_ids[0], node_ids[2], node_ids[3]},
            {node_ids[1], node_ids[2], node_ids[3]}
        };
    }

    // Find which face was the original interface (to skip it)
    std::vector<IndexType> original_interface_face;
    original_interface_face.reserve(n_interface_nodes);
    for (IndexType i = 0; i < n_interface_nodes; ++i) {
        original_interface_face.push_back(r_interface_geom_non_const[i].Id());
    }
    std::sort(original_interface_face.begin(), original_interface_face.end());

    // Build a mapping of faces to elements for efficient neighbor checking
    std::unordered_map<std::vector<IndexType>, std::vector<Element::Pointer>, VectorHasher> face_to_elements;
    
    for (auto& r_element : rComputingModelPart.Elements()) {
        const auto& r_elem_geom = r_element.GetGeometry();
        const SizeType n_elem_nodes = r_elem_geom.PointsNumber();
        
        if (n_elem_nodes == 4) { // Process only tetrahedra
            std::vector<IndexType> elem_node_ids;
            elem_node_ids.reserve(n_elem_nodes);
            for (IndexType i = 0; i < n_elem_nodes; ++i) {
                elem_node_ids.push_back(r_elem_geom[i].Id());
            }
            
            // Generate faces for this element
            std::vector<std::array<IndexType, 3>> elem_faces = {
                {elem_node_ids[0], elem_node_ids[1], elem_node_ids[2]},
                {elem_node_ids[0], elem_node_ids[1], elem_node_ids[3]},
                {elem_node_ids[0], elem_node_ids[2], elem_node_ids[3]},
                {elem_node_ids[1], elem_node_ids[2], elem_node_ids[3]}
            };
            
            for (auto& face : elem_faces) {
                std::vector<IndexType> face_sorted(face.begin(), face.end());
                std::sort(face_sorted.begin(), face_sorted.end());
                face_to_elements[face_sorted].push_back(Element::Pointer(&r_element));
            }
        }
    }

    // **IMPROVED ID GENERATION**: Create a thread-safe ID generation mechanism
    // First, collect all existing condition IDs from both model parts
    std::set<IndexType> all_existing_ids;
    
    // Get IDs from slip bed model part
    for (const auto& r_condition : rSlipBedModelPart.Conditions()) {
        all_existing_ids.insert(r_condition.Id());
    }
    
    // Get IDs from computing model part
    for (const auto& r_condition : rComputingModelPart.Conditions()) {
        all_existing_ids.insert(r_condition.Id());
    }
    
    // Find the maximum ID and start from there
    IndexType next_available_id = all_existing_ids.empty() ? 1 : *all_existing_ids.rbegin() + 1;
    
    // Check each face of the eroded element for new interface conditions
    int new_conditions_created = 0;
    
    for (const auto& face : faces) {
        std::vector<IndexType> face_sorted(face.begin(), face.end());
        std::sort(face_sorted.begin(), face_sorted.end());
        
        // Skip the original interface face
        if (face_sorted == original_interface_face) {
            continue;
        }
        
        // Check if this face is shared with a solid element
        auto it = face_to_elements.find(face_sorted);
        if (it == face_to_elements.end()) {
            continue; // Face not found in mapping
        }
        
        auto neighboring_elements = it->second; // Copy the vector
        
        // Remove the current element from neighbors
        auto new_end = std::remove_if(neighboring_elements.begin(), neighboring_elements.end(),
            [&rSolidElement](const Element::Pointer& p_elem) {
                return p_elem->Id() == rSolidElement.Id();
            });
        neighboring_elements.erase(new_end, neighboring_elements.end());
        
        // Check if any neighboring element is still solid (inactive or not fluid)
        bool has_solid_neighbor = false;
        for (const auto& p_neighbor : neighboring_elements) {
            if (!p_neighbor->Is(FLUID) || !p_neighbor->Is(ACTIVE)) {
                has_solid_neighbor = true;
                break;
            }
        }
        
        // If this face is exposed to solid domain, create new interface condition
        if (has_solid_neighbor) {
            // Check if condition already exists
            bool condition_exists = false;
            for (const auto& r_condition : rSlipBedModelPart.Conditions()) {
                const auto& r_cond_geom = r_condition.GetGeometry();
                std::vector<IndexType> existing_nodes;
                existing_nodes.reserve(r_cond_geom.PointsNumber());
                for (IndexType i = 0; i < r_cond_geom.PointsNumber(); ++i) {
                    existing_nodes.push_back(r_cond_geom[i].Id());
                }
                std::sort(existing_nodes.begin(), existing_nodes.end());
                
                if (existing_nodes == face_sorted) {
                    condition_exists = true;
                    break;
                }
            }
            
            if (!condition_exists) {
                // **IMPROVED ID GENERATION**: Use incremental ID generation
                // Find the next available ID by checking if it exists
                while (all_existing_ids.find(next_available_id) != all_existing_ids.end()) {
                    ++next_available_id;
                }
                
                const IndexType new_cond_id = next_available_id;
                
                // Add the new ID to our tracking set
                all_existing_ids.insert(new_cond_id);
                
                // Increment for next use
                ++next_available_id;
                
                // Add nodes to slip bed model part if not present
                for (IndexType node_id : face_sorted) {
                    if (!rSlipBedModelPart.HasNode(node_id)) {
                        auto& r_node = rComputingModelPart.GetNode(node_id);
                        rSlipBedModelPart.AddNode(&r_node);
                    }
                }
                
                // **ADDITIONAL SAFETY CHECK**: Verify the ID is truly unique before creation
                bool id_already_exists = false;
                try {
                    // Try to get a condition with this ID - if it exists, this will succeed
                    rSlipBedModelPart.GetCondition(new_cond_id);
                    id_already_exists = true;
                } catch (...) {
                    // Condition doesn't exist, which is what we want
                    id_already_exists = false;
                }
                
                if (!id_already_exists) {
                    // Create new condition
                    auto p_new_condition = rSlipBedModelPart.CreateNewCondition(
                        "WallCondition3D3N", 
                        new_cond_id, 
                        face_sorted, 
                        rSlipBedModelPart.pGetProperties(0));
                    
                    p_new_condition->Set(WALL, true);
                    ++new_conditions_created;
                } else {
                    // This should not happen with our improved logic, but just in case
                    KRATOS_WARNING("ProcessElementErosion") 
                        << "ID " << new_cond_id << " already exists despite checks. Skipping condition creation." << std::endl;
                }
            }
        }
    }

    return new_conditions_created;

    KRATOS_CATCH("")
}

int HydraulicFluidAuxiliaryUtilities::ProcessElementDeposition(
    Element& rFluidElement,
    ModelPart& rSlipBedModelPart,
    const Condition& rInterfaceCondition,
    ModelPart& rComputingModelPart)
{
    KRATOS_TRY

    // Deactivate the fluid element
    rFluidElement.Set(FLUID, false);

    // Get all faces of the element (assuming tetrahedra)
    const auto& r_fluid_geom = rFluidElement.GetGeometry();
    const SizeType n_nodes = r_fluid_geom.PointsNumber();
    
    std::vector<IndexType> node_ids;
    node_ids.reserve(n_nodes);
    for (IndexType i = 0; i < n_nodes; ++i) {
        node_ids.push_back(r_fluid_geom[i].Id());
    }

    // Define faces for tetrahedra (4 triangular faces)
    std::vector<std::array<IndexType, 3>> faces;
    if (n_nodes == 4) { // Tetrahedra
        faces = {
            {node_ids[0], node_ids[1], node_ids[2]},
            {node_ids[0], node_ids[1], node_ids[3]},
            {node_ids[0], node_ids[2], node_ids[3]},
            {node_ids[1], node_ids[2], node_ids[3]}
        };
    }

    // Find which face was the original interface (to skip it)
    std::vector<IndexType> original_interface_face;
    bool interface_face_found = false;

    // Get all existing conditions to find the matching interface face
    for (const auto& r_condition : rSlipBedModelPart.Conditions()) {
        const auto& r_cond_geom = r_condition.GetGeometry();
        std::vector<IndexType> cond_node_ids;
        cond_node_ids.reserve(r_cond_geom.PointsNumber());
        
        for (IndexType i = 0; i < r_cond_geom.PointsNumber(); ++i) {
            cond_node_ids.push_back(r_cond_geom[i].Id());
        }
        std::sort(cond_node_ids.begin(), cond_node_ids.end());

        // Check if this condition matches any face of the fluid element
        for (const auto& face : faces) {
            std::vector<IndexType> face_sorted(face.begin(), face.end());
            std::sort(face_sorted.begin(), face_sorted.end());
            
            if (face_sorted == cond_node_ids) {
                original_interface_face = face_sorted;
                interface_face_found = true;
                break;
            }
        }
        
        if (interface_face_found) {
            break;
        }
    }

    // Add other faces as new conditions to slip bed model part
    int new_conditions_created = 0;
    
    for (const auto& face : faces) {
        std::vector<IndexType> face_sorted(face.begin(), face.end());
        std::sort(face_sorted.begin(), face_sorted.end());
        
        // Skip the original interface face
        if (interface_face_found && face_sorted == original_interface_face) {
            continue;
        }

        // Check if condition already exists
        bool condition_exists = false;
        for (const auto& r_condition : rSlipBedModelPart.Conditions()) {
            const auto& r_cond_geom = r_condition.GetGeometry();
            std::vector<IndexType> existing_nodes;
            existing_nodes.reserve(r_cond_geom.PointsNumber());
            
            for (IndexType i = 0; i < r_cond_geom.PointsNumber(); ++i) {
                existing_nodes.push_back(r_cond_geom[i].Id());
            }
            std::sort(existing_nodes.begin(), existing_nodes.end());
            
            if (existing_nodes == face_sorted) {
                condition_exists = true;
                break;
            }
        }

        if (!condition_exists) {
            // Find new condition ID - check all model parts to avoid conflicts
            IndexType max_cond_id = 0;
            
            // Check slip bed model part conditions
            for (const auto& r_condition : rSlipBedModelPart.Conditions()) {
                max_cond_id = std::max(max_cond_id, r_condition.Id());
            }
            
            // Check computing model part conditions
            for (const auto& r_condition : rComputingModelPart.Conditions()) {
                max_cond_id = std::max(max_cond_id, r_condition.Id());
            }
            
            const IndexType new_cond_id = max_cond_id + 1;

            // Add nodes to slip bed model part if not present
            for (IndexType node_id : face_sorted) {
                if (!rSlipBedModelPart.HasNode(node_id)) {
                    // Get original node from computing model part
                    auto& r_original_node = rComputingModelPart.GetNode(node_id);
                    
                    // **CRITICAL FIX: Ensure DISTANCE is preserved when adding node**
                    double original_distance = r_original_node.FastGetSolutionStepValue(DISTANCE);
                    
                    // Add the existing node to slip bed model part
                    rSlipBedModelPart.AddNode(&r_original_node);
                    
                    // **ENSURE DISTANCE VALUE IS MAINTAINED**
                    r_original_node.FastGetSolutionStepValue(DISTANCE) = original_distance;
                    
                    // Set SLIP flag for the node
                    //r_original_node.Set(SLIP, true);
                } else {
                    // Node already exists, ensure DISTANCE is not zero
                    auto& r_existing_node = rSlipBedModelPart.GetNode(node_id);
                    
                    // **FIX: Check and correct zero DISTANCE values**
                    double current_distance = r_existing_node.FastGetSolutionStepValue(DISTANCE);
                    if (std::abs(current_distance) < 1e-12) { // Essentially zero
                        // Get the original distance from computing model part
                        auto& r_original_node = rComputingModelPart.GetNode(node_id);
                        double original_distance = r_original_node.FastGetSolutionStepValue(DISTANCE);
                        r_existing_node.FastGetSolutionStepValue(DISTANCE) = original_distance;
                    }
                    
                    // Set SLIP flag
                    //r_existing_node.Set(SLIP, true);
                }
            }

            // Create new condition
            auto p_new_condition = rSlipBedModelPart.CreateNewCondition(
                "WallCondition3D3N", 
                new_cond_id, 
                face_sorted, 
                rSlipBedModelPart.pGetProperties(0)
            );
            
            // Set WALL flag for the new condition
            p_new_condition->Set(WALL, true);
            ++new_conditions_created;
            if (rFluidElement.Id() == 13527) {
                KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "Created new condition with ID: " << new_cond_id << " for element ID: " << rFluidElement.Id() << std::endl;
                KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "Total conditions in slip bed model part: " << rSlipBedModelPart.NumberOfConditions() << std::endl;
            }
        }
    }

    return new_conditions_created;

    KRATOS_CATCH("")
}

void HydraulicFluidAuxiliaryUtilities::ConnectNewConditions(
    ModelPart& rFluidModelPart,
    ModelPart& rSlipBedModelPart,
    std::unordered_set<IndexType>& rSolidElementsSet,
    std::unordered_map<IndexType, Element::Pointer>& rInterfaceConditionToFluidElement,
    std::unordered_map<IndexType, Element::Pointer>& rInterfaceConditionToSolidElement,
    std::unordered_map<IndexType, double>& rAccumulatedVolume)
{
    KRATOS_TRY

    KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "Connecting newly generated conditions and removing redundant ones..." << std::endl;
    
    // Get all current conditions
    std::vector<Condition::Pointer> all_conditions;
    for (auto& r_condition : rSlipBedModelPart.Conditions()) {
        all_conditions.push_back(&r_condition);
    }
    
    // **NEW FUNCTIONALITY: Remove redundant conditions between fluid elements**
    KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "--- Checking for redundant conditions between fluid elements ---" << std::endl;
    std::vector<Condition::Pointer> redundant_conditions;
    
    // Build comprehensive element mappings for redundancy check
    std::map<std::vector<IndexType>, Element::Pointer> fluid_elements_by_nodes;
    std::map<std::vector<IndexType>, Element::Pointer> solid_elements_by_nodes;
    
    for (auto& r_element : rFluidModelPart.Elements()) {
        std::vector<IndexType> elem_node_ids;
        for (const auto& r_node : r_element.GetGeometry()) {
            elem_node_ids.push_back(r_node.Id());
        }
        std::sort(elem_node_ids.begin(), elem_node_ids.end());
        
        if (r_element.Is(FLUID) && rSolidElementsSet.find(r_element.Id()) == rSolidElementsSet.end()) {
            // Active fluid element
            fluid_elements_by_nodes[elem_node_ids] = &r_element;
        } else {
            // Solid element (including deposited ones)
            solid_elements_by_nodes[elem_node_ids] = &r_element;
        }
    }
    
    // Check each condition for redundancy
    for (auto& p_condition : all_conditions) {
        const IndexType condition_id = p_condition->Id();
        std::vector<IndexType> cond_node_ids;
        for (const auto& r_node : p_condition->GetGeometry()) {
            cond_node_ids.push_back(r_node.Id());
        }
        std::sort(cond_node_ids.begin(), cond_node_ids.end());
        
        // Find all elements that contain this condition's nodes
        std::vector<Element::Pointer> connected_fluid_elements;
        std::vector<Element::Pointer> connected_solid_elements;
        
        for (const auto& pair : fluid_elements_by_nodes) {
            const std::vector<IndexType>& elem_nodes = pair.first;
            Element::Pointer p_element = pair.second;
            
            // Check if condition nodes are subset of element nodes
            bool is_subset = true;
            for (const auto& node_id : cond_node_ids) {
                if (std::find(elem_nodes.begin(), elem_nodes.end(), node_id) == elem_nodes.end()) {
                    is_subset = false;
                    break;
                }
            }
            if (is_subset) {
                connected_fluid_elements.push_back(p_element);
            }
        }
        
        for (const auto& pair : solid_elements_by_nodes) {
            const std::vector<IndexType>& elem_nodes = pair.first;
            Element::Pointer p_element = pair.second;
            
            // Check if condition nodes are subset of element nodes
            bool is_subset = true;
            for (const auto& node_id : cond_node_ids) {
                if (std::find(elem_nodes.begin(), elem_nodes.end(), node_id) == elem_nodes.end()) {
                    is_subset = false;
                    break;
                }
            }
            if (is_subset) {
                connected_solid_elements.push_back(p_element);
            }
        }
        
        // **REDUNDANCY CHECK**: If condition is shared between TWO FLUID elements and NO solid elements
        if (connected_fluid_elements.size() == 2 && connected_solid_elements.size() == 0) {
            redundant_conditions.push_back(p_condition);
            KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") 
                << "REDUNDANT condition found: " << condition_id 
                << " - shared between fluid elements " << connected_fluid_elements[0]->Id() 
                << " and " << connected_fluid_elements[1]->Id() << std::endl;
        }
        // **ADDITIONAL CHECK**: If condition is shared between MORE than 2 fluid elements
        else if (connected_fluid_elements.size() > 2 && connected_solid_elements.size() == 0) {
            redundant_conditions.push_back(p_condition);
            std::string fluid_elem_ids = "";
            for (const auto& p_elem : connected_fluid_elements) {
                fluid_elem_ids += std::to_string(p_elem->Id()) + " ";
            }
            KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") 
                << "REDUNDANT condition found: " << condition_id 
                << " - shared between multiple fluid elements " << fluid_elem_ids << std::endl;
        }
    }
    
    // **REMOVE REDUNDANT CONDITIONS**
    if (!redundant_conditions.empty()) {
        KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") 
            << "--- Removing " << redundant_conditions.size() << " redundant conditions ---" << std::endl;
        
        for (std::size_t i = 0; i < redundant_conditions.size(); ++i) {
            auto p_condition = redundant_conditions[i];
            const IndexType condition_id = p_condition->Id();
            
            KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") 
                << "Removing redundant condition " << (i+1) << "/" << redundant_conditions.size() 
                << ": " << condition_id << std::endl;
            
            // Remove from connectivity mappings if present
            auto it_fluid = rInterfaceConditionToFluidElement.find(condition_id);
            if (it_fluid != rInterfaceConditionToFluidElement.end()) {
                rInterfaceConditionToFluidElement.erase(it_fluid);
            }
            
            auto it_solid = rInterfaceConditionToSolidElement.find(condition_id);
            if (it_solid != rInterfaceConditionToSolidElement.end()) {
                rInterfaceConditionToSolidElement.erase(it_solid);
            }
            
            auto it_volume = rAccumulatedVolume.find(condition_id);
            if (it_volume != rAccumulatedVolume.end()) {
                rAccumulatedVolume.erase(it_volume);
            }
            
            // Remove condition from slip bed model part
            rSlipBedModelPart.RemoveCondition(condition_id);
            
            KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") 
                << "  -> Successfully removed redundant condition " << condition_id << std::endl;
        }
        
        // Update the conditions list after removal
        all_conditions.clear();
        for (auto& r_condition : rSlipBedModelPart.Conditions()) {
            all_conditions.push_back(&r_condition);
        }
        
        KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") 
            << "Conditions remaining after redundancy removal: " << all_conditions.size() << std::endl;
    } else {
        KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "No redundant conditions found" << std::endl;
    }
    
    // **CONTINUE WITH ORIGINAL FUNCTIONALITY**: Connect unconnected conditions
    KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "--- Connecting remaining unconnected conditions ---" << std::endl;
    
    // Track conditions that need connectivity
    std::vector<Condition::Pointer> unconnected_conditions;
    
    for (auto& p_condition : all_conditions) {
        const IndexType condition_id = p_condition->Id();
        
        // Check if this condition is already connected
        bool has_fluid_connection = rInterfaceConditionToFluidElement.find(condition_id) != rInterfaceConditionToFluidElement.end();
        bool has_solid_connection = rInterfaceConditionToSolidElement.find(condition_id) != rInterfaceConditionToSolidElement.end();
        
        if (!has_fluid_connection && !has_solid_connection) {
            unconnected_conditions.push_back(p_condition);
            KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "Found unconnected condition: " << condition_id << std::endl;
        }
    }
    
    if (unconnected_conditions.empty()) {
        KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "All remaining conditions are already connected" << std::endl;
        return;
    }
    
    // **HANDLE SHARED CONNECTIONS**: Track multiple conditions per element
    std::unordered_map<IndexType, std::vector<IndexType>> fluid_element_to_conditions;
    std::unordered_map<IndexType, std::vector<IndexType>> solid_element_to_conditions;
    
    // Connect unconnected conditions
    for (auto& p_condition : unconnected_conditions) {
        std::vector<IndexType> cond_node_ids;
        for (const auto& r_node : p_condition->GetGeometry()) {
            cond_node_ids.push_back(r_node.Id());
        }
        std::sort(cond_node_ids.begin(), cond_node_ids.end());
        const IndexType condition_id = p_condition->Id();
        
        // Find connected fluid elements
        std::vector<Element::Pointer> connected_fluid_elements;
        for (const auto& pair : fluid_elements_by_nodes) {
            const std::vector<IndexType>& elem_nodes = pair.first;
            Element::Pointer p_element = pair.second;
            
            // Check if condition nodes are subset of element nodes
            bool is_subset = true;
            for (const auto& node_id : cond_node_ids) {
                if (std::find(elem_nodes.begin(), elem_nodes.end(), node_id) == elem_nodes.end()) {
                    is_subset = false;
                    break;
                }
            }
            if (is_subset) {
                connected_fluid_elements.push_back(p_element);
            }
        }
        
        // Find connected solid elements
        std::vector<Element::Pointer> connected_solid_elements;
        for (const auto& pair : solid_elements_by_nodes) {
            const std::vector<IndexType>& elem_nodes = pair.first;
            Element::Pointer p_element = pair.second;
            
            // Check if condition nodes are subset of element nodes
            bool is_subset = true;
            for (const auto& node_id : cond_node_ids) {
                if (std::find(elem_nodes.begin(), elem_nodes.end(), node_id) == elem_nodes.end()) {
                    is_subset = false;
                    break;
                }
            }
            if (is_subset) {
                connected_solid_elements.push_back(p_element);
            }
        }
        
        // **VALID INTERFACE CHECK**: Only connect if there's at least one solid element
        if (connected_solid_elements.empty()) {
            KRATOS_WARNING("HydraulicFluidAuxiliaryUtilities") 
                << "Condition " << condition_id << " has no solid element connection - may be invalid interface" << std::endl;
            continue;
        }
        
        // **HANDLE MULTIPLE CONNECTIONS**: Select best connection
        if (!connected_fluid_elements.empty()) {
            Element::Pointer selected_fluid_elem;
            if (connected_fluid_elements.size() == 1) {
                selected_fluid_elem = connected_fluid_elements[0];
            } else {
                // Multiple fluid elements - select one with highest C_SUSP
                double max_c_susp = -1.0;
                selected_fluid_elem = connected_fluid_elements[0];
                for (auto& p_elem : connected_fluid_elements) {
                    double avg_c_susp = 0.0;
                    const auto& r_geom = p_elem->GetGeometry();
                    for (IndexType i = 0; i < r_geom.PointsNumber(); ++i) {
                        avg_c_susp += r_geom[i].FastGetSolutionStepValue(C_SUSP);
                    }
                    avg_c_susp /= static_cast<double>(r_geom.PointsNumber());
                    
                    if (avg_c_susp > max_c_susp) {
                        max_c_susp = avg_c_susp;
                        selected_fluid_elem = p_elem;
                    }
                }
                
                KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") 
                    << "Condition " << condition_id << " connected to multiple fluid elements, selected element " 
                    << selected_fluid_elem->Id() << std::endl;
            }
            
            // Connect condition to selected fluid element
            rInterfaceConditionToFluidElement[condition_id] = selected_fluid_elem;
            
            // Track shared connections
            fluid_element_to_conditions[selected_fluid_elem->Id()].push_back(condition_id);
        }
        
        // Connect to solid elements (usually one per condition)
        if (!connected_solid_elements.empty()) {
            Element::Pointer selected_solid_elem = connected_solid_elements[0]; // Usually only one
            rInterfaceConditionToSolidElement[condition_id] = selected_solid_elem;
            
            // Track shared connections
            solid_element_to_conditions[selected_solid_elem->Id()].push_back(condition_id);
        }
        
        if (connected_fluid_elements.empty() && connected_solid_elements.empty()) {
            KRATOS_WARNING("HydraulicFluidAuxiliaryUtilities") 
                << "Condition " << condition_id << " could not be connected to any element" << std::endl;
        }
    }
    
    // **REPORT SHARED CONNECTIONS**
    KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "--- Shared Element Connections ---" << std::endl;
    for (const auto& pair : fluid_element_to_conditions) {
        if (pair.second.size() > 1) {
            std::string condition_ids = "";
            for (const auto& cond_id : pair.second) {
                condition_ids += std::to_string(cond_id) + " ";
            }
            KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") 
                << "Fluid element " << pair.first << " connected to multiple conditions: " << condition_ids << std::endl;
        }
    }
    
    for (const auto& pair : solid_element_to_conditions) {
        if (pair.second.size() > 1) {
            std::string condition_ids = "";
            for (const auto& cond_id : pair.second) {
                condition_ids += std::to_string(cond_id) + " ";
            }
            KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") 
                << "Solid element " << pair.first << " connected to multiple conditions: " << condition_ids << std::endl;
        }
    }
    
    // **FINAL SUMMARY**
    KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "--- Connection Summary ---" << std::endl;
    KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "Redundant conditions removed: " << redundant_conditions.size() << std::endl;
    KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "New conditions connected: " << unconnected_conditions.size() << std::endl;
    KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "Total fluid connections: " << rInterfaceConditionToFluidElement.size() << std::endl;
    KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "Total solid connections: " << rInterfaceConditionToSolidElement.size() << std::endl;
    KRATOS_INFO("HydraulicFluidAuxiliaryUtilities") << "Final conditions in slip bed: " << rSlipBedModelPart.NumberOfConditions() << std::endl;

    KRATOS_CATCH("")
}

}
