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

#pragma once

// System includes

// External includes

// Project includes
#include "geometries/geometry.h"
#include "includes/node.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "includes/global_pointer_variables.h"
#include "modified_shape_functions/modified_shape_functions.h"
#include "../../FluidDynamicsApplication/custom_utilities/fluid_auxiliary_utilities.h"

// Application includes
#include "hydraulic_fluid_auxiliary_utilities.h"

namespace Kratos
{
///@addtogroup FluidDynamicsHydraulicsApplication
///@{

///@name Kratos classes
///@{

class KRATOS_API(KRATOS_CORE) HydraulicFluidAuxiliaryUtilities
{
public:
    ///@name Type Definitions
    ///@{

    using SizeType = std::size_t;

    using IndexType = std::size_t;

    using NodeType = Node;

    using GeometryType = Geometry<NodeType>;

    using PointsArrayType = typename GeometryType::PointsArrayType;

    using ModifiedShapeFunctionsFactoryType = std::function<ModifiedShapeFunctions::UniquePointer(const GeometryType::Pointer, const Vector&)>;

    ///@}
    ///@name Static Operations
    ///@{
    /**
     * @brief This functions calculates the wetted perimeter of a given condition in order to apply the
     * corresponding inlet boundary condition.
     * @param rModelPart Inlet Model Part
     * @param rSkinFlag Flag that marks the conditions to be included in the calculation
     * @param rDistanceVariable Reference to the variable containing the distance
     * @param IsHistorical True if the distance is in the historical database, false otherwise
     * @return Wetted perimeter 
     */
    static double CalculateWettedPetimeter(
        ModelPart &rModelPart,
        const Flags &rSkinFlag,
        const Variable<double> &rDistanceVariable,
        const bool IsHistorical);
    ///@}
    ///@{
    /**
     * @brief This functions calculates the wetted area of a given condition in order to apply the
     * corresponding inlet boundary condition.
     * @param rModelPart Inlet Model Part
     * @param rSkinFlag Flag that marks the conditions to be included in the calculation
     * @param rDistanceVariable Reference to the variable containing the distance
     * @param IsHistorical True if the distance is in the historical database, false otherwise
     * @return Wetted area
     */
    static double CalculateWettedArea(
        ModelPart &rModelPart,
        const Flags &rSkinFlag,
        const Variable<double> &rDistanceVariable,
        const bool IsHistorical);

    /**
     * @brief Calculates initial water depth guess by taking the average between the maximum and minimum coordinates.
     *@param rModelPart Inlet Model Part
     * @return Initial water depth guess
     */
    static double InitialWaterDepth(ModelPart &rModelPart);

    /**
     * @brief
     * Assign the inlet velocity to all nodes that are wet in the input model part. For dry nodes(air)is assumed that the inlet velocity is null.
     * @param  rModelPart Inlet Model Part
     * @param  InletVelocity the velocity value to be assigned to wet nodes.
     * @param  rDistancesVariable Variable name of the inlet distance.
     */
    static void SetInletVelocity(ModelPart &rModelPart, double InletVelocity, const Variable<double> &rDistanceVariable);

    /**
     * @brief Free the inlet velocity in the nodes belonging to inlet model part.
     * @param  rModelPart Inlet Model Part
     */
    static void FreeInlet(ModelPart& rModelPart);

    /**
     * @brief  Set the free surface (DISTANCE) in the rModelPart equal to the water depth corresponding to Froude 1
     * @param  rModelPart Inlet Model Part
     * @param  rSkinFlag Flag that marks the conditions to be included in the calculation
     * @param  rDistancesVariable Variable name of the inlet distance.
     */
    static void SetInletFreeSurface(ModelPart &rModelPart, const Flags &rSkinFlag, const Variable<double> &rDistanceVariable);
    
    /**
     * @brief This function calculates the artificial viscosity for the given model part.
     * @param rModelPart Model Part
     * @param artificial_limiter_coefficient Coefficient for the artificial viscosity limiter
     */
    static void CalculateArtificialViscosity(ModelPart &rModelPart, double artificial_limiter_coefficient);

    /**
     * @brief Method to find elements neighboring conditions
     * @param rModelPart Model Part
     * @param rConditionFlag Flag that marks the conditions
     */
    static void FindElementsNeighbouringConditions(
        ModelPart &rModelPart,
        const Flags &rConditionFlag);

    /**
     * @brief Calculate erosion rate for a fluid element at the interface
     * @details This function calculates the erosion rate using Shields parameter
     * and threshold conditions based on Brownlie (1981) and van Rijn formulations.
     * @param rFluidElement The fluid element at the interface
     * @param rInterfaceCondition The interface condition
     * @param D50 Median particle diameter (m)
     * @param SedimentDensity Sediment density (kg/m³)
     * @return Volumetric erosion rate (m³/s)
     */
    static double CalculateErosionRate(
        const Element& rFluidElement,
        const Condition& rInterfaceCondition,
        const double D50 = 1e-04,
        const double SedimentDensity = 2650.0);

    /**
     * @brief Process erosion of a solid element
     * @details This function processes the erosion of a solid element by:
     * - Activating the solid element (converting it to fluid)
     * - Transferring averaged values from connected fluid element
     * - Calculating proper distance field values
     * - Creating new interface conditions for newly exposed solid faces
     * @param rSolidElement The solid element to be eroded
     * @param rSlipBedModelPart The slip bed model part containing interface conditions
     * @param rInterfaceCondition The interface condition being processed
     * @param rComputingModelPart The main computing model part
     * @param rConnectedFluidElement Connected fluid element for value transfer (can be null)
     * @return Number of new interface conditions created
     */
    static int ProcessElementErosion(
        Element& rSolidElement,
        ModelPart& rSlipBedModelPart,
        const Condition& rInterfaceCondition,
        ModelPart& rComputingModelPart,
        const Element* pConnectedFluidElement = nullptr);

    /**
     * @brief Process deposition of a fluid element
     * @details This function processes the deposition of a fluid element by:
     * - Deactivating the fluid element
     * - Adding all its faces (except the original interface face) as new conditions
     * @param rFluidElement The fluid element to be deactivated for deposition
     * @param rSlipBedModelPart The slip bed model part containing interface conditions
     * @param rInterfaceCondition The interface condition being processed
     * @param rComputingModelPart The main computing model part
     * @return Number of new interface conditions created
     */
    static int ProcessElementDeposition(
        Element& rFluidElement,
        ModelPart& rSlipBedModelPart,
        const Condition& rInterfaceCondition,
        ModelPart& rComputingModelPart);

    ///@}

private :

    struct EdgeDataContainer
    {
        NodeType::Pointer pNodeI = nullptr;
        NodeType::Pointer pNodeJ = nullptr;
        SizeType NumberOfRepetitions = 0;
    };

}; // namespace Kratos
}
