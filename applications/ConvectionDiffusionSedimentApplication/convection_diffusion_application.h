// KRATOS ___ ___  _  ___   __   ___ ___ ___ ___
//       / __/ _ \| \| \ \ / /__|   \_ _| __| __|
//      | (_| (_) | .` |\ V /___| |) | || _|| _|
//       \___\___/|_|\_| \_/    |___/___|_| |_|  APPLICATION
//
//  License: BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Riccardo Rossi
//

#if !defined(KRATOS_KRATOS_CONVECTION_DIFFUSION_SEDIMENT_APPLICATION_H_INCLUDED )
#define  KRATOS_KRATOS_CONVECTION_DIFFUSION_SEDIMENT_APPLICATION_H_INCLUDED


// System includes
#include <string>
#include <iostream>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "convection_diffusion_sediment_application_variables.h"

#include "custom_elements/axisymmetric_eulerian_convection_diffusion_sediment.h"
#include "custom_elements/conv_diff_sed_2d.h"
#include "custom_elements/conv_diff_sed_3d.h"
#include "custom_elements/eulerian_diff.h"
#include "custom_elements/eulerian_conv_diff_sed.h"
#include "custom_elements/laplacian_element.h"
#include "custom_elements/laplacian_shifted_boundary_element.h"
#include "custom_elements/mixed_laplacian_element.h"
#include "custom_elements/mixed_laplacian_shifted_boundary_element.h"
#include "custom_elements/embedded_laplacian_element.h"
#include "custom_elements/adjoint_diffusion_element.h"
#include "custom_elements/qs_convection_diffusion_sediment_explicit.h"
#include "custom_elements/d_convection_diffusion_sediment_explicit.h"

#include "custom_conditions/axisymmetric_thermal_face.h"
#include "custom_conditions/thermal_face.h"
#include "custom_conditions/flux_condition.h"
#include "custom_conditions/adjoint_thermal_face.h"
#include "custom_conditions/laplacian_shifted_boundary_condition.h"
#include "custom_conditions/mixed_laplacian_shifted_boundary_condition.h"

#include "includes/variables.h"
#include "includes/condition.h"


namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

///@}
///@name Kratos Classes
///@{

/**
 * @class KratosConvectionDiffusionSedimentApplication
 * @ingroup KratosConvectionDiffusionSedimentApplication
 * @brief The Convection Diffusion Sediment Application contains a series of elements and conditions and the corresponding strategies and solvers within Kratos Multiphysics necesaries in order to simulate a convection-diffusion-sediment problem
 * @details The application includes tests to check the proper functioning of the application. Features:
- A set of *Neumann* conditions:
     * Flux conditions
     * Thermal conditions
- Elements:
    * Laplacian element (both 2D/3D)
    * Laplacian embedded element (both 2D/3D)
    * Mixed Laplacian element (both 2D/3D)
    * Eulerian convection-diffusion-sediment (both 2D/3D)
    * Axisymmetric Eulerian convection-diffusion-sediment
    * Convection-diffusion-sediment (both 2D/3D)
    * Convection-diffusion-sediment with change of phase (2D)
- Strategies:
    * Non-linear/linear convection-diffusion-sediment strategy
    * Eulerian convection-diffusion-sediment strategy
    * Semi-Eulerian convection-diffusion-sediment strategy
- Utilities and others:
    * BFECC convection utility
    * BFECC elemental limiter convection utility
    * Convection particle
    * Face-heat utilities
    * Move particle utility
    * Pure convection tools
    * Pure convection (Crank-Nicolson) tools
 * @author Riccardo Rossi
 * @author Pablo Becker
 * @author Jordi Cotela
 * @author Ruben Zorrilla
 * @see KratosApplication
*/
class KRATOS_API(CONVECTION_DIFFUSION_SEDIMENT_APPLICATION) KratosConvectionDiffusionSedimentApplication : public KratosApplication
{
public:
    ///@name Type Definitions
    ///@{


    /// Pointer definition of KratosConvectionDiffusionSedimentApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosConvectionDiffusionSedimentApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosConvectionDiffusionSedimentApplication();

    /// Destructor.
    virtual ~KratosConvectionDiffusionSedimentApplication() {}


    ///@}
    ///@name Operators
    ///@{


    ///@}
    ///@name Operations
    ///@{

    void Register() override;



    ///@}
    ///@name Access
    ///@{


    ///@}
    ///@name Inquiry
    ///@{


    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    std::string Info() const override
    {
        return "KratosConvectionDiffusionSedimentApplication";
    }

    /// Print information about this object.
    void PrintInfo(std::ostream& rOStream) const override
    {
        rOStream << Info();
        PrintData(rOStream);
    }

    ///// Print object's data.
    void PrintData(std::ostream& rOStream) const override
    {
        KRATOS_WATCH("in KratosConvectionDiffusionSedimentApplication");
        KRATOS_WATCH(KratosComponents<VariableData>::GetComponents().size() );
        rOStream << "Variables:" << std::endl;
        KratosComponents<VariableData>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Elements:" << std::endl;
        KratosComponents<Element>().PrintData(rOStream);
        rOStream << std::endl;
        rOStream << "Conditions:" << std::endl;
        KratosComponents<Condition>().PrintData(rOStream);
    }


    ///@}
    ///@name Friends
    ///@{


    ///@}

protected:
    ///@name Protected static Member Variables
    ///@{


    ///@}
    ///@name Protected member Variables
    ///@{


    ///@}
    ///@name Protected Operators
    ///@{


    ///@}
    ///@name Protected Operations
    ///@{


    ///@}
    ///@name Protected  Access
    ///@{


    ///@}
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{

    ///@}
    ///@name Member Variables
    ///@{

    const AxisymmetricEulerianConvectionDiffusionSedimentElement<2,3>  mAxisymmetricEulerianConvectionDiffusionSediment2D3N;
    const AxisymmetricEulerianConvectionDiffusionSedimentElement<2,4>  mAxisymmetricEulerianConvectionDiffusionSediment2D4N;

    const EulerianConvectionDiffusionSedimentElement<2,3>  mEulerianConvDiff2D3N;
    const EulerianConvectionDiffusionSedimentElement<2,4>  mEulerianConvDiff2D4N;
    const EulerianConvectionDiffusionSedimentElement<3,4>  mEulerianConvDiff3D4N;
    const EulerianConvectionDiffusionSedimentElement<3,8>  mEulerianConvDiff3D8N;
    const EulerianDiffusionElement<2,3>  mEulerianDiffusion2D3N;
    const EulerianDiffusionElement<3,4>  mEulerianDiffusion3D4N;

    const ConvDiff2D  mConvDiff2D;
    const ConvDiff3D  mConvDiff3D;
    const LaplacianElement mLaplacian2D3N;
    const LaplacianElement mLaplacian3D4N;
    const LaplacianElement mLaplacian3D8N;
    const LaplacianElement mLaplacian3D27N;

    const LaplacianShiftedBoundaryElement<2> mLaplacianShiftedBoundary2D3N;
    const LaplacianShiftedBoundaryElement<3> mLaplacianShiftedBoundary3D4N;

    const EmbeddedLaplacianElement<2> mEmbeddedLaplacian2D3N;
    const EmbeddedLaplacianElement<3> mEmbeddedLaplacian3D4N;

    const MixedLaplacianElement<2,3> mMixedLaplacianElement2D3N;
    const MixedLaplacianElement<3,4> mMixedLaplacianElement3D4N;
    const MixedLaplacianShiftedBoundaryElement<2> mMixedLaplacianShiftedBoundary2D3N;
    const MixedLaplacianShiftedBoundaryElement<3> mMixedLaplacianShiftedBoundary3D4N;

    const AdjointDiffusionElement<LaplacianElement> mAdjointDiffusionElement2D3N;
    const AdjointDiffusionElement<LaplacianElement> mAdjointDiffusionElement3D4N;

    const AxisymmetricThermalFace mAxisymmetricThermalFace2D2N;
    const ThermalFace mThermalFace2D2N;
    const ThermalFace mThermalFace3D3N;
    const ThermalFace mThermalFace3D4N;
    const FluxCondition<2>  mFluxCondition2D2N;
    const FluxCondition<3>  mFluxCondition3D3N;
    const FluxCondition<4>  mFluxCondition3D4N;

    const LaplacianShiftedBoundaryCondition mLaplacianShiftedBoundaryCondition;
    const MixedLaplacianShiftedBoundaryCondition mMixedLaplacianShiftedBoundaryCondition;

    const AdjointThermalFace mAdjointThermalFace2D2N;
    const AdjointThermalFace mAdjointThermalFace3D3N;

    const QSConvectionDiffusionSedimentExplicit<2,3> mQSConvectionDiffusionSedimentExplicit2D3N;
    const QSConvectionDiffusionSedimentExplicit<3,4> mQSConvectionDiffusionSedimentExplicit3D4N;
    const DConvectionDiffusionSedimentExplicit<2,3> mDConvectionDiffusionSedimentExplicit2D3N;
    const DConvectionDiffusionSedimentExplicit<3,4> mDConvectionDiffusionSedimentExplicit3D4N;

    ///@}
    ///@name Private Operators
    ///@{


    ///@}
    ///@name Private Operations
    ///@{


    ///@}
    ///@name Private  Access
    ///@{


    ///@}
    ///@name Private Inquiry
    ///@{


    ///@}
    ///@name Un accessible methods
    ///@{

    /// Assignment operator.
    KratosConvectionDiffusionSedimentApplication& operator=(KratosConvectionDiffusionSedimentApplication const& rOther);

    /// Copy constructor.
    KratosConvectionDiffusionSedimentApplication(KratosConvectionDiffusionSedimentApplication const& rOther);


    ///@}

}; // Class KratosConvectionDiffusionSedimentApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_KRATOS_CONVECTION_DIFFUSION_SEDIMENT_APPLICATION_H_INCLUDED  defined
