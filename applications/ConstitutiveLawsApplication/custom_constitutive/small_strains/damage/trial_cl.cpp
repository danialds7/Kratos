// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    cmcmcmcm
//
// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_utilities/automatic_differentiation_tangent_utilities.h"
#include "constitutive_laws_application_variables.h"
#include "trial_cl.h"
#include "structural_mechanics_application_variables.h"
#include "constitutive_laws_application_variables.h"

//#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

TrialCl::TrialCl()
    : ElasticIsotropic3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

TrialCl::TrialCl(const TrialCl& rOther)
    : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer TrialCl::Clone() const
{
    return Kratos::make_shared<TrialCl>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

TrialCl::~TrialCl()
{
}

/***********************************************************************************/
/***********************************************************************************/

bool& TrialCl::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    // This Constitutive Law has been checked with Stenberg Stabilization
    if (rThisVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE) {
        rValue = true;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Matrix& TrialCl::GetValue(const Variable<Matrix>& rThisVariable, Matrix& rValue)
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

Vector& TrialCl::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/

double& TrialCl::GetValue(const Variable<double>& rThisVariable, double& rValue)
{
    return rValue;
}

//*************************CONSTITUTIVE LAW GENERAL FEATURES *************************
/***********************************************************************************/

void TrialCl::GetLawFeatures(Features& rFeatures)
{
    //Set the type of law
    rFeatures.mOptions.Set( PLANE_STRAIN_LAW);
    rFeatures.mOptions.Set( INFINITESIMAL_STRAINS );
    rFeatures.mOptions.Set( ISOTROPIC );

    //Set strain measure required by the consitutive law
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);

    //Set the strain size
    rFeatures.mStrainSize = 6;

    //Set the spacedimension
    rFeatures.mSpaceDimension = 3;
}

/***********************************************************************************/
/***********************************************************************************/

void TrialCl::CalculateElasticMatrix(VoigtSizeMatrixType& rC, ConstitutiveLaw::Parameters& rValues)
{
    const auto &r_props = rValues.GetMaterialProperties();
    const double E = r_props[YOUNG_MODULUS];
    const double NU = r_props[POISSON_RATIO];
    ConstitutiveLawUtilities<6>::CalculateElasticMatrix(rC, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

void TrialCl::CalculatePK2Stress(
    const ConstitutiveLaw::StrainVectorType& rStrainVector,
    ConstitutiveLaw::StressVectorType& rStressVector,
    ConstitutiveLaw::Parameters& rValues
    )
{
    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    ConstitutiveLawUtilities<6>::CalculatePK2StressFromStrain(rStressVector, rStrainVector, E, NU);
}

/***********************************************************************************/
/***********************************************************************************/

void TrialCl::CalculateCauchyGreenStrain(Parameters& rValues, ConstitutiveLaw::StrainVectorType& rStrainVector)
{
    ConstitutiveLawUtilities<6>::CalculateCauchyGreenStrain(rValues, rStrainVector);
}

/***********************************************************************************/
/***********************************************************************************/

void  TrialCl::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)

{
    KRATOS_TRY;

    const Properties& r_material_properties = rValues.GetMaterialProperties();
    const double E = r_material_properties[YOUNG_MODULUS];
    const double NU = r_material_properties[POISSON_RATIO];

    ConstitutiveLaw::StrainVectorType &r_strain_vector = rValues.GetStrainVector();
    ConstitutiveLaw::StressVectorType &r_stress_vector = rValues.GetStressVector();
    ConstitutiveLaw::VoigtSizeMatrixType &r_constitutive_matrix = rValues.GetConstitutiveMatrix();

    CalculateCauchyGreenStrain(rValues, r_strain_vector);
    AddInitialStrainVectorContribution<StrainVectorType>(r_strain_vector);

    const double c1 = E / ((1.0 + NU) * (1.0 - 2.0 * NU));
    const double c2 = c1 * (1 - NU);
    const double c3 = c1 * NU;
    const double c4 = c1 * 0.5 * (1 - 2 * NU);

    r_constitutive_matrix(0, 0) = c2;
    r_constitutive_matrix(0, 1) = c3;
    r_constitutive_matrix(0, 2) = c3;
    r_constitutive_matrix(1, 0) = c3;
    r_constitutive_matrix(1, 1) = c2;
    r_constitutive_matrix(1, 2) = c3;
    r_constitutive_matrix(2, 0) = c3;
    r_constitutive_matrix(2, 1) = c3;
    r_constitutive_matrix(2, 2) = c2;
    r_constitutive_matrix(3, 3) = c4;
    r_constitutive_matrix(4, 4) = c4;
    r_constitutive_matrix(5, 5) = c4;

    r_stress_vector[0] = c2 * r_strain_vector[0] + c3 * r_strain_vector[1] + c3 * r_strain_vector[2];
    r_stress_vector[1] = c3 * r_strain_vector[0] + c2 * r_strain_vector[1] + c3 * r_strain_vector[2];
    r_stress_vector[2] = c3 * r_strain_vector[0] + c3 * r_strain_vector[1] + c2 * r_strain_vector[2];
    r_stress_vector[3] = c4 * r_strain_vector[3];
    r_stress_vector[4] = c4 * r_strain_vector[4];
    r_stress_vector[5] = c4 * r_strain_vector[5];


    KRATOS_CATCH("");
}


} // Namespace Kratos
