// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    cmcmcmcm
//
//

// System includes

// External includes

// Project includes
//#include "custom_utilities/tangent_operator_calculator_utility.h"
#include "custom_utilities/automatic_differentiation_tangent_utilities.h"
#include "constitutive_laws_application_variables.h"
#include "damage_cl.h"
//#include "custom_constitutive/auxiliary_files/cl_integrators/generic_cl_integrator_damage.h"


namespace Kratos
{


//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

DamageCl::DamageCl()
  //  : ElasticIsotropic3D()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

DamageCl::DamageCl(const DamageCl& rOther)
  //  : ElasticIsotropic3D(rOther)
{
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer DamageCl::Clone() const
{
    return Kratos::make_shared<DamageCl>(*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

DamageCl::~DamageCl()
{
}


void DamageCl::CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void DamageCl::CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void DamageCl::CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    this->CalculateMaterialResponseCauchy(rValues);
}


/***********************************************************************************/
/***********************************************************************************/

void DamageCl::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // We get the strain vector
    auto& r_strain_vector = rValues.GetStrainVector();

    // We compute the stress
    auto& r_integrated_stress_vector = rValues.GetStressVector();
    // Elastic Matrix
    auto& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
    this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

    // Converged values
    double threshold = mThreshold;
    double damage =mDamage;

    // S0 = C:(E-E0) + S0
    BoundedArrayType predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector);
      
    /***********************************************************************************/
    /***********************************************************************************/
    //CalculateEquivalentStress
    const double uniaxial_stress = ConstitutiveLawUtilities<6>::CalculateVonMisesEquivalentStress(predictive_stress_vector);

    /**
    * @brief This method the uniaxial equivalent stress of Huber-VonMises
    * @param rPredictiveStressVector The predictive stress vector S = C:(E-Ep)
    * @param rStrainVector The StrainVector vector
    * @param rValues Parameters of the constitutive law
    */

    const double F = uniaxial_stress - threshold;

    if (F <= threshold_tolerance) { // Elastic case
        noalias(r_integrated_stress_vector) = (1.0 - damage) * predictive_stress_vector;
        //if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_CONSTITUTIVE_TENSOR)) {
            r_constitutive_matrix *= (1.0 - damage);
        //}

    } else { // Damage case
        
        const double characteristic_length = AdvancedConstitutiveLawUtilities<6>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
            // This routine updates the PredictiveStress to verify the yield surf
            DamageCl::IntegrateStressVector(predictive_stress_vector, uniaxial_stress, damage, threshold, rValues, characteristic_length);

            // Updated Values
            noalias(r_integrated_stress_vector) = predictive_stress_vector;
            rValues.GetConstitutiveMatrix() *= (1.0 - mDamage);
        }
    } // End CalculateMaterialResponseCauchy

/***********************************************************************************/
/***********************************************************************************/

void DamageCl::GetInitialUniaxialThreshold(
        ConstitutiveLaw::Parameters& rValues,
        double& rThreshold
        )
    {
        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double yield_tension = r_material_properties.Has(YIELD_STRESS) ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_TENSION];
		rThreshold = std::abs(yield_tension);
    }


/***********************************************************************************/
/***********************************************************************************/

    void DamageCl::CalculateDamageParameter(
                ConstitutiveLaw::Parameters& rValues,
                double& rAParameter,
                const double CharacteristicLength
                )
        {
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        const double fracture_energy = r_material_properties[FRACTURE_ENERGY];
        const double young_modulus = r_material_properties[YOUNG_MODULUS];
        const double yield_compression = r_material_properties.Has(YIELD_STRESS) ? r_material_properties[YIELD_STRESS] : r_material_properties[YIELD_STRESS_COMPRESSION];

        rAParameter = 1.00 / (fracture_energy * young_modulus / (CharacteristicLength * std::pow(yield_compression, 2)) - 0.5);

        }
    void DamageCl::CalculateExponentialDamage(
                    const double UniaxialStress,
                    const double Threshold,
                    const double DamageParameter,
                    const double CharacteristicLength,
                    ConstitutiveLaw::Parameters& rValues,
                    double& rDamage
                    )
            {
            double initial_threshold;
            DamageCl::GetInitialUniaxialThreshold(rValues, initial_threshold);
            rDamage = 1.0 - (initial_threshold / UniaxialStress) * std::exp(DamageParameter *
                    (1.0 - UniaxialStress / initial_threshold));
            }           

    void DamageCl::IntegrateStressVector(
                array_1d<double, VoigtSize>& rPredictiveStressVector,
                const double UniaxialStress,
                double& rDamage,
                double& rThreshold,
                ConstitutiveLaw::Parameters& rValues,
                const double CharacteristicLength
                )
                {   
        const Properties& r_material_properties = rValues.GetMaterialProperties();

        const int softening_type = r_material_properties[SOFTENING_TYPE];
        double damage_parameter;
        CalculateDamageParameter(rValues, damage_parameter, CharacteristicLength);

        CalculateExponentialDamage(UniaxialStress, rThreshold, damage_parameter, CharacteristicLength, rValues, rDamage);

        rDamage = (rDamage > 0.99999) ? 0.99999 : rDamage;
        rDamage = (rDamage < 0.0) ? 0.0 : rDamage;
        rPredictiveStressVector *= (1.0 - rDamage);
        }
/***********************************************************************************/
/***********************************************************************************/

void DamageCl::InitializeMaterial(
    const Properties& rMaterialProperties,
    const GeometryType& rElementGeometry,
    const Vector& rShapeFunctionsValues
    )
{
    // We construct the CL parameters
    ProcessInfo dummy_process_info;
    ConstitutiveLaw::Parameters aux_param(rElementGeometry, rMaterialProperties, dummy_process_info);

    // We call the integrator
    double initial_threshold;
    DamageCl::GetInitialUniaxialThreshold(aux_param, initial_threshold);
    this->SetThreshold(initial_threshold);
}

/***********************************************************************************/
/***********************************************************************************/

void DamageCl::FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void DamageCl::FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void DamageCl::FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues)
{
    // Small deformation so we can call the Cauchy method
    FinalizeMaterialResponseCauchy(rValues);
}

/***********************************************************************************/
/***********************************************************************************/

void DamageCl::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    // Integrate Stress Damage
    const Flags& r_constitutive_law_options = rValues.GetOptions();

    // We get the strain vector
    Vector& r_strain_vector = rValues.GetStrainVector();

    //NOTE: SINCE THE ELEMENT IS IN SMALL STRAINS WE CAN USE ANY STRAIN MEASURE. HERE EMPLOYING THE CAUCHY_GREEN
    //if (r_constitutive_law_options.IsNot( ConstitutiveLaw::USE_ELEMENT_PROVIDED_STRAIN)) {
        //::CalculateCauchyGreenStrain( rValues, r_strain_vector);
    //}

    // We compute the stress
    if (r_constitutive_law_options.Is(ConstitutiveLaw::COMPUTE_STRESS)) {
        // Elastic Matrix
        Matrix& r_constitutive_matrix = rValues.GetConstitutiveMatrix();
        this->CalculateValue(rValues, CONSTITUTIVE_MATRIX, r_constitutive_matrix);

        //this->template AddInitialStrainVectorContribution<Vector>(r_strain_vector);

        // S0 = C:(E-E0) + S0
        BoundedArrayType predictive_stress_vector = prod(r_constitutive_matrix, r_strain_vector);
        //this->template AddInitialStressVectorContribution<BoundedArrayType>(predictive_stress_vector);

        // Initialize Plastic Parameters
        double uniaxial_stress;
        //TConstLawIntegratorType::YieldSurfaceType::CalculateEquivalentStress(predictive_stress_vector, r_strain_vector, uniaxial_stress, rValues);
        ConstitutiveLawUtilities<6>::CalculateVonMisesEquivalentStress(predictive_stress_vector);
        const double F = uniaxial_stress - mThreshold;

        if (F >= threshold_tolerance) { // Damage case
            const double characteristic_length = AdvancedConstitutiveLawUtilities<VoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(rValues.GetElementGeometry());
            // This routine updates the PredictiveStress to verify the yield surf
            DamageCl::IntegrateStressVector(predictive_stress_vector, uniaxial_stress, mDamage, mThreshold, rValues, characteristic_length);
            mThreshold = uniaxial_stress;
        }
    }


}

//     int Check(
//         const Properties& rMaterialProperties,
//         const GeometryType& rElementGeometry,
//         const ProcessInfo& rCurrentProcessInfo
//         ) const 
        
//         {
//         const int check_base = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
//         //const int check_integrator = Da::Check(rMaterialProperties);
//         KRATOS_ERROR_IF_NOT(VoigtSize == this->GetStrainSize()) << "You are combining not compatible constitutive laws" << std::endl;
//         if ((check_base) > 0) return 1;
//             return 0;
// }

} // namespace Kratos
