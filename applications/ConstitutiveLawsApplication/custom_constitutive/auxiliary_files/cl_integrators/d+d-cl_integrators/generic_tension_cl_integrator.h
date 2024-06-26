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
//  Main authors:    Alejandro Cornejo 
//

#pragma once

// System includes

// Project includes
#include "includes/define.h"
#include "includes/checks.h"
#include "includes/serializer.h"
#include "includes/properties.h"
#include "utilities/math_utils.h"
#include "constitutive_laws_application_variables.h"

namespace Kratos
{

///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    // The size type definition
    typedef std::size_t SizeType;
    
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
 * @class GenericTensionConstitutiveLawIntegratorDplusDminusDamage
 * @ingroup StructuralMechanicsApplication
 * @brief: This object integrates the predictive stress using the isotropic the d+d- damage theory 
 * @details The definitions of these classes is completely static, the derivation is done in a static way
 * The damage integrator requires the definition of the following properties:
 * - SOFTENING_TYPE: The fosftening behaviour considered (linear, exponential,etc...)
 * @tparam TYieldSurfaceType The yield surface considered
 * @author Alejandro Cornejo 
 */
template <class TYieldSurfaceType>
class GenericTensionConstitutiveLawIntegratorDplusDminusDamage
{
  public:

    ///@name Type Definitions
    ///@{

    /// The type of yield surface
    typedef TYieldSurfaceType YieldSurfaceType;

    /// The define the working dimension size, already defined in the yield surface
    static constexpr SizeType Dimension = YieldSurfaceType::Dimension;

    /// The define the Voigt size, already defined in the yield surface
    static constexpr SizeType VoigtSize = YieldSurfaceType::VoigtSize;
    
    /// The type of plastic potential
    typedef typename YieldSurfaceType::PlasticPotentialType PlasticPotentialType;

    /// Counted pointer of GenericTensionConstitutiveLawIntegratorDplusDminusDamage
    KRATOS_CLASS_POINTER_DEFINITION(GenericTensionConstitutiveLawIntegratorDplusDminusDamage);

    /// Initialization constructor
    GenericTensionConstitutiveLawIntegratorDplusDminusDamage()
    {
    }

    /// Copy constructor
    GenericTensionConstitutiveLawIntegratorDplusDminusDamage(GenericTensionConstitutiveLawIntegratorDplusDminusDamage const &rOther)
    {
    }

    /// Assignment operator
    GenericTensionConstitutiveLawIntegratorDplusDminusDamage &operator=(GenericTensionConstitutiveLawIntegratorDplusDminusDamage const &rOther)
    {
        return *this;
    }

    /// Destructor
    virtual ~GenericTensionConstitutiveLawIntegratorDplusDminusDamage()
    {
    }

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method integrates the predictive stress vector with the CL using linear or exponential softening
     * @param PredictiveStressVector The predictive stress vector
     * @param UniaxialStress The equivalent uniaxial stress
     * @param Damage The internal variable of the damage model
     * @param Threshold The maximum uniaxial stress achieved previously
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void IntegrateStressVector(
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
        TYieldSurfaceType::CalculateDamageParameter(rValues, damage_parameter, CharacteristicLength);

        switch (softening_type)
        {
        case static_cast<int>(SofteningType::Linear):
            CalculateLinearDamage(UniaxialStress, rThreshold, damage_parameter, CharacteristicLength, rValues, rDamage);
            break;
        case static_cast<int>(SofteningType::Exponential):
            CalculateExponentialDamage(UniaxialStress, rThreshold, damage_parameter, CharacteristicLength, rValues, rDamage);
            break;
        default:
            KRATOS_ERROR << "SOFTENING_TYPE not defined or wrong..." << softening_type << std::endl;
            break;
        }
        rPredictiveStressVector *= (1.0 - rDamage);
    }

    /**
     * @brief This computes the damage variable according to exponential softening
     * @param UniaxialStress The equivalent uniaxial stress
     * @param Threshold The maximum uniaxial stress achieved previously
     * @param rDamage The internal variable of the damage model
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateExponentialDamage(
        const double UniaxialStress,
        const double Threshold,
        const double DamageParameter,
        const double CharacteristicLength,
        ConstitutiveLaw::Parameters& rValues,
        double& rDamage
        )
    {
        double initial_threshold;
        TYieldSurfaceType::GetInitialUniaxialThreshold(rValues, initial_threshold);
        rDamage = 1.0 - (initial_threshold / UniaxialStress) * std::exp(DamageParameter *
                 (1.0 - UniaxialStress / initial_threshold));
    }

    /**
     * @brief This computes the damage variable according to linear softening
     * @param UniaxialStress The equivalent uniaxial stress
     * @param Threshold The maximum uniaxial stress achieved previously
     * @param rDamage The internal variable of the damage model
     * @param rValues Parameters of the constitutive law
     * @param CharacteristicLength The equivalent length of the FE
     */
    static void CalculateLinearDamage(
        const double UniaxialStress,
        const double Threshold,
        const double DamageParameter,
        const double CharacteristicLength,
        ConstitutiveLaw::Parameters& rValues,
        double& rDamage
        )
    {
        double initial_threshold;
        TYieldSurfaceType::GetInitialUniaxialThreshold(rValues, initial_threshold);
        rDamage = (1.0 - initial_threshold / UniaxialStress) / (1.0 + DamageParameter);
    }

    /**
     * @brief This method returns the initial uniaxial stress threshold
     * @param rThreshold The uniaxial stress threshold
     * @param rValues Parameters of the constitutive law
     */
    static void GetInitialUniaxialThreshold(
        ConstitutiveLaw::Parameters& rValues,
        double& rThreshold
        )
    {
        TYieldSurfaceType::GetInitialUniaxialThreshold(rValues, rThreshold);
    }

    /**
     * @brief This method defines in the CL integrator
     * @return 0 if OK, 1 otherwise
     */
    static int Check(const Properties& rMaterialProperties)
    {
        KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(SOFTENING_TYPE)) << "SOFTENING_TYPE is not a defined value" << std::endl;
        return TYieldSurfaceType::Check(rMaterialProperties);
    }

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

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

    ///@}

};
} // namespace Kratos
