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
//  Main authors:    Alejandro Cornejo & Lucia Barbu
//  Collaborator:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_VISCOUS_GENERALIZED_KELVIN_H_INCLUDED)
#define KRATOS_VISCOUS_GENERALIZED_KELVIN_H_INCLUDED

// System includes

// Project includes
#include "includes/constitutive_law.h"
#include "custom_constitutive/elastic_isotropic_3d.h"

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
 * @class ViscousGeneralizedKelvin
 * @ingroup StructuralMechanicsApplication
 * @brief This is a constitutive law that reproduces the behaviour of viscous Kelvin material
 * @details The definition of the Kelvin-Voigt material can be founf in Wikipedia https://en.wikipedia.org/wiki/Kelvin%E2%80%93Voigt_material
 * The Kelvin–Voigt model, also called the Voigt model, can be represented by a purely viscous damper and purely elastic spring
 *
 *                   Spring K
 *                  |--^^^^--|
 *             -----          ------
 *                  |---[----|
 *                   Damper C
 *
 * @param TElasticBehaviourLaw Defines the elastic behaviour of the constitutive law (can be hyperelastic or just linear elastic, or any desired elastic behaviour)
 * @author Alejandro Cornejo&  Lucia Barbu
 */

template<class TElasticBehaviourLaw>
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) ViscousGeneralizedKelvin
    : public TElasticBehaviourLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// Definition of the base CL class
    typedef ConstitutiveLaw CLBaseType;

    /// Definition of the base class
    typedef TElasticBehaviourLaw BaseType;

    /// The index definition
    typedef std::size_t IndexType;

    /// The size definition
    typedef std::size_t SizeType;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = TElasticBehaviourLaw::Dimension;

    /// Static definition of the VoigtSize
    static constexpr SizeType VoigtSize = TElasticBehaviourLaw::VoigtSize;

    /// Counted pointer of GenericYieldSurface
    KRATOS_CLASS_POINTER_DEFINITION(ViscousGeneralizedKelvin);

    /// The node definition
    typedef Node<3> NodeType;

    /// The geometry definition
    typedef Geometry<NodeType> GeometryType;

    /// Definition of the machine precision tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    ViscousGeneralizedKelvin();

    /**
    * @brief Clone.
    */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
    * @brief Copy constructor.
    */
    ViscousGeneralizedKelvin(const ViscousGeneralizedKelvin& rOther);

    /**
    * @brief Destructor.
    */
    ~ViscousGeneralizedKelvin() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief Computes the material response in terms of 1st Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of 2nd Piola-Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of Kirchhoff stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Computes the material response in terms of Cauchy stresses and constitutive tensor
     * @see Parameters
     */
    void CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief To be called at the end of each solution step (e.g. from Element::FinalizeSolutionStep)
     * @param rMaterialProperties the Properties instance of the current element
     * @param rElementGeometry the geometry of the current element
     * @param rShapeFunctionsValues the shape functions values in the current integration point
     * @param the current ProcessInfo instance
     */
    void FinalizeSolutionStep(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const Vector& rShapeFunctionsValues,
        const ProcessInfo& rCurrentProcessInfo
        ) override;

    /**
     * Finalize the material response in terms of 1st Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK1(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of 2nd Piola-Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponsePK2(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of Kirchhoff stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseKirchhoff(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * Finalize the material response in terms of Cauchy stresses
     * @see Parameters
     */
    void FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief Calculates the value of a specified variable (Vector)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Vector& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Vector>& rThisVariable,
        Vector& rValue
        ) override;

    /**
     * @brief Calculates the value of a specified variable (Matrix)
     * @param rParameterValues the needed parameters for the CL calculation
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    Matrix& CalculateValue(
        ConstitutiveLaw::Parameters& rParameterValues,
        const Variable<Matrix>& rThisVariable,
        Matrix& rValue
        ) override;

    /**
     * @brief This function provides the place to perform checks on the completeness of the input.
     * @details It is designed to be called only once (or anyway, not often) typically at the beginning
     * of the calculations, so to verify that nothing is missing from the input or that no common error is found.
     * @param rMaterialProperties The properties of the material
     * @param rElementGeometry The geometry of the element
     * @param rCurrentProcessInfo The current process info instance
     * @return 0 if OK, 1 otherwise
     */
    int Check(
        const Properties& rMaterialProperties,
        const GeometryType& rElementGeometry,
        const ProcessInfo& rCurrentProcessInfo
        ) const override;

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

    // Converged values
    Vector mPrevStressVector = ZeroVector(VoigtSize);
    Vector mPrevInelasticStrainVector = ZeroVector(VoigtSize);

    ///@}
    ///@name Private Operators
    ///@{

    ///@}
    ///@name Private Operations
    ///@{

    Vector& GetPreviousStressVector() { return mPrevStressVector; }
    void SetPreviousStressVector(const Vector& PrevStressVector) { mPrevStressVector = PrevStressVector; }

    Vector& GetPreviousInelasticStrainVector() { return mPrevInelasticStrainVector; }
    void SetPreviousInelasticStrainVector(const Vector& PrevInelasticStrainVector) { mPrevInelasticStrainVector = PrevInelasticStrainVector; }

    /**
     * @brief Compute visco-elasticity
     */
    void ComputeViscoElasticity(ConstitutiveLaw::Parameters& rValues);

    ///@}
    ///@name Private  Access
    ///@{

    ///@}
    ///@name Private Inquiry
    ///@{

    ///@}
    ///@name Un accessible methods
    ///@{

    // Serialization

    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.save("PrevStressVector", mPrevStressVector);
        rSerializer.save("PrevInelasticStrainVector", mPrevInelasticStrainVector);
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, ConstitutiveLaw)
        rSerializer.load("PrevStressVector", mPrevStressVector);
        rSerializer.load("PrevInelasticStrainVector", mPrevInelasticStrainVector);
    }

    ///@}

}; // Class GenericYieldSurface

} // namespace Kratos
#endif
