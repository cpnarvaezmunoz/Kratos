// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Vahid Galavi
//

#if !defined (KRATOS_LINEAR_ELASTIC_2D_INTERFACE_LAW_GEO_H_INCLUDED)
#define  KRATOS_LINEAR_ELASTIC_2D_INTERFACE_LAW_GEO_H_INCLUDED

// System includes

// External includes

// Project includes
#include "custom_constitutive/linear_elastic_plane_strain_2D_law.h"

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
 * @class LinearElastic2DInterfaceLaw
 * @ingroup GeoMechanicsApplication
 * @brief This class defines a small deformation linear elastic constitutive model for plane strain cases
 * @details This class derives from the linear elastic case on 3D
 * @author Vahid Galavi
 */
class KRATOS_API(GEO_MECHANICS_APPLICATION) LinearElastic2DInterfaceLaw
    : public GeoLinearElasticPlaneStrain2DLaw
{
public:
    ///@name Type Definitions
    ///@{

    /// The base class ConstitutiveLaw type definition
    typedef ConstitutiveLaw       CLBaseType;

    /// The base class ElasticIsotropicK03DLaw type definition
    typedef GeoLinearElasticPlaneStrain2DLaw      BaseType;

    /// The size type definition
    typedef std::size_t             SizeType;

    /// Static definition of the dimension
    static constexpr SizeType Dimension = N_DIM_2D;

    /// Counted pointer of GeoLinearElasticPlaneStrain2DLaw
    KRATOS_CLASS_POINTER_DEFINITION( LinearElastic2DInterfaceLaw );

    ///@name Life Cycle
    ///@{

    /**
     * @brief Default constructor.
     */
    LinearElastic2DInterfaceLaw();

    /**
     * @brief The clone operation
     */
    ConstitutiveLaw::Pointer Clone() const override;

    /**
     * Copy constructor.
     */
    LinearElastic2DInterfaceLaw (const LinearElastic2DInterfaceLaw& rOther);


    /**
     * @brief Destructor.
     */
    ~LinearElastic2DInterfaceLaw() override;

    ///@}
    ///@name Operators
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This function is designed to be called once to check compatibility with element
     * @param rFeatures: The Features of the law
     */
    void GetLawFeatures(Features& rFeatures) override;

    /**
     * @brief Dimension of the law:
     * @return The dimension were the law is working
     */
    SizeType WorkingSpaceDimension() override
    {
        return Dimension;
    }

    /**
     * @brief Voigt tensor size:
     * @return The size of the strain vector in Voigt notation
     */
    SizeType GetStrainSize() const override
    {
        return VOIGT_SIZE_2D_INTERFACE;
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

    /**
     * @brief  Itreturns the value of a specified variable
     * @param rThisVariable the variable to be returned
     * @param rValue a reference to the returned value
     * @param rValue output: the value of the specified variable
     */
    bool& GetValue(const Variable<bool>& rThisVariable, bool& rValue) override;

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

    // /**
    //  * @brief It calculates the constitutive matrix C
    //  * @param C The constitutive matrix
    //  * @param rValues Parameters of the constitutive law
    //  */
    void CalculateElasticMatrix(Matrix& C, ConstitutiveLaw::Parameters& rValues) override;

    /**
     * @brief It calculates the stress vector
     * @param rStrainVector The strain vector in Voigt notation
     * @param rStressVector The stress vector in Voigt notation
     * @param rValues Parameters of the constitutive law
     */
    void CalculatePK2Stress(const Vector& rStrainVector,
                            Vector& rStressVector,
                            ConstitutiveLaw::Parameters& rValues) override;

    // /**
    //  * @brief It calculates the strain vector
    //  * @param rValues The internal values of the law
    //  * @param rStrainVector The strain vector in Voigt notation
    //  */
    // void CalculateCauchyGreenStrain(ConstitutiveLaw::Parameters& rValues, Vector& rStrainVector) override;

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

    ///@}
    ///@name Private  Access
    ///@{
    ///@}

    ///@}
    ///@name Serialization
    ///@{
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, GeoLinearElasticPlaneStrain2DLaw)
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, GeoLinearElasticPlaneStrain2DLaw)
    }

    // stress vector indices
    // const int VOIGT_INDEX_XX = 0;
    // const int VOIGT_INDEX_YY = 1;

}; // Class GeoLinearElasticPlaneStrain2DLaw
}  // namespace Kratos.
#endif // KRATOS_LINEAR_PLANE_STRAIN_K0_LAW_H_INCLUDED  defined
