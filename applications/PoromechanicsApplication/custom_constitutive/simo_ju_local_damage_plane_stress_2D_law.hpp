//   
//   Project Name:        KratosPoromechanicsApplication $
//   Last Modified by:    $Author:    Ignasi de Pouplana $
//   Date:                $Date:           February 2016 $
//   Revision:            $Revision:                 1.0 $
//

#if !defined (KRATOS_SIMO_JU_LOCAL_DAMAGE_PLANE_STRESS_2D_LAW_H_INCLUDED)
#define  KRATOS_SIMO_JU_LOCAL_DAMAGE_PLANE_STRESS_2D_LAW_H_INCLUDED

// System includes
#include <cmath>

// Project includes
#include "includes/serializer.h"

// Application includes
#include "custom_constitutive/local_damage_plane_stress_2D_law.hpp"
#include "custom_constitutive/custom_hardening_laws/exponential_damage_hardening_law.hpp"
#include "custom_constitutive/custom_yield_criteria/simo_ju_yield_criterion.hpp"
#include "custom_constitutive/custom_flow_rules/local_damage_flow_rule.hpp"
#include "poromechanics_application_variables.h"

namespace Kratos
{

class KRATOS_API(POROMECHANICS_APPLICATION) SimoJuLocalDamagePlaneStress2DLaw : public LocalDamagePlaneStress2DLaw
{

public:

    KRATOS_CLASS_POINTER_DEFINITION(SimoJuLocalDamagePlaneStress2DLaw);

    typedef FlowRule::Pointer FlowRulePointer;
    typedef YieldCriterion::Pointer YieldCriterionPointer;
    typedef HardeningLaw::Pointer HardeningLawPointer;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    /// Default Constructor
    SimoJuLocalDamagePlaneStress2DLaw();
    
    /// Second Constructor
    SimoJuLocalDamagePlaneStress2DLaw(FlowRulePointer pFlowRule, YieldCriterionPointer pYieldCriterion, HardeningLawPointer pHardeningLaw); 
    
    /// Copy Constructor
    SimoJuLocalDamagePlaneStress2DLaw (const SimoJuLocalDamagePlaneStress2DLaw& rOther);

    /// Destructor
    ~SimoJuLocalDamagePlaneStress2DLaw() override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
    int Check(const Properties& rMaterialProperties, const GeometryType& rElementGeometry, const ProcessInfo& rCurrentProcessInfo) const override;
    
    ConstitutiveLaw::Pointer Clone() const override;
        
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    
protected:

    /// Member Variables
        
///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
    void CalculateCharacteristicSize( double& rCharacteristicSize, const GeometryType& DomainGeometry ) override;

///----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

private:
    
    /// Serialization
    
    friend class Serializer;

    void save(Serializer& rSerializer) const override
    {
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

    void load(Serializer& rSerializer) override
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, ConstitutiveLaw )
    }

}; // Class SimoJuLocalDamagePlaneStress2DLaw
}  // namespace Kratos.
#endif // KRATOS_SIMO_JU_LOCAL_DAMAGE_PLANE_STRESS_2D_LAW_H_INCLUDED  defined 
