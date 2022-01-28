//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    @{KRATOS_APP_AUTHOR}
//


#if !defined(KRATOS_EHDTEST_APPLICATION_H_INCLUDED )
#define  KRATOS_EHDTEST_APPLICATION_H_INCLUDED


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "ehdtest_application_variables.h"

#include "custom_elements/electric_element.h"
#include "custom_elements/electric_distance_element.h"
#include "custom_elements/electric_enriched_element.h"
//#include "custom_elements/electric_enriched_n_element.h"
//#include "custom_elements/electric_enrichedt_element.h"

#include "custom_conditions/efield_condition.h"
#include "custom_conditions/point_charge_condition.h"

#include "includes/variables.h"
#include "includes/condition.h"


namespace Kratos {

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

/// Short class definition.
/** Detail class definition.
*/
class KRATOS_API(EHDTEST_APPLICATION) KratosEhdtestApplication : public KratosApplication {
public:
    ///@name Type Definitions
    ///@{

    /// Pointer definition of KratosEhdtestApplication
    KRATOS_CLASS_POINTER_DEFINITION(KratosEhdtestApplication);

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor.
    KratosEhdtestApplication();

    /// Destructor.
    ~KratosEhdtestApplication() override {}

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
        return "KratosEhdtestApplication";
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
          KRATOS_WATCH("in my application");
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


    ///@}julio
    ///@name Protected Inquiry
    ///@{


    ///@}
    ///@name Protected LifeCycle
    ///@{


    ///@}

private:
    ///@name Static Member Variables
    ///@{

    // static const ApplicationCondition  msApplicationCondition;

    ///@}
    ///@name Member Variables
    ///@{

    const ElectricElement mElectricElement;
    const ElectricDistanceElement mElectricDistanceElement;
    const ElectricEnrichedElement mElectricEnrichedElement;
    //const ElectricEnrichedNElement mElectricEnrichedNElement;
    //const ElectricEnrichedtElement mElectricEnrichedtElement;
    const EfieldCondition mEfieldCondition;
    const PointChargeCondition  mPointChargeCondition;

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
    KratosEhdtestApplication& operator=(KratosEhdtestApplication const& rOther);

    /// Copy constructor.
    KratosEhdtestApplication(KratosEhdtestApplication const& rOther);


    ///@}

}; // Class KratosEhdtestApplication

///@}


///@name Type Definitions
///@{


///@}
///@name Input and output
///@{

///@}


}  // namespace Kratos.

#endif // KRATOS_EHDTEST_APPLICATION_H_INCLUDED  defined
