/////////////////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
// Date: July 2024
/////////////////////////////////////////////////////////////

// Project includes
#include "custom_elements/polyhedron_contact_element.h"
#include "utilities/math_utils.h"
#include "DEM_application_variables.h"

namespace Kratos
{
//************************************************************************************
//************************************************************************************
PolyhedronContactElement::PolyhedronContactElement(IndexType NewId, 
                                                    PolyhedronParticle* PolyhedronParticle1, 
                                                    PolyhedronParticle* PolyhedronParticle2): 
mId(NewId),
mPolyhedronParticle1(PolyhedronParticle1),
mPolyhedronParticle2(PolyhedronParticle2)
{

}

PolyhedronContactElement::PolyhedronContactElement(IndexType NewId, 
                                                    PolyhedronParticle* PolyhedronParticle1, 
                                                    PolyhedronParticle* PolyhedronParticle2, 
                                                    PropertiesType::Pointer pProperties):
mId(NewId),
mPolyhedronParticle1(PolyhedronParticle1),
mPolyhedronParticle2(PolyhedronParticle2),
mpProperties(pProperties)
{

}

//create contact elements instances.
PolyhedronContactElement::Pointer PolyhedronContactElement::Create(IndexType NewId, PolyhedronParticle* PolyhedronParticle1, PolyhedronParticle* PolyhedronParticle2, PropertiesType::Pointer pProperties) const
{
    return PolyhedronContactElement::Pointer ( new PolyhedronContactElement(NewId, PolyhedronParticle1, PolyhedronParticle2, pProperties));
}

PolyhedronContactElement::~PolyhedronContactElement()
{
}

void PolyhedronContactElement::Initialize(const ProcessInfo& r_process_info) {
    KRATOS_TRY

    mFailureState = 0.0;
    mContactForce[0] = 0.0;
    mContactForce[1] = 0.0;
    mContactForce[2] = 0.0;
    mRotationalMoment[0] = 0.0;
    mRotationalMoment[1] = 0.0;
    mRotationalMoment[2] = 0.0;

    /*
    array_1d<double, 3> vector_of_zeros(3, 0.0);
    this->SetValue(LOCAL_CONTACT_FORCE, vector_of_zeros);
    this->SetValue(ELASTIC_LOCAL_ROTATIONAL_MOMENT, vector_of_zeros);
    this->SetValue(FAILURE_CRITERION_STATE, 0.0);
    this->SetValue(CONTACT_RADIUS, 0.0);*/

    KRATOS_CATCH( "" )
}

void PolyhedronContactElement::PrepareForPrinting() {
    KRATOS_TRY

    /*
    this->GetValue(LOCAL_CONTACT_FORCE)[0] = mContactForce[0];
    this->GetValue(LOCAL_CONTACT_FORCE)[1] = mContactForce[1];
    this->GetValue(LOCAL_CONTACT_FORCE)[2] = mContactForce[2];
    this->GetValue(ELASTIC_LOCAL_ROTATIONAL_MOMENT)[0] = mContactForce[0];
    this->GetValue(ELASTIC_LOCAL_ROTATIONAL_MOMENT)[1] = mRotationalMoment[1];
    this->GetValue(ELASTIC_LOCAL_ROTATIONAL_MOMENT)[2] = mRotationalMoment[2];
    this->GetValue(FAILURE_CRITERION_STATE)= mFailureCriterionState;
    this->GetValue(CONTACT_RADIUS)         = mContactRadius;*/

    KRATOS_CATCH( "" )
}

void PolyhedronContactElement::InitializeSolutionStep(const ProcessInfo& r_process_info )
{
    //
}

void PolyhedronContactElement::FinalizeSolutionStep(const ProcessInfo& r_process_info) {

}

void PolyhedronContactElement::CalculateRightHandSide(VectorType& rRightHandSideVector,
                                                    const ProcessInfo& rCurrentProcessInfo){
    KRATOS_TRY


    //Poly_particle_1


    //Poly_particle_2




    KRATOS_CATCH( "" )
    }

void PolyhedronContactElement::GJK() {
    KRATOS_TRY

    KRATOS_CATCH( "" )
}

void PolyhedronContactElement::EPA() {
    KRATOS_TRY

    KRATOS_CATCH( "" )
}

void SetId(IndexType NewId) { mId = NewId;}

} // Namespace Kratos


