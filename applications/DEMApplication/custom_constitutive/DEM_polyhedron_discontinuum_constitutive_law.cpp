/////////////////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
// Date: July 2024
/////////////////////////////////////////////////////////////

#include "DEM_polyhedron_discontinuum_constitutive_law.h"

namespace Kratos {

    DEMPolyhedronDiscontinuumConstitutiveLaw::DEMPolyhedronDiscontinuumConstitutiveLaw() {}

    //copy constructor
    DEMPolyhedronDiscontinuumConstitutiveLaw::DEMPolyhedronDiscontinuumConstitutiveLaw(const DEMPolyhedronDiscontinuumConstitutiveLaw &rReferenceDiscontinuumConstitutiveLaw) {}

    void DEMPolyhedronDiscontinuumConstitutiveLaw::SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose) {
        if (verbose) KRATOS_INFO("DEM")  << "Assigning " << pProp->GetValue(DEM_POLYHEDRON_DISCONTINUUM_CONSTITUTIVE_LAW_NAME) << " to Properties " << pProp->Id() << std::endl;
        pProp->SetValue(DEM_POLYHEDRON_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER, this->Clone());
        this->Check(pProp);
    }

    std::unique_ptr<DEMPolyhedronDiscontinuumConstitutiveLaw> DEMPolyhedronDiscontinuumConstitutiveLaw::CloneUnique() {
        KRATOS_ERROR << "This function (DEMDiscontinuumConstitutiveLaw::CloneUnique) shouldn't be accessed, use derived class instead"<<std::endl;
    }

    DEMPolyhedronDiscontinuumConstitutiveLaw::Pointer DEMPolyhedronDiscontinuumConstitutiveLaw::Clone() const {
        DEMPolyhedronDiscontinuumConstitutiveLaw::Pointer p_clone(new DEMPolyhedronDiscontinuumConstitutiveLaw(*this));
        return p_clone;
    }

    DEMPolyhedronDiscontinuumConstitutiveLaw::~DEMPolyhedronDiscontinuumConstitutiveLaw() {}

    void DEMPolyhedronDiscontinuumConstitutiveLaw::Check(Properties::Pointer pProp) const {
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::Check) shouldn't be accessed, use derived class instead"<<std::endl;
    }

    void DEMPolyhedronDiscontinuumConstitutiveLaw::Initialize(Properties::Pointer pProps) {
        mpProperties = pProps;
    }

    void DEMPolyhedronDiscontinuumConstitutiveLaw::InitializeContact() {
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::InitializeContact) shouldn't be accessed, use derived class instead"<<std::endl;
    }

    std::string DEMPolyhedronDiscontinuumConstitutiveLaw::GetTypeOfLaw() {
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::GetTypeOfLaw) shouldn't be accessed, use derived class instead"<<std::endl;
        std::string type_of_law = "";
        return type_of_law;
    }

    void DEMPolyhedronDiscontinuumConstitutiveLaw::CalculateForces(const ProcessInfo& r_process_info, 
                                                                    PolyhedronParticle* PolyhedronParticle1, 
                                                                    PolyhedronParticle* PolyhedronParticle2, 
                                                                    Vector3 mOverlapVector, 
                                                                    Vector3& contact_force) {
        KRATOS_ERROR << "This function (DEMContinuumConstitutiveLaw::CalculateForces) shouldn't be accessed, use derived class instead"<<std::endl;
    }

} // KRATOS
