/////////////////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
// Date: July 2024
/////////////////////////////////////////////////////////////

#include "DEM_P_D_Linear_viscous_Coulomb_CL.h"

namespace Kratos {

    DEMPolyhedronDiscontinuumConstitutiveLaw::Pointer DEM_P_D_Linear_viscous_Coulomb::Clone() const {
        DEMPolyhedronDiscontinuumConstitutiveLaw::Pointer p_clone(new DEM_P_D_Linear_viscous_Coulomb(*this));
        return p_clone;
    }

    std::unique_ptr<DEMPolyhedronDiscontinuumConstitutiveLaw> DEM_P_D_Linear_viscous_Coulomb::CloneUnique() {
        return Kratos::make_unique<DEM_P_D_Linear_viscous_Coulomb>();
    }

    std::string DEM_P_D_Linear_viscous_Coulomb::GetTypeOfLaw() {
        std::string type_of_law = "Linear";
        return type_of_law;
    }

    void DEM_P_D_Linear_viscous_Coulomb::Check(Properties::Pointer pProp) const {
        if(!pProp->Has(STATIC_FRICTION)) {
            if(!pProp->Has(FRICTION)) { //deprecated since April 6th, 2020
                KRATOS_WARNING("DEM")<<std::endl;
                KRATOS_WARNING("DEM")<<"WARNING: Variable STATIC_FRICTION or FRICTION should be present in the properties when using DEMDiscontinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
                KRATOS_WARNING("DEM")<<std::endl;
                pProp->GetValue(STATIC_FRICTION) = 0.0;
            }
            else {
                pProp->GetValue(STATIC_FRICTION) = pProp->GetValue(FRICTION);
            }
        }

        if(!pProp->Has(DYNAMIC_FRICTION)) {
            if(!pProp->Has(FRICTION)) { //deprecated since April 6th, 2020
                KRATOS_WARNING("DEM")<<std::endl;
                KRATOS_WARNING("DEM")<<"WARNING: Variable DYNAMIC_FRICTION or FRICTION should be present in the properties when using DEMDiscontinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
                KRATOS_WARNING("DEM")<<std::endl;
                pProp->GetValue(DYNAMIC_FRICTION) = 0.0;
            }
            else {
                pProp->GetValue(DYNAMIC_FRICTION) = pProp->GetValue(FRICTION);
            }
        }

        if(!pProp->Has(COEFFICIENT_OF_RESTITUTION)) {
            KRATOS_WARNING("DEM")<<std::endl;
            KRATOS_WARNING("DEM")<<"WARNING: Variable COEFFICIENT_OF_RESTITUTION should be present in the properties when using DEMDiscontinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(COEFFICIENT_OF_RESTITUTION) = 0.0;
        }
    }

    /////////////////////////
    // DEM-DEM INTERACTION //
    /////////////////////////

    void DEM_P_D_Linear_viscous_Coulomb::InitializeContact() {

    }

    void DEM_P_D_Linear_viscous_Coulomb::CalculateForces(const ProcessInfo& r_process_info, Vector3 mOverlapVector, Vector3& contact_force) {

        KRATOS_TRY

        double kn = 100000.0;
		contact_force = mOverlapVector * kn;

        KRATOS_CATCH( "" )
    }

} // namespace Kratos
