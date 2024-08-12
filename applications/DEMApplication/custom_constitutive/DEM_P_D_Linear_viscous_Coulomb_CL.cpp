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
            KRATOS_WARNING("DEM")<<"WARNING: Variable COEFFICIENT_OF_RESTITUTION should be present in the properties when using DEMPolyhedronDiscontinuumConstitutiveLaw. 0.0 value assigned by default."<<std::endl;
            KRATOS_WARNING("DEM")<<std::endl;
            pProp->GetValue(COEFFICIENT_OF_RESTITUTION) = 0.0;
        }

        if(!pProp->Has(CONTACT_K_N)) {

            KRATOS_ERROR << "Variable CONTACT_K_N should be present in the properties when using DEMPolyhedronDiscontinuumConstitutiveLaw."<<std::endl;

        }
    }

    /////////////////////////
    // DEM-DEM INTERACTION //
    /////////////////////////

    void DEM_P_D_Linear_viscous_Coulomb::InitializeContact() {

    }

    void DEM_P_D_Linear_viscous_Coulomb::CalculateForces(const ProcessInfo& r_process_info, 
                                                        PolyhedronParticle* PolyhedronParticle1, 
                                                        PolyhedronParticle* PolyhedronParticle2, 
                                                        Vector3 mOverlapVector, 
                                                        Vector3& contact_force) {

        KRATOS_TRY

        auto& central_node_1 = PolyhedronParticle1->GetGeometry()[0];
		auto& central_node_2 = PolyhedronParticle2->GetGeometry()[0];

        const array_1d<double, 3>& velocity_1     = central_node_1.FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3>& delta_displ_1  = central_node_1.FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        const array_1d<double, 3>& ang_velocity_1 = central_node_1.FastGetSolutionStepValue(ANGULAR_VELOCITY);
        const array_1d<double, 3>& velocity_2     = central_node_2.FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3>& delta_displ_2  = central_node_2.FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        const array_1d<double, 3>& ang_velocity_2 = central_node_2.FastGetSolutionStepValue(ANGULAR_VELOCITY);

        Vector3 relVel = Vector3(velocity_1[0], velocity_1[1], velocity_1[2]) - 
                         Vector3(velocity_2[0], velocity_2[1], velocity_2[2]);

        // Put the values into a more useful form
        Vector3 angVel1(ang_velocity_1[0], ang_velocity_1[1], ang_velocity_1[2]);
        Vector3 angVel2(ang_velocity_2[0], ang_velocity_2[1], ang_velocity_2[2]);

        // The unit vector from element 1 to the contact point
        Vector3 unitCPVect = mOverlapVector;
        unitCPVect.Normalised();

        // normal and tangential components of the relative velocity at the contact point
        Vector3 relVel_n = unitCPVect * Vector3::Dot(unitCPVect, relVel);
        Vector3 relVel_t = relVel - relVel_n;

        // Damping calculation
        const double kn = (*mpProperties)[CONTACT_K_N];
        Vector3 F_nd = unitCPVect * 2 * std::sqrt(5.0 / 6.0) * 0.2 * std::sqrt(kn * 1) *  relVel_n.length();

        // Are we in a loading situation?
        if (Vector3::Dot(relVel_n, unitCPVect) > 0.0)
        {
            F_nd = -F_nd;
        }

        const double kt = (*mpProperties)[CONTACT_K_N];
        Vector3 nOverlap_t(tangentialPhysicalOverlapX, tangentialPhysicalOverlapY, tangentialPhysicalOverlapZ);
        Vector3 F_t = -nOverlap_t * kt;

        // Damping
        Vector3 F_td, newF_t;

        if (F_t.length() > F_n.length() * (staticFriction))
        {
            newF_t = F_t * F_n.length() * (staticFriction) / F_t.length();
            nOverlap_t = -newF_t / kt; //slippage has occurred so the tangential overlap is reduced a bit

            //at this point we get energy loss from the sliding!
            F_td = newF_t;
        }
        else
        {
            //at this point we get energy loss from the damping!
            F_td = -relVel_t * 2 * sqrt(5.0 / 6.0) * 0.2 * std::sqrt(kt * nEquivMass);
            newF_t = F_t + F_td;
        }

        //double kn = 100000.0;
		contact_force = mOverlapVector * kn;

        KRATOS_CATCH( "" )
    }

} // namespace Kratos
