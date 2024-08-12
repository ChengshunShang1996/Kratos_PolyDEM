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

    void DEM_P_D_Linear_viscous_Coulomb::CalculateForces(const ProcessInfo& r_process_info, Vector3 mOverlapVector, Vector3& contact_force) {

        KRATOS_TRY

        Vector3 relVel = Vector3(elem1ContactPointVelX, elem1ContactPointVelY, elem1ContactPointVelZ) -
        Vector3(elem2ContactPointVelX, elem2ContactPointVelY, elem2ContactPointVelZ);

        // Put the values into a more useful form
        Vector3 angVel1(elem1AngVelX, elem1AngVelY, elem1AngVelZ);
        Vector3 angVel2(elem2AngVelX, elem2AngVelY, elem2AngVelZ);
        CSimple3DPoint  contactPoint(contactPointX, contactPointY, contactPointZ);

        // The unit vector from element 1 to the contact point
        Vector3 unitCPVect = contactPoint - CSimple3DPoint(elem1PosX, elem1PosY, elem1PosZ);
        unitCPVect.normalise();

        // normal and tangential components of the relative velocity at the contact point
        Vector3 relVel_n = unitCPVect * unitCPVect.dot(relVel);
        Vector3 relVel_t = relVel - relVel_n;

        // Damping calculation
        double B = 0.0;
        if (coeffRest > 0.0)
        {
            double myLog = log(coeffRest);
            B = -myLog / sqrt(myLog * myLog + PI * PI);
        }

        double S_n = 2.0 * nEquivYoungsMod * sqrt(nEquivRadius * normalPhysicalOverlap);
        Vector3 F_nd = unitCPVect * 2 * sqrt(5.0 / 6.0) * B * sqrt(S_n * nEquivMass) *  relVel_n.length();

        // Are we in a loading situation?
        if (relVel_n.dot(unitCPVect) > 0.0)
        {
            F_nd = -F_nd;
        }

        double S_t = 8.0 * nEquivShearMod * sqrt(nEquivRadius * normalPhysicalOverlap);
        Vector3 nOverlap_t(tangentialPhysicalOverlapX, tangentialPhysicalOverlapY, tangentialPhysicalOverlapZ);
        Vector3 F_t = -nOverlap_t * S_t;

        // Damping
        Vector3 F_td, newF_t;

        if (F_t.length() > F_n.length() * (staticFriction))
        {
            newF_t = F_t * F_n.length() * (staticFriction) / F_t.length();
            nOverlap_t = -newF_t / S_t; //slippage has occurred so the tangential overlap is reduced a bit

            //at this point we get energy loss from the sliding!
            F_td = newF_t;
        }
        else
        {
            //at this point we get energy loss from the damping!
            F_td = -relVel_t * 2 * sqrt(5.0 / 6.0) * B * sqrt(S_t * nEquivMass);
            newF_t = F_t + F_td;
        }

        //double kn = 100000.0;
        const double kn = (*mpProperties)[CONTACT_K_N];
		contact_force = mOverlapVector * kn;

        KRATOS_CATCH( "" )
    }

} // namespace Kratos
