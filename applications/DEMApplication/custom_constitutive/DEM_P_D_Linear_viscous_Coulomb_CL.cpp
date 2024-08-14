/////////////////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
// Date: July 2024
/////////////////////////////////////////////////////////////

#include "DEM_P_D_Linear_viscous_Coulomb_CL.h"
#include "custom_elements/polyhedron_particle.h"

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
                                                        Vector3& contact_force,
                                                        Vector3& TangentialElasticContactForce) {

        KRATOS_TRY

        auto& central_node_1 = PolyhedronParticle1->GetGeometry()[0];
		auto& central_node_2 = PolyhedronParticle2->GetGeometry()[0];

        const array_1d<double, 3>& velocity_1     = central_node_1.FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3>& delta_disp_1  = central_node_1.FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        const array_1d<double, 3>& ang_velocity_1 = central_node_1.FastGetSolutionStepValue(ANGULAR_VELOCITY);
        const array_1d<double, 3>& velocity_2     = central_node_2.FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3>& delta_disp_2  = central_node_2.FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        const array_1d<double, 3>& ang_velocity_2 = central_node_2.FastGetSolutionStepValue(ANGULAR_VELOCITY);

        Vector3 relVel = Vector3(velocity_1[0], velocity_1[1], velocity_1[2]) - 
                         Vector3(velocity_2[0], velocity_2[1], velocity_2[2]);

        Vector3 relDisp = Vector3(delta_disp_1[0], delta_disp_1[1], delta_disp_1[2]) - 
                          Vector3(delta_disp_2[0], delta_disp_2[1], delta_disp_2[2]);

        // Put the values into a more useful form
        Vector3 angVel1(ang_velocity_1[0], ang_velocity_1[1], ang_velocity_1[2]);
        Vector3 angVel2(ang_velocity_2[0], ang_velocity_2[1], ang_velocity_2[2]);

        // The unit vector from element 1 to the contact point
        Vector3 unitCPVect = mOverlapVector;
        unitCPVect.Normalised();

        // normal and tangential components of the relative velocity at the contact point
        Vector3 relVel_n = unitCPVect * Vector3::Dot(unitCPVect, relVel);
        Vector3 relVel_t = relVel - relVel_n;

        Vector3 relDisp_n = unitCPVect * Vector3::Dot(unitCPVect, relDisp);
        Vector3 relDisp_t = relDisp - relDisp_n;

        const double kn = (*mpProperties)[CONTACT_K_N];
        const double kt = (*mpProperties)[CONTACT_K_N] * 0.5;
        const double static_friction = (*mpProperties)[STATIC_FRICTION];
        const double dynamic_friction = (*mpProperties)[DYNAMIC_FRICTION];
        const double equiv_friction_decay_coefficient = (*mpProperties)[FRICTION_DECAY];

        const double my_mass    = PolyhedronParticle1->GetMass();
        const double other_mass = PolyhedronParticle2->GetMass();

        const double equiv_mass = 1.0 / (1.0/my_mass + 1.0/other_mass);

        Properties& properties_of_this_contact = PolyhedronParticle1->GetProperties().GetSubProperties(PolyhedronParticle2->GetProperties().Id());
        const double damping_gamma = properties_of_this_contact[DAMPING_GAMMA];

        const double equiv_visco_damp_coeff_normal     = 2.0 * damping_gamma * sqrt(equiv_mass * kn);
        const double equiv_visco_damp_coeff_tangential = 2.0 * damping_gamma * sqrt(equiv_mass * kt);

        // Damping calculation
        Vector3 F_n = mOverlapVector * kn;
        Vector3 F_nd = -unitCPVect * equiv_visco_damp_coeff_normal *  relVel_n.Length();

        // Are we in a loading situation?
        if (Vector3::Dot(relVel_n, unitCPVect) > 0.0)
        {
            F_nd = -F_nd;
        }

        TangentialElasticContactForce -= relDisp_t * kt;
        Vector3 F_t = TangentialElasticContactForce;

        // Damping
        Vector3 F_td, newF_t;

        if (F_t.Length() > F_n.Length() * (static_friction))
        {
            newF_t = F_t * F_n.Length() * (dynamic_friction) / F_t.Length();
            TangentialElasticContactForce = newF_t;

            //at this point we get energy loss from the sliding!
            F_td = Vector3(0.0, 0.0, 0.0);
        }
        else
        {
            //at this point we get energy loss from the damping!
            F_td = -relVel_t * equiv_visco_damp_coeff_tangential;
            newF_t = F_t + F_td;
        }

        /*
        double equiv_friction = dynamic_friction + (static_friction - dynamic_friction) * exp(-equiv_friction_decay_coefficient * relVel_t.Length());
        double maximum_admissible_shear_force = F_n.Length() * equiv_friction;

        F_td = -relVel_t * equiv_visco_damp_coeff_tangential;

        Vector3 tangential_contact_force = F_t + F_td;

        const double ActualTotalShearForce = tangential_contact_force.Length();

        if (ActualTotalShearForce > maximum_admissible_shear_force) {

            const double ActualElasticShearForce = F_t.Length();

            const double dot_product = Vector3::Dot(F_t, F_td);
            const double ViscoDampingLocalContactForceModule = F_td.Length();

            if (dot_product >= 0.0) {

                if (ActualElasticShearForce > maximum_admissible_shear_force) {
                    const double fraction = maximum_admissible_shear_force / ActualElasticShearForce;
                    F_t *= fraction;
                    F_td = Vector3(0.0, 0.0, 0.0);
                }
                else {
                    const double ActualViscousShearForce = maximum_admissible_shear_force - ActualElasticShearForce;
                    const double fraction = ActualViscousShearForce / ViscoDampingLocalContactForceModule;
                    F_td *= fraction;
                }
            }
            else {
                if (ViscoDampingLocalContactForceModule >= ActualElasticShearForce) {
                    const double fraction = (maximum_admissible_shear_force + ActualElasticShearForce) / ViscoDampingLocalContactForceModule;
                    F_td *= fraction;
                }
                else {
                    const double fraction = maximum_admissible_shear_force / ActualElasticShearForce;
                    F_t *= fraction;
                    F_td = Vector3(0.0, 0.0, 0.0);
                }
            }
        }*/

        //double kn = 100000.0;
		//contact_force = mOverlapVector * kn;
        //contact_force = F_n + F_nd + F_t + F_td;
        contact_force = F_n + F_nd + newF_t;

        KRATOS_CATCH( "" )
    }

} // namespace Kratos
