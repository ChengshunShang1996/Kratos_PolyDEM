/////////////////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
// Date: July 2024
/////////////////////////////////////////////////////////////

#if !defined(DEM_P_D_LINEAR_VISCOUS_COULOMB_CL_H_INCLUDED)
#define DEM_P_D_LINEAR_VISCOUS_COULOMB_CL_H_INCLUDED

#include <string>
#include <iostream>
#include "DEM_polyhedron_discontinuum_constitutive_law.h"

namespace Kratos {

    class SphericParticle;

    class KRATOS_API(DEM_APPLICATION) DEM_P_D_Linear_viscous_Coulomb : public DEMPolyhedronDiscontinuumConstitutiveLaw {

    public:


        KRATOS_CLASS_POINTER_DEFINITION(DEM_P_D_Linear_viscous_Coulomb);

        DEM_P_D_Linear_viscous_Coulomb() {}

        ~DEM_P_D_Linear_viscous_Coulomb() {}

        std::string GetTypeOfLaw() override;

        void Check(Properties::Pointer pProp) const override;

        DEMPolyhedronDiscontinuumConstitutiveLaw::Pointer Clone() const override;

        std::unique_ptr<DEMPolyhedronDiscontinuumConstitutiveLaw> CloneUnique() override;

        void InitializeContact() override;

        void CalculateForces(const ProcessInfo& r_process_info) override;

        Properties& GetPropertiesOfThisContact(SphericParticle* const element, SphericParticle* const neighbour);

    protected:

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, DEMPolyhedronDiscontinuumConstitutiveLaw)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, DEMPolyhedronDiscontinuumConstitutiveLaw)
                    //rSerializer.load("MyMemberName",myMember);
        }

    }; //class DEM_P_D_Linear_viscous_Coulomb

} /* namespace Kratos.*/

#endif /* DEM_P_D_LINEAR_VISCOUS_COULOMB_CL_H_INCLUDED  defined */
