/////////////////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
// Date: July 2024
/////////////////////////////////////////////////////////////

#if !defined(DEM_POLYHEDRON_DISCONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED)
#define  DEM_POLYHEDRON_DISCONTINUUM_CONSTITUTIVE_LAW_H_INCLUDED

/* Project includes */
#include "includes/define.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "includes/serializer.h"
#include "containers/flags.h"

#include "custom_elements/discrete_element.h"
#include "custom_utilities/vector3.h"


namespace Kratos {

    class Properties;
    class PolyhedronParticle; // forward declaration of polyhedron particle

    class KRATOS_API(DEM_APPLICATION) DEMPolyhedronDiscontinuumConstitutiveLaw : public Flags {
    public:

        KRATOS_CLASS_POINTER_DEFINITION(DEMPolyhedronDiscontinuumConstitutiveLaw);

        DEMPolyhedronDiscontinuumConstitutiveLaw();

        DEMPolyhedronDiscontinuumConstitutiveLaw(const DEMPolyhedronDiscontinuumConstitutiveLaw &rReferenceDiscontinuumConstitutiveLaw);

        virtual void SetConstitutiveLawInProperties(Properties::Pointer pProp, bool verbose = true);

        virtual void Check(Properties::Pointer pProp) const;

        virtual ~DEMPolyhedronDiscontinuumConstitutiveLaw();

        virtual DEMPolyhedronDiscontinuumConstitutiveLaw::Pointer Clone() const;

        virtual std::unique_ptr<DEMPolyhedronDiscontinuumConstitutiveLaw> CloneUnique();

        virtual std::string GetTypeOfLaw();

        virtual void Initialize(Properties::Pointer pProps);
        
        virtual void InitializeContact();

        virtual void CalculateForces(const ProcessInfo& r_process_info, 
                                    PolyhedronParticle* PolyhedronParticle1, 
                                    PolyhedronParticle* PolyhedronParticle2, 
                                    Vector3 OverlapVector, 
                                    Vector3 ContactPoint,
                                    Vector3& contact_force,
                                    Vector3& TangentialElasticContactForce);

    protected:

        Properties::Pointer mpProperties;

    private:

        friend class Serializer;

        virtual void save(Serializer& rSerializer) const override {
            KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, Flags)
                    //rSerializer.save("MyMemberName",myMember);
        }

        virtual void load(Serializer& rSerializer) override {
            KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, Flags)
                    //rSerializer.load("MyMemberName",myMember);
        }
    };

    //This definition is done here to avoid recursive inclusion of header files
    KRATOS_DEFINE_APPLICATION_VARIABLE(DEM_APPLICATION, DEMPolyhedronDiscontinuumConstitutiveLaw::Pointer, DEM_POLYHEDRON_DISCONTINUUM_CONSTITUTIVE_LAW_POINTER)

} /* namespace Kratos.*/

#endif /* DEM_POLYHEDRON_CONSTITUTIVE_LAW_H_INCLUDED  defined */

