/////////////////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
// Date: July 2024
/////////////////////////////////////////////////////////////

#if !defined(KRATOS_POLYHEDRON_CONTACT_ELEMENT_H_INCLUDED)
#define  KRATOS_POLYHEDRON_CONTACT_ELEMENT_H_INCLUDED

// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "includes/serializer.h"
#include "includes/element.h"
#include "discrete_element.h"
#include "includes/ublas_interface.h"
#include "includes/variables.h"

namespace Kratos
{

    class PolyhedronParticle;
    
    class KRATOS_API(DEM_APPLICATION) PolyhedronContactElement {

    public:

        using Pointer = PolyhedronContactElement*;

        typedef Properties PropertiesType;
        typedef Vector VectorType;
        
        //KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(PolyhedronContactElement);

        PolyhedronContactElement() {}
        PolyhedronContactElement(std::string ElementName) {}
        PolyhedronContactElement(IndexType NewId, PolyhedronParticle* PolyhedronParticle1, PolyhedronParticle* PolyhedronParticle2);
        PolyhedronContactElement(IndexType NewId, PolyhedronParticle* PolyhedronParticle1, PolyhedronParticle* PolyhedronParticle2, PropertiesType::Pointer pProperties);

        virtual ~PolyhedronContactElement();

        virtual Pointer Create(IndexType NewId, PolyhedronParticle* PolyhedronParticle1, PolyhedronParticle* PolyhedronParticle2, PropertiesType::Pointer pProperties) const;

        virtual void Initialize(const ProcessInfo& r_process_info);

        virtual void InitializeSolutionStep(const ProcessInfo& r_process_info);

        virtual void CalculateRightHandSide(const ProcessInfo& r_process_info, double dt, const array_1d<double, 3>& gravity);

        virtual void GJK();

        virtual void EPA();

        virtual void FinalizeSolutionStep(const ProcessInfo& r_process_info) ;

        virtual void PrepareForPrinting();

        void SetId(IndexType NewId);

        array_1d<double, 3> mContactForce;
        array_1d<double, 3> mRotationalMoment;
        double mFailureState;
        double mContactRadius;
        IndexType mId;
        PolyhedronParticle* mPolyhedronParticle1;
        PolyhedronParticle* mPolyhedronParticle2;
        Properties::Pointer mpProperties;

    protected:

    private:


    }; // Class PolyhedronContactElement

}  // namespace Kratos.

#endif // KRATOS_POLYHEDRON_CONTACT_ELEMENT_H_INCLUDED  defined
