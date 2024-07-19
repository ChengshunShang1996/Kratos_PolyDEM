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
#include "includes/ublas_interface.h"
#include "includes/variables.h"

namespace Kratos
{

    class PolyhedronParticle;
    
    class KRATOS_API(DEM_APPLICATION) PolyhedronContactElement {

    public:

        using Pointer = PolyhedronContactElement*;
        
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(PolyhedronContactElement);

        PolyhedronContactElement(IndexType NewId, GeometryType::Pointer pGeometry);
        PolyhedronContactElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties);

        virtual ~PolyhedronContactElement();

        virtual Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const;

        void Initialize(const ProcessInfo& r_process_info) override;

        void InitializeSolutionStep(const ProcessInfo& r_process_info) override;

        void CalculateRightHandSide(VectorType& rRightHandSideVector, const ProcessInfo& rCurrentProcessInfo) override;

        void GJK();

        void EPA();

        void FinalizeSolutionStep(const ProcessInfo& r_process_info) override;

        void PrepareForPrinting();

        array_1d<double,3> mLocalContactForce;
        array_1d<double,3> mElasticLocalRotationalMoment;
        double mContactSigma;
        double mContactTau;
        PolyhedronParticle* mPolyhedronParticle1;
        PolyhedronParticle* mPolyhedronParticle2;

    protected:


    private:


    }; // Class PolyhedronContactElement

}  // namespace Kratos.

#endif // KRATOS_POLYHEDRON_CONTACT_ELEMENT_H_INCLUDED  defined
