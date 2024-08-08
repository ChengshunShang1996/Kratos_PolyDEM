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
#include "custom_utilities/GeometryFunctions.h"
#include "custom_elements/polyhedron_particle.h"
#include "custom_utilities/vector3.h"

namespace Kratos
{

    class PolyhedronParticle;
    
    class PolyPlane {
    public:
        
        Vector3 normal;
        double d;

        PolyPlane(const Vector3& normal, double d) : normal(normal), d(d) {}

        static PolyPlane PlaneFromTri(const Vector3& p0, const Vector3& p1, const Vector3& p2) {
            Vector3 normal = (Vector3::Cross(p1 - p0, p2 - p0)).Normalised();
            return PolyPlane(normal, -Vector3::Dot(p0, normal));
        }

        Vector3 ProjectPointOntoPlane(const Vector3& point) const {
            double distance = Vector3::Dot(point, normal) + d;
            return point - (normal * distance);
        }
    };
    
    class KRATOS_API(DEM_APPLICATION) PolyhedronContactElement {

    public:

        using Pointer = PolyhedronContactElement*;

        typedef Properties PropertiesType;
        typedef Vector VectorType;
        struct Point {
                Vector3 p; //Conserve Minkowski Difference
                Vector3 a; //Result coordinate of object A's support function 
                Vector3 b; //Result coordinate of object B's support function 
            };

        #define GJK_MAX_NUM_ITERATIONS 100
        #define EPA_TOLERANCE 1e-6
        #define EPA_MAX_NUM_FACES 64
        #define EPA_MAX_NUM_LOOSE_EDGES 32
        #define EPA_MAX_NUM_ITERATIONS 100
        
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

        bool GJK();
        void update_simplex3(Point& a, Point& b, Point& c, Point& d, int& simp_dim, Vector3& search_dir);
        bool update_simplex4(Point& a, Point& b, Point& c, Point& d, int& simp_dim, Vector3& search_dir);
        //Expanding Polytope Algorithm. 
        void EPA(Point& a, Point& b, Point& c, Point& d);
        void CalculateSearchPoint(Point& point, Vector3& search_dir);
        void Barycentric(const Vector3 & a, const Vector3 & b, const Vector3 & c, const Vector3 & p, double& u, double& v, double& w);

        virtual void FinalizeSolutionStep(const ProcessInfo& r_process_info) ;

        virtual void PrepareForPrinting();

        void SetId(IndexType NewId);
        virtual void SetPolyElement1(PolyhedronParticle* custom_poly_element);
        virtual void SetPolyElement2(PolyhedronParticle* custom_poly_element);
        virtual PolyhedronParticle* GetPolyElement1();
        virtual PolyhedronParticle* GetPolyElement2();
        void SetDeleteFlag(bool this_flag);
        bool GetDeleteFlag();

        array_1d<double, 3> mContactForce;
        array_1d<double, 3> mRotationalMoment;
        double mFailureState;
        double mContactRadius;
        IndexType mId;
        PolyhedronParticle* mPolyhedronParticle1;
        PolyhedronParticle* mPolyhedronParticle2;
        Properties::Pointer mpProperties;
        Vector3 mOverlapVector;
        Vector3 mContactPoint1;
        Vector3 mContactPoint2;
        bool mDeleteFlag;

    protected:

    private:


    }; // Class PolyhedronContactElement

}  // namespace Kratos.

#endif // KRATOS_POLYHEDRON_CONTACT_ELEMENT_H_INCLUDED  defined
