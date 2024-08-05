/////////////////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
// Date: July 2024
/////////////////////////////////////////////////////////////

#if !defined KRATOS_POLYHEDRON_PARTICLE_H_INCLUDED
#define KRATOS_POLYHEDRON_PARTICLE_H_INCLUDED

// System includes
#include <string>

// Project includes
#include "includes/define.h"
#include "includes/node.h"
#include "includes/element.h"
#include "geometries/geometry.h"
#include "includes/properties.h"
#include "includes/indexed_object.h"
#include "containers/global_pointers_vector.h"
//#include "custom_elements/rigid_body_element.h"
#include "custom_elements/spheric_particle.h"
#include "custom_elements/polyhedron_contact_element.h"
#include "custom_utilities/vector3.h"

namespace Kratos
{
    class Element;
    class PolyhedronContactElement;
    class KRATOS_API(DEM_APPLICATION) PolyhedronParticle : public SphericParticle {

    public:
        /// Pointer definition of PolyhedronParticle
        KRATOS_CLASS_INTRUSIVE_POINTER_DEFINITION(PolyhedronParticle);

        PolyhedronParticle();
        PolyhedronParticle(IndexType NewId, GeometryType::Pointer pGeometry);
        PolyhedronParticle(IndexType NewId, NodesArrayType const& ThisNodes);
        PolyhedronParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties);
        Element::Pointer Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const override;

        /// Destructor
        virtual ~PolyhedronParticle();

        void Initialize(const ProcessInfo& r_process_info) override;
        void InitializeSolutionStep(const ProcessInfo& r_process_info) override;
        void ComputeExternalForces(const array_1d<double,3>& gravity);
        virtual void ComputeNewNeighboursHistoricalData(DenseVector<int>& temp_neighbours_ids, std::vector<array_1d<double, 3> >& temp_neighbour_elastic_contact_forces);
        double CalculateVolume() override;
        virtual double GetRadius() override;
        virtual void   SetRadius(double radius) override;
        virtual void   SetRadius() override;
        virtual double GetSearchRadius() override;
        virtual void   SetSearchRadius(const double radius) override;
        double         GetMass() override;
        void           SetMass(double real_mass) override;
        double         GetDensity() override;
        double         SlowGetDensity();
        std::vector<array_1d<double, 3>> GetListOfVertices();
        std::vector<std::vector<int>> GetListOfFaces();
        Vector3 GetFurthestPoint(Vector3 direction);
        void   SetYoungFromProperties(double* young);
        void   SetPoissonFromProperties(double* poisson);
        void   SetDensityFromProperties(double* density);
        void   SetParticleMaterialFromProperties(int* particle_material);

        // 
        double mEnginePower; 

        array_1d<double,3> mDragConstantVector;

        virtual std::string Info() const override
        {
            std::stringstream buffer;
            buffer << "Polyhedron Element #" << Id();
            return buffer.str();
        }

        /// Print information about this object.
        virtual void PrintInfo(std::ostream& rOStream) const override
        {
	        rOStream << "Polyhedron Element #" << Id();
        }

        /// Print object's data.
        virtual void PrintData(std::ostream& rOStream) const override {}

        PropertiesProxy* GetFastProperties();
        void   SetFastProperties(PropertiesProxy* pProps);
        void   SetFastProperties(std::vector<PropertiesProxy>& list_of_proxies);

        std::vector<PolyhedronContactElement*> mPolyhedronContactElements;
        std::vector<PolyhedronParticle*>               mNeighbourElements;
        std::vector<array_1d<double, 3> > mNeighbourRigidFacesTotalContactForce;
        std::vector<array_1d<double, 3> > mNeighbourRigidFacesElasticContactForce;
        std::vector<array_1d<double, 3> > mNeighbourElasticContactForces;
        std::vector<array_1d<double, 3> > mNeighbourElasticExtraContactForces;
        double mRadius;
        double mSearchRadius;
        double mRealMass;
        PropertiesProxy* mFastProperties;
        std::vector<array_1d<double, 3>> mListOfVertices;
        std::vector<std::vector<int>> mListOfFaces;

    protected:


    private:

        friend class Serializer;
        virtual void save(Serializer& rSerializer) const override{ KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, SphericParticle); }
        virtual void load(Serializer& rSerializer) override{ KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, SphericParticle); }

    }; // Class PolyhedronParticle

    /// input stream function
    inline std::istream& operator >> (std::istream& rIStream, PolyhedronParticle& rThis);

    /// output stream function
    inline std::ostream& operator << (std::ostream& rOStream, const PolyhedronParticle& rThis)
    {
        rThis.PrintInfo(rOStream);
        rOStream << std::endl;
        rThis.PrintData(rOStream);
        return rOStream;
    }

}  // namespace Kratos

#endif // KRATOS_POLYHEDRON_PARTICLE_H_INCLUDED defined
