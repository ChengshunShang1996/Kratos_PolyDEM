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
#include "custom_elements/rigid_body_element.h"
#include "custom_elements/polyhedron_contact_element.h"

namespace Kratos
{
    class Element;
    class KRATOS_API(DEM_APPLICATION) PolyhedronParticle : public RigidBodyElement3D {

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
        void CustomInitialize(ModelPart& rigid_body_element_sub_model_part) override;
        void ComputeExternalForces(const array_1d<double,3>& gravity) override;
        virtual void ComputeNewNeighboursHistoricalData(DenseVector<int>& temp_neighbours_ids, std::vector<array_1d<double, 3> >& temp_neighbour_elastic_contact_forces);
        double CalculateVolume();
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
        PropertiesProxy* mFastProperties;


    protected:


    private:

        friend class Serializer;
        virtual void save(Serializer& rSerializer) const override{ KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, RigidBodyElement3D); }
        virtual void load(Serializer& rSerializer) override{ KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, RigidBodyElement3D); }

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
