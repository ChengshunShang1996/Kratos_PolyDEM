/////////////////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
// Date: July 2024
/////////////////////////////////////////////////////////////

// System includes
#include <string>
#include <iostream>
#include <cstdlib>

// Project includes
#include "polyhedron_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/variables.h"

namespace Kratos {

    PolyhedronParticle::PolyhedronParticle() : RigidBodyElement3D() {}

    PolyhedronParticle::PolyhedronParticle(IndexType NewId, GeometryType::Pointer pGeometry)
    : RigidBodyElement3D(NewId, pGeometry) {
        mRadius = 0;
        mRealMass = 0;
    }

    PolyhedronParticle::PolyhedronParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : RigidBodyElement3D(NewId, pGeometry, pProperties) {
        mRadius = 0;
        mRealMass = 0;
    }

    PolyhedronParticle::PolyhedronParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : RigidBodyElement3D(NewId, ThisNodes) {
        mRadius = 0;
        mRealMass = 0;
    }

    Element::Pointer PolyhedronParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
        return Element::Pointer(new PolyhedronParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    PolyhedronParticle::~PolyhedronParticle() {}

    void PolyhedronParticle::Initialize(const ProcessInfo& r_process_info) {

        KRATOS_TRY

        RigidBodyElement3D::Initialize(r_process_info);

        auto& central_node = GetGeometry()[0];

        SetRadius();

	    KRATOS_ERROR_IF_NOT(GetProperties().Has(POLYHEDRON_INFORMATION))<<"Something went wrong. Properties do not contain POLYHEDRON_INFORMATION.";
        const PolyhedronInformation& poly_info = GetProperties()[POLYHEDRON_INFORMATION];
        const double reference_size = poly_info.mSize;
        const double reference_volume = poly_info.mVolume;
        const std::vector<std::vector<int>>& reference_list_of_faces = poly_info.mListOfFaces;
        const std::vector<array_1d<double,3> >& reference_list_of_vertices = poly_info.mListOfVertices;

        const unsigned int number_of_vertices = reference_list_of_vertices.size();

        mListOfVertices.resize(number_of_vertices);

        const double scaling_factor = (mRadius * 2.0) / reference_size;

        for (int i = 0; i < (int)number_of_vertices; i++) {
            mListOfVertices[i][0] = scaling_factor * reference_list_of_vertices[i][0];
            mListOfVertices[i][1] = scaling_factor * reference_list_of_vertices[i][1];
            mListOfVertices[i][2] = scaling_factor * reference_list_of_vertices[i][2];
        }

        const double particle_density = this->SlowGetDensity();
        const double polyhedron_volume = reference_volume * scaling_factor * scaling_factor * scaling_factor;
        const double polyhedron_mass = particle_density * polyhedron_volume;

        central_node.FastGetSolutionStepValue(NODAL_MASS) = polyhedron_mass;
        //central_node.FastGetSolutionStepValue(POLYHEDRON_VOLUME) = polyhedron_volume;

        //SetMass(GetDensity() * CalculateVolume());
        SetMass(polyhedron_mass);

        KRATOS_CATCH("")
    }

    void PolyhedronParticle::InitializeSolutionStep(const ProcessInfo& r_process_info)
    {
        KRATOS_TRY

        auto& central_node = GetGeometry()[0];
        mRadius = central_node.FastGetSolutionStepValue(RADIUS); //Just in case someone is overwriting the radius in Python
        central_node.FastGetSolutionStepValue(REPRESENTATIVE_VOLUME) = CalculateVolume();

        KRATOS_CATCH("")
    }
    
    void PolyhedronParticle::ComputeExternalForces(const array_1d<double,3>& gravity) {

        KRATOS_TRY

        noalias(GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES)) += GetMass() * gravity;
        
        const array_1d<double, 3> external_applied_torque = GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT);
        noalias(GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT)) += external_applied_torque;

        KRATOS_CATCH("")
    }

    void PolyhedronParticle::ComputeNewNeighboursHistoricalData(DenseVector<int>& temp_neighbours_ids,
                                                         std::vector<array_1d<double, 3> >& temp_neighbour_elastic_contact_forces)
    {
        KRATOS_TRY

        /*
        std::vector<array_1d<double, 3> > temp_neighbour_elastic_extra_contact_forces;
        unsigned int new_size = mNeighbourElements.size();
        array_1d<double, 3> vector_of_zeros = ZeroVector(3);
        temp_neighbours_ids.resize(new_size, false);
        temp_neighbour_elastic_contact_forces.resize(new_size);
        temp_neighbour_elastic_extra_contact_forces.resize(new_size);

        DenseVector<int>& vector_of_ids_of_neighbours = GetValue(NEIGHBOUR_IDS);

        for (unsigned int i = 0; i < new_size; i++) {
            noalias(temp_neighbour_elastic_contact_forces[i]) = vector_of_zeros;
            noalias(temp_neighbour_elastic_extra_contact_forces[i]) = vector_of_zeros;

            if (mNeighbourElements[i] == NULL) { // This is required by the continuum sphere which reorders the neighbors
                temp_neighbours_ids[i] = -1;
                continue;
            }

            temp_neighbours_ids[i] = mNeighbourElements[i]->Id();

            for (unsigned int j = 0; j < vector_of_ids_of_neighbours.size(); j++) {
                if (int(temp_neighbours_ids[i]) == vector_of_ids_of_neighbours[j] && vector_of_ids_of_neighbours[j] != -1) {
                    noalias(temp_neighbour_elastic_contact_forces[i]) = mNeighbourElasticContactForces[j];
                    noalias(temp_neighbour_elastic_extra_contact_forces[i]) = mNeighbourElasticExtraContactForces[j]; //TODO: remove this from discontinuum!!
                    break;
                }
            }
        }

        vector_of_ids_of_neighbours.swap(temp_neighbours_ids);
        mNeighbourElasticContactForces.swap(temp_neighbour_elastic_contact_forces);
        mNeighbourElasticExtraContactForces.swap(temp_neighbour_elastic_extra_contact_forces);
        */

        KRATOS_CATCH("")
    }

    double PolyhedronParticle::GetMass()                                                         { return mRealMass;       }
    void   PolyhedronParticle::SetMass(double real_mass)                                         { mRealMass = real_mass;  GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = real_mass;}
    //TODO: need to be updated
    double PolyhedronParticle::CalculateVolume()                                                 { return 4.0 * Globals::Pi / 3.0 * mRadius * mRadius * mRadius; }
    double PolyhedronParticle::GetRadius()                                                       { return mRadius;         }
    void   PolyhedronParticle::SetRadius(double radius)                                          { mRadius = radius;       }
    void   PolyhedronParticle::SetRadius()                                                       { mRadius = GetGeometry()[0].FastGetSolutionStepValue(RADIUS); }
    double PolyhedronParticle::GetSearchRadius()                                                 { return mSearchRadius;   }
    void   PolyhedronParticle::SetSearchRadius(const double radius)                              { mSearchRadius = radius; }
    double PolyhedronParticle::GetDensity()                                                      { return GetFastProperties()->GetDensity();}
    double PolyhedronParticle::SlowGetDensity()                                                  { return GetProperties()[PARTICLE_DENSITY];}
    std::vector<array_1d<double, 3>> PolyhedronParticle::GetListOfVertices()                     { return mListOfVertices;}

    PropertiesProxy* PolyhedronParticle::GetFastProperties()                                     { return mFastProperties;   }
    void   PolyhedronParticle::SetFastProperties(PropertiesProxy* pProps)                        { mFastProperties = pProps; }
    void   PolyhedronParticle::SetFastProperties(std::vector<PropertiesProxy>& list_of_proxies)  {
        for (unsigned int j = 0; j < list_of_proxies.size(); j++){
            if (list_of_proxies[j].GetId() == GetProperties().Id()) {
                SetFastProperties(&list_of_proxies[j]);
                return;
            }
        }
    }

}  // namespace Kratos
