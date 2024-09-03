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

    PolyhedronParticle::PolyhedronParticle() : SphericParticle() {}

    PolyhedronParticle::PolyhedronParticle(IndexType NewId, GeometryType::Pointer pGeometry)
    : SphericParticle(NewId, pGeometry) {
        mRadius = 0;
        mRealMass = 0;
    }

    PolyhedronParticle::PolyhedronParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : SphericParticle(NewId, pGeometry, pProperties) {
        mRadius = 0;
        mRealMass = 0;
    }

    PolyhedronParticle::PolyhedronParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : SphericParticle(NewId, ThisNodes) {
        mRadius = 0;
        mRealMass = 0;
    }

    Element::Pointer PolyhedronParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
        return SphericParticle::Pointer(new PolyhedronParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    PolyhedronParticle::~PolyhedronParticle() {}

    void PolyhedronParticle::Initialize(const ProcessInfo& r_process_info) {

        KRATOS_TRY

        auto& central_node = GetGeometry()[0];

        if (central_node.GetDof(VELOCITY_X).IsFixed())          central_node.Set(DEMFlags::FIXED_VEL_X, true);
        else                                                        central_node.Set(DEMFlags::FIXED_VEL_X, false);
        if (central_node.GetDof(VELOCITY_Y).IsFixed())          central_node.Set(DEMFlags::FIXED_VEL_Y, true);
        else                                                        central_node.Set(DEMFlags::FIXED_VEL_Y, false);
        if (central_node.GetDof(VELOCITY_Z).IsFixed())          central_node.Set(DEMFlags::FIXED_VEL_Z, true);
        else                                                        central_node.Set(DEMFlags::FIXED_VEL_Z, false);
        if (central_node.GetDof(ANGULAR_VELOCITY_X).IsFixed())  central_node.Set(DEMFlags::FIXED_ANG_VEL_X, true);
        else                                                        central_node.Set(DEMFlags::FIXED_ANG_VEL_X, false);
        if (central_node.GetDof(ANGULAR_VELOCITY_Y).IsFixed())  central_node.Set(DEMFlags::FIXED_ANG_VEL_Y, true);
        else                                                        central_node.Set(DEMFlags::FIXED_ANG_VEL_Y, false);
        if (central_node.GetDof(ANGULAR_VELOCITY_Z).IsFixed())  central_node.Set(DEMFlags::FIXED_ANG_VEL_Z, true);
        else                                                        central_node.Set(DEMFlags::FIXED_ANG_VEL_Z, false);

        DEMIntegrationScheme::Pointer& translational_integration_scheme = GetProperties()[DEM_TRANSLATIONAL_INTEGRATION_SCHEME_POINTER];
        DEMIntegrationScheme::Pointer& rotational_integration_scheme = GetProperties()[DEM_ROTATIONAL_INTEGRATION_SCHEME_POINTER];
        SetIntegrationScheme(translational_integration_scheme, rotational_integration_scheme);
        //SphericParticle::Initialize(r_process_info);

        //auto& central_node = GetGeometry()[0];

        ///central_node.GetSolutionStepValue(PARTICLE_MATERIAL) = GetParticleMaterial();

        SetRadius();

	    KRATOS_ERROR_IF_NOT(GetProperties().Has(POLYHEDRON_INFORMATION))<<"Something went wrong. Properties do not contain POLYHEDRON_INFORMATION.";
        const PolyhedronInformation& poly_info = GetProperties()[POLYHEDRON_INFORMATION];
        const int poly_shape_index = central_node.GetSolutionStepValue(POLYHEDRON_SHAPE_INDEX);
        const std::vector<array_1d<double,3> >& reference_list_of_vertices = poly_info.mListOfVerticesList[poly_shape_index];
        const std::vector<std::vector<int>>& reference_list_of_faces = poly_info.mListOfFacesList[poly_shape_index];
        const double reference_size = poly_info.mListOfSize[poly_shape_index];
        const double reference_volume = poly_info.mListOfVolume[poly_shape_index];

        const unsigned int number_of_vertices = reference_list_of_vertices.size();

        mListOfVertices.resize(number_of_vertices);

        const double scaling_factor = (mRadius * 2.0) / reference_size;

        for (int i = 0; i < (int)number_of_vertices; i++) {
            mListOfVertices[i][0] = scaling_factor * reference_list_of_vertices[i][0];
            mListOfVertices[i][1] = scaling_factor * reference_list_of_vertices[i][1];
            mListOfVertices[i][2] = scaling_factor * reference_list_of_vertices[i][2];
        }

        InitializeVerticesDueToRotation();

        const unsigned int number_of_faces = reference_list_of_faces.size();

        mListOfFaces.resize(number_of_faces);

        for (int i = 0; i < (int)number_of_faces; i++) {
            const unsigned int number_of_face = reference_list_of_faces[i].size();
            mListOfFaces[i].resize(number_of_face);
            for (int j = 0; j < (int)number_of_face; j++){
                mListOfFaces[i][j] = reference_list_of_faces[i][j];
            }
                
        }

        const double particle_density = this->SlowGetDensity();
        const double polyhedron_volume = reference_volume * scaling_factor * scaling_factor * scaling_factor;
        const double polyhedron_mass = particle_density * polyhedron_volume;

        central_node.FastGetSolutionStepValue(NODAL_MASS) = polyhedron_mass;
        //central_node.FastGetSolutionStepValue(POLYHEDRON_VOLUME) = polyhedron_volume;
        central_node.FastGetSolutionStepValue(REPRESENTATIVE_VOLUME) = polyhedron_volume;

        //SetMass(GetDensity() * CalculateVolume());
        SetMass(polyhedron_mass);
        SetMomentOfInertia();

        KRATOS_CATCH("")
    }

    void PolyhedronParticle::InitializeSolutionStep(const ProcessInfo& r_process_info)
    {
        KRATOS_TRY

        auto& central_node = GetGeometry()[0];
        mRadius = central_node.FastGetSolutionStepValue(RADIUS); //Just in case someone is overwriting the radius in Python

        SetMomentOfInertia();

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

    Vector3 PolyhedronParticle::GetFurthestPoint(Vector3 direction){
        
        KRATOS_TRY

        double max_distance = -1e20;
        auto& central_node = GetGeometry()[0];
        Vector3 max_point;
        for (int i = 0; i < mListOfVertices.size(); ++i) {
            Vector3 vertex_point(mListOfVertices[i][0] + central_node[0], mListOfVertices[i][1] + central_node[1], mListOfVertices[i][2] + central_node[2]);
            double distance = Vector3::Dot(vertex_point, direction);
            if (distance > max_distance){
                max_distance = distance;
                max_point = vertex_point;
            }
        }
        return max_point;

        KRATOS_CATCH("")
    }

    void PolyhedronParticle::SetMomentOfInertia(){

        KRATOS_TRY

        auto& central_node = GetGeometry()[0];
        const unsigned int number_of_vertices = mListOfVertices.size();
        double mass_per_vertex = mRealMass / number_of_vertices;
        Matrix this_moment_of_inertia(3, 3, 0.0);
        double this_identify[3][3] = {
            {1.0, 0.0, 0.0},
            {0.0, 1.0, 0.0},
            {0.0, 0.0, 1.0}
        };

        for (int i = 0; i < (int)number_of_vertices; i++) {
            
            double v_dot_v = GeometryFunctions::DotProduct(mListOfVertices[i], mListOfVertices[i]);
            double v_dot_v_identity[3][3];
            for (int j = 0; j < 3; ++j) {
                for (int k = 0; k < 3; ++k) {
                    v_dot_v_identity[j][k] = v_dot_v * this_identify[j][k];
                }
            }
            double v_outer_v[3][3];
            GeometryFunctions::OuterProduct(mListOfVertices[i], mListOfVertices[i], v_outer_v);

            std::vector<std::vector<double>> result(3, std::vector<double>(3));
            for (int m = 0; m < 3; ++m) {
                for (int n = 0; n < 3; ++n) {
                    this_moment_of_inertia(m,n) += mass_per_vertex * (v_dot_v_identity[m][n] - v_outer_v[m][n]);
                }
            }
        }

        central_node.FastGetSolutionStepValue(POLYHEDRON_MOMENT_OF_INERTIA) = this_moment_of_inertia;

        KRATOS_CATCH("")
    }

    void PolyhedronParticle::Move(const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag) {
        KRATOS_TRY

        GetTranslationalIntegrationScheme().Move(GetGeometry()[0], delta_t, force_reduction_factor, StepFlag);
        if (rotation_option) {
            GetRotationalIntegrationScheme().RotatePolyhedron(this, GetGeometry()[0], delta_t, force_reduction_factor, StepFlag);
        }

        KRATOS_CATCH("")
    }

    void PolyhedronParticle::UpdateVerticesDueToRotation(){
        auto& central_node = GetGeometry()[0];
        array_1d<double, 3>& delta_rotation = central_node.FastGetSolutionStepValue(DELTA_ROTATION);
        double modulus_square = delta_rotation[0]*delta_rotation[0] + delta_rotation[1]*delta_rotation[1] + delta_rotation[2]*delta_rotation[2];
        if (modulus_square != 0.0){
            Quaternion<double> rotation_quaternion = Quaternion<double>::FromRotationVector(delta_rotation);

            for (int k = 0; k < mListOfVertices.size(); k++) {
                rotation_quaternion.RotateVector3(mListOfVertices[k]);
            }
        }
    }

    void PolyhedronParticle::InitializeVerticesDueToRotation(){
        auto& central_node = GetGeometry()[0];
        array_1d<double, 3>& delta_rotation = central_node.FastGetSolutionStepValue(INITIAL_ROTATION_VECTOR);
        double modulus_square = delta_rotation[0]*delta_rotation[0] 
                                + delta_rotation[1]*delta_rotation[1] 
                                + delta_rotation[2]*delta_rotation[2];
        if (modulus_square != 0.0){
            
            Quaternion<double> rotation_quaternion = Quaternion<double>::FromRotationVector(delta_rotation);

            for (int k = 0; k < mListOfVertices.size(); k++) {
                rotation_quaternion.RotateVector3(mListOfVertices[k]);
            }

        }
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
    //double PolyhedronParticle::GetDensity()                                                      { return GetFastProperties()->GetDensity();}
    double PolyhedronParticle::GetDensity()                                                      { return GetProperties()[PARTICLE_DENSITY];}
    double PolyhedronParticle::SlowGetDensity()                                                  { return GetProperties()[PARTICLE_DENSITY];}
    std::vector<array_1d<double, 3>> PolyhedronParticle::GetListOfVertices()                     { return mListOfVertices;}
    std::vector<std::vector<int>> PolyhedronParticle::GetListOfFaces()                           { return mListOfFaces;}

    void   PolyhedronParticle::SetYoungFromProperties(double* young)                             { GetFastProperties()->SetYoungFromProperties( young);                                            }
    void   PolyhedronParticle::SetPoissonFromProperties(double* poisson)                         { GetFastProperties()->SetPoissonFromProperties( poisson);                                        }
    void   PolyhedronParticle::SetDensityFromProperties(double* density)                         { GetFastProperties()->SetDensityFromProperties( density);                                        }
    void   PolyhedronParticle::SetParticleMaterialFromProperties(int* particle_material)         { GetFastProperties()->SetParticleMaterialFromProperties( particle_material);                     }
    
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
