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
    : RigidBodyElement3D(NewId, pGeometry) {}

    PolyhedronParticle::PolyhedronParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : RigidBodyElement3D(NewId, pGeometry, pProperties) {}

    PolyhedronParticle::PolyhedronParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : RigidBodyElement3D(NewId, ThisNodes) {}

    Element::Pointer PolyhedronParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
        return Element::Pointer(new PolyhedronParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    PolyhedronParticle::~PolyhedronParticle() {}

    void PolyhedronParticle::InitializeSolutionStep(const ProcessInfo& r_process_info)
    {
        KRATOS_TRY

        auto& central_node = GetGeometry()[0];
        mRadius = central_node.FastGetSolutionStepValue(RADIUS); //Just in case someone is overwriting the radius in Python
        central_node.FastGetSolutionStepValue(REPRESENTATIVE_VOLUME) = CalculateVolume();

        if (this->Is(DEMFlags::HAS_ROTATION)) {
            if (this->Is(DEMFlags::HAS_ROLLING_FRICTION)) {
                if (mRollingFrictionModel != nullptr){
                    mRollingFrictionModel->InitializeSolutionStep();
                }
            }
        }

        KRATOS_CATCH("")
    }

    void PolyhedronParticle::CustomInitialize(ModelPart& rigid_body_element_sub_model_part) {
        
        KRATOS_TRY
        
        RigidBodyElement3D::CustomInitialize(rigid_body_element_sub_model_part);
        
        mEnginePower = rigid_body_element_sub_model_part[DEM_ENGINE_POWER];
        mMaxEngineForce = rigid_body_element_sub_model_part[DEM_MAX_ENGINE_FORCE];
        mThresholdVelocity = rigid_body_element_sub_model_part[DEM_THRESHOLD_VELOCITY];
        mEnginePerformance = rigid_body_element_sub_model_part[DEM_ENGINE_PERFORMANCE];
        
        mDragConstantVector = ZeroVector(3);
        mDragConstantVector[0] = rigid_body_element_sub_model_part[DEM_DRAG_CONSTANT_X];
        mDragConstantVector[1] = rigid_body_element_sub_model_part[DEM_DRAG_CONSTANT_Y];
        mDragConstantVector[2] = rigid_body_element_sub_model_part[DEM_DRAG_CONSTANT_Z];

        KRATOS_CATCH("")
    }
    
    /*
    void PolyhedronParticle::ComputeWaterDragForce() {

        KRATOS_TRY

        const double water_density = 1000;
        const double drag_coefficient = 0.75; //https://www.engineeringtoolbox.com/drag-coefficient-d_627.html
        const double water_level = 0.0;

        for (unsigned int i = 0; i != mListOfRigidFaces.size(); ++i) {

            double rigid_face_area = 0.0;
            Point rigid_face_centroid;
            array_1d<double, 3> drag_force = ZeroVector(3);
            array_1d<double, 3> rigid_body_centroid_to_rigid_face_centroid_vector = ZeroVector(3);
            array_1d<double, 3> drag_moment = ZeroVector(3);
            size_t RF_size = mListOfRigidFaces[i]->GetGeometry().size();
            array_1d<double, 3> rigid_face_z_coords_values = ZeroVector(RF_size);
            unsigned int number_of_RF_nodes_out_of_water = 0;

            for (unsigned int j = 0; j < RF_size; j++) {
                rigid_face_z_coords_values[j] = mListOfRigidFaces[i]->GetGeometry()[j].Coordinates()[2];
                if (rigid_face_z_coords_values[j] > water_level) number_of_RF_nodes_out_of_water++;
            }
            
            if (number_of_RF_nodes_out_of_water == RF_size) continue;
                        
            array_1d<double, 3> rigid_face_velocity;
            noalias(rigid_face_velocity) = mListOfRigidFaces[i]->GetVelocity();
            
            double velocity_modulus = DEM_MODULUS_3(rigid_face_velocity);
            
            if (velocity_modulus) {
                DEM_MULTIPLY_BY_SCALAR_3(rigid_face_velocity, 1.0/velocity_modulus)
            }
            
            rigid_face_centroid = mListOfRigidFaces[i]->GetGeometry().Center();
            rigid_face_area = mListOfRigidFaces[i]->GetGeometry().Area();
                        
            DEM_MULTIPLY_BY_SCALAR_3(rigid_face_velocity, -0.5 * drag_coefficient * water_density * velocity_modulus * velocity_modulus * rigid_face_area)
            DEM_COPY_SECOND_TO_FIRST_3(drag_force, rigid_face_velocity)
                                
            rigid_body_centroid_to_rigid_face_centroid_vector[0] = rigid_face_centroid.Coordinates()[0] - GetGeometry()[0].Coordinates()[0];
            rigid_body_centroid_to_rigid_face_centroid_vector[1] = rigid_face_centroid.Coordinates()[1] - GetGeometry()[0].Coordinates()[1];
            rigid_body_centroid_to_rigid_face_centroid_vector[2] = rigid_face_centroid.Coordinates()[2] - GetGeometry()[0].Coordinates()[2];

            array_1d<double, 3>& total_forces = GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
            DEM_ADD_SECOND_TO_FIRST(total_forces, drag_force)
            array_1d<double, 3>& total_moments = GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);
            GeometryFunctions::CrossProduct(rigid_body_centroid_to_rigid_face_centroid_vector, drag_force, drag_moment);
            DEM_ADD_SECOND_TO_FIRST(total_moments, drag_moment)
        }

        KRATOS_CATCH("")
    } */

    
    void PolyhedronParticle::ComputeExternalForces(const array_1d<double,3>& gravity) {

        KRATOS_TRY

        noalias(GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES)) += RigidBodyElement3D::GetMass() * gravity;
        
        ComputeBuoyancyEffects();
        ComputeEngineForce();
        ComputeWaterDragForce();
        
        const array_1d<double, 3> external_applied_torque = GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT);
        noalias(GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT)) += external_applied_torque;

        KRATOS_CATCH("")
    }

    void PolyhedronParticle::ComputeNewNeighboursHistoricalData(DenseVector<int>& temp_neighbours_ids,
                                                         std::vector<array_1d<double, 3> >& temp_neighbour_elastic_contact_forces)
    {
        KRATOS_TRY

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

        KRATOS_CATCH("")
    }

}  // namespace Kratos
