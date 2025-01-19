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
        mIsBelongingToDEMWall = false;
    }

    PolyhedronParticle::PolyhedronParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : SphericParticle(NewId, pGeometry, pProperties) {
        mRadius = 0;
        mRealMass = 0;
        mIsBelongingToDEMWall = false;
    }

    PolyhedronParticle::PolyhedronParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : SphericParticle(NewId, ThisNodes) {
        mRadius = 0;
        mRealMass = 0;
        mIsBelongingToDEMWall = false;
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

    void PolyhedronParticle::InitializeFromFEM(const ProcessInfo& r_process_info, ModelPart& r_fem_model_part){

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
        
        SetRadius();

        // Check if the first Condition type of r_fem_model_part is a tetrahedron
        const auto& first_condition = *r_fem_model_part.ConditionsBegin();
        const unsigned int number_of_vertices = first_condition.GetGeometry().size(); 
        const unsigned int number_of_face_vertices = number_of_vertices; // = first_condition.GetGeometry().size()

        mListOfVertices.resize(number_of_vertices);

        // Ensure the condition has the same id as this polyhedron particle
        for (auto it = r_fem_model_part.GetSubModelPart("SurfaceForPolyWall").ConditionsBegin(); 
             it != r_fem_model_part.GetSubModelPart("SurfaceForPolyWall").ConditionsEnd(); ++it) {
            if (it->Id() == this->Id()) {
                for (int i = 0; i < (int)number_of_vertices; i++) {
                    mListOfVertices[i][0] = it->GetGeometry()[i].Coordinates()[0];
                    mListOfVertices[i][1] = it->GetGeometry()[i].Coordinates()[1];
                    mListOfVertices[i][2] = it->GetGeometry()[i].Coordinates()[2];
                }

                array_1d<double, 3> center = ZeroVector(3);
                for (int i = 0; i < (int)mListOfVertices.size(); i++) {
                    center[0] += mListOfVertices[i][0];
                    center[1] += mListOfVertices[i][1];
                    center[2] += mListOfVertices[i][2];
                }
                center[0] /= (int)mListOfVertices.size();
                center[1] /= (int)mListOfVertices.size();
                center[2] /= (int)mListOfVertices.size();

                KRATOS_ERROR_IF_NOT(mListOfVertices.size()) << "The number of vertices is zero." << std::endl;

                // Calculate the new vertices
                for (int i = 0; i < (int)mListOfVertices.size(); i++) {
                    mListOfVertices[i][0] -= center[0];
                    mListOfVertices[i][1] -= center[1];
                    mListOfVertices[i][2] -= center[2];
                }

                for (int i = 0; i < 3; i++) {
                    central_node.Coordinates()[i] = center[i];
                }

                break;
            }
        }

        //TODO:limited to surface
        const unsigned int number_of_faces = 1;
        mListOfFaces.resize(number_of_faces);
        mListOfFaces[0].resize(number_of_face_vertices);

        if (number_of_vertices == 3) {
            mListOfFaces[0][0] = 0;
            mListOfFaces[0][1] = 1;
            mListOfFaces[0][2] = 2;
        } else if (number_of_vertices == 4) {
            mListOfFaces[0][0] = 0;
            mListOfFaces[0][1] = 1;
            mListOfFaces[0][2] = 2;
            mListOfFaces[0][3] = 3;
        } else {
            KRATOS_ERROR << "The number of vertices of the first condition is not 3 or 4." << std::endl;
            KRATOS_INFO("DEM") << "The number of vertices of the first condition is not 3 or 4." << std::endl;
        }

        //those values are not important for the wall elements
        const double polyhedron_volume = 1.0;
        const double polyhedron_mass = 1.0;

        central_node.FastGetSolutionStepValue(NODAL_MASS) = polyhedron_mass;
        central_node.FastGetSolutionStepValue(REPRESENTATIVE_VOLUME) = polyhedron_volume;

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

    void PolyhedronParticle::UpdateVerticesFromFEM(ModelPart& r_fem_model_part){

        KRATOS_TRY

        auto& central_node = GetGeometry()[0];

        //TODO: "SurfaceForPolyWall" is a temporary name, it should be changed to something more general
        for (auto it = r_fem_model_part.GetSubModelPart("SurfaceForPolyWall").ConditionsBegin(); 
             it != r_fem_model_part.GetSubModelPart("SurfaceForPolyWall").ConditionsEnd(); ++it) {
            if (it->Id() == this->Id()) {
                for (int i = 0; i < (int)mListOfVertices.size(); i++) {
                    mListOfVertices[i][0] = it->GetGeometry()[i].Coordinates()[0];
                    mListOfVertices[i][1] = it->GetGeometry()[i].Coordinates()[1];
                    mListOfVertices[i][2] = it->GetGeometry()[i].Coordinates()[2];
                }
                
                array_1d<double, 3> center = ZeroVector(3);
                for (int i = 0; i < (int)mListOfVertices.size(); i++) {
                    center[0] += mListOfVertices[i][0];
                    center[1] += mListOfVertices[i][1];
                    center[2] += mListOfVertices[i][2];
                }
                center[0] /= (int)mListOfVertices.size();
                center[1] /= (int)mListOfVertices.size();
                center[2] /= (int)mListOfVertices.size();

                KRATOS_ERROR_IF_NOT(mListOfVertices.size()) << "The number of vertices is zero." << std::endl;

                for (int i = 0; i < (int)mListOfVertices.size(); i++) {
                    mListOfVertices[i][0] -= center[0];
                    mListOfVertices[i][1] -= center[1];
                    mListOfVertices[i][2] -= center[2];
                }

                for (int i = 0; i < 3; i++) {
                    central_node[i] = center[i];
                }

                array_1d<double, 3>& velocity = it->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
                array_1d<double, 3>& delta_disp = it->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
                array_1d<double, 3>& angular_velocity = it->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
                array_1d<double, 3>& delta_rotation = it->GetGeometry()[0].FastGetSolutionStepValue(DELTA_ROTATION);
                    
                noalias(central_node.FastGetSolutionStepValue(VELOCITY)) = velocity;
                noalias(central_node.FastGetSolutionStepValue(DELTA_DISPLACEMENT)) = delta_disp;
                noalias(central_node.FastGetSolutionStepValue(ANGULAR_VELOCITY)) = angular_velocity;
                noalias(central_node.FastGetSolutionStepValue(DELTA_ROTATION)) = delta_rotation;

                break;
            }
        }

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
            Matrix rotation_matrix(3, 3, 0.0);
            rotation_quaternion.ToRotationMatrix(rotation_matrix);
            for (int k = 0; k < mListOfVertices.size(); k++) {
                array_1d<double,3> tempVector;
                GeometryFunctions::ProductMatrix3X3Vector3X1(rotation_matrix, mListOfVertices[k], tempVector);
                mListOfVertices[k] = tempVector;
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

    std::vector<Vector3> PolyhedronParticle::GetIntersectingFaceVertices(const Vector3& closestPoint, const Vector3& ContactVector, bool& find_face)
    {
        std::vector<Vector3> bestFaceVertices;
        double bestDotProduct = -1.0; 

        auto& central_node = GetGeometry()[0];

        for (int i = 0; i < mListOfFaces.size(); ++i) {
            const auto& face = mListOfFaces[i];

            if (IsPointOnFace(face, closestPoint)) {

                Vector3 faceNormal = CalculateFaceNormal(face);

                double dotProduct = std::abs(Vector3::Dot(faceNormal.Normalised(), ContactVector.Normalised()));

                if (dotProduct > bestDotProduct) {
                    bestDotProduct = dotProduct;
                    bestFaceVertices.clear();

                    for (int j = 0; j < face.size(); ++j) {
                        Vector3 vertexPoint = {
                            mListOfVertices[face[j]][0] + central_node[0],
                            mListOfVertices[face[j]][1] + central_node[1],
                            mListOfVertices[face[j]][2] + central_node[2]};
                        bestFaceVertices.push_back(vertexPoint);
                    }
                }
            }
        }

        if (bestFaceVertices.empty()) {
            if (mListOfFaces.size() == 1) {
                const auto& defaultFace = mListOfFaces[0];
                for (int j = 0; j < defaultFace.size(); ++j) {
                    Vector3 vertexPoint = {
                        mListOfVertices[defaultFace[j]][0] + central_node[0],
                        mListOfVertices[defaultFace[j]][1] + central_node[1],
                        mListOfVertices[defaultFace[j]][2] + central_node[2]};
                    bestFaceVertices.push_back(vertexPoint);
                }
            } else {
                const auto& defaultFace = mListOfFaces[0];
                for (int j = 0; j < defaultFace.size(); ++j) {
                    Vector3 vertexPoint = {
                        mListOfVertices[defaultFace[j]][0] + central_node[0],
                        mListOfVertices[defaultFace[j]][1] + central_node[1],
                        mListOfVertices[defaultFace[j]][2] + central_node[2]};
                    bestFaceVertices.push_back(vertexPoint);
                }
                find_face = false;
            }
            
        }

        return bestFaceVertices;
    }


    bool PolyhedronParticle::IsPointOnFace(const std::vector<int> face, const Vector3& point)
    {
        // Calculate the normal of the face
        Vector3 normal = CalculateFaceNormal(face);

        auto& central_node = GetGeometry()[0];
        // Check if the point is on the plane of the face
        Vector3 facePoint = {mListOfVertices[face[0]][0] + central_node[0], mListOfVertices[face[0]][1] + central_node[1], mListOfVertices[face[0]][2] + central_node[2]};
        double distance = Vector3::Dot(normal, point - facePoint);
        const double epsilon = 1e-6; // Tolerance for floating-point comparison
        if (std::abs(distance) > epsilon) {
            return false; // Point is not on the plane of the face
        }

        return true;
        /*
        // Project the face vertices and the point onto a 2D plane
        std::vector<Vector3> vertices;
        for (int i = 0; i < face.size(); ++i) {
            Vector3 VertexPoint = {mListOfVertices[face[i]][0] + central_node[0], mListOfVertices[face[i]][1] + central_node[1], mListOfVertices[face[i]][2] + central_node[2]};
            vertices.push_back(VertexPoint);
        } 
        std::vector<Vector2> projectedVertices = ProjectToPlane(vertices, normal);
        Vector2 projectedPoint = ProjectToPlane(point, normal);

        // Check if the projected point is inside the projected polygon
        return IsPointInPolygon(projectedVertices, projectedPoint);*/
        //return true;
    }

    Vector3 PolyhedronParticle::CalculateFaceNormal(const std::vector<int> face)
    {

        // Calculate the normal using the cross product of two edges of the face
        Vector3 normal;
        Vector3 edge1 = {mListOfVertices[face[1]][0] - mListOfVertices[face[0]][0], mListOfVertices[face[1]][1] - mListOfVertices[face[0]][1], mListOfVertices[face[1]][2] - mListOfVertices[face[0]][2]};
        for (size_t i = 0; i < face.size() - 2; ++i) {
            Vector3 edge2 = {mListOfVertices[face[i + 2]][0] - mListOfVertices[face[0]][0], mListOfVertices[face[i + 2]][1] - mListOfVertices[face[0]][1], mListOfVertices[face[i + 2]][2] - mListOfVertices[face[0]][2]};
            normal = Vector3::Cross(edge1, edge2);

            if (normal.Length() > 1e-6) {
                break;
            }
        }
        //Vector3 edge2 = {mListOfVertices[face[2]][0] - mListOfVertices[face[0]][0], mListOfVertices[face[2]][1] - mListOfVertices[face[0]][1], mListOfVertices[face[2]][2] - mListOfVertices[face[0]][2]};
        //Vector3 normal = Vector3::Cross(edge1, edge2);

        // Normalize the normal vector
        normal.Normalise();

        return normal;
    }

    std::vector<Vector2> PolyhedronParticle::ProjectToPlane(const std::vector<Vector3>& vertices, const Vector3& normal)
    {
        // Choose a basis for the plane
        Vector3 u;
        if (std::abs(normal[0]) > std::abs(normal[1])) {
            u = Vector3(-normal[2], 0, normal[0]).Normalised();
        } else {
            u = Vector3(0, -normal[2], normal[1]).Normalised();
        }
        //Vector3 u = normal.Perpendicular();
        Vector3 v = Vector3::Cross(normal, u);

        // Project each vertex onto the plane
        std::vector<Vector2> projectedVertices;
        for (const auto& vertex : vertices) {
            double x = Vector3::Dot(vertex, u);
            double y = Vector3::Dot(vertex, v);
            projectedVertices.push_back(Vector2(x, y));
        }

        return projectedVertices;
    }

    Vector2 PolyhedronParticle::ProjectToPlane(const Vector3& vertex, const Vector3& normal)
    {
        // Choose a basis for the plane
        Vector3 u;
        if (std::abs(normal[0]) > std::abs(normal[1])) {
            u = Vector3(-normal[2], 0, normal[0]).Normalised();
        } else {
            u = Vector3(0, -normal[2], normal[1]).Normalised();
        }
        //Vector3 u = normal.Perpendicular();
        Vector3 v = Vector3::Cross(normal, u);

        // Project the vertex onto the plane
        double x = Vector3::Dot(vertex, u);
        double y = Vector3::Dot(vertex, v);
        return Vector2(x, y);
    }

    bool PolyhedronParticle::IsPointInPolygon(const std::vector<Vector2>& polygon, const Vector2& point)
    {
        // Implement a point-in-polygon test (e.g., ray-casting algorithm)
        int intersections = 0;
        for (size_t i = 0; i < polygon.size(); ++i) {
            Vector2 v1 = polygon[i];
            Vector2 v2 = polygon[(i + 1) % polygon.size()];

            if (((v1.y > point.y) != (v2.y > point.y)) &&
                (point.x < (v2.x - v1.x) * (point.y - v1.y) / (v2.y - v1.y) + v1.x)) {
                intersections++;
            }
        }

        return (intersections % 2) != 0;
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
