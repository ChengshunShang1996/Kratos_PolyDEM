/////////////////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
// Date: July 2024
/////////////////////////////////////////////////////////////

//Important: Current contact-based explicit strategy solver is not appliable to Clusters!!!

#include "explicit_solver_contact_strategy.h"
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"

namespace Kratos {

    void ContactExplicitSolverStrategy::RebuildPropertiesProxyPointersForPolyhedron(std::vector<PolyhedronParticle*>& rCustomListOfPolyhedronParticles) {
        //This function is called for the local mesh and the ghost mesh, so mListOfSphericElements must not be used here.
        KRATOS_TRY


        std::vector<PropertiesProxy>& vector_of_properties_proxies = PropertiesProxiesManager().GetPropertiesProxies(*mpPolyhedron_model_part);

        IndexPartition<unsigned int>(rCustomListOfPolyhedronParticles.size()).for_each([&](unsigned int i){
            rCustomListOfPolyhedronParticles[i]->SetFastProperties(vector_of_properties_proxies);
        });

        return;
        KRATOS_CATCH("")
    }
    
    void ContactExplicitSolverStrategy::Initialize() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        ModelPart& fem_model_part = GetFemModelPart();
        ModelPart& polyhedron_model_part = GetPolyhedronModelPart();

        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        SendProcessInfoToClustersModelPart();

        if (polyhedron_model_part.GetCommunicator().MyPID() == 0) {
            KRATOS_INFO("DEM") << "------------------CONTACT-BASED EXPLICIT SOLVER STRATEGY---------------------" << "\n" << std::endl;
        }

        mNumberOfThreads = ParallelUtilities::GetNumThreads();
        DisplayThreadInfo();

        RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
        RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);
        RebuildListOfSphericParticles <PolyhedronParticle> (polyhedron_model_part.GetCommunicator().LocalMesh().Elements(), mListOfPolyhedronParticles);
        RebuildListOfSphericParticles <PolyhedronParticle> (polyhedron_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostPolyhedronParticles);

        mSearchControlVector.resize(mNumberOfThreads);
        for (int i = 0; i < mNumberOfThreads; i++) mSearchControlVector[i] = 0;

        PropertiesProxiesManager().CreatePropertiesProxies(*mpDem_model_part, *mpInlet_model_part, *mpCluster_model_part, *mpPolyhedron_model_part);

        bool has_mpi = false;
        Check_MPI(has_mpi);

        if (has_mpi) {
            RepairPointersToNormalProperties(mListOfSphericParticles);
            RepairPointersToNormalProperties(mListOfGhostSphericParticles);
            RepairPointersToNormalPropertiesOfPolyhedron(mListOfPolyhedronParticles);
            RepairPointersToNormalPropertiesOfPolyhedron(mListOfGhostPolyhedronParticles);
        }

        RebuildPropertiesProxyPointers(mListOfSphericParticles);
        RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);
        RebuildPropertiesProxyPointersForPolyhedron(mListOfPolyhedronParticles);
        RebuildPropertiesProxyPointersForPolyhedron(mListOfGhostPolyhedronParticles);

        GetSearchControl() = r_process_info[SEARCH_CONTROL];

        InitializeDEMElements();
        InitializeFEMElements();
        InitializeClusters(); // This adds elements to the balls modelpart
        InitializePolyhedrons();

        UpdateMaxIdOfCreatorDestructor();

        RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
        RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostSphericParticles);
        RebuildListOfSphericParticles <PolyhedronParticle> (polyhedron_model_part.GetCommunicator().LocalMesh().Elements(), mListOfPolyhedronParticles);
        RebuildListOfSphericParticles <PolyhedronParticle> (polyhedron_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostPolyhedronParticles);

        InitializeSolutionStep();
        ApplyInitialConditions();

        if (r_model_part.Nodes().size() > 0) {
            SetSearchRadiiOnAllParticles(*mpDem_model_part, mpDem_model_part->GetProcessInfo()[SEARCH_RADIUS_INCREMENT], 1.0);
            SearchNeighbours();
            ComputeNewNeighboursHistoricalData();
        }

        if (fem_model_part.Nodes().size() > 0) {
            SetSearchRadiiWithFemOnAllParticles(*mpDem_model_part, mpDem_model_part->GetProcessInfo()[SEARCH_RADIUS_INCREMENT_FOR_WALLS], 1.0);
            SearchRigidFaceNeighbours();
            ComputeNewRigidFaceNeighboursHistoricalData();
        }

        SetSearchRadiiOnAllPolyhedronParticles(*mpPolyhedron_model_part, mpDem_model_part->GetProcessInfo()[SEARCH_RADIUS_INCREMENT_FOR_WALLS], 1.0);
        SearchPolyhedronNeighbours();
        ComputePolyhedronNewNeighboursHistoricalData();

        //if (r_process_info[CONTACT_MESH_OPTION]) {
            // TODO: 
            // I comment this function temperally as there is no sphere particles actually. 
            // It is more complex to take into account of different types of particles
            // in the future, it could be improved

            // CreateContactElements(); //only for spheres
            // CreatePolyhedronContactElements();
            // InitializePolyhedronContactElements();
        //}

        CreatePolyhedronContactElements();
        //InitializePolyhedronContactElements();

        //r_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(NEIGHBOUR_IDS);
        //r_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(NEIGHBOURS_CONTACT_AREAS);

        ComputeNodalArea();

        KRATOS_CATCH("")
    }// Initialize()

    void ContactExplicitSolverStrategy::RepairPointersToNormalPropertiesOfPolyhedron(std::vector<PolyhedronParticle*>& rCustomListOfPolyhedronParticles) {

        /*
        KRATOS_TRY

        bool found = false;
        // Using IndexPartition should be fine since 'break' affects the internal for loops while the replaced continues only has an effect on the for_each loop.
        IndexPartition<unsigned int>(rCustomListOfPolyhedronParticles.size()).for_each([&](unsigned int i){

            int own_properties_id = rCustomListOfPolyhedronParticles[i]->GetProperties().Id();
            for (PropertiesIterator props_it = mpPolyhedron_model_part->GetMesh(0).PropertiesBegin(); props_it != mpPolyhedron_model_part->GetMesh(0).PropertiesEnd(); props_it++) {
                int model_part_id = props_it->GetId();
                if (own_properties_id == model_part_id) {
                    rCustomListOfPolyhedronParticles[i]->SetProperties(*(props_it.base()));
                    found = true;
                    break;
                }
            }

            KRATOS_ERROR_IF_NOT(found) << "This particle could not find its properties!!" << std::endl;
        });

        KRATOS_CATCH("")
        */
        KRATOS_TRY

        bool found = false;
        // Using IndexPartition should be fine since 'break' affects the internal for loops while the replaced continues only has an effect on the for_each loop.
        IndexPartition<unsigned int>(rCustomListOfPolyhedronParticles.size()).for_each([&](unsigned int i){

            int own_properties_id = rCustomListOfPolyhedronParticles[i]->GetProperties().Id();
            for (PropertiesIterator props_it = mpDem_model_part->GetMesh(0).PropertiesBegin(); props_it != mpDem_model_part->GetMesh(0).PropertiesEnd(); props_it++) {
                int model_part_id = props_it->GetId();
                if (own_properties_id == model_part_id) {
                    rCustomListOfPolyhedronParticles[i]->SetProperties(*(props_it.base()));
                    found = true;
                    break;
                }
            }
            if (found) return;

            for (PropertiesIterator props_it = mpInlet_model_part->GetMesh(0).PropertiesBegin(); props_it != mpInlet_model_part->GetMesh(0).PropertiesEnd(); props_it++) {
                int model_part_id = props_it->GetId();
                if (own_properties_id == model_part_id) {
                    rCustomListOfPolyhedronParticles[i]->SetProperties(*(props_it.base()));
                    found = true;
                    break;
                }
            }
            if (found) return;

            for (PropertiesIterator props_it = mpCluster_model_part->GetMesh(0).PropertiesBegin(); props_it != mpCluster_model_part->GetMesh(0).PropertiesEnd(); props_it++) {
                int model_part_id = props_it->GetId();
                if (own_properties_id == model_part_id) {
                    rCustomListOfPolyhedronParticles[i]->SetProperties(*(props_it.base()));
                    found = true;
                    break;
                }
            }
            if (found) return;

            for (PropertiesIterator props_it = mpPolyhedron_model_part->GetMesh(0).PropertiesBegin(); props_it != mpPolyhedron_model_part->GetMesh(0).PropertiesEnd(); props_it++) {
                int model_part_id = props_it->GetId();
                if (own_properties_id == model_part_id) {
                    rCustomListOfPolyhedronParticles[i]->SetProperties(*(props_it.base()));
                    found = true;
                    break;
                }
            }

            KRATOS_ERROR_IF_NOT(found) << "This particle could not find its properties!!" << std::endl;
        });

        KRATOS_CATCH("")
    }

    void ContactExplicitSolverStrategy::InitializePolyhedrons() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

        IndexPartition<unsigned int>(mListOfPolyhedronParticles.size()).for_each([&](unsigned int i){
            mListOfPolyhedronParticles[i]->Initialize(r_process_info);
        });

        KRATOS_CATCH("")
    }

    void ContactExplicitSolverStrategy::UpdateMaxIdOfCreatorDestructor() {

        KRATOS_TRY

        int max_Id = mpParticleCreatorDestructor->GetCurrentMaxNodeId();
        ModelPart& r_model_part = GetModelPart();
        int max_DEM_Id = mpParticleCreatorDestructor->FindMaxNodeIdInModelPart(r_model_part);
        int max_FEM_Id = mpParticleCreatorDestructor->FindMaxNodeIdInModelPart(*mpFem_model_part);
        int max_cluster_Id = mpParticleCreatorDestructor->FindMaxNodeIdInModelPart(*mpCluster_model_part);
        int max_polyhedron_Id = mpParticleCreatorDestructor->FindMaxNodeIdInModelPart(*mpPolyhedron_model_part);

        max_Id = std::max(max_Id, max_DEM_Id);
        max_Id = std::max(max_Id, max_FEM_Id);
        max_Id = std::max(max_Id, max_cluster_Id);
        max_Id = std::max(max_Id, max_polyhedron_Id);
        mpParticleCreatorDestructor->SetMaxNodeId(max_Id);

        KRATOS_CATCH("")
    }

    void ContactExplicitSolverStrategy::InitializeSolutionStep() {
        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        const ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        ElementsArrayType& pElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        ModelPart& r_fem_model_part = GetFemModelPart();
        const ProcessInfo& r_fem_process_info = r_fem_model_part.GetProcessInfo();
        ConditionsArrayType& pConditions = r_fem_model_part.GetCommunicator().LocalMesh().Conditions();

        ModelPart& r_polyhedron_model_part = GetPolyhedronModelPart();
        const ProcessInfo& r_polyhedron_process_info = r_polyhedron_model_part.GetProcessInfo();
        ElementsArrayType& pPolyElements = r_polyhedron_model_part.GetCommunicator().LocalMesh().Elements();

        RebuildListOfSphericParticles <SphericParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfSphericParticles);
        RebuildListOfSphericParticles <PolyhedronParticle> (r_polyhedron_model_part.GetCommunicator().LocalMesh().Elements(), mListOfPolyhedronParticles);

        SetNormalRadiiOnAllParticles(*mpDem_model_part);
        //TODO: do I need it?

        #pragma omp parallel
        {
            #pragma omp for nowait
            for (int k = 0; k < (int) pElements.size(); k++) {
                ElementsArrayType::iterator it = pElements.ptr_begin() + k;
                (it)->InitializeSolutionStep(r_process_info);
            }

            #pragma omp for nowait
            for (int k = 0; k < (int) pConditions.size(); k++) {
                ConditionsArrayType::iterator it = pConditions.ptr_begin() + k;
                (it)->InitializeSolutionStep(r_fem_process_info);
            }

            #pragma omp for nowait
            for (int k = 0; k < (int) pPolyElements.size(); k++) {
                ElementsArrayType::iterator it = pPolyElements.ptr_begin() + k;
                (it)->InitializeSolutionStep(r_polyhedron_process_info);  //TODO: what should be in this r_polyhedron_process_info
            }
        }

        ApplyPrescribedBoundaryConditions();
        KRATOS_CATCH("")
    }

    void ContactExplicitSolverStrategy::ApplyPrescribedBoundaryConditions(){
        
        KRATOS_TRY

        BaseType::ApplyPrescribedBoundaryConditions();

        ApplyPrescribedBoundaryConditionsForPolyhedron();

        KRATOS_CATCH("")
    }

    void ContactExplicitSolverStrategy::ApplyPrescribedBoundaryConditionsForPolyhedron() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        const ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();
        const double time = r_process_info[TIME];

        ModelPart& polyhedron_model_part = GetPolyhedronModelPart();

        unsigned int polyhedron_elements_counter = 0;

        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = polyhedron_model_part.SubModelPartsBegin(); sub_model_part != polyhedron_model_part.SubModelPartsEnd(); ++sub_model_part) {

            double vel_start = 0.0, vel_stop = std::numeric_limits<double>::max();
            if ((*sub_model_part).Has(VELOCITY_START_TIME)) {
                vel_start = (*sub_model_part)[VELOCITY_START_TIME];
            }
            if ((*sub_model_part).Has(VELOCITY_STOP_TIME)) {
                vel_stop = (*sub_model_part)[VELOCITY_STOP_TIME];
            }

            if (time < vel_start || time > vel_stop) continue;

            NodesArrayType& pNodes = sub_model_part->Nodes();

            if ((*sub_model_part).Has(IMPOSED_VELOCITY_X_VALUE)) {
                SetFlagAndVariableToNodes(DEMFlags::FIXED_VEL_X, VELOCITY_X, (*sub_model_part)[IMPOSED_VELOCITY_X_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(IMPOSED_VELOCITY_Y_VALUE)) {
                SetFlagAndVariableToNodes(DEMFlags::FIXED_VEL_Y, VELOCITY_Y, (*sub_model_part)[IMPOSED_VELOCITY_Y_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(IMPOSED_VELOCITY_Z_VALUE)) {
                SetFlagAndVariableToNodes(DEMFlags::FIXED_VEL_Z, VELOCITY_Z, (*sub_model_part)[IMPOSED_VELOCITY_Z_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(IMPOSED_ANGULAR_VELOCITY_X_VALUE)) {
                SetFlagAndVariableToNodes(DEMFlags::FIXED_ANG_VEL_X, ANGULAR_VELOCITY_X, (*sub_model_part)[IMPOSED_ANGULAR_VELOCITY_X_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(IMPOSED_ANGULAR_VELOCITY_Y_VALUE)) {
                SetFlagAndVariableToNodes(DEMFlags::FIXED_ANG_VEL_Y, ANGULAR_VELOCITY_Y, (*sub_model_part)[IMPOSED_ANGULAR_VELOCITY_Y_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(IMPOSED_ANGULAR_VELOCITY_Z_VALUE)) {
                SetFlagAndVariableToNodes(DEMFlags::FIXED_ANG_VEL_Z, ANGULAR_VELOCITY_Z, (*sub_model_part)[IMPOSED_ANGULAR_VELOCITY_Z_VALUE], pNodes);
            }
        } // for each mesh

        KRATOS_CATCH("")
    }

    void ContactExplicitSolverStrategy::ApplyInitialConditionsPolyhedron() {

        KRATOS_TRY
        ModelPart& polyhedron_model_part = GetPolyhedronModelPart();

        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = polyhedron_model_part.SubModelPartsBegin(); sub_model_part != polyhedron_model_part.SubModelPartsEnd(); ++sub_model_part) {

            NodesArrayType& pNodes = sub_model_part->Nodes();

            if ((*sub_model_part).Has(INITIAL_VELOCITY_X_VALUE)) {
                SetVariableToNodes(VELOCITY_X, (*sub_model_part)[INITIAL_VELOCITY_X_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(INITIAL_VELOCITY_Y_VALUE)) {
                SetVariableToNodes(VELOCITY_Y, (*sub_model_part)[INITIAL_VELOCITY_Y_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(INITIAL_VELOCITY_Z_VALUE)) {
                SetVariableToNodes(VELOCITY_Z, (*sub_model_part)[INITIAL_VELOCITY_Z_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(INITIAL_ANGULAR_VELOCITY_X_VALUE)) {
                SetVariableToNodes(ANGULAR_VELOCITY_X, (*sub_model_part)[INITIAL_ANGULAR_VELOCITY_X_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(INITIAL_ANGULAR_VELOCITY_Y_VALUE)) {
                SetVariableToNodes(ANGULAR_VELOCITY_Y, (*sub_model_part)[INITIAL_ANGULAR_VELOCITY_Y_VALUE], pNodes);
            }
            if ((*sub_model_part).Has(INITIAL_ANGULAR_VELOCITY_Z_VALUE)) {
                SetVariableToNodes(ANGULAR_VELOCITY_Z, (*sub_model_part)[INITIAL_ANGULAR_VELOCITY_Z_VALUE], pNodes);
            }
        } // for each mesh

        KRATOS_CATCH("")
    }

    void ContactExplicitSolverStrategy::ApplyInitialConditions(){
                
        KRATOS_TRY

        BaseType::ApplyInitialConditions();

        ApplyInitialConditionsPolyhedron();

        KRATOS_CATCH("")
    }

    void ContactExplicitSolverStrategy::ResetPrescribedMotionFlagsRespectingImposedDofsForPolyhedron() {
        KRATOS_TRY
        ModelPart& polyhedron_model_part = GetPolyhedronModelPart();

        NodesArrayType& polyhedron_model_part_nodes = polyhedron_model_part.Nodes();

        if (!polyhedron_model_part_nodes.size()) return;

        const unsigned int vel_x_dof_position = (polyhedron_model_part.NodesBegin())->GetDofPosition(VELOCITY_X);
        const unsigned int ang_vel_x_dof_position = (polyhedron_model_part.NodesBegin())->GetDofPosition(ANGULAR_VELOCITY_X);


        block_for_each(polyhedron_model_part_nodes, [&](ModelPart::NodeType& rNode) {

            if (rNode.Is(BLOCKED)) return;
            Node& node = rNode;

            if (node.GetDof(VELOCITY_X, vel_x_dof_position).IsFixed()) {
                node.Set(DEMFlags::FIXED_VEL_X, true);
            } else {
                node.Set(DEMFlags::FIXED_VEL_X, false);
            }
            if (node.GetDof(VELOCITY_Y, vel_x_dof_position + 1).IsFixed()) {
                node.Set(DEMFlags::FIXED_VEL_Y, true);
            } else {
                node.Set(DEMFlags::FIXED_VEL_Y, false);
            }
            if (node.GetDof(VELOCITY_Z, vel_x_dof_position + 2).IsFixed()) {
                node.Set(DEMFlags::FIXED_VEL_Z, true);
            } else {
                node.Set(DEMFlags::FIXED_VEL_Z, false);
            }
            if (node.GetDof(ANGULAR_VELOCITY_X, ang_vel_x_dof_position).IsFixed()) {
                node.Set(DEMFlags::FIXED_ANG_VEL_X, true);
            } else {
                node.Set(DEMFlags::FIXED_ANG_VEL_X, false);
            }
            if (node.GetDof(ANGULAR_VELOCITY_Y, ang_vel_x_dof_position + 1).IsFixed()) {
                node.Set(DEMFlags::FIXED_ANG_VEL_Y, true);
            } else {
                node.Set(DEMFlags::FIXED_ANG_VEL_Y, false);
            }
            if (node.GetDof(ANGULAR_VELOCITY_Z, ang_vel_x_dof_position + 2).IsFixed()) {
                node.Set(DEMFlags::FIXED_ANG_VEL_Z, true);
            } else {
                node.Set(DEMFlags::FIXED_ANG_VEL_Z, false);
            }
        });
        KRATOS_CATCH("")
    }
    
    double ContactExplicitSolverStrategy::SolveSolutionStep() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        ModelPart& r_polyhedron_model_part = GetPolyhedronModelPart();

        bool has_mpi = false;
        VariablesList r_modelpart_nodal_variables_list = r_model_part.GetNodalSolutionStepVariablesList();
        if (r_modelpart_nodal_variables_list.Has(PARTITION_INDEX)) has_mpi = true;

        SearchDEMOperations(r_model_part, has_mpi);
        SearchFEMOperations(r_model_part, has_mpi);
        SearchPolyhedronOperations(r_model_part, r_polyhedron_model_part, has_mpi);

        ForceOperations(r_model_part, r_polyhedron_model_part);
        PerformTimeIntegrationOfMotion();

        KRATOS_CATCH("")

        return 0.0;

    }//SolveSolutionStep()

    void ContactExplicitSolverStrategy::PerformTimeIntegrationOfMotion(int StepFlag) {

        KRATOS_TRY

        ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();
        double delta_t = r_process_info[DELTA_TIME];
        double virtual_mass_coeff = r_process_info[NODAL_MASS_COEFF]; //TODO: change the name of this variable to FORCE_REDUCTION_FACTOR
        bool virtual_mass_option = (bool) r_process_info[VIRTUAL_MASS_OPTION];
        double force_reduction_factor = 1.0;
        if (virtual_mass_option) {
            force_reduction_factor = virtual_mass_coeff;
            KRATOS_ERROR_IF((force_reduction_factor > 1.0) || (force_reduction_factor < 0.0)) << "The force reduction factor is either larger than 1 or negative: FORCE_REDUCTION_FACTOR= "<< virtual_mass_coeff << std::endl;
        }

        bool rotation_option = r_process_info[ROTATION_OPTION];

        const int number_of_particles       = (int) mListOfSphericParticles.size();
        const int number_of_ghost_particles = (int) mListOfGhostSphericParticles.size();
        const int number_of_polyhedron_particles       = (int) mListOfPolyhedronParticles.size();
        const int number_of_ghost_polyhedron_particles = (int) mListOfGhostPolyhedronParticles.size();

        ModelPart& r_clusters_model_part  = *mpCluster_model_part;
        ElementsArrayType& pLocalClusters = r_clusters_model_part.GetCommunicator().LocalMesh().Elements();
        ElementsArrayType& pGhostClusters = r_clusters_model_part.GetCommunicator().GhostMesh().Elements();
        ModelPart& r_fem_model_part  = *mpFem_model_part;
        ElementsArrayType& pFemElements = r_fem_model_part.GetCommunicator().LocalMesh().Elements();

        #pragma omp parallel
        {
            #pragma omp for nowait
            for (int i = 0; i < number_of_particles; i++) {
                mListOfSphericParticles[i]->Move(delta_t, rotation_option, force_reduction_factor, StepFlag);
            }

            #pragma omp for nowait
            for (int i = 0; i < number_of_ghost_particles; i++) {
                mListOfGhostSphericParticles[i]->Move(delta_t, rotation_option, force_reduction_factor, StepFlag);
            }

            #pragma omp for nowait
            for (int k = 0; k < (int) pLocalClusters.size(); k++) {
                ElementsArrayType::iterator it = pLocalClusters.ptr_begin() + k;
                Cluster3D& cluster_element = dynamic_cast<Kratos::Cluster3D&> (*it);
                cluster_element.RigidBodyElement3D::Move(delta_t, rotation_option, force_reduction_factor, StepFlag);
            }

            #pragma omp for nowait
            for (int k = 0; k < (int) pGhostClusters.size(); k++) {
                ElementsArrayType::iterator it = pGhostClusters.ptr_begin() + k;
                Cluster3D& cluster_element = dynamic_cast<Kratos::Cluster3D&> (*it);
                cluster_element.RigidBodyElement3D::Move(delta_t, rotation_option, force_reduction_factor, StepFlag);
            }

            #pragma omp for nowait
            for (int k = 0; k < (int) pFemElements.size(); k++) {
                ElementsArrayType::iterator it = pFemElements.ptr_begin() + k;
                RigidBodyElement3D& rigid_body_element = dynamic_cast<Kratos::RigidBodyElement3D&> (*it);
                rigid_body_element.Move(delta_t, rotation_option, force_reduction_factor, StepFlag);
            }

            #pragma omp for nowait
            for (int i = 0; i < number_of_polyhedron_particles; i++) {
                mListOfPolyhedronParticles[i]->Move(delta_t, rotation_option, force_reduction_factor, StepFlag);
            }

            #pragma omp for nowait
            for (int i = 0; i < number_of_ghost_polyhedron_particles; i++) {
                mListOfGhostPolyhedronParticles[i]->Move(delta_t, rotation_option, force_reduction_factor, StepFlag);
            }
        }

        //GetScheme()->Calculate(GetModelPart(), StepFlag);
        //GetScheme()->Calculate(*mpCluster_model_part, StepFlag);
        KRATOS_CATCH("")
    }

    void ContactExplicitSolverStrategy::SearchPolyhedronOperations(ModelPart& r_model_part, ModelPart& polyhedron_model_part, bool has_mpi) {

        KRATOS_TRY

        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        int time_step = r_process_info[TIME_STEPS];
        const double time = r_process_info[TIME];
        const bool is_time_to_search_neighbours = (time_step + 1) % mNStepSearch == 0 && (time_step > 0); //Neighboring search. Every N times.
        const bool is_time_to_print_results = r_process_info[IS_TIME_TO_PRINT];
        const bool is_time_to_mark_and_remove = is_time_to_search_neighbours && (r_process_info[BOUNDING_BOX_OPTION] && time >= r_process_info[BOUNDING_BOX_START_TIME] && time <= r_process_info[BOUNDING_BOX_STOP_TIME]);
        //BoundingBoxUtility(is_time_to_mark_and_remove);

        if (is_time_to_search_neighbours) {
            //if (!is_time_to_mark_and_remove) { //Just in case that some entities were marked as TO_ERASE without a bounding box (manual removal)
            //    mpParticleCreatorDestructor->DestroyParticles<Cluster3D>(*mpCluster_model_part);
            //    mpParticleCreatorDestructor->DestroyParticles<SphericParticle>(r_model_part);
            //}

            RebuildListOfSphericParticles <PolyhedronParticle> (polyhedron_model_part.GetCommunicator().LocalMesh().Elements(), mListOfPolyhedronParticles);
            RebuildListOfSphericParticles <PolyhedronParticle> (polyhedron_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostPolyhedronParticles);

            SearchPolyhedronNeighbours();

            RebuildListOfSphericParticles <PolyhedronParticle> (polyhedron_model_part.GetCommunicator().LocalMesh().Elements(), mListOfPolyhedronParticles);
            RebuildListOfSphericParticles <PolyhedronParticle> (polyhedron_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostPolyhedronParticles);
            
            RepairPointersToNormalPropertiesOfPolyhedron(mListOfPolyhedronParticles);
            RepairPointersToNormalPropertiesOfPolyhedron(mListOfGhostPolyhedronParticles);
            RebuildPropertiesProxyPointersForPolyhedron(mListOfPolyhedronParticles);
            RebuildPropertiesProxyPointersForPolyhedron(mListOfGhostPolyhedronParticles);

            ComputePolyhedronNewNeighboursHistoricalData();

            mSearchControl = 2; // Search is active and has been performed during this time step
        } else {
            mSearchControl = 1; // Search is active but no search has been done this time step;
        }

        //if (is_time_to_print_results && r_process_info[CONTACT_MESH_OPTION] == 1) {
        CreatePolyhedronContactElements();
        //InitializePolyhedronContactElements();
        //}

        KRATOS_CATCH("")
    }

    void ContactExplicitSolverStrategy::SearchPolyhedronNeighbours() {
        KRATOS_TRY

        if (!mDoSearchNeighbourElements) {
            return;
        }

        ModelPart& polyhedron_model_part = GetPolyhedronModelPart();

        int number_of_elements = polyhedron_model_part.GetCommunicator().LocalMesh().ElementsArray().end() - polyhedron_model_part.GetCommunicator().LocalMesh().ElementsArray().begin();
        if (!number_of_elements) return;

        GetResultsPoly().resize(number_of_elements);
        GetResultsDistancesPoly().resize(number_of_elements);

        mpSpSearch->SearchElementsInRadiusExclusive(polyhedron_model_part, this->GetArrayOfAmplifiedRadiiPoly(), this->GetResultsPoly(), this->GetResultsDistancesPoly());
        const int number_of_particles = (int) mListOfPolyhedronParticles.size();

        typedef std::map<PolyhedronParticle*,std::vector<PolyhedronParticle*>> ConnectivitiesMap;
        std::vector<ConnectivitiesMap> thread_maps_of_connectivities;
        thread_maps_of_connectivities.resize(ParallelUtilities::GetNumThreads());

        #pragma omp parallel for schedule(dynamic, 100)
        for (int i = 0; i < number_of_particles; i++) {
            mListOfPolyhedronParticles[i]->mNeighbourElements.clear();
            for (SpatialSearch::ResultElementsContainerType::iterator neighbour_it = this->GetResultsPoly()[i].begin(); neighbour_it != this->GetResultsPoly()[i].end(); ++neighbour_it) {
                Element* p_neighbour_element = (*neighbour_it).get();
                PolyhedronParticle* p_polyhedron_neighbour_particle = dynamic_cast<PolyhedronParticle*> (p_neighbour_element);
                //if (mListOfPolyhedronParticles[i]->Is(DEMFlags::BELONGS_TO_A_CLUSTER) && (mListOfPolyhedronParticles[i]->GetClusterId() == p_polyhedron_neighbour_particle->GetClusterId())) continue;
                //if (mListOfPolyhedronParticles[i]->Is(DEMFlags::POLYHEDRON_SKIN)) continue;
                mListOfPolyhedronParticles[i]->mNeighbourElements.push_back(p_polyhedron_neighbour_particle);
                std::vector<PolyhedronParticle*>& neighbours_of_this_neighbour_for_this_thread = thread_maps_of_connectivities[OpenMPUtils::ThisThread()][p_polyhedron_neighbour_particle];
                neighbours_of_this_neighbour_for_this_thread.push_back(mListOfPolyhedronParticles[i]);
            }
            this->GetResultsPoly()[i].clear();
            this->GetResultsDistancesPoly()[i].clear();
        }

        // the next loop ensures consistency in neighbourhood (if A is neighbour of B, B must be neighbour of A)
        #pragma omp parallel for schedule(dynamic, 100)
        for (int i = 0; i < number_of_particles; i++) {
            auto& current_neighbours = mListOfPolyhedronParticles[i]->mNeighbourElements;
            std::vector<PolyhedronParticle*> neighbours_to_add;
            for (size_t k = 0; k < thread_maps_of_connectivities.size(); k++){
                ConnectivitiesMap::iterator it = thread_maps_of_connectivities[k].find(mListOfPolyhedronParticles[i]);
                if (it != thread_maps_of_connectivities[k].end()) {
                    neighbours_to_add.insert(neighbours_to_add.end(), it->second.begin(), it->second.end());
                }
            }
            for (size_t l = 0; l < neighbours_to_add.size(); l++) {
                bool found = false;
                for (size_t m = 0; m < current_neighbours.size(); m++){
                    if (neighbours_to_add[l] == current_neighbours[m]) {
                        found = true;
                        break;
                    }
                }
                if ( found == false ) {
                    current_neighbours.push_back(neighbours_to_add[l]);
                }
            }
        }
        KRATOS_CATCH("")
    }

    void ContactExplicitSolverStrategy::ComputePolyhedronNewNeighboursHistoricalData() {

        KRATOS_TRY
        const int number_of_particles = (int) mListOfPolyhedronParticles.size();

        #pragma omp parallel
        {
            DenseVector<int> temp_neighbours_ids;
            std::vector<array_1d<double, 3> > temp_neighbour_elastic_contact_forces;

            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                mListOfPolyhedronParticles[i]->ComputeNewNeighboursHistoricalData(temp_neighbours_ids, temp_neighbour_elastic_contact_forces);
            }
        }

        KRATOS_CATCH("")
    }

    bool ContactExplicitSolverStrategy::ElementExists(const PolyhedronContactElementContainer& elements, int this_element_id, int neighbour_element_id) {
        KRATOS_TRY
        
        auto it = std::find_if(elements.begin(), elements.end(),
            [this_element_id, neighbour_element_id](const PolyhedronContactElement::Pointer& elem) {
                if (elem->GetPolyElement1()->Id() == this_element_id &&
                    elem->GetPolyElement2()->Id() == neighbour_element_id) {

                    #pragma omp critical
                    {
                        elem->SetDeleteFlag(false); 
                    }

                    return true;
                }
                return false;
            });

        return it != elements.end();

        KRATOS_CATCH("")
    }

    void ContactExplicitSolverStrategy::CreatePolyhedronContactElements() {
        KRATOS_TRY

        std::string ElementName;
        ElementName = std::string("PolyhedronContactElement");
        PolyhedronContactElement rReferenceElement(ElementName);

        //Here we are going to create contact elements when we are on a target particle and we see a neighbor whose id is higher than ours.
        //When our particle has a higher ID than the neighbor we also create a pointer to the (previously) created contact element.
        //We proceed in this way because we want to have the pointers to contact elements in a list in the same order as the initial elements order.

        const int number_of_particles = (int) mListOfPolyhedronParticles.size();
        int used_contacts_counter = 0;
        const ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();

        #pragma omp parallel
        {
 
            #pragma omp for
            for(int i=0; i<(int) mPolyhedronContactElements.size(); i++) {
                mPolyhedronContactElements[i]->SetDeleteFlag(true);;
            }
            
            int private_counter = 0;
            PolyhedronContactElement::Pointer p_new_contact_element;
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                bool add_new_bond = true;
                std::vector<PolyhedronParticle*>& neighbour_elements = mListOfPolyhedronParticles[i]->mNeighbourElements;
                unsigned int neighbors_size = mListOfPolyhedronParticles[i]->mNeighbourElements.size();

                for (unsigned int j = 0; j < neighbors_size; j++) {
                    PolyhedronParticle* neighbour_element = dynamic_cast<PolyhedronParticle*> (neighbour_elements[j]);
                    if (neighbour_element == NULL) continue; //The initial neighbor was deleted at some point in time!!
                    if (mListOfPolyhedronParticles[i]->Id() > neighbour_element->Id()) continue;

                    if (!ElementExists(mPolyhedronContactElements, mListOfPolyhedronParticles[i]->Id(), neighbour_element->Id())){
                        const Properties::Pointer& properties = mListOfPolyhedronParticles[i]->pGetProperties();
                        p_new_contact_element = rReferenceElement.Create(used_contacts_counter + 1, mListOfPolyhedronParticles[i], neighbour_element, properties);

                        #pragma omp critical
                        {
                            mPolyhedronContactElements.push_back(p_new_contact_element);
                            mPolyhedronContactElements[used_contacts_counter]->Initialize(r_process_info);
                            used_contacts_counter++;
                        }
                    }
                    
                    /*
                    #pragma omp critical
                    {
                        if (used_contacts_counter < (int) mPolyhedronContactElements.size()) {
                            add_new_bond = false;
                            private_counter = used_contacts_counter;
                            used_contacts_counter++;
                        }
                    }
                    if (!add_new_bond) {
                        PolyhedronContactElement::Pointer& p_old_contact_element = mPolyhedronContactElements[private_counter];
                        p_old_contact_element->SetPolyElement1(mListOfPolyhedronParticles[i]);
                        p_old_contact_element->SetPolyElement2(neighbour_element);
                        p_old_contact_element->SetId(used_contacts_counter);
                        //p_old_contact_element->SetProperties(mListOfPolyhedronParticles[i]->pGetProperties());
                        p_old_contact_element->Initialize(r_process_info);
                    } else {
                        const Properties::Pointer& properties = mListOfPolyhedronParticles[i]->pGetProperties();
                        p_new_contact_element = rReferenceElement.Create(used_contacts_counter + 1, mListOfPolyhedronParticles[i], neighbour_element, properties);

                        #pragma omp critical
                        {
                            mPolyhedronContactElements.push_back(p_new_contact_element);
                            mPolyhedronContactElements[used_contacts_counter]->Initialize(r_process_info);
                            used_contacts_counter++;
                        }
                    }*/

                }
            }

            /*
            #pragma omp single
            {
                if ((int) mPolyhedronContactElements.size() > used_contacts_counter) {
                    mPolyhedronContactElements.erase(mPolyhedronContactElements.begin() + used_contacts_counter, mPolyhedronContactElements.end());
                }
            }*/

            //delete useless elements
            #pragma omp for
            for (int k = 0; k < (int) mPolyhedronContactElements.size(); k++) {
                if (mPolyhedronContactElements[k]->GetDeleteFlag()) {
                    #pragma omp critical
                    {
                        delete mPolyhedronContactElements[k];
                        mPolyhedronContactElements[k] = nullptr;
                    }
                }
            }

            #pragma omp single
            {
                mPolyhedronContactElements.erase(std::remove(mPolyhedronContactElements.begin(), mPolyhedronContactElements.end(), nullptr), mPolyhedronContactElements.end());
            }

            //Renumbering the Id's of the bonds to make them unique and consecutive (otherwise the Id's are repeated)
            #pragma omp for
            for(int i=0; i<(int) mPolyhedronContactElements.size(); i++) {
                mPolyhedronContactElements[i]->SetId(i+1);
            }

        } //#pragma omp parallel
        KRATOS_CATCH("")
    } //CreateContactElements

    void ContactExplicitSolverStrategy::InitializePolyhedronContactElements() {

        KRATOS_TRY

        const ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();

        #pragma omp for
        for(int i=0; i<(int) mPolyhedronContactElements.size(); i++) {
            mPolyhedronContactElements[i]->Initialize(r_process_info);
        }

        KRATOS_CATCH("")
    }

    //TODO: update this function
    void ContactExplicitSolverStrategy::SetSearchRadiiOnAllPolyhedronParticles(ModelPart& polyhedron_model_part, const double added_search_distance, const double amplification) {
        KRATOS_TRY
        const int number_of_elements = polyhedron_model_part.GetCommunicator().LocalMesh().NumberOfElements();
        #pragma omp parallel for
        for (int i = 0; i < number_of_elements; i++) {
            mListOfPolyhedronParticles[i]->SetSearchRadius(mListOfPolyhedronParticles[i]->GetRadius());
        }

        KRATOS_CATCH("")
    }

    void ContactExplicitSolverStrategy::ForceOperations(ModelPart& r_model_part, ModelPart& r_polyhedron_model_part) {
        KRATOS_TRY

        //TODO: those force operations are still based on particle loop
        GetForce(); // Basically only calls CalculateRightHandSide()
        GetClustersForce();
        GetRigidBodyElementsForce();

        // force operations based on contact loop
        GetPolyhedronForce();

        if (r_model_part.GetProcessInfo()[COMPUTE_FEM_RESULTS_OPTION]) {
            CalculateNodalPressuresAndStressesOnWalls();
        }

        // Synchronize (should be just FORCE and TORQUE)
        SynchronizeRHS(r_model_part);
        SynchronizeRHS(r_polyhedron_model_part);

        KRATOS_CATCH("")
    }//ForceOperations;

    void ContactExplicitSolverStrategy::GetPolyhedronForce() {
        KRATOS_TRY

        ProcessInfo& r_process_info = GetModelPart().GetProcessInfo();
        double dt = r_process_info[DELTA_TIME];
        const array_1d<double, 3>& gravity = r_process_info[GRAVITY];

        ModelPart& polyhedron_model_part = GetPolyhedronModelPart();
        ElementsArrayType& pElements = polyhedron_model_part.GetCommunicator().LocalMesh().Elements();
        const int number_of_polyhedron_elements = pElements.size();

        #pragma omp parallel for schedule(dynamic, 100)
        for (int k = 0; k < number_of_polyhedron_elements; k++) {

            ElementsArrayType::iterator it = pElements.ptr_begin() + k;
            PolyhedronParticle& polyhedron_element = dynamic_cast<Kratos::PolyhedronParticle&> (*it);
            polyhedron_element.GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES).clear();
            polyhedron_element.GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT).clear();
            polyhedron_element.ComputeExternalForces(gravity);

        } 

        /*
        ElementsArrayType& pPolyhedronContactElements = GetAllElements(*mpPolyhedron_contact_model_part);
        const int number_of_polyhedron_contact_elements = (int) pPolyhedronContactElements.size();

        #pragma omp parallel for schedule(dynamic, 100)
        for (int i = 0; i < number_of_polyhedron_contact_elements; i++) {
            mPolyhedronContactElements[i]->CalculateRightHandSide(r_process_info, dt, gravity);
        }*/

        #pragma omp parallel for schedule(dynamic, 100)
        for(int i=0; i<(int) mPolyhedronContactElements.size(); i++) {
            mPolyhedronContactElements[i]->CalculateRightHandSide(r_process_info, dt, gravity);;
        }

        KRATOS_CATCH("")
    }

    void ContactExplicitSolverStrategy::FinalizeSolutionStep() {
        BaseType::FinalizeSolutionStep();
        FinalizeSolutionStepFEM();
    }

    void ContactExplicitSolverStrategy::FinalizeSolutionStepFEM() {
        KRATOS_TRY

        ConditionsArrayType& pConditions = GetFemModelPart().GetCommunicator().LocalMesh().Conditions();
        const ProcessInfo& r_process_info = GetFemModelPart().GetProcessInfo();

        block_for_each(pConditions, [&r_process_info](ModelPart::ConditionType& rCondition){
            rCondition.FinalizeSolutionStep(r_process_info);
        });
        KRATOS_CATCH("")
    }
}
