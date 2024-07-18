/////////////////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
// Date: July 2024
/////////////////////////////////////////////////////////////

#include "explicit_solver_contact_strategy.h"
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"

namespace Kratos {

    void ContactExplicitSolverStrategy::Initialize() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();
        ModelPart& fem_model_part = GetFemModelPart();
        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();
        //SendProcessInfoToClustersModelPart();

        if (r_model_part.GetCommunicator().MyPID() == 0) {
            KRATOS_INFO("DEM") << "------------------CONTACT-BASED EXPLICIT SOLVER STRATEGY---------------------" << "\n" << std::endl;
        }

        mNumberOfThreads = ParallelUtilities::GetNumThreads();
        DisplayThreadInfo();

        RebuildListOfSphericParticles <PolyhedronParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfPolyhedronParticles);
        RebuildListOfSphericParticles <PolyhedronParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostPolyhedronParticles);

        mSearchControlVector.resize(mNumberOfThreads);
        for (int i = 0; i < mNumberOfThreads; i++) mSearchControlVector[i] = 0;

        PropertiesProxiesManager().CreatePropertiesProxies(r_model_part, *mpInlet_model_part, *mpCluster_model_part);

        RepairPointersToNormalProperties(mListOfSphericParticles);
        RepairPointersToNormalProperties(mListOfGhostSphericParticles);

        RebuildPropertiesProxyPointers(mListOfSphericParticles);
        RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);

        GetSearchControl() = r_process_info[SEARCH_CONTROL];

        InitializeDEMElements();
        InitializeFEMElements();
        UpdateMaxIdOfCreatorDestructor();
        InitializeClusters(); // This adds elements to the balls modelpart

        RebuildListOfSphericParticles <PolyhedronParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfPolyhedronParticles);
        RebuildListOfSphericParticles <PolyhedronParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostPolyhedronParticles);

        InitializeSolutionStep();
        ApplyInitialConditions();

        // Search Neighbors with tolerance (after first repartition process)
        SetSearchRadiiOnAllParticles(r_model_part, r_process_info[SEARCH_RADIUS_INCREMENT_FOR_BONDS_CREATION], 1.0);
        SearchNeighbours();
        MeshRepairOperations();
        SearchNeighbours();

        RebuildListOfSphericParticles <PolyhedronParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfPolyhedronParticles);
        RebuildListOfSphericParticles <PolyhedronParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostPolyhedronParticles);

        bool has_mpi = false;
        Check_MPI(has_mpi);

        if (has_mpi) {
            RepairPointersToNormalProperties(mListOfSphericParticles);
            RepairPointersToNormalProperties(mListOfGhostSphericParticles);
        }

        RebuildPropertiesProxyPointers(mListOfSphericParticles);
        RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);

        if (has_mpi) {
            //RebuildListsOfPointersOfEachParticle(); //Serialized pointers are lost, so we rebuild them using Id's
        }

        if (fem_model_part.Nodes().size() > 0) {
            SetSearchRadiiWithFemOnAllParticles(r_model_part, mpDem_model_part->GetProcessInfo()[SEARCH_RADIUS_INCREMENT_FOR_WALLS], 1.0);
            SearchRigidFaceNeighbours();
        }

        if (r_process_info[CONTACT_MESH_OPTION] == 1) {
            CreateContactElements();
            InitializeContactElements();
        }

        r_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(NEIGHBOUR_IDS);
        r_model_part.GetCommunicator().SynchronizeElementalNonHistoricalVariable(NEIGHBOURS_CONTACT_AREAS);

        ComputeNodalArea();

        KRATOS_CATCH("")
    }// Initialize()

    double ContactExplicitSolverStrategy::SolveSolutionStep() {

        KRATOS_TRY

        ModelPart& r_model_part = GetModelPart();

        bool has_mpi = false;
        VariablesList r_modelpart_nodal_variables_list = r_model_part.GetNodalSolutionStepVariablesList();
        if (r_modelpart_nodal_variables_list.Has(PARTITION_INDEX)) has_mpi = true;

        SearchDEMOperations(r_model_part, has_mpi);
        SearchFEMOperations(r_model_part, has_mpi);
        ForceOperations(r_model_part);
        PerformTimeIntegrationOfMotion();

        KRATOS_CATCH("")

        return 0.0;

    }//SolveSolutionStep()

    void ContactExplicitSolverStrategy::SearchDEMOperations(ModelPart& r_model_part, bool has_mpi) {

        ProcessInfo& r_process_info = r_model_part.GetProcessInfo();

        if (r_process_info[SEARCH_CONTROL] == 0) {

            ElementsArrayType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();
            int some_bond_is_broken = 0;

            block_for_each(rElements, [&](ModelPart::ElementType& rElement) {

                PolyhedronParticle& r_sphere = dynamic_cast<PolyhedronParticle&>(rElement);

                for (int j=0; j<(int) r_sphere.mContinuumInitialNeighborsSize; j++) {
                    if (r_sphere.mIniNeighbourFailureId[j] != 0) {
                        AtomicAdd(some_bond_is_broken, 1);
                        break;
                    }
                }
            });

            if (some_bond_is_broken > 0) {
                r_process_info[SEARCH_CONTROL] = 1;
                KRATOS_WARNING("DEM") << "From now on, the search is activated because some failure occurred " << std::endl;
            }
        }

        const int time_step = r_process_info[TIME_STEPS];
        const double time = r_process_info[TIME];
        const bool is_time_to_search_neighbours = (time_step + 1) % GetNStepSearch() == 0 && (time_step > 0); //Neighboring search. Every N times.

        if (r_process_info[SEARCH_CONTROL] > 0) {

            if (is_time_to_search_neighbours) {

                if (r_process_info[BOUNDING_BOX_OPTION] && time >= r_process_info[BOUNDING_BOX_START_TIME] && time <= r_process_info[BOUNDING_BOX_STOP_TIME]) {
                    BoundingBoxUtility();
                } else {
                    GetParticleCreatorDestructor()->DestroyParticles<SphericParticle>(r_model_part);
                    GetParticleCreatorDestructor()->DestroyContactElements(*mpContact_model_part);
                }

                RebuildListOfSphericParticles <PolyhedronParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfPolyhedronParticles); //These lists are necessary for the loop in SearchNeighbours

                SetSearchRadiiOnAllParticles(r_model_part, r_process_info[SEARCH_RADIUS_INCREMENT], r_process_info[CONTINUUM_SEARCH_RADIUS_AMPLIFICATION_FACTOR]);

                SearchNeighbours(); //the amplification factor has been modified after the first search.

                RebuildListOfSphericParticles <PolyhedronParticle> (r_model_part.GetCommunicator().LocalMesh().Elements(), mListOfPolyhedronParticles); //These lists are necessary because the elements in this partition might have changed.
                RebuildListOfSphericParticles <PolyhedronParticle> (r_model_part.GetCommunicator().GhostMesh().Elements(), mListOfGhostPolyhedronParticles);

                if (has_mpi) {
                    RepairPointersToNormalProperties(mListOfSphericParticles);
                    RepairPointersToNormalProperties(mListOfGhostSphericParticles);
                }

                RebuildPropertiesProxyPointers(mListOfSphericParticles);
                RebuildPropertiesProxyPointers(mListOfGhostSphericParticles);

                ComputeNewNeighboursHistoricalData();

                MarkNewSkinParticles();

                r_process_info[SEARCH_CONTROL] = 2;
            } else {
                r_process_info[SEARCH_CONTROL] = 1;
            }

            if (r_process_info[CONTACT_MESH_OPTION]) {
                CreateContactElements();
                InitializeContactElements();
            }
        }
        //Synch this var.
        r_process_info[SEARCH_CONTROL] = r_model_part.GetCommunicator().GetDataCommunicator().MaxAll(r_process_info[SEARCH_CONTROL]);
    }

    void ExplicitSolverStrategy::CreateContactElements() {
        KRATOS_TRY

        std::string ElementName;
        ElementName = std::string("ParticleContactElement");
        const Element& rReferenceElement = KratosComponents<Element>::Get(ElementName);

        //Here we are going to create contact elements when we are on a target particle and we see a neighbor whose id is higher than ours.
        //We create also a pointer from the node to the element, after creating it.
        //When our particle has a higher ID than the neighbor we also create a pointer to the (previously) created contact element.
        //We proceed in this way because we want to have the pointers to contact elements in a list in the same order as the initial elements order.

        const int number_of_particles = (int) mListOfSphericParticles.size();
        int used_bonds_counter = 0;

        #pragma omp parallel
        {
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                unsigned int neighbors_size = mListOfSphericParticles[i]->mNeighbourElements.size();
                mListOfSphericParticles[i]->mBondElements.resize(neighbors_size);
                for (unsigned int j = 0; j < mListOfSphericParticles[i]->mBondElements.size(); j++) {
                    mListOfSphericParticles[i]->mBondElements[j] = NULL;
                }
            }

            int private_counter = 0;
            Element::Pointer p_new_contact_element;
            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                bool add_new_bond = true;
                std::vector<SphericParticle*>& neighbour_elements = mListOfSphericParticles[i]->mNeighbourElements;
                unsigned int neighbors_size = mListOfSphericParticles[i]->mNeighbourElements.size();

                for (unsigned int j = 0; j < neighbors_size; j++) {
                    SphericParticle* neighbour_element = dynamic_cast<SphericParticle*> (neighbour_elements[j]);
                    if (neighbour_element == NULL) continue; //The initial neighbor was deleted at some point in time!!
                    if (mListOfSphericParticles[i]->Id() > neighbour_element->Id()) continue;

                    #pragma omp critical
                    {
                        if (used_bonds_counter < (int) (*mpContact_model_part).Elements().size()) {
                            add_new_bond = false;
                            private_counter = used_bonds_counter;
                            used_bonds_counter++;
                        }
                    }
                    if (!add_new_bond) {
                        Element::Pointer& p_old_contact_element = (*mpContact_model_part).Elements().GetContainer()[private_counter];
                        p_old_contact_element->GetGeometry()(0) = mListOfSphericParticles[i]->GetGeometry()(0);
                        p_old_contact_element->GetGeometry()(1) = neighbour_element->GetGeometry()(0);
                        p_old_contact_element->SetId(used_bonds_counter);
                        p_old_contact_element->SetProperties(mListOfSphericParticles[i]->pGetProperties());
                        ParticleContactElement* p_bond = dynamic_cast<ParticleContactElement*> (p_old_contact_element.get());
                        mListOfSphericParticles[i]->mBondElements[j] = p_bond;
                    } else {
                        Geometry<Node >::PointsArrayType NodeArray(2);
                        NodeArray.GetContainer()[0] = mListOfSphericParticles[i]->GetGeometry()(0);
                        NodeArray.GetContainer()[1] = neighbour_element->GetGeometry()(0);
                        const Properties::Pointer& properties = mListOfSphericParticles[i]->pGetProperties();
                        p_new_contact_element = rReferenceElement.Create(used_bonds_counter + 1, NodeArray, properties);

                        #pragma omp critical
                        {
                            (*mpContact_model_part).Elements().push_back(p_new_contact_element);
                            used_bonds_counter++;
                        }
                        ParticleContactElement* p_bond = dynamic_cast<ParticleContactElement*> (p_new_contact_element.get());
                        mListOfSphericParticles[i]->mBondElements[j] = p_bond;
                    }

                }
            }

            #pragma omp single
            {
                if ((int) (*mpContact_model_part).Elements().size() > used_bonds_counter) {
                    (*mpContact_model_part).Elements().erase((*mpContact_model_part).Elements().ptr_begin() + used_bonds_counter, (*mpContact_model_part).Elements().ptr_end());
                }
            }

            #pragma omp for
            for (int i = 0; i < number_of_particles; i++) {
                std::vector<SphericParticle*>& neighbour_elements = mListOfSphericParticles[i]->mNeighbourElements;
                unsigned int neighbors_size = mListOfSphericParticles[i]->mNeighbourElements.size();

                for (unsigned int j = 0; j < neighbors_size; j++) {
                    SphericParticle* neighbour_element = dynamic_cast<SphericParticle*> (neighbour_elements[j]);
                    if (neighbour_element == NULL) continue; //The initial neighbor was deleted at some point in time!!
                    //ATTENTION: Ghost nodes do not have mContinuumIniNeighbourElements in general, so this bond will remain as NULL!!
                    if (mListOfSphericParticles[i]->Id() < neighbour_element->Id()) continue;
                    //In all functions using mBondElements we must check that this bond is not used.

                    for (unsigned int k = 0; k < neighbour_element->mNeighbourElements.size(); k++) {
                        //ATTENTION: Ghost nodes do not have mContinuumIniNeighbourElements in general, so this bond will remain as NULL!!
                        //In all functions using mBondElements we must check that this bond is not used.
                        if (neighbour_element->mNeighbourElements[k] == NULL) continue; //The initial neighbor was deleted at some point in time!!
                        if (neighbour_element->mNeighbourElements[k]->Id() == mListOfSphericParticles[i]->Id()) {
                            ParticleContactElement* bond = neighbour_element->mBondElements[k];
                            mListOfSphericParticles[i]->mBondElements[j] = bond;
                            break;
                        }
                    }
                }
            }

            //Renumbering the Id's of the bonds to make them unique and consecutive (otherwise the Id's are repeated)
            #pragma omp for
            for(int i=0; i<(int)(*mpContact_model_part).Elements().size(); i++) {
                (*mpContact_model_part).Elements().GetContainer()[i]->SetId(i+1);
            }

        } //#pragma omp parallel
        KRATOS_CATCH("")
    } //CreateContactElements

    void ContactExplicitSolverStrategy::SetSearchRadiiOnAllParticles(ModelPart& r_model_part, const double added_search_distance, const double amplification) {
        KRATOS_TRY
        const int number_of_elements = r_model_part.GetCommunicator().LocalMesh().NumberOfElements();
        if (GetDeltaOption() == 3){
            // In this case, the parameter "added_search_distance" is actually a multiplier for getting the added_search_distance
            const double search_radius_multiplier = added_search_distance;
            #pragma omp parallel for
            for (int i = 0; i < number_of_elements; i++) {
                mListOfPolyhedronParticles[i]->SetSearchRadius(amplification * mListOfPolyhedronParticles[i]->mLocalRadiusAmplificationFactor * ((1 + search_radius_multiplier) * mListOfPolyhedronParticles[i]->GetRadius()));
            }
        }
        else{
            #pragma omp parallel for
            for (int i = 0; i < number_of_elements; i++) {
                mListOfPolyhedronParticles[i]->SetSearchRadius(amplification * mListOfPolyhedronParticles[i]->mLocalRadiusAmplificationFactor * (added_search_distance + mListOfPolyhedronParticles[i]->GetRadius()));
            }
        }
        KRATOS_CATCH("")
    }

    void ContactExplicitSolverStrategy::MeshRepairOperations() {

        KRATOS_TRY

        const int number_of_particles = (int) mListOfPolyhedronParticles.size();
        int particle_counter = 0.0;

        #pragma omp parallel for
        for (int i = 0; i < number_of_particles; i++) {
            bool result = mListOfPolyhedronParticles[i]->OverlappedParticleRemoval();

            if (result == true) {particle_counter += 1;}
        }

        GetModelPart().GetCommunicator().SynchronizeElementalFlags();

        //KRATOS_WARNING("DEM") << "Mesh repair complete. In MPI node " <<GetModelPart().GetCommunicator().MyPID()<<". "<< particle_counter << " particles were removed. " << "\n" << std::endl;
        int total_spheres_removed = GetModelPart().GetCommunicator().GetDataCommunicator().SumAll(particle_counter);

        if(GetModelPart().GetCommunicator().MyPID() == 0 && total_spheres_removed) {
            KRATOS_WARNING("DEM") << "A total of "<<total_spheres_removed<<" spheres were removed due to excessive overlapping." << std::endl;
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
