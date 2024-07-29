/////////////////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
// Date: July 2024
/////////////////////////////////////////////////////////////


#if !defined(KRATOS_CONTACT_EXPLICIT_SOLVER_STRATEGY)
#define  KRATOS_CONTACT_EXPLICIT_SOLVER_STRATEGY
#include "custom_strategies/strategies/explicit_solver_strategy.h"
#include "custom_elements/polyhedron_particle.h"
#include "custom_utilities/create_and_destroy.h"

namespace Kratos {

    class ContactExplicitSolverSettings {
    public:
        KRATOS_CLASS_POINTER_DEFINITION(ContactExplicitSolverSettings);

        ContactExplicitSolverSettings() {
        }

        ~ContactExplicitSolverSettings() {
        }
        ModelPart* r_model_part;
        ModelPart* contact_model_part;
        ModelPart* fem_model_part;
        ModelPart* cluster_model_part;
        ModelPart* inlet_model_part;
        ModelPart* polyhedron_model_part;

    };
    
    class KRATOS_API(DEM_APPLICATION) ContactExplicitSolverStrategy : public ExplicitSolverStrategy {
    public:

        typedef ExplicitSolverStrategy BaseType;
        typedef BaseType::NodesArrayType NodesArrayType;
        typedef BaseType::ElementsArrayType ElementsArrayType;
        typedef BaseType::ElementsIterator ElementsIterator;
        typedef BaseType::ConditionsArrayType ConditionsArrayType;
        typedef GlobalPointersVector<Element> ParticleWeakVectorType;
        typedef GlobalPointersVector<Element >::iterator ParticleWeakIteratorType;
        typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;

        using BaseType::mpInlet_model_part;
        using BaseType::mpCluster_model_part;
        using BaseType::mpContact_model_part;
        using BaseType::GetModelPart;
        using BaseType::GetFemModelPart;
        using BaseType::mNumberOfThreads;
        using BaseType::mListOfSphericParticles;
        using BaseType::mListOfGhostSphericParticles;
        using BaseType::SearchNeighbours;

        using PolyhedronContactElementContainer = std::vector<PolyhedronContactElement::Pointer>;

        //Pointer definition of ExplicitSolverStrategy
        KRATOS_CLASS_POINTER_DEFINITION(ContactExplicitSolverStrategy);

        //Default constructor.
        ContactExplicitSolverStrategy() {
        }

        ContactExplicitSolverStrategy(
                ContactExplicitSolverSettings& settings,
                const double max_delta_time,
                const int n_step_search,
                const double safety_factor,
                const int delta_option,
                ParticleCreatorDestructor::Pointer p_creator_destructor,
                DEM_FEM_Search::Pointer p_dem_fem_search,
                SpatialSearch::Pointer pSpSearch,
                Parameters strategy_parameters)
        {
            mParameters = strategy_parameters;
            mDeltaOption = delta_option;
            mpParticleCreatorDestructor = p_creator_destructor;
            mpDemFemSearch = p_dem_fem_search;
            mpSpSearch = pSpSearch;

            //Also checks old flag name for backward compatibility issues.
            if(mParameters["do_search_dem_neighbours"].GetBool()) {
                mDoSearchNeighbourElements = true;
            } else mDoSearchNeighbourElements = false;
            p_creator_destructor->SetDoSearchNeighbourElements(mDoSearchNeighbourElements);

            if(mParameters["do_search_fem_neighbours"].GetBool()) mDoSearchNeighbourFEMElements = true;
            else mDoSearchNeighbourFEMElements = false;

            mMaxTimeStep = max_delta_time;
            mNStepSearch = n_step_search;
            mSafetyFactor = safety_factor;

            mpDem_model_part = &(*(settings.r_model_part));
            KRATOS_ERROR_IF(mpDem_model_part == NULL) << "Undefined settings.r_model_part in ContactExplicitSolverStrategy constructor" << std::endl;

            mpContact_model_part = &(*(settings.contact_model_part));
            KRATOS_ERROR_IF(mpContact_model_part == NULL) << "Undefined settings.contact_model_part in ContactExplicitSolverStrategy constructor" << std::endl;

            mpFem_model_part = &(*(settings.fem_model_part));
            KRATOS_ERROR_IF(mpFem_model_part == NULL) << "Undefined settings.fem_model_part in ContactExplicitSolverStrategy constructor" << std::endl;

            mpCluster_model_part = &(*(settings.cluster_model_part));
            KRATOS_ERROR_IF(mpCluster_model_part == NULL) << "Undefined settings.cluster_model_part in ContactExplicitSolverStrategy constructor" << std::endl;

            mpInlet_model_part = &(*(settings.inlet_model_part));
            KRATOS_ERROR_IF(mpInlet_model_part == NULL) << "Undefined settings.inlet_model_part in ContactExplicitSolverStrategy constructor" << std::endl;

            mpPolyhedron_model_part = &(*(settings.polyhedron_model_part));
            KRATOS_ERROR_IF(mpPolyhedron_model_part == NULL) << "Undefined settings.polyhedron_model_part in ContactExplicitSolverStrategy constructor" << std::endl;

            if(mParameters["RemoveBallsInitiallyTouchingWalls"].GetBool()) mRemoveBallsInitiallyTouchingWallsOption = true;
            else mRemoveBallsInitiallyTouchingWallsOption = false;

        }

        //Destructor.
        virtual ~ContactExplicitSolverStrategy() {}

        void Initialize() override;
        void RebuildPropertiesProxyPointersForPolyhedron(std::vector<PolyhedronParticle*>& rCustomListOfPolyhedronParticles);
        virtual void SearchPolyhedronOperations(ModelPart& r_model_part, ModelPart& polyhedron_model_part, bool has_mpi);
        virtual void SearchPolyhedronNeighbours();
        virtual void ComputePolyhedronNewNeighboursHistoricalData();
        void CreatePolyhedronContactElements();
        void CreateContactElements() override;
        double ComputeCoordinationNumber(double& standard_dev) override;
        virtual void RepairPointersToNormalPropertiesOfPolyhedron(std::vector<PolyhedronParticle*>& rCustomListOfPolyhedronParticles);
        void InitializePolyhedrons();
        virtual void InitializePolyhedronContactElements();
        void UpdateMaxIdOfCreatorDestructor() override;
        void ApplyPrescribedBoundaryConditions() override;
        void RebuildListOfPolyhedronParticles() {
            RebuildListOfSphericParticles<PolyhedronParticle>(GetModelPart().GetCommunicator().LocalMesh().Elements(), mListOfPolyhedronParticles);
        }
        virtual void SetSearchRadiiOnAllPolyhedronParticles(ModelPart& polyhedron_model_part, const double added_search_distance, const double amplification);
        virtual void MeshRepairOperations();
        void InitializeSolutionStep() override;
        void ApplyInitialConditions() override;
        double SolveSolutionStep() override;
        void FinalizeSolutionStep() override;
        void FinalizeSolutionStepFEM();

        ModelPart& GetPolyhedronModelPart() { return (*mpPolyhedron_model_part);}

        virtual void Add_As_Own(ModelPart& r_model_part, ModelPart& mcontacts_model_part, ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator, Element::Pointer p_contact_element) {
            KRATOS_TRY
            mcontacts_model_part.Elements().push_back(p_contact_element);
            KRATOS_CATCH("")
        }

        virtual void Add_As_Local(ModelPart& r_model_part, ModelPart& mcontacts_model_part, ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator, Element::Pointer p_contact_element) {
        }

        virtual void Add_As_Ghost(ModelPart& r_model_part, ModelPart& mcontacts_model_part, ParticleWeakIteratorType_ptr continuum_ini_neighbour_iterator, Element::Pointer p_contact_element) {
        }

        virtual void Sort_Contact_Modelpart(ModelPart& mcontacts_model_part) {
        }

        virtual void Reassign_Ids(ModelPart& mcontacts_model_part) {
        }

        virtual ElementsArrayType& GetElements(ModelPart& r_model_part) override {
            return r_model_part.GetCommunicator().LocalMesh().Elements();
        }

    protected:

        int mFixSwitch;
        std::vector<PolyhedronParticle*> mListOfPolyhedronParticles;
        std::vector<PolyhedronParticle*> mListOfGhostPolyhedronParticles;
        DenseVector<int> mSearchControlVector;
        ModelPart *mpPolyhedron_model_part;
        PolyhedronContactElementContainer mContactElements;

    }; //Class ContactExplicitSolverStrategy

} //namespace Kratos

#endif //KRATOS_CONTACT_EXPLICIT_SOLVER_STRATEGY  defined
