/////////////////////////////////////////////////////////////
// Main author: Chengshun Shang (CIMNE)
// Email: cshang@cimne.upc.edu; chengshun.shang1996@gmail.com
// Date: Sept 2024
/////////////////////////////////////////////////////////////

#if !defined(KRATOS_ITERATIVE_SOLVER_CONTACT_STRATEGY)
#define  KRATOS_ITERATIVE_SOLVER_CONTACT_STRATEGY
#include "custom_strategies/strategies/explicit_solver_contact_strategy.h"

namespace Kratos
{
  class ContactIterativeSolverStrategy: public ContactExplicitSolverStrategy
  {
      public:

      typedef ContactExplicitSolverStrategy  BaseType;

      typedef BaseType::NodesArrayType                             NodesArrayType;
      typedef BaseType::ElementsArrayType                          ElementsArrayType;
      typedef BaseType::ElementsIterator                           ElementsIterator;
      typedef BaseType::ConditionsArrayType                        ConditionsArrayType;

      typedef GlobalPointersVector<Element> ParticleWeakVectorType;
      typedef GlobalPointersVector<Element >::iterator ParticleWeakIteratorType;
      typedef ParticleWeakVectorType::ptr_iterator ParticleWeakIteratorType_ptr;

      /// Pointer definition
      KRATOS_CLASS_POINTER_DEFINITION(ContactIterativeSolverStrategy);

      /// Default constructor.
      ContactIterativeSolverStrategy(){}

      ContactIterativeSolverStrategy(
                             ContactExplicitSolverSettings& settings,
                             const double max_delta_time,
                             const double n_step_search,
                             const double safety_factor,
                             const int delta_option,
                             ParticleCreatorDestructor::Pointer p_creator_destructor,
                             DEM_FEM_Search::Pointer p_dem_fem_search,
                             SpatialSearch::Pointer pSpSearch,
                             Parameters strategy_parameters)
      :ContactExplicitSolverStrategy(settings, max_delta_time, n_step_search, safety_factor, delta_option, p_creator_destructor, p_dem_fem_search, pSpSearch, strategy_parameters)
      {
          
      }

      /// Destructor.
      virtual ~ContactIterativeSolverStrategy()
      {

      }


      virtual void Initialize() override
      {

        KRATOS_TRY
        ModelPart& r_model_part = BaseType::GetModelPart();
        ModelPart& r_polyhedron_model_part = BaseType::GetPolyhedronModelPart();

        BaseType::Initialize();
        BaseType::InitializeSolutionStep();
        BaseType::ForceOperations(r_model_part, r_polyhedron_model_part);
        BaseType::FinalizeSolutionStep();

        KRATOS_CATCH("")
      }// Initialize()


      void SchemeInitialize()
      {
        PerformTimeIntegrationOfMotion(0);
      }


      void SchemePredict()
      {
        PerformTimeIntegrationOfMotion(1);
      }

      void SchemeCorrect()
      {
        PerformTimeIntegrationOfMotion(2);
      }

      virtual double SolveSolutionStep() override
      {

        KRATOS_TRY

        ModelPart& r_model_part = BaseType::GetModelPart();
        ModelPart& r_polyhedron_model_part = BaseType::GetPolyhedronModelPart();

        bool has_mpi = false;
        VariablesList r_modelpart_nodal_variables_list = r_model_part.GetNodalSolutionStepVariablesList();
        if (r_modelpart_nodal_variables_list.Has(PARTITION_INDEX)) has_mpi = true;

        SchemePredict();
        BaseType::SearchDEMOperations(r_model_part, has_mpi);
        BaseType::SearchFEMOperations(r_model_part, has_mpi);
        BaseType::SearchPolyhedronOperations(r_model_part, r_polyhedron_model_part, has_mpi);
        BaseType::ForceOperations(r_model_part, r_polyhedron_model_part);
        SchemeCorrect();

        return 0.00;

        KRATOS_CATCH("")

      }//SolveSolutionStep()

  };//ClassIterativeSolverContactStrategy

}  // namespace Kratos.

#endif // KRATOS_ITERATIVE_SOLVER_CONTACT_STRATEGY  defined
