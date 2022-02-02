//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license:
// kratos/license.txt
//
//  Main authors:    Reza Najian Asl
//

#if !defined(KRATOS_HELMHOLTZ_VEC_STRATEGY)
#define KRATOS_HELMHOLTZ_VEC_STRATEGY

/* System includes */

/* External includes */

/* Project includes */
#include "helmholtz_application_variables.h"
#include "containers/model.h"
#include "solving_strategies/builder_and_solvers/residualbased_block_builder_and_solver.h"
#include "solving_strategies/schemes/residualbased_incrementalupdate_static_scheme.h"
#include "solving_strategies/strategies/residualbased_linear_strategy.h"
#include "utilities/variable_utils.h"

namespace Kratos {

/**@name Kratos Globals */
/*@{ */

/*@} */
/**@name Type Definitions */
/*@{ */

/*@} */

/**@name  Enum's */
/*@{ */

/*@} */
/**@name  Functions */
/*@{ */

/*@} */
/**@name Kratos Classes */
/*@{ */

/// Short class definition.
/**   Detail class definition.
 */
template <class TSparseSpace, class TDenseSpace, class TLinearSolver>
class HelmholtzVecStrategy
    : public ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> {
public:
  /**@name Type Definitions */
  /*@{ */

  /** Counted pointer of ClassName */
  KRATOS_CLASS_POINTER_DEFINITION(HelmholtzVecStrategy);

  typedef ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver> BaseType;
  typedef typename BaseType::TBuilderAndSolverType TBuilderAndSolverType;
  typedef Scheme<TSparseSpace, TDenseSpace> SchemeType;

  /*@} */
  /**@name Life Cycle
   */
  /*@{ */

  /** Constructor.
   */
  HelmholtzVecStrategy(ModelPart &model_part,
                               typename TLinearSolver::Pointer pNewLinearSolver,
                               bool ReformDofSetAtEachStep = false,
                               bool ComputeReactions = false,
                               int EchoLevel = 0,
                               const double PoissonRatio = 0.3)
      : ImplicitSolvingStrategy<TSparseSpace, TDenseSpace, TLinearSolver>(model_part) {
    KRATOS_TRY

    mreform_dof_set_at_each_step = ReformDofSetAtEachStep;
    mcompute_reactions = ComputeReactions;
    mecho_level = EchoLevel;
    bool calculate_norm_dx_flag = false;

    typename SchemeType::Pointer pscheme = typename SchemeType::Pointer(
        new ResidualBasedIncrementalUpdateStaticScheme<TSparseSpace,
                                                       TDenseSpace>());

    mpbulider_and_solver = typename TBuilderAndSolverType::Pointer(
        new ResidualBasedBlockBuilderAndSolver<TSparseSpace, TDenseSpace,
                                               TLinearSolver>(
            pNewLinearSolver));

    mstrategy = typename BaseType::Pointer(new ResidualBasedLinearStrategy<TSparseSpace, TDenseSpace,TLinearSolver>(
        BaseType::GetModelPart(),
        pscheme,
        mpbulider_and_solver,
        mcompute_reactions,
        mreform_dof_set_at_each_step,
        calculate_norm_dx_flag));

    mstrategy->SetEchoLevel(mecho_level);

    KRATOS_CATCH("")
  }

  virtual ~HelmholtzVecStrategy()
  {

  }

  void Initialize() override {}

  double Solve() override {
    KRATOS_TRY;

    VariableUtils().UpdateCurrentToInitialConfiguration(
        BaseType::GetModelPart().GetCommunicator().LocalMesh().Nodes());

    // Solve for the mesh movement
    mstrategy->Solve();

    // Clearing the system if needed
    if (mreform_dof_set_at_each_step == true)
      mstrategy->Clear();

    return 0.0;

    KRATOS_CATCH("");
  }

  /*@} */
  /**@name Operators
   */
  /*@{ */

  /*@} */
  /**@name Operations */
  /*@{ */

  /*@} */
  /**@name Access */
  /*@{ */

  /*@} */
  /**@name Inquiry */
  /*@{ */

  /*@} */
  /**@name Friends */
  /*@{ */

  /*@} */

protected:
  /**@name Protected static Member Variables */
  /*@{ */

  /*@} */
  /**@name Protected member Variables */
  /*@{ */

  /*@} */
  /**@name Protected Operators*/
  /*@{ */

  /*@} */
  /**@name Protected Operations*/
  /*@{ */

  /*@} */
  /**@name Protected  Access */
  /*@{ */

  /*@} */
  /**@name Protected Inquiry */
  /*@{ */

  /*@} */
  /**@name Protected LifeCycle */
  /*@{ */

  /*@} */

private:
  /**@name Static Member Variables */
  /*@{ */

  /*@} */
  /**@name Member Variables */
  /*@{ */

  typename BaseType::Pointer mstrategy;
  typename TBuilderAndSolverType::Pointer mpbulider_and_solver;

  int mecho_level;
  bool mreform_dof_set_at_each_step;
  bool mcompute_reactions;

  /*@} */
  /**@name Private Operators*/
  /*@{ */

  /*@} */
  /**@name Private Operations*/
  /*@{ */
  /*@} */
  /**@name Private  Access */
  /*@{ */

  /*@} */
  /**@name Private Inquiry */
  /*@{ */

  /*@} */
  /**@name Un accessible methods */
  /*@{ */

  /** Copy constructor.
   */
  HelmholtzVecStrategy(const HelmholtzVecStrategy &Other);

  /*@} */

}; /* Class HelmholtzVecStrategy */

/*@} */

/**@name Type Definitions */
/*@{ */

/*@} */
}
/* namespace Kratos.*/

#endif /* KRATOS_HELMHOLTZ_VEC_STRATEGY  defined */