// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "assemble/bcl_assemble_sse_factory_mc.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_sse_size.h"
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "assemble/bcl_assemble_sse_pool_insert_coil_into_sse.h"
#include "assemble/bcl_assemble_sse_pool_join_sses.h"
#include "assemble/bcl_assemble_sse_pool_move_aa.h"
#include "assemble/bcl_assemble_sse_pool_mutate_sse.h"
#include "assemble/bcl_assemble_sse_pool_split_sse.h"
#include "biol/bcl_biol_membrane.h"
#include "find/bcl_find_locator.h"
#include "fold/bcl_fold_mutate_sse_type.h"
#include "math/bcl_math_binary_function_bind_second.h"
#include "math/bcl_math_binary_sum_function.h"
#include "math/bcl_math_mutate_decision_node.h"
#include "mc/bcl_mc_approximator.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_criterion_unimproved.h"
#include "score/bcl_score_sse_pool_sses.h"
#include "score/bcl_score_sse_predictions.h"
#include "sspred/bcl_sspred_sse_factory_highest.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEFactoryMC::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEFactoryMC( sspred::GetMethods().e_Undefined, 0.0))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from sspred method
    //! @param SSMETHOD sspred method to use to generate pool of sses
    //! @param CONFIDENCE_THRESHOLD the threshold in units of z-score, above which the confidence score becomes negative
    SSEFactoryMC::SSEFactoryMC
    (
      const sspred::Method &SSMETHOD,
      const double CONFIDENCE_THRESHOLD
    ) :
      m_Method( SSMETHOD),
      m_ConfidenceThreshold( CONFIDENCE_THRESHOLD),
      m_MaxNumberIterations( s_MaxNumberIterations),
      m_NumberOptimizations( s_NumberOptimizations)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSEFactoryMC
    SSEFactoryMC *SSEFactoryMC::Clone() const
    {
      return new SSEFactoryMC( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEFactoryMC::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the thresholds to use
    //! @brief SSTYPE_THRESHOLDS the thresholds to use for the sstypes desired
    void SSEFactoryMC::SetThresholds( const storage::Map< biol::SSType, double> &SSTYPE_THRESHOLDS)
    {
      if( SSTYPE_THRESHOLDS.IsEmpty())
      {
        return;
      }

      // take first value
      m_ConfidenceThreshold = SSTYPE_THRESHOLDS.Begin()->second;
    }

    //! @brief set the scoring function
    //! @param SP_SCORING ShPtr to scoring function
    void SSEFactoryMC::SetScoringFunction( const util::ShPtr< math::BinaryFunctionInterface< SSEPool, biol::Membrane, double> > &SP_SCORING)
    {
      m_ScoringFunction = SP_SCORING;
    }

    //! @brief set the mutate object
    //! @param SP_MUTATE ShPtr to mutate
    void SSEFactoryMC::SetMutate( const util::ShPtr< math::MutateInterface< SSEPool> > &SP_MUTATE)
    {
      m_Mutate = SP_MUTATE;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate the default scoring function
    //! @return function to score the sse pool as objective function
    util::ShPtr< math::BinaryFunctionInterface< SSEPool, biol::Membrane, double> > SSEFactoryMC::DefaultScoringFunction() const
    {
      // objective function
      util::ShPtr< score::SSEPredictions> sp_sse_preds( new score::SSEPredictions( m_Method, 0.0));
      util::ShPtr< math::BinarySumFunction< SSEPool, biol::Membrane, double, double> > sp_objective
      (
        new math::BinarySumFunction< SSEPool, biol::Membrane, double, double>()
      );
      sp_objective->NewOperand
      (
        util::ShPtr< score::SSEPoolSSEs>( new score::SSEPoolSSEs( sp_sse_preds, false, false)), 1.0
      );

//      sp_objective->NewOperand( util::ShPtr< score::SSEPoolSSEs>( new score::SSEPoolSSEs( sp_sse_pred_avg2, false, false)), 1.0);
//      sp_objective->NewOperand( util::ShPtr< score::SSEPoolSSEs>( new score::SSEPoolSSEs( sp_sse_pred_avg5, false, false)), 1.0);
//      sp_objective->NewOperand( util::ShPtr< score::SSEPoolSSEs>( new score::SSEPoolSSEs( sp_sse_pred_dev, false, false)), 1.0);
//      sp_objective->NewOperand( util::ShPtr< score::SSEPoolSSEs>( new score::SSEPoolSSEs( sp_sse_pred_slope, false, false)), 1.0);
//      util::ShPtr< score::SSEPredictions> sp_sse_preds( new score::SSEPredictions( m_Method, -0.0));
//      util::ShPtr< math::SumFunctionMixin< SSEPool, double> > sp_objective( new math::SumFunctionMixin< SSEPool, double>());
//      sp_objective->NewOperand( util::ShPtr< score::SSEPoolSSEs>( new score::SSEPoolSSEs( sp_sse_pred_avg2, false, false)), 0.25);
//      sp_objective->NewOperand( util::ShPtr< score::SSEPoolSSEs>( new score::SSEPoolSSEs( sp_sse_pred_avg3, false, false)), 0.25);
//      sp_objective->NewOperand( util::ShPtr< score::SSEPoolSSEs>( new score::SSEPoolSSEs( sp_sse_pred_avg4, false, false)), 0.25);
//      sp_objective->NewOperand( util::ShPtr< score::SSEPoolSSEs>( new score::SSEPoolSSEs( sp_sse_pred_avg5, false, false)), 0.25);
//      sp_objective->NewOperand( util::ShPtr< score::SSEPoolSSEs>( new score::SSEPoolSSEs( sp_sse_pred_dev0, false, false)), 0.33);
//      sp_objective->NewOperand( util::ShPtr< score::SSEPoolSSEs>( new score::SSEPoolSSEs( sp_sse_pred_dev1, false, false)), 0.33);
//      sp_objective->NewOperand( util::ShPtr< score::SSEPoolSSEs>( new score::SSEPoolSSEs( sp_sse_pred_dev2, false, false)), 0.33);
//      sp_objective->NewOperand( util::ShPtr< score::SSEPoolSSEs>( new score::SSEPoolSSEs( sp_sse_pred_slope2, false, false)), 0.5);
//      sp_objective->NewOperand( util::ShPtr< score::SSEPoolSSEs>( new score::SSEPoolSSEs( sp_sse_pred_slope3, false, false)), 0.5);
//      sp_objective->NewOperand( util::ShPtr< score::SSEPoolSSEs>( new score::SSEPoolSSEs( sp_sse_preds, false, false)), 1.0);

      return sp_objective;
    }

    //! @brief generate the default set of mutates
    //! @return mutate to change the compositon of the sse pool
    util::ShPtr< math::MutateInterface< SSEPool> > SSEFactoryMC::DefaultMutate() const
    {
      // collector for sses size < 8
      util::ShPtr< CollectorSSESize> sp_collector_sse_size_tiny
      (
        new CollectorSSESize( math::Range< size_t>( 1, 8))
      );

      // collector for sses size >= 4
      util::ShPtr< CollectorSSESize> sp_collector_sse_size_small
      (
        new CollectorSSESize( math::Range< size_t>( 4, util::GetUndefined< size_t>()))
      );

      util::ShPtr< CollectorSSESize> sp_collector_sse_size_large
      (
        new CollectorSSESize( math::Range< size_t>( 10, util::GetUndefined< size_t>()))
      );

      // pick sse random
      util::ShPtr< PickSSERandom> sp_pick_sse_random( new PickSSERandom());

      // locate random sse below size
      util::ShPtr< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> > sp_locate_sse_tiny_random
      (
        new find::Locator< util::SiPtr< const SSE>, DomainInterface, util::SiPtrList< const SSE> >
        (
          sp_collector_sse_size_tiny,
          sp_pick_sse_random
        )
      );

      // locate random sse above size
      util::ShPtr< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> > sp_locate_sse_small_random
      (
        new find::Locator< util::SiPtr< const SSE>, DomainInterface, util::SiPtrList< const SSE> >
        (
          sp_collector_sse_size_small,
          sp_pick_sse_random
        )
      );

      util::ShPtr< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> > sp_locate_sse_large_random
      (
        new find::Locator< util::SiPtr< const SSE>, DomainInterface, util::SiPtrList< const SSE> >
        (
          sp_collector_sse_size_large,
          sp_pick_sse_random
        )
      );

      // change sse type
      util::ShPtr< fold::MutateSSEType> sp_mutate_sse_type
      (
        new fold::MutateSSEType
        (
          storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND, biol::GetSSTypes().COIL)
        )
      );

      util::ShPtr< math::MutateDecisionNode< SSEPool> > sp_mutate_node( new math::MutateDecisionNode< SSEPool>());
      util::ShPtr< SSEPoolMoveAA> sp_mutate_move_aa( new SSEPoolMoveAA( math::Range< size_t>( 1, 3)));

      // split an sse and mutate the one with the lower probability
      util::ShPtr< SSEPoolSplitSSE> sp_mutate_split_change_sse
      (
        new SSEPoolSplitSSE
        (
          m_Method,
          sp_locate_sse_small_random,
          sp_mutate_sse_type
        )
      );

      // devide an sse without mutating either resulting sses
      util::ShPtr< SSEPoolSplitSSE> sp_mutate_devide_sse
      (
        new SSEPoolSplitSSE
        (
          m_Method,
          sp_locate_sse_small_random,
          util::ShPtr< math::MutateInterface< SSE> >()
        )
      );
      util::ShPtr< SSEPoolMutateSSE> sp_mutate_pool_sse_type
      (
        new SSEPoolMutateSSE
        (
          sp_locate_sse_large_random,
          sp_mutate_sse_type
        )
      );
      util::ShPtr< SSEPoolJoinSSEs> sp_mutate_pool_join
      (
        new SSEPoolJoinSSEs( sp_locate_sse_tiny_random)
      );

      util::ShPtr< SSEPoolInsertCoilIntoSSE> sp_mutate_pool_insert_coil_1
      (
        new SSEPoolInsertCoilIntoSSE
        (
          m_Method,
          sp_locate_sse_large_random,
          1
        )
      );
      util::ShPtr< SSEPoolInsertCoilIntoSSE> sp_mutate_pool_insert_coil_3
      (
        new SSEPoolInsertCoilIntoSSE
        (
          m_Method,
          sp_locate_sse_large_random,
          3
        )
      );

      // add all mutates with a probability of 1.0 to the node
      sp_mutate_node->AddMutate( *sp_mutate_move_aa, 1.0);
      sp_mutate_node->AddMutate( *sp_mutate_split_change_sse, 1.0);
      sp_mutate_node->AddMutate( *sp_mutate_devide_sse, 1.0);
      sp_mutate_node->AddMutate( *sp_mutate_pool_sse_type, 1.0);
      sp_mutate_node->AddMutate( *sp_mutate_pool_join, 1.0);
      sp_mutate_node->AddMutate( *sp_mutate_pool_insert_coil_1, 0.6);
      sp_mutate_node->AddMutate( *sp_mutate_pool_insert_coil_3, 0.4);

      // end
      return sp_mutate_node;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that returns a set of SSEs for the given AASequence
    //! @brief SEQUENCE AASequence from which the SSEPool is going to be built
    //! @return SSEPool built from provided SEQUENCE
    SSEPool
    SSEFactoryMC::operator()( const biol::AASequence &SEQUENCE) const
    {
      // check that sequence has one aa
      if( SEQUENCE.GetSize() < 1)
      {
        return SSEPool();
      }

      // check that sequence has defined prediction method
      if( !SEQUENCE.GetFirstAA()->GetSSPrediction( m_Method).IsDefined())
      {
        BCL_MessageCrt
        (
          "cannot derive pool, since sequence does not have defined sspred method: " + m_Method.GetName()
        );
        return SSEPool();
      }

      // start pool
      SSEPool start_pool( sspred::SSEFactoryHighest( m_Method)( SEQUENCE));

      // objective function
      util::ShPtr< math::BinaryFunctionInterface< SSEPool, biol::Membrane, double> > sp_binary_objective;
      if( m_ScoringFunction.IsDefined())
      {
        sp_binary_objective = m_ScoringFunction;
      }
      else
      {
        sp_binary_objective = DefaultScoringFunction();
      }

      // bind the second argument (membrane) so that the function only takes one argument,
      // this assumes the scoring function does not care about the membrane
      util::ShPtr< math::FunctionInterfaceSerializable< SSEPool, double> > sp_objective
      (
        new math::BinaryFunctionBindSecond< SSEPool, biol::Membrane, double>
        (
          sp_binary_objective,
          biol::Membrane()
        )
      );

      // mutate
      util::ShPtr< math::MutateInterface< SSEPool> > sp_mutate_node;
      if( m_Mutate.IsDefined())
      {
        sp_mutate_node = m_Mutate;
      }
      else
      {
        sp_mutate_node = DefaultMutate();
      }

      util::ShPtr< storage::Pair< SSEPool, double> > sp_final_result
      (
        new storage::Pair< SSEPool, double>
        (
          SSEPool(), std::numeric_limits< double>::max()
        )
      );

      for( size_t pool_nr( 0); pool_nr < m_NumberOptimizations; ++pool_nr)
      {
        // metropolis
        util::ShPtr< mc::Metropolis< double> > sp_metropolis
        (
          new mc::Metropolis< double>
          (
            util::ShPtr< mc::TemperatureInterface>
            (
              new mc::TemperatureAccepted( m_MaxNumberIterations)
            )
          )
        );

        // create termination criterion
        opti::CriterionCombine< SSEPool, double> termination_criterion;
        termination_criterion.InsertCriteria( opti::CriterionUnimproved< SSEPool, double>( s_MaxNumberUnimproved));
        termination_criterion.InsertCriteria( opti::CriterionNumberIterations< SSEPool, double>( m_MaxNumberIterations));

        // create the approximator
        mc::Approximator< SSEPool, double> approximator
        (
          *sp_objective,
          *sp_mutate_node,
          *sp_metropolis,
          termination_criterion,
          start_pool
        );

        // approximate
        approximator.Approximate();
        const util::ShPtr< storage::Pair< SSEPool, double> > sp_result( approximator.GetTracker().GetBest());
        if( sp_final_result->Second() > sp_result->Second())
        {
          sp_final_result = sp_result;
        }
      }

      return sp_final_result->First();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEFactoryMC::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Method             , ISTREAM);
      io::Serialize::Read( m_ConfidenceThreshold, ISTREAM);
      io::Serialize::Read( m_ScoringFunction    , ISTREAM);
      io::Serialize::Read( m_Mutate             , ISTREAM);
      io::Serialize::Read( m_MaxNumberIterations, ISTREAM);
      io::Serialize::Read( m_NumberOptimizations, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEFactoryMC::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Method             , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ConfidenceThreshold, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ScoringFunction    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Mutate             , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MaxNumberIterations, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NumberOptimizations, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
