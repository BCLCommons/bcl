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

#ifndef BCL_MC_APPROXIMATOR_H_
#define BCL_MC_APPROXIMATOR_H_

// include the namespace header
#include "bcl_mc.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_mc_metropolis.h"
#include "assemble/bcl_assemble_printer_protein_model.h"
#include "math/bcl_math_mutate_result.h"
#include "opti/bcl_opti_approximator_modular_base.h"
#include "pdb/bcl_pdb_factory.h"
#include "util/bcl_util_logger_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace mc
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Approximator
    //! @brief This class is used for Monte Carlo minimization in the protein folding process
    //! @details Approximator knows about the mutates used, the metropolis criteria used, the tracker for the
    //! minimization as well as the termination criteria used.
    //!
    //! @tparam t_ArgumentType the type of the argument
    //! @tparam t_ResultType the type of the result
    //!
    //! @see @link example_mc_approximator.cpp @endlink
    //! @author karakam, fischea
    //! @date Dec 19, 2012
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class Approximator :
      public opti::ApproximatorModularBase< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    private:

      //! ShPtr to the objective function
      util::ShPtr< math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> > m_Objective;

      //! ShPtr to the mutate object that mutates the object
      util::ShPtr< math::MutateInterface< t_ArgumentType> > m_Mutate;

      //! ShPtr to Metropolis criteria used
      util::ShPtr< Metropolis< t_ResultType> > m_Metropolis;

      //! Upon rejected step, revert to the best model seen so far with the given probability. This is the default behavior for backwards compatibility
      //! and is primarily useful in problems where a rejected step suggests a dead-end that is difficult to recover from
      double m_RejectedStepRevertsToBest;

      //! Whether we have exceeded a prespecified maximum number of skipped steps in a row
      bool m_AllStepsSkipped;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Approximator() :
        opti::ApproximatorModularBase< t_ArgumentType, t_ResultType>( opti::e_SmallerIsBetter),
        m_Objective(),
        m_Mutate(),
        m_Metropolis(),
        m_RejectedStepRevertsToBest( 1.0),
        m_AllStepsSkipped( false)
      {
      }

      //! @brief construct from members
      //! @param OBJECTIVE the objective used to evaluate the models
      //! @param MUTATE the mutates allowed in the algorithm
      //! @param METROPOLIS the metropolis to evaluate the mutate outcome
      //! @param CRITERION the termination criterion
      Approximator
      (
        const math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> &OBJECTIVE,
        const math::MutateInterface< t_ArgumentType> &MUTATE,
        const Metropolis< t_ResultType> &METROPOLIS,
        const opti::CriterionInterface< t_ArgumentType, t_ResultType> &CRITERION,
        const double &REJECTED_STEP_REVERTS_TO_BEST_PROB = 1.0
      ) :
        opti::ApproximatorModularBase< t_ArgumentType, t_ResultType>( opti::e_SmallerIsBetter),
        m_Objective( util::CloneToShPtr( OBJECTIVE)),
        m_Mutate( util::CloneToShPtr( MUTATE)),
        m_Metropolis( util::CloneToShPtr( METROPOLIS)),
        m_RejectedStepRevertsToBest( REJECTED_STEP_REVERTS_TO_BEST_PROB),
        m_AllStepsSkipped( false)
      {
        // set the termination criterion
        this->SetCriterion( CRITERION);
      }

      //! @brief construct from members
      //! @param OBJECTIVE the objective used to evaluate the models
      //! @param MUTATE the mutates allowed in the algorithm
      //! @param METROPOLIS the metropolis to evaluate the mutate outcome
      //! @param CRITERION the termination criterion
      //! @param INITIAL_ARGUMENT starting model for the mc algorithm
      //! @param TRACKER the tracker to keep track of the best molecule or conformer
      Approximator
      (
        const math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> &OBJECTIVE,
        const math::MutateInterface< t_ArgumentType> &MUTATE,
        const Metropolis< t_ResultType> &METROPOLIS,
        const opti::CriterionInterface< t_ArgumentType, t_ResultType> &CRITERION,
        const t_ArgumentType &INITIAL_ARGUMENT,
        const opti::Tracker< t_ArgumentType, t_ResultType> &TRACKER = opti::Tracker< t_ArgumentType, t_ResultType>( opti::e_SmallerIsBetter),
        const double &REJECTED_STEP_REVERTS_TO_BEST_PROB = 1.0
      ) :
        opti::ApproximatorModularBase< t_ArgumentType, t_ResultType>( TRACKER),
        m_Objective( util::CloneToShPtr( OBJECTIVE)),
        m_Mutate( util::CloneToShPtr( MUTATE)),
        m_Metropolis( util::CloneToShPtr( METROPOLIS)),
        m_RejectedStepRevertsToBest( REJECTED_STEP_REVERTS_TO_BEST_PROB),
        m_AllStepsSkipped( false)
      {
        // set the termination criterion
        this->SetCriterion( CRITERION);

        // track the starting model
        this->Track
        (
          util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> >
          (
            new storage::Pair< t_ArgumentType, t_ResultType>
            (
              INITIAL_ARGUMENT, m_Objective->operator()( INITIAL_ARGUMENT)
            )
          )
        );
      }

      //! @brief construct from members
      //! @param OBJECTIVE the objective used to evaluate the models
      //! @param MUTATE the mutates allowed in the algorithm
      //! @param METROPOLIS the metropolis to evaluate the mutate outcome
      //! @param CRITERION the termination criterion
      //! @param INITIAL_ARGUMENT starting model for the mc algorithm
      //! @param TRACKER the tracker to keep track of the best molecule or conformer
      Approximator
      (
        const util::ShPtr< math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> > &OBJECTIVE,
        const util::ShPtr< math::MutateInterface< t_ArgumentType> > &MUTATE,
        const Metropolis< t_ResultType> &METROPOLIS,
        const opti::CriterionInterface< t_ArgumentType, t_ResultType> &CRITERION,
        const t_ArgumentType &INITIAL_ARGUMENT,
        const opti::Tracker< t_ArgumentType, t_ResultType> &TRACKER = opti::Tracker< t_ArgumentType, t_ResultType>( opti::e_SmallerIsBetter),
        const bool &REJECTED_STEP_REVERTS_TO_BEST_PROB = 1.0
      ) :
        opti::ApproximatorModularBase< t_ArgumentType, t_ResultType>( TRACKER),
        m_Objective( OBJECTIVE),
        m_Mutate( MUTATE),
        m_Metropolis( util::CloneToShPtr( METROPOLIS)),
        m_RejectedStepRevertsToBest( REJECTED_STEP_REVERTS_TO_BEST_PROB),
        m_AllStepsSkipped( false)
      {
        // set the termination criterion
        this->SetCriterion( CRITERION);

        // track the starting model
        this->Track
        (
          util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> >
          (
            new storage::Pair< t_ArgumentType, t_ResultType>
            (
              INITIAL_ARGUMENT, m_Objective->operator()( INITIAL_ARGUMENT)
            )
          )
        );
      }

      //! @brief Clone function
      //! @return pointer to a new Approximator< t_ArgumentType, t_ResultType>
      Approximator< t_ArgumentType, t_ResultType> *Clone() const
      {
        return new Approximator< t_ArgumentType, t_ResultType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief returns the objective pointer
      //! @return the objective pointer
      const util::ShPtr< math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> > &GetObjective() const
      {
        return m_Objective;
      }

      //! @brief returns the mutator
      const util::ShPtr< math::MutateInterface< t_ArgumentType> > &GetMutator() const
      {
        return m_Mutate;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns the approximation result for the given argument
      //! @param INITIAL_ARGUMENT argument to return the approximation result for
      //! @return approximation result for the given argument
      util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > ApproximateReturnCopy
      (
        const t_ArgumentType &INITIAL_ARGUMENT
      )
      {
        // track the initial argument result pair and approximate
        util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > sp_pair
        (
          new storage::Pair< t_ArgumentType, t_ResultType>
          (
            INITIAL_ARGUMENT, m_Objective->operator()( INITIAL_ARGUMENT)
          )
        );
        this->GetTracker().Track( sp_pair);
        this->Approximate();

        return this->GetTracker().GetBest();
      }

      //! @brief conducts the next approximation step and stores the approximation
      void Next()
      {
        // compute the next MC step and store the mutate result
        math::MutateResult< t_ArgumentType> mutate_result
        (
          m_Mutate->operator()( this->GetTracker().GetCurrent()->First())
        );

        // define an upper limit to the number of skipped steps
        static const size_t s_max_skipped_steps_in_a_row( 1000);

        // get a defined argument
        for( size_t n_skipped( 0); !mutate_result.GetArgument().IsDefined() && n_skipped < s_max_skipped_steps_in_a_row; ++n_skipped)
        {
          // increment skipped steps
          this->GetTracker().IncrementSkippedSteps();

          // form the mutate string
          const std::string mutate_string
          (
            ( mutate_result.GetNodes().IsEmpty() ? "_no_mutate" : mutate_result.GetNodes().LastElement()->GetScheme())
          );
          // message what status the current step has
          BCL_MessageVrb
          (
            "Skipped step: " + mutate_string
          );

          // try a different step
          mutate_result = m_Mutate->operator()( this->GetTracker().GetCurrent()->First());
        }
        if( !mutate_result.GetArgument().IsDefined())
        {
          m_AllStepsSkipped = true;
          return;
        }

        // if the mutate was successful (thus returned a defined model)
        const t_ArgumentType &arg( *mutate_result.GetArgument());

        // create a storage pair for the mutated model
        util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > sp_mutated_model
        (
          new storage::Pair< t_ArgumentType, t_ResultType>
          (
            arg,
            m_Objective->operator()( arg)
          )
        );

        // let Metropolis determine the step status
        opti::StepStatus step_status
        (
          m_Metropolis->Evaluate( this->GetTracker().GetBest()->Second(), sp_mutated_model->Second(), this->GetTracker())
        );

        // form the mutate string
        const std::string mutate_string
        (
          ( mutate_result.GetNodes().IsEmpty() ? "_no_mutate" : mutate_result.GetNodes().LastElement()->GetScheme())
        );

        // message what status the current step has
        BCL_MessageVrb
        (
          "Metropolis outcome for step " + util::Format()( this->GetTracker().GetIteration()) +
          " is " + opti::GetStepStatusName( step_status) +
          " for a score of " + util::Format()( sp_mutated_model->Second()) +
          " " + mutate_string
        );

        // let the tracker know about the new model
        this->Track( sp_mutated_model, step_status);

        // if the mutation improved the model or metropolis accepted it anyway, track the mutated model
        if( step_status == opti::e_Rejected)
        {
          if( m_RejectedStepRevertsToBest == 1.0 || random::GetGlobalRandom().Double() < m_RejectedStepRevertsToBest)
          {
            this->GetTracker().SetCurrent( this->GetTracker().GetBest());
          }
        }

        // track the step scheme
        this->GetTracker().SetStepScheme( mutate_string);

      } // function Next

      //! @brief evaluates whether the approximation can continue
      //! @detail checks if the number of steps without improvement exceed the given maximum number of steps
      //! without improvement
      //! @return true, if the approximation can continue - otherwise false
      bool CanContinue() const
      {
        return !m_AllStepsSkipped;
      }

    protected:

    //////////////////////
    // Input and output //
    //////////////////////

      //! @brief prints the start information
      void PrintStartHook()
      {
        BCL_MessageStd( "MC Minimization started");
      }

      //! @brief prints the step information
      void PrintStepHook()
      {
        util::GetLogger().LogStatus
        (
          "Iteration " + util::Format()( this->GetTracker().GetIteration())
          + " " + util::Format()( this->GetTracker().GetCurrent()->Second())
          + " " + util::Format()( this->GetTracker().GetBest()->Second())
        );
      }

      //! @brief prints the end information
      void PrintEndHook()
      {
        BCL_MessageStd( "MC Minimization ended");

        // inform user about the steps taken
        const storage::Vector< size_t> &counts( this->GetTracker().GetCounts());
        const size_t &tot_nr_steps( this->GetTracker().GetIteration());
        const size_t &nr_improved( counts( opti::e_Improved));
        const size_t &nr_accepted( counts( opti::e_Accepted));
        const size_t &nr_rejected( counts( opti::e_Rejected));
        const size_t &nr_skipped( counts( opti::e_Skipped));
        util::Format format_a, format_b;
        format_a.W( 5);
        format_b.W( 5).FFP( 2);

        BCL_MessageStd( "#MC steps: " + util::Format()( tot_nr_steps));
        BCL_MessageStd
        (
          "#MC steps improved:\t" + format_a( nr_improved) + "\t%" + format_b( 100.0 * nr_improved / tot_nr_steps)
        );
        BCL_MessageStd
        (
          "#MC steps accepted:\t" + format_a( nr_accepted) + "\t%" + format_b( 100.0 * nr_accepted / tot_nr_steps)
        );
        BCL_MessageStd
        (
          "#MC steps rejected:\t" + format_a( nr_rejected) + "\t%" + format_b( 100.0 * nr_rejected / tot_nr_steps)
        );
        BCL_MessageStd
        (
          "#MC steps skipped:\t" + format_a( nr_skipped) + "\t%" + format_b( 100.0 * nr_skipped / tot_nr_steps)
        );
      }

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // call Read of base class
        opti::ApproximatorModularBase< t_ArgumentType, t_ResultType>::Read( ISTREAM);

        // read members
        io::Serialize::Read( m_Objective, ISTREAM);
        io::Serialize::Read( m_Mutate, ISTREAM);
        io::Serialize::Read( m_Metropolis, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // call Write of base class
        opti::ApproximatorModularBase< t_ArgumentType, t_ResultType>::Write( OSTREAM, INDENT) << '\n';

        // write members
        io::Serialize::Write( m_Objective, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Mutate, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Metropolis, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class Approximator< t_ArgumentType, t_ResultType>

  } // namespace mc
} // namespace bcl

#endif // BCL_MC_APPROXIMATOR_H_
