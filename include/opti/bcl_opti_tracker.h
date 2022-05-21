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

#ifndef BCL_OPTI_TRACKER_H_
#define BCL_OPTI_TRACKER_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_improvement_type.h"
#include "bcl_opti_phase.h"
#include "bcl_opti_tracker_base.h"
#include "math/bcl_math_running_average.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_string_functions.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Tracker
    //! @brief Tracker for optimization algorithms
    //! @details Tracker is used by ApproximatorModularInterface derived classes to keep of the approximation process.
    //! Therefore Tracker stores the current model and the best model. It distinguishes between the improvement types
    //! "smaller is better", "larger is better" and "smaller absolute value is better" to evaluate the passed models.
    //! Tracker also keeps track of the number of iterations for the current model and the best model as well as the
    //! number of iterations since the last improvement.
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument
    //! @tparam t_ResultType is the type of the approximation result
    //!
    //! @see @link example_opti_tracker.cpp @endlink
    //! @author fischea
    //! @date Dec 12, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class Tracker :
      public TrackerBase
    {

    //////////
    // data //
    //////////

    protected:

      //! ShPtr to the current model
      util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > m_CurrentModel;

      //! ShPtr to the best model
      util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > m_BestModel;

      //! ShPtr to the best model
      util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > m_LastModel;

      //! Running average of best results
      math::RunningAverage< t_ResultType> m_AverageResult;

      //! Linearly weighted average result
      //! For example: e.g. past result m_ResultAverageWindow has weight 0, m_ResultAverageWindow+1 has weight 1, etc.
      math::RunningAverage< t_ResultType> m_WeightedAverageResult;

      //! Best average running result
      t_ResultType m_BestAverageResult;

      //! Window over which to compute the average result
      size_t m_ResultAverageWindow;

      //! list of past m_ResultAverageWindow results, if not <= 1
      storage::List< t_ResultType> m_PastResults;

      //! scheme of the last tracked step
      std::string m_StepScheme;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Tracker() :
        TrackerBase( e_SmallerIsBetter),
        m_CurrentModel(),
        m_BestModel(),
        m_ResultAverageWindow( 1),
        m_StepScheme( "")
      {
      }

      //! @brief construct from ImprovementType
      //! @param IMPROVEMENT_TYPE change type that indicates an improvement
      //! @param AVERAGING_WINDOW_SIZE size of averaging window
      Tracker( const ImprovementType &IMPROVEMENT_TYPE, const size_t &AVERAGING_WINDOW_SIZE = 1) :
        TrackerBase( IMPROVEMENT_TYPE),
        m_CurrentModel(),
        m_BestModel(),
        m_ResultAverageWindow( AVERAGING_WINDOW_SIZE),
        m_StepScheme( "")
      {
      }

      //! @brief Clone function
      //! @return pointer to a new Tracker< t_ArgumentType, t_ResultType>
      Tracker *Clone() const
      {
        return new Tracker< t_ArgumentType, t_ResultType>( *this);
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

      //! @brief returns a ShPtr to the current model
      //! @return ShPtr to the current model
      const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &GetCurrent() const
      {
        return m_CurrentModel;
      }

      //! @brief returns a ShPtr to the current model
      //! @return ShPtr to the current model
      void SetCurrent( const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &CURRENT)
      {
        m_CurrentModel = CURRENT;
      }

      //! @brief returns a ShPtr to the current model
      //! @return ShPtr to the current model
      const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &GetLast() const
      {
        return m_LastModel;
      }

      //! @brief returns a ShPtr to the best model
      //! @return ShPtr to the best model
      const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &GetBest() const
      {
        return m_BestModel;
      }

      //! @brief get a tag for the output from by this tracker; used in a variety of printers
      //! @return a tag
      std::string GetTag() const
      {
        return
          m_Phase == e_Start ? "start"
          : m_Phase == e_Iteration
            ? "_iteration" + util::Format().Fill( '0').W( 6)( m_IterationNumber)
            : "final"; // end
      }

      //! @brief Remove any tag that the tracker may have added to a storage key
      //! @param KEY a string that may currently start with a tag given by the tracker
      //! @return KEY without any tracker-tag prefix
      static std::string RemoveTrackerTag( const std::string &KEY)
      {
        if( util::StartsWith( KEY, "start") || util::StartsWith( KEY, "final"))
        {
          return RemoveTrackerTag( KEY.substr( 5));
        }
        else if( util::StartsWith( KEY, "_iteration") && KEY.size() > size_t( 16))
        {
          return RemoveTrackerTag( KEY.substr( 16));
        }
        return KEY;
      }

      //! @brief returns the size of the averaging window
      //! @return size of averaging window
      const size_t &GetAveragingWindowSize() const
      {
        return m_ResultAverageWindow;
      }

      //! @brief returns the size of the averaging window
      //! @param AVERAGING_WINDOW size of averaging window
      //! @return size of averaging window
      void SetAveragingWindowSize( const size_t &AVERAGING_WINDOW)
      {
        m_ResultAverageWindow = AVERAGING_WINDOW;
      }

      //! @brief returns the scheme of the last tracked step
      //! @return scheme of the last tracked step
      const std::string &GetStepScheme() const
      {
        return m_StepScheme;
      }

      //! @brief set the scheme for the last tracked step
      //! @param SCHEME scheme for the last tracked step
      void SetStepScheme( const std::string &STEP_SCHEME)
      {
        m_StepScheme = STEP_SCHEME;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief resets the tracker_WeightedAverageResult
      virtual void Reset()
      {
        // reset members
        TrackerBase::Reset();
        m_CurrentModel.Reset();
        m_BestModel.Reset();
        m_AverageResult.Reset();
        m_PastResults.Reset();
      }

      //! @brief tracks the given model
      //! @param SP_MODEL ShPtr to the model to be tracked
      //! @param STEP_STATUS; if the tracker should determine this, leave as s_NumberStepStatus
      virtual void Track
      (
        const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &SP_MODEL,
        StepStatus STEP_STATUS = s_NumberStepStatus
      )
      {
        // update result average
        if( m_ResultAverageWindow > 1)
        {
          if( this->m_IterationNumber >= m_ResultAverageWindow)
          {
            m_WeightedAverageResult -= m_AverageResult.GetAverage();
            m_AverageResult -= m_PastResults.FirstElement();
            m_PastResults.PopFront();
          }
          m_AverageResult += SP_MODEL->Second();
          m_PastResults.PushBack( SP_MODEL->Second());
          m_WeightedAverageResult.AddWeightedObservation
          (
            SP_MODEL->Second(),
            std::min( float( this->m_IterationNumber + 1), float( m_ResultAverageWindow))
          );
        }

        // determine step status if it was unknown
        if( STEP_STATUS == s_NumberStepStatus)
        {
          STEP_STATUS = e_Rejected;
          // no previous best model or no window -> last result is automatically the best
          if( !m_BestModel.IsDefined() || !m_ResultAverageWindow)
          {
            STEP_STATUS = e_Improved;
          }
          else if( DoesImprove( SP_MODEL->Second(), m_BestModel->Second(), m_ImprovementType))
          {
            if( m_ResultAverageWindow == size_t( 1))
            {
              STEP_STATUS = e_Improved;
            }
            else if( DoesImprove( m_WeightedAverageResult.GetAverage(), m_BestAverageResult, m_ImprovementType))
            {
              STEP_STATUS = e_Improved;
              m_BestAverageResult = m_WeightedAverageResult.GetAverage();
            }
            else
            {
              // no improvement in running average
              STEP_STATUS = e_Accepted;
            }
          }
          else
          {
            // an accepted step, since the calling code did not indicate otherwise
            STEP_STATUS = e_Accepted;
          }
        }

        this->Update( STEP_STATUS);

        m_LastModel = SP_MODEL;

        if( STEP_STATUS == e_Improved || STEP_STATUS == e_Accepted)
        {
          m_CurrentModel = SP_MODEL;
        }

        // if there is no current best model, or this model improves on the best
        if( STEP_STATUS == e_Improved)
        {
          m_BestModel = SP_MODEL;
        }
      }

      //! @brief note that a step was skipped; this does not increase the iteration number
      void IncrementSkippedSteps()
      {
        this->Update( e_Skipped);
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        TrackerBase::Read( ISTREAM);
        io::Serialize::Read( m_CurrentModel, ISTREAM);
        io::Serialize::Read( m_BestModel, ISTREAM);
        io::Serialize::Read( m_StepScheme, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        TrackerBase::Write( OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_CurrentModel, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_BestModel, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_StepScheme, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class Tracker< t_ArgumentType, t_ResultType>

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_TRACKER_H_
