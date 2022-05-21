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

#ifndef BCL_OPTI_TRACKER_WITH_HISTORY_H_
#define BCL_OPTI_TRACKER_WITH_HISTORY_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_improvement_type.h"
#include "bcl_opti_phase.h"
#include "bcl_opti_tracker.h"
#include "bcl_opti_tracker_base.h"
#include "math/bcl_math_running_average.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_message.h"
#include "util/bcl_util_sh_ptr_list.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TrackerWithHistory
    //! @brief Tracker for optimization algorithms
    //! @details TrackerWithHistory is used by Approximator classes to keep of the approximation process when a full or
    //! @details partial history is needed. This implementation of it keeps track of everything by default, as well as
    //! @details the initial argument
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument
    //! @tparam t_ResultType is the type of the approximation result
    //!
    //! @see @link example_opti_tracker_with_history.cpp @endlink
    //! @author geanesar
    //! @date Sept 26 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class TrackerWithHistory :
      public Tracker< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    protected:

      //! The number of previous steps that should be kept
      size_t m_NumberStepsToKeep;

      //! The initial member
      util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > m_Initial;

      //! History of previous states
      util::ShPtrList< storage::Pair< t_ArgumentType, t_ResultType> > m_History;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      TrackerWithHistory( const size_t STEPS_TO_KEEP = 1) :
        Tracker< t_ArgumentType, t_ResultType>(),
        m_NumberStepsToKeep( STEPS_TO_KEEP),
        m_Initial(),
        m_History()
      {
      }

      //! @brief construct from ImprovementType
      //! @param STEPS_TO_KEEP the number of steps to keep track of (most recent steps)
      //! @param IMPROVEMENT_TYPE change type that indicates an improvement
      //! @param AVERAGING_WINDOW_SIZE size of averaging window
      TrackerWithHistory( const ImprovementType &IMPROVEMENT_TYPE, const size_t &AVERAGING_WINDOW_SIZE = 1, const size_t &STEPS_TO_KEEP = 1) :
        Tracker< t_ArgumentType, t_ResultType>( IMPROVEMENT_TYPE, AVERAGING_WINDOW_SIZE),
        m_NumberStepsToKeep( STEPS_TO_KEEP),
        m_Initial(),
        m_History()
      {
      }

      //! @brief Clone function
      //! @return pointer to a new Tracker< t_ArgumentType, t_ResultType>
      TrackerWithHistory *Clone() const
      {
        return new TrackerWithHistory< t_ArgumentType, t_ResultType>( *this);
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

      //! @brief returns a ShPtrList of all the recorded history
      //! @return ShPtrList to all past models
      const util::ShPtrList< storage::Pair< t_ArgumentType, t_ResultType> > &GetHistory() const
      {
        return m_History;
      }

      //! @brief Return the initial argument
      //! @return ShPtr to the initial argument and score
      const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &GetInitial() const
      {
        return m_Initial;
      }

      //! @brief resets the tracker
      void Reset()
      {
        // reset members
        Tracker< t_ArgumentType, t_ResultType>::Reset();
        m_Initial.Reset();
        m_History.Reset();
      }

      //! @brief tracks the given model.  stores more recently added models at the front (queue format)
      //! @param SP_MODEL ShPtr to the model to be tracked
      //! @param STEP_STATUS; if the tracker should determine this, leave as s_NumberStepStatus
      void Track
      (
        const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &SP_MODEL,
        StepStatus STEP_STATUS = s_NumberStepStatus
      )
      {
        // Call the tracker
        Tracker< t_ArgumentType, t_ResultType>::Track( SP_MODEL, STEP_STATUS);

        // Set the initial member
        if( Tracker< t_ArgumentType, t_ResultType>::GetPhase() == e_Start)
        {
          m_Initial = SP_MODEL;
          m_History.Reset();
          m_History.PushFront( SP_MODEL);
        }
        else
        {
          // Add the model to the history
          m_History.PushFront( SP_MODEL);

          // Remove previous iterations if we don't want to keep all of them (StepsToKeep=0)
          size_t history_size( m_History.GetSize());
          while( m_NumberStepsToKeep && history_size > m_NumberStepsToKeep)
          {
            m_History.PopBack();
            --history_size;
          }
        }
      }

      //! @brief sets the number of steps that should be tracked
      //! @param NUMBER the number to keep in the history
      void SetMaxHistorySize( const size_t &NUMBER)
      {
        m_NumberStepsToKeep = NUMBER;
        // Remove previous iterations if we have resized to a smaller size
        size_t history_size( m_History.GetSize());
        while( m_NumberStepsToKeep && history_size > m_NumberStepsToKeep)
        {
          m_History.PopBack();
          --history_size;
        }
      }

      //! @brief returns the number of steps that will be kept in the history
      size_t GetMaxHistorySize() const
      {
        return m_NumberStepsToKeep;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        Tracker< t_ArgumentType, t_ResultType>::Read( ISTREAM);
        io::Serialize::Read( m_History, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        Tracker< t_ArgumentType, t_ResultType>::Write( OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_History, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }
    }; // class Tracker< t_ArgumentType, t_ResultType>
  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_TRACKER_WITH_HISTORY_H_
