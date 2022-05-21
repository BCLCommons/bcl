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
#include "opti/bcl_opti_tracker_base.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    TrackerBase::TrackerBase() :
      m_IterationNumber( 0),
      m_BestModelIterationNumber( 0),
      m_ImprovementType( e_SmallerIsBetter),
      m_Phase( e_Start),
      m_Counts( s_NumberStepStatus, size_t( 0)),
      m_IntervalStart( s_NumberStepStatus, size_t( 0)),
      m_IntervalEnd( s_NumberStepStatus, size_t( 0)),
      m_StatusOfLastStep()
    {
    }

    //! @brief construct from ImprovementType
    //! @param IMPROVEMENT_TYPE change type that indicates an improvement
    TrackerBase::TrackerBase( const ImprovementType &IMPROVEMENT_TYPE) :
      m_IterationNumber( 0),
      m_BestModelIterationNumber( 0),
      m_ImprovementType( IMPROVEMENT_TYPE),
      m_Phase( e_Start),
      m_Counts( s_NumberStepStatus, size_t( 0)),
      m_IntervalStart( s_NumberStepStatus, size_t( 0)),
      m_IntervalEnd( s_NumberStepStatus, size_t( 0)),
      m_StatusOfLastStep()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the number of iterations since the last improvement
    //! @return number of iterations since the last improvement
    size_t TrackerBase::GetIterationsSinceLastImprovement() const
    {
      return m_IterationNumber - m_BestModelIterationNumber;
    }

    //! @brief resets the tracker
    void TrackerBase::Reset()
    {
      // reset members
      m_IterationNumber = 0;
      m_BestModelIterationNumber = 0;
      m_Phase = e_Start;
      m_Counts.SetAllElements( 0);
      m_IntervalStart.SetAllElements( 0);
      m_IntervalEnd.SetAllElements( 0);
      m_StatusOfLastStep = e_Rejected;
    }

    //! @brief sets the phase of the iteration
    //! @PHASE the now current phase of the iteration
    void TrackerBase::SetPhase( const Phase &PHASE)
    {
      m_Phase = PHASE;
    }

    //! @brief return number of steps in a row of specified STEP_STATUS
    //! @param STEP_STATUS StepStatus of interest
    //! @return number of steps in a row of specified STEP_STATUS
    size_t TrackerBase::GetNumberStepsInARow( const StepStatus STEP_STATUS) const
    {
      // if this step status was observed in the last step
      if( m_IntervalEnd( STEP_STATUS) == m_IterationNumber)
      {
        return m_IterationNumber - m_IntervalStart( STEP_STATUS);
      }

      // otherwise it means some other step status type was observed last, therefore return 0
      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &TrackerBase::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_IterationNumber, ISTREAM);
      io::Serialize::Read( m_BestModelIterationNumber, ISTREAM);
      io::Serialize::Read( m_ImprovementType, ISTREAM);
      io::Serialize::Read( m_Phase, ISTREAM);
      io::Serialize::Read( m_Counts, ISTREAM);
      io::Serialize::Read( m_IntervalStart, ISTREAM);
      io::Serialize::Read( m_IntervalEnd, ISTREAM);
      io::Serialize::Read( m_StatusOfLastStep, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &TrackerBase::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write memebers
      io::Serialize::Write( m_IterationNumber, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BestModelIterationNumber, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ImprovementType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Phase, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Counts, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IntervalStart, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_IntervalEnd, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_StatusOfLastStep, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief Update the status of this tracker
    void TrackerBase::Update( const StepStatus STEP_STATUS)
    {
      // increment counts for the provided StepStatus
      ++m_Counts( STEP_STATUS);

      if( STEP_STATUS == e_Skipped)
      {
        // otherwise ignore skipped steps
        return;
      }

      // if previous iteration was not of this step status
      // in which case the end of interval would be not equal to last iteration number
      if( m_IntervalEnd( STEP_STATUS) < m_IterationNumber - 1)
      {
        // update the begin of the interval to this iteration
        m_IntervalStart( STEP_STATUS) = m_IterationNumber;
      }

      // update the end of interval to current iteration
      m_IntervalEnd( STEP_STATUS) = m_IterationNumber;

      // update last step
      m_StatusOfLastStep = STEP_STATUS;

      ++m_IterationNumber;

      if( STEP_STATUS == e_Improved)
      {
        m_BestModelIterationNumber = m_IterationNumber;
      }
    }

  } // namespace opti

} // namespace bcl
