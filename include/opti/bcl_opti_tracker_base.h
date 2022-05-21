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

#ifndef BCL_OPTI_TRACKER_BASE_H_
#define BCL_OPTI_TRACKER_BASE_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_improvement_type.h"
#include "bcl_opti_phase.h"
#include "bcl_opti_step_status.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TrackerBase
    //! @brief Base class for all trackers (template-arg independent)
    //!
    //! @see @link example_opti_tracker_base.cpp @endlink
    //! @author mendenjl
    //! @date Sep 03, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API TrackerBase :
      public util::ObjectInterface
    {

    protected:

      //! number of the current iteration
      size_t m_IterationNumber;

      //! number of iterations for the best model
      size_t m_BestModelIterationNumber;

      //! type of improvement
      ImprovementTypeEnum m_ImprovementType;

      //! Current phase
      PhaseEnum m_Phase;

      //! Stores for each status, how many times it was observed since the beginning of minimization
      storage::Vector< size_t> m_Counts;

      //! Stores for each step status the begin and end iteration numbers for the last interval it was observed
      storage::Vector< size_t> m_IntervalStart;
      storage::Vector< size_t> m_IntervalEnd;

      //! Last step
      StepStatusEnum m_StatusOfLastStep;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      TrackerBase();

      //! @brief construct from ImprovementType
      //! @param IMPROVEMENT_TYPE change type that indicates an improvement
      TrackerBase( const ImprovementType &IMPROVEMENT_TYPE);

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the current number of iterations
      //! @return current number of iterations
      size_t GetIteration() const
      {
        return m_IterationNumber;
      }

      //! @brief returns the number of iterations for the best model
      //! @return number of iterations for the best model
      size_t GetBestIteration() const
      {
        return m_BestModelIterationNumber;
      }

      //! @brief get the current phase of the approximation
      //! @return current phase of the approximation
      const PhaseEnum &GetPhase() const
      {
        return m_Phase;
      }

      //! @brief returns the enum indicating the means by which the result is considered improved
      //! @return the means by which the result is considered improved
      const ImprovementTypeEnum &GetImprovementType() const
      {
        return m_ImprovementType;
      }

      //! @brief returns step status counts
      //! @return step status counts
      const storage::Vector< size_t> &GetCounts() const
      {
        return m_Counts;
      }

      //! @brief return last step status
      //! @return last step status
      const StepStatusEnum &GetStatusOfLastStep() const
      {
        return m_StatusOfLastStep;
      }

      //! @brief return number of steps in a row of specified STEP_STATUS
      //! @param STEP_STATUS StepStatus of interest
      //! @return number of steps in a row of specified STEP_STATUS
      size_t GetNumberStepsInARow( const StepStatus STEP_STATUS) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns the number of iterations since the last improvement
      //! @return number of iterations since the last improvement
      size_t GetIterationsSinceLastImprovement() const;

      //! @brief resets the tracker
      virtual void Reset();

      //! @brief sets the phase of the iteration
      //! @PHASE the now current phase of the iteration
      void SetPhase( const Phase &PHASE);

      //! @brief Update the status of this tracker
      void Update( const StepStatus STEP_STATUS);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // template class TrackerBase< t_ArgumentType, t_ResultType>

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_TRACKER_BASE_H_
