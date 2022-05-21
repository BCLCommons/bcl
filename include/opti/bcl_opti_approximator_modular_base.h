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

#ifndef BCL_OPTI_APPROXIMATOR_MODULAR_BASE_H_
#define BCL_OPTI_APPROXIMATOR_MODULAR_BASE_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_approximator_modular_interface.h"
#include "bcl_opti_criterion_interface.h"
#include "bcl_opti_print_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorModularBase
    //! @brief This class provides basic functionality for a general approximation scheme that allows specialization
    //! for various problems and methods
    //! @detail ApproximatorModularBase provides a tracker that stores the results of the
    //! approximation process and controls printing of results.
    //!
    //! @tparam t_ArgumentType is the type of the approximation argument
    //! @tparam t_ResultType is the type of the approximation result
    //!
    //! @remarks example unnecessary
    //! @author fischea
    //! @date Dec 11, 2012
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class ApproximatorModularBase :
      public ApproximatorModularInterface< t_ArgumentType, t_ResultType>
    {

    //////////
    // data //
    //////////

    protected:

      //! shared pointer to the step printer
      util::ShPtr< PrintInterface< t_ArgumentType, t_ResultType> > m_Printer;

      //! shared pointer to a criterion which indicates whether the approximation should terminate
      util::ShPtr< CriterionInterface< t_ArgumentType, t_ResultType> > m_ShouldTerminate;

      //! the tracker used to track the approximation process
      util::OwnPtr< Tracker< t_ArgumentType, t_ResultType> > m_Tracker;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      ApproximatorModularBase() :
        m_Tracker( new Tracker< t_ArgumentType, t_ResultType>(), true)
      {
      }

      //! @brief construct from improvement type
      //! @param IMPROVEMENT_TYPE type of change to be considered as improvement
      ApproximatorModularBase( const ImprovementType &IMPROVEMENT_TYPE) :
        m_Tracker( new Tracker< t_ArgumentType, t_ResultType>( IMPROVEMENT_TYPE), true)
      {
      }

      ApproximatorModularBase( const Tracker< t_ArgumentType, t_ResultType> &TRACKER) :
        m_Tracker( TRACKER.Clone(), true)
      {
        m_Tracker->Reset();
      }

      //! @brief Clone function
      //! @return pointer to a new ApproximatorModularBase< t_ArgumentType, t_ResultType>
      virtual ApproximatorModularBase *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the tracker used
      //! @return the tracker used
      const Tracker< t_ArgumentType, t_ResultType> &GetTracker() const
      {
        return *m_Tracker;
      }

      //! @brief returns the tracker used
      //! @return the tracker used
      Tracker< t_ArgumentType, t_ResultType> &GetTracker()
      {
        return *m_Tracker;
      }

      //! @brief sets the printer
      //! @param PRINTER the printer to be used
      void SetPrinter( const PrintInterface< t_ArgumentType, t_ResultType> &PRINTER)
      {
        m_Printer = util::CloneToShPtr( PRINTER);
      }

      //! @brief sets the termination criterion
      //! @param CRITERION criterion to be used for termination
      void SetCriterion( const CriterionInterface< t_ArgumentType, t_ResultType> &CRITERION)
      {
        m_ShouldTerminate = util::CloneToShPtr( CRITERION);
      }

      //! @brief returns a shared pointer to the current approximation
      //! @return shared pointer to  the current argument result pair
      virtual const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > GetCurrentApproximation() const
      {
        return this->GetTracker().GetCurrent();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief conducts approximation
      virtual void Approximate()
      {
        this->GetTracker().SetPhase( e_Start);
        if( m_Printer.IsDefined())
        {
          m_Printer->Print( this->GetTracker());
        }
        PrintStartHook();

        // conduct approximation as long as the given termination criteria are not met
        this->GetTracker().SetPhase( e_Iteration);
        while( this->CanContinue() && this->ShouldContinue())
        {
          this->Next();
          if( m_Printer.IsDefined())
          {
            m_Printer->Print( this->GetTracker());
          }
          PrintStepHook();
        }

        this->GetTracker().SetPhase( e_End);
        if( m_Printer.IsDefined())
        {
          m_Printer->Print( this->GetTracker());
        }
        PrintEndHook();
      }

      //! @brief returns whether any of the termination criteria are met
      //! @return true, if any of the termination criteria are met
      bool ShouldContinue() const
      {
        return !m_ShouldTerminate.IsDefined() || !m_ShouldTerminate->CriteriaMet( this->GetTracker());
      }

    protected:

      //! @brief tracks the given model
      //! @param shared pointer to the model to keep track of
      //! @param STEP_STATUS status of the current step; leave as s_NumberStepStatus to determine from MODEL
      void Track
      (
        const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &MODEL,
        StepStatus STEP_STATUS = s_NumberStepStatus
      )
      {
        this->GetTracker().Track( MODEL, STEP_STATUS);
      }

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief prints start information concerning the approximation process
      //! @detail gets called at the beginning of the approximation process
      virtual void PrintStartHook()
      {
      }

      //! @brief prints step information concerning the approximation process
      //! @detail gets called at every approximation step
      virtual void PrintStepHook()
      {
      }

      //! @brief prints end information concerning the approximation process
      //! @detail gets called at the end of the approximation process
      virtual void PrintEndHook()
      {
      }

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_Printer, ISTREAM);
        io::Serialize::Read( m_ShouldTerminate, ISTREAM);
        io::Serialize::Read( this->GetTracker(), ISTREAM);

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
        io::Serialize::Write( m_Printer, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_ShouldTerminate, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( this->GetTracker(), OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class ApproximatorModularBase< t_ArgumentType, t_ResultType>

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_APPROXIMATOR_MODULAR_BASE_H_
