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

#ifndef BCL_OPTI_APPROXIMATOR_MODULAR_INTERFACE_H_
#define BCL_OPTI_APPROXIMATOR_MODULAR_INTERFACE_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_approximator_interface.h"
#include "bcl_opti_step_status.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorModularInterface
    //! @brief This class provides an interface for a general approximation scheme that allows specialization for
    //! various problems and methods. ApproximatorModularInterface provides an interface to keep track of the
    //! approximation process and control the printing of results.
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
    class ApproximatorModularInterface :
      public ApproximatorInterface< t_ArgumentType, t_ResultType>
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to a new ApproximatorModularInterface< t_ArgumentType, t_ResultType>
      virtual ApproximatorModularInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the tracker used
      //! @return the tracker used
      virtual const Tracker< t_ArgumentType, t_ResultType> &GetTracker() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief conducts approximation
      virtual void Approximate()
      {
        PrintStartHook();
        while( this->CanContinue() && this->ShouldContinue())
        {
          this->Next();
          Track( this->GetCurrentApproximation());
          PrintStepHook();
        }
        PrintEndHook();
      }

      //! @brief evaluates whether the approximation should continue
      //! @return true, if the approximation should continue - otherwise false
      virtual bool ShouldContinue() const = 0;

      //! @brief tracks the given model
      //! @param ShPtr to the model to keep track of
      virtual void Track
      (
        const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > &ARGUMENT,
        StepStatus STEP_STATUS = s_NumberStepStatus
      ) = 0;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief prints the start information
      virtual void PrintStartHook() = 0;

      //! @brief prints the step information
      virtual void PrintStepHook() = 0;

      //! @brief prints the end information
      virtual void PrintEndHook() = 0;

    }; // template class Tracker< t_ArgumentType, t_ResultType>

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_APPROXIMATOR_MODULAR_INTERFACE_H_
