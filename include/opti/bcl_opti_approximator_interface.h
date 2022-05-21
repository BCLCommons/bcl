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

#ifndef BCL_OPTI_APPROXIMATOR_INTERFACE_H_
#define BCL_OPTI_APPROXIMATOR_INTERFACE_H_

// include the namespace header
#include "bcl_opti.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opti_tracker.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opti
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorInterface
    //! @brief This class provides an interface for a general approximation scheme that allows specialization for
    //! various problems and methods.
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
    class ApproximatorInterface :
      public virtual util::ObjectInterface
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief Clone function
      //! @return pointer to a new ApproximatorInterface< t_ArgumentType, t_ResultType>
      virtual ApproximatorInterface *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns the current approximation
      //! @return current argument result pair
      virtual const util::ShPtr< storage::Pair< t_ArgumentType, t_ResultType> > GetCurrentApproximation() const = 0;

      //! @brief conducts approximation
      virtual void Approximate()
      {
        while( CanContinue())
        {
          Next();
        }
      }

      //! @brief evaluates whether the approximation can continue
      //! @return true, if the approximation can continue - otherwise false
      //! This function needs to be public to allow for optimization wrappers
      virtual bool CanContinue() const = 0;

      //! @brief conducts the next approximation step and stores the approximation
      //! This function needs to be public to allow for optimization wrappers
      virtual void Next() = 0;

    }; // template class ApproximatorInterface< t_ArgumentType, t_ResultType>

  } // namespace opti
} // namespace bcl

#endif // BCL_OPTI_APPROXIMATOR_INTERFACE_H_
