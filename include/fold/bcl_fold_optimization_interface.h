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

#ifndef BCL_FOLD_OPTIMIZATION_INTERFACE_H_
#define BCL_FOLD_OPTIMIZATION_INTERFACE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{ namespace fold
  {

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OptimizationInterface
    //! @brief Interface for optimization methods working of protein models
    //! @details Implementing classes encapsulate specific, logically coherent  parts of the approximation process,
    //! like topology sampling and loop construction.
    //!
    //! @remarks example unnecessary
    //! @author fischea
    //! @date Mar 08, 2016
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API OptimizationInterface :
      public util::ObjectInterface
    {

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

    public:

      //! @brief virtual copy constructor
      //! @return pointer to a new OptimizationInterface
      virtual OptimizationInterface *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief initiates the optimization and returns a shared pointer to the optimization result
      //! @param MODEL protein model to be optimized
      //! @return shared pointer to the optimization result
      virtual util::ShPtr< assemble::ProteinModel> Optimize( const assemble::ProteinModel &MODEL) const = 0;

    }; // class OptimizationInterface

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_OPTIMIZATION_INTERFACE_H_
