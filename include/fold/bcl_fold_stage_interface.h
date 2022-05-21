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

#ifndef BCL_FOLD_STAGE_INTERFACE_H_
#define BCL_FOLD_STAGE_INTERFACE_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{ namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StageInterface
    //! @brief Interface for stages
    //! @details Stages encapsulate specific parts of the approximation process.
    //!
    //! @remarks example unnecessary
    //! @author pinojc
    //! @date 12/10/2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StageInterface :
      public util::ObjectInterface
    {

    ////////////////////////////////////
    //  construction and destruction  //
    ////////////////////////////////////

    public:

      //! @brief virtual copy constructor
      //! @return pointer to a new StageInterface
      virtual StageInterface *Clone() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief initiates the approximation and returns a shared pointer to the approximation result
      //! @param SP_ARGUMENT shared pointer to the approximation argument
      //! @return shared pointer to the approximation result
      virtual util::ShPtr< assemble::ProteinModel> Approximate( util::ShPtr< assemble::ProteinModel> &SP_ARGUMENT) const = 0;

      //! @brief initializes the stage with the current round number and stage number
      //! @param ROUND_NUMBER number of the model currently created
      //! @param STAGE_NUMBER number of the stage
      //! @param IS_LAST_STAGE indicates if it is the last stage
      virtual void Initialize( const size_t &ROUND_NUMBER, const size_t &STAGE_NUMBER, const bool &IS_LAST_STAGE) = 0;

    }; // class StageInterface

  } // namespace fold
} // namespace bcl
  
#endif // BCL_FOLD_STAGE_INTERFACE_H_
