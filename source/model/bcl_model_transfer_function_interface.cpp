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
#include "model/bcl_model_transfer_function_interface.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //! @brief overloaded F for vectors, calls F(float)
    //! @param ARGUMENT_X the vector of arguments for the derivative of the transfer function
    //! @return function value of the transfer function for a vector of arguments
    linal::Vector< float> TransferFunctionInterface::F( const linal::VectorConstInterface< float> &ARGUMENT_X) const
    {
      linal::Vector< float> storage( ARGUMENT_X.GetSize());
      this->F( storage, ARGUMENT_X);
      return storage;
    }

    //! overloaded F for vectors, calls F(float), uses storage provided by user
    void TransferFunctionInterface::F
    (
      linal::VectorInterface< float> &STORAGE,
      const linal::VectorConstInterface< float> &ARGUMENT_X
    ) const
    {
      BCL_Assert( STORAGE.GetSize() == ARGUMENT_X.GetSize(), "storage and argument vectors must be the same size!");
      const float *x( ARGUMENT_X.Begin());
      for
      (
        float *itr_storage( STORAGE.Begin()), *itr_storage_end( STORAGE.End());
        itr_storage != itr_storage_end;
        ++itr_storage, ++x
      )
      {
        *itr_storage = this->F( *x);
      }
    }

  } // namespace model
} // namespace bcl
