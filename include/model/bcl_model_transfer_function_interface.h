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

#ifndef BCL_MODEL_TRANSFER_FUNCTION_INTERFACE_H_
#define BCL_MODEL_TRANSFER_FUNCTION_INTERFACE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TransferFunctionInterface
    //! @brief interface for neural network transfer functions.
    //! @details It provides a function, a dF function and information on the input and output range, so that data is normalized
    //! properly.
    //!
    //! @remarks example unnecessary
    //! @author woetzen
    //! @date Apr 9, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API TransferFunctionInterface :
      public util::SerializableInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new TransferFunctionInterface
      virtual TransferFunctionInterface *Clone() const = 0;

      //! @brief get the output y range of that function
      //! @return y Range this function returns
      virtual const math::Range< float> &GetOutputRange() const = 0;

      //! @brief get the range between which the function has a significant derivative
      //! @return x math::Range this function works on
      //! @note this is set to the range at which dF(y) is == 1/25 max(dF(y))
      virtual const math::Range< float> &GetDynamicOutputRange() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! function for one argument and result - depends on x ( ARGUMENT( 0)) and F( x) ( ARGUMENT( 1))
      virtual float F( const float &ARGUMENT_X) const = 0;

      //! derivative for one argument and result - depends on x ( ARGUMENT( 0)) and F( x) ( ARGUMENT( 1))
      virtual float dF( const float &ARGUMENT_X, const float &ARGUMENT_Y) const = 0;

      //! overloaded F for vectors, calls F(float)
      linal::Vector< float> F( const linal::VectorConstInterface< float> &ARGUMENT_X) const;

      //! overloaded F for vectors, calls F(float), uses storage provided by user
      void F
      (
        linal::VectorInterface< float> &STORAGE,
        const linal::VectorConstInterface< float> &ARGUMENT_X
      ) const;

      //! multiply the storage vector by dF; used to compute error derivative at prior hidden layers
      virtual void MultiplyBydF
      (
        linal::VectorInterface< float> &STORAGE,
        const linal::VectorConstInterface< float> &ARGUMENT_X,
        const linal::VectorConstInterface< float> &ARGUMENT_Y
      ) const = 0;

    ///////////////
    // operators //
    ///////////////

    }; // class TransferFunctionInterface

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_TRANSFER_FUNCTION_INTERFACE_H_
