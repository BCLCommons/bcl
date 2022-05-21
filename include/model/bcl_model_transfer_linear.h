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

#ifndef BCL_MODEL_TRANSFER_LINEAR_H_
#define BCL_MODEL_TRANSFER_LINEAR_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_transfer_function_interface.h"
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TransferLinear
    //! @brief This class is a Gaussian function suitable as transfer Function for ANN objects only
    //! @details The derivative given in dF is not the actual derivative but holds the condition dF( F( x)) == true
    //! derivative( x). This is required for all functions used as transfer Function in ANN.
    //!
    //! @see @link example_model_transfer_linear.cpp @endlink
    //! @author mendenjl
    //! @date Nov 18, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API TransferLinear :
      public TransferFunctionInterface
    {

    public:

    //////////
    // data //
    //////////

      //! static instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      TransferLinear();

      //! copy constructor
      TransferLinear *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the working x range of that function
      //! @return x math::Range this function works on
      const math::Range< float> &GetOutputRange() const;

      //! @brief get the range between which the function has a significant derivative
      //! @return x math::Range this function works on
      //! @note this is set to the range at which dF(y) is == 1/25 max(dF(y))
      const math::Range< float> &GetDynamicOutputRange() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief function for one argument and result - depends on x ( ARGUMENT_X)
      //! @param ARGUMENT_X the argument for the transfer function
      //! @return function value of sigmoid transfer function
      float F( const float &ARGUMENT_X) const;

      //! @brief derivative for one argument and result - depends on x ( ARGUMENT_X) and F( x) ( ARGUMENT_Y)
      //! @param ARGUMENT_X the argument for the derivative of the transfer function
      //! @param ARGUMENT_Y the result of F( ARGUMENT_X), can simplify calculating the derivative
      //! @return function value of the derivative of sigmoid transfer function
      float dF( const float &ARGUMENT_X, const float &ARGUMENT_Y) const;

      //! @brief overloaded F for vectors, calls F(vector) of the base class
      //! @param ARGUMENT_X the vector of arguments for the transfer function
      //! @return function value of sigmoid transfer function for a vector of arguments
      linal::Vector< float> F( const linal::VectorConstInterface< float> &ARGUMENT_X) const;

      //! @brief overloaded F for vectors, calls F(vector,vector) of the base class
      //! @param STORAGE the vector that will hold the result
      //! @param ARGUMENT_X the vector of arguments for the transfer function
      void F
      (
        linal::VectorInterface< float> &STORAGE,
        const linal::VectorConstInterface< float> &ARGUMENT_X
      ) const;

      //! multiply the storage vector by dF; used to compute error derivative at prior hidden layers
      //! @param STORAGE the vector that will hold the result
      //! @param ARGUMENT_X the vector of arguments for the derivative of the transfer function
      //! @param ARGUMENT_Y the vector of results of F( ARGUMENT_X), can simplify calculating the derivative
      void MultiplyBydF
      (
        linal::VectorInterface< float> &STORAGE,
        const linal::VectorConstInterface< float> &ARGUMENT_X,
        const linal::VectorConstInterface< float> &ARGUMENT_Y
      ) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class TransferLinear

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_TRANSFER_LINEAR_H_
