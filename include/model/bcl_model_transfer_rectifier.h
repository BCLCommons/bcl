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

#ifndef BCL_MODEL_TRANSFER_RECTIFIER_H_
#define BCL_MODEL_TRANSFER_RECTIFIER_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_transfer_function_interface.h"
#include "linal/bcl_linal_vector_const_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TransferRectifier
    //! @brief This class is a linear rectifier above 0.  Below 0, it returns 0 identically
    //!
    //! @see @link example_model_transfer_rectifier.cpp @endlink
    //! @author mendenjl
    //! @date Apr 30, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API TransferRectifier :
      public TransferFunctionInterface
    {

    private:

    //////////
    // data //
    //////////

      //! indicates whether this class was registered with the Enumerated
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! alpha - coefficient multiplied by X if it is below zero
      float m_Alpha;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      TransferRectifier() :
        m_Alpha( 0.0)
      {
      }

      //! @brief copy constructor
      TransferRectifier *Clone() const;

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

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class TransferRectifier

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_TRANSFER_RECTIFIER_H_
