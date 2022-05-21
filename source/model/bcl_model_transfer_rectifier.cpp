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
#include "model/bcl_model_transfer_rectifier.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> TransferRectifier::s_Instance
    (
      util::Enumerated< TransferFunctionInterface>::AddInstance( new TransferRectifier())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief copy constructor
    TransferRectifier *TransferRectifier::Clone() const
    {
      return new TransferRectifier( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TransferRectifier::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the working x range of that function
    //! @return x math::Range this function works on
    const math::Range< float> &TransferRectifier::GetOutputRange() const
    {
      static const math::Range< float> s_OutputRange( 0, 1);
      return s_OutputRange;
    }

    //! @brief get the range between which the function has a significant derivative
    //! @return x math::Range this function works on
    //! @note this is set to the range at which dF(y) is == 1/25 max(dF(y))
    const math::Range< float> &TransferRectifier::GetDynamicOutputRange() const
    {
      static const math::Range< float> s_OutputRange( 0.0, 1.0);
      return s_OutputRange;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &TransferRectifier::GetAlias() const
    {
      static const std::string s_Name( "Rectifier");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief function for one argument and result - depends on x ( ARGUMENT_X)
    //! @param ARGUMENT_X the argument for the transfer function
    //! @return function value of sigmoid transfer function
    float TransferRectifier::F( const float &ARGUMENT_X) const
    {
      return ARGUMENT_X >= 0.0 ? ARGUMENT_X : m_Alpha * ARGUMENT_X;
    }

    //! @brief derivative for one argument and result - depends on x ( ARGUMENT_X) and F( x) ( ARGUMENT_Y)
    //! @param ARGUMENT_X the argument for the derivative of the transfer function
    //! @param ARGUMENT_Y the result of F( ARGUMENT_X), can simplify calculating the derivative
    //! @return function value of the derivative of sigmoid transfer function
    float TransferRectifier::dF( const float &ARGUMENT_X, const float &ARGUMENT_Y) const
    {
      return ARGUMENT_X > 0.0 ? 1.0 : m_Alpha;
    }

    //! @brief overloaded F for vectors, calls F(vector) of the base class
    //! @param ARGUMENT_X the vector of arguments for the transfer function
    //! @return function value of sigmoid transfer function for a vector of arguments
    linal::Vector< float> TransferRectifier::F( const linal::VectorConstInterface< float> &ARGUMENT_X) const
    {
      return TransferFunctionInterface::F( ARGUMENT_X);
    }

    //! @brief overloaded F for vectors, calls F(vector,vector) of the base class
    //! @param STORAGE the vector that will hold the result
    //! @param ARGUMENT_X the vector of arguments for the transfer function
    void TransferRectifier::F
    (
      linal::VectorInterface< float> &STORAGE,
      const linal::VectorConstInterface< float> &ARGUMENT_X
    ) const
    {
      return TransferFunctionInterface::F( STORAGE, ARGUMENT_X);
    }

    //! multiply the storage vector by dF; used to compute error derivative at prior hidden layers
    //! @param STORAGE the vector that will hold the result
    //! @param ARGUMENT_X the vector of arguments for the derivative of the transfer function
    //! @param ARGUMENT_Y the vector of results of F( ARGUMENT_X), can simplify calculating the derivative
    void TransferRectifier::MultiplyBydF
    (
      linal::VectorInterface< float> &STORAGE,
      const linal::VectorConstInterface< float> &ARGUMENT_X,
      const linal::VectorConstInterface< float> &ARGUMENT_Y
    ) const
    {
      BCL_Assert
      (
        STORAGE.GetSize() == ARGUMENT_X.GetSize()
        && STORAGE.GetSize() == ARGUMENT_Y.GetSize(),
        "storage and argument vectors must be the same size!"
      );
      const float *x( ARGUMENT_X.Begin());

      for
      (
        float *itr_storage( STORAGE.Begin()), *itr_storage_end( STORAGE.End());
        itr_storage != itr_storage_end;
        ++itr_storage, ++x
      )
      {
        *itr_storage *= dF( *x, 0.0);
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer TransferRectifier::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "f(x) = x * (x > 0 ? 1 : alpha)");
      parameters.AddInitializer
      (
        "",
        "multiplier for X below 0",
        io::Serialization::GetAgentWithRange( &m_Alpha, 0.0, 1.0)
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
