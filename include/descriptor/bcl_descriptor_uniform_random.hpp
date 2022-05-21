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

// include header of this class
#include "bcl_descriptor_uniform_random.h"
// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief constructor from bool; whether all numbers returned in a given storage should be 1
    //! @param SAME true if all the numbers returned per call to operator() should be the same
    //! @param SIZE number of uniform random numbers to return per call
    template< typename t_DataType>
    UniformRandom< t_DataType>::UniformRandom( const bool &SAME, const size_t &SIZE) :
      m_Size( SIZE),
      m_Same( SAME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new UniformRandom
    template< typename t_DataType>
    UniformRandom< t_DataType> *UniformRandom< t_DataType>::Clone() const
    {
      return new UniformRandom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &UniformRandom< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &UniformRandom< t_DataType>::GetAlias() const
    {
      static const std::string s_name( "UniformRandom"), s_copy_name( "UniformRandomDuplicated");
      return m_Same ? s_copy_name : s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t UniformRandom< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return m_Size;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer UniformRandom< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Returns a uniform-randomly chosen "
        + std::string( m_Same ? " value, copied {size} times" : "set of values")
      );
      parameters.AddInitializer
      (
        "size",
        "number of values to return per " + Base< t_DataType, float>::GetObjectName(),
        io::Serialization::GetAgent( &m_Size)
      );

      return parameters;
    }

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType>
    void UniformRandom< t_DataType>::Calculate( linal::VectorReference< float> &STORAGE)
    {
      if( m_Same)
      {
        STORAGE = random::GetGlobalRandom().Double();
      }
      else
      {
        for( size_t i( 0); i < m_Size; ++i)
        {
          STORAGE( i) = random::GetGlobalRandom().Double();
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace descriptor
} // namespace bcl
