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
#include "bcl_descriptor_constants.h"
// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief default constructor
    template< typename t_DataType, typename t_ReturnType>
    Constants< t_DataType, t_ReturnType>::Constants() :
      m_Description()
    {
    }

    //! @brief constructor from a description
    //! @param DESCRIPTION storage::Vector< t_ReturnType> that holds the constant description
    template< typename t_DataType, typename t_ReturnType>
    Constants< t_DataType, t_ReturnType>::Constants( const t_ReturnType &DESCRIPTION) :
      m_Description( size_t( 1), DESCRIPTION)
    {
    }

    //! @brief constructor from a description
    //! @param DESCRIPTION storage::Vector< t_ReturnType> that holds the constant description
    template< typename t_DataType, typename t_ReturnType>
    Constants< t_DataType, t_ReturnType>::Constants( const linal::Vector< t_ReturnType> &DESCRIPTION) :
      m_Description( DESCRIPTION)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Constants
    template< typename t_DataType, typename t_ReturnType>
    Constants< t_DataType, t_ReturnType> *Constants< t_DataType, t_ReturnType>::Clone() const
    {
      return new Constants( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType, typename t_ReturnType>
    const std::string &Constants< t_DataType, t_ReturnType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType, typename t_ReturnType>
    const std::string &Constants< t_DataType, t_ReturnType>::GetAlias() const
    {
      static const std::string s_name( "Constant");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType, typename t_ReturnType>
    size_t Constants< t_DataType, t_ReturnType>::GetNormalSizeOfFeatures() const
    {
      return m_Description.GetSize();
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
    template< typename t_DataType, typename t_ReturnType>
    io::Serializer Constants< t_DataType, t_ReturnType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Returns a constant set of values");
      parameters.AddInitializer( "", "", io::Serialization::GetAgent( &m_Description));

      return parameters;
    }

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType, typename t_ReturnType>
    void Constants< t_DataType, t_ReturnType>::Calculate( linal::VectorReference< t_ReturnType> &STORAGE)
    {
      STORAGE.CopyValues( m_Description);
    }

    //! @brief Clone function
    //! @return pointer to new Constants
    template< typename t_DataType>
    Constants< t_DataType, char> *Constants< t_DataType, char>::Clone() const
    {
      return new Constants( *this);
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &Constants< t_DataType, char>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &Constants< t_DataType, char>::GetAlias() const
    {
      static const std::string s_name( "String");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t Constants< t_DataType, char>::GetNormalSizeOfFeatures() const
    {
      return m_Description.size();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer Constants< t_DataType, char>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "A constant string");
      parameters.AddInitializer( "", "", io::Serialization::GetAgent( &m_Description));

      return parameters;
    }

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    //! @return true, if the recalculation was completed
    template< typename t_DataType>
    void Constants< t_DataType, char>::Calculate( linal::VectorReference< char> &STORAGE)
    {
      STORAGE.CopyValues( linal::VectorConstReference< char>( m_Description.size(), &m_Description[ 0]));
    }

  } // namespace descriptor
} // namespace bcl
