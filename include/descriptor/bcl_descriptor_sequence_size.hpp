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
#include "bcl_descriptor_sequence_size.h"
// includes from bcl - sorted alphabetically
#include "io/bcl_io_serializer.h"
#include "linal/bcl_linal_vector_reference.h"
#include "type/bcl_type_compare.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief Clone function
    //! @return pointer to new SequenceSize
    template< typename t_DataType, typename t_ReturnType>
    SequenceSize< t_DataType, t_ReturnType> *SequenceSize< t_DataType, t_ReturnType>::Clone() const
    {
      return new SequenceSize( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType, typename t_ReturnType>
    const std::string &SequenceSize< t_DataType, t_ReturnType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType, typename t_ReturnType>
    const std::string &SequenceSize< t_DataType, t_ReturnType>::GetAlias() const
    {
      static const std::string s_name( "NElements");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    //! @note only one value is needed if this is a numeric descriptor, for char descriptors, assume that 99999 is the
    //! @note max, so 5 characters is sufficient
    template< typename t_DataType, typename t_ReturnType>
    size_t SequenceSize< t_DataType, t_ReturnType>::GetNormalSizeOfFeatures() const
    {
      return type::Compare< t_ReturnType, float>::e_Same ? 1 : 5;
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
    io::Serializer SequenceSize< t_DataType, t_ReturnType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "The number of " + this->GetElementName() + " in the " + this->GetObjectName());

      return parameters;
    }

    //! @brief implementation of the calculate function for floating point numbers
    //! @param STORAGE storage for the result
    //! @param SIZE actual size of the sequence
    void CalculateImpl( linal::VectorReference< float> &STORAGE, const size_t &SIZE)
    {
      STORAGE( 0) = SIZE;
    }

    //! @brief implementation of the calculate function for characters
    //! @param STORAGE storage for the result
    //! @param SIZE actual size of the sequence
    void CalculateImpl( linal::VectorReference< char> &STORAGE, const size_t &SIZE)
    {
      const std::string size_str( util::Format()( SIZE));

      STORAGE = ' ';
      if( size_str.size() > STORAGE.GetSize())
      {
        BCL_MessageCrt
        (
          "Sequence size string " + size_str + " was too long to fit into the given field,"
          "increase the constant in GetNormalSizeOfFeatures()"
        );
        std::copy( size_str.begin(), size_str.begin() + STORAGE.GetSize(), STORAGE.Begin());
      }
      else
      {
        std::copy( size_str.begin(), size_str.end(), STORAGE.Begin());
      }
    }

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType, typename t_ReturnType>
    void SequenceSize< t_DataType, t_ReturnType>::Calculate( linal::VectorReference< t_ReturnType> &STORAGE)
    {
      CalculateImpl( STORAGE, Base< t_DataType, t_ReturnType>::GetCurrentObject()->GetSize());
    }

  } // namespace descriptor
} // namespace bcl
