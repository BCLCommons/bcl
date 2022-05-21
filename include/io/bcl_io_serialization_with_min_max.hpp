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

#ifndef BCL_IO_SERIALIZATION_WITH_MIN_MAX_HPP_
#define BCL_IO_SERIALIZATION_WITH_MIN_MAX_HPP_

// include the header of this class
#include "bcl_io_serialization_with_min_max.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "type/bcl_type_enable_if.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_string_numeric_conversion.h"

// external includes - sorted alphabetically
#include <iostream>
#include <limits>

namespace bcl
{
  namespace io
  {

    //! @brief helper function to write the help for floating point types
    //! @param OSTREAM the stream to which the help is written to
    //! @note specialization for floating point t_Data types
    template< typename t_DataType>
    typename type::EnableIf< !std::numeric_limits< t_DataType>::is_integer, void>::Type
     WriteHelpImpl( const t_DataType &MIN, const t_DataType &MAX, std::ostream &OSTREAM)
    {
      // simple test for -infinity -> compare against the lowest number of that type
      if( MIN > -std::numeric_limits< t_DataType>::max())
      {
        OSTREAM << " >= " << MIN;
        if( MAX < std::numeric_limits< t_DataType>::max())
        {
          OSTREAM << " <= " << MAX;
        }
      }
      else if( MAX < std::numeric_limits< t_DataType>::max())
      {
        OSTREAM << " <= " << MAX;
      }
      else
      {
        OSTREAM << ", must be defined (not nan)";
      }
    }

    //! @brief helper function to write the help for integral types
    //! @param OSTREAM the stream to which the help is written to
    //! @note specialization for integral t_Data types
    template< typename t_DataType>
    typename type::EnableIf< std::numeric_limits< t_DataType>::is_integer, void>::Type
     WriteHelpImpl( const t_DataType &MIN, const t_DataType &MAX, std::ostream &OSTREAM)
    {
      // do not print out limits if they are already at the limits of the type
      if( MIN != std::numeric_limits< t_DataType>::min())
      {
        OSTREAM << " >= " << MIN;
      }
      if( MAX != std::numeric_limits< t_DataType>::max())
      {
        OSTREAM << " <= " << MAX;
      }
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from members
    //! @param MIN the minimum allowed character
    //! @param MAX the maximum allowed character
    template< typename t_DataType>
    SerializationWithMinMax< t_DataType>::SerializationWithMinMax
    (
      const t_DataType &MIN,
      const t_DataType &MAX
    ) :
      SerializationBuiltin< t_DataType>( NULL),
      m_MinValue( MIN),
      m_MaxValue( MAX)
    {
    }

    //! @brief constructor from members
    //! @param DATA reference to the associated member variable
    //! @param MIN the minimum allowed character
    //! @param MAX the maximum allowed character
    template< typename t_DataType>
    SerializationWithMinMax< t_DataType>::SerializationWithMinMax
    (
      const t_DataType *DATA,
      const t_DataType &MIN,
      const t_DataType &MAX
    ) :
      SerializationBuiltin< t_DataType>( DATA),
      m_MinValue( MIN),
      m_MaxValue( MAX)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SerializationWithMinMax
    template< typename t_DataType>
    SerializationWithMinMax< t_DataType> *SerializationWithMinMax< t_DataType>::Clone() const
    {
      return new SerializationWithMinMax< t_DataType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &SerializationWithMinMax< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief set the given object using the label
    //! @param OBJECT object to read
    //! @param LABEL label that is used to set the string
    //! @param ERROR_STREAM stream to write errors to
    //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
    template< typename t_DataType>
    bool SerializationWithMinMax< t_DataType>::TryReadObject
    (
      t_DataType &OBJECT,
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    ) const
    {
      // try to convert the value from a string, accounting for the allowed values
      return util::TryConvertFromString( OBJECT, m_MinValue, m_MaxValue, LABEL.GetValue(), ERR_STREAM);
    }

    //! @brief writes the help for the label
    //! @param OSTREAM the stream to which the help is written to
    //! @param INDENT the amount of indent
    //! @return the given stream to which the help was written to
    template< typename t_DataType>
    std::ostream &SerializationWithMinMax< t_DataType>::WriteHelp
    (
      std::ostream &OSTREAM,
      const size_t INDENT
    ) const
    {
      SerializationBuiltin< t_DataType>::WriteHelp( OSTREAM);
      WriteHelpImpl( m_MinValue, m_MaxValue, OSTREAM);
      return OSTREAM;
    }

  } // namespace io
} // namespace bcl

#endif // BCL_IO_SERIALIZATION_WITH_MIN_MAX_HPP_

