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

#ifndef BCL_IO_SERIALIZATION_BUILTIN_HPP_
#define BCL_IO_SERIALIZATION_BUILTIN_HPP_

// include the header of this class
#include "bcl_io_serialization_builtin.h"

// includes from bcl - sorted alphabetically
#include "bcl_io_serialize.h"
#include "type/bcl_type_enable_if.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_string_numeric_conversion.h"

// external includes - sorted alphabetically
#include <limits>

namespace bcl
{
  namespace io
  {

    //! @brief get label for a numeric type
    //! @param DATA the number to get as a label
    //! @return an object data label
    template< typename t_DataType>
    util::ObjectDataLabel GetBasicLabel( const t_DataType &DATA)
    {
      std::ostringstream writer;
      Serialize::Write( DATA, writer);
      return util::ObjectDataLabel( "", writer.str());
    }

    //! @brief get label for a character type
    //! @param DATA the char to get as a label
    //! @return an object data label
    util::ObjectDataLabel GetBasicLabel( const char &DATA)
    {
      return util::ObjectDataLabel( "", std::string( 1, DATA));
    }

    //! @brief get label for a string type
    //! @param DATA the string to get as a label
    //! @return an object data label
    util::ObjectDataLabel GetBasicLabel( const std::string &DATA)
    {
      return util::ObjectDataLabel( "", DATA);
    }

    //! @brief get label for a string type
    //! @param DATA the string to get as a label
    //! @return an object data label
    util::ObjectDataLabel GetBasicLabel( const util::ObjectDataLabel &DATA)
    {
      return DATA;
    }

    //! @brief get the type of object represented by DATA
    //! @param DATA the data type of interest
    //! @return the corresponding DataType::Type
    template< typename t_DataType>
    typename type::EnableIf
    <
      !std::numeric_limits< t_DataType>::is_integer,
      util::DataType::Type
    >::Type GetType( const t_DataType *DATA)
    {
      return util::DataType::e_Float;
    }

    //! @brief get the type of object represented by DATA
    //! @param DATA the data type of interest
    //! @return the corresponding DataType::Type
    template< typename t_DataType>
    typename type::EnableIf
    <
      std::numeric_limits< t_DataType>::is_integer
      && std::numeric_limits< t_DataType>::is_signed,
      util::DataType::Type
    >::Type GetType( const t_DataType *DATA)
    {
      return util::DataType::e_Int;
    }

    //! @brief get the type of object represented by DATA
    //! @param DATA the data type of interest
    //! @return the corresponding DataType::Type
    template< typename t_DataType>
    typename type::EnableIf
    <
      std::numeric_limits< t_DataType>::is_integer
      && !std::numeric_limits< t_DataType>::is_signed,
      util::DataType::Type
    >::Type GetType( const t_DataType *DATA)
    {
      return util::DataType::e_UnsignedInt;
    }

    //! @brief overload of GetType for boolean types
    util::DataType::Type GetType( const bool *DATA)
    {
      return util::DataType::e_Bool;
    }

    //! @brief overload of GetType for character types
    util::DataType::Type GetType( const char *DATA)
    {
      return util::DataType::e_Char;
    }

    //! @brief overload of GetType for string
    util::DataType::Type GetType( const std::string *DATA)
    {
      return util::DataType::e_String;
    }

    //! @brief overload of GetType for serialization tree
    util::DataType::Type GetType( const util::ObjectDataLabel *DATA)
    {
      return util::DataType::e_DynamicObject;
    }

    //! @brief get help for a type
    //! @param DATA the type to get help for
    //! @param OSTREAM output stream for the help
    //! @note overload for signed integral types
    template< typename t_DataType>
    typename type::EnableIf
    <
      std::numeric_limits< t_DataType>::is_integer
      && std::numeric_limits< t_DataType>::is_signed,
      void
    >::Type
    WriteHelpImpl( const t_DataType *DATA, std::ostream &OSTREAM)
    {
      OSTREAM << "Any integer";
    }

    //! @brief get help for a type
    //! @param DATA the type to get help for
    //! @param OSTREAM output stream for the help
    //! @note overload for unsigned integral types
    template< typename t_DataType>
    typename type::EnableIf
    <
      std::numeric_limits< t_DataType>::is_integer
      && !std::numeric_limits< t_DataType>::is_signed,
      void
    >::Type
    WriteHelpImpl( const t_DataType *DATA, std::ostream &OSTREAM)
    {
      OSTREAM << "Any non-negative integer";
    }

    //! @brief get help for a type
    //! @param DATA the type to get help for
    //! @param OSTREAM output stream for the help
    //! @note overload for floating point types
    template< typename t_DataType>
    typename type::EnableIf< !std::numeric_limits< t_DataType>::is_integer, void>::Type
    WriteHelpImpl( const t_DataType *DATA, std::ostream &OSTREAM)
    {
      OSTREAM << "Any decimal (floating-point) value";
    }

    //! @brief overload of WriteHelp for character type
    void WriteHelpImpl( const char *DATA, std::ostream &OSTREAM)
    {
      OSTREAM << "any letter or character";
    }

    //! @brief overload of WriteHelp for string type
    void WriteHelpImpl( const std::string *DATA, std::ostream &OSTREAM)
    {
      OSTREAM << "any string";
    }

    //! @brief overload of WriteHelp for object data label type
    void WriteHelpImpl( const util::ObjectDataLabel *DATA, std::ostream &OSTREAM)
    {
      OSTREAM << "any object serialization tree";
    }

    //! @brief set the value of the corresponding member based on the label
    //! @param DATA the data to try to read
    //! @param LABEL label that is used to set the string
    //! @param ERR_STREAM stream to write errors to
    //! @note overload for numeric types
    template< typename t_DataType>
    typename type::EnableIf< std::numeric_limits< t_DataType>::is_specialized, bool>::Type
      TryReadImpl( t_DataType &DATA, const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      // try to convert the value from a string
      return util::TryConvertFromString( DATA, LABEL.GetValue(), ERR_STREAM);
    }

    //! @brief set the value of the corresponding member based on the label
    //! @param DATA the data to try to read
    //! @param LABEL label that is used to set the string
    //! @param ERR_STREAM stream to write errors to
    //! @note overload for strings
    bool TryReadImpl( std::string &DATA, const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      // only call ToString if necessary since it might quote the values, which is undesireable
      DATA = LABEL.IsScalar() ? LABEL.GetValue() : LABEL.ToString();
      return true;
    }

    //! @brief set the value of the corresponding member based on the label
    //! @param DATA the data to try to read
    //! @param LABEL label that is used to set the string
    //! @param ERR_STREAM stream to write errors to
    //! @note overload for object data labels
    bool TryReadImpl( util::ObjectDataLabel &DATA, const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      DATA = LABEL;
      return true;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType>
    SerializationBuiltin< t_DataType>::SerializationBuiltin( const t_DataType *DATA) :
      SerializationBase< t_DataType>( DATA)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SerializationBuiltin
    template< typename t_DataType>
    SerializationBuiltin< t_DataType> *SerializationBuiltin< t_DataType>::Clone() const
    {
      return new SerializationBuiltin( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &SerializationBuiltin< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Get the label for an object of the given type
    //! @param OBJECT the object to get the label for
    //! @param WITH_DATA whether to include any data members, else, only include initialization members
    template< typename t_DataType>
    util::ObjectDataLabel SerializationBuiltin< t_DataType>::GetLabelForObject
    (
      const t_DataType &OBJECT,
      const bool &WITH_DATA
    ) const
    {
      return GetBasicLabel( OBJECT);
    }

    //! @brief determine the type of value that the handler expects to parse
    template< typename t_DataType>
    util::DataType::Type SerializationBuiltin< t_DataType>::GetSerializedType() const
    {
      return GetType( ( t_DataType *)( NULL));
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
    //! @param ERR_STREAM stream to write errors to
    template< typename t_DataType>
    bool SerializationBuiltin< t_DataType>::TryReadObject
    (
      t_DataType &OBJECT,
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    ) const
    {
      // try to convert the value from a string
      return TryReadImpl( OBJECT, LABEL, ERR_STREAM);
    }

    //! @brief writes the help for the label
    //! @param OSTREAM the stream to which the help is written to
    //! @param INDENT the amount of indent
    //! @return the given stream to which the help was written to
    template< typename t_DataType>
    std::ostream &SerializationBuiltin< t_DataType>::WriteHelp( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // even though t_DataType has a defined range of values, it is exceedingly unlikely that the actual range
      // is useful to the user at this point.  If a value outside the desired range is entered, it will be recognized
      // and reported to the user later on
      WriteHelpImpl( ( t_DataType *)( NULL), OSTREAM);
      return OSTREAM;
    }

  } // namespace io
} // namespace bcl

#endif // BCL_IO_SERIALIZATION_BUILTIN_HPP_

