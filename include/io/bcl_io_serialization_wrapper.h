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

#ifndef BCL_IO_SERIALIZATION_WRAPPER_H_
#define BCL_IO_SERIALIZATION_WRAPPER_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_io_serialization_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SerializationWrapper
    //! @brief Wraps a type that derives from SerializationInterface
    //!
    //! @see @link example_io_serialization_wrapper.cpp @endlink
    //! @author mendenjl
    //! @date Nov 01, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_DataType>
    class SerializationWrapper :
      public SerializationBase< t_DataType>
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new SerializationWrapper
      SerializationWrapper *Clone() const
      {
        return new SerializationWrapper( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief Get the label for an object of the given type
      //! @param OBJECT the object to get the label for
      //! @param WITH_DATA whether to include any data members, else, only include initialization members
      util::ObjectDataLabel GetLabelForObject( const t_DataType &OBJECT, const bool &WITH_DATA) const
      {
        return OBJECT.GetLabel( WITH_DATA);
      }

      //! @brief determine the type of value that the handler expects to parse
      util::DataType::Type GetSerializedType() const
      {
        // cache the serialization type for faster access
        static util::DataType::Type s_type( t_DataType().GetSerializedType());
        return s_type;
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief writes the help for the label
      //! @param OSTREAM the stream to which the help is written to
      //! @param INDENT the amount of indent
      //! @return the given stream to which the help was written to
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t INDENT = 0) const
      {
        return t_DataType().WriteHelp( OSTREAM, INDENT);
      }

      //! @brief set the given object using the label
      //! @param OBJECT object to read
      //! @param LABEL label that is used to set the string
      //! @param ERROR_STREAM stream to write errors to
      //! @return bool - if the parameter was set (parameters are not set if they are not allowed parameters)
      bool TryReadObject( t_DataType &OBJECT, const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM) const
      {
        return OBJECT.TryRead( LABEL, ERR_STREAM);
      }

      //! @brief Get a set of all class names used by the serializer. Useful for introspection
      //! @param TYPES set to insert type names into
      //! @param INCLUDE_OPTIONAL true to also count optional members
      //! @param INCLUDE_DATA true to also include data-containing members
      virtual void InsertDataTypesForObject
      (
        const t_DataType &OBJECT,
        storage::Map< std::string, size_t> &TYPES,
        const bool &INCLUDE_OPTIONAL,
        const bool &INCLUDE_DATA,
        const size_t &MAX_DEPTH
      ) const
      {
        OBJECT.InsertDataTypes( TYPES, INCLUDE_OPTIONAL, INCLUDE_DATA, MAX_DEPTH);
      }

    }; // class SerializationWrapper

  } // namespace io
} // namespace bcl

#endif // BCL_IO_SERIALIZATION_WRAPPER_H_
