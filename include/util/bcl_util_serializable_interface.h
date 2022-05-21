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

#ifndef BCL_UTIL_SERIALIZABLE_INTERFACE_H_
#define BCL_UTIL_SERIALIZABLE_INTERFACE_H_

// include the namespace header
#include "bcl_util.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_util_data_type.h"
#include "bcl_util_functional_type.h"
#include "bcl_util_object_interface.h"
#include "io/bcl_io_serialization_interface.h"
#include "io/bcl_io_serializer.h"

// external includes - sorted alphabetically
#include <iostream>

namespace bcl
{
  namespace util
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SerializableInterface
    //! @brief an object with I/O carried out via parameters
    //! @details Parameterized objects need to provide functions GetAlias(), GetClassDescription(), and GetParameters()
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Oct 29, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SerializableInterface :
      public io::SerializationInterface
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone the SerializableInterface
      //! @return pointer to new SerializableInterface
      virtual SerializableInterface *Clone() const = 0;

      //! @brief virtual destructor
      virtual ~SerializableInterface();

    /////////////////
    // data access //
    /////////////////

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      virtual const std::string &GetAlias() const
      {
        return GetClassIdentifier();
      }

      //! @brief return the label for this object
      //! @param INCLUDE_DATA whether to include members defined as data, otherwise, include only initialization info
      //! @return a label for this object
      virtual ObjectDataLabel GetLabel( const bool &INCLUDE_DATA = false) const;

      //! @brief return the string for this object
      //! @param INCLUDE_DATA whether to include members defined as data, otherwise, include only initialization info
      //! @return the initialization string for this object
      virtual std::string GetString( const bool &INCLUDE_DATA = false) const;

      //! @brief determine the type that this type serializes to
      virtual DataType::Type GetSerializedType() const;

    //////////////////////
    // helper functions //
    //////////////////////

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief Try to read the members of this object from the given label
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      //! @note this function is overrideable (for the moment) due to legacy reasons, for compatibility with
      //!       util::Implementation, which must call SetMembersFromDataLabel
      bool TryRead( const ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

      //! @brief set the value of the corresponding member based on the label
      //! @param LABEL label that is used to set the string
      //! @param ERROR_STREAM stream to write errors to
      //! @return result of validation
      io::ValidationResult ValidateRead( const ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

      //! @brief Read the members of this object from the given label, assert on failure, error messages sent to logger
      //! @param LABEL the label containing members that should be read of this class
      void AssertRead( const ObjectDataLabel &LABEL);

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief write this help for this object
      //! @param WRITER a fixed line width writer
      virtual io::FixedLineWidthWriter &WriteHelp( io::FixedLineWidthWriter &WRITER) const;

      //! @brief write this help for this object
      //! @param STREAM stream to write to
      //! @param INDENT indentation to give each line
      virtual std::ostream &WriteHelp( std::ostream &STREAM, const size_t INDENT = 0) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      virtual io::Serializer GetSerializer() const = 0;

      //! @brief get the complete serializer, including command line identifier and type
      //! @return the complete serializer, including command line identifier and type
      virtual io::Serializer GetCompleteSerializer() const;

      //! These functions are sometimes overridden to initialize data members from initialization variables setup during read

      //! @brief a function that derived classes can override to perform some action on a class before any data members
      //!        are read, e.g. resetting certain data members so that a post-read command can update parameters that
      //!        were not set
      //! @param SERIALIZER The serializer that will be used
      //! @param ERR_STREAM stream to write any errors encountered to
      //! @return result of any validation performed internally
      virtual io::ValidationResult PreReadHook( const ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM);

      //! @brief hook that derived classes can override, called after TryRead successfully reads a serializer containing only initializer info (no data variables)
      //! @param SERIALIZER the serializer object with initialization information
      //! @param ERR_STREAM stream to write out errors to
      //! @return true, unless there were new errors
      virtual bool ReadInitializerSuccessHook( const ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM);

      //! @brief a function that derived classes can override to take additional actions whenever Read is called successfully
      //!        AND any data members are specified for the given class
      //! @param SERIALIZER the serializer object with initialization information
      //! @param ERR_STREAM stream to write out errors to
      //! @return true, unless there were new errors
      virtual bool ReadDataSuccessHook( const ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM);

      //! @brief Get a set of all class names used by the serializer. Useful for introspection
      //! @param TYPES set to insert type names into
      //! @param INCLUDE_OPTIONAL true to also count optional members
      //! @param INCLUDE_DATA true to also include data-containing members
      //! @param RECURSIVE whether to recurse into sub-objects
      void InsertDataTypes
      (
        storage::Map< std::string, size_t> &TYPES,
        const bool &INCLUDE_OPTIONAL = true,
        const bool &INCLUDE_DATA = false,
        const size_t &MAX_DEPTH = size_t( -1)
      ) const;

      //! @brief Get addresses of all objects serialized as part of this
      //! @notes used to ensure uniqueness of serialized objects
      std::vector< size_t> GetSerializationAddresses() const;

      //! @brief infer functional type of the serializable interface, relative to a specific interface
      //! @param INTERFACE_STR the type name string of the interface of interest (by default, any interface)
      //! @return the functional type of the serializable interface, relative to a specific interface
      FunctionalType::Type InferFunctionalType( const std::string &INTERFACE_STR = "util::Implementation<") const;

    }; // class SerializableInterface

    //! @brief operator >> to read ReadWrite derived Object from std::istream
    //! @param ISTREAM input stream the object is read from
    //! @param OBJECT the object derived from ReadWrite
    //! @return the input stream where the object was read from
    BCL_API std::istream &operator >>( std::istream &ISTREAM, SerializableInterface &OBJECT);

    //! @brief operator >> to read const ReadWrite derived Object from std::istream which will assert, since this is not possible
    //! @param ISTREAM input stream the object is read from
    //! @param OBJECT the object derived from ReadWrite
    //! @return the input stream where the object was read from
    BCL_API std::istream &operator >>( std::istream &ISTREAM, const SerializableInterface &OBJECT);

    //! @brief write ReadWrite derived Object to std::ostream
    //! @param OSTREAM output stream the object is written to
    //! @param OBJECT object derived from ReadWrite
    //! @return output stream the object was written to
    BCL_API std::ostream &operator <<( std::ostream &OSTREAM, const SerializableInterface &OBJECT);

  } // namespace util
} // namespace bcl

#endif // BCL_UTIL_SERIALIZABLE_INTERFACE_H_

