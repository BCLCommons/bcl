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
#include "util/bcl_util_serializable_interface.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual destructor
    SerializableInterface::~SerializableInterface()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief return the label for this object
    //! @param INCLUDE_DATA whether to include members defined as data, otherwise, include only initialization info
    //! @return a label for this object
    ObjectDataLabel SerializableInterface::GetLabel( const bool &INCLUDE_DATA) const
    {
      return GetCompleteSerializer().GetLabel( INCLUDE_DATA);
    }

    //! @brief return the string for this object
    //! @param INCLUDE_DATA whether to include members defined as data, otherwise, include only initialization info
    //! @return the initialization string for this object
    std::string SerializableInterface::GetString( const bool &INCLUDE_DATA) const
    {
      return GetLabel( INCLUDE_DATA).ToString();
    }

    //! @brief determine the type of value that can be parsed
    DataType::Type SerializableInterface::GetSerializedType() const
    {
      return DataType::e_StaticObject;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief Try to read the members of this object from the given label
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool SerializableInterface::TryRead( const ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      return ValidateRead( LABEL, ERROR_STREAM).IsAllowed();
    }

    //! @brief Read the members of this object from the given label, assert on failure, error messages sent to logger
    //! @param LABEL the label containing members that should be read of this class
    void SerializableInterface::AssertRead( const ObjectDataLabel &LABEL)
    {
      BCL_MessageDbg( "Trying to read " + LABEL.ToString() + " into " + this->GetClassIdentifier());
      std::stringstream err_stream;
      io::ValidationResult result( ValidateRead( LABEL, err_stream));
      if( !result.IsAllowed())
      {
        if( result.IsHelp())
        {
          BCL_ExitWithoutCallstack( "Requested help: " + err_stream.str(), 0);
        }
        else if( !result.IsComplete())
        {
          BCL_Exit( "Failed to read " + this->GetClassIdentifier() + ", error was: " + err_stream.str() + ", exiting", -1);
        }
      }
    }

    //! @brief set the value of the corresponding member based on the label
    //! @param LABEL label that is used to set the string
    //! @param ERROR_STREAM stream to write errors to
    //! @return result of validation
    io::ValidationResult SerializableInterface::ValidateRead
    (
      const ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
    )
    {
      if( LABEL.GetName() == io::ValidationResult::GetHelpString())
      {
        WriteHelp( ERROR_STREAM);
        return io::ValidationResult( io::e_Help);
      }

      // call PreReadHook: overridable function to reset any necessary members
      io::ValidationResult result( PreReadHook( LABEL, ERROR_STREAM));
      if( result.IsInvalid())
      {
        ERROR_STREAM << "Failed at pre-read hook";
        return result;
      }
      else if( result.IsComplete())
      {
        return io::ValidationResult( io::e_Allowed);
      }

      // get the serializer for this object
      io::Serializer serializer( GetCompleteSerializer());

      // read in arguments
      result = serializer.ReadArguments( LABEL, ERROR_STREAM);
      if( !result)
      {
        ERROR_STREAM << "Failed at reading arguments";
        // failed to read object correctly
        return result;
      }
      else if( serializer.GetWasDataRead())
      {
        // call function to finalize object if data members were read
        result = io::ValidationResult( ReadDataSuccessHook( LABEL, ERROR_STREAM));
        if( result.IsInvalid())
        {
          ERROR_STREAM << "Failed at ReadDataSuccessHook for " << this->GetClassIdentifier();
        }
        return result;
      }

      // only read initializers
      result = io::ValidationResult( ReadInitializerSuccessHook( LABEL, ERROR_STREAM));
      if( result.IsInvalid())
      {
        ERROR_STREAM << "Failed at ReadInitializerSuccessHook for " << this->GetClassIdentifier();
      }
      return result;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SerializableInterface::Read( std::istream &ISTREAM)
    {
      AssertRead( ObjectDataLabel( ISTREAM));
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &SerializableInterface::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return GetCompleteSerializer().Write( OSTREAM, INDENT);
    }

    //! @brief write this help for this object
    //! @param STREAM stream to write to
    //! @param INDENT indentation to give each line
    std::ostream &SerializableInterface::WriteHelp( std::ostream &STREAM, const size_t INDENT) const
    {
      return GetCompleteSerializer().WriteHelp( STREAM, INDENT);
    }

    //! @brief write this help for this object
    //! @param WRITER a fixed line width writer
    io::FixedLineWidthWriter &SerializableInterface::WriteHelp( io::FixedLineWidthWriter &WRITER) const
    {
      return GetCompleteSerializer().WriteHelp( WRITER);
    }

    //! @brief get the complete serializer, including command line identifier and type
    //! @return the complete serializer, including command line identifier and type
    io::Serializer SerializableInterface::GetCompleteSerializer() const
    {
      return GetSerializer().SetCommandLineIdentifier( GetAlias()).SetTypeAndFinalize( GetSerializedType());
    }

    //! These functions are sometimes overridden to initialize data members from initialization variables setup during read

    //! @brief a function that derived classes can override to perform some action on a class before any data members
    //!        are read, e.g. resetting certain data members so that a post-read command can update parameters that
    //!        were not set
    //! @param SERIALIZER The serializer that will be used
    //! @param ERR_STREAM stream to write any errors encountered to
    //! @return result of any validation performed internally
    io::ValidationResult SerializableInterface::PreReadHook( const ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
    {
      return io::ValidationResult( io::e_Allowed);
    }

    //! @brief a function that derived classes can override to set up additional data members whenever Read is called successfully
    //! @param SERIALIZER the serializer object with initialization information
    //! @param ERR_STREAM stream to write out errors to
    //! @return true, unless there were new errors
    bool SerializableInterface::ReadInitializerSuccessHook( const ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
    {
      return true;
    }

    //! @brief a function that derived classes can override to take additional actions whenever Read is called successfully
    //!        AND any data members are specified for the given class
    //! @param SERIALIZER the serializer object with initialization information
    //! @param ERR_STREAM stream to write out errors to
    //! @return true, unless there were new errors
    bool SerializableInterface::ReadDataSuccessHook( const ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM)
    {
      return true;
    }

    //! @brief Get a set of all class names used by the serializer. Useful for introspection
    //! @param TYPES set to insert type names into
    //! @param INCLUDE_OPTIONAL true to also count optional members
    //! @param INCLUDE_DATA true to also include data-containing members
    void SerializableInterface::InsertDataTypes
    (
      storage::Map< std::string, size_t> &TYPES,
      const bool &INCLUDE_OPTIONAL,
      const bool &INCLUDE_DATA,
      const size_t &MAX_DEPTH
    ) const
    {
      ++TYPES[ this->GetClassIdentifier()];
      this->GetCompleteSerializer().InsertDataTypes( TYPES, INCLUDE_OPTIONAL, INCLUDE_DATA, MAX_DEPTH);
    }

    //! @brief Get addresses of all objects serialized as part of this
    //! @notes used to ensure uniqueness of serialized objects
    std::vector< size_t> SerializableInterface::GetSerializationAddresses() const
    {
      if( this->GetSerializedType() == DataType::e_DynamicObject)
      {
        return io::SerializationInterface::GetSerializationAddresses();
      }
      return this->GetCompleteSerializer().GetSerializationAddresses();
    }

    //! @brief infer functional type of the serializable interface, relative to a specific interface
    //! @param INTERFACE_STR the type name string of the interface of interest (by default, any interface)
    //! @return the functional type of the serializable interface, relative to a specific interface
    FunctionalType::Type SerializableInterface::InferFunctionalType( const std::string &INTERFACE_STR) const
    {
      storage::Map< std::string, size_t> type_counts;
      this->InsertDataTypes( type_counts, true, false, 2);
      size_t count_match_iface( 0), count_unmatched_iface( 0);
      for
      (
        storage::Map< std::string, size_t>::const_iterator itr( type_counts.Begin()), itr_end( type_counts.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->first != this->GetClassIdentifier())
        {
          ( StartsWith( itr->first, INTERFACE_STR) ? count_match_iface : count_unmatched_iface) += itr->second;
        }
      }
      return FunctionalType::GetType( count_match_iface, count_unmatched_iface);
    }

    //! @brief operator >> to read ReadWrite derived Object from std::istream
    //! @param ISTREAM input stream the object is read from
    //! @param OBJECT the object derived from ReadWrite
    //! @return the input stream where the object was read from
    std::istream &operator >>( std::istream &ISTREAM, SerializableInterface &OBJECT)
    {
      return OBJECT.Read( ISTREAM);
    }

    //! @brief operator >> to read const ReadWrite derived Object from std::istream which will assert, since this is not possible
    //! @param ISTREAM input stream the object is read from
    //! @param OBJECT the object derived from ReadWrite
    //! @return the input stream where the object was read from
    std::istream &operator >>( std::istream &ISTREAM, const SerializableInterface &OBJECT)
    {
      BCL_Exit( "cannot read a const OBJECT", -1);
      return ISTREAM;
    }

    //! @brief write ReadWrite derived Object to std::ostream
    //! @param OSTREAM output stream the object is written to
    //! @param OBJECT object derived from ReadWrite
    //! @return output stream the object was written to
    std::ostream &operator <<( std::ostream &OSTREAM, const SerializableInterface &OBJECT)
    {
      return OBJECT.Write( OSTREAM, 0);
    }

  } // namespace util
} // namespace bcl
