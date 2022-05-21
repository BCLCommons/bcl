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

#ifndef BCL_IO_SERIALIZER_H_
#define BCL_IO_SERIALIZER_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_io_serialization_interface.h"
#include "util/bcl_util_own_ptr.h"

// external includes - sorted alphabetically
#include <list>
#include <map>

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Serializer
    //! @brief class that holds objects to initialize members of a class
    //!
    //! @see @link example_io_initializer.cpp @endlink
    //! @author mendenjl
    //! @date Aug 10, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Serializer :
      public util::ObjectInterface
    {

    public:
    //////////
    // data //
    //////////

      //! 2-value enum to distinguish between initialization vs. data members of a class
      //! Initialization members:
      //!   - are constant until the object is re-initialized or goes out of scope
      //!   - are always written out via an object's GetString function
      //!   - Must be given when reading in an object OR have a default
      //! Data members:
      //!   - usually change during an object's lifetime, e.g. they are dynamic state variables as opposed to constants
      //!   - are only written out when calling GetString( true)
      //!   - May be completely omitted when initializing an object, but, if any of them are given, all of them must be
      //!   - be given
      enum Type
      {
        e_Initializer, //!< See enum comment
        e_Data         //!< See enum comment
      };

    private:

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class MemberInfo
      //! @brief data class that holds information about each member
      //!
      //! @remarks example unnecessary
      //! @author mendenjl
      //! @date Aug 10, 2012
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct MemberInfo
      {
        enum State //!< current state of the member's initialization
        {
          e_WasSet,    //!< parameter was attempted to be set
          e_Optional,  //!< parameter is optional and has not been set
          e_Default,   //!< parameter has a default value given but has not been set
          e_Required   //!< parameter that is still required
        };

      private:

        util::OwnPtr< SerializationInterface> m_Handler;     //!< Object that handles getting / setting the member
        std::string                           m_Description; //!< Member description
        std::string                           m_Default;     //!< Default string for object
        Type                                  m_Type;        //!< Type of parameter
        State                                 m_State;       //!< True if the parameter is neither required nor has default

      public:

        //! @brief default constructor, can optionally take handler, description, and type
        MemberInfo
        (
          const util::OwnPtr< SerializationInterface> &HANDLER = util::OwnPtr< SerializationInterface>(),
          const Type &TYPE = e_Initializer,
          const std::string &DESCRIPTION = std::string()
        );

        //! @brief get the description
        //! @return the description for this parameter
        const std::string &GetDescription() const
        {
          return m_Description;
        }

        //! @brief get the current state
        //! @return the state of this parameter
        State GetState() const
        {
          return m_State;
        }

        //! @brief get the type
        //! @return the type of this parameter
        Type GetType() const
        {
          return m_Type;
        }

        //! @brief get the label handler
        //! @return the label handler
        const SerializationInterface &GetHandler() const
        {
          return *m_Handler;
        }

        //! @brief get the label handler
        //! @return the label handler
        const util::OwnPtr< SerializationInterface> &GetHandlerPtr() const
        {
          return m_Handler;
        }

        //! @brief get the default value (if applicable)
        //! @return the default value (if applicable)
        const std::string &GetDefault() const
        {
          return m_Default;
        }

        //! @brief set the default value
        //! @param DEFAULT the default value to use
        void SetDefault( const std::string &DEFAULT);

        //! @brief change a required parameter into being optional
        void SetAsOptional();

        //! @brief setup with a tree
        //! @param NAME name of the parameter attempted to be set
        //! @param TREE tree of values to use to set the member
        //! @param ERR_STREAM error stream for output
        //! @return validation result; allowed on success, invalid if not possible,
        //!         e.g. if the parameter was non-optional or had no default
        ValidationResult Set( const std::string &NAME, const util::ObjectDataLabel &TREE, std::ostream &ERR_STREAM);

        //! @brief get the label for the underlying object
        //! @param DETAIL whether to output detailed labels
        //! @return the label for the underlying object
        util::ObjectDataLabel GetLabel( const bool &DETAIL) const;

        //! @brief finalize the parameter -> set default values up if they exist
        //! @param NAME name of the parameter attempted to be set
        //! @param ERR_STREAM error stream for output
        //! @return true if the final state is e_WasSet or e_IsOptional
        bool Finalize( const std::string &NAME, std::ostream &ERR_STREAM);

        //! @brief writes the help for the ioline
        //! @param OSTREAM the stream to which the help is written to
        //! @return the given stream to which the help was written to
        FixedLineWidthWriter &WriteHelp( FixedLineWidthWriter &OSTREAM) const;

      };

      std::string                         m_CmdLineIdentifier; //!< Class id when used in polymorphic context
      std::string                         m_ClassDescription;  //!< Class description, used for help
      std::map< std::string, MemberInfo>  m_HandlerMap;        //!< Map from name to member info
      size_t                              m_NumberInits;       //!< Count of the total # of initializers
      size_t                              m_NumberData;        //!< Count of the total # of data members
      bool                                m_ReadData;          //!< True if any data members were read
      bool                                m_InlineSelf;        //!< True if the serializer is inlined
      bool                                m_InlineArgs;        //!< True if there is a single parameter to this object that should be inlined
      util::DataType::Type                m_Type;              //!< Type of data serialized

      //!< List of members in the order that they were added
      std::list< std::map< std::string, MemberInfo>::const_iterator> m_MemberItrs;

      //! @brief undefined assignment operator (undefined to discourage use to circumvent const-correctness)
      Serializer &operator=( const Serializer &);

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Serializer();

      //! @brief copy constructor
      Serializer( const Serializer &SERIALIZER);

      //! @brief virtual copy constructor
      //! @return pointer to Serializer object
      Serializer *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of this class as is used on the command line
      //! @return the name of this class as is used on the command line
      const std::string &GetCommandLineIdentifier() const;

      //! @brief sets the command line identifier
      //! @param ID the new identifier
      //! @return a reference to this
      Serializer &SetCommandLineIdentifier( const std::string &ID);

      //! @brief get the current classes' description
      //! @return the current classes description
      const std::string &GetClassDescription();

      //! @brief sets class description, for use in help
      //! @param DESCRIPTION the new description
      Serializer &SetClassDescription( const std::string &DESCRIPTION);

      //! @brief get whether data was read during the last ReadArguments operation
      //! @return true if data was read during the last ReadArguments operation
      bool GetWasDataRead() const
      {
        return m_ReadData;
      }

      //! @brief get the type of data serialized
      //! @return the type of data serialized
      util::DataType::Type GetType() const
      {
        return m_Type;
      }

      //! @brief Get a set of all class names used by the serializer. Useful for introspection
      //! @param TYPES set to insert type names into
      //! @param INCLUDE_OPTIONAL true to also count optional members
      //! @param INCLUDE_DATA true to also include data-containing members
      void InsertDataTypes
      (
        storage::Map< std::string, size_t> &TYPES,
        const bool &INCLUDE_OPTIONAL = true,
        const bool &INCLUDE_DATA = false,
        const size_t &MAX_DEPTH = size_t( -1)
      ) const;

      //! @brief Test whether this serializer has already been finalized, such that no further initializers can be added
      //! @return true if this serializer has already been finalized, such that no further initializers can be added
      bool IsFinalized() const
      {
        return m_Type != util::DataType::s_NumberTypes;
      }

      //! @brief sets the type and finalizes this object (e.g. such that no further parameters can be added)
      //! @param TYPE the actual type
      Serializer &SetTypeAndFinalize( const util::DataType::Type &TYPE);

    ////////////////
    // operations //
    ////////////////

      //! @brief add a member initializer
      //! @param NAME name for the data member
      //! @param DESCRIPTION description for the parameter
      //! @param INIT Serializer to push back
      void AddInitializer
      (
        const std::string &NAME,
        const std::string &DESCRIPTION,
        const util::OwnPtr< SerializationInterface> &INIT
      );

      //! @brief add an optional member initializer
      //! @param NAME name for the data member
      //! @param DESCRIPTION description for the parameter
      //! @param INIT Serializer to push back
      void AddOptionalInitializer
      (
        const std::string &NAME,
        const std::string &DESCRIPTION,
        const util::OwnPtr< SerializationInterface> &INIT
      );

      //! @brief add a member initializer with a default value, to be used if the parameter was omitted
      //! @param NAME name for the data member
      //! @param DESCRIPTION description for the parameter
      //! @param INIT member handler
      //! @param DEFAULT default value for the data member
      void AddInitializer
      (
        const std::string &NAME,
        const std::string &DESCRIPTION,
        const util::OwnPtr< SerializationInterface> &INIT,
        const std::string &DEFAULT
      );

      //! @brief add a data member
      //! @param NAME name for the data member
      //! @param DATA data member to push back
      void AddDataMember
      (
        const std::string &NAME,
        const util::OwnPtr< SerializationInterface> &DATA
      );

      //! @brief add a data member with a default value that will to be used if any other data members were given
      //! @param NAME name for the data member
      //! @param DATA data member handler
      //! @param DEFAULT default value for the data member
      void AddDataMember
      (
        const std::string &NAME,
        const util::OwnPtr< SerializationInterface> &DATA,
        const std::string &DEFAULT
      );

      //! @brief set a particular initializer's default. Useful in subclasses for which certain parameters are unnecessary
      //! @param NAME name of the parameter
      //! @param VALUE value to set the default parameter with
      void SetDefault( const std::string &NAME, const std::string &VALUE);

      //! @brief merge another serializer into this serializer. Useful only for statically-typed, non-container classes
      //! @param SERIALIZER the serializer for a statically-typed object; is consumed by this operation
      //! @note asserts if the serializer's type is not static object or if there was a name collision
      void Merge( const Serializer &SERIALIZER);

      //! @brief return the initialization label for this object
      //! @param WITH_DATA whether to include data members, otherwise, only initialization members should be included
      //! @return the initialization label for this object
      util::ObjectDataLabel GetLabel( const bool &WITH_DATA = false) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from an object data label
      //! @param LABEL - the list of arguments
      //! @param ERROR_STREAM - the error stream to which errors should be written
      //! @return result of all internal validation
      ValidationResult ReadArguments( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

      //! @brief writes the help for the Serializer
      //! @param OSTREAM - the outstream to write to
      //! @param INDENT - the indent to use
      //! @return std::ostream &OSTREAM - return the stream after writing
      std::ostream &WriteHelp( std::ostream &OSTREAM, const size_t &INDENT) const;

      //! @brief writes the help for the Serializer
      //! @param OSTREAM - the outstream to write to
      //! @return return the stream after writing
      FixedLineWidthWriter &WriteHelp( FixedLineWidthWriter &OSTREAM) const;

      //! @brief writes the brief help for the Serializer
      //! @param OSTREAM - the outstream to write to
      //! @return return the stream after writing
      FixedLineWidthWriter &WriteBriefHelp( FixedLineWidthWriter &OSTREAM) const;

      //! @brief read Serializer for std::istream ISTREAM
      //! @param ISTREAM - the instream to read from
      //! @return std::istream &ISTREAM - return the stream after reading
      std::istream &Read( std::istream &ISTREAM);

      //! @brief writes Serializer to ostream OSTREAM
      //! @param OSTREAM - the outstream to write to
      //! @param INDENT indentation
      //! @return std::ostream &OSTREAM - return the stream after writing
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief Get addresses of all objects serialized as part of this
      //! @notes used to ensure uniqueness of serialized objects
      std::vector< size_t> GetSerializationAddresses() const;

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief add an additional data or member initializer
      //! @param NAME name of the member initializer
      //! @param MEMBER_INFO information about tha
      //! @param DESCRIPTION description for the parameter
      //! @param HANDLER object to handle parsing / setting this member
      //! @param TYPE the type of object
      void AddMember( const std::string &NAME, const MemberInfo &MEMBER_INFO);

      //! @brief function to get the keys from a map as a storage::Vector
      //! @param MAP the map from which to retrieve the keys
      //! @return get the keys from a map as a storage::Vector
      friend storage::Vector< std::string> GetKeys( const std::map< std::string, Serializer::MemberInfo> &MAP);

    }; //class Serializer

  } // namespace io
} // namespace bcl

#endif //BCL_IO_SERIALIZER_H_
