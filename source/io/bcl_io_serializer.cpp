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
#include "io/bcl_io_serializer.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_command_state.h"
#include "command/bcl_command_guesser.h"
#include "io/bcl_io_fixed_line_width_writer.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Serializer::s_Instance( GetObjectInstances().AddInstance( new Serializer()));

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Serializer::Serializer() :
      m_NumberInits( 0),
      m_NumberData( 0),
      m_ReadData( false),
      m_InlineSelf( false),
      m_InlineArgs( false),
      m_Type( util::DataType::s_NumberTypes)
    {
    }

    //! @brief copy constructor
    Serializer::Serializer( const Serializer &SERIALIZER) :
      m_CmdLineIdentifier( SERIALIZER.m_CmdLineIdentifier),
      m_ClassDescription ( SERIALIZER.m_ClassDescription ),
      m_HandlerMap       ( SERIALIZER.m_HandlerMap       ),
      m_NumberInits      ( SERIALIZER.m_NumberInits      ),
      m_NumberData       ( SERIALIZER.m_NumberData       ),
      m_ReadData         ( SERIALIZER.m_ReadData         ),
      m_InlineSelf       ( SERIALIZER.m_InlineSelf       ),
      m_InlineArgs       ( SERIALIZER.m_InlineArgs       ),
      m_Type             ( SERIALIZER.m_Type             )
    {
      // recreate the iterator list
      for
      (
        std::list< std::map< std::string, MemberInfo>::const_iterator>::const_iterator
          itr( SERIALIZER.m_MemberItrs.begin()), itr_end( SERIALIZER.m_MemberItrs.end());
        itr != itr_end;
        ++itr
      )
      {
        const std::map< std::string, MemberInfo>::const_iterator &itr_map( *itr);
        m_MemberItrs.push_back( m_HandlerMap.find( itr_map->first));
      }
    }

    //! @brief default constructor, can optionally take handler, description, and type
    Serializer::MemberInfo::MemberInfo
    (
      const util::OwnPtr< SerializationInterface> &HANDLER,
      const Type &TYPE,
      const std::string &DESCRIPTION
    ) :
      m_Handler( HANDLER),
      m_Description( DESCRIPTION),
      m_Default(),
      m_Type( TYPE),
      m_State( e_Required)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to Serializer object
    Serializer *Serializer::Clone() const
    {
      return new Serializer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Serializer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the default value
    //! @param DEFAULT the default value to use
    void Serializer::MemberInfo::SetDefault( const std::string &DEFAULT)
    {
      BCL_Assert( m_State == e_Required || m_State == e_Default, "Can only set a default value up for required parameters");
      m_Default = DEFAULT;
      m_State = e_Default;
    }

    //! @brief change a required parameter into being optional
    void Serializer::MemberInfo::SetAsOptional()
    {
      BCL_Assert( m_State == e_Required, "Can only set required parameters up to be optional");
      m_State = e_Optional;
    }

    //! @brief writes the help for the ioline
    //! @param OSTREAM the stream to which the help is written to
    //! @return the given stream to which the help was written to
    FixedLineWidthWriter &Serializer::MemberInfo::WriteHelp( FixedLineWidthWriter &OSTREAM) const
    {
      if( !m_Description.empty())
      {
        // write description
        OSTREAM << m_Description << ", ";
      }

      // if default parameter was given, also write this
      if( m_State == e_Default)
      {
        OSTREAM.WriteOnOneLine( "default: \"" + m_Default + "\", ");
      }
      else if( m_State == e_Optional)
      {
        OSTREAM << "optional, ";
      }
      std::ostringstream parameter_check_stream;
      const size_t actual_indent( OSTREAM.GetIndent());
      const size_t xtra_indent( OSTREAM.GetExtraIndent());
      m_Handler->WriteHelp( parameter_check_stream, actual_indent / 2);
      const std::string help_str( util::RStrip( parameter_check_stream.str(), ", \n\r"));
      const size_t first_line_length( std::min( help_str.size(), help_str.find( '\n')));
      if( OSTREAM.GetRemainingSpaceOnLine() < first_line_length)
      {
        OSTREAM.AddIndent( 2);
        OSTREAM.NewLineIndent();
        OSTREAM.PopIndent();
      }
      OSTREAM.SetIndent( 0);
      OSTREAM.SetAutoIndentExtra( xtra_indent + actual_indent);
      OSTREAM << help_str;
      OSTREAM.PopIndent();
      OSTREAM.SetAutoIndentExtra( xtra_indent);
      return OSTREAM;
    }

    //! @brief finalize the parameter -> set default values up if they exist
    //! @param NAME name of the parameter attempted to be set
    //! @param ERR_STREAM error stream for output
    //! @return validation result; allowed on success, invalid if not possible,
    //!         e.g. if the parameter was non-optional or had no default
    ValidationResult Serializer::MemberInfo::Set
    (
      const std::string &NAME,
      const util::ObjectDataLabel &TREE,
      std::ostream &ERR_STREAM
    )
    {
      m_State = e_WasSet;
      if
      (
        TREE.GetNumberArguments()
        && util::DataType::TestMustBeScalar( m_Handler->GetSerializedType())
      )
      {
        // should not be any arguments to a data/scalar value
        ERR_STREAM << "Arguments not allowed for scalar types [ in data label " << TREE.ToString() << "]; handled by "
                   << m_Handler->GetClassIdentifier();
        return ValidationResult( e_Invalid);
      }
      ValidationResult result( m_Handler->ValidateRead( TREE, ERR_STREAM));
      if( result.IsInvalid())
      {
        ERR_STREAM
          << " [ in "
          << util::DataType::GetTypeName( m_Handler->GetSerializedType())
          << " label for parameter: " << NAME
          << " handled by " << m_Handler->GetClassIdentifier()
          << " label was: "
          << TREE.ToNamedString()
          << "]";
      }
      return result;
    }

    //! @brief get the label for the underlying object
    //! @param DETAIL whether to output detailed labels
    //! @return the label for the underlying object
    util::ObjectDataLabel Serializer::MemberInfo::GetLabel( const bool &DETAIL) const
    {
      return m_Handler->GetLabel( DETAIL);
    }

    //! @brief set with the default value, if there was one
    //! @param NAME name of the parameter attempted to be set
    //! @param ERR_STREAM error stream for output
    //! @return true on success or if the parameter was already set, false if not possible because the state was not e_Default
    bool Serializer::MemberInfo::Finalize( const std::string &NAME, std::ostream &ERR_STREAM)
    {
      // optional and already set parameters are fine
      if( m_State == e_WasSet || m_State == e_Optional)
      {
        return true;
      }

      // setup default parameters
      if( m_State == e_Default)
      {
        // set the default parameter
        const bool status( Set( NAME, util::ObjectDataLabel( NAME, util::ObjectDataLabel( m_Default)), ERR_STREAM));
        m_State = e_Default; //< reset the state back to the proper value
        return status;
      }

      // error, not set
      ERR_STREAM << "\nparameter \"" << NAME << "\" was not given!";
      return false;
    }

    //! @brief get the name of this class as is used on the command line
    //! @return the name of this class as is used on the command line
    const std::string &Serializer::GetCommandLineIdentifier() const
    {
      return m_CmdLineIdentifier;
    }

    //! @brief sets the command line identifier
    //! @param ID the new identifier
    //! @return a reference to this
    Serializer &Serializer::SetCommandLineIdentifier( const std::string &ID)
    {
      m_CmdLineIdentifier = ID;
      return *this;
    }

    //! @brief get the current classes' description
    //! @return the current classes description
    const std::string &Serializer::GetClassDescription()
    {
      return m_ClassDescription;
    }

    //! @brief sets class description, for use in help
    //! @param DESCRIPTION the new description
    Serializer &Serializer::SetClassDescription( const std::string &DESCRIPTION)
    {
      m_ClassDescription = DESCRIPTION;
      return *this;
    }

    //! @brief Get a set of all class names used by the serializer. Useful for introspection
    //! @param TYPES set to insert type names into
    //! @param INCLUDE_OPTIONAL true to also count optional members
    //! @param INCLUDE_DATA true to also include data-containing members
    void Serializer::InsertDataTypes
    (
      storage::Map< std::string, size_t> &TYPES,
      const bool &INCLUDE_OPTIONAL,
      const bool &INCLUDE_DATA,
      const size_t &MAX_DEPTH
    ) const
    {
      if( !MAX_DEPTH)
      {
        return;
      }
      for
      (
        std::map< std::string, MemberInfo>::const_iterator itr( m_HandlerMap.begin()), itr_end( m_HandlerMap.end());
        itr != itr_end;
        ++itr
      )
      {
        if
        (
          !INCLUDE_OPTIONAL
          && ( itr->second.GetState() == MemberInfo::e_Default || itr->second.GetState() == MemberInfo::e_Optional)
        )
        {
          continue;
        }
        if( !INCLUDE_DATA && itr->second.GetType() == e_Data)
        {
          continue;
        }
        itr->second.GetHandler().InsertDataTypes
        (
          TYPES,
          INCLUDE_OPTIONAL,
          INCLUDE_DATA,
          MAX_DEPTH - 1
        );
      }
    }

    //! @brief sets the type and finalizes this object (e.g. such that no further parameters can be added)
    //! @param TYPE the actual type
    Serializer &Serializer::SetTypeAndFinalize( const util::DataType::Type &TYPE)
    {
      if( TYPE == util::DataType::e_StaticObject)
      {
        // it is always possible to inline calls to static objects with single parameters
        if( m_NumberInits == size_t( 1) && m_NumberData == size_t( 0))
        {
          // auto-inline
          m_Type = m_HandlerMap.begin()->second.GetHandler().GetSerializedType();
          m_InlineSelf = true;
          m_InlineArgs = true;
        }
        else
        {
          // no inlining possible, multiple parameters
          m_Type = TYPE;
          m_InlineSelf = false;
          m_InlineArgs = false;
        }
      }
      else if( !util::DataType::IsUnderlyingTypeKnown( TYPE))
      {
        if( m_InlineArgs)
        {
          // scalar types cannot be inlined into dynamic types
          if( util::DataType::TestMustBeScalar( m_Type) || !util::DataType::IsUnderlyingTypeKnown( m_Type))
          {
            // the sole parameter was not a type which is known, so it cannot be inlined
            m_InlineArgs = false;
          }
        }
        m_InlineSelf = false;
        m_Type = TYPE;
      }
      else
      {
        m_Type = TYPE;
        if( m_NumberInits == size_t( 1) && m_NumberData == size_t( 0))
        {
          m_InlineSelf = true;
          m_InlineArgs = true;
        }
      }
      return *this;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add a member initializer
    //! @param NAME name for the data member
    //! @param DESCRIPTION description for the parameter
    //! @param INIT Serializer to push back
    void Serializer::AddInitializer
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const util::OwnPtr< SerializationInterface> &INIT
    )
    {
      BCL_Assert( INIT.IsDefined(), "Tried to add undefined initializer with name: " + NAME);
      AddMember( NAME, MemberInfo( INIT, e_Initializer, DESCRIPTION));
      ++m_NumberInits;
    }

    //! @brief add an optional member initializer
    //! @param NAME name for the data member
    //! @param DESCRIPTION description for the parameter
    //! @param INIT Serializer to push back
    void Serializer::AddOptionalInitializer
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const util::OwnPtr< SerializationInterface> &INIT
    )
    {
      BCL_Assert( INIT.IsDefined(), "Tried to add undefined initializer with name: " + NAME);
      MemberInfo mem( INIT, e_Initializer, DESCRIPTION);
      mem.SetAsOptional();
      AddMember( NAME, mem);
      ++m_NumberInits;
    }

    //! @brief add a member initializer with a default value, to be used if the parameter was omitted
    //! @param NAME name for the data member
    //! @param DESCRIPTION description for the parameter
    //! @param INIT member handler
    //! @param DEFAULT default value for the data member
    void Serializer::AddInitializer
    (
      const std::string &NAME,
      const std::string &DESCRIPTION,
      const util::OwnPtr< SerializationInterface> &INIT,
      const std::string &DEFAULT
    )
    {
      BCL_Assert( INIT.IsDefined(), "Tried to add undefined initializer with name: " + NAME);
      MemberInfo mem( INIT, e_Initializer, DESCRIPTION);
      mem.SetDefault( DEFAULT);
      AddMember( NAME, mem);
      ++m_NumberInits;
    }

    //! @brief add a data member
    //! @param NAME name for the data member
    //! @param DATA data member to push back
    void Serializer::AddDataMember
    (
      const std::string &NAME,
      const util::OwnPtr< SerializationInterface> &DATA
    )
    {
      BCL_Assert( DATA.IsDefined(), "Tried to add undefined initializer with name: " + NAME);
      AddMember( NAME, MemberInfo( DATA, e_Data, ""));
      ++m_NumberData;
    }

    //! @brief add a data member with a default value that will to be used if any other data members were given
    //! @param NAME name for the data member
    //! @param DATA data member handler
    //! @param DEFAULT default value for the data member
    void Serializer::AddDataMember
    (
      const std::string &NAME,
      const util::OwnPtr< SerializationInterface> &DATA,
      const std::string &DEFAULT
    )
    {
      BCL_Assert( DATA.IsDefined(), "Tried to add undefined initializer with name: " + NAME);
      MemberInfo mem( DATA, e_Data, "");
      mem.SetDefault( DEFAULT);
      ++m_NumberData;
      AddMember( NAME, mem);
    }

    //! @brief merge another serializer into this serializer. Useful only for statically-typed, non-container classes
    //! @param SERIALIZER the serializer for a statically-typed object; is consumed by this operation
    //! @note asserts if the serializer's type is not static object or if there was a name collision
    void Serializer::Merge( const Serializer &SERIALIZER)
    {
      BCL_Assert
      (
        SERIALIZER.GetType() == util::DataType::e_StaticObject
        || SERIALIZER.GetType() == util::DataType::s_NumberTypes,
        "Cannot merge a serializer for a " + util::DataType::GetTypeName( SERIALIZER.GetType())
        + " with another serializer; only static objects may be merged, because only they have multiple, known arguments. "
        " Scalars are merged automatically"
      );

      // iterate over initializers in the order they were added
      for
      (
        std::list< std::map< std::string, MemberInfo>::const_iterator>::const_iterator
          itr( SERIALIZER.m_MemberItrs.begin()), itr_end( SERIALIZER.m_MemberItrs.end());
        itr != itr_end;
        ++itr
      )
      {
        const std::string &name( ( *itr)->first);
        const MemberInfo &member( ( *itr)->second);
        BCL_Assert
        (
          m_HandlerMap.find( name) == m_HandlerMap.end(),
          "Cannot merge fields with name: \"" + name + "\""
        );

        BCL_Assert
        (
          member.GetState() != MemberInfo::e_WasSet,
          "Cannot merge serializer with a serializer that has already been processed!"
        );
        if( member.GetType() == e_Initializer)
        {
          if( member.GetState() == MemberInfo::e_Optional)
          {
            AddOptionalInitializer( name, member.GetDescription(), member.GetHandlerPtr());
          }
          else if( member.GetState() == MemberInfo::e_Default)
          {
            AddInitializer( name, member.GetDescription(), member.GetHandlerPtr(), member.GetDefault());
          }
          else if( member.GetState() == MemberInfo::e_Required)
          {
            AddInitializer( name, member.GetDescription(), member.GetHandlerPtr());
          }
        }
        else // if( member.GetType() == e_Data)
        {
          if( member.GetState() == MemberInfo::e_Default)
          {
            AddDataMember( name, member.GetHandlerPtr(), member.GetDefault());
          }
          else if( member.GetState() == MemberInfo::e_Required)
          {
            AddDataMember( name, member.GetHandlerPtr());
          }
        }
      }
    }

    //! @brief return the initialization label for this object
    //! @param WITH_DATA whether to include data members, otherwise, only initialization members should be included
    //! @return the initialization label for this object
    util::ObjectDataLabel Serializer::GetLabel( const bool &WITH_DATA) const
    {
      // handle special case of inlined tree
      if( m_InlineSelf)
      {
        util::ObjectDataLabel basic_label( m_HandlerMap.begin()->second.GetLabel( WITH_DATA));
        return basic_label;
      }
      else if( m_InlineArgs)
      {
        util::ObjectDataLabel basic_label
        (
          m_CmdLineIdentifier, m_HandlerMap.begin()->second.GetLabel( WITH_DATA).GetArguments()
        );
        return basic_label;
      }
      else if
      (
        m_CmdLineIdentifier == "" && m_NumberInits == size_t( 1) && m_NumberData == size_t( 0)
        && util::DataType::TestMustBeScalar( m_HandlerMap.begin()->second.GetHandler().GetSerializedType())
      )
      {
        util::ObjectDataLabel basic_label( m_HandlerMap.begin()->second.GetLabel( WITH_DATA));
        return basic_label;
      }

      std::vector< util::ObjectDataLabel> labels;
      labels.reserve( WITH_DATA ? m_NumberData + m_NumberInits : m_NumberInits);

      // walk over the members
      for
      (
        std::list< std::map< std::string, MemberInfo>::const_iterator>::const_iterator
          itr( m_MemberItrs.begin()), itr_end( m_MemberItrs.end());
        itr != itr_end;
        ++itr
      )
      {
        // get the map iterator
        const std::map< std::string, MemberInfo>::const_iterator &itr_map( *itr);
        // test whether the initializer is an initializer
        if( itr_map->second.GetType() == e_Initializer)
        {
          labels.push_back( itr_map->second.GetLabel( WITH_DATA));
          if( itr_map->second.GetState() != MemberInfo::e_Optional || !labels.back().IsEmpty())
          {
            labels.back().SetName( itr_map->first);
          }
          else
          {
            // remove the last label since it was empty and optional
            labels.pop_back();
          }
        }
        else if( WITH_DATA)
        {
          labels.push_back( itr_map->second.GetLabel( WITH_DATA));
          labels.back().SetName( itr_map->first);
        }
      }

      return util::ObjectDataLabel( m_CmdLineIdentifier, labels);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief function to get the keys from a map as a storage::Vector
    //! @param MAP the map from which to retrieve the keys
    //! @return get the keys from a map as a storage::Vector
    storage::Vector< std::string> GetKeys( const std::map< std::string, Serializer::MemberInfo> &MAP)
    {
      storage::Vector< std::string> keys;
      for
      (
        std::map< std::string, Serializer::MemberInfo>::const_iterator itr( MAP.begin()), itr_end( MAP.end());
        itr != itr_end;
        ++itr
      )
      {
        keys.PushBack( itr->first);
      }
      return keys;
    }

    //! @brief read from an object data label
    //! @param LABEL - the list of arguments
    //! @param ERROR_STREAM - the error stream to which errors should be written
    //! @return result of all internal validation
    ValidationResult Serializer::ReadArguments( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      BCL_MessageDbg
      (
        "Setting " + m_CmdLineIdentifier + " up using label " + LABEL.ToNamedString()
        + std::string( m_InlineArgs || m_InlineSelf ? " inline" : " ")
      );

      m_ReadData = false;

      if
      (
        LABEL.GetNumberArguments() == size_t( 1)
        && LABEL.GetArgument( 0).GetName().empty()
        && LABEL.GetArgument( 0).GetValue() == ValidationResult::GetHelpString()
        && LABEL.GetArgument( 0).IsScalar()
      )
      {
        // user asked for help
        WriteHelp( ERROR_STREAM, 0);
        return ValidationResult( e_Help);
      }

      ValidationResult result( true);

      // handle inlined serializers
      if( m_InlineSelf || m_InlineArgs)
      {
        // inlined parameter
        result = m_HandlerMap.begin()->second.Set( m_HandlerMap.begin()->first, LABEL, ERROR_STREAM);
      }
      else
      {
        // normal parameter
        util::ObjectDataLabel::const_iterator itr_label( LABEL.Begin()), itr_label_end( LABEL.End());

        // iterate over the named arguments; these are always at the end of the list due to sorting
        while( itr_label != itr_label_end)
        {
          // keep track of the result for setting this parameter
          ValidationResult result_this_label( e_Allowed);
          std::map< std::string, MemberInfo>::iterator initializer_itr( m_HandlerMap.end());
          if( !itr_label->GetName().empty())
          {
            initializer_itr = m_HandlerMap.find( itr_label->GetName());

            if( initializer_itr == m_HandlerMap.end())
            {
              ERROR_STREAM << '\n';
              command::Guesser::GetDefaultGuesser().WriteGuesses
              (
                itr_label->GetName(),
                GetKeys( m_HandlerMap),
                ERROR_STREAM,
                "parameter of " + m_CmdLineIdentifier
              );
            }
          }
          else
          {
            // look for a parameter with the same value
            initializer_itr = m_HandlerMap.find( itr_label->GetValue());
            if( initializer_itr == m_HandlerMap.end() && !itr_label->GetValue().empty())
            {
              // the label value is itself important
              initializer_itr = m_HandlerMap.find( "");
              if( initializer_itr == m_HandlerMap.end() && itr_label->GetValue() == "help")
              {
                // user asked for help
                WriteHelp( ERROR_STREAM, 0);
                result = ValidationResult( e_Help);
                ++itr_label;
                continue;
              }
            }
            if( initializer_itr == m_HandlerMap.end())
            {
              ERROR_STREAM << '\n';
              command::Guesser::GetDefaultGuesser().WriteGuesses
              (
                itr_label->GetValue(),
                GetKeys( m_HandlerMap),
                ERROR_STREAM,
                "parameter of " + m_CmdLineIdentifier
              );
            }
          }

          if( initializer_itr != m_HandlerMap.end() && initializer_itr->second.GetState() != MemberInfo::e_WasSet)
          {
            result_this_label = initializer_itr->second.Set( initializer_itr->first, *itr_label, ERROR_STREAM);

            // test whether data members of the class were read in
            if( initializer_itr->second.GetType() == e_Data)
            {
              m_ReadData = true;
            }
          }
          else
          {
            result_this_label = ValidationResult( e_Invalid);
            if( initializer_itr != m_HandlerMap.end() && initializer_itr->second.GetState() == MemberInfo::e_WasSet)
            {
              ERROR_STREAM << "\nFound duplicate key for " << initializer_itr->first
                           << " in object " << m_CmdLineIdentifier;
            }
          }
          if( result_this_label.IsInvalid())
          {
            if( !result.IsHelp())
            {
              result = result_this_label;
            }
          }
          else if( result_this_label.IsHelp())
          {
            result = result_this_label;
          }
          ++itr_label;
        }
      }

      // make sure that all parameters were set
      // unless help was requested
      if( !result.IsHelp())
      {
        for
        (
          std::map< std::string, MemberInfo>::iterator
            itr_param( m_HandlerMap.begin()), itr_param_end( m_HandlerMap.end());
          itr_param != itr_param_end;
          ++itr_param
        )
        {
          if( itr_param->second.GetState() == MemberInfo::e_WasSet)
          {
            // skip parameters that were set
            continue;
          }

          // skip data parameters unless at least one data parameter was given
          if( itr_param->second.GetType() == e_Data && !m_ReadData)
          {
            continue;
          }

          if( !itr_param->second.Finalize( itr_param->first, ERROR_STREAM))
          {
            result = ValidationResult( e_Invalid);
          }
        }
      }

      BCL_MessageDbg
      (
        m_CmdLineIdentifier + " given label: " + LABEL.ToString() + " set up label was: " + GetLabel( m_ReadData).ToString()
      );
      if( result.IsInvalid())
      {
        ERROR_STREAM << "\n[ in label: " << LABEL.ToString() << ']';
      }
      return result;
    }

    //! @brief set a particular initializer's default. Useful in subclasses for which certain parameters are unnecessary
    //! @param NAME name of the parameter
    //! @param VALUE value to set the default parameter with
    void Serializer::SetDefault( const std::string &NAME, const std::string &VALUE)
    {
      util::ObjectDataLabel label( VALUE);
      std::map< std::string, MemberInfo>::iterator initializer_itr( m_HandlerMap.end());
      // handle inlined serializers
      if( m_InlineSelf || m_InlineArgs)
      {
        // inlined parameter
        initializer_itr = m_HandlerMap.begin();
      }
      else
      {
        // keep track of the result for setting this parameter
        initializer_itr = m_HandlerMap.find( NAME);
      }
      BCL_Assert( initializer_itr != m_HandlerMap.end(), "Could not find parameter with name " + NAME + " to fix");
      initializer_itr->second.SetDefault( label);
    }

    //! @brief writes the help for the Serializer
    //! @param OSTREAM - the outstream to write to
    //! @param INDENT - the indent to use
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &Serializer::WriteHelp( std::ostream &OSTREAM, const size_t &INDENT) const
    {
      FixedLineWidthWriter writer;
      writer.SetAutoIndentExtra( 2);
      writer.SetBclIndent( INDENT);
      this->WriteHelp( writer);
      OSTREAM << writer.String().substr( writer.GetIndent());
      return OSTREAM;
    }

    //! @brief writes the help for the Serializer
    //! @param OSTREAM - the outstream to write to
    //! @return return the stream after writing
    FixedLineWidthWriter &Serializer::WriteHelp( FixedLineWidthWriter &OSTREAM) const
    {
      if( !m_NumberInits)
      {
        if( m_ClassDescription.size())
        {
          OSTREAM << m_ClassDescription;
        }
        OSTREAM.NewLine();
      }
      else
      {
        bool added_indent( false);
        if( m_ClassDescription.size())
        {
          OSTREAM << m_ClassDescription << '\n';
          OSTREAM.AddIndent( 2);
          added_indent = true;
        }
        if( !m_InlineSelf)
        {
          OSTREAM << "Default label : " << this->GetLabel().ToString() << '\n';
          if( !added_indent)
          {
            OSTREAM.AddIndent( 2);
            added_indent = true;
          }
        }
        if( m_InlineArgs)
        {
          if( !m_HandlerMap.begin()->first.empty())
          {
            // write name
            OSTREAM << "Parameter: <" << m_HandlerMap.begin()->first << "> ";
          }
          else
          {
            OSTREAM << "(anonymous) parameter: ";
          }
          m_HandlerMap.begin()->second.WriteHelp( OSTREAM) << '\n';
        }
        else
        {
          OSTREAM << "Parameters:\n";
          // write help for the sub initializers
          for
          (
            std::list< std::map< std::string, MemberInfo>::const_iterator>::const_iterator
              itr( m_MemberItrs.begin()), itr_end( m_MemberItrs.end());
            itr != itr_end;
            ++itr
          )
          {
            const std::map< std::string, MemberInfo>::const_iterator &itr_map( *itr);
            if( itr_map->second.GetType() == e_Initializer)
            {
              if( !itr_map->first.empty())
              {
                // write name
                OSTREAM << "<" << itr_map->first << "> ";
              }
              itr_map->second.WriteHelp( OSTREAM) << '\n';
            }
          }
        }
        if( added_indent)
        {
          OSTREAM.PopIndent();
        }
      }
      return OSTREAM;
    }

    //! @brief writes the brief help for the Serializer
    //! @param OSTREAM - the outstream to write to
    //! @return return the stream after writing
    FixedLineWidthWriter &Serializer::WriteBriefHelp( FixedLineWidthWriter &OSTREAM) const
    {
      OSTREAM << "* " << m_CmdLineIdentifier << " : " << m_ClassDescription;
      if( m_NumberInits)
      {
        OSTREAM << "; ";
        OSTREAM.WriteOnOneLine( "\"" + m_CmdLineIdentifier + "(help)\" shows internal options");
      }
      return OSTREAM;
    }

    //! @brief read Serializer for std::istream ISTREAM
    //! @param ISTREAM - the instream to read from
    //! @return std::istream &ISTREAM - return the stream after reading
    std::istream &Serializer::Read( std::istream &ISTREAM)
    {
      std::ostringstream err_stream;
      // assert that there were no errors
      BCL_Assert
      (
        ReadArguments( util::ObjectDataLabel( ISTREAM), err_stream) && err_stream.tellp() == std::streampos( 0),
        "Error while reading " + m_CmdLineIdentifier + ": " + err_stream.str()
      );
      return ISTREAM;
    }

    //! @brief writes Serializer to ostream OSTREAM
    //! @param OSTREAM - the outstream to write to
    //! @param INDENT indentation
    //! @return std::ostream &OSTREAM - return the stream after writing
    std::ostream &Serializer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // Write out the object with the given indent
      OSTREAM << GetLabel( true).ToNamedStringDefaultWidth( INDENT * Serialize::s_Number_spaces_per_indent);
      return OSTREAM;
    }

    //! @brief Get addresses of all objects serialized as part of this
    //! @notes used to ensure uniqueness of serialized objects
    std::vector< size_t> Serializer::GetSerializationAddresses() const
    {
      std::vector< size_t> addresses;
      for
      (
        std::map< std::string, MemberInfo>::const_iterator
          itr_param( m_HandlerMap.begin()), itr_param_end( m_HandlerMap.end());
        itr_param != itr_param_end;
        ++itr_param
      )
      {
        std::vector< size_t> to_append( itr_param->second.GetHandler().GetSerializationAddresses());
        addresses.insert( addresses.end(), to_append.begin(), to_append.end());
      }
      return addresses;
    }

    //! @brief add an additional data or member initializer
    //! @param NAME name of the member initializer
    //! @param MEMBER_INFO information about tha
    //! @param DESCRIPTION description for the parameter
    //! @param HANDLER object to handle parsing / setting this member
    //! @param TYPE the type of object
    void Serializer::AddMember
    (
      const std::string &NAME,
      const MemberInfo &MEMBER_INFO
    )
    {
      // ensure that the serializer has not already been finalized
      BCL_Assert
      (
        !IsFinalized(),
        "Tried to add member with name: " + NAME
        + " to a serializer after calling SetTypeAndFinalize"
      );

      // try to insert the member with the given name into the map
      std::pair< std::map< std::string, MemberInfo>::const_iterator, bool>
        itr_success( m_HandlerMap.insert( std::make_pair( NAME, MEMBER_INFO)));

      // test that a member with that name was not already present
      BCL_Assert
      (
        itr_success.second,
        m_CmdLineIdentifier + " Serializer tried to add parameter with duplicate name: " + NAME
      );

      if( command::CommandState::IsInStaticInitialization())
      {
        std::vector< size_t> new_addresses( MEMBER_INFO.GetHandler().GetSerializationAddresses());
        const std::vector< size_t>::const_iterator itr_new( new_addresses.begin()), itr_new_end( new_addresses.end());
        // test to see whether the new addresses overlapped the old
        for
        (
          std::map< std::string, MemberInfo>::const_iterator
            itr( m_HandlerMap.begin()), itr_end( m_HandlerMap.end());
          itr != itr_end;
          ++itr
        )
        {
          if( itr == itr_success.first)
          {
            continue;
          }
          std::vector< size_t> existing_addresses( itr->second.GetHandler().GetSerializationAddresses());
          for
          (
            std::vector< size_t>::const_iterator
              itr_e( existing_addresses.begin()), itr_e_end( existing_addresses.end());
            itr_e != itr_e_end;
            ++itr_e
          )
          {
            BCL_Assert
            (
              std::find( itr_new, itr_new_end, *itr_e) == itr_new_end,
              itr->first + " serializes at least one of the same things as newly added member: " + NAME
              + ". This is most likely due to io::Serialization::GetAgent() being called on same class member with "
              " different parameter name for object " + m_CmdLineIdentifier
            );
          }
        }
      }

      // add the iterator to the ordered list for later use
      m_MemberItrs.push_back( itr_success.first);
    }

  } // namespace io
} // namespace bcl
