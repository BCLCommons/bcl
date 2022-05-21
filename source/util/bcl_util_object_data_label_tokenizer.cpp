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
#include "util/bcl_util_object_data_label_tokenizer.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_object_instances.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const SiPtr< const ObjectInterface> ObjectDataLabelTokenizer::s_Instance
    (
      GetObjectInstances().AddInstance( new ObjectDataLabelTokenizer( ""))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from string
    ObjectDataLabelTokenizer::ObjectDataLabelTokenizer( const std::string &STR, const size_t &POSITION) :
      m_String( STR),
      m_Size( STR.size()),
      m_Position( POSITION),
      m_ScopeDepth( 0),
      m_NextType( e_Start)
    {
      DetermineNextTokenType();
    }

    //! @brief Clone function
    //! @return pointer to new ObjectDataLabelTokenizer
    ObjectDataLabelTokenizer *ObjectDataLabelTokenizer::Clone() const
    {
      return new ObjectDataLabelTokenizer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ObjectDataLabelTokenizer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get token type
    //! @return the current token type
    const ObjectDataLabelTokenizer::Type &ObjectDataLabelTokenizer::GetLastTokenType() const
    {
      return m_LastType;
    }

    //! @brief get type of token that will next be returned by Pop
    //! @return the type of token that will next be returned by Pop
    const ObjectDataLabelTokenizer::Type &ObjectDataLabelTokenizer::GetNextTokenType() const
    {
      return m_NextType;
    }

    //! @brief get token type
    //! @return the current token type
    const size_t &ObjectDataLabelTokenizer::GetScopeDepth() const
    {
      return m_ScopeDepth;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get the string for the next token
    std::string ObjectDataLabelTokenizer::Pop()
    {
      m_LastType = m_NextType;
      if( m_LastType == e_End)
      {
        return std::string();
      }

      // handle scalars
      if( m_LastType == e_Scalar)
      {
        size_t start_position( m_Position);
        const char first_char( m_String[ m_Position]);

        // handle quoted strings
        if( first_char == '"')
        {
          ++start_position;
          m_Position = ObjectDataLabelTokenizer::FindEndQuote( m_String, m_Position);
          if( m_Position == std::string::npos)
          {
            // no end quote, return complete string
            return m_String.substr( start_position);
          }

          const size_t end_position( m_Position++);

          // determine next type
          DetermineNextTokenType();

          return m_String.substr( start_position, end_position - start_position);
        }

        // handle normal scalars, keep track of the last position that was not a space

        // find the next delimiter
        m_Position = m_String.find_first_of( ObjectDataLabel::GetAllDelimiters(), m_Position);

        // if there are no remaining delimiters, reset to the end of the string
        if( m_Position == std::string::npos)
        {
          m_Position = m_Size;
        }

        // determine the next token type
        DetermineNextTokenType();

        // find the last non-space character, since the scalar does not include leading/trailing spaces
        size_t last_non_space_pos( m_Position - 1);
        while( isspace( m_String[ last_non_space_pos]))
        {
          --last_non_space_pos;
        }

        // return the substring containing only the scalar
        return m_String.substr( start_position, last_non_space_pos + 1 - start_position);
      }

      // handle delimiters
      // get the current character and advance the position
      const char character( m_String[ m_Position++]);

      // determine the next token type
      DetermineNextTokenType();

      // handle scope depth changes
      if( m_LastType == e_ScopeOpen)
      {
        ++m_ScopeDepth;
      }
      else if( m_LastType == e_ScopeClose)
      {
        --m_ScopeDepth;
      }

      // return the string with the character
      return std::string( 1, character);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief find the end of the quoted string that starts at a given position, ignoring escaped quotes
    //! @param STR string of interest
    //! @param START_QUOTE_POS the position of the starting quote
    //! @return the end quote position, or std::string::npos if no such quote exists
    size_t ObjectDataLabelTokenizer::FindEndQuote( const std::string &STR, const size_t &START_QUOTE_POS)
    {
      const size_t last_char_pos( STR.size() - 1);
      const char quote( STR[ START_QUOTE_POS]);

      // check for quotes at the end of the string
      if( START_QUOTE_POS == last_char_pos)
      {
        return std::string::npos;
      }

      // find the next quote character
      size_t pos( STR.find( quote, START_QUOTE_POS + 1));

      // continue while the quote character is preceded by an escape character
      while( pos < last_char_pos && STR[ pos - 1] == '\\')
      {
        pos = STR.find( quote, pos + 1);
      }

      // if pos is exactly at the last character, make sure that it was not preceded by an escape
      if( pos == last_char_pos && STR[ pos - 1] == '\\')
      {
        return std::string::npos;
      }

      // return position
      return pos;
    }

    //! @brief  validate a string, write out errors to a given stream
    //! @param  STR the string to validate
    //! @param  ERR the stream to write errors to
    //! @return true if the string was valid
    bool ObjectDataLabelTokenizer::Validate( const std::string &STR, std::ostream &ERR)
    {
      ObjectDataLabelTokenizer validator( STR);
      Type previous_type( e_Start);
      while( validator.GetNextTokenType() != e_End)
      {
        // check for invalid combinations of last/next tokens
        if( validator.GetScopeDepth() == 0)
        {
          if( validator.GetNextTokenType() == e_ArgDelimiter)
          {
            ERR << "Arguments are not allowed without an associated object, that is, ',' cannot appear outside\n"
                << "parenthetical scope for object data labels, but did in given string:\n"
                << STR << std::endl;
            return false;
          }
          else if( validator.GetNextTokenType() == e_ScopeClose)
          {
            ERR << "Excessive end-scope parenthesis in given string:\n"
                << STR << std::endl;
            return false;
          }
        }
        if( validator.GetLastTokenType() == e_Scalar && validator.GetNextTokenType() == e_Scalar)
        {
          ERR << "Variable next to quoted string is forbidden, but was found in given string:\n" << STR << std::endl;
          return false;
        }
        if( validator.GetLastTokenType() == e_ScopeClose)
        {
          if
          (
            validator.GetNextTokenType() != e_ArgDelimiter
            && validator.GetNextTokenType() != e_ScopeClose
            && validator.GetNextTokenType() != e_End
          )
          {
            ERR << "Missing comma near " << validator.Pop() << " in given string:\n"
                << STR << std::endl;
            return false;
          }
        }

        // check that the string =X= never appears
        // people with programming might try to use this to set multiple variables equal to some value,
        // but this currently is not supported
        if
        (
          previous_type == e_TagDelimiter
          && validator.GetLastTokenType() == e_Scalar
          && validator.GetNextTokenType() == e_TagDelimiter
        )
        {
          ERR << "Multiple parameter assignment is forbidden but was found in given string:\n"
              << STR << std::endl;
          return false;
        }

        previous_type = validator.GetLastTokenType();
        validator.Pop();
      }
      if( validator.GetScopeDepth())
      {
        ERR << validator.GetScopeDepth() << " unclosed parenthesis in given string:\n"
            << STR << std::endl;
        return false;
      }
      return true;
    }

    //! @brief update the next token type, after skipping spaces
    void ObjectDataLabelTokenizer::DetermineNextTokenType()
    {
      while( m_Position < m_Size && isspace( m_String[ m_Position]))
      {
        ++m_Position;
      }

      m_LastType = m_NextType;
      if( m_Position >= m_Size)
      {
        m_NextType = e_End;
      }
      else
      {
        // switch to find delimiter-based enums
        bool is_file( false);
        switch( m_String[ m_Position])
        {
          case '(':
            m_NextType = e_ScopeOpen;
            break;
          case ')':
            m_NextType = e_ScopeClose;
            break;
          case ',':
            m_NextType = e_ArgDelimiter;
            break;
          case '=':
            m_NextType = e_TagDelimiter;
            break;
          case '@':
            is_file = true;
            break;
          default:
            // all other letters indicate scalars
            m_NextType = e_Scalar;
            break;
        };
        if( is_file)
        {
          const size_t end_filename_pos( m_String.find_first_of( ObjectDataLabel::GetAllDelimiters(), m_Position + 1));
          const std::string filename( RStrip( m_String.substr( m_Position + 1, end_filename_pos - m_Position - 1), " \t\n\r"));

          io::IFStream ist;
          io::File::MustOpenIFStream( ist, filename);
          const std::string complete_file
          (
            Join( "\n", util::StringListFromIStream( ist))
          );
          m_String = m_String.substr( 0, m_Position) + complete_file + m_String.substr( end_filename_pos);
          m_Size = m_String.size();
          DetermineNextTokenType();
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ObjectDataLabelTokenizer::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_String, ISTREAM);
      m_Size = m_String.size();

      size_t new_position( 0);
      io::Serialize::Read( new_position, ISTREAM);

      // get back to the previous state
      m_Position = 0;
      m_LastType = m_NextType = e_Start;
      while( m_Position < new_position)
      {
        Pop();
      }

      DetermineNextTokenType();
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ObjectDataLabelTokenizer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_String, OSTREAM, INDENT) << '\t' << m_Position;
      return OSTREAM;
    }

  } // namespace util
} // namespace bcl
