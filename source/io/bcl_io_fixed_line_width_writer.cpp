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
#include "io/bcl_io_fixed_line_width_writer.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_guesser.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace io
  {

    //! @brief default constructor, uses logger's width - 1
    FixedLineWidthWriter::FixedLineWidthWriter() :
      m_CurrentLinePosition( 0),
      m_DesiredLineWidth( util::GetLogger().GetMaxLineWidth() - 1),
      m_AutoExtraIndent( 0)
    {
    }

    //! @brief default constructor
    //! @param AUTO_INDENT_EXTRA # spaces extra to indent continued lines
    //! @param DESIRED_LINE_WIDTH desired width of lines
    FixedLineWidthWriter::FixedLineWidthWriter( const size_t &AUTO_INDENT_EXTRA, const size_t &DESIRED_LINE_WIDTH) :
      m_CurrentLinePosition( 0),
      m_DesiredLineWidth( DESIRED_LINE_WIDTH),
      m_AutoExtraIndent( AUTO_INDENT_EXTRA)
    {
    }

    //! @brief virtual desctructor
    FixedLineWidthWriter::~FixedLineWidthWriter()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &FixedLineWidthWriter::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the current indent level in actual spaces
    size_t FixedLineWidthWriter::GetIndent() const
    {
      return m_CurrentIndent.size();
    }

    //! @brief return the amount each newly-created line will be indented
    size_t FixedLineWidthWriter::GetExtraIndent() const
    {
      return m_AutoExtraIndent;
    }

    //! @brief get the position on the current line
    //! @return the position on the current line
    size_t FixedLineWidthWriter::GetLinePosition()
    {
      return m_CurrentLinePosition + m_CurrentLine.tellp();
    }

    //! @brief get the amount of space remaining on the current line
    //! @return the amount of space remaining on the current line
    size_t FixedLineWidthWriter::GetRemainingSpaceOnLine()
    {
      const size_t pos( GetLinePosition());
      if( !pos)
      {
        return m_DesiredLineWidth - m_CurrentIndent.size();
      }
      return std::max( pos, m_DesiredLineWidth) - pos;
    }

    //! @brief get the maximum # of non-space characters on a given line currently
    //! @return the maximum # of non-space characters on a given line currently
    size_t FixedLineWidthWriter::GetEffectiveLineWidth()
    {
      return std::max( m_DesiredLineWidth, m_CurrentIndent.size()) - m_CurrentIndent.size();
    }

    //! @brief set autoindent extra to a different value
    //! @param AUTO_INDENT_EXTRA # spaces extra to indent continued lines
    void FixedLineWidthWriter::SetAutoIndentExtra( const size_t &AUTO_INDENT_EXTRA)
    {
      m_AutoExtraIndent = AUTO_INDENT_EXTRA;
    }

    //! @brief add a certain amount to the indent
    //! @param INDENT amount extra to indent
    void FixedLineWidthWriter::AddIndent( const size_t &INDENT)
    {
      if( m_IndentLevels.empty())
      {
        m_IndentLevels.push_back( INDENT);
      }
      else
      {
        m_IndentLevels.push_back( m_IndentLevels.back() + INDENT);
      }
      if( INDENT)
      {
        UpdateIndentString();
      }
    }

    //! @brief pop the last indent level off the stack
    void FixedLineWidthWriter::PopIndent()
    {
      BCL_Assert( !m_IndentLevels.empty(), "Tried to pop indent, but no indentation exists");
      m_IndentLevels.pop_back();
      UpdateIndentString();
    }

    //! @brief specifically set the indent, without changing previous indent levels
    //! @param INDENT the new indentation level
    void FixedLineWidthWriter::SetIndent( const size_t &INDENT)
    {
      m_IndentLevels.push_back( INDENT);
      UpdateIndentString();
    }

    //! @brief specifically set the indent in terms of bcl-indent levels (io::Serialize::s_Number_spaces_per_indent)
    //! @param INDENT the new indentation level
    void FixedLineWidthWriter::SetBclIndent( const size_t &INDENT)
    {
      SetIndent( Serialize::s_Number_spaces_per_indent * INDENT);
    }

    //! @brief clear all indentation
    void FixedLineWidthWriter::ClearIndent( const size_t &INDENT)
    {
      m_IndentLevels.clear();
      UpdateIndentString();
    }

    //! @brief write a given value on a single line (the current line, if it will fit, the next if it will not)
    //! @param VALUE the string to write
    void FixedLineWidthWriter::WriteOnOneLine( const std::string &VALUE)
    {
      if( GetLinePosition() == 0)
      {
        m_CurrentLine << m_CurrentIndent;
      }

      const size_t end_line_pos( std::min( VALUE.find( '\n'), VALUE.size()));
      if( end_line_pos <= m_DesiredLineWidth - GetLinePosition())
      {
        m_CurrentLine << VALUE;
      }
      else
      {
        NewLineIndent();
        m_CurrentLine << std::string( m_AutoExtraIndent, ' ');
        m_CurrentLine << VALUE;
        FlushLineAutoNewline();
      }
    }

    //! @brief Write a heading centered on one or more lines, optionally with a bordering line (top and bottom)
    //! @param HEADING the heading to write centered.  May be divided up into more than one line, if necessary.
    //!        HEADING should not contain any new lines though
    //! @param FILLER the character to use to fill the space around the heading. If not ' ', one ' ' will be
    //!        automatically added around the heading
    //! @param BORDER if true, write one line before and after the heading, completely full of FILLER
    void FixedLineWidthWriter::WriteHeading( const std::string &HEADING, const char &FILLER, const bool &BORDER)
    {
      if( !GetLinePosition())
      {
        m_CurrentLine << m_CurrentIndent;
      }
      else if( GetLinePosition() > GetIndent())
      {
        NewLineIndent();
      }
      // effective line width
      const size_t eff_line_width( GetEffectiveLineWidth());
      if( BORDER)
      {
        m_CurrentLine << std::string( eff_line_width, FILLER);
        NewLineIndent();
      }
      if( HEADING.size())
      {
        const size_t eff_heading_size( FILLER == ' ' ? HEADING.size() + 2 : HEADING.size() + 4);
        size_t n_lines( ( eff_heading_size - 1) / eff_line_width + 1);
        const size_t n_filler( ( n_lines * eff_line_width - eff_heading_size) / 2);
        operator <<( std::string( n_filler, FILLER));
        if( FILLER != ' ')
        {
          operator <<( ' ');
        }
        operator <<( HEADING);
        if( FILLER != ' ')
        {
          operator <<( ' ');
        }
        m_CurrentLine << std::string( GetRemainingSpaceOnLine(), FILLER);
        NewLine();
      }
      if( BORDER)
      {
        m_CurrentLine << m_CurrentIndent << std::string( eff_line_width, FILLER);
        NewLine();
      }
    }

    //! @brief shortcut for operator << '\n'
    void FixedLineWidthWriter::NewLine()
    {
      const std::string strn( m_CurrentLine.str());
      m_Complete.write( strn.c_str(), strn.size());
      m_Complete << '\n';
      m_CurrentLinePosition = 0;
      m_CurrentLine.str( std::string());
    }

    //! @brief newline, followed by indent
    void FixedLineWidthWriter::NewLineIndent()
    {
      NewLine();
      m_CurrentLine << m_CurrentIndent;
    }

    //! @brief get the current string
    std::string FixedLineWidthWriter::String()
    {
      FlushLine();
      return m_Complete.str();
    }

    //! @brief flush the current line to the complete string stream and write a newline
    void FixedLineWidthWriter::Endl()
    {
      NewLine();
      FlushLine();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief overload specifically for strings, which should be indented properly
    //! @param VALUE value to write
    FixedLineWidthWriter &FixedLineWidthWriter::operator <<( const std::string &VALUE)
    {
      if( VALUE.empty())
      {
        // just return
        return *this;
      }
      const size_t new_line_pos( VALUE.rfind( '\n'));
      if( new_line_pos != std::string::npos)
      {
        // initialize results
        storage::Vector< std::string> output_strings;

        // find all \n; do not use SplitString, as it would remove blank lines
        {
          std::string::size_type last_pos( 0);
          std::string::size_type pos( VALUE.find( '\n'));

          // as long as something was found
          while( pos != std::string::npos)
          {
            // Found a token, add it to the vector.
            output_strings.PushBack( VALUE.substr( last_pos, pos - last_pos));
            last_pos = pos + 1;
            // Find next "non-delimiter"
            pos = VALUE.find( '\n', last_pos);
          }
          if( last_pos < VALUE.size())
          {
            output_strings.PushBack( VALUE.substr( last_pos));
          }
        }

        this->operator <<( output_strings.FirstElement());
        for
        (
          storage::Vector< std::string>::const_iterator
            itr( output_strings.Begin() + 1), itr_end( output_strings.End());
          itr != itr_end;
          ++itr
        )
        {
          NewLine();
          this->operator <<( *itr);
        }
        if( VALUE[ VALUE.size() - 1] == '\n')
        {
          NewLine();
        }
        return *this;
      }

      if( GetLinePosition() == 0)
      {
        m_CurrentLine << m_CurrentIndent;
      }

      // main case, VALUE contains no newlines
      std::ostringstream output;
      output << VALUE;
      const std::string output_str( output.str());

      size_t total_line_length( GetLinePosition());

      if( total_line_length + output_str.size() <= m_DesiredLineWidth)
      {
        m_CurrentLine << output_str;
        return *this;
      }
      else
      {
        // handle excessive indentation by just writing out whatever was desired, no word wrap
        if( m_CurrentIndent.size() >= m_DesiredLineWidth / 2)
        {
          m_CurrentLine << output_str;
          NewLine();
          return *this;
        }

        size_t additional_indent( output_str.find_first_not_of( " \t"));
        if( additional_indent == std::string::npos)
        {
          additional_indent = 0;
        }
        else
        {
          additional_indent += GetLinePosition();
        }

        // line that requires some wrapping
        size_t previous_position( 0);

        // handle the first line separately, as it is subject to START_POSITION, unlike the others
        {
          // line length limit accounting for indent and first line size
          const size_t start_limit( m_DesiredLineWidth - GetLinePosition());

          // find a space that is ideally at LINE_LENGTH_LIMIT away from the current position (accounting for indent)
          // or the size of the current string, whichever comes first
          size_t position( output_str.rfind( ' ', start_limit));

          // handle the case where no spaces are available in output_str before start_limit
          if( position == std::string::npos)
          {
            size_t first_space_position( output_str.find( ' ', start_limit));
            if( first_space_position + m_CurrentIndent.size() < m_DesiredLineWidth)
            {
              NewLine();
              m_CurrentLine << m_CurrentIndent;
              position = output_str.rfind( ' ', m_DesiredLineWidth - GetLinePosition());
            }
            else
            {
              // no spaces for more than m_DesiredLineWidth - current indent (e.g. a very very long word), just write it
              // on this line
              position = first_space_position;
            }
          }
          m_CurrentLine << output_str.substr( 0, position);

          // no spaces found
          if( position == std::string::npos)
          {
            return *this;
          }
          previous_position = output_str.find_first_not_of( ' ', position);
          if( previous_position == std::string::npos)
          {
            return *this;
          }
          NewLine();
        }

        // update the indent, if necessary
        AddIndent( m_AutoExtraIndent);

        // line length limit accounting for indent
        // indent continued lines more if there was more than one block
        const size_t effective_line_length_limit( m_DesiredLineWidth - m_CurrentIndent.size());

        while( previous_position + effective_line_length_limit < output_str.size())
        {
          // find a space that is ideally at LINE_LENGTH_LIMIT away from the current position (accounting for indent)
          // or the size of the current string, whichever comes first
          const size_t next_line_ideal_end( previous_position + effective_line_length_limit);
          size_t position( output_str.rfind( ' ', next_line_ideal_end));
          if( position == std::string::npos || position <= previous_position)
          {
            position = output_str.find( ' ', next_line_ideal_end);
          }
          m_CurrentLine << m_CurrentIndent;
          m_CurrentLine << output_str.substr( previous_position, position - previous_position);
          previous_position = position;
          if( position == std::string::npos)
          {
            PopIndent();
            return *this;
          }
          previous_position = output_str.find_first_not_of( ' ', position);
          if( previous_position == std::string::npos)
          {
            PopIndent();
            return *this;
          }
          NewLine();
        }
        if( previous_position < output_str.size())
        {
          m_CurrentLine << m_CurrentIndent;
          m_CurrentLine << output_str.substr( previous_position);
        }
        FlushLineAutoNewline();
        PopIndent();
      }
      return *this;
    }

    //! @brief overload specifically for characters
    //! @param VALUE value to write
    FixedLineWidthWriter &FixedLineWidthWriter::operator <<( const char &VALUE)
    {
      if( VALUE == '\n')
      {
        NewLine();
      }
      else if( VALUE != ' ' || GetLinePosition() > 0)
      {
        std::ostringstream output;
        if( GetLinePosition() == 0)
        {
          m_CurrentLine << m_CurrentIndent;
        }
        output << VALUE;
        if( m_CurrentLinePosition + m_CurrentIndent.size() + output.tellp() >= m_DesiredLineWidth)
        {
          NewLine();
        }
        m_CurrentLine << output.str();
      }
      return *this;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief flush the current line to the complete string stream
    void FixedLineWidthWriter::FlushLine()
    {
      m_CurrentLinePosition = m_CurrentLine.tellp();
      if( m_CurrentLinePosition)
      {
        const std::string strn( m_CurrentLine.str());
        m_Complete.write( strn.c_str(), strn.size());
        m_CurrentLine.str( std::string());
      }
    }

    //! @brief flush the current line to the complete string stream; carefully evaluate where the previous new line was; insert a
    //! newline automatically if the last new line was too long ago
    //! This is necessary if it is possible that an output string contained a new line
    void FixedLineWidthWriter::FlushLineAutoNewline()
    {
      if( m_CurrentLine.tellp())
      {
        const std::string output( m_CurrentLine.str());
        size_t last_newline( std::min( output.rfind( '\n'), output.size() - size_t( 1)));
        m_CurrentLinePosition = output.size() - last_newline;
        if( m_CurrentLinePosition >= m_DesiredLineWidth)
        {
          NewLine();
        }
        else
        {
          m_Complete.write( output.c_str(), output.size());
          m_CurrentLine.str( std::string());
        }
      }
    }

    //! @brief update the current indentation string
    void FixedLineWidthWriter::UpdateIndentString()
    {
      if( m_IndentLevels.empty())
      {
        m_CurrentIndent.erase();
      }
      else if( m_CurrentIndent.size() != m_IndentLevels.back())
      {
        m_CurrentIndent = std::string( m_IndentLevels.back(), ' ');
      }
    }

  } // namespace io
} // namespace bcl
