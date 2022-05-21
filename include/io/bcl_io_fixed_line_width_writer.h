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

#ifndef BCL_IO_FIXED_LINE_WIDTH_WRITER_H_
#define BCL_IO_FIXED_LINE_WIDTH_WRITER_H_

// include the namespace header
#include "bcl_io.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically
#include <iostream>
#include <list>
#include <sstream>

namespace bcl
{
  namespace io
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FixedLineWidthWriter
    //! @brief String stream that handles new lines and indentation to maintain max line width and desired indentation
    //!
    //! @see @link example_io_fixed_line_width_writer.cpp @endlink
    //! @author mendenjl
    //! @date Feb 28, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FixedLineWidthWriter
    {

    //////////
    // data //
    //////////

      std::ostringstream  m_Complete;            //!< Stream containing all completed lines
      std::ostringstream  m_CurrentLine;         //!< Stream for current line
      size_t              m_CurrentLinePosition; //!< Position on the current line
      std::list< size_t>  m_IndentLevels;        //!< Indentation levels current used
      std::string         m_CurrentIndent;       //!< Current indentation string
      size_t              m_DesiredLineWidth;    //!< Desired line width
      size_t              m_AutoExtraIndent;     //!< # spaces extra to indent continued lines

    public:

      //! @brief default constructor, uses logger's width - 1
      FixedLineWidthWriter();

      //! @brief constructor from desired line width and amount extra to indent continued lines
      //! @param AUTO_INDENT_EXTRA # spaces extra to indent continued lines
      //! @param DESIRED_LINE_WIDTH desired width of lines
      FixedLineWidthWriter
      (
        const size_t &AUTO_INDENT_EXTRA,
        const size_t &DESIRED_LINE_WIDTH
      );

      //! @brief virtual desctructor
      virtual ~FixedLineWidthWriter();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the current indent level in actual spaces
      size_t GetIndent() const;

      //! @brief return the amount each newly-created line will be indented
      size_t GetExtraIndent() const;

      //! @brief get the position on the current line
      //! @return the position on the current line
      size_t GetLinePosition();

      //! @brief get the amount of space remaining on the current line
      //! @return the amount of space remaining on the current line
      size_t GetRemainingSpaceOnLine();

      //! @brief get the maximum # of non-space characters on a given line currently
      //! @return the maximum # of non-space characters on a given line currently
      size_t GetEffectiveLineWidth();

      //! @brief set autoindent extra to a different value
      //! @param AUTO_INDENT_EXTRA # spaces extra to indent continued lines
      void SetAutoIndentExtra( const size_t &AUTO_INDENT_EXTRA);

    ////////////////
    // operations //
    ////////////////

      //! @brief add a certain amount to the indent
      //! @param INDENT amount extra to indent
      void AddIndent( const size_t &INDENT);

      //! @brief pop the last indent level off the stack
      void PopIndent();

      //! @brief specifically set the indent, without changing previous indent levels
      //! @param INDENT the new indentation level
      void SetIndent( const size_t &INDENT);

      //! @brief specifically set the indent in terms of bcl-indent levels (io::Serialize::s_Number_spaces_per_indent)
      //! @param INDENT the new indentation level
      void SetBclIndent( const size_t &INDENT);

      //! @brief clear all indentation
      void ClearIndent( const size_t &INDENT);

      //! @brief write a given value on a single line (the current line, if it will fit, the next if it will not)
      //! @param VALUE the string to write
      void WriteOnOneLine( const std::string &VALUE);

      //! @brief Write a heading centered on one or more lines, optionally with a bordering line (top and bottom)
      //! @param HEADING the heading to write centered.  May be divided up into more than one line, if necessary.
      //!        HEADING should not contain any new lines though
      //! @param FILLER the character to use to fill the space around the heading. If not ' ', one ' ' will be
      //!        automatically added around the heading
      //! @param BORDER if true, write one line before and after the heading, completely full of FILLER
      void WriteHeading( const std::string &HEADING, const char &FILLER, const bool &BORDER);

      //! @brief shortcut for operator << '\n'
      void NewLine();

      //! @brief newline, followed by indent
      void NewLineIndent();

      //! @brief get the current string
      std::string String();

      //! @brief flush the current line to the complete string stream and write a newline
      void Endl();

    ///////////////
    // operators //
    ///////////////

      //! @brief stream insertion operator
      //! @param VALUE value to write
      template< class t_DataType>
      FixedLineWidthWriter &operator <<( const t_DataType &VALUE)
      {
        std::ostringstream output;
        output << VALUE;
        return *this << output.str();
      }

      //! @brief overload specifically for strings, which should be indented properly
      //! @param VALUE value to write
      FixedLineWidthWriter &operator <<( const std::string &VALUE);

      //! @brief overload specifically for characters
      //! @param VALUE value to write
      FixedLineWidthWriter &operator <<( const char &VALUE);

    //////////////////////
    // input and output //
    //////////////////////

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief flush the current line to the complete string stream
      void FlushLine();

      //! @brief flush the current line to the complete string stream; carefully evaluate where the previous new line was; insert a
      //! newline automatically if the last new line was too long ago
      //! This is necessary if it is possible that an output string contained a new line
      void FlushLineAutoNewline();

      //! @brief update the current indentation string
      void UpdateIndentString();

    }; //class FixedLineWidthWriter

  } // namespace io
} // namespace bcl

#endif //BCL_IO_FIXED_LINE_WIDTH_WRITER_H_
