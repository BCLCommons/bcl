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
#include "fold/bcl_fold_handler_locator_loop_segment.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_locator_loop_segment.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! number of columns that the line containing the information should have
    const size_t HandlerLocatorLoopSegment::s_ExpectedLineColumns( 4);

    //! the column number where the chain id should be in the line
    const size_t HandlerLocatorLoopSegment::s_ChainIDColumn( 0);

    //! the column number where the starting residue of the sse should be in the line
    const size_t HandlerLocatorLoopSegment::s_SSEStartSeqIDColumn( 1);

    //! the column number where the ending residue of the sse should be in the line
    const size_t HandlerLocatorLoopSegment::s_SSEEndSeqIDColumn( 2);

    //! the column number where the boolean indicating rigid or not should be in the line
    const size_t HandlerLocatorLoopSegment::s_SSEIsRigidIDColumn( 3);

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> HandlerLocatorLoopSegment::s_Instance
    (
      GetObjectInstances().AddInstance( new HandlerLocatorLoopSegment())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new HandlerLocatorLoopSegment
    HandlerLocatorLoopSegment *HandlerLocatorLoopSegment::Clone() const
    {
      return new HandlerLocatorLoopSegment( *this);
    }

    //! @brief virtual destructor
    HandlerLocatorLoopSegment::~HandlerLocatorLoopSegment()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &HandlerLocatorLoopSegment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief HandleRead creates a LocatorLoopSegment from a stream
    //! @param ISTREAM the stream from which the LocatorLoopSegment will be created
    //! @return LocatorLoopSegment which was created from the contents of ISTREAM
    LocatorLoopSegment HandlerLocatorLoopSegment::HandleRead( std::istream &ISTREAM) const
    {
      // create string to hold a line from ISTREAM
      std::string line;

      // get the current line from ISTREAM and put it into "line"
      std::getline( ISTREAM, line);

      // create LocatorLoopSegment from "line"
      return HandleRead( line);
    }

    //! @brief HandleRead creates a LocatorLoopSegment from a string
    //! @param STRING the string from which the LocatorLoopSegment will be created
    //! @return LocatorLoopSegment which was created from the contents of STRING
    LocatorLoopSegment HandlerLocatorLoopSegment::HandleRead( const std::string &STRING) const
    {
      // trim and split "STRING" and put into vector "split_line"
      const storage::Vector< std::string> split_line( util::SplitString( util::TrimString( STRING)));

      // make sure "split_line" has the correct size
      BCL_Assert
      (
        split_line.GetSize() == s_ExpectedLineColumns,
        "split line should have " + util::Format()( s_ExpectedLineColumns) + " but has " +
        util::Format()( split_line.GetSize()) + "\n" + GetFormat()
      );

      // get the chain id
      const std::string &chain_id_string( split_line( s_ChainIDColumn));

      // make sure chain id consists of three characters
      BCL_Assert
      (
        chain_id_string.size() == 3,
        "chain id column " + util::Format()( s_ChainIDColumn) +
        " should have length 3 with chain id (e.g. 'A') but is " + chain_id_string + "\n" + GetFormat()
      );

      // make sure the chain id has the correct formatting
      BCL_Assert
      (
        chain_id_string[ 0] == '\'' && chain_id_string[ 2] == '\'',
        "chain id column " + util::Format()( s_ChainIDColumn) +
        " should have length 3 with chain id (e.g. 'A') but is " + chain_id_string + "\n" + GetFormat()
      );

      // get the chain id character
      const char chain_id( chain_id_string[ 1]);

      // get the starting residue of the sse
      const std::string &sse_start_seq_id_string( split_line( s_SSEStartSeqIDColumn));

      // get the numerical value seq id of the starting residue of the sse
      const int sse_start_seq_id( util::ConvertStringToNumericalValue< int>( sse_start_seq_id_string));

      // make sure "sse_start_seq_id" is defined
      BCL_Assert
      (
        util::IsDefined( sse_start_seq_id),
        "starting seq id of SSE is undefined. Read in value is  \"" + sse_start_seq_id_string + "\""
      );

      // get the ending residue of the sse
      const std::string sse_end_seq_id_string( split_line( s_SSEEndSeqIDColumn));

      // get the numerical value seq id of the ending residue of the sse
      const int sse_end_seq_id( util::ConvertStringToNumericalValue< int>( sse_end_seq_id_string));

      // make sure "sse_end_seq_id" is defined
      BCL_Assert
      (
        util::IsDefined( sse_end_seq_id),
        "ending seq id of SSE is undefined. Read in value is  \"" + sse_end_seq_id_string + "\""
      );

      // get the boolean indicating rigid or not
      const std::string &sse_is_rigid_string( split_line( s_SSEIsRigidIDColumn));

      // convert the is_rigid string to a boolean value
      const bool sse_is_rigid( util::ConvertStringToBoolean( sse_is_rigid_string));

      // create a LocatorLoopSegment from "chain_id", "sse_start_seq_id", "sse_end_seq_id", and "sse_is_rigid"
      const LocatorLoopSegment loop_segment_locator
      (
        assemble::LocatorSSE( chain_id, sse_start_seq_id, sse_end_seq_id), sse_is_rigid
      );

      // return the created LocatorLoopSegment
      return loop_segment_locator;
    }

    //! @brief GetFormat gives the format that this handler needs in order to work
    //! @return string which describes the format needed by this handler in order for it to work
    std::string HandlerLocatorLoopSegment::GetFormat() const
    {
      return "<chain_id> <sse_start_seq_id> <sse_end_seq_id> <true/false>\n"
        "where \"true\" or \"false\" indicate whether the sse is rigid (true) or not (false)\n"
        "chain id should be surrounded by '', e.g. 'B'"
        "example\n'A' 32 34 true\nexample\n'A' 32 34 false\n";
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &HandlerLocatorLoopSegment::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &HandlerLocatorLoopSegment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace fold
} // namespace bcl
