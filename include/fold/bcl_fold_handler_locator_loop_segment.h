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

#ifndef BCL_FOLD_HANDLER_LOCATOR_LOOP_SEGMENT_H_
#define BCL_FOLD_HANDLER_LOCATOR_LOOP_SEGMENT_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class HandlerLocatorLoopSegment
    //! @brief for creating LocatorLoopSegment objects from a file.
    //!
    //! @see @link example_fold_handler_locator_loop_segment.cpp @endlink
    //! @author alexanns
    //! @date Sep 7, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API HandlerLocatorLoopSegment :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! number of columns that the line containing the information should have
      static const size_t s_ExpectedLineColumns;

      //! the column number where the chain id should be in the line
      static const size_t s_ChainIDColumn;

      //! the column number where the starting residue of the sse should be in the line
      static const size_t s_SSEStartSeqIDColumn;

      //! the column number where the ending residue of the sse should be in the line
      static const size_t s_SSEEndSeqIDColumn;

      //! the column number where the boolean indicating rigid or not should be in the line
      static const size_t s_SSEIsRigidIDColumn;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new HandlerLocatorLoopSegment
      virtual HandlerLocatorLoopSegment *Clone() const;

      //! @brief virtual destructor
      virtual ~HandlerLocatorLoopSegment();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief HandleRead creates a LocatorLoopSegment from a stream
      //! @param ISTREAM the stream from which the LocatorLoopSegment will be created
      //! @return LocatorLoopSegment which was created from the contents of ISTREAM
      virtual LocatorLoopSegment HandleRead( std::istream &ISTREAM) const;

      //! @brief HandleRead creates a LocatorLoopSegment from a string
      //! @param STRING the string from which the LocatorLoopSegment will be created
      //! @return LocatorLoopSegment which was created from the contents of STRING
      LocatorLoopSegment HandleRead( const std::string &STRING) const;

      //! @brief GetFormat gives the format that this handler needs in order to work
      //! @return string which describes the format needed by this handler in order for it to work
      virtual std::string GetFormat() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

    private:

    }; // class HandlerLocatorLoopSegment

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_HANDLER_LOCATOR_LOOP_SEGMENT_H_ 
