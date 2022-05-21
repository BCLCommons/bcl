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

#ifndef BCL_FOLD_LOOP_SEGMENT_H_
#define BCL_FOLD_LOOP_SEGMENT_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LoopSegment
    //! @brief is for representing a segment of a loop.
    //! @details  The segment has an sse and a can be specified as being rigid or not. If the LoopSegment is rigid,
    //! the dihedral angles of the sse will not be changed. Otherwise, the dihedral angles can be changed.
    //!
    //! @see @link example_fold_loop_segment.cpp @endlink
    //! @author alexanns
    //! @date Aug 28, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LoopSegment :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the sse which the loop segment refers to
      util::ShPtr< assemble::SSE> m_SSE;

      //! boolean indicating if the dihedral angles of the LoopSegment can be changed or not
      bool m_IsRigid;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LoopSegment();

      //! @brief constructor taking member variable parameters
      //! @param SSE the sse which the loop segment refers to
      //! @param IS_RIGID boolean indicating if the dihedral angles of the LoopSegment can be changed or not
      LoopSegment( const util::ShPtr< assemble::SSE> &SSE, const bool IS_RIGID);

      //! @brief copy constructor
      //! @param OTHER loop segment to be copied
      LoopSegment( const LoopSegment &OTHER);

      //! @brief Clone function
      //! @return pointer to new LoopSegment
      LoopSegment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief GetSSE gives the sse which the loop segment refers to
      //! @return "m_SSE" the sse which the loop segment refers to
      const util::ShPtr< assemble::SSE> &GetSSE() const;

      //! @brief GetSSE gives the sse which the loop segment refers to
      //! @return "m_SSE" the sse which the loop segment refers to
      util::ShPtr< assemble::SSE> &GetSSE();

      //! @brief GetSSE gives the sse which the loop segment refers to
      //! @return "m_SSE" the sse which the loop segment refers to
      assemble::SSE &GetSSEReference();

      //! @brief GetSSE gives the sse which the loop segment refers to
      //! @return "m_SSE" the sse which the loop segment refers to
      const assemble::SSE &GetConstSSEReference() const;

      //! @brief IsRigid gives boolean indicating if the dihedral angles of the LoopSegment can be changed or not
      //! @return "m_IsRigid" the boolean indicating if the dihedral angles of the LoopSegment can be changed or not
      bool IsRigid() const;

      //! @brief gives iterator to residue provided
      //! @param AA_BASE residue for which an iterator will be given
      //! @return iterator to the provided residue
      biol::AASequence::iterator GetAAIterator( const biol::AABase &AA_BASE);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class LoopSegment

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_LOOP_SEGMENT_H_ 
