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

#ifndef BCL_FOLD_LOCATOR_LOOP_SEGMENT_H_
#define BCL_FOLD_LOCATOR_LOOP_SEGMENT_H_

// include the namespace header
#include "bcl_fold.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorLoopSegment
    //! @brief is for finding and creating a LoopSegment from a domain.
    //!
    //! @see @link example_fold_locator_loop_segment.cpp @endlink
    //! @author alexanns
    //! @date Aug 28, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorLoopSegment :
      public find::LocatorInterface< LoopSegment, assemble::DomainInterface>
    {

    private:

    //////////
    // data //
    //////////

      //! locator to locate the SSE that will be a part of the loop segment
      assemble::LocatorSSE m_SSELocator;

      //! indicates if the dihedral angles of the SSE must be kept rigid (true) or if they can be changed (false)
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
      LocatorLoopSegment();

      //! @brief constructor taking member variable parameters
      //! @param SSE_LOCATOR locator to locate the SSE that will be a part of the loop segment
      //! @param IS_RIGID indicates if the dihedral angles of the SSE must be kept rigid (true), false otherwise
      LocatorLoopSegment( const assemble::LocatorSSE &SSE_LOCATOR, const bool IS_RIGID);

      //! @brief Clone function
      //! @return pointer to new LoopSegment
      LocatorLoopSegment *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief GetIdentification
      //! @return GetIdentification
      const std::string GetIdentification() const;

      //! @brief GetLocatorSSE gives the locator to locate the SSE that will be a part of the loop segment
      //! @return the locator to locate the SSE that will be a part of the loop segment
      const assemble::LocatorSSE &GetLocatorSSE() const;

      //! @brief IsRigid indicates if the dihedral angles of the SSE must be kept rigid or if they can be changed
      //! @return true if the dihedral angles of the sse must be kept rigid, false otherwise
      const bool IsRigid() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Locate locates and creates the loop segment in a domain
      //! @param SSE_DOMAIN is the domain from which the loop segment will be located and created
      //! @return a loop segment that has been created from SSE_DOMAIN
      LoopSegment Locate( const assemble::DomainInterface &SSE_DOMAIN) const;

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

    }; // class LocatorLoopSegment

  } // namespace fold
} // namespace bcl

#endif // BCL_FOLD_LOCATOR_LOOP_SEGMENT_H_ 
