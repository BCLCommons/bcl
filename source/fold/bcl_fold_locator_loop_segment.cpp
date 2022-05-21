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
#include "fold/bcl_fold_locator_loop_segment.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_loop_segment.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LocatorLoopSegment::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorLoopSegment())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LocatorLoopSegment::LocatorLoopSegment() :
      m_SSELocator(),
      m_IsRigid()
    {
    }

    //! @brief constructor taking member variable parameters
    //! @param SSE_LOCATOR locator to locate the SSE that will be a part of the loop segment
    //! @param IS_RIGID indicates if the dihedral angles of the SSE must be kept rigid (true), false otherwise
    LocatorLoopSegment::LocatorLoopSegment( const assemble::LocatorSSE &SSE_LOCATOR, const bool IS_RIGID) :
      m_SSELocator( SSE_LOCATOR),
      m_IsRigid( IS_RIGID)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LoopSegment
    LocatorLoopSegment *LocatorLoopSegment::Clone() const
    {
      return new LocatorLoopSegment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &LocatorLoopSegment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetIdentification
    //! @return GetIdentification
    const std::string LocatorLoopSegment::GetIdentification() const
    {
      return m_SSELocator.GetIdentification() + " is rigid? " + util::Format()( m_IsRigid);
    }

    //! @brief GetLocatorSSE gives the locator to locate the SSE that will be a part of the loop segment
    //! @return the locator to locate the SSE that will be a part of the loop segment
    const assemble::LocatorSSE &LocatorLoopSegment::GetLocatorSSE() const
    {
      return m_SSELocator;
    }

    //! @brief IsRigid indicates if the dihedral angles of the SSE must be kept rigid or if they can be changed
    //! @return true if the dihedral angles of the sse must be kept rigid, false otherwise
    const bool LocatorLoopSegment::IsRigid() const
    {
      return m_IsRigid;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Locate locates and creates the loop segment in a domain
    //! @param SSE_DOMAIN is the domain from which the loop segment will be located and created
    //! @return a loop segment that has been created from SSE_DOMAIN
    LoopSegment LocatorLoopSegment::Locate( const assemble::DomainInterface &SSE_DOMAIN) const
    {
      // find the sse denoted by "m_SSELocator" in "PROTEIN_MODEL"
      const util::SiPtr< const assemble::SSE> located_sse( m_SSELocator.Locate( SSE_DOMAIN));

      // make sure that the sse was able to be located
      BCL_Assert
      (
        located_sse.IsDefined(), "sse could not be located with chain " +
        util::Format()( m_SSELocator.GetChainID()) + " starting at " +
        util::Format()( m_SSELocator.GetSSEID().First()) + " and ending at " +
        util::Format()( m_SSELocator.GetSSEID().Second())
      );

      // get a shptr to a clone of "located_sse"
      const util::ShPtr< assemble::SSE> new_sse( located_sse->Clone());

      // make sure "new_sse" is defined
      BCL_Assert( new_sse.IsDefined(), "new_sse is not defined");

      // create a LoopSegment out of "new_sse" and "m_IsRigid"
      const LoopSegment loop_segment( new_sse, m_IsRigid);

      // make sure the sse in "loop_segment" is defined after construction
      BCL_Assert( loop_segment.GetSSE().IsDefined(), "sse not defined");

      // return the created LoopSegment
      return loop_segment;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LocatorLoopSegment::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSELocator, ISTREAM);
      io::Serialize::Read( m_IsRigid, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LocatorLoopSegment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSELocator, OSTREAM, INDENT);
      io::Serialize::Write( m_IsRigid, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  } // namespace fold
} // namespace bcl
