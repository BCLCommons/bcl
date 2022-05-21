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
#include "assemble/bcl_assemble_locator_sse_terminus_residue.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LocatorSSETerminusResidue::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorSSETerminusResidue())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LocatorSSETerminusResidue::LocatorSSETerminusResidue() :
      m_ChainID(),
      m_StartSSE(),
      m_Terminus()
    {
    }

    //! @brief construct from a ChainID, and the SeqID of the amino acid from which you want to grow
    //! @param CHAINID char which indicates the chain
    //! @param SSE_START residue from which to be grown
    //! @param GROW_C_TO_N if true, grow from C to N terminus
    LocatorSSETerminusResidue::LocatorSSETerminusResidue
    (
      const char CHAINID, const int SSE_START, const biol::AASequenceFlexibility::SequenceDirection &TERMINUS
    ) :
      m_ChainID( CHAINID),
      m_StartSSE( SSE_START),
      m_Terminus( TERMINUS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LocatorSSETerminusResidue
    LocatorSSETerminusResidue *LocatorSSETerminusResidue::Clone() const
    {
      return new LocatorSSETerminusResidue( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorSSETerminusResidue::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetIdentification
    //! @return the GetIdentification
    const std::string LocatorSSETerminusResidue::GetIdentification() const
    {
      return util::Format()( m_ChainID) + " " + util::Format()( m_StartSSE) + " " +
        util::Format()( m_Terminus.GetString());
    }

    //! @brief gives the residue that identifies the sse of interest
    //! @return returns the locator
    int LocatorSSETerminusResidue::GetSSEID() const
    {
      return m_StartSSE;
    }

    //! @brief if the N or the C terminus is to be located
    //! @return true if the cterminus of the sse is to be located, false if n-terminus
    biol::AASequenceFlexibility::SequenceDirection LocatorSSETerminusResidue::GetLocateCTerminus() const
    {
      return m_Terminus;
    }

    //! @brief gives an identifying string for this locator
    //! @return returns an identifying string for this locator
    std::string LocatorSSETerminusResidue::GetSSEIDString() const
    {
      return
          std::string
          (
            util::Format()(             m_ChainID) + " " +
            util::Format().W( 4).R()( m_StartSSE) + " " +
            util::Format()( m_Terminus.GetString())
          );
    }

    //! returns the chain locator
    //! @return returns the const reference to the chain locator
    const char &LocatorSSETerminusResidue::GetChainID() const
    {
      return m_ChainID;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Locate a sse in the domain
    //! @param SSE_DOMAIN domain which the LocatorSSE refers to
    //! @return returns SiPtr to the SSE denoted by the LocatorSSE
    util::SiPtr< const SSE> LocatorSSETerminusResidue::Locate( const DomainInterface &SSE_DOMAIN) const
    {
      // get all sses
      const util::SiPtrVector< const SSE> sses( SSE_DOMAIN.GetSSEs());

      // iterate over the SSEs
      for
      (
        util::SiPtrVector< const SSE>::const_iterator itr( sses.Begin()), itr_end( sses.End());
        itr != itr_end;
        ++itr
      )
      {
        // check to see if the SeqIDs of the first and last amino acids of current SSE match m_Start_SSE and m_End_SSE
        if
        (
          ( *itr)->GetChainID()             == m_ChainID &&
          ( *itr)->GetFirstAA()->GetSeqID() == m_StartSSE &&
          ( m_Terminus == biol::AASequenceFlexibility::e_NTerminal)
        )
        {
          // if so then return this SSE
          return *itr;
        }
        else if
        (
          ( *itr)->GetChainID()             == m_ChainID &&
          ( *itr)->GetLastAA()->GetSeqID() == m_StartSSE &&
          ( m_Terminus == biol::AASequenceFlexibility::e_CTerminal)
        )
        {
          // if so then return this SSE
          return *itr;
        }
      }

      // if not found, return empty SiPtr< const SSE>
      return util::SiPtr< const SSE>();
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
    std::istream &LocatorSSETerminusResidue::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ChainID        , ISTREAM);
      io::Serialize::Read( m_StartSSE       , ISTREAM);
      io::Serialize::Read( m_Terminus, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &LocatorSSETerminusResidue::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ChainID        , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_StartSSE       , OSTREAM)         << '\t';
      io::Serialize::Write( m_Terminus, OSTREAM);
      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace assemble
} // namespace bcl
