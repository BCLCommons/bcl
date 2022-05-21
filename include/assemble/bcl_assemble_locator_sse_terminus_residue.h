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

#ifndef BCL_ASSEMBLE_LOCATOR_SSE_TERMINUS_RESIDUE_H_
#define BCL_ASSEMBLE_LOCATOR_SSE_TERMINUS_RESIDUE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_sequence_flexibility.h"
#include "find/bcl_find_locator_interface.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorSSETerminusResidue
    //! @brief method for locating the teminal residue of an sse based on a single residue
    //! @details
    //!
    //! @see @link example_assemble_locator_sse_terminus_residue.cpp @endlink
    //! @author bitterd
    //! @date Aug 4, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorSSETerminusResidue :
      public find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface>
    {

    private:

    //////////
    // data //
    //////////

      char m_ChainID;         //!< Chain id
      int  m_StartSSE;        //!< the SeqID of the first or last amino acid in the SSE

      //! indicates if the residue is at the c or n terminus of the sse
      biol::AASequenceFlexibility::DirectionEnum m_Terminus;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LocatorSSETerminusResidue();

      //! @brief construct from a ChainID, and the SeqID of the amino acid from which you want to grow
      //! @param CHAINID char which indicates the chain
      //! @param SSE_START residue from which to be grown
      //! @param TERMINUS if true, grow from C to N terminus
      LocatorSSETerminusResidue
      (
        const char CHAINID, const int SSE_START, const biol::AASequenceFlexibility::SequenceDirection &TERMINUS
      );

      //! @brief Clone function
      //! @return pointer to new LocatorSSETerminusResidue
      LocatorSSETerminusResidue *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief GetIdentification
      //! @return the GetIdentification
      const std::string GetIdentification() const;

      //! @brief gives the residue that identifies the sse of interest
      //! @return returns the locator
      int GetSSEID() const;

      //! @brief if the N or the C terminus is to be located
      //! @return true if the cterminus of the sse is to be located, false if n-terminus
      biol::AASequenceFlexibility::SequenceDirection GetLocateCTerminus() const;

      //! @brief gives an identifying string for this locator
      //! @return returns an identifying string for this locator
      std::string GetSSEIDString() const;

      //! returns the chain locator
      //! @return returns the const reference to the chain locator
      const char &GetChainID() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Locate a sse in the domain
      //! @param SSE_DOMAIN domain which the LocatorSSE refers to
      //! @return returns SiPtr to the SSE denoted by the LocatorSSE
      util::SiPtr< const SSE> Locate( const DomainInterface &SSE_DOMAIN) const;

    ///////////////
    // operators //
    ///////////////

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

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class LocatorSSETerminusResidue

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_LOCATOR_SSE_TERMINUS_RESIDUE_H_ 
