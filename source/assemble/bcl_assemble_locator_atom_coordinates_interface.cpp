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
#include "assemble/bcl_assemble_locator_atom_coordinates_interface.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_atom.h"
#include "assemble/bcl_assemble_sse.h"
#include "util/bcl_util_string_replacement.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  /////////////////
  // data access //
  /////////////////

    //! @brief gives formatted string describing the locator
    //! @return formatted string describing the locator
    std::string LocatorAtomCoordinatesInterface::GetIdentification() const
    {
      return util::Format()( GetChainID()) + " " +
        util::Format()( GetSeqID()) + " " + GetAtomType().GetName() + " " + GetAAType().GetName();
    }

    //! @brief gives identifier that can be used for selections in pymol
    //! @return gives string which is an identifier that can be used for selections in pymol
    const std::string LocatorAtomCoordinatesInterface::GetPymolName() const
    {
      return util::Format()( GetChainID()) + "_" +
        util::Format()( GetSeqID()) + "_" + GetAtomType().GetName();
    }

    //! @brief gives identifier that can be used for selections in pymol
    //! @return gives string which is an identifier that can be used for selections in pymol
    const std::string LocatorAtomCoordinatesInterface::GetPymolAtomSelection() const
    {
      return " chain " + util::Format()( GetChainID()) + " and resi " +
        util::Format()( GetSeqID()) + " and name " + GetAtomType().GetName();
    }

    //! @brief gives identifier that can be used for selections in pymol
    //! @return gives string which is an identifier that can be used for selections in pymol
    const std::string LocatorAtomCoordinatesInterface::GetPymolResidueSelection() const
    {
      return " chain " + util::Format()( GetChainID()) + " and resi " +
        util::Format()( GetSeqID());
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief locates an atom from a protein model
    //! @param PROTEIN_MODEL model from which the atom will be located
    //! @return siptr to atom that has been located from PROTEIN_MODEL
    util::SiPtr< const biol::Atom>
    LocatorAtomCoordinatesInterface::LocateAtom( const ProteinModel &PROTEIN_MODEL) const
    {
        return LocatorAtom::LocateAtomFromModel( PROTEIN_MODEL, GetChainID(), GetSeqID(), GetAtomType());
    }

    //! @brief locates an atom from a protein model
    //! @param PROTEIN_MODEL model from which the atom will be located
    //! @return copy of atom that has been located from PROTEIN_MODEL
    biol::Atom LocatorAtomCoordinatesInterface::LocateAtomCopy( const ProteinModel &PROTEIN_MODEL) const
    {
      // try to locate atom
      const util::SiPtr< const biol::Atom> atom
      (
        LocatorAtom::LocateAtomFromModel( PROTEIN_MODEL, GetChainID(), GetSeqID(), GetAtomType())
      );

      // true if atom could not be found
      if( !atom.IsDefined())
      {
        // return empty atom
        return biol::Atom();
      }

      // atom is defined so return atom behind pointer
      return *atom;
    }

    //! @brief locates an aa from a protein model
    //! @param PROTEIN_MODEL model from which the aa will be located
    //! @return siptr to aabase that has been located from PROTEIN_MODEL
    util::SiPtr< const biol::AABase>
    LocatorAtomCoordinatesInterface::LocateAA( const ProteinModel &PROTEIN_MODEL) const
    {
      return LocatorAA( GetChainID(), GetSeqID()).Locate( PROTEIN_MODEL);
    }

    //! @brief locates an aa from an sse
    //! @param SSE the sse from which the aa will be located
    //! @return siptr to aabase that has been located from SSE
    util::SiPtr< const biol::AABase> LocatorAtomCoordinatesInterface::LocateAA( const SSE &SSE)
    {
      return LocateAA( SSE.GetMembers());
    }

    //! @brief locates an aa from a list of residues
    //! @param RESIDUE_LIST the list of residues which will be searched for the residue of interest
    //! @return siptr to aabase that has been located from RESIDUE_LIST
    util::SiPtr< const biol::AABase>
    LocatorAtomCoordinatesInterface::LocateAA( const util::SiPtrVector< const biol::AABase> &RESIDUE_LIST) const
    {
      return LocatorAA( GetChainID(), GetSeqID()).Locate( RESIDUE_LIST);
    }

    //! @brief locates an sse from a protein model
    //! @param PROTEIN_MODEL model from which the sse will be found
    //! @return siptr to aa base that has been located from PROTEIN_MODEL
    util::SiPtr< const SSE>
    LocatorAtomCoordinatesInterface::LocateSSE( const ProteinModel &PROTEIN_MODEL) const
    {
      return LocatorAA( GetChainID(), GetSeqID()).LocateSSE( PROTEIN_MODEL);
    }

    //! @brief locates the SSE this locator is in from a list of SSEs
    //! @param SSES the list of SSEs from which the sse this locator is within will be found
    //! @return siptr to sse that contains this locator
    util::SiPtr< const SSE> LocatorAtomCoordinatesInterface::LocateSSE
    (
      const util::SiPtrVector< const SSE> &SSES
    ) const
    {
      // iterate through the sses
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( SSES.Begin()), sse_itr_end( SSES.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // true if the locator indicates a part of the current sse
        if( IsWithin( **sse_itr))
        {
          // return siptr to current sse
          return *sse_itr;
        }
      }

      // locator does not refer to any part of the provided sses
      return util::SiPtr< const SSE>();
    }

    //! @brief locates the SSE this locator is in from a set of SSEs
    //! @param SSES the list of SSEs from which the sse this locator is within will be found
    //! @return siptr to sse that contains this locator
    util::SiPtr< const SSE> LocatorAtomCoordinatesInterface::LocateSSE
    (
      const storage::Set< util::SiPtr< const SSE>, SSELessThanNoOverlap> &SSES
    ) const
    {
      // iterate through the sses
      for
      (
        storage::Set< util::SiPtr< const SSE>, SSELessThanNoOverlap>::const_iterator
          sse_itr( SSES.Begin()), sse_itr_end( SSES.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // true if the locator indicates a part of the current sse
        if( IsWithin( **sse_itr))
        {
          // return siptr to current sse
          return *sse_itr;
        }
      }

      // locator does not refer to any part of the provided sses
      return util::SiPtr< const SSE>();
    }

    //! @brief determines if this locator is within a given SSE
    //! @param SSE the sse which will be checked to see if this lcoator is within it
    //! @return bool true if this locator is within the SSE - false otherwise
    bool LocatorAtomCoordinatesInterface::IsWithin( const SSE &SSE) const
    {
      // get information about this locator and the SSE
      const char chain_id( GetChainID());
      const int seq_id( GetSeqID());
      const int sse_start( SSE.GetFirstAA()->GetSeqID());
      const int sse_end( SSE.GetLastAA()->GetSeqID());
      const char sse_chain_id( SSE.GetChainID());

      return chain_id == sse_chain_id && seq_id >= sse_start && seq_id <= sse_end;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief gives a name from an atom locator with no spaces
    //! @param LOCATOR_A first locator to get name from
    //! @param LOCATOR_B second locator to get name from
    //! @return string that represents the pair of locators
    std::string LocatorAtomCoordinatesInterface::GetNameFromPair
    (
      const LocatorAtomCoordinatesInterface &LOCATOR_A,
      const LocatorAtomCoordinatesInterface &LOCATOR_B
    )
    {
      // string replacer to remove undefined atom names
      static const util::StringReplacement s_undefined_replacer( util::StringReplacement::e_Any, "Undefined", "");

      // first locator name
      const std::string name_a
      (
        util::Format()( LOCATOR_A.GetChainID()) + util::Format()( LOCATOR_A.GetSeqID())
        + LOCATOR_A.GetAtomType().GetName()
      );

      // second locator name
      const std::string name_b
      (
        util::Format()( LOCATOR_B.GetChainID()) + util::Format()( LOCATOR_B.GetSeqID())
        + LOCATOR_B.GetAtomType().GetName()
      );

      // name of pair
      std::string name( name_a + name_b);

      // remove undefined atom names
      s_undefined_replacer.ReplaceEachIn( name);

      // return name
      return name;
    }

    //! @brief less than operator for comparing two LocatorAtomCoordinatesInterface
    //! @param LHS the first LocatorAtomCoordinatesInterface
    //! @param RHS the second LocatorAtomCoordinatesInterface
    //! @return bool true if LHS is less than RHS - false otherwise
    bool operator <( const LocatorAtomCoordinatesInterface &LHS, const LocatorAtomCoordinatesInterface &RHS)
    {
      return
        LocatorAtom( LHS.GetChainID(), LHS.GetSeqID(), LHS.GetAtomType()) <
        LocatorAtom( RHS.GetChainID(), RHS.GetSeqID(), RHS.GetAtomType());
    }

    //! @brief not equal operator for comparing two LocatorAtomCoordinatesInterface
    //! @param LHS the first LocatorAtomCoordinatesInterface
    //! @param RHS the second LocatorAtomCoordinatesInterface
    //! @return bool true if LHS not equal to RHS - false otherwise
    bool operator !=( const LocatorAtomCoordinatesInterface &LHS, const LocatorAtomCoordinatesInterface &RHS)
    {
      return
        LocatorAtom( LHS.GetChainID(), LHS.GetSeqID(), LHS.GetAtomType()) !=
        LocatorAtom( RHS.GetChainID(), RHS.GetSeqID(), RHS.GetAtomType());
    }

    //! @brief helper binary operator struct for comparing objects behind pointers to this interface
    //! @param LHS first  pointer to LocatorAtomCoordinatesInterface
    //! @param RHS second pointer to LocatorAtomCoordinatesInterface
    //! @return bool true if LHS object behind pointer is less than RHS object behind pointer - false otherwise
    bool LocatorAtomCoordinatesInterface::PtrLessThan::operator()
    (
      const util::PtrInterface< LocatorAtomCoordinatesInterface> &LHS,
      const util::PtrInterface< LocatorAtomCoordinatesInterface> &RHS
    ) const
    {
      return *LHS < *RHS;
    }

    //! @brief helper binary operator struct for comparing objects behind pointers to this interface
    //! @param LHS first  pointer to LocatorAtomCoordinatesInterface
    //! @param RHS second pointer to LocatorAtomCoordinatesInterface
    //! @return bool true if LHS object behind pointer is less than RHS object behind pointer - false otherwise
    bool LocatorAtomCoordinatesInterface::PtrResidueLessThan::operator()
    (
      const util::PtrInterface< LocatorAtomCoordinatesInterface> &LHS,
      const util::PtrInterface< LocatorAtomCoordinatesInterface> &RHS
    ) const
    {
      return LocatorAA( LHS->GetChainID(), LHS->GetSeqID()) < LocatorAA( RHS->GetChainID(), RHS->GetSeqID());
    }

  } // namespace assemble
  
} // namespace bcl
