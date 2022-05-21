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

#ifndef BCL_ASSEMBLE_LOCATOR_ATOM_COORDINATES_INTERFACE_H_
#define BCL_ASSEMBLE_LOCATOR_ATOM_COORDINATES_INTERFACE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_compare.h"
#include "biol/bcl_biol_aa_sequence.h"
#include "biol/bcl_biol_atom.h"
#include "find/bcl_find_locator_coordinates_interface.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorAtomCoordinatesInterface
    //! @brief Interface for an object that can locate coordinates as well as other structures the coordinates are in
    //! @details provides functionality for locating coordinates and higher structural levels
    //!
    //! @see @link example_assemble_locator_atom_coordinates_interface.cpp @endlink
    //! @author alexanns
    //! @date May 6, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorAtomCoordinatesInterface :
      public find::LocatorCoordinatesInterface< ProteinModel>
    {

    public:

      //! @brief Clone function
      //! @return pointer to new LocatorCoordinatesInterface
      virtual LocatorAtomCoordinatesInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief gives formatted string describing the locator
      //! @return formatted string describing the locator
      virtual std::string GetIdentification() const;

      //! @brief reads formatted string describing the locator
      //! @return formatted string describing the locator
      virtual std::istream &ReadIdentification( std::istream &ISTREAM) = 0;

      //! @brief gives identifier that can be used for selections in pymol
      //! @return gives string which is an identifier that can be used for selections in pymol
      virtual const std::string GetPymolName() const;

      //! @brief gives identifier that can be used for selections in pymol
      //! @return gives string which is an identifier that can be used for selections in pymol
      const std::string GetPymolAtomSelection() const;

      //! @brief gives identifier that can be used for selections in pymol
      //! @return gives string which is an identifier that can be used for selections in pymol
      const std::string GetPymolResidueSelection() const;

      //! @brief gives the chain id the locator corresponds to
      //! @return the chain id the locator corresponds to
      virtual char GetChainID() const = 0;

      //! @brief gives the seq id the locator corresponds to
      //! @return the seq id the locator corresponds to
      virtual int GetSeqID() const = 0;

      //! @brief gives the atom type the locator corresponds to
      //! @return the atom type the locator corresponds to
      virtual const biol::AtomType &GetAtomType() const = 0;

      //! @brief gives the aa type the locator corresponds to
      //! @return the aa type the locator corresponds to
      virtual const biol::AAType &GetAAType() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief locates an atom from a protein model
      //! @param PROTEIN_MODEL model from which the atom will be located
      //! @return siptr to atom that has been located from PROTEIN_MODEL
      virtual util::SiPtr< const biol::Atom> LocateAtom( const ProteinModel &PROTEIN_MODEL) const;

      //! @brief locates an atom from a protein model
      //! @param PROTEIN_MODEL model from which the atom will be located
      //! @return copy of atom that has been located from PROTEIN_MODEL
      virtual biol::Atom LocateAtomCopy( const ProteinModel &PROTEIN_MODEL) const;

      //! @brief locates an aa from a protein model
      //! @param PROTEIN_MODEL model from which the aa will be located
      //! @return siptr to aabase that has been located from PROTEIN_MODEL
      virtual util::SiPtr< const biol::AABase> LocateAA( const ProteinModel &PROTEIN_MODEL) const;

      //! @brief locates an aa from an sse
      //! @param SSE the sse from which the aa will be located
      //! @return siptr to aabase that has been located from SSE
      virtual util::SiPtr< const biol::AABase> LocateAA( const SSE &SSE);

      //! @brief locates an aa from a list of residues
      //! @param RESIDUE_LIST the list of residues which will be searched for the residue of interest
      //! @return siptr to aabase that has been located from RESIDUE_LIST
      virtual util::SiPtr< const biol::AABase>
      LocateAA( const util::SiPtrVector< const biol::AABase> &RESIDUE_LIST) const;

      //! @brief locates an sse from a protein model
      //! @param PROTEIN_MODEL model from which the sse will be found
      //! @return siptr to aa base that has been located from PROTEIN_MODEL
      virtual util::SiPtr< const SSE> LocateSSE( const ProteinModel &PROTEIN_MODEL) const;

      //! @brief locates the SSE this locator is in from a list of SSEs
      //! @param SSES the list of SSEs from which the sse this locator is within will be found
      //! @return siptr to sse that contains this locator
      virtual util::SiPtr< const SSE> LocateSSE( const util::SiPtrVector< const SSE> &SSES) const;

      //! @brief locates the SSE this locator is in from a set of SSEs
      //! @param SSES the list of SSEs from which the sse this locator is within will be found
      //! @return siptr to sse that contains this locator
      virtual util::SiPtr< const SSE> LocateSSE
      (
        const storage::Set< util::SiPtr< const SSE>, SSELessThanNoOverlap> &SSES
      ) const;

      //! @brief determines if this locator is within a given SSE
      //! @param SSE the sse which will be checked to see if this lcoator is within it
      //! @return bool true if this locator is within the SSE - false otherwise
      virtual bool IsWithin( const SSE &SSE) const;

    public:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! @class PtrLessThan
      //! @brief helper binary operator struct for comparing objects behind pointers to this interface
      //! @author alexanns
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct PtrLessThan
      {
        //! @brief binary operator for comparing objects behind pointers to this interface
        //! @param LHS first  pointer to LocatorAtomCoordinatesInterface
        //! @param RHS second pointer to LocatorAtomCoordinatesInterface
        //! @return bool true if LHS object behind pointer is less than RHS object behind pointer - false otherwise
        bool operator()
        (
          const util::PtrInterface< LocatorAtomCoordinatesInterface> &LHS,
          const util::PtrInterface< LocatorAtomCoordinatesInterface> &RHS
        ) const;
      };

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //! @class PtrResidueLessThan
      //! @brief helper binary operator struct for comparing objects behind pointers to this interface
      //!        compares the chain id and residue seq id
      //! @author alexanns
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct PtrResidueLessThan
      {
        //! @brief binary operator for comparing objects behind pointers to this interface
        //! @param LHS first  pointer to LocatorAtomCoordinatesInterface
        //! @param RHS second pointer to LocatorAtomCoordinatesInterface
        //! @return bool true if LHS object behind pointer is less than RHS object behind pointer - false otherwise
        bool operator()
        (
          const util::PtrInterface< LocatorAtomCoordinatesInterface> &LHS,
          const util::PtrInterface< LocatorAtomCoordinatesInterface> &RHS
        ) const;
      };

      //! @brief creates a set of pointers to interface locators from a sequence
      //! @tparam t_LocatorType the type of locator that should be behind the pointer
      //! @param SEQUENCE the sequence for which a set of locators will be created
      //! @return set which has locators created from the sequence
      template< typename t_LocatorType> static storage::Set
      <
        util::ShPtr< LocatorAtomCoordinatesInterface>,
        LocatorAtomCoordinatesInterface::PtrLessThan
      > CreateLocators( const biol::AASequence &SEQUENCE);

      //! @brief creates a pointer to a locator interface from an aa base
      //! @tparam t_LocatorType the type of locator that should be behind the pointer
      //! @param AA_BASE aa base from which the locator will be created
      //! @return shptr to locator with type t_LocatorType behind it
      template< typename t_LocatorType> static util::ShPtr< LocatorAtomCoordinatesInterface>
      CreateLocator( const biol::AABase &AA_BASE);

      //! @brief gives a name from an atom locator with no spaces
      //! @param LOCATOR_A first locator to get name from
      //! @param LOCATOR_B second locator to get name from
      //! @return string that represents the pair of locators
      static std::string GetNameFromPair
      (
        const LocatorAtomCoordinatesInterface &LOCATOR_A,
        const LocatorAtomCoordinatesInterface &LOCATOR_B
      );

    }; // class LocatorAtomCoordinatesInterface

    //! @brief less than operator for comparing two LocatorAtomCoordinatesInterface
    //! @param LHS the first LocatorAtomCoordinatesInterface
    //! @param RHS the second LocatorAtomCoordinatesInterface
    //! @return bool true if LHS is less than RHS - false otherwise
    bool operator <( const LocatorAtomCoordinatesInterface &LHS, const LocatorAtomCoordinatesInterface &RHS);

    //! @brief not equal operator for comparing two LocatorAtomCoordinatesInterface
    //! @param LHS the first LocatorAtomCoordinatesInterface
    //! @param RHS the second LocatorAtomCoordinatesInterface
    //! @return bool true if LHS not equal to RHS - false otherwise
    bool operator !=( const LocatorAtomCoordinatesInterface &LHS, const LocatorAtomCoordinatesInterface &RHS);

    //! @brief creates a set of pointers to interface locators from a sequence
    //! @tparam t_LocatorType the type of locator that should be behind the pointer
    //! @param SEQUENCE the sequence for which a set of locators will be created
    //! @return set which has locators created from the sequence
    template< typename t_LocatorType>
    storage::Set
    <
      util::ShPtr< LocatorAtomCoordinatesInterface>,
      LocatorAtomCoordinatesInterface::PtrLessThan
    > LocatorAtomCoordinatesInterface::CreateLocators( const biol::AASequence &SEQUENCE)
    {
      // set to hold created locators
      storage::Set
      <
        util::ShPtr< LocatorAtomCoordinatesInterface>,
        LocatorAtomCoordinatesInterface::PtrLessThan
      > locators;

      // iterate through the sequence
      for
      (
        biol::AASequence::const_iterator aa_itr_a( SEQUENCE.Begin()), aa_itr_end( SEQUENCE.End());
        aa_itr_a != aa_itr_end;
        ++aa_itr_a
      )
      {
        // insert created locator
        locators.Insert( CreateLocator< t_LocatorType>( **aa_itr_a));
      }

      // return the set of locators
      return locators;
    }

    //! @brief creates a pointer to a locator interface from an aa base
    //! @tparam t_LocatorType the type of locator that should be behind the pointer
    //! @param AA_BASE aa base from which the locator will be created
    //! @return shptr to locator with type t_LocatorType behind it
    template< typename t_LocatorType> util::ShPtr< LocatorAtomCoordinatesInterface>
    LocatorAtomCoordinatesInterface::CreateLocator( const biol::AABase &AA_BASE)
    {
      // create locator with desired type
      const util::ShPtr< LocatorAtomCoordinatesInterface> locator
      (
        new t_LocatorType
        (
          AA_BASE.GetChainID(), AA_BASE.GetSeqID(), AA_BASE.GetFirstSidechainAtom().GetType(), AA_BASE.GetType()
        )
      );

      // return locator with desired type behind the pointer
      return locator;
    }

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_LOCATOR_ATOM_COORDINATES_INTERFACE_H_
