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

#ifndef BCL_CHEMISTRY_POSSIBLE_ATOM_TYPES_FOR_ATOM_H_
#define BCL_CHEMISTRY_POSSIBLE_ATOM_TYPES_FOR_ATOM_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_conformational_interface.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PossibleAtomTypesForAtom
    //! @brief A helper class by which AtomTypes can return all possible atom types for a given atom in a structure
    //!        that is easily accessed by orbital type
    //!
    //! @see @link example_chemistry_possible_atom_types_for_atom.cpp @endlink
    //! @author mendenjl
    //! @date Aug 27, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PossibleAtomTypesForAtom :
      public util::ObjectInterface
    {
    private:

    //////////
    // data //
    //////////

      storage::Vector< size_t>   m_NumberAtomTypesWithHybridization; //!< Number of atom types with each hybridization
      storage::List< AtomType>   m_AtomTypesByDecreasingStability;   //!< Most stable types first
      size_t                     m_NumberConjugatedTypes;            //!< Number of conjugated types in the list

      //! Function used to resolve the final atom type; only used in a few ambiguous cases
      //! such as trigonal vs. tetrahedral nitrogen with 3 bonds, otherwise NULL
      //! Parameters are
      //! 1. the atom
      //! 2. the size of the smallest ring that this atom is part of
      void ( PossibleAtomTypesForAtom::*m_FinalizeFunction)( const AtomConformationalInterface &, const size_t &);

      //! @brief Create the map from atom environment string to possible atom types
      //! @param IN_AROMATIC_RING whether to only include types that could be in an aromatic ring
      //! @param EXPLICIT_CHARGE true iff when set to 0 the expected charge must be interpreted literally as a neutral atom;
      //! default is false, which preserves backwards compatibility with old behavior where 0 allows searching of all atom types
      static storage::Map< std::string, PossibleAtomTypesForAtom> CreateAtomicEnvironmentToTypesMap
      (
        const bool IN_AROMATIC_RING,
        const bool EXPLICIT_CHARGE = false
      );

    public:

      //! @brief write out the atom typing scheme
      //! @param OSTREAM stream to write the atom typing scheme to
      static std::ostream &WriteDetailedScheme( std::ostream &OSTREAM);

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PossibleAtomTypesForAtom();

      //! @brief constructor from the known information about the atom
      //! @param ELEMENT element type,
      //! @param NUMBER_ELECTRONS_IN_BONDS number of electrons in bonds for the atom type
      //! @param NUMBER_BONDS number of bonds for the atom
      //! @param SUSPECTED_CHARGE; expected charge, ignored if no atom type matching the other criteria if found
      //! @param IN_AROMATIC_RING true iff the atom has bonds of the aromatic unspecified type
      //! @param EXPLICIT_CHARGE true iff when set to 0 the SUSPECTED_CHARGE must be interpreted literally as a neutral atom;
      //! default is false, which preserves backwards compatibility with old behavior where 0 allows searching of all atom types
      PossibleAtomTypesForAtom
      (
        const ElementType &ELEMENT,
        const size_t NUMBER_ELECTRONS_IN_BONDS,
        const size_t NUMBER_BONDS,
        const short SUSPECTED_CHARGE,
        const bool IN_AROMATIC_RING,
        const bool EXPLICIT_CHARGE = false
      );

      //! @brief Clone function
      //! @return pointer to new PossibleAtomTypesForAtom
      PossibleAtomTypesForAtom *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief tell whether a particular hybrid orbital type is possible given what we know about this atom
      //! @param HYBRID the type of hybrid orbital
      //! @return true iff there is a possible atom type for that hybrid orbital
      bool CouldHaveHybridization( const HybridOrbitalType &HYBRID) const;

      //! @brief return the number of types that the atom has the potential to become
      //! @return the number of types that the atom has the potential to become
      size_t GetNumberPossibleTypes() const;

      //! @brief return the most stable type
      //! @return the most stable type
      AtomType GetMostStableType() const;

      //! @brief get the alternate atom types
      //! @return the alternative atom types
      storage::Vector< AtomType> GetAlternateTypes() const;

      //! @brief get the alternate atom type with the given charge
      //! @param CHARGE the charge desired
      //! @return an alternative atom type
      AtomType GetAlternateTypeWithCharge( const short &CHARGE) const;

      //! @brief returns true if any of the candidate types are conjugated
      //! @return true if any of the candidate types are conjugated
      bool CouldBeConjugated() const;

      //! @brief returns true if all of the candidate types are conjugated
      //! @return true if all of the candidate types are conjugated
      bool MustBeConjugated() const;

      //! @brief determine the maximal # of pi-electrons in the pi-electron system
      //! @return the maximal # of pi-electrons in the pi-electron system
      size_t GetMaxElectronsParticipatingInPiSystem() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief add an atom type to be considered
      //! @param ATOM_TYPE the type of atom to consider
      void AddAtomType( const AtomType &ATOM_TYPE);

      //! @brief set this object to only consider the given atom type
      //! @param ATOM_TYPE the atom type desired
      void SetToType( const AtomType &ATOM_TYPE);

      //! @brief set the final type based on the given atom and smallest ring size
      //! @param ATOM the atom of interest
      //! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
      void Finalize( const AtomConformationalInterface &ATOM, const size_t &SMALLEST_RING_SIZE);

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

      //! @brief remove a particular hybrid orbital type from the possible types, unless that would remove all possibilities
      //! @param HYBRID the type of hybrid orbital to remove
      void RemoveHybridization( const HybridOrbitalType &HYBRID);

      //! @brief add an atom type to the search using a set of rules for atom types in aromatic rings
      //! @param ATOM_TYPE the type of atom to consider
      //! @param DESIRED_CHARGE the charge desired
      //! The atom type will be ordered using the distance from the desired charge as the first sort key, second by
      //! the stability.  Unlike AddAtomType, AddAromaticAtomType always adds the type to the considered list
      void AddAromaticAtomType( const AtomType &ATOM_TYPE, const short &DESIRED_CHARGE);

      //! @brief Select the best choice for the atom type wherever possible
      //! @see @link https://structbio.vanderbilt.edu:8443/display/MeilerLab/RethinkingAtomTypeDetection @endlink
      //! @details the link above contains the statistics and models used to select the current set of rules
      void Finalize();

      //! @brief choose the preferred atom type (using VSEPR theory) assuming that the orbitals do not hybridize
      //! @details This is used for elements in group 1, 2, & 7, which do hybridize in the gasteiger scheme
      void FinalizeUnhybridized();

      //! @brief only keep the most stable types for the atom that span the set of desired pi orbital electrons (0-2)
      //! @param DESIRED_CHARGE the preferred charge
      //! used during construction of the maps when there is no part of standardization that
      //! should edit this class
      void FinalizeAromatic( const short &DESIRED_CHARGE);

      //! @brief choose the final atom type for Nitrogen with two single bonds
      //! @param ATOM the atom of interest
      //! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
      void FinalizeNitrogenTwoSingle( const AtomConformationalInterface &ATOM, const size_t &SMALLEST_RING_SIZE);

      //! @brief choose the final atom type for a nitrogen with a single and a double bond
      //! @param ATOM the atom of interest
      //! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
      void FinalizeNitrogenSingleDouble( const AtomConformationalInterface &ATOM, const size_t &SMALLEST_RING_SIZE);

      //! @brief choose the final atom type for Nitrogen with three single bonds
      //! @param ATOM the atom of interest
      //! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
      void FinalizeNitrogenThreeSingle( const AtomConformationalInterface &ATOM, const size_t &SMALLEST_RING_SIZE);

      //! @brief choose the final atom type for Oxygen with two single bonds
      //! @param ATOM the atom of interest
      //! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
      void FinalizeOxygenTwoSingle( const AtomConformationalInterface &ATOM, const size_t &SMALLEST_RING_SIZE);

      //! @brief choose the final atom type for Oxygen with a single and a double bond
      //! @param ATOM the atom of interest
      //! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
      void FinalizeOxygenSingleDouble( const AtomConformationalInterface &ATOM, const size_t &SMALLEST_RING_SIZE);

      //! @brief choose the final atom type for Oxygen with three single bonds
      //! @param ATOM the atom of interest
      //! @param SMALLEST_RING_SIZE the size of the smallest ring this atom is part of
      void FinalizeOxygenThreeSingle( const AtomConformationalInterface &ATOM, const size_t &SMALLEST_RING_SIZE);

      //! @brief get connected element types
      //! @param ATOM the atom of interest
      //! @return a set of the connected element types
      static storage::Set< ElementType> GetConnectedElementTypes( const AtomConformationalInterface &ATOM);

      //! @brief test whether a particular atom is unsaturated, without relying on atom types having already been set
      //! @param ATOM the atom of interest
      //! @return true if atom has no A. unsaturated bonds or B. is part of an aromatic ring or C. has empty orbitals
      static bool IsUnsaturated( const AtomConformationalInterface &ATOM);

      //! @brief count unsaturated neighbors
      //! @param ATOM the atom of interest
      //! @return the number of unsaturated neighbors around ATOM
      static size_t CountUnsaturatedNeighbors( const AtomConformationalInterface &ATOM);

      //! @brief test whether atom is bonded to any halogens
      //! @param ATOM the atom of interest
      //! @return true if the atom is bonded to any halogens
      static bool IsBondedToAHalogen( const AtomConformationalInterface &ATOM);

    }; // class PossibleAtomTypesForAtom

  } // namespace chemistry
} // namespace bcl

#endif //BCL_CHEMISTRY_POSSIBLE_ATOM_TYPES_FOR_ATOM_H_

