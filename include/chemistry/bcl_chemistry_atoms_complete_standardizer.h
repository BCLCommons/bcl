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

#ifndef BCL_CHEMISTRY_ATOMS_COMPLETE_STANDARDIZER_H_
#define BCL_CHEMISTRY_ATOMS_COMPLETE_STANDARDIZER_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_complete.h"
#include "bcl_chemistry_atom_vector.h"
#include "bcl_chemistry_possible_atom_types_for_atom.h"
#include "graph/bcl_graph_const_graph.h"
#include "graph/bcl_graph_ring.h"
#include "sdf/bcl_sdf_molecule_reading_pref.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomsCompleteStandardizer
    //! @brief Standardizes FragmentComplete object
    //!
    //! @see @link example_fragment_complete.cpp @endlink
    //! @author mendenjl
    //! @date Jan 30, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API AtomsCompleteStandardizer :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      AtomVector< AtomComplete>         &m_Atoms;                  //!< Atoms in the small molecule to standardize
      const std::string                 &m_ID;                     //!< ID of this molecule
      graph::ConstGraph< size_t, size_t> m_Graph;                  //!< The associated graph
      storage::List< graph::Ring>        m_Rings;                  //!< The edge-cover ring perception of the graph
      storage::List< graph::Ring>        m_AromaticRings;          //!< The aromatic rings that were detected
      storage::List< graph::Ring>        m_ConjugatedRings;        //!< The conjugated, but not aromatic, rings that were detected
      storage::Vector< size_t>           m_NumberRings;            //!< Number of rings this atom is part of
      storage::Vector< size_t>           m_NumberDeclaredAromaticBonds; //!< Number of rings with aromatic bond type that this atom is part of
      storage::Vector< size_t>           m_SmallestRingSize;       //!< Smallest ring size for each atom
      storage::Vector< size_t>           m_IsONDoubleBondedToRing; //!< 1 if the atom is an O or N double bonded to a dissimilarly electronegative atom in the ring
      bool                               m_NeedToRecomputeOtherAtoms; //!< 1 if other atom types may need to be recomputed based on the current state

      //! possible types for each atom
      storage::Vector< PossibleAtomTypesForAtom> m_PossibleAtomTypes;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a small molecule to be standardized
      //! @param ATOMS the atoms to standardize
      //! @param ID identifier for the molecule
      //! @param FORCE_RECALCULATION true if all atom types should be recalculated
      AtomsCompleteStandardizer
      (
        AtomVector< AtomComplete> &ATOMS,
        const std::string &ID,
        const bool &FORCE_RECALCULATION
      );

      //! @brief Clone function
      //! @return pointer to new AtomsCompleteStandardizer
      AtomsCompleteStandardizer *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief set bond type conjugations to conjugated/aromatic/or non-conjugated
      static void SetConjugationOfBondTypes( AtomVector< AtomComplete> &ATOMS);

      //! @brief try to neutralize a given molecule
      static void TryNeutralize( AtomVector< AtomComplete> &ATOMS, const sdf::NeutralizationPref &PREF);

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
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief add ring information to the bond types
      void AddRingInformationToBondTypes();

      //! @brief initialize the standardization
      //! @param FORCE whether the force recalculation of atom types
      //! @return true if all types could be assigned an atom type
      bool Initialize( const bool &FORCE);

      //! @brief split rings by conjugation (non-conjugated, conjugated, aromatic)
      void SplitRingsByAromaticity();

      //! @brief check a ring for aromaticity
      //! @param RING ring which has to be detected for aromaticity
      //! @param IS_CONJUGATED_ATOM_IN_RING  0 for atoms that are not conjugated atoms in rings, 1 otherwise
      //! @return true if ring is aromatic; false otherwise
      bool IsAromatic
      (
        const graph::Ring &RING,
        const storage::Vector< size_t> &IS_CONJUGATED_ATOM_IN_RING
      );

      //! @brief test whether a particular path contains any non-conjugatable atoms
      //! @param RING the path to check
      //! @return true if the path contains any non-conjugatable atoms
      bool RingContainsNonConjugableAtoms( const graph::Ring &RING);

      //! @brief test whether any atoms in a path contain double bonds to atoms outside the conjugated atoms list
      //! @param RING the path whose atoms to check
      //! @param IS_CONJUGATED_ATOM_IN_RING 0 for atoms that are not conjugated atoms in rings
      bool DoubleBondsAreToAtomsInSet( const graph::Ring &RING, const storage::Vector< size_t> &IS_CONJUGATED_ATOM_IN_RING) const;

      //! @brief test whether a ring appears to have purely aromatic bond types
      //! @param RING the path whose atoms to check
      bool RingHasAllAromaticBondTypes( const graph::Ring &RING) const;

      //! @brief refine atom types for a ring with all aromatic declared bonds
      //! @param RING the ring, for which RingHasAllAromaticBondTypes must have returned true
      void RefineAtomTypesForAromaticRing( const graph::Ring &RING);

      //! @brief test whether a ring appears to have purely aromatic bond types
      //! @param RING the path whose atoms to check
      bool CheckAllAtomTypeValid( const graph::Ring &RING);

      //! @brief count the electrons in a pi-system along atoms in RING
      //! @param RING the path to examine
      //! @return the number of electrons in a pi-system along atoms in RING
      size_t CountNominalEInPiSystem( const graph::Ring &RING);

      //! @brief test whether a given ring contains an atom with a lone pair
      //! @param RING the path to examine
      //! @return the number of atoms that can donate a lone pair to the pi system
      size_t CountAtomsWithLonePairForPiSystem( const graph::Ring &RING);

      //! Test whether all the bond types in a ring were already set to being aromatic
      bool TestBondTypesAlreadySetToAromatic( const graph::Ring &RING) const;

      //! Set bond types of all bonds in a ring to include aromatic character
      void SetAromaticBondTypes( const graph::Ring &RING);

      //! @brief set the selected atom types for all atoms in the molecule
      //! @return true if all types were defined
      bool SetAtomTypes();

      //! @brief set bond type conjugations to conjugated/aromatic/or non-conjugated
      void SetConjugationOfBondTypes();

      //! @brief get all atom types matching a given atom considering its bonds
      //! @param ATOM atom of interest
      //! @return PossibleAtomTypesForAtom
      PossibleAtomTypesForAtom GetPossibleTypesForAtom( AtomComplete &ATOM);

      //! @brief using known atom types, determine a possible bond type for bonds that have unknown order
      void DetermineUnknownBondOrders();

      //! @brief remove all bonds to / from group 1 elements if the # of bonds from this atom exceeds 1
      //! In configuration compounds, sometimes group 1 elements form a configuration center; the group 1
      //! atom becomes charged and centered around other atoms with a partial charge
      //! The bcl only supports normal covalent bonds at this point, so these bonds are thrown out
      //! @return reference to ATOMS that was passed in
      static AtomVector< AtomComplete> &RemoveObviousIonicBonds( AtomVector< AtomComplete> &ATOMS);

    }; // class AtomsCompleteStandardizer

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_ATOMS_COMPLETE_STANDARDIZER_H_

