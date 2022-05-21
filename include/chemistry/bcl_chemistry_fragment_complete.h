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

#ifndef BCL_CHEMISTRY_FRAGMENT_COMPLETE_H_
#define BCL_CHEMISTRY_FRAGMENT_COMPLETE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_complete.h"
#include "bcl_chemistry_atom_vector.h"
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_has_properties.h"
#include "sdf/bcl_sdf_molecule_reading_pref.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentComplete
    //! @brief Class that contains molecular configuration data
    //! @details Models stereochemistry, isomeric fragments, chemical adjacency, and aromatic and ring structures
    //!
    //! @see @link example_chemistry_fragment_complete.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Jan 11, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentComplete :
      public HasProperties< ConformationInterface>
    {

    protected:

    //////////
    // data //
    //////////

      AtomVector< AtomComplete> m_Atoms; //!< vector of atoms with configuration
      std::string               m_Name;  //!< Name of this constitution

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      FragmentComplete();

      //! @brief constructor from members
      //! @param ATOMS atom vector of atom completes
      //! @param NAME name of the fragment
      FragmentComplete
      (
        const AtomVector< AtomComplete> &ATOMS,
        const std::string &NAME
      );

      //! @brief constructor from members
      //! @param ATOMS atom vector of atom completes
      //! @param NAME name of the fragment
      //! @param STORED_PROPERTIES stored properties of the fragment
      FragmentComplete
      (
        const AtomVector< AtomComplete> &ATOMS,
        const std::string &NAME,
        const storage::Map< std::string, std::string> &STORED_PROPERTIES
      );

      //! @brief constructor with a constitution interface
      //! @params CONSTITUTION constitution interface
      FragmentComplete( const ConstitutionInterface &CONSTITUTION);

      //! @brief constructor with a configuration interface
      //! @params CONFIGURATION configuration interface
      FragmentComplete( const ConfigurationInterface &CONFIGURATION);

      //! @brief constructor with a conformation interface
      //! @params CONFORMATION conformation interface
      FragmentComplete( const ConformationInterface &CONFORMATION);

      //! @brief Clone function
      //! @return pointer to new FragmentComplete
      FragmentComplete *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns an iterator of atom conformations
      //! @return a generic iterator for atom conformations
      iterate::Generic< const AtomConformationalInterface> GetAtomsIterator() const;

      //! @brief get name of small molecule
      const std::string &GetName() const;

      //! @brief set name of small molecule
      void SetName( const std::string &NAME);

    ////////////////
    // operations //
    ////////////////

      //! @brief return the adjacency list
      //! @param BOND_SCHEME how to represent bond data as a size_t
      //! @return the adjacency list
      storage::Vector< graph::UndirectedEdge< size_t> >
        GetAdjacencyList( const ConfigurationalBondTypeData::Data &BOND_SCHEME) const;

      //! @brief return the number of atoms
      //! @return the number of atoms
      const size_t &GetNumberAtoms() const;

      //! @brief return the number of bonds
      //! @return the number of bonds
      const size_t &GetNumberBonds() const;

      //! @brief get atom info for the molecule
      //! @return atom info for the molecule
      storage::Vector< sdf::AtomInfo> GetAtomInfo() const;

      //! @brief get bonds of molecule
      //! @return all the connectivities of a molecule
      storage::Vector< sdf::BondInfo> GetBondInfo() const;

      //! @brief get the index of a particular atom conformational interface
      //! @param ATOM the atom of interest
      //! @return the index of a particular atom conformational interface (undefined if atom is not in this molecule)
      size_t GetAtomIndex( const AtomConformationalInterface &ATOM) const;

      //! @brief get the number of hydrogen atoms that preceed (i.e. occur sequentially prior to) the given atom
      //! @param ATOM the atom of interest
      //! @return the number of preceeding hydrogen atoms
      size_t GetNumberPreceedingHydrogenAtoms( const AtomConformationalInterface &ATOM) const;

      //! @brief get the number of hydrogen atoms that preceed (i.e. occur sequentially prior to) the given atom
      //! @param ATOM the atom index of interest
      //! @return the number of preceeding hydrogen atoms
      size_t GetNumberPreceedingHydrogenAtoms( const size_t ATOM) const;

      //! @brief test whether this fragment has valences
      bool HasValences() const
      {
        return GetNumberValences() > size_t( 0);
      }

      //! @brief get the conformational properties merged with the constitutional and configurational properties
      //! if a properties are in multiple layers, the conformation's value will be used
      SmallMoleculeMiscProperties GetMergedProperties() const
      {
        return GetStoredProperties();
      }

      const AtomVector< AtomComplete> &GetAtomVector() const
      {
        return m_Atoms;
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief Remove all hydrogens from the molecule
      void RemoveH();

      //! @brief Saturate all single bond valences with H
      void SaturateWithH();

      //! @brief Update all hydrogen positions
      void UpdateH();

      //! @brief Canonicalize the order of atoms in the molecule, such that the atoms are arranged in descending CIP priority
      void Canonicalize();

      //! @brief Encode resonance of terminal atoms, such that e.g. carboxylic acid groups have the same atom type
      //!        Unfortunately there isn't a special gasteiger resonant O/S/etc type. This function just attempts to make
      //!        terminal species that exist in a resonance state equivalent, primarily for comparison using graph-based
      //!        methods
      void EncodeTerminalResonance();

      //! @brief Standardize non-ring bond lengths according to covalent radii
      //! @param MOBILE_ATOMS atoms that can be perturbed during standardization
      void StandardizeBondLengths( const storage::Vector< size_t> &MOBILE_ATOMS = storage::Vector< size_t>());

      //! @brief set coordinates of the molecule
      //! @param COORDINATES vector of coordinates of atoms( indices) of molecule
      //! @param INDICES vector of indices of molecule
      void SetAtomCoordinates
      (
        const storage::Vector< linal::Vector3D> &COORDINATES,
        const storage::Vector< size_t> &INDICES = storage::Vector< size_t>()
      );

      //! @brief try to generate an initial idealized geometry. This will fail miserably for rings but generate something
      //!        reasonable (but potentially clashing) otherwise. Assumes only a single molecule (all atoms reachable from all
      //!        other atoms via covalent bonds, not a complex)
      void IdealizeGeometry();

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

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief merge a given fragment with this fragment
      //! @param FRAGMENT the given fragment
      //! @param BONDS_BETWEEN interconnectivity between this and the given fragment where first member of triplet is
      //!        is index of atom of this fragment while third member is index of atom in the given fragment.
      void MergeWithFragment
      (
        const FragmentComplete &FRAGMENT,
        const storage::Vector< sdf::BondInfo> &BONDS_BETWEEN
      );

      //! @brief get changable atoms conformational interface
      //! @return iterator to changable atoms conformational interface, which allows the base class to call SetPosition
      iterate::Generic< AtomConformationalInterface> GetAtomsIteratorNonConst();

    }; // class FragmentComplete

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_COMPLETE_H_
