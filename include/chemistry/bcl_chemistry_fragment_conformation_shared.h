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

#ifndef BCL_CHEMISTRY_FRAGMENT_CONFORMATION_SHARED_H_
#define BCL_CHEMISTRY_FRAGMENT_CONFORMATION_SHARED_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_vector.h"
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_fragment_configuration_shared.h"
#include "bcl_chemistry_has_properties.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentConformationShared
    //! @brief Class that contains molecular configuration data
    //! @details Models stereochemistry, isomeric fragments, chemical adjacency, and aromatic and ring structures
    //!
    //! @see @link example_chemistry_fragment_conformation_shared.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Dec 20, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentConformationShared :
      public HasProperties< ConformationInterface>
    {

    private:

    //////////
    // data //
    //////////

      util::ShPtr< FragmentConfigurationShared>  m_Configuration; //!< pointer to molecule configuration object
      AtomVector< AtomConformationalShared>      m_Atoms;         //!< vector of atoms
      std::string                                m_Name;          //!< Name of this atom collection

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
      FragmentConformationShared();

      //! @brief constructor given atoms with conformation and molecule configuration
      //! @params CONFIGURATION molecule configuration
      //! @params ATOMS_CONFORMATIONAL atom conformation objects
      FragmentConformationShared
      (
        const util::ShPtr< FragmentConfigurationShared> &CONFIGURATION,
        const AtomVector< AtomConformationalShared> &ATOMS_CONFORMATIONAL
      );

      //! @brief constructor with a conformation interface
      //! @params CONFIGURATION configuration interface
      FragmentConformationShared( const ConfigurationInterface &CONFIGURATION);

      //! @brief constructor with a conformation interface
      //! @params CONFORMATION conformation interface
      FragmentConformationShared( const ConformationInterface &CONFORMATION);

      //! @brief copy constructor
      //! @params CONFORMATION the parent conformation
      FragmentConformationShared( const FragmentConformationShared &CONFORMATION);

      //! @brief Clone function
      //! @return pointer to new FragmentConformationShared
      FragmentConformationShared *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return an iterator that iterates through atom conformation
      //! @return an iterator of atom conformation
      iterate::Generic< const AtomConformationalInterface> GetAtomsIterator() const;

      //! @brief return the number of atoms
      //! @return the number of atoms
      const size_t &GetNumberAtoms() const;

      //! @brief return the number of bonds
      //! @return the number of bonds
      const size_t &GetNumberBonds() const;

      //! @brief return constitution interface-derived object
      //! @return a reference to the associated constitution interface-derived object
      const ConstitutionInterface &GetConstitution() const;

      //! @brief return configuration interface-derived object
      //! @return a reference to the associated configuration interface-derived object
      const ConfigurationInterface &GetConfiguration() const;

      //! @brief return the adjacency list
      //! @param BOND_SCHEME how to represent bond data as a size_t
      //! @return the adjacency list
      storage::Vector< graph::UndirectedEdge< size_t> >
        GetAdjacencyList( const ConfigurationalBondTypeData::Data &BOND_SCHEME) const;

      //! @brief get info on all atoms in the molecule
      //! @return all the atom info about the molecule
      storage::Vector< sdf::AtomInfo> GetAtomInfo() const;

      //! @brief get bonds of molecule
      //! @return all the connectivities of a molecule
      storage::Vector< sdf::BondInfo> GetBondInfo() const;

      //! @brief get the index of a particular atom conformational interface
      //! @param ATOM the atom of interest
      //! @return the index of a particular atom conformational interface (undefined if atom is not in this molecule)
      size_t GetAtomIndex( const AtomConformationalInterface &ATOM) const;

      //! @brief get name of small molecule
      const std::string &GetName() const;

      //! @brief set name of small molecule
      void SetName( const std::string &NAME);

      //! @brief get the conformational properties merged with the constitutional and configurational properties
      //! if a properties are in multiple layers, the conformation's value will be used
      SmallMoleculeMiscProperties GetMergedProperties() const;

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

      //! @brief get changable atoms conformational interface
      //! @return iterator to changable atoms conformational interface, which allows this class to call SetPosition
      iterate::Generic< AtomConformationalInterface> GetAtomsIteratorNonConst();

    }; // class FragmentConformationShared

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_CONFORMATION_SHARED_H_
