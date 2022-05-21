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

#ifndef BCL_CHEMISTRY_FRAGMENT_CONFIGURATION_SHARED_H_
#define BCL_CHEMISTRY_FRAGMENT_CONFIGURATION_SHARED_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_configurational_shared.h"
#include "bcl_chemistry_atom_vector.h"
#include "bcl_chemistry_configuration_interface.h"
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_constitution_interface.h"
#include "bcl_chemistry_has_properties.h"
#include "bcl_chemistry_molecular_constitution_shared.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentConfigurationShared
    //! @brief Class that contains molecular configuration data
    //! @details Models stereochemistry, isomeric fragments, chemical adjacency, and aromatic and ring structures
    //!
    //! @see @link example_chemistry_fragment_configuration_shared.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Dec 20, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentConfigurationShared :
      public HasProperties< ConfigurationInterface>
    {

    private:

    //////////
    // data //
    //////////

      util::ShPtr< FragmentConstitutionShared>  m_Constitution; //!< pointer to constitution object
      AtomVector< AtomConfigurationalShared>    m_Atoms;        //!< vector of atoms
      std::string                               m_Name;         //!< Name of this molecule or fragment

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
      FragmentConfigurationShared();

      //! @brief constructor given atoms with configuration and molecule constitution
      //! @params CONSTITUTION fragment constitution
      //! @params ATOMS_CONFIGURATIONAL atom configuration objects
      FragmentConfigurationShared
      (
        const util::ShPtr< FragmentConstitutionShared> &CONSTITUTION,
        const AtomVector< AtomConfigurationalShared> &ATOMS_CONFIGURATIONAL
      );

      //! @brief constructor with a configuration interface
      //! @params CONFIGURATION configuration interface
      FragmentConfigurationShared( const ConfigurationInterface &CONFIGURATION);

      //! @brief constructor with a conformation interface
      //! @params CONFORMATION conformation interface
      FragmentConfigurationShared( const ConformationInterface &CONFORMATION);

      //! @brief copy constructor
      //! @params CONFIGURATION the parent configuration
      FragmentConfigurationShared( const FragmentConfigurationShared &CONFORMATION);

      //! @brief Clone function
      //! @return pointer to new FragmentConfigurationShared
      FragmentConfigurationShared *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return an iterator that iterates through atom configuration.
      //! @return an iterator of atom configuration
      iterate::Generic< const AtomConfigurationalInterface> GetAtomsIterator() const;

      //! @brief return the number of atoms
      //! @return the number of atoms
      const size_t &GetNumberAtoms() const;

      //! @brief return the number of bonds
      //! @return the number of bonds
      const size_t &GetNumberBonds() const;

      //! @brief return constitution interface-derived object
      //! @return a reference to the associated constitution interface-derived object
      const ConstitutionInterface &GetConstitution() const;

      //! @brief return the adjacency list
      //! @param BOND_SCHEME how to represent bond data as a size_t
      //! @return the adjacency list
      storage::Vector< graph::UndirectedEdge< size_t> >
        GetAdjacencyList( const ConfigurationalBondTypeData::Data &BOND_SCHEME) const;

      //! @brief get bonds of molecule
      //! @return all the connectivities of a molecule
      storage::Vector< sdf::BondInfo> GetBondInfo() const;

      //! @brief get atom info molecule
      //! @return atom info of a molecule
      storage::Vector< sdf::AtomInfo> GetAtomInfo() const;

      //! @brief get name of small molecule
      const std::string &GetName() const;

      //! @brief set name of small molecule
      void SetName( const std::string &NAME);

      //! @brief get the configurational properties merged with the constitutional properties
      //! if a properties are in both configuration and constitution, the configuration's value will be used
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

    }; // class FragmentConfigurationShared

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_CONFIGURATION_SHARED_H_
