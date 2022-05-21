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
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> FragmentConfigurationShared::s_Instance
    (
      GetObjectInstances().AddInstance( new FragmentConfigurationShared())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    FragmentConfigurationShared::FragmentConfigurationShared() :
      m_Constitution( new FragmentConstitutionShared)
    {
    }

    //! @brief constructor given atoms with configuration and molecule constitution
    //! @params CONSTITUTION fragment constitution
    //! @params ATOMS_CONFIGURATIONAL atom configuration objects
    FragmentConfigurationShared::FragmentConfigurationShared
    (
      const util::ShPtr< FragmentConstitutionShared> &CONSTITUTION,
      const AtomVector< AtomConfigurationalShared> &ATOMS_CONFIGURATIONAL
    ) :
      m_Constitution( CONSTITUTION),
      m_Atoms( ATOMS_CONFIGURATIONAL)
    {
      m_Atoms.LinkToLayer( m_Constitution->GetAtomsIterator());
    }

    //! @brief constructor with a configuration interface
    //! @params CONFIGURATION configuration interface
    FragmentConfigurationShared::FragmentConfigurationShared( const ConfigurationInterface &CONFIGURATION) :
      HasProperties< ConfigurationInterface>( CONFIGURATION),
      m_Constitution( new FragmentConstitutionShared( CONFIGURATION)),
      m_Name( CONFIGURATION.GetName())
    {
      // create atoms from
      m_Atoms = AtomVector< AtomConfigurationalShared>( CONFIGURATION.GetAtomInfo(), CONFIGURATION.GetBondInfo());
      m_Atoms.LinkToLayer( m_Constitution->GetAtomsIterator());
    }

    //! @brief constructor with a conformation interface
    //! @params CONFORMATION conformation interface
    FragmentConfigurationShared::FragmentConfigurationShared( const ConformationInterface &CONFORMATION) :
      HasProperties< ConfigurationInterface>( CONFORMATION),
      m_Constitution( new FragmentConstitutionShared( CONFORMATION)),
      m_Name( CONFORMATION.GetName())
    {
      // create atoms from
      m_Atoms = AtomVector< AtomConfigurationalShared>( CONFORMATION.GetAtomInfo(), CONFORMATION.GetBondInfo());
      m_Atoms.LinkToLayer( m_Constitution->GetAtomsIterator());
      CacheNumeric( GetStoredProperties().GetMDLProperties());
    }

    //! @brief copy constructor
    //! @params CONFIGURATION the parent configuration
    FragmentConfigurationShared::FragmentConfigurationShared( const FragmentConfigurationShared &CONFIGURATION) :
      HasProperties< ConfigurationInterface>( CONFIGURATION),
      m_Constitution( CONFIGURATION.m_Constitution),
      m_Atoms( CONFIGURATION.m_Atoms),
      m_Name( CONFIGURATION.m_Name)
    {
      // link atoms to the constitutional layer
      m_Atoms.LinkToLayer( m_Constitution->GetAtomsIterator());
      CacheNumeric( GetStoredProperties().GetMDLProperties());
    }

    //! @brief Clone function
    //! @return pointer to new FragmentConfigurationShared
    FragmentConfigurationShared *FragmentConfigurationShared::Clone() const
    {
      return new FragmentConfigurationShared( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &FragmentConfigurationShared::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return an iterator that iterates through atom configuration.
    //! @return an iterator of atom configuration
    iterate::Generic< const AtomConfigurationalInterface> FragmentConfigurationShared::GetAtomsIterator() const
    {
      return iterate::Generic< const AtomConfigurationalInterface>( m_Atoms.Begin(), m_Atoms.End());
    }

    //! @brief return the number of atoms
    //! @return the number of atoms
    const size_t &FragmentConfigurationShared::GetNumberAtoms() const
    {
      return m_Atoms.GetSize();
    }

    //! @brief return the number of bonds
    //! @return the number of bonds
    const size_t &FragmentConfigurationShared::GetNumberBonds() const
    {
      return m_Atoms.GetNumberBonds();
    }

    //! @brief return constitution interface-derived object
    //! @return a reference to the associated constitution interface-derived object
    const ConstitutionInterface &FragmentConfigurationShared::GetConstitution() const
    {
      return *m_Constitution;
    }

    //! @brief return the adjacency list
    //! @param BOND_SCHEME how to represent bond data as a size_t
    //! @return the adjacency list
    storage::Vector< graph::UndirectedEdge< size_t> >
      FragmentConfigurationShared::GetAdjacencyList( const ConfigurationalBondTypeData::Data &BOND_SCHEME) const
    {
      return m_Atoms.GetAdjacencyList( BOND_SCHEME);
    }

    //! @brief get bonds of molecule
    //! @return all the connectivities of a molecule
    storage::Vector< sdf::BondInfo> FragmentConfigurationShared::GetBondInfo() const
    {
      return m_Atoms.GetBondInfo();
    }

    //! @brief get atom info molecule
    //! @return atom info of a molecule
    storage::Vector< sdf::AtomInfo> FragmentConfigurationShared::GetAtomInfo() const
    {
      return m_Atoms.GetAtomInfo();
    }

    //! @brief get name of small molecule
    const std::string &FragmentConfigurationShared::GetName() const
    {
      return m_Name;
    }

    //! @brief set name of small molecule
    void FragmentConfigurationShared::SetName( const std::string &NAME)
    {
      m_Name = NAME;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentConfigurationShared::Read( std::istream &ISTREAM)
    {

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &FragmentConfigurationShared::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the configurational properties merged with the constitutional properties
    //! if a properties are in both configuration and constitution, the configuration's value will be used
    SmallMoleculeMiscProperties FragmentConfigurationShared::GetMergedProperties() const
    {
      return SmallMoleculeMiscProperties( GetStoredProperties(), GetConstitution().GetStoredProperties());
    }

  } // namespace chemistry
} // namespace bcl
