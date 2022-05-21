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
#include "chemistry/bcl_chemistry_fragment_conformation_shared.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> FragmentConformationShared::s_Instance
    (
      GetObjectInstances().AddInstance( new FragmentConformationShared())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    FragmentConformationShared::FragmentConformationShared() :
      m_Configuration( new FragmentConfigurationShared)
    {
    }

    //! @brief constructor given atoms with conformation and molecule configuration
    //! @params CONFIGURATION molecule configuration
    //! @params ATOMS_CONFORMATIONAL atom conformation objects
    FragmentConformationShared::FragmentConformationShared
    (
      const util::ShPtr< FragmentConfigurationShared> &CONFIGURATION,
      const AtomVector< AtomConformationalShared> &ATOMS_CONFORMATIONAL
    ) :
      HasProperties< ConformationInterface>( *CONFIGURATION),
      m_Configuration( CONFIGURATION),
      m_Atoms( ATOMS_CONFORMATIONAL),
      m_Name( CONFIGURATION->GetName())
    {
      m_Atoms.LinkToLayer( m_Configuration->GetAtomsIterator());
      CacheNumeric( GetStoredProperties().GetMDLProperties());
    }

    //! @brief constructor with a conformation interface
    //! @params CONFIGURATION configuration interface
    FragmentConformationShared::FragmentConformationShared( const ConfigurationInterface &CONFIGURATION) :
      HasProperties< ConformationInterface>( CONFIGURATION),
      m_Configuration( new FragmentConfigurationShared( CONFIGURATION)),
      m_Name( CONFIGURATION.GetName())
    {
      m_Atoms = AtomVector< AtomConformationalShared>( CONFIGURATION.GetAtomInfo(), CONFIGURATION.GetBondInfo());
      m_Atoms.LinkToLayer( m_Configuration->GetAtomsIterator());
      CacheNumeric( GetStoredProperties().GetMDLProperties());
    }

    //! @brief constructor with a conformation interface
    //! @params CONFORMATION conformation interface
    FragmentConformationShared::FragmentConformationShared( const ConformationInterface &CONFORMATION) :
      HasProperties< ConformationInterface>( CONFORMATION),
      m_Configuration( new FragmentConfigurationShared( CONFORMATION)),
      m_Name( CONFORMATION.GetName())
    {
      m_Atoms = AtomVector< AtomConformationalShared>( CONFORMATION.GetAtomInfo(), CONFORMATION.GetBondInfo());
      m_Atoms.LinkToLayer( m_Configuration->GetAtomsIterator());
      // no need to copy stored properties since the conformation interface's cache was already copied
    }

    //! @brief copy constructor
    //! @params CONFORMATION the parent conformation
    FragmentConformationShared::FragmentConformationShared( const FragmentConformationShared &CONFORMATION) :
      HasProperties< ConformationInterface>( CONFORMATION),
      m_Configuration( CONFORMATION.m_Configuration),
      m_Atoms( CONFORMATION.m_Atoms),
      m_Name( CONFORMATION.m_Name)
    {
      // link atoms to the configurational layer
      m_Atoms.LinkToLayer( m_Configuration->GetAtomsIterator());
      // no need to copy stored properties since the conformation interface's cache was already copied
    }

    //! @brief Clone function
    //! @return pointer to new FragmentConformationShared
    FragmentConformationShared *FragmentConformationShared::Clone() const
    {
      return new FragmentConformationShared( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &FragmentConformationShared::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return an iterator that iterates through atom Conformation.
    //! @return an iterator of atoms with conformation
    iterate::Generic< const AtomConformationalInterface> FragmentConformationShared::GetAtomsIterator() const
    {
      return iterate::Generic< const AtomConformationalInterface>( m_Atoms.Begin(), m_Atoms.End());
    }

    //! @brief return the number of atoms
    //! @return the number of atoms
    const size_t &FragmentConformationShared::GetNumberAtoms() const
    {
      return m_Atoms.GetSize();
    }

    //! @brief return the number of bonds
    //! @return the number of bonds
    const size_t &FragmentConformationShared::GetNumberBonds() const
    {
      return m_Atoms.GetNumberBonds();
    }

    //! @brief return constitution interface-derived object
    //! @return a reference to the associated constitution interface-derived object
    const ConstitutionInterface &FragmentConformationShared::GetConstitution() const
    {
      return m_Configuration->GetConstitution();
    }

    //! @brief return configuration interface-derived object
    //! @return a reference to the associated configuration interface-derived object
    const ConfigurationInterface &FragmentConformationShared::GetConfiguration() const
    {
      return *m_Configuration;
    }

    //! @brief return the adjacency list
    //! @param BOND_SCHEME how to represent bond data as a size_t
    //! @return the adjacency list
    storage::Vector< graph::UndirectedEdge< size_t> >
      FragmentConformationShared::GetAdjacencyList( const ConfigurationalBondTypeData::Data &BOND_SCHEME) const
    {
      return m_Atoms.GetAdjacencyList( BOND_SCHEME);
    }

    //! @brief get info on all atoms in the molecule
    //! @return all the atom info about the molecule
    storage::Vector< sdf::AtomInfo> FragmentConformationShared::GetAtomInfo() const
    {
      return m_Atoms.GetAtomInfo();
    }

    //! @brief get bonds of molecule
    //! @return all the connectivities of a molecule
    storage::Vector< sdf::BondInfo> FragmentConformationShared::GetBondInfo() const
    {
      return m_Atoms.GetBondInfo();
    }

    //! @brief get the index of a particular atom conformational interface
    //! @param ATOM the atom of interest
    //! @return the index of a particular atom conformational interface (undefined if atom is not in this molecule)
    size_t FragmentConformationShared::GetAtomIndex( const AtomConformationalInterface &ATOM) const
    {
      return m_Atoms.GetAtomIndex( ATOM);
    }

    //! @brief get name of small molecule
    const std::string &FragmentConformationShared::GetName() const
    {
      return m_Name;
    }

    //! @brief set name of small molecule
    void FragmentConformationShared::SetName( const std::string &NAME)
    {
      m_Name = NAME;
    }

    //! @brief get the conformational properties merged with the constitutional and configurational properties
    //! if a properties are in multiple layers, the conformation's value will be used
    SmallMoleculeMiscProperties FragmentConformationShared::GetMergedProperties() const
    {
      return SmallMoleculeMiscProperties( GetStoredProperties(), m_Configuration->GetMergedProperties());
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentConformationShared::Read( std::istream &ISTREAM)
    {

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &FragmentConformationShared::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get changable atoms conformational interface
    //! @return iterator to changable atoms conformational interface, which allows this class to call SetPosition
    iterate::Generic< AtomConformationalInterface> FragmentConformationShared::GetAtomsIteratorNonConst()
    {
      return iterate::Generic< AtomConformationalInterface>( m_Atoms.Begin(), m_Atoms.End());
    }

  } // namespace chemistry
} // namespace bcl
