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
#include "chemistry/bcl_chemistry_fragment_constitution_shared.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentConstitutionShared::FragmentConstitutionShared()
    {
    }

    //! @brief constructor given atoms with constitution
    //! @params CONSTITUTION atom constitution objects
    FragmentConstitutionShared::FragmentConstitutionShared
    (
      const AtomVector< AtomConstitutionalShared> &CONSTITUTION
    ) :
      m_Atoms( CONSTITUTION)
    {
    }

    //! @brief constructor with a constitution interface
    //! @params CONSTITUTION constitution interface
    FragmentConstitutionShared::FragmentConstitutionShared( const ConstitutionInterface &CONSTITUTION) :
      HasProperties< ConstitutionInterface>( CONSTITUTION),
      m_Name( CONSTITUTION.GetName())
    {
      // create atoms from initializers and bonds
      m_Atoms = AtomVector< AtomConstitutionalShared>( CONSTITUTION.GetAtomInfo(), CONSTITUTION.GetBondInfo());
    }

    //! @brief constructor with a configuration interface
    //! @params CONFIGURATION configuration interface
    FragmentConstitutionShared::FragmentConstitutionShared( const ConfigurationInterface &CONFIGURATION) :
      HasProperties< ConstitutionInterface>( CONFIGURATION),
      m_Name( CONFIGURATION.GetName())
    {
      // create atoms from initializers and bonds
      m_Atoms = AtomVector< AtomConstitutionalShared>( CONFIGURATION.GetAtomInfo(), CONFIGURATION.GetBondInfo());
    }

    //! @brief constructor with a conformation interface
    //! @params CONFORMATION conformation interface
    FragmentConstitutionShared::FragmentConstitutionShared( const ConformationInterface &CONFORMATION) :
      HasProperties< ConstitutionInterface>( CONFORMATION),
      m_Name( CONFORMATION.GetName())
    {
      // create atoms from initializers and bonds
      m_Atoms = AtomVector< AtomConstitutionalShared>( CONFORMATION.GetAtomInfo(), CONFORMATION.GetBondInfo());
    }

    //! @brief copy constructor
    //! @params CONSTITUTION parent constitution
    FragmentConstitutionShared::FragmentConstitutionShared( const FragmentConstitutionShared &CONSTITUTION) :
      HasProperties< ConstitutionInterface>( CONSTITUTION),
      m_Atoms( CONSTITUTION.m_Atoms),
      m_Name( CONSTITUTION.m_Name)
    {
    }

    //! @brief Clone function
    //! @return pointer to new FragmentConstitutionShared
    FragmentConstitutionShared *FragmentConstitutionShared::Clone() const
    {
      return new FragmentConstitutionShared( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &FragmentConstitutionShared::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return an iterator that iterates through atom constitution.
    //! @return an iterator of atom constitutions
    iterate::Generic< const AtomConstitutionalInterface>
      FragmentConstitutionShared::GetAtomsIterator() const
    {
      return iterate::Generic< const AtomConstitutionalInterface>( m_Atoms.Begin(), m_Atoms.End());
    }

    //! @brief return the number of atoms
    //! @return the number of atoms
    const size_t &FragmentConstitutionShared::GetNumberAtoms() const
    {
      return m_Atoms.GetSize();
    }

    //! @brief return the number of bonds
    //! @return the number of bonds
    const size_t &FragmentConstitutionShared::GetNumberBonds() const
    {
      return m_Atoms.GetNumberBonds();
    }

    //! @brief checks if a given atom is in the molecule
    //! @param ATOM atom of interest
    //! @return bool whether atom is in molecule
    bool FragmentConstitutionShared::IsAtomInMolecule( const AtomConstitutionalInterface &ATOM) const
    {
      return &ATOM >= m_Atoms.Begin() && &ATOM < m_Atoms.End();
    }

    //! @brief return the adjacency list
    //! @param BOND_SCHEME how to represent bond data as a size_t
    //! @return the adjacency list
    storage::Vector< graph::UndirectedEdge< size_t> >
      FragmentConstitutionShared::GetAdjacencyList( const ConstitutionalBondTypeData::Data &BOND_SCHEME) const
    {
      return m_Atoms.GetAdjacencyList( BOND_SCHEME);
    }

    //! @brief get bonds of molecule
    //! @return all the connectivities of a molecule
    storage::Vector< sdf::BondInfo>
      FragmentConstitutionShared::GetBondInfo() const
    {
      return m_Atoms.GetBondInfo();
    }

    //! @brief get atom info molecule
    //! @return atom info of a molecule
    storage::Vector< sdf::AtomInfo> FragmentConstitutionShared::GetAtomInfo() const
    {
      return m_Atoms.GetAtomInfo();
    }

    //! @brief get name of small molecule
    const std::string &FragmentConstitutionShared::GetName() const
    {
      return m_Name;
    }

    //! @brief set name of small molecule
    void FragmentConstitutionShared::SetName( const std::string &NAME)
    {
      m_Name = NAME;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentConstitutionShared::Read( std::istream &ISTREAM)
    {

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &FragmentConstitutionShared::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {

      return OSTREAM;
    }

    // define s_Instance
    const util::SiPtr< const util::ObjectInterface> FragmentConstitutionShared::s_Instance
    (
      GetObjectInstances().AddInstance( new FragmentConstitutionShared())
    );

  } // namespace chemistry
} // namespace bcl
