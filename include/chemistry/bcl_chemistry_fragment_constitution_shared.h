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

#ifndef BCL_CHEMISTRY_FRAGMENT_CONSTITUTION_SHARED_H_
#define BCL_CHEMISTRY_FRAGMENT_CONSTITUTION_SHARED_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_constitutional_interface.h"
#include "bcl_chemistry_atom_constitutional_shared.h"
#include "bcl_chemistry_atom_vector.h"
#include "bcl_chemistry_conformation_interface.h"
#include "bcl_chemistry_constitution_interface.h"
#include "bcl_chemistry_has_properties.h"
#include "iterate/bcl_iterate_generic.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FragmentConstitutionShared
    //! @brief Models molecular data at the constitutional level
    //! Features basic connectivity and atom information
    //!
    //! @see @link example_chemistry_fragment_constitution_shared.cpp @endlink
    //! @author kothiwsk
    //! @date Dec 20, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FragmentConstitutionShared :
      public HasProperties< ConstitutionInterface>
    {

    private:

    //////////
    // data //
    //////////

      AtomVector< AtomConstitutionalShared> m_Atoms; //!< vector of atoms with constitution
      std::string                           m_Name;  //!< Name of this constitution

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FragmentConstitutionShared();

      //! @brief constructor given atoms with constitution
      //! @params CONSTITUTION atom constitution objects
      FragmentConstitutionShared( const AtomVector< AtomConstitutionalShared> &CONSTITUTION);

      //! @brief constructor with a constitution interface
      //! @params CONSTITUTION constitution interface
      FragmentConstitutionShared( const ConstitutionInterface &CONSTITUTION);

      //! @brief constructor with a configuration interface
      //! @params CONFIGURATION configuration interface
      FragmentConstitutionShared( const ConfigurationInterface &CONFIGURATION);

      //! @brief constructor with a conformation interface
      //! @params CONFORMATION conformation interface
      FragmentConstitutionShared( const ConformationInterface &CONFORMATION);

      //! @brief copy constructor
      //! @params CONSTITUTION parent constitution
      FragmentConstitutionShared( const FragmentConstitutionShared &CONSTITUTION);

      //! @brief Clone function
      //! @return pointer to new FragmentConstitutionShared
      FragmentConstitutionShared *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return an iterator that iterates through atom constitution.
      //! @return an iterator of atom constitutions
      iterate::Generic< const AtomConstitutionalInterface> GetAtomsIterator() const;

      //! @brief return the number of atoms
      //! @return the number of atoms
      const size_t &GetNumberAtoms() const;

      //! @brief return the number of bonds
      //! @return the number of bonds
      const size_t &GetNumberBonds() const;

      //! @brief checks if a given atom is in the molecule
      //! @param ATOM atom of interest
      //! @return bool whether atom is in molecule
      bool IsAtomInMolecule( const AtomConstitutionalInterface &ATOM) const;

      //! @brief return the adjacency list
      //! @param BOND_SCHEME how to represent bond data as a size_t
      //! @return the adjacency list
      storage::Vector< graph::UndirectedEdge< size_t> >
        GetAdjacencyList( const ConstitutionalBondTypeData::Data &BOND_SCHEME) const;

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

    ////////////////
    // operations //
    ////////////////

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

    }; // class FragmentConstitutionShared

  } // namespace chemistry

} // namespace bcl

#endif // BCL_CHEMISTRY_FRAGMENT_CONSTITUTION_SHARED_H_
