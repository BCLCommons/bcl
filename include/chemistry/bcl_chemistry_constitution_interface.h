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

#ifndef BCL_CHEMISTRY_CONSTITUTION_INTERFACE_H_
#define BCL_CHEMISTRY_CONSTITUTION_INTERFACE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_constitutional_shared.h"
#include "bcl_chemistry_bond_constitutional.h"
#include "bcl_chemistry_has_properties_interface.h"
#include "bcl_chemistry_small_molecule_misc_properties.h"
#include "graph/bcl_graph_undirected_edge.h"
#include "iterate/bcl_iterate_generic.h"
#include "sdf/bcl_sdf_bond_info.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConstitutionInterface
    //! @brief Interface class for molecule/fragment constitution
    //! @details Interface class for constitution of molecule or fragment. Handles molecular/fragment constitution data.
    //!
    //! @remarks example unnecessary
    //! @author kothiwsk, mendenjl
    //! @date Dec 08, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConstitutionInterface :
      public HasPropertiesInterface< util::ObjectInterface>
    {

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief return an iterator that iterates through atom constitution.
      //! @return an iterator of atom constitutions
      virtual iterate::Generic< const AtomConstitutionalInterface> GetAtomsIterator() const = 0;

      //! @brief return the number of atoms
      //! @return the number of atoms
      virtual const size_t &GetNumberAtoms() const = 0;

      //! @brief return the number of bonds
      //! @return the number of bonds
      virtual const size_t &GetNumberBonds() const = 0;

      //! @brief return a number of hydrogens
      //! @return the number of hydrogens
      size_t GetNumberHydrogens() const;

      //! @brief get the number of valences for the entire molecule or fragment
      //! @return the number of valences on the entire molecule or fragment
      size_t GetNumberValences() const;

      //! @return a string with each atom type in the molecule
      std::string GetAtomTypesString() const;

      //! @return a string with each bond type in the molecule
      std::string GetBondTypesString() const;

      //! @brief checks if a given atom is in the molecule
      //! @param ATOM atom of interest
      //! @return bool whether atom is in molecule
      virtual bool IsAtomInMolecule( const AtomConstitutionalInterface &ATOM) const = 0;

      //! @brief return the adjacency list
      //! @param BOND_SCHEME how to represent bond data as a size_t
      //! @return the adjacency list
      virtual storage::Vector< graph::UndirectedEdge< size_t> >
        GetAdjacencyList( const ConstitutionalBondTypeData::Data &BOND_SCHEME) const = 0;

      //! @brief get bonds of molecule
      //! @return all the connectivities of a molecule
      virtual storage::Vector< sdf::BondInfo> GetBondInfo() const = 0;

      //! @brief get info on all atoms in the molecule
      //! @return all the atom info about the molecule
      virtual storage::Vector< sdf::AtomInfo> GetAtomInfo() const = 0;

      //! @brief get name of small molecule
      virtual const std::string &GetName() const = 0;

      //! @brief set name of small molecule
      virtual void SetName( const std::string &NAME) = 0;

      //! @brief returns a string with the sum formula
      //! @return a string with the sum formula
      std::string GetSumFormula() const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &WriteMDL( std::ostream &OSTREAM) const;

    }; // class ConstitutionInterface

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONSTITUTION_INTERFACE_H_
