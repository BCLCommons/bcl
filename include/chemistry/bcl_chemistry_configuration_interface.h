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

#ifndef BCL_CHEMISTRY_CONFIGURATION_INTERFACE_H_
#define BCL_CHEMISTRY_CONFIGURATION_INTERFACE_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_atom_configurational_interface.h"
#include "bcl_chemistry_atom_configurational_shared.h"
#include "bcl_chemistry_constitution_interface.h"
#include "descriptor/bcl_descriptor_has_cache.h"
#include "descriptor/bcl_descriptor_sequence_interface.h"
#include "graph/bcl_graph_undirected_edge.h"
#include "iterate/bcl_iterate_generic.h"
#include "sdf/bcl_sdf_bond_info.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConfigurationInterface
    //! @brief Interface class for molecule/fragment configuration
    //! @details Interface class for configuration of molecule or fragment. Handles molecular/fragment configuration data.
    //!
    //! @remarks example unnecessary
    //! @author kothiwsk, mendenjl
    //! @date Dec 08, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConfigurationInterface :
      public descriptor::SequenceInterface< AtomConfigurationalInterface>,
      public HasPropertiesInterface< util::ObjectInterface>
    {
    /////////////
    // friends //
    /////////////

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief return an iterator that iterates through atom configuration.
      //! @return an iterator of atom configuration
      virtual iterate::Generic< const AtomConfigurationalInterface> GetAtomsIterator() const = 0;

      //! @brief returns an iterator of atom configuration
      //! @return a generic iterator for atom configuration
      virtual iterate::Generic< const AtomConfigurationalInterface> GetIterator() const
      {
        return GetAtomsIterator();
      }

      //! @brief return the number of atoms
      //! @return the number of atoms
      virtual const size_t &GetNumberAtoms() const = 0;

      //! @brief return the number of bonds
      //! @return the number of bonds
      virtual const size_t &GetNumberBonds() const = 0;

      //! @brief return the adjacency list
      //! @param BOND_SCHEME how to represent bond data as a size_t
      //! @return the adjacency list
      virtual storage::Vector< graph::UndirectedEdge< size_t> >
        GetAdjacencyList( const ConfigurationalBondTypeData::Data &BOND_SCHEME) const = 0;

      //! @brief get bonds of molecule
      //! @return all the connectivities of a molecule
      virtual storage::Vector< sdf::BondInfo> GetBondInfo() const = 0;

      //! @brief get info on all atoms in the molecule
      //! @return all the atom info about the molecule
      virtual storage::Vector< sdf::AtomInfo> GetAtomInfo() const = 0;

      //! @brief get the number of valences for the entire molecule or fragment
      //! @return the number of valences on the entire molecule or fragment
      size_t GetNumberValences() const;

      //! @return a string with each atom type in the molecule
      std::string GetAtomTypesString() const;

      //! @return a string with each bond type in the molecule
      std::string GetBondTypesString() const;

      //! @brief return a string containing the chirality of the atoms of this molecule as a string
      std::string GetChiralityString() const;

      //! @brief get name of small molecule
      virtual const std::string &GetName() const = 0;

      //! @brief set name of small molecule
      virtual void SetName( const std::string &NAME) = 0;

      //! @brief get the configurational properties merged with the constitutional properties
      //! if a properties are in both configuration and constitution, the configuration's value will be used
      virtual SmallMoleculeMiscProperties GetMergedProperties() const = 0;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @return output stream which was written to
      std::ostream &WriteMDL( std::ostream &OSTREAM) const;

    }; // class ConfigurationInterface

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFIGURATION_INTERFACE_H_

