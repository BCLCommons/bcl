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

#ifndef BCL_CHEMISTRY_CONFIGURATION_SET_SAME_CONSTITUTION_H_
#define BCL_CHEMISTRY_CONFIGURATION_SET_SAME_CONSTITUTION_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_fragment_configuration_shared.h"
#include "graph/bcl_graph_const_graph.h"
#include "graph/bcl_graph_undirected_edge.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConfigurationSetSameConstitution
    //! @brief Container class for FragmentConstitutionShared objects
    //!
    //! @see @link example_chemistry_configuration_set_same_constitution.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Mar 07, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConfigurationSetSameConstitution :
      public util::ObjectInterface
    {
    public:

    //////////////
    // typedefs //
    //////////////

      typedef util::ShPtrList< FragmentConfigurationShared>::const_iterator const_iterator;

    private:

    //////////
    // data //
    //////////

      //! base constitution; all configurations given will link with this
      util::ShPtr< FragmentConstitutionShared> m_Constitution;

      //! list of FragmentConfigurationShared
      util::ShPtrList< FragmentConfigurationShared> m_Configurations;

      //! whether the configuration is guaranteed to be unique for the given constitution;
      //! this is determined the first time Insert is called
      bool m_IsUnique;

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class Node
      //! @brief Storage class for information about a molecule in the set
      //!
      //! @remarks example unnecessary
      //! @author mendenjl
      //! @date Oct 24, 2012
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct Node
      {
        const_iterator                     m_Itr;   //!< iterator on m_Configurations
        graph::ConstGraph< size_t, size_t> m_Graph; //!< ConstGraph of the molecule

        //! @brief test whether the graph is cached
        //! @return true if the graph is cached on this node
        bool IsGraphCached() const
        {
          return m_Graph.GetSize() > size_t( 0);
        }

        //! @brief cache the graph for the molecule using the iterator
        void CacheGraph();

      };

      //! map of hash string to itr and graph::ConstGraph
      std::map< std::string, std::list< Node> > m_HashToNodeMap;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, provided only so that Map's can use this class
      ConfigurationSetSameConstitution()
      {
      }

      //! @brief constructor from constitution
      //! @param CONSTITUTION the constitution to use
      ConfigurationSetSameConstitution( const util::ShPtr< FragmentConstitutionShared> &CONSTITUTION);

      //! @brief Clone function
      //! @return pointer to new ConfigurationSetSameConstitution
      ConfigurationSetSameConstitution *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief set the constitution for this configuration
      //! @param CONSTITUTION shared pointer to the new constitution
      void SetConstitution( const util::ShPtr< FragmentConstitutionShared> &CONSTITUTION);

      //! @brief return the sole constitution
      //! @return the sole constitution
      const util::ShPtr< FragmentConstitutionShared> &GetConstitution() const
      {
        return m_Constitution;
      }

      //! @brief return a list containing the configurations
      //! @return a list containing the configurations
      const util::ShPtrList< FragmentConfigurationShared> &GetConfigurations() const
      {
        return m_Configurations;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns the number of configurations of a given constitution
      //! @return the number of configurations of a given constitution
      size_t GetSize() const
      {
        return m_Configurations.GetSize();
      }

      //! @brief return const iterator begin to configurations
      //! @return const iterator begin to configurations
      const_iterator Begin() const
      {
        return m_Configurations.Begin();
      }

      //! @brief return const iterator end to configurations
      //! @return const iterator end to configurations
      const_iterator End() const
      {
        return m_Configurations.End();
      }

      //! @brief find a particular configuration
      //! @param FRAGMENT a configuration to search for
      //! @param ISOMORPHISM isomorphism vector to store isomorphism
      //! @return an iterator to the configuration in this set
      const_iterator Find
      (
        const ConfigurationInterface &FRAGMENT,
        storage::Vector< size_t> &ISOMORPHISM
      ) const;

      //! @brief append a molecule to SetConstitution
      //! @param FRAGMENT fragment constitution shared that needs to be added to set constitution
      //! @param CONSTITUTION shared pointer to the constitution for the given fragment
      //! @param ISOMORPHISM  isomorphism vector to store isomorphism
      std::pair< const_iterator, bool> Insert
      (
        const ConfigurationInterface &FRAGMENT,
        const util::ShPtr< FragmentConstitutionShared> &CONSTITUTION,
        storage::Vector< size_t> &ISOMORPHISM
      );

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

    private:

      //! @brief reorder a given fragment configuration shared according to the given isomoprhism
      //! @param CONFIGURATION the configuration to consider
      //! @param REORDERING the mapping indices, usually from an isomorphism
      //! @return a shptr to a newly generated configuration, which is linked to m_Constitution
      util::ShPtr< FragmentConfigurationShared> Reorder
      (
        const ConfigurationInterface &CONFIGURATION,
        const storage::Vector< size_t> &REORDERING
      ) const;

      //! @brief get the atom types given an iterator
      //! @param ITR iterator to the atom types
      template< typename t_AtomType>
      static storage::Vector< AtomType> GetAtomTypes( iterate::Generic< t_AtomType> ITR)
      {
        storage::Vector< AtomType> atom_types;
        atom_types.AllocateMemory( ITR.GetSize());
        for( ; ITR.NotAtEnd(); ++ITR)
        {
          atom_types.PushBack( ITR->GetAtomType());
        }
        return atom_types;
      }

      //! @brief validate that a given configuration matches *m_Constitution, including atom ordering
      //! @param CONFIGURATION atom vector for the configuration of interest
      //! @param REORDERING any reordering to be applied to the configuration first
      void ValidateConfigurationAlignment
      (
        const ConfigurationInterface &CONFIGURATION,
        const storage::Vector< size_t> &REORDERING
      ) const;

      //! @brief function used internally to search for a given fragment within member data structures
      //! @param FRAGMENT a configuration to search for
      //! @param NODES molecules to compare FRAGMENT against
      //! @param ISOMORPHISM isomorphism vector to store isomorphism if it is calculated
      //! @param GRAPH reference to a const graph, will be set to the graph for FRAGMENT if a graph is created
      //! @param FRAG_PTR pointer to a configuration that will be set if one is created internally
      const_iterator InternalFind
      (
        const ConfigurationInterface &FRAGMENT,
        const std::list< Node> &NODES,
        storage::Vector< size_t> &ISOMORPHISM,
        graph::ConstGraph< size_t, size_t> &GRAPH,
        util::ShPtr< FragmentConfigurationShared> &FRAG_PTR
      ) const;

      //! @brief Compose two isomorphisms
      //! @param CONFIGURATION atom vector for the configuration of interest
      //! @param CONSTITUTIONAL_ISOMORPHISM if present, the constitutional isomorphism vector
      //! @param CONFIGURATIONAL_ISOMORPHISM the configurational isomorphism
      //! @return the composed isomorphism
      static storage::Vector< size_t> ComposeIsomorphisms
      (
        const storage::Vector< size_t> &CONSTITUTIONAL_ISOMORPHISM,
        const storage::Vector< size_t> &CONFIGURATIONAL_ISOMORPHISM
      );

      //! @brief make a hash string of the atom types and bond types (order independent)
      //! @param FRAGMENT a fragment configuration
      //! @return the hash string of atom types and bond types for the given fragment
      static std::string MakeBasicHashString( const ConfigurationInterface &FRAGMENT);

      //! @brief test whether a hash string yields a unique isomer
      //! @param HASH the hash string
      //! @return true if only one isomer is possible for the given hash string
      static bool HashStringYieldsUniqueIsomer( const std::string &HASH);

      //! @brief test whether a configuration has no stereocenters or bond isometry
      //! @param FRAGMENT a fragment configuration
      //! @return true if the configuration could have no stereoisomers
      static bool IsUniqueConfigurationForConstitution( const ConfigurationInterface &FRAGMENT);

    }; // class ConfigurationSetSameConstitution

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFIGURATION_SET_SAME_CONSTITUTION_H_
