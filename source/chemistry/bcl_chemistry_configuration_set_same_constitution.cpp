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
#include "chemistry/bcl_chemistry_configuration_set_same_constitution.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_configuration_graph_converter.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief cache the graph for the given molecule
    void ConfigurationSetSameConstitution::Node::CacheGraph()
    {
      if( !IsGraphCached())
      {
        static ConfigurationGraphConverter make_graph
        (
          ConfigurationGraphConverter::e_AtomTypeAndChirality,
          ConfigurationalBondTypeData::e_BondOrderAmideWithIsometryOrAromaticWithRingness
        );
        m_Graph = make_graph( **m_Itr);
      }
    }

    //! @brief constructor from constitution
    //! @param CONSTITUTION the constitution to use
    ConfigurationSetSameConstitution::ConfigurationSetSameConstitution
    (
      const util::ShPtr< FragmentConstitutionShared> &CONSTITUTION
    )
    {
      if( CONSTITUTION.IsDefined())
      {
        SetConstitution( CONSTITUTION);
      }
    }

    //! @brief Clone function
    //! @return pointer to new ConfigurationSetSameConstitution
    ConfigurationSetSameConstitution *ConfigurationSetSameConstitution::Clone() const
    {
      return new ConfigurationSetSameConstitution( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ConfigurationSetSameConstitution::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the constitution for this configuration
    //! @param CONSTITUTION shared pointer to the new constitution
    void ConfigurationSetSameConstitution::SetConstitution( const util::ShPtr< FragmentConstitutionShared> &CONSTITUTION)
    {
      BCL_Assert
      (
        m_Configurations.IsEmpty(),
        "Cannot set constitution on " + GetClassIdentifier() + " when configurations have already been given"
      );
      BCL_Assert( CONSTITUTION.IsDefined(), GetClassIdentifier() + " cannot be given an undefined constitution");
      m_Constitution = CONSTITUTION;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ConfigurationSetSameConstitution::Read( std::istream &ISTREAM)
    {
      BCL_Exit( "This class cannot be read; it contains iterators internally", -1);
//      *this = ConfigurationSetSameConstitution( ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ConfigurationSetSameConstitution::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
//      SmallMoleculeFactory::WriteToMDLFile( *this, OSTREAM);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief find a particular configuration
    //! @param FRAGMENT a configuration to search for
    //! @param ISOMORPHISM isomorphism vector to store isomorphism
    //! @return an iterator to the configuration in this set
    ConfigurationSetSameConstitution::const_iterator
      ConfigurationSetSameConstitution::Find
      (
        const ConfigurationInterface &FRAGMENT,
        storage::Vector< size_t> &ISOMORPHISM
      ) const
    {
      // handle the case that this is the first configuration seen
      if( m_Configurations.IsEmpty())
      {
        return End();
      }

      // ensure that the right sized molecule was passed in
      ValidateConfigurationAlignment( FRAGMENT, ISOMORPHISM);

      // if there is only one configuration possible for the constitution, return begin
      if( m_IsUnique)
      {
        return Begin();
      }

      // create the hash string
      const std::string hash( MakeBasicHashString( FRAGMENT));

      // look for the nodes
      std::map< std::string, std::list< Node> >::const_iterator itr_node_map( m_HashToNodeMap.find( hash));
      if( itr_node_map == m_HashToNodeMap.end())
      {
        return End();
      }

      const std::list< Node> &nodes( itr_node_map->second);

      // next, check whether the hash must yield a unique isomer
      if( !nodes.empty() && HashStringYieldsUniqueIsomer( hash))
      {
        // unique isomer for this hash; return the front element
        return nodes.front().m_Itr;
      }

      util::ShPtr< FragmentConfigurationShared> reordered_mol;
      graph::ConstGraph< size_t, size_t> fragment_graph;

      // call the internal find command
      return InternalFind( FRAGMENT, nodes, ISOMORPHISM, fragment_graph, reordered_mol);
    }

    //! @brief append a molecule to SetConstitution
    //! @param FRAGMENT fragment constitution shared that needs to be added to set constitution
    //! @param CONSTITUTION shared pointer to the constitution for the given fragment
    //! @param ISOMORPHISM isomorphism which if it was determined at constitution layer
    //! @return a pair containing iterator to configuration layer of FRAGMENT and whether configuration has been seen
    //!         earlier (false) or seen for the first time (true)
    std::pair< ConfigurationSetSameConstitution::const_iterator, bool>
      ConfigurationSetSameConstitution::Insert
      (
        const ConfigurationInterface &FRAGMENT,
        const util::ShPtr< FragmentConstitutionShared> &CONSTITUTION,
        storage::Vector< size_t> &ISOMORPHISM
      )
    {
      if( !m_Constitution.IsDefined())
      {
        BCL_Assert
        (
          ISOMORPHISM.IsEmpty(),
          "Insert should only be given an isomorphism if it already has a constitution"
        );
        SetConstitution( CONSTITUTION);
      }
      else
      {
        // validate the isomorphism
        ValidateConfigurationAlignment( FRAGMENT, ISOMORPHISM);
      }

      // handle the case that this is the first configuration seen
      if( m_Configurations.IsEmpty())
      {
        // test whether the configuration was trivial/unique for the given constitution
        m_IsUnique = IsUniqueConfigurationForConstitution( FRAGMENT);
        m_Configurations.PushBack( Reorder( FRAGMENT, ISOMORPHISM));

        if( !m_IsUnique)
        {
          // non-unique constitution, use hash strings for later lookups of configurations
          Node new_node;
          new_node.m_Itr = m_Configurations.Last();
          m_HashToNodeMap[ MakeBasicHashString( FRAGMENT)].push_back( new_node);
        }
        return std::make_pair( m_Configurations.Last(), true);
      }

      // handle the case that the configuration is unique given the constitution
      if( m_IsUnique)
      {
        // if only one configuration is possible for this constitution (a common case) just return it
        return std::make_pair( m_Configurations.Begin(), false);
      }

      // make hash string for the configuration
      const std::string hash_fragment( MakeBasicHashString( FRAGMENT));

      // non-unique constitution, use hash strings for later lookups of configurations
      std::list< Node> &nodes_for_hash( m_HashToNodeMap[ hash_fragment]);

      // if configuration hash string was seen before and it specifies a unique isomer
      if( !nodes_for_hash.empty())
      {
        if( HashStringYieldsUniqueIsomer( hash_fragment))
        {
          return std::make_pair( nodes_for_hash.front().m_Itr, false);
        }
        else
        {
          // cache the first graph after the first time that it was seen
          nodes_for_hash.front().CacheGraph();
        }
      }

      util::ShPtr< FragmentConfigurationShared> reordered_mol;
      graph::ConstGraph< size_t, size_t> fragment_graph;

      // look for the new fragment
      const_iterator itr( InternalFind( FRAGMENT, nodes_for_hash, ISOMORPHISM, fragment_graph, reordered_mol));

      // check whether the fragment was found
      if( itr != End())
      {
        // yep, was found, so return it
        return std::make_pair( itr, false);
      }

      // new fragment
      if( !reordered_mol.IsDefined())
      {
        // reordering was not necessary for Find, but we need a real molecule to insert
        reordered_mol = Reorder( FRAGMENT, ISOMORPHISM);
      }

      // add the molecule
      m_Configurations.PushBack( reordered_mol);

      // add the new node
      nodes_for_hash.push_back( Node());
      nodes_for_hash.back().m_Itr = m_Configurations.Last();
      nodes_for_hash.back().m_Graph = fragment_graph;
      nodes_for_hash.back().CacheGraph();

      return std::make_pair( m_Configurations.Last(), true);
    }

    //! @brief reorder a given fragment configuration shared according to the given isomoprhism
    //! @param CONFIGURATION the configuration to consider
    //! @param REORDERING the mapping indices, usually from an isomorphism
    //! @return a shptr to a newly generated configuration, which is linked to m_Constitution
    util::ShPtr< FragmentConfigurationShared> ConfigurationSetSameConstitution::Reorder
    (
      const ConfigurationInterface &CONFIGURATION,
      const storage::Vector< size_t> &REORDERING
    ) const
    {
      AtomVector< AtomConfigurationalShared> atoms_reordered( CONFIGURATION.GetAtomInfo(), CONFIGURATION.GetBondInfo());
      if( REORDERING.GetSize() == CONFIGURATION.GetNumberAtoms())
      {
        atoms_reordered.Reorder( REORDERING);
      }
      return
        util::ShPtr< FragmentConfigurationShared>
        (
          new FragmentConfigurationShared( m_Constitution, atoms_reordered)
        );
    }

    //! @brief validate that a given configuration matches *m_Constitution, including atom ordering
    //! @param CONFIGURATION atom vector for the configuration of interest
    //! @param REORDERING any reordering to be applied to the configuration first
    void ConfigurationSetSameConstitution::ValidateConfigurationAlignment
    (
      const ConfigurationInterface &CONFIGURATION,
      const storage::Vector< size_t> &REORDERING
    ) const
    {
      BCL_Assert
      (
        REORDERING.IsEmpty()
        || REORDERING.GetSize() == CONFIGURATION.GetNumberAtoms(),
        "Bad isomorphism given to Matches; # atoms: " + util::Format()( CONFIGURATION.GetNumberAtoms())
        + " but isomorphism had size: " + util::Format()( REORDERING.GetSize())
      );

      // check # of atoms and bonds
      BCL_Assert
      (
        CONFIGURATION.GetNumberAtoms() == m_Constitution->GetNumberAtoms()
        && CONFIGURATION.GetNumberBonds() == m_Constitution->GetNumberBonds(),
        "Wrong size of molecule given"
      );

      // check atom types
      {
        storage::Vector< AtomType> atom_types( GetAtomTypes( CONFIGURATION.GetAtomsIterator()));

        // handle reordering
        if( !REORDERING.IsEmpty())
        {
          atom_types.Reorder( REORDERING);
        }

        // check that the atom types are identical
        BCL_Assert( atom_types == GetAtomTypes( m_Constitution->GetAtomsIterator()), "Atom types did not match");
      }

      // check bond types
      storage::Vector< graph::UndirectedEdge< size_t> > bond_info
      (
        CONFIGURATION.GetAdjacencyList( ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness)
      );

      // sort bond info if necessary due to non-trivial inverse ordering
      if( !REORDERING.IsEmpty())
      {
        // handle reordering
        const size_t n_atoms( CONFIGURATION.GetNumberAtoms());
        storage::Vector< size_t> new_order_inverse( n_atoms, util::GetUndefined< size_t>());
        for( size_t i( 0); i < n_atoms; ++i)
        {
          BCL_Assert( REORDERING( i) < n_atoms, "Bad order vector (index out of range)");
          BCL_Assert
          (
            !util::IsDefined( new_order_inverse( REORDERING( i))),
            "Bad order vector (indices non-unique)"
          );
          new_order_inverse( REORDERING( i)) = i;
        }
        for
        (
          storage::Vector< graph::UndirectedEdge< size_t> >::iterator
            itr( bond_info.Begin()), itr_end( bond_info.End());
          itr != itr_end;
          ++itr
        )
        {
          // update indices too, if reorder was performed
          *itr =
            graph::UndirectedEdge< size_t>
            (
              new_order_inverse( itr->GetIndexLow()),
              new_order_inverse( itr->GetIndexHigh()),
              itr->GetEdgeData()
            );
        }
        bond_info.Sort( std::less< graph::UndirectedEdge< size_t> >());
      }
      const storage::Vector< graph::UndirectedEdge< size_t> > constitutional_bond_info
      (
        m_Constitution->GetAdjacencyList( ConstitutionalBondTypeData::e_BondOrderAmideOrAromaticWithRingness)
      );
      BCL_Assert
      (
        bond_info == constitutional_bond_info,
        "Bond types did not match " + util::Format()( bond_info) + " " + util::Format()( constitutional_bond_info)
      );
    }

    //! @brief function used internally to search for a given fragment within member data structures
    //! @param FRAGMENT a configuration to search for
    //! @param NODES molecules to compare FRAGMENT against
    //! @param ISOMORPHISM isomorphism vector to store isomorphism if it is calculated
    //! @param GRAPH reference to a const graph, will be set to the graph for FRAGMENT if a graph is created
    //! @param FRAG_PTR pointer to a configuration that will be set if one is created internally
    ConfigurationSetSameConstitution::const_iterator ConfigurationSetSameConstitution::InternalFind
    (
      const ConfigurationInterface &FRAGMENT,
      const std::list< Node> &NODES,
      storage::Vector< size_t> &ISOMORPHISM,
      graph::ConstGraph< size_t, size_t> &GRAPH,
      util::ShPtr< FragmentConfigurationShared> &FRAG_PTR
    ) const
    {
      // handle the simple case where no configurations are present
      if( NODES.empty())
      {
        return End();
      }

      static ConfigurationGraphConverter make_graph
      (
        ConfigurationGraphConverter::e_AtomTypeAndChirality,
        ConfigurationalBondTypeData::e_BondOrderAmideWithIsometryOrAromaticWithRingness
      );

      // create the reordered molecule, if the isomorphism was given
      if( !ISOMORPHISM.IsEmpty())
      {
        FRAG_PTR = Reorder( FRAGMENT, ISOMORPHISM);
      }
      GRAPH = make_graph( ISOMORPHISM.IsEmpty() ? FRAGMENT : *FRAG_PTR);

      // perform isomorphism matching
      graph::SubgraphIsomorphism< size_t, size_t> isomorphism;
      isomorphism.SetSubgraphExternalOwnership( GRAPH);
      for
      (
        std::list< Node>::const_iterator itr_node( NODES.begin()), itr_node_end( NODES.end());
        itr_node != itr_node_end;
        ++itr_node
      )
      {
        if( itr_node->IsGraphCached())
        {
          isomorphism.SetGraphExternalOwnership( itr_node->m_Graph);
        }
        else
        {
          // graph was not cached, create it
          isomorphism.SetGraph( make_graph( **itr_node->m_Itr));
        }
        // test whether the molecules really are equal
        if( isomorphism.FindIsomorphism())
        {
          // graphs are equal
          ISOMORPHISM = ComposeIsomorphisms( ISOMORPHISM, isomorphism.GetInverseIsomorphism());
          return itr_node->m_Itr;
        }
      }

      return End();
    }

    //! @brief Compose two isomorphisms
    //! @param CONSTITUTIONAL_ISOMORPHISM if present, the constitutional isomorphism vector
    //! @param CONFIGURATIONAL_ISOMORPHISM the configurational isomorphism
    //! @return the composed isomorphism
    storage::Vector< size_t> ConfigurationSetSameConstitution::ComposeIsomorphisms
    (
      const storage::Vector< size_t> &CONSTITUTIONAL_ISOMORPHISM,
      const storage::Vector< size_t> &CONFIGURATIONAL_ISOMORPHISM
    )
    {
      // if an isomorphism was performed at the constitutional level
      if( CONSTITUTIONAL_ISOMORPHISM.GetSize())
      {
        // update the isomorphism mapping to maintain the invariant that ISOMORPHISM maps a fragment outside to the
        // set's fragment
        const size_t n_atoms( CONFIGURATIONAL_ISOMORPHISM.GetSize());
        storage::Vector< size_t> composition( n_atoms);
        for( size_t i( 0); i < n_atoms; ++i)
        {
          composition( i) = CONSTITUTIONAL_ISOMORPHISM( CONFIGURATIONAL_ISOMORPHISM( i));
        }
        return composition;
      }
      // no constitutional isomorphism given, so no remapping is necessary
      return CONFIGURATIONAL_ISOMORPHISM;
    }

    //! @brief make a hash string of the isometry
    //! @param a fragment constitution
    //! @return the hash string of isometry
    std::string ConfigurationSetSameConstitution::MakeBasicHashString( const ConfigurationInterface &FRAGMENT)
    {
      size_t r_chiral_atoms( 0);
      size_t s_chiral_atoms( 0);
      size_t c_chiral_atoms( 0);
      size_t t_chiral_atoms( 0);
      size_t e_isometric_bonds( 0);
      size_t z_isometric_bonds( 0);
      storage::Map< std::string, std::pair< size_t, size_t> > rs_atom_type_counts;
      storage::Map< std::string, std::pair< size_t, size_t> > ez_atom_type_counts;
      storage::Map< std::string, std::pair< size_t, size_t> > ct_atom_type_counts;

      // iterate through all atoms of the molecule and count the appearance of each confrigurational element
      for
      (
        iterate::Generic< const AtomConfigurationalInterface> itr_atom( FRAGMENT.GetAtomsIterator());
        itr_atom.NotAtEnd();
        ++itr_atom
      )
      {
        // accumulate chiral center counts
        const std::string &atom_type_name( itr_atom->GetAtomType().GetName());
        if( itr_atom->GetChirality() == e_RChirality)
        {
          ++r_chiral_atoms;
          ++rs_atom_type_counts[ atom_type_name].first;
        }
        else if( itr_atom->GetChirality() == e_SChirality)
        {
          ++s_chiral_atoms;
          ++rs_atom_type_counts[ atom_type_name].second;
        }
        else if( itr_atom->GetChirality() == e_CisRingChirality)
        {
          ++c_chiral_atoms;
          ++ct_atom_type_counts[ atom_type_name].first;     
        }
        else if( itr_atom->GetChirality() == e_TransRingChirality)
        {
          ++t_chiral_atoms;
          ++ct_atom_type_counts[ atom_type_name].second;
        }

        // get all the bonds
        const storage::Vector< BondConfigurational> &bonds( itr_atom->GetBonds());

        // count e-z in bonds
        for
        (
          storage::Vector< BondConfigurational>::const_iterator itr_bond( bonds.Begin()), itr_bond_end( bonds.End());
          itr_bond != itr_bond_end;
          ++itr_bond
        )
        {
          const ConfigurationalBondType &bond_type( itr_bond->GetBondType());
          // skip the bond going one direction to avoid double-counting
          if( &*itr_atom > &itr_bond->GetTargetAtom() || bond_type->GetIsometry() == e_NonIsometric)
          {
            continue;
          }
          std::string key;
          // form the key for this bond using its constitent atom types in sorted order
          if( atom_type_name <= itr_bond->GetTargetAtom().GetAtomType().GetName())
          {
            key = atom_type_name + "-" + itr_bond->GetTargetAtom().GetAtomType().GetName();
          }
          else
          {
            key = itr_bond->GetTargetAtom().GetAtomType().GetName() + '-' + atom_type_name;
          }
          if( bond_type->GetIsometry() == e_EIsometry)
          {
            ++e_isometric_bonds;
            ++ez_atom_type_counts[ key].first;
          }
          else if( bond_type->GetIsometry() == e_ZIsometry)
          {
            ++z_isometric_bonds;
            ++ez_atom_type_counts[ key].second;
          }
        }
      }

      // stringstream for the hash
      std::stringstream hash;

      // count # of options for this hash string
      // This will be prepended to the start of the hash
      size_t permutations( 1);

      // add atom type count information for chiral centers if necessary
      if( r_chiral_atoms && s_chiral_atoms)
      {
        // if there are both R and S chiral centers, output additional information
        if( rs_atom_type_counts.GetSize() == size_t( 1))
        {
          // common case: just one atom type (usually C_TeTeTeTe) yielded stereocenters
          // Thus, no need to break down the stereocenter hash by atom types
          hash << 'R' << r_chiral_atoms << 'S' << s_chiral_atoms;
          // just compute the permutations
          const size_t r( rs_atom_type_counts.Begin()->second.first);
          const size_t s( rs_atom_type_counts.Begin()->second.second);
          permutations *= math::BinomialCoefficient( r + s, r);
        }
        else
        {
          // Several atom types yielded stereocenters
          // Use atom type names to improve the chances of declaring this to be a unique stereocenter
          for
          (
            storage::Map< std::string, std::pair< size_t, size_t> >::const_iterator
              itr( rs_atom_type_counts.Begin()), itr_end( rs_atom_type_counts.End());
            itr != itr_end;
            ++itr
          )
          {
            const size_t r( itr->second.first);
            const size_t s( itr->second.second);
            hash << ' ' << itr->first << '-';
            if( r && s)
            {
              hash << 'R' << r << 'S' << s;
              // if there are both R and Z types for this atom type, then the hash string is impure
              // to get an accurate  multiply by binomial_coefficient( s, r+s)
              permutations *= math::BinomialCoefficient( r + s, r);
            }
            else
            {
              hash << ( r ? 'R' : 'S');
            }
          }
        }
      }
      // handle case of pure R and S chirality
      else if( r_chiral_atoms)
      {
        hash << 'R';
      }
      else if( s_chiral_atoms)
      {
        hash << 'S';
      }

      // add atom type count information for ring chiral centers if necessary
      if( c_chiral_atoms && t_chiral_atoms)
      {
        // if there are both R and S chiral centers, output additional information
        if( ct_atom_type_counts.GetSize() == size_t( 1))
        {
          // common case: just one atom type (usually C_TeTeTeTe) yielded stereocenters
          // Thus, no need to break down the stereocenter hash by atom types
          hash << 'C' << c_chiral_atoms << 'T' << t_chiral_atoms;
          // just compute the permutations
          const size_t c( ct_atom_type_counts.Begin()->second.first);
          const size_t t( ct_atom_type_counts.Begin()->second.second);
          permutations *= math::BinomialCoefficient( c + t, c);
        }
        else
        {
          // Several atom types yielded stereocenters
          // Use atom type names to improve the chances of declaring this to be a unique stereocenter
          for
          (
            storage::Map< std::string, std::pair< size_t, size_t> >::const_iterator
              itr( ct_atom_type_counts.Begin()), itr_end( ct_atom_type_counts.End());
            itr != itr_end;
            ++itr
          )
          {
            const size_t c( itr->second.first);
            const size_t t( itr->second.second);
            hash << ' ' << itr->first << '-';
            if( c && t)
            {
              hash << 'C' << c << 'T' << t;
              // if there are both R and Z types for this atom type, then the hash string is impure
              // to get an accurate  multiply by binomial_coefficient( s, r+s)
              permutations *= math::BinomialCoefficient( c + t, c);
            }
            else
            {
              hash << ( c ? 'C' : 'T');
            }
          }
        }
      }
      // handle case of pure R and S chirality
      else if( c_chiral_atoms)
      {
        hash << 'C';
      }
      else if( t_chiral_atoms)
      {
        hash << 'T';
      } 

      // add atom type count information for bond isomers if necessary
      if( e_isometric_bonds && z_isometric_bonds)
      {
        if( ez_atom_type_counts.GetSize() != size_t( 1))
        {
          for
          (
            storage::Map< std::string, std::pair< size_t, size_t> >::const_iterator
              itr( ez_atom_type_counts.Begin()), itr_end( ez_atom_type_counts.End());
            itr != itr_end;
            ++itr
          )
          {
            const size_t e( itr->second.first);
            const size_t z( itr->second.second);
            hash << ' ' << itr->first << '-';
            if( e && z)
            {
              hash << 'E' << e << 'Z' << z;
              // if there are both E and Z types for this atom type pair, then the hash string is impure
              // to get an accurate  multiply by binomial_coefficient( s, r+s)
              permutations *= math::BinomialCoefficient( e + z, e);
            }
            else
            {
              // all double bond isomers for this atom type are either E or Z, just write out which they are
              hash << ( e ? 'E' : 'Z');
            }
          }
        }
        else
        {
          // just one atom type pair for e-z isomers; no need to clutter the hash with its identity since it would
          // not reduce the # of hash collisions.  Just write out the number of E and Z bonds
          hash << 'E' << e_isometric_bonds << 'Z' << z_isometric_bonds;
          // just compute the permutations
          const size_t e( ez_atom_type_counts.Begin()->second.first);
          const size_t z( ez_atom_type_counts.Begin()->second.second);
          permutations *= math::BinomialCoefficient( e + z, e);
        }
      }
      // handle case of pure R and S chirality
      else if( e_isometric_bonds)
      {
        hash << 'E';
      }
      else if( z_isometric_bonds)
      {
        hash << 'Z';
      }

      // return the isometry hash
      return util::Format()( permutations) + " " + hash.str();
    }

    //! @brief test whether a hash string yields a unique isomer
    //! @param HASH the hash string
    //! @return true if only one isomer is possible for the given hash string
    bool ConfigurationSetSameConstitution::HashStringYieldsUniqueIsomer( const std::string &HASH)
    {
      // a hash string gives rise to a pure configuration if the number of permutations was just 1
      return HASH.size() >= 2 && HASH[ 0] == '1' && HASH[ 1] == ' ';
    }

    //! @brief test whether a configuration has no stereocenters or bond isometry
    //! @param FRAGMENT a fragment configuration
    //! @return true if the configuration could have no stereoisomers
    bool ConfigurationSetSameConstitution::IsUniqueConfigurationForConstitution( const ConfigurationInterface &FRAGMENT)
    {
      // test whether any stereocenters are present
      for
      (
        iterate::Generic< const AtomConfigurationalInterface> itr( FRAGMENT.GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        if( itr->GetChirality() != e_NonChiral)
        {
          return false;
        }

        const storage::Vector< BondConfigurational> &bonds( itr->GetBonds());

        // skip atoms with excessive valence bonds, they could not be isometric
        if( itr->GetAtomType()->GetNumberBonds() != size_t( 3) || bonds.GetSize() <= size_t( 1))
        {
          continue;
        }

        // check for isometry on the bonds that are present
        if
        (
          bonds( 0).GetBondType()->GetIsometry() != e_NonIsometric
          || bonds( 1).GetBondType()->GetIsometry() != e_NonIsometric
          || ( bonds.GetSize() == size_t( 3) && bonds( 2).GetBondType()->GetIsometry() != e_NonIsometric)
        )
        {
          return false;
        }
      }

      // no stereo isometry is possible for this molecule
      return true;
    }

  } // namespace chemistry
} // namespace bcl
