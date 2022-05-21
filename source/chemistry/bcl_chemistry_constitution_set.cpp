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
#include "chemistry/bcl_chemistry_constitution_set.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_constitution_graph_converter.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief construct from iterator of ConstitutionInterface
    //! @param MOLECULES an iterator to configuration of molecules
    ConstitutionSet::ConstitutionSet( iterate::Generic< const ConstitutionInterface> MOLECULES)
    {
      for( ; MOLECULES.NotAtEnd(); ++MOLECULES)
      {
        Insert( *MOLECULES);
      }
    }

    //! @brief Clone function
    //! @return pointer to new ConstitutionSet
    ConstitutionSet *ConstitutionSet::Clone() const
    {
      return new ConstitutionSet( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief cache the graph for the given molecule
    void ConstitutionSet::Node::CacheGraph()
    {
      if( !IsGraphCached())
      {
        static ConstitutionGraphConverter make_graph;
        m_Graph = make_graph( **m_Itr);
      }
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ConstitutionSet::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ConstitutionSet::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ConstitutionSet::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief find a particular constitution
    //! @param FRAGMENT a constitution to search for
    //! @param ISOMORPHISM isomorphism vector to store isomorphism
    //! @return an iterator to the constitution in this set
    ConstitutionSet::const_iterator ConstitutionSet::Find
    (
      const ConstitutionInterface &FRAGMENT,
      util::SiPtr< storage::Vector< size_t> > ISOMORPHISM
    ) const
    {
      // handle the case that this is the first configuration seen
      if( m_Constitutions.IsEmpty())
      {
        return End();
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

      graph::ConstGraph< size_t, size_t> fragment_graph;

      storage::Vector< size_t> inverse_isomorphism;

      if( !ISOMORPHISM.IsDefined())
      {
        ISOMORPHISM = util::ToSiPtrNonConst( inverse_isomorphism);
      }

      // call the internal find command
      return InternalFind( FRAGMENT, nodes, *ISOMORPHISM, fragment_graph);
    }

    //! @brief append a molecule to SetConstitution
    //! @param FRAGMENT fragment constitution shared that needs to be added to set constitution
    //! @param ISOMORPHISM isomorphism vector to store isomorphism
    //! @return a pair containing iterator to constitution layer of FRAGMENT and whether constitution has been seen
    //!         earlier (false) or seen for the first time (true)
    std::pair< ConstitutionSet::const_iterator, bool> ConstitutionSet::Insert
    (
      const ConstitutionInterface &FRAGMENT,
      util::SiPtr< storage::Vector< size_t> > ISOMORPHISM
    )
    {
      std::string hash_fragment( MakeBasicHashString( FRAGMENT));

      storage::Vector< size_t> inverse_isomorphism;

      if( !ISOMORPHISM.IsDefined())
      {
        ISOMORPHISM = util::ToSiPtrNonConst( inverse_isomorphism);
      }

      std::list< Node> &nodes_for_hash( m_HashToNodeMap[ hash_fragment]);

      graph::ConstGraph< size_t, size_t> fragment_graph;

      // if configuration hash string was seen before and it specifies a unique isomer
      if( !nodes_for_hash.empty())
      {
        // cache the first graph after the first time that it was seen
        nodes_for_hash.front().CacheGraph();

        // look for the new fragment
        const_iterator itr( InternalFind( FRAGMENT, nodes_for_hash, *ISOMORPHISM, fragment_graph));

        // check whether the fragment was found
        if( itr != End())
        {
          // yep, was found, so return it
          return std::make_pair( itr, false);
        }
      }

      // add the molecule
      m_Constitutions.PushBack( util::ShPtr< FragmentConstitutionShared>( new FragmentConstitutionShared( FRAGMENT)));

      // add the new node
      nodes_for_hash.push_back( Node());
      nodes_for_hash.back().m_Itr = m_Constitutions.Last();
      nodes_for_hash.back().m_Graph = fragment_graph;

      return std::make_pair( m_Constitutions.Last(), true);
    }

    //! @brief make a hash string of the atom types and bond types
    //! @param a fragment constitution
    //! @return the hash string of atom types and bond types for the given fragment
    std::string ConstitutionSet::MakeBasicHashString( const ConstitutionInterface &FRAGMENT)
    {
      storage::Vector< size_t> atom_types( GetAtomTypes().GetEnumCount() + 1, size_t( 0));
      storage::Vector< size_t> bond_types( 10, size_t( 0));

      // iterate through all atoms of the molecule and count their appearance
      for
      (
        iterate::Generic< const AtomConstitutionalInterface> itr_atoms( FRAGMENT.GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        if( itr_atoms->GetAtomType().IsDefined())
        {
          ++atom_types( itr_atoms->GetAtomType().GetIndex());
        }
        else
        {
          ++atom_types.LastElement();
        }
      }

      storage::Vector< sdf::BondInfo> bond_info( FRAGMENT.GetBondInfo());
      for
      (
        storage::Vector< sdf::BondInfo>::const_iterator
          itr_bonds( bond_info.Begin()), itr_bonds_end( bond_info.End());
        itr_bonds != itr_bonds_end;
        ++itr_bonds
      )
      {
        if( itr_bonds->GetConstitutionalBondType().IsDefined())
        {
          ++bond_types( itr_bonds->GetConstitutionalBondType()->GetBondData( ConstitutionalBondTypeData::e_BondOrderAmideOrAromaticWithRingness));
        }
        else
        {
          ++bond_types.LastElement();
        }
      }

      // stringstream for the final sum formula
      std::stringstream sum_formula;

      // assemble sum formula with delimiter "_"
      for( size_t atom_type_id( 0), n_atom_types( atom_types.GetSize()); atom_type_id < n_atom_types; ++atom_type_id)
      {
        if( atom_types( atom_type_id) != 0)
        {
          sum_formula << atom_type_id;

          // only print the number of atoms if it is greater than 1
          if( atom_types( atom_type_id) > 1)
          {
            sum_formula << 'x' << atom_types( atom_type_id);
          }

          sum_formula << ',';
        }
      }

      // append bond types
      sum_formula << "BT,";
      for( size_t bond_type_id( 0), n_bond_types( bond_types.GetSize()); bond_type_id < n_bond_types; ++bond_type_id)
      {
        if( bond_types( bond_type_id) != 0)
        {
          sum_formula << bond_type_id;

          // only print the number of bonds if it is greater than 1
          if( bond_types( bond_type_id) > 1)
          {
            sum_formula << 'x' << bond_types( bond_type_id);
          }

          sum_formula << ',';
        }
      }

      // return hash
      return sum_formula.str();
    }

    //! @brief function used internally to search for a given fragment within member data structures
    //! @param FRAGMENT a constitution to search for
    //! @param NODES molecules to compare FRAGMENT against
    //! @param ISOMORPHISM isomorphism vector to store isomorphism if it is calculated
    //! @param GRAPH reference to a const graph, will be set to the graph for FRAGMENT if a graph is created
    ConstitutionSet::const_iterator ConstitutionSet::InternalFind
    (
      const ConstitutionInterface &FRAGMENT,
      const std::list< Node> &NODES,
      storage::Vector< size_t> &ISOMORPHISM,
      graph::ConstGraph< size_t, size_t> &GRAPH
    ) const
    {
      // handle the simple case where no configurations are present
      if( NODES.empty())
      {
        return End();
      }

      static ConstitutionGraphConverter make_graph;

      GRAPH = make_graph( FRAGMENT);
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
          ISOMORPHISM = isomorphism.GetInverseIsomorphism();
          return itr_node->m_Itr;
        }
      }

      return End();
    }

  } // namespace chemistry
} // namespace bcl
