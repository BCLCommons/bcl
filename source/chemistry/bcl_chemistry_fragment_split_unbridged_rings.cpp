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
#include "chemistry/bcl_chemistry_fragment_split_unbridged_rings.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_connectivity.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentSplitUnbridgedRings::s_InstanceChain
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new FragmentSplitUnbridgedRings( false))
    );
    const util::SiPtr< const util::ObjectInterface> FragmentSplitUnbridgedRings::s_InstanceRing
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new FragmentSplitUnbridgedRings( true))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! @param GET_RINGS true if rings are desired, false if chains are desired
    //! @param MIN_SIZE get the minimum size of ring or chain that is desired
    FragmentSplitUnbridgedRings::FragmentSplitUnbridgedRings( bool GET_AROMATIC, const size_t MAX_SIZE) :
      m_Aromatic( GET_AROMATIC),
      m_MaxSize( MAX_SIZE)
    {
    }

    //! virtual copy constructor
    FragmentSplitUnbridgedRings *FragmentSplitUnbridgedRings::Clone() const
    {
      return new FragmentSplitUnbridgedRings( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &FragmentSplitUnbridgedRings::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentSplitUnbridgedRings::GetAlias() const
    {
      static const std::string s_name_ring( "UnbridgedRings");
      static const std::string s_name_aro( "UnbridgedAromaticRings");
      return m_Aromatic ? s_name_aro : s_name_ring;
    }

    //! @brief Get a description for what this class does (used when writing help)
    //! @return a description for what this class does (used when writing help)
    const std::string &FragmentSplitUnbridgedRings::GetClassDescription() const
    {
      static const std::string s_desc_ring( "Gets all aromatic systems that do not contain bridges (aromatic or otherwise)");
      static const std::string s_desc_chain( "Gets all ring systems that do not contain bridges");
      return m_Aromatic ? s_desc_ring : s_desc_chain;
    }

    //! get the minimum size of a component of interest
    const size_t FragmentSplitUnbridgedRings::GetMaxSize() const
    {
      return m_MaxSize;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief returns list of ring or chain fragments
    //! @param MOLECULE molecule of interest
    //! @param MOLECULE_GRAPH graph of molecule of interest
    //! @return list of ring or chain fragments
    storage::List< storage::Vector< size_t> > FragmentSplitUnbridgedRings::GetComponentVertices
    (
      const ConformationInterface &MOLECULE,
      ConformationGraphConverter::t_AtomGraph &MOLECULE_GRAPH
    ) const
    {
      // declare a list for storing atom indices of rings/chain
      storage::List< storage::Vector< size_t> > indices;

      ConfigurationalBondTypeData::Data data_to_check
      (
        m_Aromatic ? ConfigurationalBondTypeData::e_IsAromatic : ConfigurationalBondTypeData::e_IsInRing
      );

      // count ring bonds or chain bonds
      const size_t number_ring_bonds
      (
        MOLECULE.CountNonValenceBondsWithProperty( data_to_check, size_t( 1))
      );
      // skip molecules that have no ring bonds
      if( number_ring_bonds == size_t( 0))
      {
        // molecule has no ring bonds, so there are no rings to isolate
        return storage::List< storage::Vector< size_t> >();
      }

      // determine whether each bond type exists in a ring
      storage::Vector< sdf::BondInfo> bond_info( MOLECULE.GetBondInfo());

      // remove all bonds outside rings
      for
      (
        storage::Vector< sdf::BondInfo>::const_iterator
          itr_bond( bond_info.Begin()), itr_bond_end( bond_info.End());
        itr_bond != itr_bond_end;
        ++itr_bond
      )
      {
        if( itr_bond->GetConfigurationalBondType()->GetBondData( data_to_check) != size_t( 1))
        {
          MOLECULE_GRAPH.RemoveEdge( itr_bond->GetAtomIndexLow(), itr_bond->GetAtomIndexHigh());
        }
      }
      size_t n_bridgeheads( 0);
      for( iterate::Generic< const AtomConformationalInterface> itr( MOLECULE.GetAtomsIterator()); itr.NotAtEnd(); ++itr)
      {
        if( IsBridgeHead( *itr))
        {
          ++n_bridgeheads;
          MOLECULE_GRAPH.RemoveAllEdges( itr.GetPosition());
        }
      }
      // remove all atoms in the system with the bridgehead so that we don't need to check the number of bonds later
      bool removed_edges( n_bridgeheads);
      while( removed_edges)
      {
        removed_edges = false;
        for( size_t i( 0), n_atoms( MOLECULE_GRAPH.GetSize()); i < n_atoms; ++i)
        {
          if( MOLECULE_GRAPH.GetNeighborIndices( i).GetSize() == size_t( 1))
          {
            removed_edges = true;
            MOLECULE_GRAPH.RemoveAllEdges( i);
          }
        }
      }

      // get components of graph
      storage::List< storage::Vector< size_t> > components( graph::Connectivity::GetComponents( MOLECULE_GRAPH));

      for
      (
        storage::List< storage::Vector< size_t> >::const_iterator itr( components.Begin()), itr_end( components.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->GetSize() > size_t( 1) && itr->GetSize() < m_MaxSize)
        {
          indices.PushBack( *itr);
        }
      }
      return indices;
    }

    bool FragmentSplitUnbridgedRings::IsBridgeHead( const AtomConformationalInterface &ATOM)
    {
      size_t number_bonds( ATOM.GetBonds().GetSize());
      if( number_bonds < 3)
      {
        return false;
      }

      size_t ring_bonds( ATOM.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1)));
      if( ring_bonds < 3)
      {
        return false;
      }
      if( ring_bonds >= 4)
      {
        return true;
      }
      for
      (
        storage::Vector< BondConformational>::const_iterator itr_bonds( ATOM.GetBonds().Begin()), itr_bonds_end( ATOM.GetBonds().End());
          itr_bonds != itr_bonds_end;
        ++itr_bonds
      )
      {
        if( !itr_bonds->GetBondType()->IsBondInRing())
        {
          continue;
        }
        if( itr_bonds->GetTargetAtom().CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1)) == 3)
        {
          return false;
        }
      }
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer FragmentSplitUnbridgedRings::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "returns either ring or chans of a molecule");
      return parameters;
    }

  } // namespace chemistry
} // namespace bcl
