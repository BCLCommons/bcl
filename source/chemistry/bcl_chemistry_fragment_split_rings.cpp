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
#include "chemistry/bcl_chemistry_fragment_split_rings.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_connectivity.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentSplitRings::s_InstanceChain
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new FragmentSplitRings( false, 1))
    );
    const util::SiPtr< const util::ObjectInterface> FragmentSplitRings::s_InstanceRing
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new FragmentSplitRings( true))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! @param GET_RINGS true if rings are desired, false if chains are desired
    //! @param MIN_SIZE get the minimum size of ring or chain that is desired
    FragmentSplitRings::FragmentSplitRings( bool GET_RINGS, const size_t MIN_SIZE) :
      m_GetRings( GET_RINGS),
      m_MinSize( MIN_SIZE)
    {
    }

    //! virtual copy constructor
    FragmentSplitRings *FragmentSplitRings::Clone() const
    {
      return new FragmentSplitRings( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &FragmentSplitRings::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentSplitRings::GetAlias() const
    {
      static const std::string s_name_ring( "Rings");
      static const std::string s_name_chain( "Chains");
      return m_GetRings ? s_name_ring : s_name_chain;
    }

    //! @brief Get a description for what this class does (used when writing help)
    //! @return a description for what this class does (used when writing help)
    const std::string &FragmentSplitRings::GetClassDescription() const
    {
      static const std::string s_desc_ring( "Gets all complete ring systems");
      static const std::string s_desc_chain( "Gets all complete chains");
      return m_GetRings ? s_desc_ring : s_desc_chain;
    }

    //! get the minimum size of a component of interest
    const size_t FragmentSplitRings::GetMinSize() const
    {
      return m_MinSize;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief returns list of ring or chain fragments
    //! @param MOLECULE molecule of interest
    //! @param MOLECULE_GRAPH graph of molecule of interest
    //! @return list of ring or chain fragments
    storage::List< storage::Vector< size_t> > FragmentSplitRings::GetComponentVertices
    (
      const ConformationInterface &MOLECULE,
      ConformationGraphConverter::t_AtomGraph &MOLECULE_GRAPH
    ) const
    {
      // declare a list for storing atom indices of rings/chain
      storage::List< storage::Vector< size_t> > indices;

      // count ring bonds or chain bonds
      const size_t number_ring_bonds
      (
        MOLECULE.CountNonValenceBondsWithProperty
        (
          ConfigurationalBondTypeData::e_IsInRing,
          size_t( m_GetRings)
        )
      );
      // skip molecules that have no ring bonds
      if( number_ring_bonds == size_t( 0))
      {
        // molecule has no ring bonds, so there are no rings to isolate
        return storage::List< storage::Vector< size_t> >();
      }

      if( number_ring_bonds == MOLECULE.GetNumberBonds())
      {
        // keep track of which molecule this ring came from
        linal::Vector< size_t> molecule_atoms( linal::FillVector< size_t>( MOLECULE.GetNumberAtoms(), size_t( 0), size_t( 1)));
        return storage::List< storage::Vector< size_t> >
        (
          1,
          storage::Vector< size_t>( molecule_atoms.Begin(), molecule_atoms.End())
        );
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
        if( itr_bond->GetConfigurationalBondType()->IsBondInRing() != m_GetRings)
        {
          MOLECULE_GRAPH.RemoveEdge( itr_bond->GetAtomIndexLow(), itr_bond->GetAtomIndexHigh());
        }
      }

      // get components of graph
      storage::List< storage::Vector< size_t> > components( graph::Connectivity::GetComponents( MOLECULE_GRAPH));

      // add the components to the list that will be returned from the function
      for
      (
        storage::List< storage::Vector< size_t> >::const_iterator itr( components.Begin()), itr_end( components.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->GetSize() > GetMinSize())
        {
          indices.PushBack( *itr);
        }
      }
      return indices;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer FragmentSplitRings::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "returns either ring or chains of a molecule");
      return parameters;
    }

  } // namespace chemistry
} // namespace bcl
