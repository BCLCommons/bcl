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
#include "chemistry/bcl_chemistry_fragment_split_scaffolds.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_connectivity.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentSplitScaffolds::s_Instance
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new FragmentSplitScaffolds( 2, false))
    );
    const util::SiPtr< const util::ObjectInterface> FragmentSplitScaffolds::s_InverseInstance
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new FragmentSplitScaffolds( 1, true))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! @param MIN_SIZE get the minimum size of scaffold that is desired
    FragmentSplitScaffolds::FragmentSplitScaffolds( const size_t MIN_SIZE, const bool &INVERT) :
      m_MinSize( MIN_SIZE),
      m_Invert( INVERT)
    {
    }

    //! virtual copy constructor
    FragmentSplitScaffolds *FragmentSplitScaffolds::Clone() const
    {
      return new FragmentSplitScaffolds( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &FragmentSplitScaffolds::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentSplitScaffolds::GetAlias() const
    {
      static const std::string s_name_scaffold( "Scaffolds"), s_name_inverse( "InverseScaffold");
      return m_Invert ? s_name_inverse : s_name_scaffold;
    }

    //! get the minimum size of a component of interest
    const size_t FragmentSplitScaffolds::GetMinSize() const
    {
      return m_MinSize;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief returns list of scaffolds that are contained in the molecule of interest
    //! @param MOLECULE molecule of interest
    //! @param MOLECULE_GRAPH graph of molecule of interest
    //! @return list of scaffolds that are contained in the molecule of interest
    storage::List< storage::Vector< size_t> > FragmentSplitScaffolds::GetComponentVertices
    (
      const ConformationInterface &MOLECULE,
      ConformationGraphConverter::t_AtomGraph &MOLECULE_GRAPH
    ) const
    {
      // atom status is one of
      // 0 - Terminal or in terminal chain
      // 1 - In chain, possibly terminal
      // 2 - In ring
      storage::Vector< size_t> atom_status( MOLECULE.GetNumberAtoms(), size_t( 0));
      storage::Vector< size_t>::iterator itr_atom_status( atom_status.Begin());
      size_t status_unknown( 0); // number of atoms in chain for which it is unknown whether they are part of a terminal chain
      size_t n_ring_atoms( 0);
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr( MOLECULE.GetAtomsIterator());
          itr.NotAtEnd();
        ++itr, ++itr_atom_status
      )
      {
        if( itr->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, 1) > 0)
        {
          // in ring
          *itr_atom_status = 2;
          ++n_ring_atoms;
        }
        else if( itr->GetBonds().GetSize() > 1)
        {
          // non-terminal atom
          *itr_atom_status = 1;
          ++status_unknown;
        }
      }

      // detect no ring atoms -> no scaffold
      if( n_ring_atoms == size_t( 0))
      {
        if( !m_Invert)
        {
          return storage::List< storage::Vector< size_t> >();
        }
        // inverse scaffold, return whole molecule
        linal::Vector< size_t> molecule_atoms( linal::FillVector< size_t>( MOLECULE.GetNumberAtoms(), size_t( 0), size_t( 1)));
        return storage::List< storage::Vector< size_t> >
        (
          1,
          storage::Vector< size_t>( molecule_atoms.Begin(), molecule_atoms.End())
        );
      }

      // initialize previous count of status unknown to setup finding terminal chains
      size_t old_status_unknown( status_unknown + 1);
      while( old_status_unknown > status_unknown && status_unknown)
      {
        old_status_unknown = status_unknown;
        status_unknown = 0;

        size_t id( 0); // atom id
        for
        (
          iterate::Generic< const AtomConformationalInterface> itr( MOLECULE.GetAtomsIterator());
          itr.NotAtEnd();
          ++itr, ++id
        )
        {
          if( atom_status( id) != size_t( 1))
          {
            continue;
          }

          size_t neighbors_not_in_terminal_chain( 0);
          for
          (
            storage::Vector< BondConformational>::const_iterator
              itr_bond( itr->GetBonds().Begin()), itr_bond_end( itr->GetBonds().End());
            itr_bond != itr_bond_end;
            ++itr_bond
          )
          {
            // get the atom index
            const size_t atom_index( MOLECULE.GetAtomIndex( itr_bond->GetTargetAtom()));

            // get atom status
            const size_t bonded_atom_status( atom_status( atom_index));

            if( bonded_atom_status != size_t( 0))
            {
              ++neighbors_not_in_terminal_chain;
              if( neighbors_not_in_terminal_chain == 2)
              {
                break;
              }
            }
          }

          if( neighbors_not_in_terminal_chain < size_t( 2))
          {
            atom_status( id) = size_t( 0);
          }
          else
          {
            // non-terminal atom
            ++status_unknown;
          }
        }
      }

      if( !m_Invert)
      {
        // find atoms that have been colored as not being in terminal rings
        storage::Vector< size_t> non_terminal_chain_atom_ids;
        non_terminal_chain_atom_ids.AllocateMemory( atom_status.GetSize());
        for( size_t atom_id( 0), n_atoms( atom_status.GetSize()); atom_id < n_atoms; ++atom_id)
        {
          if( atom_status( atom_id))
          {
            non_terminal_chain_atom_ids.PushBack( atom_id);
          }
        }
        if( non_terminal_chain_atom_ids.GetSize() > GetMinSize())
        {
          // get connected components of the graph, which are either isolated atoms or rings
          return storage::List< storage::Vector< size_t> >( 1, non_terminal_chain_atom_ids);
        }
        else
        {
          return storage::List< storage::Vector< size_t> >();
        }
      }

      // change the atom status for any atom adjacent to a terminal-chain directed atom
      storage::Vector< size_t> new_atom_status( atom_status);
      for( size_t atom_id( 0), n_atoms( atom_status.GetSize()); atom_id < n_atoms; ++atom_id)
      {
        if( !atom_status( atom_id))
        {
          for
          (
            auto itrb( MOLECULE_GRAPH.GetNeighborIndices( atom_id).Begin()),
                 itrb_end( MOLECULE_GRAPH.GetNeighborIndices( atom_id).End());
            itrb != itrb_end;
            ++itrb
          )
          {
            new_atom_status( *itrb) = 0;
          }
        }
      }

      storage::Vector< size_t> terminal_chain_atom_ids;
      terminal_chain_atom_ids.AllocateMemory( atom_status.GetSize());
      for( size_t atom_id( 0), n_atoms( atom_status.GetSize()); atom_id < n_atoms; ++atom_id)
      {
        if( !new_atom_status( atom_id))
        {
          terminal_chain_atom_ids.PushBack( atom_id);
        }
      }
      auto subgraph( MOLECULE_GRAPH.GetSubgraph( terminal_chain_atom_ids));
      auto components( graph::Connectivity::GetComponents( subgraph));
      for( auto itr( components.Begin()), itr_end( components.End()); itr != itr_end;)
      {
        if( itr->GetSize() > GetMinSize())
        {
          // update atom indices, since they relate to the subgraph induced by terminal atom ids
          for( auto itrc( itr->Begin()), itrc_end( itr->End()); itrc != itrc_end; ++itrc)
          {
            *itrc = terminal_chain_atom_ids( *itrc);
          }
          ++itr;
        }
        else
        {
          // component too small to be of interest; remove
          auto itrp( itr);
          ++itr;
          components.Remove( itrp);
        }
      }
      return components;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer FragmentSplitScaffolds::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        m_Invert ?
        "returns the remaining components of a molecule after the murcko scaffold is removed"
         : "returns scaffolds that are present in the molecule of interest"
      );
      return parameters;
    }

  } // namespace chemistry
} // namespace bcl
