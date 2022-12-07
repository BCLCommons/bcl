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
#include "chemistry/bcl_chemistry_fragment_mutate_ring_swap.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_align_to_scaffold.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_mutate_bond_angles.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "find/bcl_find_collector_interface.h"
#include "graph/bcl_graph_edge_cover_ring_perception.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentMutateRingSwap::s_Instance
    (
      util::Enumerated< FragmentMutateInterface>::AddInstance( new FragmentMutateRingSwap())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentMutateRingSwap::FragmentMutateRingSwap() :
      m_RotamerLibrarySearcher( util::ShPtr< SearchFragmentLibraryFromTree>()),
      m_RingInitiationProbability( 0.1),
      m_FixGeometry( true),
      m_Neutralize( false),
      m_RestrictToNoMoreThanOneRingSizeChange( true),
      m_AllowLargeRingCollapse( true),
      m_AlignRings( false),
      m_ExtendAdjacentAtoms( size_t( 1)),
      m_ChooseBestAlignedConf( false),
      m_BondComparisonType( ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness),
      m_AtomComparisonType( ConformationGraphConverter::e_ElementType)
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief pose-dependent constructor
    FragmentMutateRingSwap::FragmentMutateRingSwap
    (
      const util::ShPtr< SearchFragmentLibraryFromTree> &FRAGMENT_LIBRARY,
      const std::string &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const std::string &MDL,
      const descriptor::CheminfoProperty &PROPERTY_SCORER,
      const bool &RESOLVE_CLASHES,
      const storage::Vector< float> &BFACTORS,
      const bool &CORINA,
      const bool &FIX_GEOMETRY,
      const bool &NEUTRALIZE,
      const double &RING_INITIATION_PROBABILITY,
      const bool &PREVENT_MORE_THAN_ONE_RING_FROM_CHANGING_SIZE,
      const bool &ALLOW_LARGE_RING_COLLAPSE
    ) :
      m_RotamerLibrarySearcher( FRAGMENT_LIBRARY),
      m_RingInitiationProbability( RING_INITIATION_PROBABILITY),
      m_FixGeometry( FIX_GEOMETRY),
      m_Neutralize( NEUTRALIZE),
      m_RestrictToNoMoreThanOneRingSizeChange( PREVENT_MORE_THAN_ONE_RING_FROM_CHANGING_SIZE),
      m_AllowLargeRingCollapse( ALLOW_LARGE_RING_COLLAPSE),
      m_AlignRings( false),
      m_ExtendAdjacentAtoms( size_t( 1)),
      m_ChooseBestAlignedConf( false),
      m_BondComparisonType( ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness),
      m_AtomComparisonType( ConformationGraphConverter::e_ElementType)
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_ScaffoldFragment = SCAFFOLD_FRAGMENT;
      m_MutableFragments = MUTABLE_FRAGMENTS;
      m_MutableAtomIndices = MUTABLE_ATOM_INDICES;
      m_MDL = MDL;
      m_PropertyScorer = PROPERTY_SCORER;
      m_ResolveClashes = RESOLVE_CLASHES;
      m_BFactors = BFACTORS;
      m_Corina = CORINA;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief pose-independent constructor
    FragmentMutateRingSwap::FragmentMutateRingSwap
    (
      const util::ShPtr< SearchFragmentLibraryFromTree> &FRAGMENT_LIBRARY,
      const std::string &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const bool &CORINA,
      const bool &FIX_GEOMETRY,
      const bool &NEUTRALIZE,
      const double &RING_INITIATION_PROBABILITY,
      const bool &PREVENT_MORE_THAN_ONE_RING_FROM_CHANGING_SIZE,
      const bool &ALLOW_LARGE_RING_COLLAPSE
    ) :
      m_RotamerLibrarySearcher( FRAGMENT_LIBRARY),
      m_RingInitiationProbability( RING_INITIATION_PROBABILITY),
      m_FixGeometry( FIX_GEOMETRY),
      m_Neutralize( NEUTRALIZE),
      m_RestrictToNoMoreThanOneRingSizeChange( PREVENT_MORE_THAN_ONE_RING_FROM_CHANGING_SIZE),
      m_AllowLargeRingCollapse( ALLOW_LARGE_RING_COLLAPSE),
      m_AlignRings( false),
      m_ExtendAdjacentAtoms( size_t( 1)),
      m_ChooseBestAlignedConf( false),
      m_BondComparisonType( ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness),
      m_AtomComparisonType( ConformationGraphConverter::e_ElementType)
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_ScaffoldFragment = SCAFFOLD_FRAGMENT;
      m_MutableFragments = MUTABLE_FRAGMENTS;
      m_MutableAtomIndices = MUTABLE_ATOM_INDICES;
      m_Corina = CORINA;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief clone constructor
    FragmentMutateRingSwap *FragmentMutateRingSwap::Clone() const
    {
      return new FragmentMutateRingSwap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentMutateRingSwap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentMutateRingSwap::GetAlias() const
    {
      static const std::string s_name( "RingSwap");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @return MutateResult with Constitution after the mutate
    math::MutateResult< FragmentComplete> FragmentMutateRingSwap::operator()( const FragmentComplete &FRAGMENT) const
    {
      BCL_MessageStd( "RingSwap!");

      // choose whether to create/remove a ring
      AtomVector< AtomComplete> atoms( FRAGMENT.GetAtomVector());
      if( random::GetGlobalRandom().Double() < m_RingInitiationProbability)
      {
        // expand or remove a ring
        if( random::GetGlobalRandom().Double() < 0.5)
        {
          // pick random atom to transform
          util::SiPtr< const AtomConformationalInterface> picked_atom;
          size_t chosen_atom( util::GetUndefinedSize_t());
          if( m_MutableAtomIndices.GetSize() || m_MutableElements.GetSize() || m_MutableFragments.GetSize())
          {
            picked_atom = this->PickAtom( FRAGMENT, false);
            size_t picked_atom_index( FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom));
            chosen_atom = GetRandomNonringAtom( FRAGMENT, picked_atom_index);
          }
          else
          {
            picked_atom = this->PickAtom( FRAGMENT, true);
            chosen_atom = GetGlobalRandomNonringAtom( FRAGMENT);
          }

          // create a new ring
          if( !util::IsDefined( chosen_atom))
          {
            // signal no-op
            return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
          }
          size_t n_double_bonds
          (
            atoms( chosen_atom).GetAtomType()->GetNumberElectronsInBonds() - atoms( chosen_atom).GetAtomType()->GetNumberBonds()
          );
          if( n_double_bonds >= m_FragmentPool.GetSize())
          {
            // signal no-op
            return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
          }

          // collect substituents of the atom
          auto substituents( GetSubstitutents( FRAGMENT, storage::Vector< size_t>( size_t( 1), chosen_atom), false)( 0));

          // select a ring
          bool found_a_good_one( false);
          auto itr_ring_select( m_FragmentPool( n_double_bonds).Begin()), itr_ring_select_end( m_FragmentPool( n_double_bonds).End());
          while( !found_a_good_one)
          {
            long frag_value( random::GetGlobalRandom().Random< size_t>( 0, m_FragmentPoolScaffoldSums( n_double_bonds) - 1));
            itr_ring_select = m_FragmentPool( n_double_bonds).Begin();
            while( itr_ring_select != itr_ring_select_end)
            {
              frag_value -= GetCounts( *itr_ring_select);
              if( frag_value < 0)
              {
                break;
              }
              ++itr_ring_select;
            }
            if( itr_ring_select == itr_ring_select_end)
            {
              --itr_ring_select;
            }
            found_a_good_one = true;
            if( m_RestrictToNoMoreThanOneRingSizeChange)
            {
              if( GetRingSizeInformation( *itr_ring_select).GetSize() != size_t( 1))
              {
                found_a_good_one = false;
              }
            }
          }
          // determine where substituents will go on the new rings by collecting valences on the new ring
          storage::Vector< storage::Vector< size_t> > ring_valence_atom_indices( size_t( 2));
          for( auto itr_ring( itr_ring_select->GetAtomsIterator()); itr_ring.NotAtEnd(); ++itr_ring)
          {
            size_t n_single_valence_bonds( itr_ring->GetNumberofValenceBondsWithOrder( 1));
            size_t n_double_valence_bonds( itr_ring->GetNumberofValenceBondsWithOrder( 2));
            if( n_single_valence_bonds)
            {
              ring_valence_atom_indices( 0).Append
                  (
                    storage::Vector< size_t>( n_single_valence_bonds, itr_ring.GetPosition() + atoms.GetSize())
                  );
            }
            if( n_double_valence_bonds)
            {
              ring_valence_atom_indices( 1).Append
                  (
                    storage::Vector< size_t>( n_double_valence_bonds, itr_ring.GetPosition() + atoms.GetSize())
                  );
            }
          }
          storage::Vector< sdf::BondInfo> new_bonds;

          // add atoms with connectivity to the new ring. and substituents to the old atom
          ring_valence_atom_indices( 0).Shuffle();
          ring_valence_atom_indices( 1).Shuffle();
          new_bonds.Reset();
          new_bonds.AllocateMemory( substituents.GetSize());
          size_t single_valence_count( 0), double_valence_count( 0);

          for( auto itr( substituents.Begin()), itr_end( substituents.End()); itr != itr_end; ++itr)
          {

            // avoid trying to add too many substituents to a smaller ring
            if( single_valence_count >= itr_ring_select->GetSize())
            {
              break;
            }

            if( itr->Second()->GetNumberOfElectrons() == size_t( 2))
            {
              if( single_valence_count == ring_valence_atom_indices( 0).GetSize())
              {
                continue; // too many double bond valences
              }
              new_bonds.PushBack( sdf::BondInfo( itr->First(), ring_valence_atom_indices( 0)( single_valence_count), itr->Second()));
              ++single_valence_count;
            }
            else
            {
              if( double_valence_count == ring_valence_atom_indices( 1).GetSize())
              {
                continue; // too many double bond valences
              }
              new_bonds.PushBack( sdf::BondInfo( itr->First(), ring_valence_atom_indices( 1)( double_valence_count), itr->Second()));
              ++double_valence_count;
            }
          }
          atoms.AddAtomsWithConnectivity( itr_ring_select->GetAtomVector(), new_bonds);

          // Reorder to remove the old atom
          storage::Vector< size_t> indices_to_keep( storage::CreateIndexVector( atoms.GetSize()));
          indices_to_keep.RemoveElements( chosen_atom, 1);
          atoms.Reorder( indices_to_keep);
        }
        else
        {
          // collapse an old ring system, if possible
          // find a collapsible ring
          // pick random atom to transform
          util::SiPtr< const AtomConformationalInterface> picked_atom;
          storage::Vector< size_t> ring_to_collapse;
          if( m_MutableAtomIndices.GetSize() || m_MutableElements.GetSize() || m_MutableFragments.GetSize())
          {
            picked_atom = this->PickAtom( FRAGMENT, false);
            size_t picked_atom_index( FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom));
            ring_to_collapse = GetRandomRingSystem( FRAGMENT, picked_atom_index, m_AllowLargeRingCollapse);
          }
          else
          {
            picked_atom = this->PickAtom( FRAGMENT, true);
            ring_to_collapse = GetGlobalRandomRingSystem( FRAGMENT, m_AllowLargeRingCollapse);
          }

          if( ring_to_collapse.IsEmpty())
          {
            return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
          }

          // ring center
          math::RunningAverage< linal::Vector3D> ring_center;
          for( size_t ring_index( 0), ring_sz( ring_to_collapse.GetSize()); ring_index < ring_sz; ++ring_index)
          {
            ring_center += FRAGMENT.GetAtomVector()( ring_to_collapse( ring_index)).GetPosition();
          }
          storage::Vector< sdf::AtomInfo> single_atom
          (
            size_t( 1),
            sdf::AtomInfo
            (
              GetAtomTypes().GetAtomType( GetElementTypes().e_Carbon, 0),
              e_NonChiral,
              ring_center.GetAverage()
            )
          );

          // collect substituents
          storage::Vector< sdf::BondInfo> bonds;
          auto substituents( GetSubstitutents( FRAGMENT, ring_to_collapse, true));
          substituents( 0).Shuffle();
          size_t n_electrons( 0);
          for
          (
              auto itr_substituents( substituents( 0).Begin()), itr_substituents_end( substituents( 0).End());
              itr_substituents != itr_substituents_end;
              ++itr_substituents
          )
          {
            if( ( n_electrons += itr_substituents->Second()->GetNumberOfElectrons() / 2) > 4)
            {
              break;
            }
            bonds.PushBack( sdf::BondInfo( atoms.GetSize(), itr_substituents->First(), itr_substituents->Second()));
          }
          AtomVector< AtomComplete> new_atoms( single_atom, storage::Vector< sdf::BondInfo>());
          atoms.AddAtomsWithConnectivity( new_atoms, bonds);
          storage::Set< size_t> indices_to_remove( ring_to_collapse.Begin(), ring_to_collapse.End());
          storage::Vector< size_t> reorder;
          reorder.AllocateMemory( atoms.GetSize() - ring_to_collapse.GetSize());
          for( size_t i( 0), na( atoms.GetSize()); i < na; ++i)
          {
            if( !indices_to_remove.Contains( i))
            {
              reorder.PushBack( i);
            }
          }
          atoms.Reorder( reorder);
        }
      }
      else
      {
        // swap a ring
        // determine where substituents will go on the new rings by collecting valences on the new ring
        // add atoms with connectivity to the new ring. and substituents to the old atom
        // Reorder to remove the old atom
        // create a new ring
        util::SiPtr< const AtomConformationalInterface> picked_atom;
        storage::Vector< size_t> chosen_ring_system;

        if( m_MutableAtomIndices.GetSize() || m_MutableElements.GetSize() || m_MutableFragments.GetSize())
        {
          picked_atom = this->PickAtom( FRAGMENT, false);
          size_t picked_atom_index( FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom));
          chosen_ring_system = GetRandomRingSystem( FRAGMENT, picked_atom_index, m_AllowLargeRingCollapse);
        }
        else
        {
          picked_atom = this->PickAtom( FRAGMENT, true);
          chosen_ring_system = GetGlobalRandomRingSystem( FRAGMENT, m_AllowLargeRingCollapse);
        }
        if( chosen_ring_system.IsEmpty())
        {
          // signal no-op
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }

        // double bonds
        size_t n_double_bonds( GetNumberDoubleBondValences( FRAGMENT, chosen_ring_system));
        if( n_double_bonds >= m_FragmentPool.GetSize())
        {
          // signal no-op
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }
        // collect substituents of the atom
        auto substituents( GetSubstitutents( FRAGMENT, chosen_ring_system, true)( 0));

        // select a ring
        bool found_a_good_one( false);
        auto ring_first_info
        (
          m_RestrictToNoMoreThanOneRingSizeChange
          ? GetRingSizeInformation( FRAGMENT, chosen_ring_system)
              : linal::Vector< float>()
        );
        size_t trial( 0);
        auto itr_ring_select( m_FragmentPool( n_double_bonds).Begin()), itr_ring_select_end( m_FragmentPool( n_double_bonds).End());
        while( !found_a_good_one && trial < FRAGMENT.GetSize())
        {
          long frag_value( random::GetGlobalRandom().Random< size_t>( 0, m_FragmentPoolScaffoldSums( n_double_bonds) - 1));
          itr_ring_select = m_FragmentPool( n_double_bonds).Begin();
          while( itr_ring_select != itr_ring_select_end)
          {
            frag_value -= GetCounts( *itr_ring_select);
            if( frag_value < 0)
            {
              break;
            }
            ++itr_ring_select;
          }
          if( itr_ring_select == itr_ring_select_end)
          {
            --itr_ring_select;
          }
          found_a_good_one = true;
          if( m_RestrictToNoMoreThanOneRingSizeChange)
          {
            if
            (
                SortedVectorsDifferByNoMoreThanOneElement
                (
                  GetRingSizeInformation( *itr_ring_select),
                  GetRingSizeInformation( FRAGMENT, chosen_ring_system)
                )
            )
            {
              found_a_good_one = false;
              ++trial;
            }
          }
        }

        // try to preserve substituent relative distances when placed on new ring
        storage::Vector< sdf::BondInfo> new_bonds;
        new_bonds.Reset();
        new_bonds.AllocateMemory( substituents.GetSize());
        if( m_AlignRings)
        {
          // identify which atoms in the selected ring currently attach to the substituents
          storage::Map< size_t, size_t> substituents_to_ringatoms;

          // go over our substituent atoms
          for
          (
              auto substituent_itr( substituents.Begin()), substituents_itr_end( substituents.End());
              substituent_itr != substituents_itr_end;
              ++substituent_itr
          )
          {
            // find the bonds made by the substituent atom
            auto &atom( FRAGMENT.GetAtomVector()( substituent_itr->First()));
            for
            (
                auto bond_itr( atom.GetBonds().Begin()), bond_itr_end( atom.GetBonds().End());
                bond_itr != bond_itr_end;
                ++bond_itr
            )
            {
              // check if the bond is attached to one of our ring atoms
              for( size_t ring_atom_i( 0), n_ring_atoms( chosen_ring_system.GetSize()); ring_atom_i < n_ring_atoms; ++ring_atom_i)
              {
                if( FRAGMENT.GetAtomIndex( bond_itr->GetTargetAtom()) == chosen_ring_system( ring_atom_i))
                {
                  // map our substituent to its connecting atom in the selected ring
                  substituents_to_ringatoms.Insert
                  (
                    std::make_pair
                    (
                      FRAGMENT.GetAtomIndex( atom),
                      chosen_ring_system( ring_atom_i)
                    )
                  );
                }
              }
            }
          }

          // get new fragment from ring indices
          AtomVector< AtomComplete> chosen_ring_atom_vec( FRAGMENT.GetAtomVector());
          chosen_ring_atom_vec.Reorder( chosen_ring_system);
          FragmentComplete chosen_ring_frag( chosen_ring_atom_vec, "");
          FragmentComplete ring_select( *itr_ring_select);

          FragmentAlignToScaffold align_to_scaffold
          (
            m_AtomComparisonType,
            m_BondComparisonType,
            size_t( 3)
          );

          bool aligned( align_to_scaffold.AlignToScaffold( ring_select, chosen_ring_frag));

          // use a special graph converter to check for substitutable atoms
          ConformationGraphConverter graph_converter
          (
            m_AtomComparisonType,
            m_BondComparisonType
          );

          // build graphs of current and new rings
          graph::ConstGraph< size_t, size_t> current_ring_graph( graph_converter( chosen_ring_frag));
          graph::ConstGraph< size_t, size_t> new_ring_graph( graph_converter( ring_select));
          util::OwnPtr< graph::ConstGraph< size_t, size_t> > current_ring_graph_ptr( &current_ring_graph, false);
          util::OwnPtr< graph::ConstGraph< size_t, size_t> > new_ring_graph_ptr( &new_ring_graph, false);

          // set graph positions for common iso
          graph::CommonSubgraphIsomorphism< size_t, size_t> common_subgraph_iso( graph::CommonSubgraphIsomorphismBase::e_GreedyUnconnected);
          common_subgraph_iso.SetGraphA( current_ring_graph_ptr);
          common_subgraph_iso.SetGraphB( new_ring_graph_ptr);
  
          // MCS comparison of current and new rings; get the mapping between atoms
          common_subgraph_iso.FindIsomorphism( common_subgraph_iso.EstimateUpperBounds());

          graph::Subgraph< size_t, size_t> current_ring_subgraph( common_subgraph_iso.GetSubgraphIsomorphismsOfGraphA().FirstElement());
          graph::Subgraph< size_t, size_t> new_ring_subgraph( common_subgraph_iso.GetSubgraphIsomorphismsOfGraphB().FirstElement());

          // iterate over our substituents mapped to our original current mol indices
          size_t pair_i( 0);
          for
          (
              auto pair_itr( substituents_to_ringatoms.Begin()), pair_itr_end( substituents_to_ringatoms.End());
              pair_itr != pair_itr_end;
              ++pair_itr, ++pair_i
          )
          {
            // key is substituent, value is original index of attached atom in current ring
            // so let's loop over our current mol subgraph and grab those indices
            for
            (
                size_t subgraph_i( 0), subgraph_i_sz( current_ring_subgraph.GetVertexIndices().GetSize());
                subgraph_i < subgraph_i_sz;
                ++subgraph_i
            )
            {
              // map the original indices of the current ring to the MCS atoms in the new ring
              size_t original_i( chosen_ring_system( current_ring_subgraph.GetVertexIndices()( subgraph_i)));
              if( pair_itr->second == original_i)
              {
                new_bonds.PushBack
                (
                  sdf::BondInfo
                  (
                    pair_itr->first,
                    new_ring_subgraph.GetVertexIndices()( subgraph_i) + FRAGMENT.GetSize(),
                    substituents( pair_i).Second()
                  )
                );
                break;
              }
            }
          }
          atoms.AddAtomsWithConnectivity( ring_select.GetAtomVector(), new_bonds);

          // Reorder to remove the old ring
          storage::Set< size_t> indices_to_remove( chosen_ring_system.Begin(), chosen_ring_system.End());
          storage::Vector< size_t> reorder;
          reorder.AllocateMemory( atoms.GetSize() - ring_select.GetSize());
          for( size_t i( 0), na( atoms.GetSize()); i < na; ++i)
          {
            if( !indices_to_remove.Contains( i))
            {
              reorder.PushBack( i);
            }
          }
          atoms.Reorder( reorder);
        }

        // randomly assign valid substituent connections to new ring
        else
        {
          // determine where substituents will go on the new rings by collecting valences on the new ring
          storage::Vector< storage::Vector< size_t> > ring_valence_atom_indices( size_t( 2));
          for( auto itr_ring( itr_ring_select->GetAtomsIterator()); itr_ring.NotAtEnd(); ++itr_ring)
          {
            size_t n_single_valence_bonds( itr_ring->GetNumberofValenceBondsWithOrder( 1));
            size_t n_double_valence_bonds( itr_ring->GetNumberofValenceBondsWithOrder( 2));
            if( n_single_valence_bonds)
            {
              ring_valence_atom_indices( 0).Append
                  (
                    storage::Vector< size_t>( n_single_valence_bonds, itr_ring.GetPosition() + atoms.GetSize())
                  );
            }
            if( n_double_valence_bonds)
            {
              ring_valence_atom_indices( 1).Append
                  (
                    storage::Vector< size_t>( n_double_valence_bonds, itr_ring.GetPosition() + atoms.GetSize())
                  );
            }
          }
          ring_valence_atom_indices( 0).Shuffle();
          ring_valence_atom_indices( 1).Shuffle();
          size_t single_valence_count( 0), double_valence_count( 0);

          // add atoms with connectivity to the new ring. and substituents to the old atom
          for( auto itr( substituents.Begin()), itr_end( substituents.End()); itr != itr_end; ++itr)
          {
            // avoid trying to add too many substituents to a smaller ring
            if( single_valence_count >= itr_ring_select->GetSize())
            {
              break;
            }

            if( itr->Second()->GetNumberOfElectrons() == size_t( 2))
            {
              if( single_valence_count == ring_valence_atom_indices( 0).GetSize())
              {
                continue; // too many double bond valences
              }
              new_bonds.PushBack( sdf::BondInfo( itr->First(), ring_valence_atom_indices( 0)( single_valence_count), itr->Second()));
              ++single_valence_count;
            }
            else
            {
              if( double_valence_count == ring_valence_atom_indices( 1).GetSize())
              {
                continue; // too many double bond valences
              }
              new_bonds.PushBack( sdf::BondInfo( itr->First(), ring_valence_atom_indices( 1)( double_valence_count), itr->Second()));
              ++double_valence_count;
            }
          }
          atoms.AddAtomsWithConnectivity( itr_ring_select->GetAtomVector(), new_bonds);

          // Reorder to remove the old ring
          storage::Set< size_t> indices_to_remove( chosen_ring_system.Begin(), chosen_ring_system.End());
          storage::Vector< size_t> reorder;
          reorder.AllocateMemory( atoms.GetSize() - itr_ring_select->GetSize());
          for( size_t i( 0), na( atoms.GetSize()); i < na; ++i)
          {
            if( !indices_to_remove.Contains( i))
            {
              reorder.PushBack( i);
            }
          }
          atoms.Reorder( reorder);
        }
      } // end swap ring

      // for cleaning and optimizing the new molecule conformer
      FragmentMapConformer cleaner
      (
        m_DrugLikenessType,
        m_MDL,
        FRAGMENT.GetMDLProperty( m_MDL),
        m_PropertyScorer,
        m_ResolveClashes,
        m_BFactors,
        m_Corina,
        storage::Vector< size_t>(),
        m_ChooseBestAlignedConf,
        m_FixGeometry,
        m_ExtendAdjacentAtoms
      );

      // clean the molecule
      util::ShPtr< FragmentComplete> frag( util::ShPtr< FragmentComplete>( new FragmentComplete( atoms, "")));
      HydrogensHandler::Remove( atoms);
      m_ScaffoldFragment.GetSize()
          ? frag = cleaner.Clean( atoms, m_ScaffoldFragment, m_DrugLikenessType)
          : frag = cleaner.Clean( atoms, FRAGMENT, m_DrugLikenessType);
      // return the new constitution
      return math::MutateResult< FragmentComplete>( frag, *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief compute whether two vectors differ by at most one element
    bool FragmentMutateRingSwap::SortedVectorsDifferByNoMoreThanOneElement
    (
      const linal::Vector< float> &VEC_A,
      const linal::Vector< float> &VEC_B
    )
    {
      auto itr_a( VEC_A.Begin()), itr_a_end( VEC_A.End());
      auto itr_b( VEC_B.Begin()), itr_b_end( VEC_B.End());

      if( VEC_A.GetSize() == VEC_B.GetSize())
      {
        for( size_t n_diff( 0); itr_a != itr_a_end; ++itr_a, ++itr_b)
        {
          if( *itr_a != *itr_b && ++n_diff > size_t( 1))
          {
            return false;
          }
        }
        return true;
      }
      else if( VEC_A.GetSize() == VEC_B.GetSize() + size_t( 1))
      {
        for( ; itr_a != itr_a_end; ++itr_a, ++itr_b)
        {
          if( *itr_a != *itr_b && ( ++itr_a == itr_a_end || *itr_a != *itr_b))
          {
            return false;
          }
        }
        return itr_b == itr_b_end;
      }
      else if( VEC_A.GetSize() + 1 == VEC_B.GetSize())
      {
        for( ; itr_b != itr_b_end; ++itr_a, ++itr_b)
        {
          if( *itr_a != *itr_b && ( ++itr_b == itr_b_end || *itr_a != *itr_b))
          {
            return false;
          }
        }
        return itr_a == itr_a_end;
      }
      return false;
    }

    //! @brief add ring size information as a property to the molecule
    typename descriptor::CacheMap::value_type FragmentMutateRingSwap::GetRingSizeInformation( const FragmentComplete &FRAG)
    {
      static util::ObjectDataLabel s_label( "RingSizes");
      auto cached_val( FRAG.GetFromCache( s_label));
      if( !cached_val.IsEmpty())
      {
        return cached_val;
      }

      ConformationGraphConverter graph_maker;
      auto rings( graph::EdgeCoverRingPerception( graph_maker( FRAG)).GetRings());
      storage::Vector< size_t> ring_sizes;
      for( auto itr_rings( rings.Begin()), itr_rings_end( rings.End()); itr_rings != itr_rings_end; ++itr_rings)
      {
        ring_sizes.PushBack( itr_rings->GetSize());
      }
      ring_sizes.Sort( std::less< size_t>());
      linal::Vector< float> ring_sizes_flt( ring_sizes.Begin(), ring_sizes.End());
      return FRAG.Cache( s_label, ring_sizes_flt);
    }

    //! @brief add ring size information as a property to the molecule
    typename descriptor::CacheMap::value_type FragmentMutateRingSwap::GetRingSizeInformation
    (
      const FragmentComplete &FRAG,
      const storage::Vector< size_t> &RING
    )
    {
      auto grap( ConformationGraphConverter()( FRAG));
      graph::Subgraph< size_t, size_t> subgraph( util::OwnPtr< const graph::ConstGraph< size_t, size_t> >( &grap, false), RING);
      auto rings( graph::EdgeCoverRingPerception( subgraph.ToGraph()).GetRings());
      storage::Vector< size_t> ring_sizes;
      for( auto itr_rings( rings.Begin()), itr_rings_end( rings.End()); itr_rings != itr_rings_end; ++itr_rings)
      {
        ring_sizes.PushBack( itr_rings->GetSize());
      }
      ring_sizes.Sort( std::less< size_t>());
      linal::Vector< float> ring_sizes_flt( ring_sizes.Begin(), ring_sizes.End());
      return ring_sizes_flt;
    }

    //! @brief add ring size information as a property to the molecule
    size_t FragmentMutateRingSwap::GetNumberDoubleBondValences( const FragmentComplete &FRAG)
    {
      static util::ObjectDataLabel s_label( "DoubleBondValences");
      auto cached_val( FRAG.GetFromCache( s_label));
      if( !cached_val.IsEmpty())
      {
        return size_t( cached_val( 0));
      }
      size_t value( 0);
      for( auto itr( FRAG.GetAtomsIterator()); itr.NotAtEnd(); ++itr)
      {
        value += itr->GetNumberofValenceBondsWithOrder( size_t( 2));
      }
      FRAG.Cache( s_label, linal::Vector< float>( size_t( 1), float( value)));
      return value;
    }

    //! @brief add ring size information as a property to the molecule
    size_t FragmentMutateRingSwap::GetNumberDoubleBondValences( const FragmentComplete &FRAG, const storage::Vector< size_t> &RING)
    {
      size_t value( 0);
      for( auto itr( RING.Begin()), itr_end( RING.End()); itr != itr_end; ++itr)
      {
        value += FRAG.GetAtomVector()( *itr).GetNumberofValenceBondsWithOrder( size_t( 2));
      }
      return value;
    }

    //! @brief Get the counts for the given fragment
    size_t FragmentMutateRingSwap::GetCounts( const FragmentComplete &FRAG)
    {
      static util::ObjectDataLabel s_label( "ScaffoldCount");
      if( !FRAG.IsPropertyStored( s_label))
      {
        return size_t( 1);
      }
      return FRAG.GetFromCache( s_label)( 0);
    }

    //! @brief get a random ring system of the given fragment
    storage::Vector< size_t> FragmentMutateRingSwap::GetRandomRingSystem( const FragmentComplete &FRAG, const size_t ATOM_INDEX, bool ALLOW_LARGE_RING_COLLAPSE)
    {
      // create a ring splitter
      FragmentSplitRings ringsplit( size_t( 2), true);

      // create a graph of our molecule of interest
      ConformationGraphConverter::t_AtomGraph mol_graph( ConformationGraphConverter::CreateGraphWithAtoms( FRAG));

      // get atom indices involved in each ring, check if our picked atom index is in it
      auto ring_systems( ringsplit.GetComponentVertices( FRAG, mol_graph));
      size_t ring_index( 0);
      bool this_ring( false);
      for
      (
          auto ring_itr( ring_systems.Begin()), ring_itr_end( ring_systems.End());
          ring_itr != ring_itr_end;
          ++ring_itr, ++ring_index
      )
      {
        for
        (
            auto atom_itr( ring_itr->Begin()), atom_itr_end( ring_itr->End());
            atom_itr != atom_itr_end;
            ++atom_itr
        )
        {
          if( *atom_itr == ATOM_INDEX)
          {
            this_ring = true;
          }
        }
        if( this_ring)
        {
          break;
        }
      }

      // quit if our atom is not in any of the ring systems
      if( !this_ring)
      {
        return storage::Vector< size_t>();
      }

      // start walking over all identified rings
      size_t n_rings( ring_systems.GetSize());
      auto itr_rings( ring_systems.Begin());
      std::advance( itr_rings, ring_index);

      // when true, do not allow rings greater than 4 atoms to collapse
      if( ALLOW_LARGE_RING_COLLAPSE)
      {
        while( !IsCollapsible( FRAG, *itr_rings))
        {
          if( !--n_rings)
          {
            return storage::Vector< size_t>();
          }
          ring_systems.Remove( itr_rings);
          itr_rings = ring_systems.Begin();
          std::advance( itr_rings, ring_index);
        }
      }
      return *itr_rings;
    }

    //! @brief get a random atom that is not in a ring system for the given fragment
    size_t FragmentMutateRingSwap::GetRandomNonringAtom( const FragmentComplete &FRAG, const size_t ATOM_INDEX)
    {
      // enforce that the picked atom is not in a ring or triple bond
      if
      (
          !FRAG.GetAtomVector()( ATOM_INDEX).CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1)) &&
          !FRAG.GetAtomVector()( ATOM_INDEX).CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_BondOrder, size_t( 3))
      )
      {
        return ATOM_INDEX;
      }
        return util::GetUndefinedSize_t();
    }

    //! @brief get a random ring system of the given fragment
    storage::Vector< size_t> FragmentMutateRingSwap::GetGlobalRandomRingSystem( const FragmentComplete &FRAG, bool ALLOW_LARGE_RING_COLLAPSE)
    {
      FragmentSplitRings ringsplit( size_t( 2), true);
      // create a graph
      ConformationGraphConverter::t_AtomGraph mol_graph( ConformationGraphConverter::CreateGraphWithAtoms( FRAG));
      auto ring_systems( ringsplit.GetComponentVertices( FRAG, mol_graph));
      if( ring_systems.IsEmpty())
      {
        return storage::Vector< size_t>();
      }
      size_t n_rings( ring_systems.GetSize());
      auto itr_rings( ring_systems.Begin());
      std::advance( itr_rings, random::GetGlobalRandom().Random< size_t>( 0, n_rings - 1));
      if( ALLOW_LARGE_RING_COLLAPSE)
      {
        while( !IsCollapsible( FRAG, *itr_rings))
        {
          if( !--n_rings)
          {
            return storage::Vector< size_t>();
          }
          ring_systems.Remove( itr_rings);
          itr_rings = ring_systems.Begin();
          std::advance( itr_rings, random::GetGlobalRandom().Random< size_t>( 0, n_rings - 1));
        }
      }
      return *itr_rings;
    }

    //! @brief get a random atom that is not in a ring system for the given fragment
    size_t FragmentMutateRingSwap::GetGlobalRandomNonringAtom( const FragmentComplete &FRAG)
    {
      storage::Vector< size_t> available_chain_atoms;
      for( auto itr( FRAG.GetAtomsIterator()); itr.NotAtEnd(); ++itr)
      {
        // skip triple bonds and ring bonded atoms
        if( !itr->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1)))
        {
          if( !itr->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_BondOrder, size_t( 3)))
          {
            available_chain_atoms.PushBack( itr.GetPosition());
          }
        }
      }
      if( available_chain_atoms.IsEmpty())
      {
        return util::GetUndefinedSize_t();
      }
      size_t n_chain_atms( available_chain_atoms.GetSize());
      return available_chain_atoms( random::GetGlobalRandom().Random< size_t>( 0, n_chain_atms - 1));
    }

    //! @brief return the immediate substituent indices for all atoms of a given ring
    storage::Vector< storage::Vector< storage::Pair< size_t, ConfigurationalBondType> > > FragmentMutateRingSwap::GetSubstitutents
    (
      const FragmentComplete &FRAG,
      const storage::Vector< size_t> &RING,
      const bool &COLLAPSE
    )
    {
      storage::Vector< storage::Vector< storage::Pair< size_t, ConfigurationalBondType> > > substituents( COLLAPSE ? size_t( 1) : RING.GetSize());
      for( size_t ring_index( 0), ring_sz( RING.GetSize()); ring_index < ring_sz; ++ring_index)
      {
        for
        (
          auto itr_bnd( FRAG.GetAtomVector()( RING( ring_index)).GetBonds().Begin()),
               itr_bnd_end( FRAG.GetAtomVector()( RING( ring_index)).GetBonds().End());
          itr_bnd != itr_bnd_end;
          ++itr_bnd
        )
        {
          if( itr_bnd->GetBondType()->IsBondInRing() || itr_bnd->GetTargetAtom().GetAtomType() == GetAtomTypes().H_S)
          {
            continue;
          }
          substituents( COLLAPSE ? size_t( 0) : ring_index).PushBack
          (
            storage::Pair< size_t, ConfigurationalBondType>
            (
              FRAG.GetAtomIndex( itr_bnd->GetTargetAtom()),
              itr_bnd->GetBondType()->GetNumberOfElectrons() < size_t( 4)
              ? GetConfigurationalBondTypes().e_NonConjugatedSingleBond
              : GetConfigurationalBondTypes().e_ConjugatedDoubleBond
            )
          );
        }
      }
      return substituents;
    }

    //! @brief return whether the ring can be collapsed to a single atom
    bool FragmentMutateRingSwap::IsCollapsible
    (
      const FragmentComplete &FRAG,
      const storage::Vector< size_t> &RING
    )
    {
      size_t n_bonds( 0), n_e_in_valence_bonds( 0);
      for( size_t ring_index( 0), ring_sz( RING.GetSize()); ring_index < ring_sz; ++ring_index)
      {
        for
        (
          auto itr_bnd( FRAG.GetAtomVector()( RING( ring_index)).GetBonds().Begin()),
               itr_bnd_end( FRAG.GetAtomVector()( RING( ring_index)).GetBonds().End());
          itr_bnd != itr_bnd_end;
          ++itr_bnd
        )
        {
          if( itr_bnd->GetBondType()->IsBondInRing() || itr_bnd->GetTargetAtom().GetAtomType() == GetAtomTypes().H_S)
          {
            continue;
          }
          if
          (
            ++n_bonds > size_t( 4)
            ||
            ( n_e_in_valence_bonds += itr_bnd->GetBondType()->GetNumberOfElectrons() / size_t( 2)) > size_t( 4)
          )
          {
            return false;
          }
        }
      }
      return true;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    io::Serializer FragmentMutateRingSwap::GetSerializer() const
    {
      io::Serializer parameters( FragmentMutateInterface::GetSerializer());
      parameters.SetClassDescription( "Swaps ring systems; can also create and delete ring systems");
      parameters.AddInitializer
      (
        "restricted",
        "restricts the types of ring swaps that are allowed to those that change at most one ring's size. This also implies "
        "that if rings are grown from single atoms, it will be a simple closed ring, not a system. Also that at most one "
        "ring in a ring system will change size (or be removed/expanded",
        io::Serialization::GetAgent( &m_RestrictToNoMoreThanOneRingSizeChange),
        "false"
      );

      parameters.AddInitializer
      (
        "conservative",
        "aligns new ring to current ring prior to substitution to preserve topological distances between substituents "
        "as best as possible; useful for focused re-design of an existing scaffold. Set to False to increase the diversity "
        "of new scaffolds created.",
        io::Serialization::GetAgent( &m_AlignRings),
        "true"
      );

      parameters.AddInitializer
      (
        "fix_geometry",
        "If True, then any atom/bonds with bad geometry is included for conformational sampling. If False, "
        "then atoms with bad geometry will not be included unless they are also one of the perturbed atoms or "
        "included as adjacent to the perturbed atoms.",
        io::Serialization::GetAgent( &m_FixGeometry),
        "true"
      );

      parameters.AddInitializer
      (
        "refine_alignment",
        "If True, then choose the returned conformer based on a flexible substructure-based alignment scored with ChargeRMSD. "
        "This method generates a conformational ensemble, performs a greedy disconnected substructure alignment of each conformer, "
        "and then chooses the best one by ChargeRMSD score. If False, select the best conformer based on BCL::Conf score. "
        "This option will reduce the speed of the mutate and is mostly recommended for pose-dependent replacement of ring "
        "structures at the core of the molecule.",
        io::Serialization::GetAgent( &m_ChooseBestAlignedConf),
        "false"
      );

      parameters.AddInitializer
      (
        "ring_initiation_probability",
        "Probability of initiating / removing a ring vs just swapping rings by default. So at 0.1, new rings will be "
        "created from individual atoms or removed 10% of the time, and just swapped 90% of the time",
        io::Serialization::GetAgent( &m_RingInitiationProbability),
        "0.1"
      );

      parameters.AddInitializer
      (
        "atom_comparison",
        "atom data that is compared to determine whether atoms are equivalent",
        io::Serialization::GetAgent( &m_AtomComparisonType),
        "CouldHaveSubstituents"
      );
      parameters.AddInitializer
      (
        "bond_comparison",
        "bond data that is compared",
        io::Serialization::GetAgent( &m_BondComparisonType),
        "Identity"
      );

      parameters.AddInitializer
      (
        "ring_library",
        "path to the ring library file",
        io::Serialization::GetAgent( &m_RingLibraryFilename),
        RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "ring_libraries/drug_ring_database.simple.sdf.gz"
      );

      parameters.AddInitializer
      (
        "extend_adjacent_atoms",
        "include adjacent atoms out this many bonds from any perturbed atom when generating a new 3D conformer",
        io::Serialization::GetAgent( &m_ExtendAdjacentAtoms),
        "1"
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentMutateRingSwap::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      // static initialization check
      if( command::CommandState::IsInStaticInitialization())
      {
        return true;
      }

      // call RISH function of the base class
      if( !FragmentMutateInterface::ReadInitializerSuccessHook( LABEL, ERROR_STREAM))
      {
        return false;
      }

      // read in ring library
      io::IFStream input;
      io::File::MustOpenIFStream
      (
        input,
        !m_RingLibraryFilename.empty() ?
        m_RingLibraryFilename : RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "ring_libraries/drug_ring_database.simple.sdf.gz"
      );
      FragmentEnsemble ensemble( input, sdf::e_Remove);
      io::File::CloseClearFStream( input);
      for
      (
         auto itr_ensemble( ensemble.Begin()), itr_ensemble_end( ensemble.End());
         itr_ensemble != itr_ensemble_end;
         ++itr_ensemble
      )
      {
        size_t n_dbv( GetNumberDoubleBondValences( *itr_ensemble));
        if( n_dbv + size_t( 1) > m_FragmentPool.GetSize())
        {
          m_FragmentPool.Resize( n_dbv + size_t( 1));
          m_FragmentPoolScaffoldSums.Resize( n_dbv + size_t( 1), size_t( 0));
        }
        m_FragmentPool( n_dbv).PushBack( *itr_ensemble);
        m_FragmentPoolScaffoldSums( n_dbv) += GetCounts( *itr_ensemble);
      }

      // we require a ring library for this mutate to function
      BCL_Assert( m_FragmentPool( 0).GetSize(), "No valid rings in database! Check number of double bond valences!");

      // done
      return true;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace chemistry
} // namespace bcl
