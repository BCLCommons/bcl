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
#include "chemistry/bcl_chemistry_small_molecule_fragment_mapping.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_sh_ptr_vector.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SmallMoleculeFragmentMapping::s_Instance
    (
      GetObjectInstances().AddInstance( new SmallMoleculeFragmentMapping())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new SmallMoleculeFragmentMapping
    SmallMoleculeFragmentMapping *SmallMoleculeFragmentMapping::Clone() const
    {
      return new SmallMoleculeFragmentMapping( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SmallMoleculeFragmentMapping::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operators //
  ////////////////

    //! @brief get all possible mappings between fragments identified and the molecule of interest
    //! @param FRAGMENTS fragment isomorphisms of fragments contained in a molecule of interest
    //! @return map of number of rotatable bond of fragment to its RotateDihedralBondData
    util::ShPtrVector< RotamerDihedralBondData> SmallMoleculeFragmentMapping::MapFragmentIsomorphisms
    (
      const FragmentComplete &MOLECULE,
      const util::ShPtrVector< SmallMoleculeFragmentIsomorphism> &FRAGMENTS,
      const storage::Set< size_t> &SAMPLE_PARTS
    )
    {
      size_t n_skipped( 0), n_skipped_better( 0);
      // initialize a vector to store RotamerDihedralBondData objects
      util::ShPtrVector< RotamerDihedralBondData> fragment_interface_vector;

      // get bond info of the molecule of interest
      storage::Vector< sdf::BondInfo> bond_info( MOLECULE.GetBondInfo());

      storage::Set< graph::UndirectedEdge< size_t> > ring_bonds;
      storage::Set< graph::UndirectedEdge< size_t> > all_bonds;
      for( auto itr( bond_info.Begin()), itr_end( bond_info.End()); itr != itr_end; ++itr)
      {
        if( itr->GetConfigurationalBondType()->IsBondInRing())
        {
          ring_bonds.Insert( graph::UndirectedEdge< size_t>( itr->GetAtomIndexLow(), itr->GetAtomIndexHigh(), size_t( 0)));
        }
        all_bonds.Insert( graph::UndirectedEdge< size_t>( itr->GetAtomIndexLow(), itr->GetAtomIndexHigh(), size_t( 0)));
      }

      // counts of different central bonds
      storage::Map
      <
        util::ShPtr< RotamerDihedralBondData>,
        util::SiPtr< const storage::Set< graph::UndirectedEdge< size_t> > >
      > central_bond_infos;
      storage::Map
      <
        storage::Set< graph::UndirectedEdge< size_t> >,
        storage::Vector
        <
          storage::Pair< util::ShPtr< SmallMoleculeFragmentIsomorphism>, storage::List< storage::Vector< size_t> > >
        >
      > central_bonds_to_best_isomorphisms;
      storage::Map
      <
        storage::Set< graph::UndirectedEdge< size_t> >,
        storage::Vector
        <
          storage::Pair< util::ShPtr< SmallMoleculeFragmentIsomorphism>, storage::List< storage::Vector< size_t> > >
        >
      > central_bonds_to_best_isomorphisms_flexible_rings;

      // create appropriate RotamerDihedralBondData by iterating over all fragments of the molecule
      storage::Vector< double> energies;
      for
      (
        util::ShPtrVector< SmallMoleculeFragmentIsomorphism>::const_iterator
          itr_iso( FRAGMENTS.Begin()), itr_iso_end( FRAGMENTS.End());
        itr_iso != itr_iso_end;
        ++itr_iso
      )
      {
        const util::ShPtr< SmallMoleculeFragmentIsomorphism> &sp_iso( *itr_iso);
//        if( sp_iso->ContainsRings() && ( *itr_iso)->GetNumberChainBonds() && !( *itr_iso)->GetNumberDihedralChainBonds())
//        {
//          continue;
//        }
        if( !sp_iso->ContainsRings() && !( *itr_iso)->GetNumberDihedralChainBonds())
        {
          continue;
        }
//        if( sp_iso->GetFragment().GetNumberHydrogens())
//        {
//          continue;
//        }
//        if( sp_iso->GetFragmentCounts() < double( 4.0))
//        {
//          continue;
//        }

        // > 2 ring systems - too many degrees of freedom. Assume that the signal is pairwise decomposable
        if( ( *itr_iso)->GetNumberRingSystems() > size_t( 2))
        {
          continue;
        }
        // non-aromatic rings don't tend to interact strongly (e.g. no pi-stacking, etc.) so
        // ignore extended fragments with them
        if
        (
          ( *itr_iso)->GetNumberDihedralChainBonds() > size_t( 1)
          && ( *itr_iso)->GetNumberAromaticRingSystems() != ( *itr_iso)->GetNumberRingSystems()
        )
        {
          continue;
        }

        // excessive degrees of freedom, likely no additional information over smaller fragments
        if( ( *itr_iso)->GetNumberDihedralChainBonds() > size_t( 4))
        {
          continue;
        }
        if( ( *itr_iso)->GetNumberDihedralChainBonds() > size_t( 1))
        {
          bool has_terminal_h( false);
          for( auto itr_atom( ( *itr_iso)->GetFragment().GetAtomsIterator()); itr_atom.NotAtEnd(); ++itr_atom)
          {
            const size_t nh( itr_atom->GetNumberCovalentlyBoundHydrogens());
            if( nh && itr_atom->GetBonds().GetSize() <= nh + size_t( 1))
            {
              has_terminal_h = true;
              break;
            }
          }
          if( has_terminal_h)
          {
            continue;
          }
        }

//        // chain systems separated by a ring; typically only correlated to avoid clashes or due to chance
//        if( ( *itr_iso)->GetNumberChainSystems() > size_t( 2))
//        {
//          continue;
//        }
//        // 4 dihedral chain bonds and no aromatic rings - no possibility for pi stacking
//        if( ( *itr_iso)->GetNumberDihedralChainBonds() > size_t( 1) && !( *itr_iso)->GetNumberAromaticRingSystems())
//        {
//          continue;
//        }

        // get priority dihedral angles of the fragment of interest
        const storage::Vector< storage::VectorND< 4, size_t> > &
          priority_angle_info( sp_iso->GetDihedralAngleIndices());
        if( priority_angle_info.IsEmpty())
        {
          continue;
        }

        // get reference to isomorphisms of the fragment
        const storage::Vector< storage::Vector< size_t> > &isomorphisms( sp_iso->GetFragmentCSI());
        size_t iso_number( 0);
        for
        (
          storage::Vector< storage::Vector< size_t> >::const_iterator
            itr_isomorphisms( isomorphisms.Begin()), itr_isomorphisms_end( isomorphisms.End());
          itr_isomorphisms != itr_isomorphisms_end;
          ++itr_isomorphisms, ++iso_number
        )
        {
          // get isomorphism between fragment and the molecule of interest
          const storage::Vector< size_t> &isomorphism( *itr_isomorphisms);
          if( !SAMPLE_PARTS.IsEmpty())
          {
            const storage::Set< size_t> isomorphism_set( isomorphism.Begin(), isomorphism.End());
            bool should_sample( true);
            // check for whether all dihedral atoms are contained in SAMPLE_PARTS
            for
            (
                auto itr_iso_set( isomorphism_set.Begin()), itr_iso_set_end( isomorphism_set.End());
                itr_iso_set != itr_iso_set_end;
                ++itr_iso_set
            )
            {
              if
              (
                  !SAMPLE_PARTS.Contains( *itr_iso_set)
                  && MOLECULE.GetAtomVector()( *itr_iso_set).GetBonds().GetSize() > size_t( 1)
              )
              {
                should_sample = false;
                break;
              }
            }
            if( !should_sample)
            {
              continue;
            }
          }

          if( sp_iso->ContainsIncompleteRings( iso_number) && !( *itr_iso)->GetNumberDihedralChainBonds())
          {
            continue;
          }

          storage::Vector< graph::UndirectedEdge< size_t> > chain_dihedrals;
          for( auto itr_b( priority_angle_info.Begin()), itr_b_end( priority_angle_info.End()); itr_b != itr_b_end; ++itr_b)
          {
            graph::UndirectedEdge< size_t> n_edge( isomorphism( itr_b->Second()), isomorphism( itr_b->Third()), size_t( 0));
            if( ring_bonds.Contains( n_edge))
            {
              continue;
            }
            chain_dihedrals.PushBack( n_edge);
          }
//          if( chain_dihedrals.GetSize() == size_t( 2))
//          {
//            if
//            (
//              all_bonds.Contains( graph::UndirectedEdge< size_t>( chain_dihedrals( 0).GetIndexLow(), chain_dihedrals( 1).GetIndexLow(), size_t( 0)))
//              || all_bonds.Contains( graph::UndirectedEdge< size_t>( chain_dihedrals( 0).GetIndexLow(), chain_dihedrals( 1).GetIndexHigh(), size_t( 0)))
//              || all_bonds.Contains( graph::UndirectedEdge< size_t>( chain_dihedrals( 0).GetIndexHigh(), chain_dihedrals( 1).GetIndexLow(), size_t( 0)))
//              || all_bonds.Contains( graph::UndirectedEdge< size_t>( chain_dihedrals( 0).GetIndexHigh(), chain_dihedrals( 1).GetIndexHigh(), size_t( 0)))
//            )
//            {
//            }
//            else
//            {
//              continue;
//            }
//          }
//          else if( chain_dihedrals.GetSize() > size_t( 2))
//          {
//            continue;
//          }

          // create a set to store center bonds of rotatable bonds which fragment represents in a molecule
          storage::Set< graph::UndirectedEdge< size_t> > center_bonds;

          // get atom indices from the frame of reference of the molecule of interest
          for
          (
            storage::Vector< storage::VectorND< 4, size_t> >::const_iterator
              itr( priority_angle_info.Begin()), itr_end( priority_angle_info.End());
            itr != itr_end;
            ++itr
          )
          {
            // get center bonds of the fragment in terms of molecule atom indices to get mapping between fragment and molecule
            graph::UndirectedEdge< size_t> new_edge( isomorphism( itr->Second()), isomorphism( itr->Third()), size_t( 0));
            center_bonds.InsertElement( new_edge);
          }

          auto &cbbi_vec
          (
            ( *itr_iso)->ContainsRings() && ( *itr_iso)->GetNumberDihedralChainBonds() == size_t( 0)
            ? central_bonds_to_best_isomorphisms_flexible_rings[ center_bonds]
            : central_bonds_to_best_isomorphisms[ center_bonds]
          );
//          if( !( *itr_iso)->ContainsRings() && !cbbi_vec.IsEmpty())
//          {
//            const size_t cbbi_chain( cbbi_vec.FirstElement().First()->GetNumberChainBonds());
//            const size_t cbbi_ring( cbbi_vec.FirstElement().First()->GetNumberRingBonds());
//            if( cbbi_vec.FirstElement().First()->GetExpectedFreeEnergy( false).GetStandardDeviation() > sp_iso->GetExpectedFreeEnergy( false).GetStandardDeviation())
//            {
//              BCL_Assert
//              (
//                !( *itr_iso)->ContainsRings() || ( *itr_iso)->GetNumberChainBonds(),
//                "Multiple disparate isomorphisms found for the same ring system. Fragment library corruption?"
//              );
//              BCL_MessageStd( "Discarding molecule due to molecule with better energy: ");
//              ( *itr_iso)->GetFragment().WriteMDL( util::GetLogger());
//              ++n_skipped_better;
//              continue;
//            }
//            else if
//            (
//              cbbi_vec.FirstElement().First()->GetExpectedFreeEnergy( false).GetStandardDeviation() < sp_iso->GetExpectedFreeEnergy( false).GetStandardDeviation()
//            )
//            {
//              BCL_Assert
//              (
//                !( *itr_iso)->ContainsRings() || ( *itr_iso)->GetNumberChainBonds(),
//                "Multiple disparate isomorphisms found for the same ring system. Fragment library corruption?"
//              );
//              n_skipped_better += cbbi_vec.GetSize();
//              BCL_MessageStd( "Discarding molecule due to molecule with better energy: ");
//              cbbi_vec.FirstElement().First()->GetFragment().WriteMDL( util::GetLogger());
//              cbbi_vec.Reset();
//            }
//          }
          bool found( false);
          for( auto itr_cbbi( cbbi_vec.Begin()), itr_cbbi_end( cbbi_vec.End()); itr_cbbi != itr_cbbi_end; ++itr_cbbi)
          {
            if( itr_cbbi->First() == *itr_iso)
            {
              found = true;
              if( std::find( itr_cbbi->Second().Begin(), itr_cbbi->Second().End(), isomorphism) == itr_cbbi->Second().End())
              {
                itr_cbbi->Second().PushBack( isomorphism);
              }
              break;
            }
          }
          if( !found)
          {
            cbbi_vec.PushBack
            (
              storage::Pair< util::ShPtr< SmallMoleculeFragmentIsomorphism>, storage::List< storage::Vector< size_t> > >
              (
                *itr_iso,
                storage::List< storage::Vector< size_t> >( size_t( 1), isomorphism)
              )
            );

          }
        }
      }

      for
      (
        auto itr_best( central_bonds_to_best_isomorphisms.Begin()), itr_best_end( central_bonds_to_best_isomorphisms.End());
        itr_best != itr_best_end;
        ++itr_best
      )
      {
        // create a set to store center bonds of rotatable bonds which fragment represents in a molecule
        const storage::Set< graph::UndirectedEdge< size_t> > &center_bonds( itr_best->first);
        // going over each isomorphism with the same central bonds
        for( auto itr_iso( itr_best->second.Begin()), itr_iso_end( itr_best->second.End()); itr_iso != itr_iso_end; ++itr_iso)
        {
          fragment_interface_vector.PushBack
          (
            util::ShPtr< RotamerDihedralBondData>
            (
              new RotamerDihedralBondData( MOLECULE, *itr_iso->First(), itr_iso->Second())
            )
          );
          central_bond_infos[ fragment_interface_vector.LastElement()] = util::ToSiPtr( center_bonds);
        }
      }

      for
      (
        auto itr_best( central_bonds_to_best_isomorphisms_flexible_rings.Begin()), itr_best_end( central_bonds_to_best_isomorphisms_flexible_rings.End());
        itr_best != itr_best_end;
        ++itr_best
      )
      {
        // create a set to store center bonds of rotatable bonds which fragment represents in a molecule
        const storage::Set< graph::UndirectedEdge< size_t> > &center_bonds( itr_best->first);

        // going over each isomorphism with the same central bonds
        for( auto itr_iso( itr_best->second.Begin()), itr_iso_end( itr_best->second.End()); itr_iso != itr_iso_end; ++itr_iso)
        {
          fragment_interface_vector.PushBack
          (
            util::ShPtr< RotamerDihedralBondData>
            (
              new RotamerDihedralBondData( MOLECULE, *itr_iso->First(), itr_iso->Second())
            )
          );
          central_bond_infos[ fragment_interface_vector.LastElement()] = util::ToSiPtr( center_bonds);
        }
      }

      // pre-filter useless fragments
      auto itr_place( fragment_interface_vector.Begin());
      for
      (
        auto itr_next( fragment_interface_vector.Begin()), itr_end( fragment_interface_vector.End());
        itr_next != itr_end;
        ++itr_next
      )
      {
        if( !( *itr_next)->ContainsRings() && !( *itr_next)->GetFragment().GetNumberDihedralChainBonds())
        {
          continue;
        }
        if( ( *itr_next)->HasIncompleteRings() && !( *itr_next)->GetFragment().GetNumberDihedralChainBonds())
        {
          continue;
        }
        if( itr_place != itr_next)
        {
          *itr_place = *itr_next;
        }
        ++itr_place;
      }
      const size_t new_rot_bond_data_size( size_t( std::distance( fragment_interface_vector.Begin(), itr_place)));
      if( new_rot_bond_data_size != fragment_interface_vector.GetSize())
      {
        fragment_interface_vector.Resize( new_rot_bond_data_size);
      }

      return fragment_interface_vector;
    }

    //! @brief Filter out RDBD, choosing the best one to represent each case
    //! @param FRAGMENTS fragment isomorphisms of fragments contained in a molecule of interest
    //! @return map of number of rotatable bond of fragment to its RotateDihedralBondData
    util::ShPtrVector< RotamerDihedralBondData> SmallMoleculeFragmentMapping::Filter
    (
      const util::ShPtrVector< RotamerDihedralBondData> &ROTAMERS
    )
    {
      util::ShPtrVector< RotamerDihedralBondData> rdbd;

      storage::Map< storage::Set< size_t>, util::SiPtr< const RotamerDihedralBondData> > all_atoms_to_rdbd;
      storage::Set< util::SiPtr< const RotamerDihedralBondData> > to_ignore;
      for( auto itr( ROTAMERS.Begin()), itr_end( ROTAMERS.End()); itr != itr_end; ++itr)
      {
        storage::Set< size_t> atom_ids;
        for
        (
          auto itr_iso( ( *itr)->GetIsomorphisms().Begin()), itr_iso_end( ( *itr)->GetIsomorphisms().End());
          itr_iso != itr_iso_end;
          ++itr_iso
        )
        {
          atom_ids = atom_ids + storage::Set< size_t>( itr_iso->Begin(), itr_iso->End());
        }
        util::SiPtr< const RotamerDihedralBondData> &match( all_atoms_to_rdbd[ atom_ids]);
        if( match.IsDefined())
        {
          BCL_MessageStd
          (
            "Test: " + util::Format()( match->GetFragment().GetFragment().GetSize())
            + " " + util::Format()( ( *itr)->GetFragment().GetFragment().GetSize())
            + " " + util::Format()( match->GetIsomorphisms().GetSize())
            + " " + util::Format()( ( *itr)->GetIsomorphisms().GetSize())
            + " " + util::Format()( match->GetRotamers().GetSize())
            + " " + util::Format()( ( *itr)->GetRotamers().GetSize())
            + " " + util::Format()( match->GetFragment().GetExpectedFreeEnergy( false).GetStandardDeviation())
            + " " + util::Format()( ( *itr)->GetFragment().GetExpectedFreeEnergy( false).GetStandardDeviation())
          );
          if( match->GetFragment().GetFragment().GetSize() > ( *itr)->GetFragment().GetFragment().GetSize())
          {
            to_ignore.Insert( *match);
            match = *itr;
          }
          else
          {
            to_ignore.Insert( **itr);
          }
        }
        else
        {
          match = *itr;
        }
      }
      if( to_ignore.GetSize())
      {
        BCL_MessageStd( "Number of fragments filtered: " + util::Format()( to_ignore.GetSize()));
      }

      for( auto itrb( ROTAMERS.Begin()), itr_end( ROTAMERS.End()); itrb != itr_end; ++itrb)
      {
        if( !to_ignore.Contains( *itrb))
        {
          rdbd.PushBack( *itrb);
        }
      }

      return ROTAMERS;
    }

    //! @brief get all possible rings contained in the molecule of  interest
    //! @param MOLECULE molecule whose rings have to be found
    //! @return list of atom vertices that are contained in rings
    const storage::List< storage::Vector< size_t> > SmallMoleculeFragmentMapping::GetAllRings
    (
      const ConformationInterface &MOLECULE
    )
    {
      ConformationGraphConverter graph_maker( ConformationGraphConverter::e_AtomType, ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness);

      // create a graph
      graph::ConstGraph< size_t, size_t> mol_graph( graph_maker( MOLECULE));

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
        if( !itr_bond->GetConstitutionalBondType()->IsBondInRing())
        {
          mol_graph.RemoveEdge( itr_bond->GetAtomIndexLow(), itr_bond->GetAtomIndexHigh());
        }
      }

      // get connected components of the graph, which are either isolated atoms or rings
      storage::List< storage::Vector< size_t> > components( graph::Connectivity::GetComponents( mol_graph));

      storage::List< storage::Vector< size_t> > ring_vertices;
      for
      (
        storage::List< storage::Vector< size_t> >::const_iterator itr( components.Begin()), itr_end( components.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->GetSize() > 2)
        {
          ring_vertices.PushBack( *itr);
        }
      }
      return ring_vertices;
    }

    //! @brief get all possible rings contained in the molecule of  interest
    //! @param MOLECULE molecule whose rings have to be found
    //! @return list of atom vertices that are contained in rings
    const storage::List< storage::Vector< size_t> > SmallMoleculeFragmentMapping::GetRingVertices
    (
      const ConformationInterface &MOLECULE
    )
    {
      storage::List< storage::Vector< size_t> > ring_vertices;
      // count ring bonds
      const size_t number_ring_bonds
      (
        MOLECULE.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1))
      );

      // skip molecules that have no ring bonds
      if( number_ring_bonds >= size_t( 3))
      {
        if( number_ring_bonds == MOLECULE.GetNumberBonds())
        {
          ring_vertices.PushBack( storage::CreateIndexVector( MOLECULE.GetNumberAtoms()));
        }
        else
        {
          ring_vertices = GetAllRings( MOLECULE);
        }
      }
      return ring_vertices;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SmallMoleculeFragmentMapping::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &SmallMoleculeFragmentMapping::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl

