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
#include "chemistry/bcl_chemistry_fragment_probability_score.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_comparison_by_dihedral_bins.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "chemistry/bcl_chemistry_priority_dihedral_angles.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_running_min_max.h"

// external includes - sorted alphabetically

//#define BCL_PROFILE_FragmentProbabilityScore
#ifdef BCL_PROFILE_FragmentProbabilityScore
#include "util/bcl_util_stopwatch.h"
#endif

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> FragmentProbabilityScore::s_Instance
    (
      GetObjectInstances().AddInstance( new FragmentProbabilityScore())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////
    FragmentProbabilityScore::FragmentProbabilityScore() :
      m_NumberDihedralBonds( 0),
      m_DihedralBinComparer( PriorityDihedralAngles::GetWrappingAngle() * 2.0)
    {
    }

    //! @brief constructor
    //! @param FRAGMENT_DATA vector of fragments that are part of the molecule of interest
    //! @param PRIORITY_ANGLES priority dihedral angle object for molecule whose conformations are being sampled
    //! @param FRAGMENT_SIZE_WT weight to be given to fragment size
    //! @param FRAGMENT_COUNTS_WT weight to be given to fragment counts
    //! @param ROTAMER_PROPENSITY_WT weight to be given to rotamer propensity
    FragmentProbabilityScore::FragmentProbabilityScore
    (
      const util::ShPtrVector< RotamerDihedralBondData> &FRAGMENT_DATA,
      const bool &CONSIDER_ISOMETRY_CHANGE
    )
    :
      m_RotamerData( FRAGMENT_DATA),
      m_NumberDihedralBonds( 0),
      m_ConsiderIsometryChange( CONSIDER_ISOMETRY_CHANGE)
    {
      // iterate through each fragment
      storage::Vector< graph::UndirectedEdge< ConfigurationalBondType> > dihedral_edges;
      storage::List< storage::Vector< size_t> > ring_list;
      size_t n_atoms( 0);
      for
      (
        util::ShPtrVector< RotamerDihedralBondData>::const_iterator itr( m_RotamerData.Begin()), itr_end( m_RotamerData.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr == m_RotamerData.Begin())
        {
          m_DihedralBinComparer = ConformationComparisonByDihedralBins( ( *itr)->GetBinSize());
          auto graph( ConformationGraphConverter::CreateGraphWithAtoms( ( *itr)->GetFragment().GetMolecule()));
          ring_list = FragmentSplitRings( true, 3).GetComponentVertices( ( *itr)->GetFragment().GetMolecule(), graph);
          dihedral_edges = PriorityDihedralAngles::GetDihedralEdges( ( *itr)->GetFragment().GetMolecule());
          m_NumberDihedralBonds = dihedral_edges.GetSize();
          n_atoms = ( *itr)->GetFragment().GetMolecule().GetSize();
        }
        else
        {
          BCL_Assert
          (
            m_DihedralBinComparer.GetBinSize() == ( *itr)->GetBinSize(),
            "Different bin sizes used in fragment library"
          );
        }
      }
      m_DihedralComponentSize.Resize( m_NumberDihedralBonds);
      m_DihedralComponentSize.SetAllElements( size_t( 1));
      storage::Vector< size_t> main_component_size( n_atoms, size_t( 1));
      for( auto itr_ring_list( ring_list.Begin()), itr_ring_list_end( ring_list.End()); itr_ring_list != itr_ring_list_end; ++itr_ring_list)
      {
        for( auto itr_ring( itr_ring_list->Begin()), itr_ring_end( itr_ring_list->End()); itr_ring != itr_ring_end; ++itr_ring)
        {
          main_component_size( *itr_ring) = itr_ring_list->GetSize();
        }
      }
      size_t dihedral_n( 0);
      for( auto itr_dihedral( dihedral_edges.Begin()), itr_dihedral_end( dihedral_edges.End()); itr_dihedral != itr_dihedral_end; ++itr_dihedral, ++dihedral_n)
      {
        if( itr_dihedral->GetEdgeData()->IsBondInRing())
        {
          m_DihedralComponentSize( dihedral_n) = main_component_size( itr_dihedral->GetIndexLow());
        }
      }
      // remove irrelevant rotamers
      auto itr_place( m_RotamerData.Begin());
      for
      (
        util::ShPtrVector< RotamerDihedralBondData>::iterator
          itr_next( m_RotamerData.Begin()), itr_end( m_RotamerData.End());
        itr_next != itr_end;
        ++itr_next
      )
      {
        // for pure ring fragments
        if( !( *itr_next)->GetFragment().ContainsRingConformations() && !( *itr_next)->GetFragment().GetNumberDihedralChainBonds())
        {
          // no ring conformations and no dihedrals -> nothing to score
          continue;
        }
        if( ( *itr_next)->GetFragment().ContainsOnlyIncompleteHydrogenation())
        {
          // lacking full set of hydrogens. Most likely there are many isomorphisms too (due to hydrogen equivalence)
          // Generally these fragments are very expensive to score and provide no additional information that isn't available
          // from the fully-hydrogenated variants, so better to skip them
          continue;
        }

        if( itr_place != itr_next)
        {
          *itr_place = *itr_next;
        }
        ++itr_place;
      }
      m_RotamerData.Resize( std::distance( m_RotamerData.Begin(), itr_place));
      BCL_MessageStd
      (
        "Number of fragments in rotamer library that map to this molecule: "
        + util::Format()( m_RotamerData.GetSize())
      );
    }

    //! @brief Clone function
    //! @return pointer to new FragmentProbabilityScore
    FragmentProbabilityScore *FragmentProbabilityScore::Clone() const
    {
      return new FragmentProbabilityScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentProbabilityScore::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const util::ShPtrVector< RotamerDihedralBondData> &FragmentProbabilityScore::GetRotamerData() const
    {
      return m_RotamerData;
    }

    //! @brief get the dihedral component sizes used to normalize per-dihedral scores
    const storage::Vector< size_t> &FragmentProbabilityScore::GetDihedralComponentSize() const
    {
      return m_DihedralComponentSize;
    }

    //! @brief get the number of dihedrals
    const size_t FragmentProbabilityScore::GetNumberDihedralBonds() const
    {
      return m_NumberDihedralBonds;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief evaluate clashes for given atom pair
    //! @param MOLECULE molecule that needs to scored
    //! @return propensity score for observing rotamers that exist in conformation
    double FragmentProbabilityScore::operator()
    (
      const FragmentComplete &MOLECULE
    ) const
    {

      #ifdef BCL_PROFILE_FragmentProbabilityScore
      static util::Stopwatch s_timer( "fragment_probability_score", util::Message::e_Standard, true);
      s_timer.Start();
      #endif

      // get the individual dihedral scores
      storage::Vector< double> individual_dihedral_scores( GetDihedralScores( MOLECULE));

      // take the average dihedral score
      math::RunningAverage< double> val;
      for( size_t dhb_index( 0); dhb_index < m_NumberDihedralBonds; ++dhb_index)
      {
        // checks for 0-weighted positions to ignore them
        if( individual_dihedral_scores( dhb_index))
        {
          val.AddWeightedObservation( individual_dihedral_scores( dhb_index), 1.0 / double( m_DihedralComponentSize( dhb_index)));
        }
      }
      double temp_score( val.GetAverage());
      temp_score += MOLECULE.GetTotalAmideBondNonPlanarity( 10.0) / 100.0
                    + MOLECULE.GetTotalAmideBondNonPlanarity( 15.0) / 10.0
                    + MOLECULE.GetTotalAmideBondNonPlanarity( 25.0);

      #ifdef BCL_PROFILE_FragmentProbabilityScore
      s_timer.Stop();
      #endif
      return temp_score;
    }

    //! @brief evaluate 1-4 interaction score
    //! @param MOLECULE molecule that needs to scored
    //! @return 1-4 interaction score
    double FragmentProbabilityScore::Get14InteractionScore
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      return m_14Score( MOLECULE);
    }

    //! @brief get back each component of the score as a separate molecule with full details/energies/etc
    FragmentEnsemble FragmentProbabilityScore::GetScoreComponents
    (
      const FragmentComplete &FRAG_COMP
    ) const
    {
      storage::Vector< math::RunningAverage< double> > dihedral_scores( m_NumberDihedralBonds);
      FragmentEnsemble ensemble;

      t_Map dihedral_bin_hash_to_bin;
      t_Map dihedral_bin_hash_to_value;

      // iterate through each fragment of input molecule conformer
      for
      (
        util::ShPtrVector< RotamerDihedralBondData>::const_iterator itr( m_RotamerData.Begin()), itr_end( m_RotamerData.End());
        itr != itr_end;
        ++itr
      )
      {
        if( ( *itr)->GetRotamerBonds().IsEmpty())
        {
          continue;
        }
        //       const double best_score
        //       (
        //         m_ConsiderIsometryChange
        //         ? ( *itr)->GetFragment().GetBestFreeEnergy()
        //         : ( *itr)->GetFragment().GetBestFreeEnergyNoIsometry()
        //       );

        // calculate background probability of observing a rotamer
        // number of rotamers of current fragment divided by the number of times this fragment appears in CSD/PDB
        const storage::Vector< storage::Vector< size_t> > &center_bond_indices_v( ( *itr)->GetCenterBondIsomorphisms());

        FragmentComplete copy( ( *itr)->GetFragment().GetFragment());
        copy.StoreProperty( "NIsomorphisms", util::Format()( ( *itr)->GetIsomorphisms().GetSize()));
        copy.StoreProperty
        (
          "NDihedralChains",
          util::Format()( ( *itr)->GetFragment().GetNumberDihedralChainBonds())
        );
        copy.StoreProperty
        (
          "NIsomorphisms",
          util::Format()( ( *itr)->GetIsomorphisms().GetSize())
        );
        copy.StoreProperty
        (
          "NRingBonds",
          util::Format()( ( *itr)->GetFragment().GetNumberRingBonds())
        );
        copy.StoreProperty
        (
          "EstimatedMaxRotamers",
          util::Format()( ( *itr)->GetFragment().GetEstimatedMaxRotamers())
        );
        copy.StoreProperty
        (
          "NRotamers",
          util::Format()( ( *itr)->GetFragment().GetRotamerNumbers())
        );
        copy.StoreProperty
        (
          "NRotBond",
          util::Format()( ( *itr)->GetFragment().GetRotatableBondNumber())
        );
        copy.StoreProperty( "UpdatedCounts", ( *itr)->GetRotamerCounts());
        storage::Vector< int> bins;
        for
        (
          auto itr_bins( ( *itr)->GetFragment().GetRotamerBins().Begin()),
               itr_bins_end( ( *itr)->GetFragment().GetRotamerBins().End());
          itr_bins != itr_bins_end;
          ++itr_bins
        )
        {
          bins.Append( storage::Vector< int>( itr_bins->Begin(), itr_bins->End()));
        }
        copy.StoreProperty( "UpdatedBins", bins);

        math::RunningAverage< double> score_avg;
        // now find the fragment rotamers that exist in the given conformation
        auto itr_center_bond_indices( center_bond_indices_v.Begin());
        auto itr_center_bond_indices_best( itr_center_bond_indices);
        auto itr_weights( ( *itr)->GetIsomorphismWeights().Begin());
        size_t isonumber( 0);
        for
        (
          auto itr_signature( ( *itr)->GetRotamerBonds().Begin()), itr_signature_end( ( *itr)->GetRotamerBonds().End());
          itr_signature != itr_signature_end;
          ++itr_signature, ++itr_center_bond_indices, ++itr_weights, ++isonumber
        )
        {
          auto dihedral_bins_and_angles
          (
            GetMolecularBinFromFragment
            (
              FRAG_COMP,
              *itr_signature,
              dihedral_bin_hash_to_bin
            )
          );
          double local_score
          (
            ( *itr)->GetFragment().GetRotamerFreeEnergy
            (
              dihedral_bins_and_angles.first,
              m_ConsiderIsometryChange
            )
          );
          copy.StoreProperty( "DihedralScore" + util::Format()( isonumber), util::Format()( local_score));
          copy.StoreProperty( "CenterBondIndices" + util::Format()( isonumber), *itr_center_bond_indices);
          copy.StoreProperty
          (
            "DihedralSignature" + util::Format()( isonumber),
            GetMolecularBinFromFragment
            (
              FRAG_COMP,
              *itr_signature,
              dihedral_bin_hash_to_bin
            ).first
          );
          copy.StoreProperty
          (
            "DihedralAngles" + util::Format()( isonumber),
            GetMolecularBinFromFragment
            (
              FRAG_COMP,
              *itr_signature,
              dihedral_bin_hash_to_bin
            ).second
          );
          copy.StoreProperty
          (
            "IsomerWeight" + util::Format()( isonumber),
            util::Format()( *itr_weights)
          );

          score_avg.AddWeightedObservation( local_score, *itr_weights);
        }
        const auto &center_bond_indices( *itr_center_bond_indices_best);
        if( ( ( *itr)->ContainsRings() && !( *itr)->GetFragment().GetNumberDihedralChainBonds()) || !( *itr)->ContainsRings())
        {
          for
          (
            size_t bond_nr_frag( 0), n_bonds_frag( center_bond_indices.GetSize());
            bond_nr_frag < n_bonds_frag;
            ++bond_nr_frag
          )
          {
            dihedral_scores( center_bond_indices( bond_nr_frag)).AddWeightedObservation( score_avg.GetAverage(), 1.0);
          }
        }
        else
        {
          for
          (
            size_t bond_nr_frag( 0), n_bonds_frag( center_bond_indices.GetSize());
            bond_nr_frag < n_bonds_frag;
            ++bond_nr_frag
          )
          {
            if( m_DihedralComponentSize( center_bond_indices( bond_nr_frag)) == size_t( 1))
            {
              dihedral_scores( center_bond_indices( bond_nr_frag)).AddWeightedObservation( score_avg.GetAverage(), 1.0);
            }
          }
        }
        copy.StoreProperty( "FinalDihedralScore", util::Format()( score_avg.GetAverage()));
        ensemble.PushBack( copy);
      }

      return ensemble;
    }

    //! @brief get back the weighted average dihedral scores prior to aggregation
    storage::Vector< double> FragmentProbabilityScore::GetDihedralScores
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      // initialize score objects
      storage::Vector< math::RunningAverage< double> > dihedral_scores( m_NumberDihedralBonds);
      storage::Vector< double> final_dihedral_scores( m_NumberDihedralBonds);
      t_Map dihedral_bin_hash_to_bin;
      t_Map dihedral_bin_hash_to_value;

      // iterate through each fragment of input molecule conformer
      for
      (
          util::ShPtrVector< RotamerDihedralBondData>::const_iterator itr( m_RotamerData.Begin()), itr_end( m_RotamerData.End());
          itr != itr_end;
          ++itr
      )
      {
        if( ( *itr)->GetRotamerBonds().IsEmpty())
        {
          continue;
        }

        const double unseen_score
        (
          ( *itr)->GetFragment().GetRotamerFreeEnergy( linal::Vector< int>(), m_ConsiderIsometryChange)
        );

        // calculate background probability of observing a rotamer
        // number of rotamers of current fragment divided by the number of times this fragment appears in CSD/PDB
        const storage::Vector< storage::Vector< size_t> > &center_bond_indices_v( ( *itr)->GetCenterBondIsomorphisms());
        double local_dihedral_score( 1000.0);

        math::RunningAverage< double> score_avg;
        // now find the fragment rotamers that exist in the given conformation
        auto itr_center_bond_indices( center_bond_indices_v.Begin());
        auto itr_center_bond_indices_best( itr_center_bond_indices);
        auto itr_weights( ( *itr)->GetIsomorphismWeights().Begin());
        size_t iso_number( 0), tested_number( 0);
        for
        (
            auto itr_signature( ( *itr)->GetRotamerBonds().Begin()), itr_signature_end( ( *itr)->GetRotamerBonds().End());
            itr_signature != itr_signature_end;
            ++itr_signature, ++itr_center_bond_indices, ++itr_weights, ++iso_number
        )
        {
          if( ( *itr)->GetFragment().ContainsIncompleteHydrogenation( iso_number))
          {
            continue;
          }
          if( ++tested_number > size_t( 16))
          {
            break;
          }
          double local_score
          (
            ( *itr)->GetFragment().GetRotamerFreeEnergy
            (
              GetMolecularBinFromFragment
              (
                MOLECULE,
                *itr_signature,
                dihedral_bin_hash_to_bin
              ).first,
              m_ConsiderIsometryChange
            )
          );
          if( local_score < local_dihedral_score)
          {
            local_dihedral_score = local_score;
            itr_center_bond_indices_best = itr_center_bond_indices;
          }
          if( local_score != unseen_score)
          {
            score_avg.AddWeightedObservation( local_score, *itr_weights);
          }
        }
        if( !score_avg.GetWeight())
        {
          score_avg.AddWeightedObservation( unseen_score, 1.0);
        }
        const auto &center_bond_indices( *itr_center_bond_indices_best);
        // dihedral scores for rings are given on a per-ring basis
        if( ( ( *itr)->ContainsRings() && !( *itr)->GetFragment().GetNumberDihedralChainBonds()) || !( *itr)->ContainsRings())
        {
          for
          (
              size_t bond_nr_frag( 0), n_bonds_frag( center_bond_indices.GetSize());
              bond_nr_frag < n_bonds_frag;
              ++bond_nr_frag
          )
          {
            dihedral_scores( center_bond_indices( bond_nr_frag)).AddWeightedObservation( score_avg.GetAverage(), 1.0 / double( n_bonds_frag));
          }
        }
        // dihedrals for non-rings are given on a per-dihedral basis
        else
        {
          for
          (
              size_t bond_nr_frag( 0), n_bonds_frag( center_bond_indices.GetSize());
              bond_nr_frag < n_bonds_frag;
              ++bond_nr_frag
          )
          {
            if( m_DihedralComponentSize( center_bond_indices( bond_nr_frag)) == size_t( 1))
            {
              dihedral_scores( center_bond_indices( bond_nr_frag)).AddWeightedObservation( score_avg.GetAverage(), 1.0 / double( n_bonds_frag));
            }
          }
        }
      }
      for( size_t dhb_index( 0); dhb_index < m_NumberDihedralBonds; ++dhb_index)
      {
        if( dihedral_scores( dhb_index).GetWeight())
        {
          if( final_dihedral_scores.GetSize())
          {
            final_dihedral_scores( dhb_index) = dihedral_scores( dhb_index).GetAverage();
          }
        }
      }
      return final_dihedral_scores;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentProbabilityScore::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FragmentProbabilityScore::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns the returns the bin signature for molecule derived from dihedral angles of fragments
    //! @param MOLECULE molecule of interest
    //! @param FRAGMENT_ISOMORPHISM fragment that will be used to calculate the bin signature
    //! @oaran BIN_SIZE bin size being used for determining dihedral keys
    //! @return all possible dihedral bin signatures for part of molecule corresponding to the fragment being provided
    std::pair< linal::Vector< int>, linal::Vector< double> > FragmentProbabilityScore::GetMolecularBinFromFragment //all of the dihedral bins that make a conformer, hence signature
    (
      const FragmentComplete &MOLECULE,
      const storage::Vector< storage::VectorND< 4, size_t> > &FRAGMENT_ISOMORPHISM,
      typename FragmentProbabilityScore::t_Map &BIN_HASH_TO_BIN
    ) const
    {
      std::pair< linal::Vector< int>, linal::Vector< double> > conf_bin_angles
      (
        linal::Vector< int>( FRAGMENT_ISOMORPHISM.GetSize()),
        linal::Vector< double>( FRAGMENT_ISOMORPHISM.GetSize())
      );
      for( size_t dihedral_index( 0), sz( FRAGMENT_ISOMORPHISM.GetSize()); dihedral_index < sz; ++dihedral_index)
      {
        auto bin_angle_pair( GetDihedralBin( MOLECULE, FRAGMENT_ISOMORPHISM( dihedral_index), BIN_HASH_TO_BIN));
        conf_bin_angles.first( dihedral_index) = bin_angle_pair.first;
        conf_bin_angles.second( dihedral_index) = bin_angle_pair.second;
      }
      return conf_bin_angles;
    }

    //! @brief returns the dihedral bin in which dihedral bond belongs
    //! @param MOLECULE molecule of interest
    //! @param DIHEDRAL_BOND atom indices that make up the dihedral bond
    //! @oaran BIN_SIZE bin size being used for determining dihedral keys
    //! @return dihedral bin in which dihedral bond belongs
    std::pair< int, double> FragmentProbabilityScore::GetDihedralBin // dihedral bin
    (
      const FragmentComplete &MOLECULE,
      const storage::VectorND< 4, size_t> &DIHEDRAL_BOND,
      typename FragmentProbabilityScore::t_Map &BIN_HASH_TO_BIN
    ) const
    {
      uint64_t hash;
      if( sizeof( uint64_t) == sizeof( size_t))
      {
        // 64 bit machine; hash with bit shifts. This works because molecules are anyway limited to 1000 atoms, and shifting
        // by 10 gives us 1024 space per atom index
        hash = ( DIHEDRAL_BOND( 0)) + ( DIHEDRAL_BOND( 1) << 10) + ( DIHEDRAL_BOND( 2) << 20) + ( DIHEDRAL_BOND( 3) << 30);
      }
      else
      {
        hash = uint64_t( DIHEDRAL_BOND( 0)) + ( uint64_t( DIHEDRAL_BOND( 1)) << 10)
               + ( uint64_t( DIHEDRAL_BOND( 2)) << 20) + ( uint64_t( DIHEDRAL_BOND( 3)) << 30);
      }
      auto itr_hashed( BIN_HASH_TO_BIN.Insert( std::make_pair( hash, std::make_pair( 0, 30.0))));
      if( !itr_hashed.second)
      {
        // dihedral was already calculated
        return itr_hashed.first->second;
      }
      // get the existing angle of the dihedral bond, unless the bond has to be linear because either atom type is
      // linear and the bond is not in a ring
      const AtomComplete &second_atom( MOLECULE.GetAtomVector()( DIHEDRAL_BOND.Second()));
      const AtomComplete &third_atom( MOLECULE.GetAtomVector()( DIHEDRAL_BOND.Third()));
      double existing_angle
      (
        ( second_atom.GetAtomType()->GetFormsOnlyLinearBonds() || third_atom.GetAtomType()->GetFormsOnlyLinearBonds())
        && !second_atom.GetBondTypeTo( third_atom)->IsBondInRing()
        ? math::g_Pi
        : linal::Dihedral
          (
            MOLECULE.GetAtomVector()( DIHEDRAL_BOND.First()).GetPosition(),
            second_atom.GetPosition(),
            third_atom.GetPosition(),
            MOLECULE.GetAtomVector()( DIHEDRAL_BOND.Fourth()).GetPosition()
          )
      );
      existing_angle = m_DihedralBinComparer.DetermineWrappedAroundAngle( existing_angle, true);

      // determine the bin from the angle measure
      return
        itr_hashed.first->second
          = std::make_pair( m_DihedralBinComparer.DetermineDihedralKey( existing_angle, true, true), existing_angle);
    }

  } // namespace chemistry
} // namespace bcl
