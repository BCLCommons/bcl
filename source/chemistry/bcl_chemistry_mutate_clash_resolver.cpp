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
#include "chemistry/bcl_chemistry_mutate_clash_resolver.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_angle_assignment.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_conformation_comparison_by_dihedral_bins.h"
#include "chemistry/bcl_chemistry_constitution_graph_converter.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_mutate_chirality.h"
#include "chemistry/bcl_chemistry_ring_fragment_map.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "io/bcl_io_ofstream.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_mutate_result.h"
#include "quality/bcl_quality_rmsd.h"

// external includes - sorted alphabetically

#define BCL_PROFILE_MutateClashResolver
#ifdef BCL_PROFILE_MutateClashResolver
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
    const util::SiPtr< const util::ObjectInterface> MutateClashResolver::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateClashResolver())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new MutateClashResolver
    MutateClashResolver *MutateClashResolver::Clone() const
    {
      return new MutateClashResolver( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateClashResolver::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateClashResolver::GetScheme() const
    {
      static const std::string s_name( "MutateMolecule");
      return s_name;
    }

    //! @brief setup this mutate to handle a new molecule
    void MutateClashResolver::Setup
    (
      const FragmentComplete &CONF,
      const util::ShPtrVector< MutateBondAngles> &MUT_BA,
      const util::ShPtr< AtomClashScore> &CLASH
    )
    {
      // Get indices that we are allowed to sample
      linal::Vector< size_t> sample_by_parts( CONF.GetStoredProperties().GetMDLPropertyAsVector( "SampleByParts"));
      storage::Set< size_t> parts_to_sample;
      for
      (
        linal::Vector< size_t>::const_iterator itr_indices( sample_by_parts.Begin()), itr_indices_end( sample_by_parts.End());
        itr_indices != itr_indices_end;
        ++itr_indices
      )
      {
        parts_to_sample.Insert( *itr_indices - size_t( 0));
      }

      m_ClashScore = CLASH;
      m_BondAngleMutators = MUT_BA;
      m_DihedralMutators.Reset();
      m_Map.Reset();
      // graph. It is critical that the edges be labelled by identity here because the FindMinimalPath function used later
      // considers the edge color to be a weighting or cost of traversing the edge
      m_Graph = ConformationGraphConverter
                (
                  ConformationGraphConverter::e_Identity,
                  ConfigurationalBondTypeData::e_Identity
                )( CONF);
      m_DihedralMutators.Resize( CONF.GetSize());
      m_NumberDihedrals = 0;

      auto all_edges( PriorityDihedralAngles::GetDihedralEdges( CONF));
      for( auto itr_edges( all_edges.Begin()), itr_edges_end( all_edges.End()); itr_edges != itr_edges_end; ++itr_edges)
      {
        if
        (
          !parts_to_sample.IsEmpty()
          && ( !parts_to_sample.Contains( itr_edges->GetIndexLow()) || !parts_to_sample.Contains( itr_edges->GetIndexHigh()))
        )
        {
          continue;
        }
        if( itr_edges->GetEdgeData()->IsBondInRing())
        {
          continue;
        }

        util::ShPtr< MutateDihedralBond> mut
        (
          new MutateDihedralBond
          (
            CONF.GetAtomVector()( itr_edges->GetIndexLow()),
            CONF.GetAtomVector()( itr_edges->GetIndexHigh()),
            m_Graph,
            CONF,
            false,
            true
          )
        );
        m_DihedralMutators( itr_edges->GetIndexLow()).PushBack( mut);
        ++m_NumberDihedrals;
      }
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutateClashResolver::GetSerializer() const
    {
      io::Serializer serial;
      serial.SetClassDescription( "Resolves clashes");
      serial.AddInitializer
      (
        "max cycles",
        "maximum number of times clashes are redected and removed",
        io::Serialization::GetAgent( &m_MaxCycles),
        "100"
      );
      return serial;
    }

    //! @brief cache the resolution path for a particular clash
    const storage::Pair< util::ShPtrVector< MutateDihedralBond>, util::ShPtrVector< MutateBondAngles> > &
    MutateClashResolver::GetResolutionPaths
    (
      const size_t &INDEX_A,
      const size_t &INDEX_B
    ) const
    {
      storage::Pair< util::ShPtrVector< MutateDihedralBond>, util::ShPtrVector< MutateBondAngles> > &ret
      (
        m_Map[ std::make_pair( INDEX_A, INDEX_B)]
      );
      if( ret.First().IsEmpty() && ret.Second().IsEmpty())
      {
        graph::Path shortest_path( graph::Connectivity::FindMinimalPath( m_Graph, INDEX_A, INDEX_B));
        util::ShPtrVector< MutateBondAngles> should_look_for_in_path( m_Graph.GetSize());

        // add the bond angle mutates. Arguably we may want to only use the bond angles mutates that are in rings, because
        // otherwise there's probably a better way to resolve
        for( auto itr_ba( m_BondAngleMutators.Begin()), itr_ba_end( m_BondAngleMutators.End()); itr_ba != itr_ba_end; ++itr_ba)
        {
          if( shortest_path.Contains( ( *itr_ba)->GetFragmentData()->GetCentralAtomIndex()))
          {
            if
            (
              ( *itr_ba)->GetFragmentData()->GetCentralAtomIndex() != INDEX_A
              && ( *itr_ba)->GetFragmentData()->GetCentralAtomIndex() != INDEX_B
            )
            {
              // skip bonds angles that involve non-mobile atoms
              if( m_MobileAtoms.Find( ( *itr_ba)->GetFragmentData()->GetCentralAtomIndex()) < m_MobileAtoms.GetSize())
              {
                continue;
              }
              // skip bond angles that involve rings. We'll add them in the next loop if they're necessary
              else if( !( *itr_ba)->GetFragmentData()->GetNumberRingBonds())
              {
                ret.Second().PushBack( *itr_ba);
              }
              // if the bond angle is in a ring, check if the number of ring bonds is less than the total
              // number of neighbor atoms of the central atom in the bond angle
              else if
              (
                ( *itr_ba)->GetFragmentData()->GetNumberRingBonds()
                < m_Graph.GetNeighborIndices( ( *itr_ba)->GetFragmentData()->GetCentralAtomIndex()).GetSize()
              )
              {
                // note that we should check for this vertex when we go over the path to add the dihedral mutates.
                // in general we want to only add bond angle mutate for atom A if the shortest path from the clashing atoms
                // traverses no more than one ring bond of A
                should_look_for_in_path( ( *itr_ba)->GetFragmentData()->GetCentralAtomIndex()) = *itr_ba;
              }
            }
          }
        }

        if( shortest_path.GetSize() > size_t( 3))
        {
          graph::Path::const_iterator itr( shortest_path.Begin()), itr_end( shortest_path.End());
          ++itr;
          --itr_end;
          graph::Path::const_iterator itr_next( itr);
          int position( 0), center( ( shortest_path.GetSize() - 2) / 2);
          util::ShPtrVector< MutateDihedralBond> first_half, second_half;
          for( ++itr_next; itr_next != itr_end; ++itr, ++itr_next, ++position)
          {
            const size_t low_index( std::min( *itr, *itr_next));
            const size_t high_index( std::max( *itr, *itr_next));
            const auto &dihedral_mutates( m_DihedralMutators( low_index));
            bool found( false);
            for( auto itr_mut( dihedral_mutates.Begin()), itr_mut_end( dihedral_mutates.End()); itr_mut != itr_mut_end; ++itr_mut)
            {
              if( ( *itr_mut)->GetAtomIndexC() == high_index || ( *itr_mut)->GetAtomIndexB() == high_index)
              {
                if( position < center)
                {
                  first_half.PushBack( *itr_mut);
                }
                else
                {
                  second_half.PushBack( *itr_mut);
                }
                found = true;
                break;
              }
            }
            if( found && should_look_for_in_path( low_index).IsDefined())
            {
              ret.Second().PushBack( should_look_for_in_path( low_index));
              should_look_for_in_path( low_index).Reset();
            }
            if( found && should_look_for_in_path( high_index).IsDefined())
            {
              ret.Second().PushBack( should_look_for_in_path( high_index));
              should_look_for_in_path( high_index).Reset();
            }
          }
          while( first_half.GetSize() > second_half.GetSize())
          {
            ret.First().PushBack( first_half.LastElement());
            first_half.PopBack();
          }
          while( second_half.GetSize() > first_half.GetSize())
          {
            ret.First().PushBack( second_half.LastElement());
            second_half.PopBack();
          }
          size_t sh_pos( 0);
          while( !first_half.IsEmpty())
          {
            ret.First().PushBack( second_half( sh_pos++));
            ret.First().PushBack( first_half.LastElement());
            first_half.PopBack();
          }
        }
      }
      return ret;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief operator taking a conformation and returns a mutated conformation
    //! @param MOLECULE conformation of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< FragmentComplete> MutateClashResolver::operator()
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      double clash_score( ( *m_ClashScore)( MOLECULE));
      if( clash_score < 1.0e-8)
      {
        // return empty; mutate is skipped
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
      }
      #ifdef BCL_PROFILE_MutateClashResolver
      static util::Stopwatch s_timer( "Clash resolution (not counting clash scoring)", util::Message::e_Standard, true);
      s_timer.Start();
      #endif
      util::ShPtr< FragmentComplete> current_sp( new FragmentComplete( MOLECULE.GetAtomVector(), MOLECULE.GetName()));
      util::ShPtr< FragmentComplete> best_mol_sp( current_sp->Clone());
      size_t n_mutates( 0);
      size_t max_mutates( m_MaxCycles * ( m_NumberDihedrals + MOLECULE.GetNumberAtoms() / 2));
      for( size_t cycle_n( 0); n_mutates < max_mutates && clash_score; ++cycle_n, ++n_mutates)
      {
        double old_cs( clash_score);
        #ifdef BCL_PROFILE_MutateClashResolver
        s_timer.Stop();
        #endif
        clash_score = ( ( *m_ClashScore)( *current_sp));
        #ifdef BCL_PROFILE_MutateClashResolver
        s_timer.Start();
        #endif
        if( old_cs > clash_score)
        {
          *best_mol_sp = *current_sp;
        }

        if( old_cs < clash_score)
        {
          clash_score = old_cs;
        }
        auto clashed_triplets( m_ClashScore->GetClashingPairs( *current_sp));
        for
        (
          auto itr_clashed( clashed_triplets.Begin()), itr_clashed_end( clashed_triplets.End());
          itr_clashed != itr_clashed_end && n_mutates < max_mutates;
          ++itr_clashed
        )
        {
          const double non_clashing_distance
          (
            m_ClashScore->GetMinimumNonClashingDistance( itr_clashed->First(), itr_clashed->Second())
          );
          itr_clashed->Third() = non_clashing_distance -
            linal::Distance
            (
              current_sp->GetAtomVector()( itr_clashed->First()).GetPosition(),
              current_sp->GetAtomVector()( itr_clashed->Second()).GetPosition()
            );
          if( itr_clashed->Third() < 1.0e-8)
          {
            continue;
          }
          const auto &resolution_path( GetResolutionPaths( itr_clashed->First(), itr_clashed->Second()));
//          storage::Vector< size_t> index_visit( storage::CreateIndexVector( resolution_path.First().GetSize()));
//          index_visit.Shuffle();
          for
          (
            size_t mutate_visit( 0), sz_mutates( resolution_path.First().GetSize());
            mutate_visit < sz_mutates && n_mutates < max_mutates;
            ++mutate_visit
          )
          {
            itr_clashed->Third() = non_clashing_distance -
              linal::Distance
              (
                current_sp->GetAtomVector()( itr_clashed->First()).GetPosition(),
                current_sp->GetAtomVector()( itr_clashed->Second()).GetPosition()
              );
            if( itr_clashed->Third() < 1.0e-8)
            {
              break;
            }
            //BCL_MessageStd( "Mutation # " + util::Format()( mutate_visit) + " preclash: " + util::Format()( itr_clashed->Third()));
            const auto &mutate( resolution_path.First()( mutate_visit));
            auto nr
            (
              mutate->RemoveClash( *current_sp, itr_clashed->First(), itr_clashed->Second(), non_clashing_distance)
            );
            ++n_mutates;
            if( nr.GetArgument().IsDefined())
            {
              current_sp = nr.GetArgument();
              itr_clashed->Third() = non_clashing_distance -
                linal::Distance
                (
                  current_sp->GetAtomVector()( itr_clashed->First()).GetPosition(),
                  current_sp->GetAtomVector()( itr_clashed->Second()).GetPosition()
                );
              if( itr_clashed->Third() < 1.0e-8)
              {
                break;
              }
            }
          }
          if( itr_clashed->Third() >= 1.0e-8)
          {
            storage::Vector< size_t> index_visit_ba
            (
              storage::CreateIndexVector( resolution_path.Second().GetSize())
            );
            index_visit_ba.Shuffle();
            for
            (
              size_t ba_choice( 0), ba_mx_choice( index_visit_ba.GetSize());
              ba_choice < ba_mx_choice;
              ++ba_choice
            )
            {
              const auto &mutate( resolution_path.Second()( index_visit_ba( ba_choice)));
              auto nr( ( *mutate)( *current_sp));
              if( nr.GetArgument().IsDefined())
              {
                const double new_clash
                (
                  non_clashing_distance - linal::Distance
                  (
                    nr.GetArgument()->GetAtomVector()( itr_clashed->First()).GetPosition(),
                    nr.GetArgument()->GetAtomVector()( itr_clashed->Second()).GetPosition()
                  )
                );
                if( new_clash < itr_clashed->Third())
                {
                  current_sp = nr.GetArgument();
                  itr_clashed->Third() = new_clash;
                  if( new_clash < 1.0e-8)
                  {
                    break;
                  }
                }
              }
            }
          }
        }
      }
      #ifdef BCL_PROFILE_MutateClashResolver
      s_timer.Stop();
      #endif
      double old_cs( clash_score);
      clash_score = ( ( *m_ClashScore)( *current_sp));
      BCL_MessageDbg( "Done...Clash: " + util::Format()( std::min( clash_score, old_cs)));
      if( old_cs < clash_score)
      {
        return math::MutateResult< FragmentComplete>( best_mol_sp, *this);
      }

      return math::MutateResult< FragmentComplete>( current_sp, *this);
    }

  ///////////////
  // operators //
  ///////////////

  } // namespace chemistry
} // namespace bcl
