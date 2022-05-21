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
#include "assemble/bcl_assemble_voxel_grid_aa.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "biol/bcl_biol_aa_base.h"
#include "graph/bcl_graph_const_graph.h"
#include "graph/bcl_graph_undirected_edge.h"
#include "math/bcl_math_running_average_sd.h"
#include "score/bcl_score_aa_pair_clash.h"
#include "score/bcl_score_aa_pair_hi_res_clash.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! @brief calculates AA pair clash score for a vector of AAs via Slicelist usage
    //! @note NOTE: This might have to be moved into a scoring class derived from protein model
    //! @param AAS vector of pointers to the AAs
    //! @return score of the AA pair clash
    double VoxelGridAA::GetAAClashScore( const util::SiPtrVector< const biol::AABase> &AAS)
    {
      // amino acid clash score
      const score::AAPairClash &clash_score( score::AAPairClash::GetInstance());

      double score( 0.0);

      for( size_t i( 0); i < AAS.GetSize(); ++i)
      {
        storage::Vector< storage::Pair< util::SiPtr< const biol::AABase>, double> >
          res( util::VoxelGrid< biol::AABase>::GetNeighbors( *AAS( i), clash_score.GetDistanceCutoff()));
        for( size_t j( 0), r_size( res.GetSize()); j < r_size; ++j)
        {
          const storage::Pair< util::SiPtr< const biol::AABase>, double> &pair( res( j));
          const double closest_distance
          (
            score::AAPairClash::GetInstance().GetClosestDistance( AAS( i)->GetType(), pair.First()->GetType())
          );
          score += score::AAPairClash::GetInstance().CalculateRepulsiveTerm( pair.Second(), closest_distance);
        }
      }

      return score;
    }

    //! @brief calculates AA pair clash score for a vector of AAs via Slicelist usage
    //! @note NOTE: This might have to be moved into a scoring class derived from protein model
    //! @param AAS vector of pointers to the AAs
    //! @return score of the AA pair clash
    double VoxelGridAA::GetAAClashScore()
    {
      // amino acid clash score
      const score::AAPairClash &clash_score( score::AAPairClash::GetInstance());

      double score( 0.0);

      storage::Vector< storage::Triplet< util::SiPtr< const biol::AABase>, util::SiPtr< const biol::AABase>, double> >
        res( util::VoxelGrid< biol::AABase>::GetNeighbors( clash_score.GetDistanceCutoff()));
      for( size_t j( 0), r_size( res.GetSize()); j < r_size; ++j)
      {
        const storage::Triplet< util::SiPtr< const biol::AABase>, util::SiPtr< const biol::AABase>, double>
          &triplet( res( j));
        const double closest_distance
        (
          clash_score.GetClosestDistance( triplet.First()->GetType(), triplet.Second()->GetType())
        );
        score += clash_score.CalculateRepulsiveTerm( triplet.Third(), closest_distance);
      }

      return score;
    }

    //! @brief Get a matrix with counts of AAs interacting at particular distance between each SSE
    linal::Matrix< float> VoxelGridAA::GetSSEInteractionMatrix
    (
      const util::SiPtrVector< const SSE> &SSES,
      const util::SiPtrVector< const biol::AABase> &AAS,
      const size_t &SEQ_EXCLUSION,
      const double &RESOLUTION,
      const bool   &CONSIDER_LOOPS,
      const double &MIN_CONTACT_P,
      const bool   &CONSIDER_POINT_CONTACTS,
      const bool   &MUTUALLY_CLOSEST_DIFF_SSE_CONNECTIONS_ONLY
    )
    {
      bool prev_prefer_ca( m_PreferCA);
      m_PreferCA = true;
      this->SetObjects( AAS);
      storage::Vector< int> first_ids( size_t( 128), std::numeric_limits< int>::max());
      storage::Vector< int> last_ids( size_t( 128), -std::numeric_limits< int>::max());
      for( auto itr( AAS.Begin()), itr_end( AAS.End()); itr != itr_end; ++itr)
      {
        size_t chain_id( ( *itr)->GetChainID());
        if( first_ids( chain_id) > ( *itr)->GetSeqID())
        {
          first_ids( chain_id) = ( *itr)->GetSeqID();
        }
        if( last_ids( chain_id) < ( *itr)->GetSeqID())
        {
          last_ids( chain_id) = ( *itr)->GetSeqID();
        }
      }

      // hash sses for all amino acid ids
      storage::Vector< storage::Vector< size_t> > sse_id( size_t( 128));
      storage::Vector< storage::Vector< size_t> > have_used( size_t( 128));
      for( size_t chain_id( 0), mx_chain_id( first_ids.GetSize()); chain_id < mx_chain_id; ++chain_id)
      {
        if( last_ids( chain_id) >= first_ids( chain_id))
        {
          sse_id( chain_id).Resize( last_ids( chain_id) - first_ids( chain_id) + 1, util::GetUndefined< size_t>());
          if( MUTUALLY_CLOSEST_DIFF_SSE_CONNECTIONS_ONLY)
          {
            have_used( chain_id).Resize( last_ids( chain_id) - first_ids( chain_id) + 1, size_t( 0));
          }
        }
      }

      size_t current_sse_id( 0);
      for
      (
        auto itr_sses( SSES.Begin()), itr_sses_end( SSES.End());
        itr_sses != itr_sses_end;
        ++itr_sses, ++current_sse_id
      )
      {
        // skip loops unless they're desired
        if( !CONSIDER_LOOPS && !( *itr_sses)->GetType()->IsStructured())
        {
          --current_sse_id;
          continue;
        }
        size_t chain_id( ( *itr_sses)->GetChainID());
        const int first_id( first_ids( chain_id));
        storage::Vector< size_t> &chain_vec( sse_id( chain_id));

        int seq_id_start( ( *itr_sses)->GetFirstAA()->GetSeqID() - first_id);
        int seq_id_end( ( *itr_sses)->GetLastAA()->GetSeqID() - first_id);
        for( int pos( seq_id_start); pos <= seq_id_end; ++pos)
        {
          chain_vec( pos) = current_sse_id;
        }
      }

      const size_t n_sses( current_sse_id);
      linal::Matrix< float> interactions( n_sses, n_sses, float( 0));
      score::AAPairHiResClash clash;
      typedef storage::Vector< storage::Triplet< util::SiPtr< const biol::AABase>, util::SiPtr< const biol::AABase>, double> >
              t_Vector;
      t_Vector neighbors( GetNeighbors( RESOLUTION));
      storage::Vector< storage::Vector< math::RunningMinMax< int> > > aas_involved
      (
        n_sses,
        storage::Vector< math::RunningMinMax< int> >( n_sses)
      );
      // if only looking for mutually closest
      if( MUTUALLY_CLOSEST_DIFF_SSE_CONNECTIONS_ONLY)
      {
        // first, remove intra-sse connections and other connections that are simply too close
        size_t n_used( 0);
        for
        (
          auto itr_valid( neighbors.Begin()), itr_placement( neighbors.Begin()), itr_end( neighbors.End());
          itr_valid != itr_end;
          ++itr_valid
        )
        {
          const auto &triplet( *itr_valid);
          const size_t chain_id_a( triplet.First()->GetChainID()), chain_id_b( triplet.Second()->GetChainID());
          const size_t sse1( sse_id( chain_id_a)( triplet.First()->GetSeqID() - first_ids( chain_id_a)));
          const size_t sse2( sse_id( chain_id_b)( triplet.Second()->GetSeqID() - first_ids( chain_id_b)));
          if( sse1 == sse2)
          {
            continue;
          }
          if( biol::SequenceSeparation( *triplet.First(), *triplet.Second()) < SEQ_EXCLUSION)
          {
            // ignore adjacent aas
            continue;
          }
          if( !CONSIDER_LOOPS && !util::IsDefined( std::max( sse1, sse2)))
          {
            continue;
          }
          if( itr_placement != itr_valid)
          {
            *itr_placement = *itr_valid;
          }
          ++itr_placement;
          ++n_used;
        }
        neighbors.Resize( n_used);

        // Now sort connections by distance
        neighbors.Sort( storage::LessThanThird());

        // go through connections between sses in order of ascending distance
        for
        (
          auto itr_valid( neighbors.Begin()), itr_end( neighbors.End());
          itr_valid != itr_end;
          ++itr_valid
        )
        {
          const auto &triplet( *itr_valid);
          const size_t chain_id_a( triplet.First()->GetChainID()), chain_id_b( triplet.Second()->GetChainID());
//            BCL_MessageStd
//            (
//              triplet.First()->GetIdentification()
//              + " " + triplet.Second()->GetIdentification()
//              + " " + util::Format()( triplet.Third())
//            );
          size_t &have_used_a( have_used( chain_id_a)( triplet.First()->GetSeqID() - first_ids( chain_id_a)));
          size_t &have_used_b( have_used( chain_id_b)( triplet.Second()->GetSeqID() - first_ids( chain_id_b)));
          if( have_used_a && have_used_b)
          {
            continue;
          }
          have_used_a = have_used_b = 1;
          const size_t sse1( sse_id( chain_id_a)( triplet.First()->GetSeqID() - first_ids( chain_id_a)));
          const size_t sse2( sse_id( chain_id_b)( triplet.Second()->GetSeqID() - first_ids( chain_id_b)));
          if( !CONSIDER_LOOPS && ( !util::IsDefined( sse1) || !util::IsDefined( sse2)))
          {
            continue;
          }
          aas_involved( sse1)( sse2) += ( triplet.First()->GetSeqID());
          aas_involved( sse2)( sse1) += ( triplet.Second()->GetSeqID());
        }
      }
      else
      {
        for( auto itr_valid( neighbors.Begin()), itr_end( neighbors.End()); itr_valid != itr_end; ++itr_valid)
        {
          const auto &triplet( *itr_valid);
          const size_t chain_id_a( triplet.First()->GetChainID()), chain_id_b( triplet.Second()->GetChainID());
          const size_t sse1( sse_id( chain_id_a)( triplet.First()->GetSeqID() - first_ids( chain_id_a)));
          const size_t sse2( sse_id( chain_id_b)( triplet.Second()->GetSeqID() - first_ids( chain_id_b)));
          if( !util::IsDefined( sse1) || !util::IsDefined( sse2))
          {
            continue;
          }
          if( sse1 == sse2)
          {
            continue;
          }
          if( biol::SequenceSeparation( *triplet.First(), *triplet.Second()) < SEQ_EXCLUSION)
          {
            // ignore adjacent aas
            continue;
          }
          const double p( clash.GetContactProbability( *triplet.First(), *triplet.Second(), triplet.Third()));
          //          BCL_MessageStd
          //          (
          //            triplet.First()->GetIdentification()
          //            + " " + triplet.Second()->GetIdentification()
          //            + " " + util::Format()( p)
          //            + " " + util::Format()( triplet.Third())
          //          );
          if( p > MIN_CONTACT_P)
          {
            aas_involved( sse1)( sse2) += ( triplet.First()->GetSeqID());
            aas_involved( sse2)( sse1) += ( triplet.Second()->GetSeqID());
          }
        }
      }
      for( size_t i( 0); i < n_sses; ++i)
      {
        for( size_t j( i + 1); j < n_sses; ++j)
        {
          const auto &set_ij( aas_involved( i)( j));
          if( set_ij.GetMin() > set_ij.GetMax())
          {
            continue;
          }
          const size_t ij_range( set_ij.GetRange() + 1);
          const auto &set_ji( aas_involved( j)( i));
          const size_t ji_range( set_ji.GetRange() + 1);
          const size_t min_sz( std::min( ij_range, ji_range));
          const size_t interaction_size
          (
            MUTUALLY_CLOSEST_DIFF_SSE_CONNECTIONS_ONLY
            ? std::max( ( ij_range + ji_range) / 2, min_sz)
            : min_sz
          );
          if( CONSIDER_POINT_CONTACTS || min_sz > size_t( 2))
          {
            interactions( i, j) = interactions( j, i) = interaction_size;
          }
        }
      }
      m_PreferCA = prev_prefer_ca;

      return interactions;
    }

    //! @brief Get the SSE connections
    storage::Vector< storage::Triplet< util::SiPtr< const biol::AABase>, util::SiPtr< const biol::AABase>, double> >
    VoxelGridAA::GetSSEConnections
    (
      const util::SiPtrVector< const SSE> &SSES,
      const util::SiPtrVector< const biol::AABase> &AAS,
      const size_t &SEQ_EXCLUSION,
      const double &RESOLUTION,
      const bool   &CONSIDER_LOOPS,
      const bool   &MUTUALLY_CLOSEST_DIFF_SSE_CONNECTIONS_ONLY,
      const bool   &CONSIDER_INTRA_LOOP_CLASHES
    )
    {
      this->SetObjects( AAS);
      storage::Vector< int> first_ids( size_t( 128), std::numeric_limits< int>::max());
      storage::Vector< int> last_ids( size_t( 128), -std::numeric_limits< int>::max());
      for( auto itr( AAS.Begin()), itr_end( AAS.End()); itr != itr_end; ++itr)
      {
        size_t chain_id( ( *itr)->GetChainID());
        if( first_ids( chain_id) > ( *itr)->GetSeqID())
        {
          first_ids( chain_id) = ( *itr)->GetSeqID();
        }
        if( last_ids( chain_id) < ( *itr)->GetSeqID())
        {
          last_ids( chain_id) = ( *itr)->GetSeqID();
        }
      }

      // hash sses for all amino acid ids
      storage::Vector< storage::Vector< size_t> > sse_id( size_t( 128));
      storage::Vector< storage::Vector< double> > min_distance( size_t( 128));
      for( size_t chain_id( 0), mx_chain_id( first_ids.GetSize()); chain_id < mx_chain_id; ++chain_id)
      {
        if( last_ids( chain_id) >= first_ids( chain_id))
        {
          sse_id( chain_id).Resize( last_ids( chain_id) - first_ids( chain_id) + 1, util::GetUndefined< size_t>());
          if( MUTUALLY_CLOSEST_DIFF_SSE_CONNECTIONS_ONLY)
          {
            min_distance( chain_id).Resize( last_ids( chain_id) - first_ids( chain_id) + 1, math::GetHighestBoundedValue< double>());
          }
        }
      }

      size_t current_sse_id( 0);
      for
      (
        auto itr_sses( SSES.Begin()), itr_sses_end( SSES.End());
        itr_sses != itr_sses_end;
        ++itr_sses, ++current_sse_id
      )
      {
        // skip loops (leave them as undefined, since they are not part of a real sse)
        if( !( *itr_sses)->GetType()->IsStructured() && ( !CONSIDER_LOOPS || CONSIDER_INTRA_LOOP_CLASHES))
        {
          continue;
        }
        size_t chain_id( ( *itr_sses)->GetChainID());
        const int first_id( first_ids( chain_id));
        storage::Vector< size_t> &chain_vec( sse_id( chain_id));

        int seq_id_start( ( *itr_sses)->GetFirstAA()->GetSeqID() - first_id);
        int seq_id_end( ( *itr_sses)->GetLastAA()->GetSeqID() - first_id);
        for( int pos( seq_id_start); pos <= seq_id_end; ++pos)
        {
          chain_vec( pos) = current_sse_id;
        }
      }

      typedef storage::Vector< storage::Triplet< util::SiPtr< const biol::AABase>, util::SiPtr< const biol::AABase>, double> >
              t_Vector;
      t_Vector neighbors( GetNeighbors( RESOLUTION));
      typename t_Vector::iterator itr_valid( neighbors.Begin()),
                                  itr_placement( neighbors.Begin()),
                                  itr_end( neighbors.End());
      for( ; itr_valid != itr_end; ++itr_valid)
      {
        const auto &triplet( *itr_valid);
        const size_t chain_id_a( triplet.First()->GetChainID()), chain_id_b( triplet.Second()->GetChainID());
        const size_t sse1( sse_id( chain_id_a)( triplet.First()->GetSeqID() - first_ids( chain_id_a)));
        const size_t sse2( sse_id( chain_id_b)( triplet.Second()->GetSeqID() - first_ids( chain_id_b)));
        if( !CONSIDER_LOOPS && ( !util::IsDefined( sse1) || !util::IsDefined( sse2)))
        {
          continue;
        }
        if( sse1 == sse2 && ( !CONSIDER_LOOPS || util::IsDefined( sse1)))
        {
          continue;
        }
        if( biol::SequenceSeparation( *triplet.First(), *triplet.Second()) < SEQ_EXCLUSION)
        {
          // ignore adjacent aas
          continue;
        }
        if( itr_placement != itr_valid)
        {
          *itr_placement = *itr_valid;
        }
        ++itr_placement;
      }
      neighbors.Resize( std::distance( neighbors.Begin(), itr_placement));
      if( MUTUALLY_CLOSEST_DIFF_SSE_CONNECTIONS_ONLY)
      {
        for( auto itr( neighbors.Begin()), itr_end( neighbors.End()); itr != itr_end; ++itr)
        {
          const auto &triplet( *itr_valid);
          const size_t chain_id_a( triplet.First()->GetChainID()), chain_id_b( triplet.Second()->GetChainID());
          double &min_distance_a( min_distance( chain_id_a)( triplet.First()->GetSeqID() - first_ids( chain_id_a)));
          double &min_distance_b( min_distance( chain_id_b)( triplet.Second()->GetSeqID() - first_ids( chain_id_b)));
          min_distance_a = std::min( min_distance_a, triplet.Third());
          min_distance_b = std::min( min_distance_b, triplet.Third());
        }
        itr_valid = neighbors.Begin();
        itr_placement = neighbors.Begin();
        itr_end = neighbors.End();
        for( ; itr_valid != itr_end; ++itr_valid)
        {
          const auto &triplet( *itr_valid);
          const size_t chain_id_a( triplet.First()->GetChainID()), chain_id_b( triplet.Second()->GetChainID());
          const double &min_distance_a( min_distance( chain_id_a)( triplet.First()->GetSeqID() - first_ids( chain_id_a)));
          const double &min_distance_b( min_distance( chain_id_b)( triplet.Second()->GetSeqID() - first_ids( chain_id_b)));
          if( triplet.Third() <= min_distance_a && triplet.Third() <= min_distance_b)
          {
            if( itr_placement != itr_valid)
            {
              *itr_placement = *itr_valid;
            }
            ++itr_placement;
          }
        }
        neighbors.Resize( std::distance( neighbors.Begin(), itr_placement));
      }
      return neighbors;
    }

    //! @brief Get the SSE connections
    size_t VoxelGridAA::GetMinSSEMovesToRemoveClashes
    (
      const util::SiPtrVector< const SSE> &SSES,
      const util::SiPtrVector< const biol::AABase> &AAS,
      const bool   &CONSIDER_LOOPS
    )
    {
      this->SetObjects( AAS);
      storage::Vector< int> first_ids( size_t( 128), std::numeric_limits< int>::max());
      storage::Vector< int> last_ids( size_t( 128), -std::numeric_limits< int>::max());
      for( auto itr( AAS.Begin()), itr_end( AAS.End()); itr != itr_end; ++itr)
      {
        size_t chain_id( ( *itr)->GetChainID());
        if( first_ids( chain_id) > ( *itr)->GetSeqID())
        {
          first_ids( chain_id) = ( *itr)->GetSeqID();
        }
        if( last_ids( chain_id) < ( *itr)->GetSeqID())
        {
          last_ids( chain_id) = ( *itr)->GetSeqID();
        }
      }

      // hash sses for all amino acid ids
      storage::Vector< storage::Vector< size_t> > sse_id( size_t( 128));
      for( size_t chain_id( 0), mx_chain_id( first_ids.GetSize()); chain_id < mx_chain_id; ++chain_id)
      {
        if( last_ids( chain_id) >= first_ids( chain_id))
        {
          sse_id( chain_id).Resize( last_ids( chain_id) - first_ids( chain_id) + 1, util::GetUndefined< size_t>());
        }
      }

      size_t current_sse_id( 0);
      for
      (
        auto itr_sses( SSES.Begin()), itr_sses_end( SSES.End());
        itr_sses != itr_sses_end;
        ++itr_sses, ++current_sse_id
      )
      {
        // skip loops unless they're desired
        if( !CONSIDER_LOOPS && !( *itr_sses)->GetType()->IsStructured())
        {
          continue;
        }
        size_t chain_id( ( *itr_sses)->GetChainID());
        const int first_id( first_ids( chain_id));
        storage::Vector< size_t> &chain_vec( sse_id( chain_id));

        int seq_id_start( ( *itr_sses)->GetFirstAA()->GetSeqID() - first_id);
        int seq_id_end( ( *itr_sses)->GetLastAA()->GetSeqID() - first_id);
        for( int pos( seq_id_start); pos <= seq_id_end; ++pos)
        {
          chain_vec( pos) = current_sse_id;
        }
      }

      auto neighbors( GetNeighbors( 4.5));
      storage::Set< graph::UndirectedEdge< size_t> > edges;
      score::AAPairHiResClash clash;
      for( auto itr_valid( neighbors.Begin()), itr_end( neighbors.End()); itr_valid != itr_end; ++itr_valid)
      {
        const auto &triplet( *itr_valid);
        const size_t chain_id_a( triplet.First()->GetChainID()), chain_id_b( triplet.Second()->GetChainID());
        const size_t sse1( sse_id( chain_id_a)( triplet.First()->GetSeqID() - first_ids( chain_id_a)));
        const size_t sse2( sse_id( chain_id_b)( triplet.Second()->GetSeqID() - first_ids( chain_id_b)));
        if( !util::IsDefined( sse1) || !util::IsDefined( sse2))
        {
          continue;
        }
        if( sse1 == sse2)
        {
          continue;
        }
        if( biol::SequenceSeparation( *triplet.First(), *triplet.Second()) < 1)
        {
          // ignore adjacent aas
          continue;
        }
        if( clash( *triplet.First(), *triplet.Second(), triplet.Third()) > 0.05)
        {
          edges.InsertElement( graph::UndirectedEdge< size_t>( sse1, sse2, size_t( 1)));
        }
      }
      graph::ConstGraph< size_t, size_t> clash_graph
      (
        storage::CreateIndexVector( SSES.GetSize()),
        storage::Vector< graph::UndirectedEdge< size_t> >( edges.Begin(), edges.End()),
        0
      );
      if( clash_graph.NumEdges() < 2)
      {
        return clash_graph.NumEdges();
      }
      size_t vertex_cover_size( 0);
      while( clash_graph.NumEdges())
      {
        bool had_tree_nodes( false);
        size_t max_node_size( 0), max_node_id( 0);
        // run through the graph; prune any leaf nodes from any sub-trees in the graph
        for( size_t i( 0), n_sses( SSES.GetSize()); i < n_sses; ++i)
        {
          const size_t n_neigh( clash_graph.GetNeighborIndices( i).GetSize());
          if( n_neigh == size_t( 1))
          {
            const size_t neighbor( clash_graph.GetNeighborIndices( i).FirstElement());
            ++vertex_cover_size;
            clash_graph.RemoveAllEdges( neighbor);
            had_tree_nodes = true;
          }
          else if( n_neigh > max_node_size)
          {
            max_node_size = n_neigh;
            max_node_id = i;
          }
        }
        if( !had_tree_nodes)
        {
          // no tree nodes. Approximation : select the node with the highest number of edges and add it to the cover
          ++vertex_cover_size;
          clash_graph.RemoveAllEdges( max_node_id);
        }
      }
      return vertex_cover_size;
    }

    //! @brief Get the SSE connections
    storage::Vector< storage::Pair< util::SiPtr< const SSE>, linal::Vector3D> >
    VoxelGridAA::GetMinSSEMoveIDsToRemoveClashes
    (
      const util::SiPtrVector< const SSE> &SSES,
      const util::SiPtrVector< const biol::AABase> &AAS,
      const bool   &CONSIDER_LOOPS
    )
    {
      // options that control how the optimal movement from one sse away from another is computed.
      static bool s_use_sse_level_moves( false);
      this->SetObjects( AAS);
      storage::Vector< int> first_ids( size_t( 128), std::numeric_limits< int>::max());
      storage::Vector< int> last_ids( size_t( 128), -std::numeric_limits< int>::max());
      for( auto itr( AAS.Begin()), itr_end( AAS.End()); itr != itr_end; ++itr)
      {
        size_t chain_id( ( *itr)->GetChainID());
        if( first_ids( chain_id) > ( *itr)->GetSeqID())
        {
          first_ids( chain_id) = ( *itr)->GetSeqID();
        }
        if( last_ids( chain_id) < ( *itr)->GetSeqID())
        {
          last_ids( chain_id) = ( *itr)->GetSeqID();
        }
      }

      // hash sses for all amino acid ids
      storage::Vector< storage::Vector< size_t> > sse_id( size_t( 128));
      for( size_t chain_id( 0), mx_chain_id( first_ids.GetSize()); chain_id < mx_chain_id; ++chain_id)
      {
        if( last_ids( chain_id) >= first_ids( chain_id))
        {
          sse_id( chain_id).Resize( last_ids( chain_id) - first_ids( chain_id) + 1, util::GetUndefined< size_t>());
        }
      }

      size_t current_sse_id( 0);
      for
      (
        auto itr_sses( SSES.Begin()), itr_sses_end( SSES.End());
        itr_sses != itr_sses_end;
        ++itr_sses, ++current_sse_id
      )
      {
        // skip loops unless they're desired
        if( !CONSIDER_LOOPS && !( *itr_sses)->GetType()->IsStructured())
        {
          continue;
        }
        size_t chain_id( ( *itr_sses)->GetChainID());
        const int first_id( first_ids( chain_id));
        storage::Vector< size_t> &chain_vec( sse_id( chain_id));

        int seq_id_start( ( *itr_sses)->GetFirstAA()->GetSeqID() - first_id);
        int seq_id_end( ( *itr_sses)->GetLastAA()->GetSeqID() - first_id);
        for( int pos( seq_id_start); pos <= seq_id_end; ++pos)
        {
          chain_vec( pos) = current_sse_id;
        }
      }

      const double allowance( 0.05);

      static score::AAPairHiResClash clash;
      auto neighbors( GetNeighbors( clash.GetDistanceCutoff()));
      storage::Set< graph::UndirectedEdge< double> > edges, interactions;
      storage::Map< storage::Pair< size_t, size_t>, math::RunningAverageSD< linal::Vector3D> > clashed_sse_pairs_to_directions;
      for( auto itr_valid( neighbors.Begin()), itr_end( neighbors.End()); itr_valid != itr_end; ++itr_valid)
      {
        auto &triplet( *itr_valid);
        const size_t chain_id_a( triplet.First()->GetChainID()), chain_id_b( triplet.Second()->GetChainID());
        size_t sse1( sse_id( chain_id_a)( triplet.First()->GetSeqID() - first_ids( chain_id_a)));
        size_t sse2( sse_id( chain_id_b)( triplet.Second()->GetSeqID() - first_ids( chain_id_b)));
        if( !util::IsDefined( sse1) || !util::IsDefined( sse2))
        {
          continue;
        }
        if( sse1 == sse2)
        {
          continue;
        }
        if( biol::SequenceSeparation( *triplet.First(), *triplet.Second()) < 1)
        {
          // ignore adjacent aas
          continue;
        }

        if( double clash_amt = clash( *triplet.First(), *triplet.Second(), triplet.Third()))
        {
          clash_amt += allowance;
          linal::Vector3D translation( triplet.Second()->GetCA().GetCoordinates() - triplet.First()->GetCA().GetCoordinates());
          translation.Normalize();
          translation *= clash_amt;
          clashed_sse_pairs_to_directions[ storage::Pair< size_t, size_t>( sse1, sse2)] += translation;
          clashed_sse_pairs_to_directions[ storage::Pair< size_t, size_t>( sse2, sse1)] += -translation;
          auto itr_clash( edges.InternalData().insert( graph::UndirectedEdge< double>( sse1, sse2, clash_amt)));
          if( !itr_clash.second && itr_clash.first->GetEdgeData() < clash_amt)
          {
            edges.RemoveElement( itr_clash.first);
            edges.Insert( graph::UndirectedEdge< double>( sse1, sse2, clash_amt));
          }
        }
        interactions.InsertElement( graph::UndirectedEdge< double>( sse1, sse2, triplet.Third()));
      }
      graph::ConstGraph< size_t, double> clash_graph
      (
        storage::CreateIndexVector( SSES.GetSize()),
        storage::Vector< graph::UndirectedEdge< double> >( edges.Begin(), edges.End()),
        0
      );
      graph::ConstGraph< size_t, double> interaction_graph
      (
        storage::CreateIndexVector( SSES.GetSize()),
        storage::Vector< graph::UndirectedEdge< double> >( interactions.Begin(), interactions.End()),
        0
      );

      storage::Vector< storage::Pair< util::SiPtr< const SSE>, linal::Vector3D> > sses_to_move;
      if( !clash_graph.NumEdges())
      {
        return sses_to_move;
      }

      // Determine translation vectors
      storage::Map< storage::Pair< size_t, size_t>, linal::Vector3D> clash_removal_translations;
      for
      (
        auto itr_trans( clashed_sse_pairs_to_directions.Begin()), itr_trans_end( clashed_sse_pairs_to_directions.End());
        itr_trans != itr_trans_end;
        ++itr_trans
      )
      {
        linal::Vector3D &trans( clash_removal_translations[ itr_trans->first]);
        if
        (
          !s_use_sse_level_moves
          || itr_trans->second.GetAverage().Norm() * 1.5 > itr_trans->second.GetStandardDeviation().Norm()
        )
        {
          // side chains or just one side of backbone touching backbone of other SSE. In this case, just move
          // the minimal amount in the average direction to resolve the clashes
          trans = itr_trans->second.GetAverage();
          trans.Normalize();

          // multiply by the clashed radius
          trans *=
            std::max
            (
              edges.Find
              (
                graph::UndirectedEdge< double>( itr_trans->first.First(), itr_trans->first.Second(), double( 0.0))
              )->GetEdgeData(),
              0.125
            );
        }
        else
        {
          // severe overlap, e.g. one sse going partly through the core of another. In this case, we treat them as cylinders
          // and move them apart according to their radii
          SSEGeometryPacking packing( *SSES( itr_trans->first.Second()), *SSES( itr_trans->first.First()));
          // highly intersecting SSEs. Have to take more drastic steps to remove the clash
          const linal::Vector3D &direction_a( SSES( itr_trans->first.First())->GetMainAxis().GetDirection()),
                                &direction_b( SSES( itr_trans->first.Second())->GetMainAxis().GetDirection());
          linal::Vector3D norm( linal::CrossProduct( direction_a, direction_b));
          double dist( packing.GetDistance());
          double nrm( norm.Norm());
          const double half_radius
          (
            (
              SSES( itr_trans->first.Second())->GetType()->GetRadialExtent()
              + SSES( itr_trans->first.First())->GetType()->GetRadialExtent()
            ) / 2.0
          );
          if( nrm < 1.0 || packing.GetDistance() > 0.5)
          {
            trans = packing.GetShortestConnection().GetDirection();
          }
          else
          {
            trans = -norm;
          }
          trans.Normalize();
          trans *= -edges.Find
                     (
                       graph::UndirectedEdge< double>( itr_trans->first.First(), itr_trans->first.Second(), double( 0.0))
                     )->GetEdgeData();
        }
      }
      size_t vertex_cover_size( 0);
      while( clash_graph.NumEdges())
      {
        bool had_tree_nodes( false);
        size_t min_node_size( clash_graph.NumEdges() + 1), min_node_id( 0);
        // run through the graph; prune any leaf nodes from any sub-trees in the graph
        for( size_t i( 0), n_sses( SSES.GetSize()); i < n_sses; ++i)
        {
          const size_t n_neigh( clash_graph.GetNeighborIndices( i).GetSize());
          if( !n_neigh)
          {
            continue;
          }
          if( n_neigh == size_t( 1))
          {
            const size_t neighbor( clash_graph.GetNeighborIndices( i).FirstElement());
            const size_t n_interaction_i( interaction_graph.GetNeighborIndices( i).GetSize());
            const size_t n_interaction_neigh( interaction_graph.GetNeighborIndices( neighbor).GetSize());
            if( n_interaction_i <= n_interaction_neigh || clash_graph.GetNeighborIndices( neighbor).GetSize() > size_t( 1))
            {
              sses_to_move.PushBack
              (
                storage::Pair< util::SiPtr< const SSE>, linal::Vector3D>
                (
                  SSES( i),
                  clash_removal_translations[ storage::Pair< size_t, size_t>( i, neighbor)]
                )
              );
              clash_graph.RemoveAllEdges( i);
            }
            else
            {
              sses_to_move.PushBack
              (
                storage::Pair< util::SiPtr< const SSE>, linal::Vector3D>
                (
                  SSES( neighbor),
                  clash_removal_translations[ storage::Pair< size_t, size_t>( neighbor, i)]
                )
              );
              clash_graph.RemoveAllEdges( neighbor);
            }
            ++vertex_cover_size;
            had_tree_nodes = true;
          }
          else if
          (
            n_neigh < min_node_size
            ||
            (
              n_neigh == min_node_size
              && interaction_graph.GetNeighborIndices( min_node_id).GetSize()
                 > interaction_graph.GetNeighborIndices( i).GetSize()
            )
          )
          {
            min_node_size = n_neigh;
            min_node_id = i;
          }
        }
        if( !had_tree_nodes)
        {

          // no tree nodes. Approximation : select the node with the highest number of edges and add it to the cover
          ++vertex_cover_size;
          math::RunningAverageSD< linal::Vector3D> ave_direction;
          double max_distance( 0.0);
          for( size_t j( 0); j < min_node_size; ++j)
          {
            const linal::Vector3D &trans
            (
              clash_removal_translations
              [
                 storage::Pair< size_t, size_t>( min_node_id, clash_graph.GetNeighborIndices( min_node_id)( j))
              ]
            );
            ave_direction += trans;
            max_distance = std::max( clash_graph.GetNeighborData( min_node_id)( j), max_distance);
          }
          BCL_MessageDbg
          (
            "Multiple connections on all nodes, taking node with the lowest degree "
            + util::Format()( ave_direction.GetAverage().Norm())
            + " "
            + util::Format()( ave_direction.GetStandardDeviation().Norm())
          );
          linal::Vector3D translation( ave_direction.GetAverage());
          translation.Normalize();
          translation *= max_distance;
          sses_to_move.PushBack
          (
            storage::Pair< util::SiPtr< const SSE>, linal::Vector3D>
            (
              SSES( min_node_id),
              translation
            )
          );
          clash_graph.RemoveAllEdges( min_node_id);
        }
      }
      return sses_to_move;
    }

  } // namespace assemble
} // namespace bcl
