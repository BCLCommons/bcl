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
#include "score/bcl_score_protein_model_sse_linear_loop_proximity.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_si_ptr_list.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelSSELinearLoopProximity::s_Instance
    (
      util::Enumerated< ProteinModel>::AddInstance( new ProteinModelSSELinearLoopProximity())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct parameters
    //! @param NORMALIZE if true, final score will be normalized by the number of sses/sequences in the protein model
    //! @param CONSIDER_DISTANCE_ALONG_SSE onsider the distance along the sse if true; usually useful, see comment
    //!        for m_ConsiderDistanceAlongSSE
    //! @param FOOTPOINT_OFFSET sets m_FootPointOffset
    //! @param MAXIMUM_LINEAR_LOOP_RESIDUES maximum number of residues between adjacent sses for virtual loops; useful
    //!                                     so that incomplete models are not penalized in cases where the linear loop
    //!                                     approximation breaks down
    ProteinModelSSELinearLoopProximity::ProteinModelSSELinearLoopProximity
    (
      const bool NORMALIZE,
      const bool CONSIDER_DISTANCE_ALONG_SSE,
      const double &FOOTPOINT_OFFSET,
      const bool &CONSIDER_VIRTUAL_LOOP_CLASHES,
      const size_t &MAXIMUM_LINEAR_LOOP_RESIDUES
    ) :
      m_Normalize( NORMALIZE),
      m_ConsiderDistanceAlongSSE( CONSIDER_DISTANCE_ALONG_SSE),
      m_FootPointOffset( FOOTPOINT_OFFSET),
      m_ConsiderVirtualLoopClashes( CONSIDER_VIRTUAL_LOOP_CLASHES),
      m_MaximumSequenceSeparation( MAXIMUM_LINEAR_LOOP_RESIDUES)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ProteinModelSSELinearLoopProximity
    ProteinModelSSELinearLoopProximity *ProteinModelSSELinearLoopProximity::Clone() const
    {
      return new ProteinModelSSELinearLoopProximity( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelSSELinearLoopProximity::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ProteinModelSSELinearLoopProximity::GetScheme() const
    {
      static const std::string s_name( "sse_linear_loop_proximity"), s_aa_name( "sse_linear_loop_proximity_aa");
      return m_Normalize ? s_aa_name : s_name;
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &ProteinModelSSELinearLoopProximity::GetAlias() const
    {
      static const std::string s_name( "SSELinearLoopProximity");
      return s_name;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinModelSSELinearLoopProximity::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      ( "score for proximity of imaginary linear loops between all neighborings sses to other sses in the model");
      serializer.AddInitializer
      (
        "normalize",
        "Normalizes the final score to # of aas in each associated sse pair",
        io::Serialization::GetAgent( &m_Normalize),
        "false"
      );
      serializer.AddInitializer
      (
        "consider distance alongSSE",
        "Weights each close contact by (x+y+2*footpoint_offset)/(1+2*footpoint_offset), "
        "where x = closest point on sse norm distance from end (0-0.5), "
        "y = closest point on loop norm distance from end (0-0.5)",
        io::Serialization::GetAgent( &m_ConsiderDistanceAlongSSE),
        "true"
      );
      serializer.AddInitializer
      (
        "foot point offset",
        "See equation for m_ConsiderDistanceAlongSSE",
        io::Serialization::GetAgent( &m_FootPointOffset),
        "0.125"
      );
      serializer.AddInitializer
      (
        "consider virtual loop clashes",
        "whether to consider virtual loop clashes",
        io::Serialization::GetAgent( &m_ConsiderVirtualLoopClashes),
        "true"
      );
      serializer.AddInitializer
      (
        "maximum sequence separation",
        "Maximum sequence distance between sses to consider; useful during folding and with incomplete models, "
        "where linear loop approximation breaks down, potentially in such a way that it prevents non-consecutive insertion of"
        " SSEs.  Set to undefined or similarly high number to consider all loops",
        io::Serialization::GetAgent( &m_MaximumSequenceSeparation),
        "40"
      );
      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that scores the chain
    //! @param CHAIN the chain for which all neighbor scores are calculated
    //! @param MISSING_SSES sse pool of missing sses
    //! @return score
    double ProteinModelSSELinearLoopProximity::operator()
    (
      const assemble::Chain &CHAIN,
      const assemble::SSEPool &MISSING_SSES
    ) const
    {
      std::stringstream dummy_stream;
      return ScoreLoops( CHAIN, MISSING_SSES, dummy_stream, false);
    }

    //! @brief operator that scores the Protein model
    //! @param PROTEIN_MODEL the protein model for which all neighbor scores are calculated
    //! @return score
    double ProteinModelSSELinearLoopProximity::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // score sum
      double score( 0.0);

      // get the pool from the given ProteinModel and make sure it is valid
      const util::ShPtr< assemble::SSEPool> sp_pool
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool)
      );

      // find all sse's not yet in the model. do not assume a linear loop between sses if another SSE could exist
      // between them, because that would invalidate the linear loop assumption
      util::SiPtrList< const assemble::SSE> non_overlapping
      (
        sp_pool.IsDefined()
        ? sp_pool->GetNonOverlappingSSEs( PROTEIN_MODEL)
        : util::SiPtrList< const assemble::SSE>()
      );

      // create an sse-pool with the missing SSEs
      assemble::SSEPool sse_pool( non_overlapping, true, false);

      // iterate through the chains of PROTEIN_MODEL
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        score += operator()( **chain_itr, sse_pool);
      }

      // return score
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ProteinModelSSELinearLoopProximity::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Normalize, ISTREAM);
      io::Serialize::Read( m_ConsiderDistanceAlongSSE, ISTREAM);
      io::Serialize::Read( m_FootPointOffset, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ProteinModelSSELinearLoopProximity::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Normalize, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_ConsiderDistanceAlongSSE, OSTREAM, 0) << '\t';
      io::Serialize::Write( m_FootPointOffset, OSTREAM, 0);

      // return the stream
      return OSTREAM;
    }

    //! @brief write detailed scheme and values to OSTREAM
    //! @param PROTEIN_MODEL ProteinModel to be used to evaluate the function
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &ProteinModelSSELinearLoopProximity::WriteDetailedSchemeAndValues
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      std::ostream &OSTREAM
    ) const
    {
      // get the pool from the given ProteinModel and make sure it is valid
      const util::ShPtr< assemble::SSEPool> sp_pool
      (
        PROTEIN_MODEL.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool)
      );

      // find all sse's not yet in the model. do not assume a linear loop between sses if another SSE could exist
      // between them, because that would invalidate the linear loop assumption
      util::SiPtrList< const assemble::SSE> non_overlapping
      (
        sp_pool.IsDefined()
        ? sp_pool->GetNonOverlappingSSEs( PROTEIN_MODEL)
        : util::SiPtrList< const assemble::SSE>()
      );

      // iterate through the chains of PROTEIN_MODEL
      for
      (
        util::ShPtrVector< assemble::Chain>::const_iterator
          chain_itr( PROTEIN_MODEL.GetChains().Begin()), chain_itr_end( PROTEIN_MODEL.GetChains().End());
        chain_itr != chain_itr_end;
        ++chain_itr
      )
      {
        ScoreLoops( **chain_itr, non_overlapping, OSTREAM, true);
      }
      // end
      return OSTREAM;
    }

    //! @brief helper function called by WriteDetailedSchemeAndValues and operator() so that the code remains in-sync
    //! @param CHAIN the chain of interest
    //! @param MISSING_SSES missing sses from the pool
    //! @param OSTREAM the output stream to write the detailed scheme to for this chain
    //! @param DO_WRITE set to true to actually write to the output stream; otherwise, nothing will be written
    //! @return the final score
    double ProteinModelSSELinearLoopProximity::ScoreLoops
    (
      const assemble::Chain &CHAIN,
      const assemble::SSEPool &MISSING_SSES,
      std::ostream &OSTREAM,
      const bool &DO_WRITE
    ) const
    {
      // sum of all scores
      double score( 0.0);

      // need at least three sses in the chain
      if( CHAIN.GetNumberSSEs() < 3)
      {
        return score;
      }

      // collect all non-coil sse's
      const util::SiPtrVector< const assemble::SSE> structured_sses
      (
        CHAIN.GetSSEs( storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND))
      );

      // need at least three sses in the chain
      if( structured_sses.GetSize() < 3)
      {
        return score;
      }

      // determine the number of residues between each structured sse
      storage::List< size_t> residues_between_sses;
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          itr( structured_sses.Begin()), itr_next( structured_sses.Begin() + 1), itr_end( structured_sses.End());
        itr_next != itr_end;
        ++itr, ++itr_next
      )
      {
        const int end_last_sse_seq_id( ( *itr)->GetLastAA()->GetPdbID());
        const int start_next_sse_seq_id( ( *itr_next)->GetFirstAA()->GetPdbID());
        if( size_t( start_next_sse_seq_id - end_last_sse_seq_id) >= m_MaximumSequenceSeparation)
        {
          residues_between_sses.PushBack( start_next_sse_seq_id - end_last_sse_seq_id);
        }
        else
        {
          bool found_potential_sse_between( false);
          // check for any overlapping sses in the pool
          for
          (
            auto itr_pool( MISSING_SSES.Begin()), itr_pool_end( MISSING_SSES.End());
            itr_pool != itr_pool_end;
            ++itr_pool
          )
          {
            const assemble::SSE &potential_sse( **itr_pool);
            if
            (
              potential_sse.GetFirstAA()->GetSeqID() > end_last_sse_seq_id
              && potential_sse.GetLastAA()->GetSeqID() < start_next_sse_seq_id
            )
            {
              // overlapping sse. inflate the number of residues between these sses to ensure that we don't check
              // the distance between the pair
              found_potential_sse_between = true;
              break;
            }
          }
          residues_between_sses.PushBack
          (
            found_potential_sse_between
            ? m_MaximumSequenceSeparation + 1
            : start_next_sse_seq_id - end_last_sse_seq_id
          );
        }
      }

      // map from sse pointer to pair containing: center of the sse, and vector of lines making up the sse geometry
      // This is considerably faster than always converting the individual fragments to lines every time through the
      // inner loop of this algorithm
      storage::Map
      <
        util::SiPtr< const assemble::SSEGeometryInterface>,
        storage::Pair< linal::Vector3D, storage::Vector< coord::LineSegment3D> >
      > sses_to_center_and_segments;

      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          sse_itr( structured_sses.Begin()), sse_itr_end( structured_sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // insert this sse into the map
        storage::Pair< linal::Vector3D, storage::Vector< coord::LineSegment3D> > &
          this_sse_center_and_segments( sses_to_center_and_segments[ **sse_itr]);

        // set the center
        this_sse_center_and_segments.First() = ( *sse_itr)->GetCenter();

        // insert all the line segments
        storage::Vector< coord::LineSegment3D> &segments( this_sse_center_and_segments.Second());
        util::SiPtrVector< const assemble::SSEGeometryInterface> geometries( ( *sse_itr)->GetSSEGeometries());
        segments.AllocateMemory( geometries.GetSize());
        for
        (
          util::SiPtrVector< const assemble::SSEGeometryInterface>::const_iterator
            itr_geo( geometries.Begin()), itr_geo_end( geometries.End());
          itr_geo != itr_geo_end;
          ++itr_geo
        )
        {
          // get the line segment for this fragment
          segments.PushBack( ( *itr_geo)->GetMainAxis());
        }
      }

      // create virtual linear loops between all adjacent SSE elements
      storage::List< assemble::SSEGeometry> linear_loops;

      // iterate through vector of loop sizes
      storage::List< size_t>::const_iterator itr_loop_size( residues_between_sses.Begin());

      // iterate through consecutive sses of "CHAIN" to create all necessary line segments
      // DO NOT call GetCenter() on any SSEGeometryInterface in sses_to_center_and_segments, since the center is not
      // defined for all of them (e.g. coils)
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          sse_itr_a( structured_sses.Begin()),
          sse_itr_b( ++structured_sses.Begin()),
          sse_itr_end( structured_sses.End());
        sse_itr_b != sse_itr_end;
        ++sse_itr_a, ++sse_itr_b, ++itr_loop_size
      )
      {
        // cast the pointers down to geometry interfaces for easier comparison later on
        util::SiPtr< const assemble::SSEGeometryInterface> sse_a_ptr( **sse_itr_a), sse_b_ptr( **sse_itr_b);

        // test for excessive residues between the two ends of the loop; if so, it is likely that the linear loop
        // approximation will be invalid
        if( *itr_loop_size > m_MaximumSequenceSeparation)
        {
          continue;
        }

        // get the line from the end of the last structured element to the beginning of this one
        const coord::LineSegment3D coil( ( *sse_itr_a)->EndOfZ(), ( *sse_itr_b)->BeginOfZ());

        // walk over the map
        for
        (
          storage::Map
          <
            util::SiPtr< const assemble::SSEGeometryInterface>,
            storage::Pair< linal::Vector3D, storage::Vector< coord::LineSegment3D> >
          >::const_iterator
            itr_map( sses_to_center_and_segments.Begin()),
            itr_map_end( sses_to_center_and_segments.End());
          itr_map != itr_map_end;
          ++itr_map
        )
        {
          const assemble::SSEGeometryInterface &sse_geometry( *itr_map->first);
          // skip sse's that are immediately adjacent to this sse (e.g. the sses bordering this coil)
          if( sse_a_ptr == itr_map->first || sse_b_ptr == itr_map->first)
          {
            continue;
          }

          // for each segment of the SSE
          const linal::Vector3D &sse_center( itr_map->second.First());
          const storage::Vector< coord::LineSegment3D> &sse_segments( itr_map->second.Second());

          // find the highest (most-clashing) score for any segment of this sse
          // One could add up all the scores for all the segments of the sse, but then the score would increase (gets worse)
          // whenever the SSEs are bent, so it is therefore desireable to only take the highest-scoring
          double highest_score( 0), best_interaction_weight( 0);
          double highest_score_distance( std::numeric_limits< double>::infinity());
          const double &radial_extent( sse_geometry.GetType()->GetRadialExtent());
          const double sse_length( sse_geometry.GetLength());

          for
          (
            storage::Vector< coord::LineSegment3D>::const_iterator
              itr_segment( sse_segments.Begin()), itr_segment_end( sse_segments.End());
            itr_segment != itr_segment_end;
            ++itr_segment
          )
          {
            const coord::LineSegment3D &sse_segment( *itr_segment);

            // test for closest distance and whether it was orthogonal
            storage::Pair< coord::LineSegment3D, bool>
              line_orthog( coord::ShortestConnectionBetweenLineSegments3D( coil, sse_segment));

            // if the foot-point of the connection is at the end of either line, skip it, because this indicates that the
            // sse is spatially staggered from the linear coil (as opposed to being in its way or nearly in its way)
            if( !line_orthog.Second())
            {
              continue;
            }

            // compute the distance weight
            // if the loop is < 1/2 the radial extent of the SSE, leave it at 1,
            // otherwise decrease linearly to 0
            // TODO: make this continuously diff; maybe cosine
            const double distance( line_orthog.First().GetLength());
            double distance_weight( 0.0);
            if( distance < 0.5 * radial_extent)
            {
              distance_weight = 1.0;
            }
            else if( distance < radial_extent)
            {
              distance_weight = 1.0 - 2.0 * distance / radial_extent;
            }
            else
            {
              continue;
            }

            double interaction_weight( 1);
            if( m_ConsiderDistanceAlongSSE)
            {
              // compute the interaction weight based on how far along the near-intersection point is along the sse and
              // the loop
              const coord::LineSegment3D &shortest_line_loop_to_sse( line_orthog.First());
              const linal::Vector3D &footpoint_loop( shortest_line_loop_to_sse.GetStartPoint());
              const linal::Vector3D &footpoint_sse( shortest_line_loop_to_sse.GetEndPoint());

              // get the distance along the loop
              const double norm_dist_along_loop
              (
                linal::Distance( footpoint_loop, coil.GetStartPoint()) / coil.GetLength()
              );
              // normalized distances to loop end (0 - 0.5)
              const double norm_dist_to_loop_end( std::min( norm_dist_along_loop, 1.0 - norm_dist_along_loop));

              // get the distance from the center of the sse to the closest connection point
              const double distance_along_sse
              (
                0.5 - std::min( linal::Distance( footpoint_sse, sse_center) / sse_length, 0.5)
              );

              // compute the weighting
              // average the normalized distances; this is because it is bad for either sse to be intersecting anywhere
              // other than at the ends (which are much more commonly near-intersected at in native models)
              interaction_weight = ( distance_along_sse + norm_dist_to_loop_end + 2.0 * m_FootPointOffset)
                                   /
                                   ( 1.0 + 2.0 * m_FootPointOffset);
            }

            const double local_score( interaction_weight * distance_weight);

            // track the min distance for this sse segment
            if( local_score > highest_score)
            {
              highest_score_distance = distance;
              highest_score = local_score;
              best_interaction_weight = interaction_weight;
            }
          }

          if( highest_score)
          {
            double sse_length_av( double( ( *sse_itr_a)->GetSize() + ( *sse_itr_b)->GetSize()) / 20.0);
            score += ( m_Normalize ? sse_length_av : 1.0) * highest_score;
            if( DO_WRITE)
            {
              OSTREAM << "Loop between " << ( *sse_itr_a)->GetIdentification()
                << " - "
                << ( *sse_itr_b)->GetIdentification()
                << " closest approach from "
                << sse_geometry.GetIdentification()
                << " was: " << highest_score_distance
                << " Angstrom; score contribution was: " << highest_score << " unweighted; weighted: " << sse_length_av * highest_score
                << " interaction weight was: " << best_interaction_weight
                << '\n';
            }
          }
        }

        if( m_ConsiderVirtualLoopClashes)
        {
          // add this linear coil to the map
          // Give it type helix so that it has a defined radial extent.  Coils have larger extents than either helices or
          // strands typically, but helix is the better approximation
          linear_loops.PushBack
          (
            assemble::SSEGeometry
            (
              biol::GetSSTypes().HELIX,
              "COIL " + ( *sse_itr_a)->GetLastAA()->GetIdentification()
              + " - " + ( *sse_itr_b)->GetFirstAA()->GetIdentification(),
              coil.GetLength()
            )
          );

          // insert this sse into the map
          storage::Pair< linal::Vector3D, storage::Vector< coord::LineSegment3D> > &
            this_sse_center_and_segments( sses_to_center_and_segments[ linear_loops.LastElement()]);

          // set the center
          this_sse_center_and_segments.First() = 0.5 * ( coil.GetStartPoint() + coil.GetEndPoint());

          // insert all the single line segment
          storage::Vector< coord::LineSegment3D> &segments( this_sse_center_and_segments.Second());
          segments.PushBack( coil);

          // insert the loops size
          residues_between_sses.PushBack( *itr_loop_size);
        }
      }

      // normalize
      if( DO_WRITE)
      {
        OSTREAM << "Final " << GetScheme() << " score: " << score << '\n';
      }

      // return the score sum
      return score;
    }

  } // namespace score
} // namespace bcl
