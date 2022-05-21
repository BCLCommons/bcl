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
#include "score/bcl_score_protein_model_sse_chirality.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "assemble/bcl_assemble_voxel_grid_aa.h"
#include "command/bcl_command_command_state.h"
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_const_graph.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serializer.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_si_ptr_list.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelSSEChirality::s_Instance
    (
      util::Enumerated< ProteinModel>::AddInstance( new ProteinModelSSEChirality())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    //! @param PSEUDOCOUNT pseudocount to use
    ProteinModelSSEChirality::ProteinModelSSEChirality
    (
      const double &PSEUDOCOUNT,
      const bool   &CONTACT_SPECIFIC,
      const bool   &SCALE_BY_AAS
    ) :
      m_Pseudocount( PSEUDOCOUNT),
      m_UseContactSpecificWeights( CONTACT_SPECIFIC),
      m_ScaleByAAs( SCALE_BY_AAS)
    {
      std::stringstream err_str;
      BCL_Assert
      (
        this->ReadInitializerSuccessHook( util::ObjectDataLabel(), err_str),
        "Couldn't read histogram file: " + err_str.str()
      );
    }

    //! @brief Clone function
    //! @return pointer to new ProteinModelSSEChirality
    ProteinModelSSEChirality *ProteinModelSSEChirality::Clone() const
    {
      return new ProteinModelSSEChirality( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelSSEChirality::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ProteinModelSSEChirality::GetScheme() const
    {
      static const std::string
        s_contact_name( "all_sse_chirality"),
        s_contact_sp_name( "all_sse_chirality_contact"),
        s_contact_name_aa( "all_sse_chirality_aa"),
        s_contact_sp_name_aa( "all_sse_chirality_contact_aa");
      return m_ScaleByAAs
             ? (
                 m_UseContactSpecificWeights
                 ? s_contact_sp_name_aa
                 : s_contact_name_aa
               )
             : (
                 m_UseContactSpecificWeights
                 ? s_contact_sp_name
                 : s_contact_name
               );
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &ProteinModelSSEChirality::GetAlias() const
    {
      static const std::string s_combined( "SSEChirality"), s_ann_name( "SSEChiralityEff"), s_con_name( "SSEChiralityConSpec");
      return m_UseContactSpecificWeights ? s_con_name : s_combined;
    }

    //! @brief get the propensities vector
    const linal::Vector< double> &ProteinModelSSEChirality::GetPropensities() const
    {
      return m_Propensities;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinModelSSEChirality::GetSerializer() const
    {
      io::Serializer serializer;

      serializer.SetClassDescription
      (
        "Scores interaction energy of all SSE pairs and chirality of all triplets"
      );
      serializer.AddInitializer
      (
        "pseudocount",
        "pseudocount to add to each count",
        io::Serialization::GetAgent( &m_Pseudocount),
        "20"
      );
      serializer.AddInitializer
      (
        "contact specific",
        "use hi-res weights that further break down which SSEs are in contact",
        io::Serialization::GetAgent( &m_UseContactSpecificWeights),
        "0"
      );
      serializer.AddInitializer
      (
        "scale by aas",
        "whether to scale by the number of aas in the protein",
        io::Serialization::GetAgent( &m_ScaleByAAs),
        "True"
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
    double ProteinModelSSEChirality::operator()
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
    double ProteinModelSSEChirality::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
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

    //! @brief write detailed scheme and values to OSTREAM
    //! @param PROTEIN_MODEL ProteinModel to be used to evaluate the function
    //! @param OSTREAM std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &ProteinModelSSEChirality::WriteDetailedSchemeAndValues
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
      assemble::SSEPool non_overlapping
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
    double ProteinModelSSEChirality::ScoreLoops
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
      if( CHAIN.GetNumberSSEs() < 2)
      {
        return score;
      }

      // min probability to score two residues as being an almost-certain contact
      const double min_p( 0.04);

      // Create a matrix that will hold the x,y, and z coordinates for terminii of the first & second strands
      linal::Matrix3x3< double> xyz_coordinates( 0.0); // make a matrix of size 3 X 3

      // collect all non-coil sse's
      const util::SiPtrVector< const assemble::SSE> structured_sses
      (
        CHAIN.GetSSEs( storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND))
      );
      const storage::Set< biol::AtomType> types_of_interest
      (
        biol::GetAtomTypes().GetFirstSidechainAtomTypes()
      );

      util::SiPtrVector< const biol::Atom> si_atoms;
      math::RunningAverage< double> score_ave;
      assemble::VoxelGridAA interactions_detector( m_HashMaker.GetInteractionDistance());
      const size_t min_atoms_helix( m_HashMaker.GetMinAtomsInContactHelix());
      const size_t min_atoms_strand( m_HashMaker.GetMinAtomsInContactStrand());
      linal::Matrix< float> interactions_matrix
      (
        interactions_detector.GetSSEInteractionMatrix
        (
          structured_sses,
          CHAIN.GetAminoAcids(),
          2,
          m_HashMaker.GetInteractionDistance(),
          false,
          min_p,
          false
        )
      );
      // cache of the packing objects computed so far
      storage::Vector< storage::Vector< assemble::SSEGeometryPacking::OrientationEnum> > orientations
      (
        structured_sses.GetSize(),
        storage::Vector< assemble::SSEGeometryPacking::OrientationEnum>
        (
          structured_sses.GetSize(),
          assemble::SSEGeometryPacking::OrientationEnum( assemble::SSEGeometryPacking::s_NumberOrientations)
        )
      );

      storage::Vector< math::RunningAverage< double> > sse_average_chirality_score( structured_sses.GetSize());
      storage::Vector< math::RunningAverage< double> > sse_average_adj_chirality_score( structured_sses.GetSize());

      // handle SSE-pair interaction-weight and parallel/anti-parallel bias
      size_t a( 0), atoms_in_sses( 0);
      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          itr_a( structured_sses.Begin()),
          itr_end( structured_sses.End());
        itr_a != itr_end;
        ++itr_a, ++a
      )
      {
        size_t b( a + 1);
        atoms_in_sses += ( *itr_a)->GetSize();
        const size_t min_atoms_a( ( *itr_a)->GetType() == biol::GetSSTypes().STRAND ? min_atoms_strand : min_atoms_helix);
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator itr_b( itr_a + 1);
          itr_b != itr_end;
          ++itr_b, ++b
        )
        {
          // for every pair of SSEs
          const assemble::SSEGeometryPacking::OrientationEnum &packing_ab
          (
            m_HashMaker.GetCacheOrientation( orientations, a, b, **itr_a, **itr_b)
          );
          const size_t min_atoms_b( ( *itr_b)->GetType() == biol::GetSSTypes().STRAND ? min_atoms_strand : min_atoms_helix);

          util::SiPtrVector< const assemble::SSE>::const_iterator itr_c( itr_b + 1);
          if( itr_c == itr_end)
          {
            continue;
          }
          size_t c( b + 1);
          const size_t interactions_ab( interactions_matrix( a, b));
          const size_t min_atoms_c( ( *itr_c)->GetType() == biol::GetSSTypes().STRAND ? min_atoms_strand : min_atoms_helix);

          if
          (
            a + 1 == b
            &&
            ( interactions_ab >= min_atoms_a * min_atoms_b ? 1 : 0)
            + ( interactions_matrix( a, c) >= min_atoms_a * min_atoms_c ? 1 : 0)
            + ( interactions_matrix( b, c) >= min_atoms_b * min_atoms_c ? 1 : 0)
            > 1
          )
          {
            const int end_last_sse_seq_id( ( *itr_a)->GetLastAA()->GetPdbID());
            const int start_next_sse_seq_id( ( *itr_c)->GetFirstAA()->GetPdbID());

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
            if( !found_potential_sse_between)
            {
              const assemble::SSEGeometryPacking::OrientationEnum &packing_ac
              (
                m_HashMaker.GetCacheOrientation( orientations, a, c, **itr_a, **itr_c)
              );
              const assemble::SSEGeometryPacking::OrientationEnum &packing_bc
              (
                m_HashMaker.GetCacheOrientation( orientations, b, c, **itr_b, **itr_c)
              );
              const linal::Vector3D root_position( ( *itr_b)->GetCenter());

              xyz_coordinates.GetRow( 0).CopyValues( ( *itr_a)->GetMainAxis().GetStartPoint() - root_position);
              xyz_coordinates.GetRow( 1).CopyValues( ( *itr_a)->GetMainAxis().GetEndPoint() - root_position);
              xyz_coordinates.GetRow( 2).CopyValues( ( *itr_c)->GetCenter() - root_position);

              // Calculate the determinant of a 3 X 3 matrix that has rows sorted in descending order of priority.
              // Opposite orders will have opposite signs. The sign is assigned to a value of 1 and returned as the
              // value of the stereocenter.
              const float determinant( xyz_coordinates.Determinant());
              const double triplet_adj_score
              (
                m_Propensities
                (
                  m_HashMaker.GetPackingTripletNumber
                  (
                    ( *itr_a)->GetType(),
                    ( *itr_b)->GetType(),
                    ( *itr_c)->GetType(),
                    packing_ab,
                    packing_bc,
                    packing_ac,
                    interactions_ab,
                    interactions_matrix( b, c),
                    interactions_matrix( a, c),
                    true,
                    determinant >= 0.0
                  )
                )
              );
              sse_average_adj_chirality_score( a).AddWeightedObservation( triplet_adj_score, interactions_ab + interactions_matrix( a, c));
              sse_average_adj_chirality_score( b).AddWeightedObservation( triplet_adj_score, interactions_ab + interactions_matrix( b, c));
              sse_average_adj_chirality_score( c).AddWeightedObservation( triplet_adj_score, interactions_matrix( b, c) + interactions_matrix( a, c));
              if( DO_WRITE)
              {
                OSTREAM << "Triplet combined_energy adjacent: "
                        << ( *itr_a)->GetIdentification()
                        << " " << ( *itr_b)->GetIdentification()
                        << " " << ( *itr_c)->GetIdentification()
                        << " " << packing_ab.GetString()
                        << " " << packing_bc.GetString()
                        << " " << packing_ac.GetString()
                        << " " << bool( determinant >= 0)
                        << " " << interactions_ab
                        << " " << interactions_matrix( b, c)
                        << " " << interactions_matrix( a, c)
                        << " " << triplet_adj_score << std::endl;
              }
              ++itr_c;
              ++c;
            }
          }
          for( ; itr_c != itr_end; ++itr_c, ++c)
          {
            const size_t min_atoms_c( ( *itr_c)->GetType() == biol::GetSSTypes().STRAND ? min_atoms_strand : min_atoms_helix);

            if
            (
              ( interactions_ab >= min_atoms_a * min_atoms_b ? 1 : 0)
              + ( interactions_matrix( a, c) >= min_atoms_a * min_atoms_c ? 1 : 0)
              + ( interactions_matrix( b, c) >= min_atoms_b * min_atoms_c ? 1 : 0)
              <= 1
            )
            {
              continue;
            }
            const assemble::SSEGeometryPacking::OrientationEnum &packing_ac
            (
              m_HashMaker.GetCacheOrientation( orientations, a, c, **itr_a, **itr_c)
            );
            const assemble::SSEGeometryPacking::OrientationEnum &packing_bc
            (
              m_HashMaker.GetCacheOrientation( orientations, b, c, **itr_b, **itr_c)
            );
            const linal::Vector3D root_position( ( *itr_b)->GetCenter());

            xyz_coordinates.GetRow( 0).CopyValues( ( *itr_a)->GetMainAxis().GetStartPoint() - root_position);
            xyz_coordinates.GetRow( 1).CopyValues( ( *itr_a)->GetMainAxis().GetEndPoint() - root_position);
            xyz_coordinates.GetRow( 2).CopyValues( ( *itr_c)->GetCenter() - root_position);

            // Calculate the determinant of a 3 X 3 matrix that has rows sorted in descending order of priority.
            // Opposite orders will have opposite signs. The sign is assigned to a value of 1 and returned as the
            // value of the stereocenter.
            const float determinant( xyz_coordinates.Determinant());
            const double triplet_nonadj_score
            (
              m_Propensities
              (
                m_HashMaker.GetPackingTripletNumber
                (
                  ( *itr_a)->GetType(),
                  ( *itr_b)->GetType(),
                  ( *itr_c)->GetType(),
                  packing_ab,
                  packing_bc,
                  packing_ac,
                  interactions_ab,
                  interactions_matrix( b, c),
                  interactions_matrix( a, c),
                  false,
                  determinant >= 0.0
                )
              )
            );
            sse_average_chirality_score( a).AddWeightedObservation( triplet_nonadj_score, interactions_ab + interactions_matrix( a, c));
            sse_average_chirality_score( b).AddWeightedObservation( triplet_nonadj_score, interactions_ab + interactions_matrix( b, c));
            sse_average_chirality_score( c).AddWeightedObservation( triplet_nonadj_score, interactions_matrix( b, c) + interactions_matrix( a, c));
            if( DO_WRITE)
            {
              OSTREAM << "Triplet combined_energy non-adjacent: "
                      << ( *itr_a)->GetIdentification()
                      << " " << ( *itr_b)->GetIdentification()
                      << " " << ( *itr_c)->GetIdentification()
                      << " " << packing_ab.GetString()
                      << " " << packing_bc.GetString()
                      << " " << packing_ac.GetString()
                      << " " << bool( determinant >= 0)
                      << " " << interactions_ab
                      << " " << interactions_matrix( b, c)
                      << " " << interactions_matrix( a, c)
                      << " " << triplet_nonadj_score << std::endl;
            }
          }
        }
      }

      for( size_t i( 0), n_sses( structured_sses.GetSize()); i < n_sses; ++i)
      {
        const double ave_adj( sse_average_adj_chirality_score( i).GetAverage());
        const double ave_nonadj( sse_average_chirality_score( i).GetAverage());

        // if there are both adjacent and non-adjacent terms that agree in sign,
        // take whichever is larger in magnitude. This avoids that a highly-favorable adjacent chirality term becomes
        // less favorable when paired with other favorable non-adjacent terms
        if
        (
          sse_average_adj_chirality_score( i).GetWeight()
          && sse_average_chirality_score( i).GetWeight()
          && ( ( ave_adj > 0.0) == ( ave_nonadj > 0.0))
        )
        {
          score_ave.AddWeightedObservation
          (
            ( ave_adj < 0.0 ? std::min( ave_adj, ave_nonadj) : std::max( ave_adj, ave_nonadj)),
            structured_sses( i)->GetSize()
          );
        }
        else
        {
          score_ave.AddWeightedObservation
          (
            ave_adj + ave_nonadj,
            structured_sses( i)->GetSize()
          );
        }
      }
      // return the score sum
      return score_ave.GetAverage() * ( m_ScaleByAAs ? atoms_in_sses : 1.0);
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @return ERROR_STREAM stream with which to write errors
    bool ProteinModelSSEChirality::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      if( !command::CommandState::GetGlobalCommandState().IsInStaticInitialization())
      {
        io::IFStream input;
        io::File::MustOpenIFStream
        (
          input,
          Score::AddHistogramPath
          (
            m_UseContactSpecificWeights
            ? "sse_triplet_chirality_contact.histogram"
            : "sse_triplet_chirality.histogram"
          )
        );
        util::ChopHeader( input);
        bool read_success( m_HashMaker.TryRead( util::ObjectDataLabel( input), ERROR_STREAM));
        if( !read_success)
        {
          return false;
        }
        util::ChopHeader( input);
        m_Propensities = linal::Vector< double>( m_HashMaker.GetNumberHashes(), 0.0);
        storage::Map< std::string, size_t> hash_to_index;
        for( size_t i( 0); i < m_HashMaker.GetNumberHashes(); ++i)
        {
          hash_to_index[ m_HashMaker.GetPackingTripletString( size_t( i))] = i;
        }
        storage::Vector< storage::Vector< std::string> > tokens
        (
          util::SplittedStringLineListFromIStream( input, " \n\t")
        );
        if( !tokens.IsEmpty() && tokens.LastElement().IsEmpty())
        {
          // remove empty lines at the end
          tokens.PopBack();
        }
        if( tokens.GetSize() != m_HashMaker.GetNumberHashes())
        {
          ERROR_STREAM << "Wrong number of hashes in file, should be: "
                       << m_HashMaker.GetNumberHashes() << " but found " << tokens.GetSize() << '\n';
          return false;
        }
        for( auto itr_line( tokens.Begin()), itr_line_end( tokens.End()); itr_line != itr_line_end; ++itr_line)
        {
          const storage::Vector< std::string> &line( *itr_line);
          if( line.GetSize() != 4 && line.GetSize() != 2)
          {
            ERROR_STREAM << "Every line should have two fields (ID and propensity) or four fields (ID, weighted SSEs, raw SSEs, expected SSEs)"
                         << " but one field had: "
                         << line.GetSize()
                         << " : " << util::Join( " ", line) << '\n';
            return false;
          }
          auto itr_hash_to_index( hash_to_index.Find( line( 0)));
          if( itr_hash_to_index == hash_to_index.End())
          {
            ERROR_STREAM << "Hash: " << line( 0) << " is not recognized.\n";
            return false;
          }
          const size_t hash_index( itr_hash_to_index->second);
          if( line.GetSize() == size_t( 2))
          {
            const double counts( util::ConvertStringToNumericalValue< double>( line( 1)));
            m_Propensities( hash_index) = counts;
          }
          else
          {
            const double counts( util::ConvertStringToNumericalValue< double>( line( 2))),
                         expected_counts( util::ConvertStringToNumericalValue< double>( line( 3)));
            m_Propensities( hash_index) = -std::log( ( counts + m_Pseudocount) / ( expected_counts + m_Pseudocount));
          }
        }
        io::File::CloseClearFStream( input);
      }
      return true;
    }

  } // namespace score
} // namespace bcl
