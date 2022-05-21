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
// include header for this class
#include "scorestat/bcl_scorestat_protein_model_sse_triplet_chirality.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_voxel_grid_aa.h"
#include "assemble/bcl_assemble_voxel_grid_atom.h"
#include "biol/bcl_biol_aa_back_bone.h"
#include "chemistry/bcl_chemistry_aa_fragment_complete.h"
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "descriptor/bcl_descriptor_atom_sigma_charge.h"
#include "descriptor/bcl_descriptor_coulombic_force.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_2d.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_histogram_3d.h"
#include "math/bcl_math_running_average.h"
#include "score/bcl_score_aa_pair_hi_res_clash.h"
#include "score/bcl_score_aa_pair_sidechain_interaction.h"
#include "util/bcl_util_string_functions.h"
#include "util/bcl_util_wrapper.h"
namespace bcl
{
  namespace scorestat
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinModelSSETripletChirality::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new ProteinModelSSETripletChirality())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    ProteinModelSSETripletChirality *ProteinModelSSETripletChirality::Clone() const
    {
      return new ProteinModelSSETripletChirality( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelSSETripletChirality::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &ProteinModelSSETripletChirality::GetOutFilePostfix() const
    {
      static const std::string s_name( "sse_triplet_chirality"), s_con_name( "sse_triplet_chirality_contact"), s_rows_name( "sse_triplet_chirality_per_protein");
      return m_OutputPerProtein ? s_rows_name : ( m_ConsiderWhichSSEsInContact ? s_con_name : s_name);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &ProteinModelSSETripletChirality::GetAlias() const
    {
      static const std::string s_name( "SSETripletChirality");
      return s_name;
    }

    //! @brief helper function to cache orientational information
    //! @param CACHE the cache to retrieve/store data in
    //! @param SSE_A_ID index of SSE_A
    //! @param SSE_B_ID index of SSE_B
    //! @param SSE_A, SSE_B the two sses of interest
    //! @return the orientation
    const assemble::SSEGeometryPacking::OrientationEnum &ProteinModelSSETripletChirality::GetCacheOrientation
    (
      storage::Vector< storage::Vector< assemble::SSEGeometryPacking::OrientationEnum> > &CACHE,
      const size_t &SSE_A_ID,
      const size_t &SSE_B_ID,
      const assemble::SSE &SSE_A,
      const assemble::SSE &SSE_B
    ) const
    {
      assemble::SSEGeometryPacking::OrientationEnum &ori( CACHE( SSE_A_ID)( SSE_B_ID));
      if( ori == assemble::SSEGeometryPacking::s_NumberOrientations)
      {
        ori = assemble::SSEGeometryPacking::OrientationEnum
              (
                assemble::SSEGeometryPacking::OrientationFromSSEs( SSE_A, SSE_B)
              );
      }
      return ori;
    }

    //! @brief hash a given packing. The given string is intended to be fast to hash, not necessarily easily readable
    std::string ProteinModelSSETripletChirality::GetPackingTripletHash
    (
      const biol::SSType &TYPE_A,
      const biol::SSType &TYPE_B,
      const biol::SSType &TYPE_C,
      const assemble::SSEGeometryPacking::OrientationEnum &PACKING_AB,
      const assemble::SSEGeometryPacking::OrientationEnum &PACKING_BC,
      const assemble::SSEGeometryPacking::OrientationEnum &PACKING_AC,
      const size_t &COUNTS_PACK_A,
      const size_t &COUNTS_PACK_B,
      const size_t &COUNTS_PACK_C,
      const bool   &ADJACENT,
      const bool   &RHS
    ) const
    {
      const size_t min_atoms_a
      (
        TYPE_A == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
      );
      const size_t min_atoms_b
      (
        TYPE_B == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
      );
      const size_t min_atoms_c
      (
        TYPE_C == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
      );
      return std::string
      (
        std::string( size_t( 1), TYPE_A.GetName()[ 0])
        + std::string( size_t( 1), TYPE_B.GetName()[ 0])
        + std::string( size_t( 1), TYPE_C.GetName()[ 0])
        + std::string( size_t( 1), PACKING_AB.GetString()[ 0])
        + std::string( size_t( 1), PACKING_BC.GetString()[ 0])
        + std::string( size_t( 1), PACKING_AC.GetString()[ 0])
        + (
            m_ConsiderWhichSSEsInContact
            ? util::Format()( COUNTS_PACK_A >= min_atoms_a * min_atoms_b)
              + util::Format()( COUNTS_PACK_B >= min_atoms_b * min_atoms_c)
              + util::Format()( COUNTS_PACK_C >= min_atoms_a * min_atoms_c)
            : std::string()
          )
        + std::string( ADJACENT ? "Adj" : "NAdj")
        + std::string( size_t( 1), RHS ? 'R' : 'L')
      );
    }

    //! @brief hash a given packing into a size_t for speed. The size_t will have the range 0 - 223 ( 7 * 8 * 4)
    size_t ProteinModelSSETripletChirality::GetPackingTripletNumber
    (
      const biol::SSType &TYPE_A,
      const biol::SSType &TYPE_B,
      const biol::SSType &TYPE_C,
      const assemble::SSEGeometryPacking::OrientationEnum &PACKING_AB,
      const assemble::SSEGeometryPacking::OrientationEnum &PACKING_BC,
      const assemble::SSEGeometryPacking::OrientationEnum &PACKING_AC,
      const size_t &COUNTS_PACK_A,
      const size_t &COUNTS_PACK_B,
      const size_t &COUNTS_PACK_C,
      const bool   &ADJACENT,
      const bool   &RHS
    ) const
    {
      const size_t min_atoms_a
      (
        TYPE_A == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
      );
      const size_t min_atoms_b
      (
        TYPE_B == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
      );
      const size_t min_atoms_c
      (
        TYPE_C == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
      );
      if( m_ConsiderWhichSSEsInContact)
      {
        bool ncontact_ab( COUNTS_PACK_A < min_atoms_a * min_atoms_b);
        bool ncontact_ac( COUNTS_PACK_C < min_atoms_a * min_atoms_c);
        bool ncontact_bc( COUNTS_PACK_B < min_atoms_b * min_atoms_c);
        if( ( ncontact_ab && ( ncontact_ac || ncontact_bc)) || ( ncontact_ac && ncontact_bc))
        {
          ncontact_ab = !ncontact_ab;
          ncontact_ac = !ncontact_ac;
          ncontact_bc = !ncontact_bc;
        }
        return
           ( TYPE_A == biol::GetSSTypes().STRAND ? 512 : 0)
         | ( TYPE_B == biol::GetSSTypes().STRAND ? 256 : 0)
         | ( TYPE_C == biol::GetSSTypes().STRAND ? 128 : 0)
         | ( PACKING_AB == assemble::SSEGeometryPacking::e_Parallel ? 64 : 0)
         | ( PACKING_BC == assemble::SSEGeometryPacking::e_Parallel ? 32 : 0)
         | ( PACKING_AC == assemble::SSEGeometryPacking::e_Parallel ? 16 : 0)
         | ( ADJACENT ? 8 : 0)
         | ( ncontact_ab ? 2 : ( ncontact_bc ? 4 : ( ncontact_ac ? 6 : 0)))
         | ( RHS ? 1 : 0);
      }
      return ( TYPE_A == biol::GetSSTypes().STRAND ? 128 : 0)
             | ( TYPE_B == biol::GetSSTypes().STRAND ? 64 : 0)
             | ( TYPE_C == biol::GetSSTypes().STRAND ? 32 : 0)
             | ( PACKING_AB == assemble::SSEGeometryPacking::e_Parallel ? 16 : 0)
             | ( PACKING_BC == assemble::SSEGeometryPacking::e_Parallel ? 8 : 0)
             | ( PACKING_AC == assemble::SSEGeometryPacking::e_Parallel ? 4 : 0)
             | ( ADJACENT ? 2 : 0)
             | ( RHS ? 1 : 0);
    }

    //! @brief convert a packing triplet number to a string
    std::string ProteinModelSSETripletChirality::GetPackingTripletString( const size_t &HASH_NUMBER) const
    {
      static storage::VectorND< 2, storage::Vector< std::string> > s_all_hashes;
      storage::Vector< std::string> &hashes( s_all_hashes( m_ConsiderWhichSSEsInContact ? 1 : 0));
      if( hashes.IsEmpty())
      {
        hashes.Resize( m_ConsiderWhichSSEsInContact ? s_NumberHashesContactSplit : s_NumberHashesNoContactSplit);
        const storage::Vector< biol::SSType> sstypes
        (
          storage::Vector< biol::SSType>::Create( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND)
        );
        const storage::Vector< assemble::SSEGeometryPacking::OrientationEnum> oris
        (
          storage::Vector< assemble::SSEGeometryPacking::OrientationEnum>::Create
          (
            assemble::SSEGeometryPacking::e_AntiParallel,
            assemble::SSEGeometryPacking::e_Parallel
          )
        );
        storage::Vector< storage::VectorND< 3, size_t> > counts
        (
          storage::Vector< storage::VectorND< 3, size_t> >::Create
          (
            storage::VectorND< 3, size_t>( 1000, 1000, 1000),
            storage::VectorND< 3, size_t>( 0, 1000, 1000),
            storage::VectorND< 3, size_t>( 1000, 0, 1000),
            storage::VectorND< 3, size_t>( 1000, 1000, 0)
          )
        );
        if( !m_ConsiderWhichSSEsInContact)
        {
          counts.Resize( 1);
        }
        // iterate through each allowed combination of ss types, orientation, adjaceny, and counts
        for( auto itr_a( sstypes.Begin()), ss_end( sstypes.End()); itr_a != ss_end; ++itr_a)
        {
          for( auto itr_b( sstypes.Begin()); itr_b != ss_end; ++itr_b)
          {
            for( auto itr_c( sstypes.Begin()); itr_c != ss_end; ++itr_c)
            {
              for( auto itr_d( oris.Begin()), ori_end( oris.End()); itr_d != ori_end; ++itr_d)
              {
                for( auto itr_e( oris.Begin()); itr_e != ori_end; ++itr_e)
                {
                  for( auto itr_f( oris.Begin()); itr_f != ori_end; ++itr_f)
                  {
                    for( auto itr_g( counts.Begin()), itr_g_end( counts.End()); itr_g != itr_g_end; ++itr_g)
                    {
                      for( int adjacent( 0); adjacent < 2; ++adjacent)
                      {
                        for( int rhs( 0); rhs < 2; ++rhs)
                        {
                          // associate the hash number with the string
                          hashes
                          (
                            GetPackingTripletNumber
                            (
                              *itr_a, *itr_b, *itr_c,
                              *itr_d, *itr_e, *itr_f,
                              itr_g->First(), itr_g->Second(), itr_g->Third(),
                              bool( adjacent),
                              bool( rhs)
                            )
                          ) = GetPackingTripletHash
                              (
                                *itr_a, *itr_b, *itr_c,
                                *itr_d, *itr_e, *itr_f,
                                itr_g->First(), itr_g->Second(), itr_g->Third(),
                                bool( adjacent),
                                bool( rhs)
                              );
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
      return hashes( HASH_NUMBER);
    }

    //! @brief convert a packing triplet number to a string
    linal::Vector< double> ProteinModelSSETripletChirality::ComputeExpectedCounts
    (
      const linal::Matrix< double> &ADJ_CONTACT_SEEN, // Rows - ss triplet-type, Cols - contact type
      const linal::Matrix< double> &DIST_CONTACT_SEEN, // Rows - ss triplet-type, Cols - contact type
      const linal::Matrix< double> &PARALLEL_PROBS, // Rows - ss pair-type, Cols - SSE distance (1,2,3+)
      const linal::Vector< size_t> &SSE_TRIPLET_COUNTS_ADJACENT, // size - 7, sse triplet type, lowest value indicates first sse
      const linal::Vector< size_t> &SSE_TRIPLET_COUNTS // size - 7, sse triplet type, lowest value indicates first sse, 0 - helix, 1 - strand,
    ) const
    {
      linal::Vector< double> output( GetNumberHashes(), double( 0.0));
      const storage::Vector< biol::SSType> sstypes
      (
        storage::Vector< biol::SSType>::Create( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND)
      );
      const storage::Vector< assemble::SSEGeometryPacking::OrientationEnum> oris
      (
        storage::Vector< assemble::SSEGeometryPacking::OrientationEnum>::Create
        (
          assemble::SSEGeometryPacking::e_AntiParallel,
          assemble::SSEGeometryPacking::e_Parallel
        )
      );
      // compute probabilities for each orientation triplet. This arises from the integral
      // sin(phi)*(1/2-phi/(2*pi))/2 from 0 to pi/2 (in the case of odd # of parallel elements, the more common case)
      // and integral( sin(phi)*phi/(2*pi)/2 from 0 to pi/2 for the even # of parallel elements
      // So the index is the number of parallel sses mod two
      const double orientation_prob[ 2] = { 0.25 / math::g_Pi, 0.25 * ( math::g_Pi - 1.0) / math::g_Pi};

      storage::Vector< storage::VectorND< 3, size_t> > counts
      (
        storage::Vector< storage::VectorND< 3, size_t> >::Create
        (
          storage::VectorND< 3, size_t>( 1000, 1000, 1000),
          storage::VectorND< 3, size_t>( 0, 1000, 1000),
          storage::VectorND< 3, size_t>( 1000, 0, 1000),
          storage::VectorND< 3, size_t>( 1000, 1000, 0)
        )
      );
      if( !m_ConsiderWhichSSEsInContact)
      {
        counts.Resize( 1);
      }
      linal::Vector< double> expected_contact_probs( 8);
      linal::Vector< double> expected_ori_probs( 8);
      linal::Vector< double> relevant_contact_probs( counts.GetSize());
      for( int adjacent( 0); adjacent < 2; ++adjacent)
      {
        const linal::Matrix< double> &seen_matrix( adjacent ? ADJ_CONTACT_SEEN : DIST_CONTACT_SEEN);
        const int adj_index( adjacent ? 0 : 2);
        const int para_index( adjacent ? 1 : 2);
        const linal::Vector< size_t> &counts_vector( adjacent ? SSE_TRIPLET_COUNTS_ADJACENT : SSE_TRIPLET_COUNTS);
        // iterate through each allowed combination of ss types, orientation, adjaceny, and counts
        size_t int_type_a( 0);
        for( auto itr_a( sstypes.Begin()), ss_end( sstypes.End()); itr_a != ss_end; ++itr_a, ++int_type_a)
        {
          size_t int_type_b( 0);
          for( auto itr_b( sstypes.Begin()); itr_b != ss_end; ++itr_b, ++int_type_b)
          {
            const size_t int_type_ab( int_type_a + int_type_b * 2);
            size_t int_type_c( 0);
            for( auto itr_c( sstypes.Begin()); itr_c != ss_end; ++itr_c, ++int_type_c)
            {
              const size_t int_type_ac( int_type_a + int_type_c * 2);
              const size_t int_type_bc( int_type_b + int_type_c * 2);
              const size_t int_type_abc( int_type_a + 2 * ( int_type_b + 2 * int_type_c));
              expected_contact_probs = seen_matrix.GetRow( int_type_abc);
              expected_contact_probs.SetToSum( 1.0);
              for( size_t c( 0); c < size_t( 8); ++c)
              {;
                expected_ori_probs( c) =
                  ( c & 4 ? PARALLEL_PROBS( int_type_ab, adj_index) : 1.0 - PARALLEL_PROBS( int_type_ab, adj_index))
                  * ( c & 2 ? PARALLEL_PROBS( int_type_bc, adj_index) : 1.0 - PARALLEL_PROBS( int_type_bc, adj_index))
                  * ( c & 1 ? PARALLEL_PROBS( int_type_ac, para_index) : 1.0 - PARALLEL_PROBS( int_type_ac, para_index));
                expected_ori_probs( c) *= orientation_prob[ ( ( ( c & 4) >> 2) + ( ( c & 2) >> 1) + ( c & 1)) & 1];
              }
              expected_ori_probs.SetToSum( 1.0);
              if( !m_ConsiderWhichSSEsInContact)
              {
                relevant_contact_probs( 0) =
                  expected_contact_probs( 3) + expected_contact_probs( 5)
                  + expected_contact_probs( 6) + expected_contact_probs( 7);
              }
              else
              {
                relevant_contact_probs( 0) = expected_contact_probs( 7);
                relevant_contact_probs( 1) = expected_contact_probs( 3);
                relevant_contact_probs( 2) = expected_contact_probs( 5);
                relevant_contact_probs( 3) = expected_contact_probs( 6);
              }

              const double count( counts_vector( int_type_abc));
              size_t ori_index( 0);
              for( auto itr_d( oris.Begin()), ori_end( oris.End()); itr_d != ori_end; ++itr_d)
              {
                for( auto itr_e( oris.Begin()); itr_e != ori_end; ++itr_e)
                {
                  for( auto itr_f( oris.Begin()); itr_f != ori_end; ++itr_f, ++ori_index)
                  {
                    for
                    (
                      size_t count_vec_i( 0), n_counts( relevant_contact_probs.GetSize());
                      count_vec_i < n_counts;
                      ++count_vec_i
                    )
                    {
                      const storage::VectorND< 3, size_t> &cvec( counts( count_vec_i));
                      double expected_val( relevant_contact_probs( count_vec_i) * count / 2.0 * expected_ori_probs( ori_index));
                      output
                      (
                        GetPackingTripletNumber
                        (
                          *itr_a, *itr_b, *itr_c,
                          *itr_d, *itr_e, *itr_f,
                          cvec.First(), cvec.Second(), cvec.Third(),
                          bool( adjacent),
                          false
                        )
                      ) = expected_val;
                      output
                      (
                        GetPackingTripletNumber
                        (
                          *itr_a, *itr_b, *itr_c,
                          *itr_d, *itr_e, *itr_f,
                          cvec.First(), cvec.Second(), cvec.Third(),
                          bool( adjacent),
                          true
                        )
                      ) = expected_val;
                    }
                  }
                }
              }
            }
          }
        }
      }
      return output;
    }

  //////////////
  // operator //
  //////////////

    //! @brief add statistics about the given protein and chains to the internal table and histograms
    //! @param ENSEMBLE the protein ensemble the statistics of which to be computed
    std::string ProteinModelSSETripletChirality::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the protein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

      // write statistics
      std::ostringstream stream;

      // min probability to score two residues as being an almost-certain contact
      const double min_p( 0.04);

      const size_t n_hashes( m_ConsiderWhichSSEsInContact ? s_NumberHashesContactSplit : s_NumberHashesNoContactSplit);
      // initializes maps for storing histograms of distances between all types of amino acid pairs
      linal::Vector< double> triplet_packing_counts( n_hashes, 0.0);
      storage::Vector< size_t> triplet_packing_sse_counts( n_hashes, size_t( 0));
      storage::Vector< std::string> triplet_strings( n_hashes);
      linal::Vector< double> triplet_packing_counts_nic( n_hashes, 0.0);

      if( m_OutputPerProtein)
      {
        stream << "storage::Table<double>\tIsNative\t";
        for( size_t hash_id( 0); hash_id < n_hashes; ++hash_id)
        {
          triplet_strings( hash_id) = GetPackingTripletString( hash_id);
          stream << triplet_strings( hash_id) << '\t';
        }
        stream << '\n';
      }
      linal::Matrix< double> adj_seen( 8, 8, double( 0));
      linal::Matrix< double> dist_seen( 8, 8, double( 0));

      storage::Vector< math::RunningAverage< double> > adj_fract_parallel( 4);
      storage::Vector< math::RunningAverage< double> > semiadj_fract_parallel( 4);
      storage::Vector< math::RunningAverage< double> > dist_fract_parallel( 4);
      linal::Vector< size_t> total_sses_seen( 8, size_t( 0));
      linal::Vector< size_t> total_sses_seen_adj( 8, size_t( 0));

      storage::Set< biol::AtomType> types_of_interest( biol::GetAtomTypes().GetFirstSidechainAtomTypes());
      // Create a matrix that will hold the x,y, and z coordinates for terminii of the first & second strands
      linal::Matrix3x3< double> xyz_coordinates( 0.0); // make a matrix of size 3 X 3
      for( auto itr_ensemble( ENSEMBLE.Begin()), itr_ensemble_end( ENSEMBLE.End()); itr_ensemble != itr_ensemble_end; ++itr_ensemble)
      {
        // get current protein model
        util::ShPtr< assemble::ProteinModel> sp_current_model( *itr_ensemble);
        const assemble::ProteinModel &protein_model( *sp_current_model);

        // get current pdb name before all the loops start
        const util::ShPtr< util::Wrapper< std::string> > &model_name_ptr(
          protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile));
        const std::string model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");

        // instantiate sequences CaCb from pdb for each chain one
        BCL_MessageStd( "computing chirality and packing for " + model_name);

        for
        (
          auto itr_chain( ( *itr_ensemble)->GetChains().Begin()), itr_chain_end( ( *itr_ensemble)->GetChains().End());
          itr_chain != itr_chain_end;
          ++itr_chain
        )
        {
          // need at least three sses in the chain
          if( ( *itr_chain)->GetNumberSSEs() < 2)
          {
            continue;
          }
          // collect all non-coil sse's
          const util::SiPtrVector< const assemble::SSE> structured_sses
          (
            ( *itr_chain)->GetSSEs( storage::Set< biol::SSType>( biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND))
          );

          if( structured_sses.GetSize() < 2)
          {
            continue;
          }
          storage::Vector< size_t> atoms_in_sse;
          for( auto itr_sse( structured_sses.Begin()), itr_sse_end( structured_sses.End()); itr_sse != itr_sse_end; ++itr_sse)
          {
            auto new_atoms( ( *itr_sse)->GetAtoms( types_of_interest));
            atoms_in_sse.PushBack( new_atoms.GetSize());
          }
          assemble::VoxelGridAA interactions_detector( m_InteractionDistance);
          linal::Matrix< float> interactions_matrix
          (
            interactions_detector.GetSSEInteractionMatrix
            (
              structured_sses,
              ( *itr_chain)->GetAminoAcids(),
              2,
              m_InteractionDistance,
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

          // handle SSE-pair interaction-weight and parallel/anti-parallel bias
          size_t a( 0);
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              itr_a( structured_sses.Begin()),
              itr_end( structured_sses.End());
            itr_a != itr_end;
            ++itr_a, ++a
          )
          {
            const biol::SSType &type_a( ( *itr_a)->GetType());
            const size_t a_type( type_a == biol::GetSSTypes().STRAND ? 1 : 0);
            size_t b( a + 1);
            const size_t min_atoms_a
            (
              type_a == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
            );
            for
            (
              util::SiPtrVector< const assemble::SSE>::const_iterator itr_b( itr_a + 1);
              itr_b != itr_end;
              ++itr_b, ++b
            )
            {
              const size_t b_type( ( *itr_b)->GetType() == biol::GetSSTypes().STRAND ? 2 : 0);
              const biol::SSType &type_b( ( *itr_b)->GetType());
              const size_t min_atoms_b
              (
                type_b == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
              );

              const size_t interactions_ab( interactions_matrix( a, b));
              // for every pair of SSEs
              const assemble::SSEGeometryPacking::OrientationEnum &packing_ab
              (
                GetCacheOrientation( orientations, a, b, **itr_a, **itr_b)
              );

              const double is_parallel( packing_ab == assemble::SSEGeometryPacking::e_Parallel ? 1.0 : 0.0);
              if( itr_b == itr_a + 1)
              {
                adj_fract_parallel( b_type | a_type) += is_parallel;
              }
              else
              {
                if( itr_b == itr_a + 2)
                {
                  semiadj_fract_parallel( b_type | a_type) += is_parallel;
                }
                dist_fract_parallel( b_type | a_type) += is_parallel;
              }

              util::SiPtrVector< const assemble::SSE>::const_iterator itr_c( itr_b + 1);
              if( itr_c == itr_end)
              {
                continue;
              }

              size_t c( b + 1);

              for( ; itr_c != itr_end; ++itr_c, ++c)
              {
                const biol::SSType &type_c( ( *itr_c)->GetType());
                const size_t c_type( type_c == biol::GetSSTypes().STRAND ? 4 : 0);
                const size_t min_atoms_c
                (
                  type_c == biol::GetSSTypes().STRAND ? m_MinAtomsInContactStrand : m_MinAtomsInContactHelix
                );
                const size_t interactions_ac( interactions_matrix( a, c));
                const size_t interactions_bc( interactions_matrix( b, c));
                const size_t adj_type
                (
                  ( interactions_ac < min_atoms_a * min_atoms_c ? 0 : 1)
                  + ( interactions_ab < min_atoms_a * min_atoms_b ? 0 : 4)
                  + ( interactions_bc < min_atoms_b * min_atoms_c ? 0 : 2)
                );
                // adjacent
                if( c == a + 2)
                {
                  ++total_sses_seen_adj( a_type | b_type | c_type);
                  adj_seen( a_type | b_type | c_type, adj_type) += 1.0;
                }
                else
                {
                  ++total_sses_seen( a_type | b_type | c_type);
                  dist_seen( a_type | b_type | c_type, adj_type) += 1.0;
                }

                const assemble::SSEGeometryPacking::OrientationEnum &packing_ac
                (
                  GetCacheOrientation( orientations, a, c, **itr_a, **itr_c)
                );
                const assemble::SSEGeometryPacking::OrientationEnum &packing_bc
                (
                  GetCacheOrientation( orientations, b, c, **itr_b, **itr_c)
                );
                const linal::Vector3D root_position( ( *itr_b)->GetCenter());

                xyz_coordinates.GetRow( 0).CopyValues( ( *itr_a)->GetMainAxis().GetStartPoint() - root_position);
                xyz_coordinates.GetRow( 1).CopyValues( ( *itr_a)->GetMainAxis().GetEndPoint() - root_position);
                xyz_coordinates.GetRow( 2).CopyValues( ( *itr_c)->GetCenter() - root_position);
                // Calculate the determinant of a 3 X 3 matrix that has rows sorted in descending order of priority.
                // Opposite orders will have opposite signs. The sign is assigned to a value of 1 and returned as the
                // value of the stereocenter.
                const float determinant( xyz_coordinates.Determinant());

                size_t triplet_number
                (
                  GetPackingTripletNumber
                  (
                    type_a,
                    type_b,
                    type_c,
                    packing_ab,
                    packing_bc,
                    packing_ac,
                    interactions_ab,
                    interactions_bc,
                    interactions_ac,
                    itr_a + 1 == itr_b && itr_b + 1 == itr_c,
                    determinant >= 0
                  )
                );

                if
                (
                  ( interactions_ac < min_atoms_a * min_atoms_c ? 1 : 0)
                  + ( interactions_ab < min_atoms_a * min_atoms_b ? 1 : 0)
                  + ( interactions_bc < min_atoms_b * min_atoms_c ? 1 : 0)
                  > 1
                )
                {
                  ++triplet_packing_counts_nic( triplet_number);
                  continue;
                }

                BCL_MessageVrb
                (
                  GetPackingTripletString( triplet_number)
                  + ( *itr_a)->GetIdentification() + " " + ( *itr_b)->GetIdentification() + " " + ( *itr_c)->GetIdentification()
                  + " " + util::Format()( interactions_ab)
                  + " " + util::Format()( interactions_ac)
                  + " " + util::Format()( interactions_bc)
                );
                triplet_packing_counts( triplet_number)
                    += double( interactions_ac + interactions_ab + interactions_bc) / double( atoms_in_sse( a) + atoms_in_sse( b) + atoms_in_sse( c));
                triplet_packing_sse_counts( triplet_number) += 1;
              }
            }
          }
        }
        if( m_OutputPerProtein && triplet_packing_counts.Sum())
        {
          stream << model_name << "\t1\t";
          for( size_t hash_id( 0); hash_id < n_hashes; ++hash_id)
          {
            stream << triplet_packing_counts( hash_id) << '\t';
          }
          stream << '\n';
          stream << model_name << "Mirror\t-1\t";
          for( size_t hash_id( 0); hash_id < n_hashes; ++hash_id)
          {
            stream << triplet_packing_counts( hash_id ^ size_t( 1)) << '\t';
          }
          stream << '\n';
          triplet_packing_counts = 0.0;
          triplet_packing_sse_counts.SetAllElements( 0);
        }
      }

      if( !m_OutputPerProtein)
      {
        stream << this->GetLabel().ToString() << '\n';
        stream << "# Triplet-Hash\tWeightedSSEs\tSSEs\tExpectedSSEsInContact\n";
        auto itr_sse( triplet_packing_sse_counts.Begin());
        size_t hash_id( 0);
        linal::Matrix< double> parallel_probs( 4, 3);
        for( int i( 0); i < 4; ++i)
        {
          parallel_probs( i, 0) = adj_fract_parallel( i).GetAverage();
          parallel_probs( i, 1) = semiadj_fract_parallel( i).GetAverage();
          parallel_probs( i, 2) = dist_fract_parallel( i).GetAverage();
        }

        linal::Vector< double> expected_counts
        (
          ComputeExpectedCounts( adj_seen, dist_seen, parallel_probs, total_sses_seen_adj, total_sses_seen)
        );
        for
        (
          auto itr( triplet_packing_counts.Begin()),
               itr_end( triplet_packing_counts.End());
          itr != itr_end;
          ++itr, ++itr_sse, ++hash_id
        )
        {
          stream << GetPackingTripletString( hash_id) << ' '
                 << *itr << ' '
                 << *itr_sse << ' '
                 << expected_counts( hash_id)
                 << '\n';
        }
      }

      return stream.str();

    } // end of operator

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinModelSSETripletChirality::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Computes chirality-based entropies");
      parameters.AddInitializer
      (
        "distance",
        "max distance between atoms to be considered in contact",
        io::Serialization::GetAgent( &m_InteractionDistance),
        util::Format()( ProteinModelSSETripletChirality().m_InteractionDistance)
      );
      parameters.AddInitializer
      (
        "min sc atoms strand",
        "minimum side chains atoms betwen contacting SSEs",
        io::Serialization::GetAgent( &m_MinAtomsInContactStrand),
        util::Format()( ProteinModelSSETripletChirality().m_MinAtomsInContactStrand)
      );
      parameters.AddInitializer
      (
        "min sc atoms helix",
        "minimum side chains atoms betwen contacting SSEs",
        io::Serialization::GetAgent( &m_MinAtomsInContactHelix),
        util::Format()( ProteinModelSSETripletChirality().m_MinAtomsInContactHelix)
      );
      parameters.AddInitializer
      (
        "consider which sses contact",
        "whether to consider which SSEs are in contact. This is rather important for several SSE types, but only rarely"
        "affects the preferred handedness",
        io::Serialization::GetAgent( &m_ConsiderWhichSSEsInContact),
        util::Format()( ProteinModelSSETripletChirality().m_ConsiderWhichSSEsInContact)
      );
      parameters.AddInitializer
      (
        "output per protein",
        "If true, output will be raw counts per protein, suitable for use in MinimizeScoreWeightset. Two rows will be "
        "output per protein model. One contains the real chirality, the other the mirror image. IsNative score will be "
        "output as 1 for the native model, -1 for the mirror",
        io::Serialization::GetAgent( &m_OutputPerProtein),
        util::Format()( ProteinModelSSETripletChirality().m_OutputPerProtein)
      );
      return parameters;
    }
  } // namespace scorestat
} // namespace bcl

