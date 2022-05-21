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
#include "scorestat/bcl_scorestat_protein_model_packing.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "assemble/bcl_assemble_voxel_grid_aa.h"
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
#include "storage/bcl_storage_table.h"
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
    const util::SiPtr< const util::ObjectInterface> ProteinModelPacking::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new ProteinModelPacking())
    );

    //! @brief Category as string
    //! @param CATEGORY the desired category
    //! @return the string for the category
    const std::string &ProteinModelPacking::GetCategoryName( const Category &CATEGORY)
    {
      static std::string s_names[ s_NumberCategories + 1] =
      {
        "AdjacentInContact",
        "AdjacentNotInContact",
        "AdjacentParallel",
        "AdjacentAntiParallel",
        "OneSSEApartParallel",
        "OneSSEApartAntiParallel",
        "TwoSSEApartParallel",
        "TwoSSEApartAntiParallel",
        "FarApartParallel",
        "FarApartAntiParallel",
        "WeakInteraction",
        "ModerateInteraction",
        "StrongInteraction",
        "BackgroundWeakInteraction",
        "BackgroundModerateInteraction",
        "BackgroundStrongInteraction",
        "BackgroundParallel",
        "BackgroundAntiParallel",
        GetStaticClassName< Category>()
      };
      return s_names[ CATEGORY];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ProteinModelPacking::ProteinModelPacking
    (
      const double &MIN_INTERACTION_DISTANCE,
      const size_t &MIN_ATOMS_IN_CONTACT
    ) :
      m_InteractionDistance( MIN_INTERACTION_DISTANCE),
      m_MinAtomsInContact( MIN_ATOMS_IN_CONTACT)
    {
    }

    //! @brief virtual copy constructor
    ProteinModelPacking *ProteinModelPacking::Clone() const
    {
      return new ProteinModelPacking( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinModelPacking::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &ProteinModelPacking::GetOutFilePostfix() const
    {
      static const std::string s_name( "sse_packing_type");
      return s_name;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &ProteinModelPacking::GetAlias() const
    {
      static const std::string s_name( "SSEPackingType");
      return s_name;
    }

    //! @brief hash a given packing. The given string is intended to be fast to hash, not necessarily easily readable
    void ProteinModelPacking::AddPackingType
    (
      const assemble::SSEGeometryPacking &PACKING,
      const bool &IS_IN_CONTACT,
      const size_t &N_SSES_APART,
      const bool &IS_BACKGROUND,
      const bool &ORIENTATION_COULD_BE_OPPOSITE,
      linal::VectorInterface< size_t> &COUNTS
    )
    {
      if( !IS_BACKGROUND)
      {
        // adjacency contacts. Bins 0 - 1
        // Orientation contacts. Bins 2 - 5 (Adjacent Parallel, Adjacent Anti-Parallel, Non-Adj Parallel, Non-Adj Anti-Parallel)
        if( !N_SSES_APART)
        {
          ++COUNTS( IS_IN_CONTACT ? e_AdjacentInContact : e_AdjacentNotInContact);
        }
        if( !IS_IN_CONTACT)
        {
          return;
        }
        if( ORIENTATION_COULD_BE_OPPOSITE)
        {
          ++COUNTS
          (
            std::min( N_SSES_APART, size_t( 3)) * 2 + size_t( e_AdjacentParallel)
            + ( PACKING.GetOrientation() == assemble::SSEGeometryPacking::e_AntiParallel ? 1 : 0)
          );
        }
        ++COUNTS
        (
          PACKING.GetInteractionWeight() < 0.5
          ? e_WeakInteraction
          : (
              PACKING.GetInteractionWeight() < 0.95
              ? e_ModerateInteraction
              : e_StrongInteraction
            )
        );
      }
      else
      {
        ++COUNTS
        (
          PACKING.GetInteractionWeight() < 0.5
          ? e_BackgroundWeakInteraction
          : (
              PACKING.GetInteractionWeight() < 0.95
              ? e_BackgroundModerateInteraction
              : e_BackgroundStrongInteraction
            )
        );
        ++COUNTS
        (
          PACKING.GetOrientation() == assemble::SSEGeometryPacking::e_Parallel
          ? e_BackgroundParallel
          : e_BackgroundAntiParallel
        );
      }
    }

  //////////////
  // operator //
  //////////////

    //! @brief add statistics about the given protein and chains to the internal table and histograms
    //! @param ENSEMBLE the protein ensemble the statistics of which to be computed
    std::string ProteinModelPacking::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // make sure the protein ensemble is not empty
      BCL_Assert( ENSEMBLE.GetSize() != 0, "protein ensemble is empty");

      // min probability to score two residues as being an almost-certain contact
      const double min_p( 0.04);

      // initializes maps for storing histograms of distances between all types of amino acid pairs
      storage::Vector< linal::Vector< size_t> > packing_info_by_contact_type
      (
        contact::GetTypes().GetEnumCount(),
        linal::Vector< size_t>( size_t( s_NumberCategories), size_t( 0))
      );
      storage::Vector< linal::Vector< size_t> > contact_type_by_ss_types
      (
        size_t( 4),
        linal::Vector< size_t>( size_t( contact::GetTypes().GetEnumCount()), size_t( 0))
      );
      storage::Vector< linal::Vector< size_t> > contact_type_by_ss_types_bg
      (
        size_t( 4),
        linal::Vector< size_t>( size_t( contact::GetTypes().GetEnumCount()), size_t( 0))
      );
      const long max_helix_helix_contacts( 8), max_strand_strand_contacts( 6);
      const long max_helix_strand_contacts( 6), max_strand_helix_contacts( 3);

      linal::Matrix< double> helix_adjacency_counts( max_helix_helix_contacts + 1, max_helix_strand_contacts + 1, double( 0.0));
      linal::Matrix< double> strand_adjacency_counts( max_strand_helix_contacts + 1, max_strand_strand_contacts + 1, double( 0.0));
      linal::Matrix< double> helix_adjacency_plaus( helix_adjacency_counts), strand_adjacency_plaus( strand_adjacency_counts);
      storage::Vector< double> counts_contacts_possible( size_t( 4));
      storage::Vector< double> counts_contacts_actual( size_t( 4));
      // for each sse, compute average contact type of contacting SSEs, weighted by # of AAs in contact
      // Also compute P(Helix|Helix) In contact, P(Strand|Strand)
      // # of helix-strand contacts / (min( # of helix * min(# of strands,4), min(# of helices,2)*# of strands))
      // # of helix-helix contacts * 2 / (# of helices * min(# of helices - 1,4))
      // # of strand-strand contacts * 2 / (# of strands * min(# of strand - 1,4))

      storage::Set< biol::AtomType> types_of_interest( biol::GetAtomTypes().GetFirstSidechainAtomTypes());

      const std::string sse_name[ 2] = { "HELIX", "STRAND"};

      size_t n_runs( 10000);
      score::AAPairHiResClash clash;
      for( int config( 0); config < 4; ++config)
      {
        const size_t hs_a( config & 1 ? 1 : 0);
        const size_t hs_b( config & 2 ? 1 : 0);
        BCL_MessageStd( "Computing background for " + sse_name[ hs_a] + " and " + sse_name[ hs_b]);

        for( size_t run_n( 0), n_runs_l( n_runs); run_n < n_runs_l; ++run_n)
        {
          const size_t strand_l( 6), helix_l( 10);
          storage::VectorND< 2, storage::VectorND< 3, assemble::SSE> > helix_strand_res;
          int o_id( 1);
          int pdb_atom_id( 1);
          for( size_t o( 0); o < size_t( 3); ++o, o_id += 2)
          {
            util::ShPtrVector< biol::AABase> strand_aas( strand_l);
            for( size_t i( 0); i < strand_l; ++i, ++o_id)
            {
              double ran_choice( random::GetGlobalRandom().Random( 1.0));
              biol::AAType type( biol::GetAATypes().ALA);
              ran_choice -= type->GetAAProperty( biol::AATypeData::e_NaturalPrevalence);
              for( int type_x( 1); ran_choice > 0.0; type_x += 1)
              {
                type = biol::AAType( type_x);
                ran_choice -= type->GetAAProperty( biol::AATypeData::e_NaturalPrevalence);
              }
              util::ShPtr< biol::AAData> aadata( new biol::AAData( type, o_id, o_id, 'A'));
              strand_aas( i) = util::ShPtr< biol::AABase>( new biol::AABackBone( aadata));
              strand_aas( i)->SetAtom( biol::Atom( biol::GetAtomTypes().N, pdb_atom_id++));
              strand_aas( i)->SetAtom( biol::Atom( biol::GetAtomTypes().C, pdb_atom_id++));
              strand_aas( i)->SetAtom( biol::Atom( biol::GetAtomTypes().CA, pdb_atom_id++));
              strand_aas( i)->SetAtom( biol::Atom( biol::GetAtomTypes().O, pdb_atom_id++));
              strand_aas( i)->SetAtom( biol::Atom( type->GetFirstSidechainAtomType(), pdb_atom_id++));
            }
            o_id += 4;
            util::ShPtrVector< biol::AABase> helix_aas( helix_l);
            for( size_t j( 0); j < helix_l; ++j, ++o_id)
            {
              double ran_choice( random::GetGlobalRandom().Random( 1.0));
              biol::AAType type( biol::GetAATypes().ALA);
              ran_choice -= type->GetAAProperty( biol::AATypeData::e_NaturalPrevalence);
              for( int type_x( 1); ran_choice > 0.0; type_x += 1)
              {
                type = biol::AAType( type_x);
                ran_choice -= type->GetAAProperty( biol::AATypeData::e_NaturalPrevalence);
              }
              util::ShPtr< biol::AAData> aadata( new biol::AAData( type, o_id, o_id, 'A'));
              helix_aas( j) = util::ShPtr< biol::AABase>( new biol::AABackBone( aadata));
              helix_aas( j)->SetAtom( biol::Atom( biol::GetAtomTypes().N, pdb_atom_id++));
              helix_aas( j)->SetAtom( biol::Atom( biol::GetAtomTypes().C, pdb_atom_id++));
              helix_aas( j)->SetAtom( biol::Atom( biol::GetAtomTypes().CA, pdb_atom_id++));
              helix_aas( j)->SetAtom( biol::Atom( biol::GetAtomTypes().O, pdb_atom_id++));
              helix_aas( j)->SetAtom( biol::Atom( type->GetFirstSidechainAtomType(), pdb_atom_id++));
            }
            o_id += 4;
            helix_strand_res( 0)( o) = assemble::SSE( biol::AASequence( helix_aas), biol::GetSSTypes().HELIX);
            helix_strand_res( 1)( o) = assemble::SSE( biol::AASequence( strand_aas), biol::GetSSTypes().STRAND);
            helix_strand_res( 0)( o).SetToIdealConformationAtOrigin();
            helix_strand_res( 1)( o).SetToIdealConformationAtOrigin();
          }
          assemble::SSE &sse_a( helix_strand_res( hs_a)( 0));
          assemble::SSE &sse_b( helix_strand_res( hs_b)( 1));

          sse_a.SetToIdealConformationAtOrigin();
          sse_b.SetToIdealConformationAtOrigin();
          sse_a.Rotate( math::RotationMatrix3D().SetRand());
          sse_b.Rotate( math::RotationMatrix3D().SetRand());
          sse_a.Translate( linal::Vector3D().SetRandomTranslation( 20.0));
          sse_b.Translate( linal::Vector3D().SetRandomTranslation( 20.0));
          if( clash( sse_a, sse_b))
          {
            --run_n;
            continue;
          }

          assemble::VoxelGridAA interactions_detector( m_InteractionDistance);
          util::SiPtrVector< const assemble::SSE> structured_sses( util::SiPtrVector< const assemble::SSE>::Create( sse_a, sse_b));
          util::SiPtrVector< const biol::AABase> aas( sse_a.Begin(), sse_a.End());
          aas.Append( util::SiPtrVector< const biol::AABase>( sse_b.Begin(), sse_b.End()));
          linal::Matrix< float> interactions_matrix
          (
            interactions_detector.GetSSEInteractionMatrix( structured_sses, aas, 2, m_InteractionDistance, false, min_p, false)
          );

          if( interactions_matrix( 0, 1) < m_MinAtomsInContact)
          {
            --run_n;
            continue;
          }
          // at this point we know that we have a Beta-Alpha-Beta with no intervening SSEs. Next, check that the strands
          // point in the same general direction
          assemble::SSEGeometryPacking packing_b( sse_a, sse_b, 0.0);
          if( !packing_b.GetContactType().IsDefined())
          {
            continue;
          }

          AddPackingType( packing_b, true, 1, true, true, packing_info_by_contact_type( packing_b.GetContactType()));
          AddPackingType( packing_b, true, 1, true, true, packing_info_by_contact_type( contact::GetTypes().Reverse( packing_b.GetContactType())));
          ++contact_type_by_ss_types_bg( config)( packing_b.GetContactType());
          if( !( run_n % ( n_runs_l / 100)))
          {
            util::GetLogger().LogStatus( util::Format()( run_n * 100 / n_runs_l) + "% complete");
          }
        }
      }

      // Create a matrix that will hold the x,y, and z coordinates for terminii of the first & second strands
      for( auto itr_ensemble( ENSEMBLE.Begin()), itr_ensemble_end( ENSEMBLE.End()); itr_ensemble != itr_ensemble_end; ++itr_ensemble)
      {
        // get current protein model
        util::ShPtr< assemble::ProteinModel> sp_current_model( *itr_ensemble);
        const assemble::ProteinModel &protein_model( *sp_current_model);

        // get current pdb name before all the loops start
        const util::ShPtr< util::Wrapper< std::string> > &model_name_ptr
        (
          protein_model.GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile)
        );
        const std::string model_name( model_name_ptr.IsDefined() ? model_name_ptr->GetData() : "");

        // instantiate sequences CaCb from pdb for each chain one
        BCL_MessageStd( "computing packing for " + model_name);

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
          const long n_helices( ( *itr_chain)->GetSSEs( biol::GetSSTypes().HELIX).GetSize());
          const long n_strands( ( *itr_chain)->GetSSEs( biol::GetSSTypes().STRAND).GetSize());

          const double helix_exp
          (
            double( n_helices)
            / std::max
              (
                1.0,
                double
                (
                  std::min( n_helices, long( max_helix_helix_contacts + 1)) *
                  std::min( n_strands + 1, long( max_helix_strand_contacts + 1))
                )
              )
          );
          const double strand_exp
          (
            double( n_strands)
            / std::max
              (
                1.0,
                double
                (
                  std::min( n_strands, long( max_strand_strand_contacts + 1)) *
                  std::min( n_helices + 1, long( max_strand_helix_contacts + 1))
                )
              )
          );

          for( long i( 0); i <= n_helices; ++i)
          {
            for( long j( 0); j <= n_strands; ++j)
            {
              if( i < n_helices && i <= max_helix_helix_contacts && j <= max_helix_strand_contacts)
              {
                helix_adjacency_plaus( i, j) += helix_exp;
              }
              if( j < n_strands && i <= max_strand_helix_contacts && j <= max_strand_strand_contacts)
              {
                strand_adjacency_plaus( i, j) += strand_exp;
              }
            }
          }

          counts_contacts_possible( 0) += n_helices * double( std::min( long( 6), std::max( long( n_helices - 1), long( 0))));
          counts_contacts_possible( 3) += n_strands * double( std::min( long( 4), std::max( long( n_strands - 1), long( 0))));
          const long expected_n_helix_strand
          (
            std::max( double( n_helices * std::min( n_strands, long( 4))), double( std::min( n_helices, long( 4)) * n_strands))
          );
          counts_contacts_possible( 1) += expected_n_helix_strand;
          counts_contacts_possible( 2) += expected_n_helix_strand;

          if( structured_sses.GetSize() < 2)
          {
            continue;
          }

          util::SiPtrVector< const biol::Atom> si_atoms;

          for( auto itr_sse( structured_sses.Begin()), itr_sse_end( structured_sses.End()); itr_sse != itr_sse_end; ++itr_sse)
          {
            si_atoms.Append( ( *itr_sse)->GetAtoms( types_of_interest));
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
          for( size_t i( 0), nr_sse( interactions_matrix.GetNumberCols()); i < nr_sse; ++i)
          {
            size_t strand_count( 0), helix_count( 0);
            for( size_t j( 0); j < nr_sse; ++j)
            {
              if( interactions_matrix( i, j) >= m_MinAtomsInContact)
              {
                ++( structured_sses( j)->GetType() == biol::GetSSTypes().STRAND ? strand_count : helix_count);
              }
            }
            linal::Matrix< double> &ref_mat( structured_sses( i)->GetType() == biol::GetSSTypes().STRAND ? strand_adjacency_counts : helix_adjacency_counts);
            if( helix_count < ref_mat.GetNumberRows() && strand_count < ref_mat.GetNumberCols())
            {
              ref_mat( helix_count, strand_count) += 1.0;
            }
          }

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
            size_t b( a + 1);
            const size_t a_type( ( *itr_a)->GetType() == biol::GetSSTypes().HELIX ? 0 : 1);
            for
            (
              util::SiPtrVector< const assemble::SSE>::const_iterator itr_b( itr_a + 1);
              itr_b != itr_end;
              ++itr_b, ++b
            )
            {
              const size_t b_type( ( *itr_b)->GetType() == biol::GetSSTypes().HELIX ? 0 : 2);
              const bool is_in_contact( interactions_matrix( a, b) >= m_MinAtomsInContact);
              const assemble::SSEGeometryPacking::Orientation ori
              (
                assemble::SSEGeometryPacking::OrientationFromSSEs( ( **itr_a), ( **itr_b))
              );
              // for every pair of SSEs
              assemble::SSEGeometryPacking packing_ab( **itr_a, **itr_b, 0.0);
              if( !packing_ab.GetContactType().IsDefined())
              {
                continue;
              }
              const contact::Type reverse_type( contact::GetTypes().Reverse( packing_ab.GetContactType()));
              if( !is_in_contact && itr_b != itr_a + 1)
              {
                continue;
              }

              const double max_loop_length
              (
                ( ( *itr_b)->GetFirstAA()->GetPdbID() - ( *itr_a)->GetLastAA()->GetPdbID()) * 2.56 + 2.11
              );
              const bool could_be_opposite_ori
              (
                ( *itr_b)->GetChainID() != ( *itr_a)->GetChainID()
                ||
                std::min
                (
                  linal::Distance( ( *itr_b)->GetLastAA()->GetCA().GetCoordinates(), ( *itr_a)->GetLastAA()->GetCA().GetCoordinates()),
                  std::min
                  (
                    linal::Distance( ( *itr_b)->GetLastAA()->GetCA().GetCoordinates(), ( *itr_b)->GetFirstAA()->GetCA().GetCoordinates()),
                    linal::Distance( ( *itr_a)->GetLastAA()->GetCA().GetCoordinates(), ( *itr_a)->GetFirstAA()->GetCA().GetCoordinates())
                  ) + 2.34
                ) < max_loop_length
              );
//              if( is_in_contact)
//              {
//                BCL_MessageStd
//                (
//                  std::string( itr_b == itr_a + 1 ? "Adj " : "NonAdj ")
//                  + std::string( could_be_opposite_ori ? " CouldSwitch " : " Constrained ")
//                  + std::string( is_parallel ? " Parallel " : " Antiparallel ")
//                  + util::Format()( max_loop_length) + " "
//                  + util::Format()( linal::Distance( ( *itr_b)->GetLastAA()->GetCA().GetCoordinates(), ( *itr_a)->GetLastAA()->GetCA().GetCoordinates()))
//                  + " " + util::Format()( linal::Distance( ( *itr_b)->GetLastAA()->GetCA().GetCoordinates(), ( *itr_b)->GetFirstAA()->GetCA().GetCoordinates()) + 2.34)
//                  + " " + util::Format()( linal::Distance( ( *itr_a)->GetLastAA()->GetCA().GetCoordinates(), ( *itr_a)->GetFirstAA()->GetCA().GetCoordinates()) + 2.34)
//                  + " " + ( *itr_a)->GetIdentification() + " " + ( *itr_b)->GetIdentification()
//                );
//              }
              const size_t sse_separation( itr_b - itr_a - 1);
              AddPackingType
              (
                packing_ab,
                is_in_contact,
                sse_separation,
                false,
                could_be_opposite_ori,
                packing_info_by_contact_type( packing_ab.GetContactType())
              );
              AddPackingType
              (
                packing_ab,
                is_in_contact,
                sse_separation,
                false,
                could_be_opposite_ori,
                packing_info_by_contact_type( reverse_type)
              );
              if( is_in_contact)
              {
                counts_contacts_actual( b_type | a_type) += 1;
                if( a_type != b_type)
                {
                  counts_contacts_actual( ( a_type ? 2 : 0) + ( b_type ? 1 : 0)) += 1;
                }
              }
              ++contact_type_by_ss_types( a_type + b_type)( packing_ab.GetContactType());
              ++contact_type_by_ss_types( ( a_type ? 2 : 0) + ( b_type ? 1 : 0))( reverse_type);
            }
          }
        }
      }

      // write statistics
      std::ostringstream stream;

      stream << this->GetString() << '\n';

      const size_t n_contact_types( contact_type_by_ss_types( 0).GetSize());
      stream << "Type\t";
      for
      (
        auto itr_ct( contact::GetTypes().Begin()), itr_ct_end( contact::GetTypes().End());
        itr_ct != itr_ct_end;
        ++itr_ct
      )
      {
        stream << itr_ct->GetName() << '\t';
      }

      stream << '\n';
      stream << "LogOdds\t";
      for( size_t ct( 0); ct < n_contact_types; ++ct)
      {
        int config( 0);
        while( config < 4 && !contact_type_by_ss_types_bg( config)( ct))
        {
          ++config;
        }
        const linal::Vector< size_t> &counts_real( contact_type_by_ss_types( config));
        const linal::Vector< size_t> &counts_bg( contact_type_by_ss_types_bg( config));
        const size_t sum_real( counts_real.Sum());
        const size_t sum_bg( counts_bg.Sum());
        const double counts_ratio( double( sum_real) / double( sum_bg));
        stream << -std::log( ( double( counts_real( ct)) / double( counts_bg( ct))) / counts_ratio) << '\t';
      }
      stream << '\n';
      stream << "\nContactType\t";
      for( size_t i( 0); i < s_NumberNaturalCategories; ++i)
      {
        stream << CategoryEnum( Category( i)).GetString() << '\t';
      }
      stream << '\n';
      for( size_t ct( 0); ct < n_contact_types; ++ct)
      {
        stream << contact::Type( ct).GetName() << '\t';
        const linal::Vector< size_t> &contact_counts( packing_info_by_contact_type( ct));
        const double adjacent_contact( contact_counts( e_AdjacentParallel) + contact_counts( e_AdjacentAntiParallel));
        const double one_sse_contact( contact_counts( e_OneSSEApartParallel) + contact_counts( e_OneSSEApartAntiParallel));
        const double two_sse_contact( contact_counts( e_TwoSSEApartParallel) + contact_counts( e_TwoSSEApartAntiParallel));
        const double three_sse_contact( contact_counts( e_ThreeOrMoreSSEApartParallel) + contact_counts( e_ThreeOrMoreSSEApartAntiParallel));
        const double real_total
        (
          contact_counts( e_WeakInteraction)
          + contact_counts( e_ModerateInteraction)
          + contact_counts( e_StrongInteraction)
        );
        const double bg_total
        (
          contact_counts( e_BackgroundWeakInteraction)
          + contact_counts( e_BackgroundModerateInteraction)
          + contact_counts( e_BackgroundStrongInteraction)
        );

        const double bg_parallel_propensity
        (
          double( contact_counts( e_BackgroundParallel)) / bg_total
        );
        const double bg_antiparallel_propensity
        (
          double( contact_counts( e_BackgroundAntiParallel)) / bg_total
        );
        const double real_bg_ratio( real_total / bg_total);
        const double adjacent_total( contact_counts( e_AdjacentInContact) + contact_counts( e_AdjacentNotInContact));
        stream << -std::log( double( contact_counts( e_AdjacentInContact) * 2) / adjacent_total) << '\t'
               << -std::log( double( contact_counts( e_AdjacentNotInContact) * 2) / adjacent_total) << '\t'
               << -std::log( double( contact_counts( e_AdjacentParallel)) / adjacent_contact / bg_parallel_propensity) << '\t'
               << -std::log( double( contact_counts( e_AdjacentAntiParallel)) / adjacent_contact / bg_antiparallel_propensity) << '\t'
               << -std::log( double( contact_counts( e_OneSSEApartParallel) / one_sse_contact / bg_parallel_propensity)) << '\t'
               << -std::log( double( contact_counts( e_OneSSEApartAntiParallel) / one_sse_contact / bg_antiparallel_propensity)) << '\t'
               << -std::log( double( contact_counts( e_TwoSSEApartParallel) / two_sse_contact / bg_parallel_propensity)) << '\t'
               << -std::log( double( contact_counts( e_TwoSSEApartAntiParallel) / two_sse_contact / bg_antiparallel_propensity)) << '\t'
               << -std::log( double( contact_counts( e_ThreeOrMoreSSEApartParallel) / three_sse_contact / bg_parallel_propensity)) << '\t'
               << -std::log( double( contact_counts( e_ThreeOrMoreSSEApartAntiParallel) / three_sse_contact / bg_antiparallel_propensity)) << '\t'
               << -std::log( double( contact_counts( e_WeakInteraction) + real_bg_ratio) / double( contact_counts( e_BackgroundWeakInteraction) + 1) / real_bg_ratio) << '\t'
               << -std::log( double( contact_counts( e_ModerateInteraction) + real_bg_ratio) / double( contact_counts( e_BackgroundModerateInteraction) + 1) / real_bg_ratio) << '\t'
               << -std::log( double( contact_counts( e_StrongInteraction) + real_bg_ratio) / double( contact_counts( e_BackgroundStrongInteraction) + 1) / real_bg_ratio) << '\n';
      }
      stream << '\n';
      stream << "\n\n# The following information is not currently use but is collected to aid in the development of further scores\n"
             << "SS-PairType\tMaxPairingsPossible\tActualNumber\n"
             << "HelixHelix\t" << counts_contacts_possible( 0) << '\t' << counts_contacts_actual( 0) << '\n'
             << "HelixStrand\t" << counts_contacts_possible( 1) << '\t' << counts_contacts_actual( 1) << '\n'
             << "StrandStrand\t" << counts_contacts_possible( 3) << '\t' << counts_contacts_actual( 3) << '\n';

      stream << "\n# Adjacency Counts for strands. Rows - # of helices in contact. Columns - # of strands in contact\n";
      stream << "\n" << strand_adjacency_counts;
      stream << "\nAdjacency Counts for helices. Rows - # of helices in contact. Columns - # of strands in contact\n"
             << "\n" << helix_adjacency_counts;
      auto itr_plaus( strand_adjacency_plaus.Begin());
      for( auto itr( strand_adjacency_counts.Begin()), itr_end( strand_adjacency_counts.End()); itr != itr_end; ++itr, ++itr_plaus)
      {
        *itr /= *itr_plaus;
      }
      itr_plaus = helix_adjacency_plaus.Begin();
      for( auto itr( helix_adjacency_counts.Begin()), itr_end( helix_adjacency_counts.End()); itr != itr_end; ++itr, ++itr_plaus)
      {
        *itr /= *itr_plaus;
      }
      stream << "\n#Adjacency Propensity for strands (# times / # times expected from background w/ same helix-strand distribution) for helices.\n"
             << "\n#Rows - # of helices in contact. Columns - # of strands in contact\n";
      stream << "\n" << strand_adjacency_counts << "\n";
      stream << "\n#Adjacency propensity for helices (# times / # times expected from background w/ same helix-strand distribution) for helices.\n"
              << "\n#Rows - # of helices in contact. Columns - # of strands in contact\n";
      stream << helix_adjacency_counts;

      return stream.str();
    } // end of operator

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinModelPacking::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Computes packing type, distance, interaction weight entropies");

      parameters.AddInitializer
      (
        "distance",
        "max distance between atoms to be considered in contact",
        io::Serialization::GetAgent( &m_InteractionDistance),
        util::Format()( ProteinModelPacking().m_InteractionDistance)
      );
      parameters.AddInitializer
      (
        "min sc atoms",
        "minimum side chains atoms betwen contacting SSEs",
        io::Serialization::GetAgent( &m_MinAtomsInContact),
        util::Format()( ProteinModelPacking().m_MinAtomsInContact)
      );
      return parameters;
    }

    //! @brief read in the SSPair to contact entropy table into a vector with
    //! @param STREAM Input stream to read from
    linal::Vector< double> ProteinModelPacking::ReadSSPairToContactEntropies( std::istream &STREAM) const
    {
      storage::Table< double> table;
      table.ReadFormatted( STREAM);
      BCL_Assert( table.GetNumberRows() == size_t( 1), "Table should have size 4!");
      linal::Vector< double> contact_type_entropies( contact::GetTypes().GetEnumCount());
      const size_t n_cols( table.GetNumberColumns());
      auto itr_col( table.Begin()->Second().GetData().Begin());
      for( size_t col_number( 0); col_number < n_cols; ++col_number, ++itr_col)
      {
        contact_type_entropies( contact::Type( table.GetHeader()( col_number)).GetIndex()) = *itr_col;
      }
      return contact_type_entropies;
    }

    //! @brief read in the entropy table into a vec of vecs. Rows are indexed by contact type, columns by Categories Enum
    //! @param STREAM Input stream to read from
    storage::Vector< linal::Vector< double> > ProteinModelPacking::ReadContactTypeEntropies( std::istream &STREAM) const
    {
      storage::Table< double> table;
      table.ReadFormatted( STREAM);
      storage::Vector< linal::Vector< double> > contact_entropies
      (
        contact::GetTypes().GetEnumCount(),
        linal::Vector< double>( size_t( s_NumberNaturalCategories), double( 0.0))
      );
      const size_t n_cols( table.GetNumberColumns());
      for( auto itr( table.Begin()), itr_end( table.End()); itr != itr_end; ++itr)
      {
        const size_t row_number( contact::Type( itr->First()).GetIndex());
        auto itr_col( itr->Second().GetData().Begin());
        for( size_t col_number( 0); col_number < n_cols; ++col_number, ++itr_col)
        {
          contact_entropies( row_number)( CategoryEnum( table.GetHeader()( col_number))) = *itr_col;
        }
      }
      return contact_entropies;
    }
  } // namespace scorestat
} // namespace bcl

