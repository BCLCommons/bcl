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
#include "score/bcl_score_aa_pair_contact_energy.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_voxel_grid_aa.h"
#include "biol/bcl_biol_aa_base.h"
#include "biol/bcl_biol_atom.h"
#include "command/bcl_command_command_state.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_histogram_2d.h"
#include "math/bcl_math_histogram_3d.h"
#include "score/bcl_score_aa_pair_hi_res_clash.h"
#include "score/bcl_score_energy_distribution.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_logger_interface.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    const util::SiPtr< const util::ObjectInterface> AAPairContactEnergy::s_InteractionInstance
    (
      util::Enumerated< ProteinModel>::AddInstance( new AAPairContactEnergy())
    );

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &AAPairContactEnergy::GetAlias() const
    {
      static const std::string s_name( "AAPairInteractionE");
      return s_name;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a specified histogram file
    AAPairContactEnergy::AAPairContactEnergy() :
      m_HistogramFileName( "aapair_contact_energies.histograms3D"),
      m_Scheme( "aa_pair_interaction"),
      m_Histograms(),
      m_DistanceCutoff( 0.0),
      m_ContactProbabilityCutoff( 0.01),
      m_ConsiderLoops( true),
      m_InterfaceOnly( true)
    {
      // read the histogram file and store the energy functions
      ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new AAPairContactEnergy object that is copied from this one
    AAPairContactEnergy *AAPairContactEnergy::Clone() const
    {
      return new AAPairContactEnergy( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AAPairContactEnergy::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate amino acid pairing potential for given amino acid pair
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @param DISTANCE distance between the amino acid pair
    //! @param TYPE_A, TYPE_B SS-type for AMINO_ACID_A/B
    //! @return amino acid pairing potential for given amino acid pair
    double AAPairContactEnergy::operator()
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const double DISTANCE,
      const biol::SSType &TYPE_A,
      const biol::SSType &TYPE_B
    ) const
    {
      if( AMINO_ACID_A.GetType() > AMINO_ACID_B.GetType())
      {
        return operator()( AMINO_ACID_B, AMINO_ACID_A, DISTANCE, TYPE_B, TYPE_A);
      }

      // calculate the angles distance
      const linal::Vector3D &ca_a( AMINO_ACID_A.GetCA().GetCoordinates());
      const linal::Vector3D &cb_a( AMINO_ACID_A.GetFirstSidechainAtom().GetCoordinates());
      const linal::Vector3D &ca_b( AMINO_ACID_B.GetCA().GetCoordinates());
      const linal::Vector3D &cb_b( AMINO_ACID_B.GetFirstSidechainAtom().GetCoordinates());

      // check that distance is defined
      if
      (
        !util::IsDefined( DISTANCE)
        || !ca_a.IsDefined()
        || !ca_b.IsDefined()
        || !cb_a.IsDefined()
        || !cb_b.IsDefined()
        || !AMINO_ACID_A.GetType()->IsNaturalAminoAcid()
        || !AMINO_ACID_B.GetType()->IsNaturalAminoAcid()
        || TYPE_A->GetIndex() > size_t( 2)
        || TYPE_B->GetIndex() > size_t( 2)
      )
      {
        return 0.0;
      }

      const double angle_b1_a1_a2( linal::ProjAngleCosinus( ca_a, cb_a, cb_b));
      const double angle_b2_a2_a1( linal::ProjAngleCosinus( ca_b, cb_b, cb_a));

      double f
      (
        ( *m_Histograms)( TYPE_A->GetIndex())( TYPE_B->GetIndex())
                        ( AMINO_ACID_A.GetType())( AMINO_ACID_B.GetType())
                        ->Value( DISTANCE, angle_b2_a2_a1, angle_b1_a1_a2)
      );

      // return the score
      return f;
    }

    //! @brief calculate clash score for a protein model
    //! @param MODEL the protein model of interest
    //! @return amino acid pairing potential for given protein
    double AAPairContactEnergy::operator()( const assemble::ProteinModel &MODEL) const
    {
      assemble::VoxelGridAA s_voxel_grid( m_DistanceCutoff);

      s_voxel_grid.SetObjects( MODEL.GetAminoAcids());
      double score( 0.0);

      // hash sses for all amino acid ids
      storage::Vector< storage::Vector< size_t> > sse_id( size_t( 128));
      storage::Vector< storage::Vector< size_t> > sse_type( size_t( 128));
      size_t current_sse_id( 0);
      for
      (
        auto itr_chain( MODEL.GetChains().Begin()), itr_chain_end( MODEL.GetChains().End());
        itr_chain != itr_chain_end;
        ++itr_chain
      )
      {
        storage::Vector< size_t> &chain_vec( sse_id( size_t( ( *itr_chain)->GetChainID())));
        chain_vec.Resize
        (
          size_t
          (
            ( *itr_chain)->GetSequence()->GetLastMember()->GetSeqID()
            + 1
          ),
          util::GetUndefined< size_t>()
        );
        storage::Vector< size_t> &chain_type_vec( sse_type( size_t( ( *itr_chain)->GetChainID())));
        chain_type_vec.Resize
        (
          size_t
          (
            ( *itr_chain)->GetSequence()->GetLastMember()->GetSeqID()
            + 1
          ),
          biol::GetSSTypes().COIL
        );
        util::SiPtrVector< const assemble::SSE> sses( ( *itr_chain)->GetSSEs());
        for( auto itr_sses( sses.Begin()), itr_sses_end( sses.End()); itr_sses != itr_sses_end; ++itr_sses)
        {
          // whether to consider loops
          if( !m_ConsiderLoops)
          {
            // for now, we aren't considering coils, though the data for this exists in the histogram file
            if( !( *itr_sses)->GetType()->IsStructured())
            {
              continue;
            }
          }
          int seq_id_start( ( *itr_sses)->GetFirstAA()->GetSeqID());
          int seq_id_end( ( *itr_sses)->GetLastAA()->GetSeqID());
          for( int pos( seq_id_start); pos <= seq_id_end; ++pos)
          {
            chain_vec( pos) = current_sse_id;
            chain_type_vec( pos) = ( *itr_sses)->GetType();
          }
          ++current_sse_id;
        }
      }
      storage::Vector< storage::Triplet< util::SiPtr< const biol::AABase>, util::SiPtr< const biol::AABase>, double> >
        res( s_voxel_grid.GetNeighbors( m_DistanceCutoff));
      for( auto itr( res.Begin()), itr_end( res.End()); itr != itr_end; ++itr)
      {
        const auto &triplet( *itr);
        const size_t sse1( sse_id( size_t( triplet.First()->GetChainID()))( triplet.First()->GetSeqID()));
        const size_t sse2( sse_id( size_t( triplet.Second()->GetChainID()))( triplet.Second()->GetSeqID()));
        const biol::SSType sse1_type( sse_type( size_t( triplet.First()->GetChainID()))( triplet.First()->GetSeqID()));
        const biol::SSType sse2_type( sse_type( size_t( triplet.Second()->GetChainID()))( triplet.Second()->GetSeqID()));
        if( !util::IsDefined( sse1) || !util::IsDefined( sse2))
        {
          continue;
        }
        if( sse1 == sse2 && util::IsDefined( sse1))
        {
          continue;
        }
        if( biol::SequenceSeparation( *triplet.First(), *triplet.Second()) < 2)
        {
          // ignore adjacent aas
          continue;
        }
        if( m_InterfaceOnly)
        {
          // ignore residue pairs in the same chain
          if( triplet.First()->GetChainID() == triplet.Second()->GetChainID())
          {
            continue;
          }
        }
        score += this->operator()( *triplet.First(), *triplet.Second(), triplet.Third(), sse1_type, sse2_type);
      }
      return score;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param AMINO_ACID_A first amino acid of interest
    //! @param AMINO_ACID_B second amino acid of interest
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    AAPairContactEnergy::WriteDetailedSchemeAndValues
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const biol::SSType &TYPE_A,
      const biol::SSType &TYPE_B,
      std::ostream &OSTREAM
    ) const
    {
      // calculate the CB Distance
      const double cb_distance( biol::FirstSidechainAtomDistance( AMINO_ACID_A, AMINO_ACID_B));

      // if distance is not defined return 0
      if( !util::IsDefined( cb_distance))
      {
        return OSTREAM;
      }
      // calculate the angles distance
      const linal::Vector3D &ca_a( AMINO_ACID_A.GetCA().GetCoordinates());
      const linal::Vector3D &cb_a( AMINO_ACID_A.GetFirstSidechainAtom().GetCoordinates());
      const linal::Vector3D &ca_b( AMINO_ACID_B.GetCA().GetCoordinates());
      const linal::Vector3D &cb_b( AMINO_ACID_B.GetFirstSidechainAtom().GetCoordinates());

      // check that distance is defined
      if
      (
        !ca_a.IsDefined()
        || !ca_b.IsDefined()
        || !cb_a.IsDefined()
        || !cb_b.IsDefined()
        || !AMINO_ACID_A.GetType()->IsNaturalAminoAcid()
        || !AMINO_ACID_B.GetType()->IsNaturalAminoAcid()
      )
      {
        return OSTREAM;
      }

      const double angle_b1_a1_a2( linal::ProjAngleCosinus( ca_a, cb_a, cb_b));
      const double angle_b2_a2_a1( linal::ProjAngleCosinus( ca_b, cb_b, cb_a));

      // write Scheme
      OSTREAM << AMINO_ACID_A.GetSeqID() << '\t'
              << AMINO_ACID_A.GetType()->GetThreeLetterCode() << '\t'
              << AMINO_ACID_B.GetSeqID() << '\t'
              << AMINO_ACID_B.GetType()->GetThreeLetterCode() << '\t'
              << cb_distance << '\t'
              << angle_b2_a2_a1 << '\t'
              << angle_b1_a1_a2 << '\t'
              << operator()( AMINO_ACID_A, AMINO_ACID_B, cb_distance, TYPE_A, TYPE_B) << '\n';

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param MODEL model of interest
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    AAPairContactEnergy::WriteDetailedSchemeAndValues
    (
      const assemble::ProteinModel &MODEL,
      std::ostream &OSTREAM
    ) const
    {
      assemble::VoxelGridAA s_voxel_grid( m_DistanceCutoff);

      s_voxel_grid.SetObjects( MODEL.GetAminoAcids());

      // hash sses for all amino acid ids
      // hash sses for all amino acid ids
      storage::Vector< storage::Vector< size_t> > sse_id( size_t( 128));
      storage::Vector< storage::Vector< size_t> > sse_type( size_t( 128));
      size_t current_sse_id( 0);
      for
      (
        auto itr_chain( MODEL.GetChains().Begin()), itr_chain_end( MODEL.GetChains().End());
        itr_chain != itr_chain_end;
        ++itr_chain
      )
      {
        storage::Vector< size_t> &chain_vec( sse_id( size_t( ( *itr_chain)->GetChainID())));
        chain_vec.Resize
        (
          size_t
          (
            ( *itr_chain)->GetSequence()->GetLastMember()->GetSeqID()
            + 1
          ),
          util::GetUndefined< size_t>()
        );
        storage::Vector< size_t> &chain_type_vec( sse_type( size_t( ( *itr_chain)->GetChainID())));
        chain_type_vec.Resize
        (
          size_t
          (
            ( *itr_chain)->GetSequence()->GetLastMember()->GetSeqID()
            + 1
          ),
          biol::GetSSTypes().COIL
        );
        util::SiPtrVector< const assemble::SSE> sses( ( *itr_chain)->GetSSEs());
        for( auto itr_sses( sses.Begin()), itr_sses_end( sses.End()); itr_sses != itr_sses_end; ++itr_sses)
        {
          // for now, we aren't considering coils, though the data for this exists in the histogram file
          if( !( *itr_sses)->GetType()->IsStructured())
          {
            continue;
          }
          int seq_id_start( ( *itr_sses)->GetFirstAA()->GetSeqID());
          int seq_id_end( ( *itr_sses)->GetLastAA()->GetSeqID());
          for( int pos( seq_id_start); pos <= seq_id_end; ++pos)
          {
            chain_vec( pos) = current_sse_id;
            chain_type_vec( pos) = ( *itr_sses)->GetType();
          }
          ++current_sse_id;
        }
      }
      storage::Vector< storage::Triplet< util::SiPtr< const biol::AABase>, util::SiPtr< const biol::AABase>, double> >
        res( s_voxel_grid.GetNeighbors( m_DistanceCutoff));
      for( auto itr( res.Begin()), itr_end( res.End()); itr != itr_end; ++itr)
      {
        const auto &triplet( *itr);
        const size_t sse1( sse_id( size_t( triplet.First()->GetChainID()))( triplet.First()->GetSeqID()));
        const size_t sse2( sse_id( size_t( triplet.Second()->GetChainID()))( triplet.Second()->GetSeqID()));
        const biol::SSType sse1_type( sse_type( size_t( triplet.First()->GetChainID()))( triplet.First()->GetSeqID()));
        const biol::SSType sse2_type( sse_type( size_t( triplet.Second()->GetChainID()))( triplet.Second()->GetSeqID()));
        if( !util::IsDefined( sse1) || !util::IsDefined( sse2))
        {
          continue;
        }
        if( sse1 == sse2 && util::IsDefined( sse1))
        {
          continue;
        }
        if( biol::SequenceSeparation( *triplet.First(), *triplet.Second()) < 2)
        {
          // ignore adjacent aas
          continue;
        }
        this->WriteDetailedSchemeAndValues( *triplet.First(), *triplet.Second(), sse1_type, sse2_type, OSTREAM);
      }
      return OSTREAM;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAPairContactEnergy::GetSerializer() const
    {
      io::Serializer serializer;

      serializer.SetClassDescription
      (
        "Scores amino acid pair distance and angles for interaction energy (neglecting clash) "
        "Uses probability of interaction (based on distance and angles) computed from the pdb and contact energies from "
        "http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-4-8"
      );
      serializer.AddInitializer
      (
        "contact p cutoff",
        "if two residues are less than this % likely to be interacting at their current distance/orientation, use 0 for "
        "this score rather than the computed propensity",
        io::Serialization::GetAgent( &m_ContactProbabilityCutoff),
        "0.01"
      );
      serializer.AddInitializer
      (
        "consider loops",
        "whether consider loops when computing contact energy and clashes",
        io::Serialization::GetAgent( &m_ConsiderLoops),
        "True"
      );
      serializer.AddInitializer
      (
        "interface only",
        "whether to score interface only",
        io::Serialization::GetAgent( &m_InterfaceOnly),
        "True"
      );

      return serializer;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief read map of amino acid pair energies based on distance from histogram files
    void AAPairContactEnergy::ReadEnergyFunctionMap()
    {
      m_Histograms = util::ToSiPtr( GetHistograms( m_HistogramFileName));

      // static distance cutoff
      // might want to shorten this. Beyond 8A, only charged or aromatic types have a significant
      // chance of interacting, and even then only if directly facing each other. Having it at
      // 12 A costs a lot of extra computation
      m_DistanceCutoff = 12.0;
    }

    //! @brief set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @return ERROR_STREAM stream with which to write errors
    bool AAPairContactEnergy::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
    )
    {
      if( !command::CommandState::IsInStaticInitialization())
      {
        ReadEnergyFunctionMap();
      }
      return true;
    }

    //! @brief Get histograms from a particular file. Caches histograms so they need only be read in once
    const AAPairContactEnergy::SSEPairAAPairHistogramType &AAPairContactEnergy::GetHistograms( const std::string &FILENAME) const
    {
      static storage::Map< std::string, SSEPairAAPairHistogramType> s_histograms;
      SSEPairAAPairHistogramType &histogram_vec( s_histograms[ FILENAME]);
      if( !histogram_vec( 0)( 0).IsEmpty() || command::CommandState::IsInStaticInitialization())
      {
        return histogram_vec;
      }

      const size_t n_sse_types( 3);
      // store map of all aa pairs with their histogram
      for( size_t i( 0); i < n_sse_types; ++i)
      {
        for( size_t j( 0); j < n_sse_types; ++j)
        {
          histogram_vec( i)( j).Reset();
          histogram_vec( i)( j).Resize
          (
            biol::AATypes::s_NumberStandardAATypes,
            util::ShPtrVector< math::Histogram3D>( biol::AATypes::s_NumberStandardAATypes)
          );
        }
      }

      AAPairHiResClash clash;
      auto clash_hist( clash.GetHistograms());

      // read file with all histograms for each pair of aa types
      io::IFStream read;
      io::File::MustOpenIFStream( read, Score::AddHistogramPath( FILENAME));

      size_t outer_ss_type( 0), inner_ss_type( 0);

      // initialize temporary strings to read
      std::string tmp_a, tmp_b;

      // while reading the aatype pair and until reaching the end of the file
      while( read >> tmp_a >> tmp_b && !read.eof())
      {
        // read the string and convert it to first aa type
        biol::AAType a( biol::GetAATypes().AATypeFromOneLetterCode( tmp_a[0]));
        biol::AAType b( biol::GetAATypes().AATypeFromOneLetterCode( tmp_b[0]));

        // read the histogram
        util::ShPtr< math::Histogram3D> histogram( new math::Histogram3D);
        read >> *histogram;
        if( histogram_vec( outer_ss_type)( inner_ss_type)( a)( b).IsDefined())
        {
          ++inner_ss_type;
          if( inner_ss_type == n_sse_types)
          {
            inner_ss_type = 0;
            ++outer_ss_type;
          }
        }
        auto binning( histogram->GetBinningXYZ());
        const size_t n_bins_x( histogram->GetNumberOfBinsX()),
                     n_bins_y( histogram->GetNumberOfBinsY()),
                     n_bins_z( histogram->GetNumberOfBinsZ());
        auto clash_hist_ab( *clash_hist( a)( b));
        for( size_t i( 0); i < n_bins_x; ++i)
        {
          for( size_t j( 0); j < n_bins_y; ++j)
          {
            for( size_t k( 0); k < n_bins_z; ++k)
            {
              if( clash_hist_ab.Interpolate( binning( 0)( i), binning( 1)( j), binning( 2)( k)) < m_ContactProbabilityCutoff)
              {
                histogram->GetChangeableHistogram()( i, j, k) = 0.0;
              }
            }
          }
        }
        histogram_vec( outer_ss_type)( inner_ss_type)( a)( b)
          = histogram_vec( outer_ss_type)( inner_ss_type)( b)( a) = histogram;
      }

      // close the stream
      io::File::CloseClearFStream( read);
      return histogram_vec;
    }

  } // namespace score
} // namespace bcl
