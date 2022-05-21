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
#include "internal/biology/bcl_app_msa2pssm.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_with_cache_storage_file.h"
#include "biol/bcl_biol_aa.h"
#include "biol/bcl_biol_aa_data.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_blast_profile.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector.h"

namespace bcl
{
  namespace app
  {
  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief default constructor
    MSA2PSSM::MSA2PSSM() :
      m_AlignmentFile( new command::Parameter( "aln", "CLUSTAL-formatted alignment file")),
      m_OutputFile( new command::Parameter( "output", "where the PSSM should be written to")),
      m_PseudocountFlag
      (
        new command::FlagStatic
        (
          "pseudocount",
          "pseudocount to add for unrepresented columns",
          command::Parameter( "pseudocount", "", command::ParameterCheckRanged< double>( 0.0, 50.0), "1.0")
        )
      ),
      m_UseMSABackgroundFlag
      (
        new command::FlagStatic
        (
          "msa_background",
          "Use aa propensities from the MSA itself + a constant pseudocount of 1. Incompatible with -aa_background"
        )
      ),
      m_FractionalPseudocountsFlag
      (
        new command::FlagStatic
        (
          "fractional_pseudocount",
          "Fraction of signal that should be pseudocount. May be better or useful in conjunction with pseudocounts if "
          "the input alignments contain regions of widely varying alignment density. Final pseudocount is computed as: "
          "alignment_count * ( 1 + fractional_pseudocount) + base_pseudocount + unique_aas_at_position * efractional_pseudocount",
          command::Parameter( "fraction", "", command::ParameterCheckRanged< double>( 0.0, 1.0), "0.0")
        )
      ),
      m_EffFractionalPseudocountsFlag
      (
        new command::FlagStatic
        (
          "efractional_pseudocount",
          "effective fraction of signal that should be pseudocount. see equation in -fractional_pseudocount description."
          " This flag controls the degree of mixing due to conservation at a position",
          command::Parameter( "fraction", "", command::ParameterCheckRanged< double>( 0.0, 10.0), "0.0")
        )
      ),
      m_LogOnlyFlag
      (
        new command::FlagStatic
        (
          "log_only",
          "If set, all pssm positions are simply log(pseudocount+real_count), without any row-based normalization"
        )
      ),
      m_Pseudocount( 1.0),
      m_FractionalPseudocount( 0.0),
      m_EffFractionalPseudocount( 0.0),
      m_NumberSeqs( 0)
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MSA2PSSM::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string MSA2PSSM::GetDescription() const
    {
      return "Converts a multiple sequence alignment into a pssm";
    }

  ////////////////
  // operations //
  ////////////////

    util::ShPtr< command::Command> MSA2PSSM::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add member parameters and flags
      sp_cmd->AddParameter( m_AlignmentFile);
      sp_cmd->AddParameter( m_OutputFile);
      sp_cmd->AddFlag( m_PseudocountFlag);
      sp_cmd->AddFlag( m_UseMSABackgroundFlag);
      sp_cmd->AddFlag( m_FractionalPseudocountsFlag);
      sp_cmd->AddFlag( m_EffFractionalPseudocountsFlag);
      sp_cmd->AddFlag( m_LogOnlyFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! Main
    int MSA2PSSM::Main() const
    {
      BCL_Assert
      (
        !m_LogOnlyFlag->GetFlag() || !m_UseMSABackgroundFlag->GetFlag(),
        "Cannot use both log-only normalization with msa background!"
      );
      m_Pseudocount = m_PseudocountFlag->GetFirstParameter()->GetNumericalValue< double>();
      m_FractionalPseudocount = m_FractionalPseudocountsFlag->GetFirstParameter()->GetNumericalValue< double>();
      m_EffFractionalPseudocount = m_EffFractionalPseudocountsFlag->GetFirstParameter()->GetNumericalValue< double>();

      io::IFStream input;
      io::File::MustOpenIFStream( input, m_AlignmentFile->GetValue());

      std::string line;
      BCL_Assert( std::getline( input, line), "Alignment file was empty");
      BCL_Assert( util::StartsWith( line, "CLUSTAL"), "Only CLUSTAL Alignments are supported at this time");
      // clustal formatted alignment file skip the next two lines
      BCL_Assert( std::getline( input, line) && std::getline( input, line), "Bad clustal file");

      // slurp in all the lines
      const storage::Vector< std::string> lines( util::StringLineListFromIStream( input));

      // vector to hold all the sequences
      storage::Vector< std::string> sequences;

      // create a shptr vector to store the sequence
      util::ShPtrVector< biol::AABase> sp_aas;

      // read the alignment file.
      for( size_t pos( 0), block_start( 0), sz( lines.GetSize()); pos < sz; ++pos)
      {
        const std::string &line( lines( pos));
        if( line.empty())
        {
          block_start = pos + 1;
          continue;
        }
        else
        {
          const storage::Vector< std::string> seqname_seq( util::SplitString( line));
          BCL_Assert( seqname_seq.GetSize() == size_t( 2), "Invalid clustal file");
          // first block, append to sequences
          if( pos - block_start >= sequences.GetSize())
          {
            sequences.PushBack( seqname_seq( 1));
          }
          else
          {
            // past first block, append to correct sequence in first block
            sequences( pos - block_start).append( seqname_seq( 1));
          }
        }
      }

      // create a mask containing 0 for each ungapped position, 1 for gapped position
      storage::Vector< size_t> gap( sequences( 0).size(), size_t( 0));

      // create the AA-Sequence
      for( size_t seq_pos( 0), seq_id( 0), seq_len( sequences( 0).size()); seq_pos < seq_len; ++seq_pos)
      {
        const char blast_aa_ch( sequences( 0)[ seq_pos]);
        if( blast_aa_ch == '-')
        {
          gap( seq_pos) = 1;
          continue;
        }
        ++seq_id;

        biol::AAType blast_aa_type( biol::GetAATypes().AATypeFromOneLetterCode( blast_aa_ch));
        sp_aas.PushBack
        (
          util::ShPtr< biol::AABase>
          (
            new biol::AA( util::ShPtr< biol::AAData>( new biol::AAData( blast_aa_type, seq_id, seq_id)))
          )
        );
      }

      // count the number of each AA-type in all the sequences
      storage::Vector< size_t> msa_aa_counts( 26, size_t( 1));
      const size_t ungapped_seq_length( gap.GetSize());
      const size_t seq_length( sp_aas.GetSize());

      // count up the number of times each type was seen
      storage::Vector< std::string> seq_transposed( seq_length);
      m_NumberSeqs = sequences.GetSize();
      for( size_t i( 0); i < m_NumberSeqs; ++i)
      {
        for( size_t j( 0), gapped( 0); j < ungapped_seq_length; ++j)
        {
          const char c( std::toupper( sequences( i)[ j]));
          if( gap( j))
          {
            continue;
          }
          if( c != '-')
          {
            BCL_Assert( c >= 'A' && c <= 'Z', "bad character in multiple sequence alignment: " + util::Format()( c));
            seq_transposed( gapped) += c;
            ++msa_aa_counts( int( c - 'A'));
          }
          ++gapped;
        }
      }

      linal::Vector< double> msa_background_p( biol::AATypes::s_NumberStandardAATypes, double( 0.0));
      for
      (
        biol::AATypes::const_iterator
          itr_aa( biol::GetAATypes().Begin()),
          itr_aa_end( biol::GetAATypes().Begin() + biol::AATypes::s_NumberStandardAATypes);
        itr_aa != itr_aa_end;
        ++itr_aa
      )
      {
        msa_background_p( itr_aa->GetIndex()) = msa_aa_counts( int( ( *itr_aa)->GetOneLetterCode() - 'A'));
      }
      // normalize
      msa_background_p.SetToSum( 1.0);

      // set blast profiles
      for( size_t seq_pos( 0), seq_len( sp_aas.GetSize()); seq_pos < seq_len; ++seq_pos)
      {
        sp_aas( seq_pos)->SetBlastProfile( ComputeBlastProfile( seq_transposed( seq_pos), msa_background_p));
      }

      // write blast profile
      io::OFStream output;
      io::File::MustOpenOFStream( output, m_OutputFile->GetValue());
      biol::BlastProfileHandler::WriteProfileForAASequence( output, biol::AASequence( sp_aas));
      io::File::CloseClearFStream( output);
      // end
      return 0;
    }

    //! @brief compute a blast profile
    //! @param ALIGNED AA-types that aligned at this position (should contain only upper case letters, no gaps)
    //! @param BACKGROUND_AA_FREQ frequency of AAs in the MSA
    //! @return BlastProfile
    biol::BlastProfile MSA2PSSM::ComputeBlastProfile
    (
      const std::string &ALIGNED,
      const linal::Vector< double> &BACKGROUND_AA_FREQ
    ) const
    {
      BCL_Assert( !ALIGNED.empty(), "Empty set of aligned residues; there should at least be the native");

      storage::Vector< size_t> aa_counts( 26, size_t( 0));
      const size_t n_seq( ALIGNED.size());
      // count up the number of times each type was seen
      for( size_t i( 1); i < n_seq; ++i)
      {
        ++aa_counts( int( ALIGNED[ i] - 'A'));
      }

      linal::Vector< double> probs( size_t( biol::AATypes::s_NumberStandardAATypes), 0.0);
      size_t n_unique_types( 0);
      for
      (
        biol::AATypes::const_iterator
          itr_aa( biol::GetAATypes().Begin()),
          itr_aa_end( biol::GetAATypes().Begin() + biol::AATypes::s_NumberStandardAATypes);
        itr_aa != itr_aa_end;
        ++itr_aa
      )
      {
        probs( itr_aa->GetIndex()) = aa_counts( int( ( *itr_aa)->GetOneLetterCode() - 'A'));
        // only consider residues seen in at least 1/100 of the seqs in the alignment in the unique count
        if( probs( itr_aa->GetIndex()) > double( n_seq) / 100.0)
        {
          ++n_unique_types;
        }
      }

      // get a vector reference to the associated row in the blosum matrix
      const char native_type_ch( ALIGNED[ 0]);
      const biol::AAType native_type( biol::GetAATypes().AATypeFromOneLetterCode( native_type_ch));

      // BLOSUM62 frequencies, derived from the original data used to construct the BLOSUM62 log-odds matrix for
      // improved accuracy
      static const double s_BLOSUM62_f[ biol::AATypes::s_NumberStandardAATypes + 1][ biol::AATypes::s_NumberStandardAATypes] =
      {
        { 2.9e-1, 3.1e-2, 2.6e-2, 3.0e-2, 2.2e-2, 2.6e-2, 4.0e-2, 7.8e-2, 1.5e-2, 4.3e-2, 5.9e-2, 4.5e-2, 1.8e-2, 2.2e-2, 3.0e-2, 8.5e-2, 5.0e-2, 5.4e-3, 1.8e-2, 6.9e-2},
        { 3.1e-2, 2.4e-1, 2.7e-2, 2.2e-2, 5.4e-3, 3.4e-2, 3.6e-2, 2.3e-2, 1.6e-2, 1.6e-2, 3.2e-2, 8.4e-2, 1.1e-2, 1.2e-2, 1.3e-2, 3.1e-2, 2.4e-2, 4.0e-3, 1.2e-2, 2.2e-2},
        { 2.6e-2, 2.7e-2, 1.9e-1, 5.0e-2, 5.4e-3, 2.0e-2, 3.0e-2, 3.9e-2, 1.9e-2, 1.3e-2, 1.9e-2, 3.2e-2, 6.7e-3, 1.1e-2, 1.2e-2, 4.2e-2, 3.0e-2, 2.7e-3, 9.4e-3, 1.6e-2},
        { 3.0e-2, 2.2e-2, 5.0e-2, 2.9e-1, 5.4e-3, 2.2e-2, 6.6e-2, 3.4e-2, 1.3e-2, 1.6e-2, 2.0e-2, 3.2e-2, 6.7e-3, 1.1e-2, 1.6e-2, 3.8e-2, 2.6e-2, 2.7e-3, 8.1e-3, 1.8e-2},
        { 2.2e-2, 5.4e-3, 5.4e-3, 5.4e-3, 1.6e-1, 4.0e-3, 5.4e-3, 1.1e-2, 2.7e-3, 1.5e-2, 2.2e-2, 6.7e-3, 5.4e-3, 6.7e-3, 5.4e-3, 1.3e-2, 1.2e-2, 1.3e-3, 4.0e-3, 1.9e-2},
        { 2.6e-2, 3.4e-2, 2.0e-2, 2.2e-2, 4.0e-3, 9.9e-2, 4.7e-2, 1.9e-2, 1.3e-2, 1.2e-2, 2.2e-2, 4.2e-2, 9.4e-3, 6.7e-3, 1.1e-2, 2.6e-2, 1.9e-2, 2.7e-3, 9.4e-3, 1.6e-2},
        { 4.0e-2, 3.6e-2, 3.0e-2, 6.6e-2, 5.4e-3, 4.7e-2, 2.2e-1, 2.6e-2, 1.9e-2, 1.6e-2, 2.7e-2, 5.5e-2, 9.4e-3, 1.2e-2, 1.9e-2, 4.0e-2, 2.7e-2, 4.0e-3, 1.2e-2, 2.3e-2},
        { 7.8e-2, 2.3e-2, 3.9e-2, 3.4e-2, 1.1e-2, 1.9e-2, 2.6e-2, 5.1e-1, 1.3e-2, 1.9e-2, 2.8e-2, 3.4e-2, 9.4e-3, 1.6e-2, 1.9e-2, 5.1e-2, 3.0e-2, 5.4e-3, 1.1e-2, 2.4e-2},
        { 1.5e-2, 1.6e-2, 1.9e-2, 1.3e-2, 2.7e-3, 1.3e-2, 1.9e-2, 1.3e-2, 1.3e-1, 8.1e-3, 1.3e-2, 1.6e-2, 5.4e-3, 1.1e-2, 6.7e-3, 1.5e-2, 9.4e-3, 2.7e-3, 2.0e-2, 8.1e-3},
        { 4.3e-2, 1.6e-2, 1.3e-2, 1.6e-2, 1.5e-2, 1.2e-2, 1.6e-2, 1.9e-2, 8.1e-3, 2.5e-1, 1.5e-1, 2.2e-2, 3.4e-2, 4.0e-2, 1.3e-2, 2.3e-2, 3.6e-2, 5.4e-3, 1.9e-2, 1.6e-1},
        { 5.9e-2, 3.2e-2, 1.9e-2, 2.0e-2, 2.2e-2, 2.2e-2, 2.7e-2, 2.8e-2, 1.3e-2, 1.5e-1, 5.0e-1, 3.4e-2, 6.6e-2, 7.3e-2, 1.9e-2, 3.2e-2, 4.5e-2, 9.4e-3, 3.0e-2, 1.3e-1},
        { 4.5e-2, 8.4e-2, 3.2e-2, 3.2e-2, 6.7e-3, 4.2e-2, 5.5e-2, 3.4e-2, 1.6e-2, 2.2e-2, 3.4e-2, 2.2e-1, 1.2e-2, 1.2e-2, 2.2e-2, 4.2e-2, 3.1e-2, 4.0e-3, 1.3e-2, 2.6e-2},
        { 1.8e-2, 1.1e-2, 6.7e-3, 6.7e-3, 5.4e-3, 9.4e-3, 9.4e-3, 9.4e-3, 5.4e-3, 3.4e-2, 6.6e-2, 1.2e-2, 5.4e-2, 1.6e-2, 5.4e-3, 1.2e-2, 1.3e-2, 2.7e-3, 8.1e-3, 3.1e-2},
        { 2.2e-2, 1.2e-2, 1.1e-2, 1.1e-2, 6.7e-3, 6.7e-3, 1.2e-2, 1.6e-2, 1.1e-2, 4.0e-2, 7.3e-2, 1.2e-2, 1.6e-2, 2.5e-1, 6.7e-3, 1.6e-2, 1.6e-2, 1.1e-2, 5.7e-2, 3.5e-2},
        { 3.0e-2, 1.3e-2, 1.2e-2, 1.6e-2, 5.4e-3, 1.1e-2, 1.9e-2, 1.9e-2, 6.7e-3, 1.3e-2, 1.9e-2, 2.2e-2, 5.4e-3, 6.7e-3, 2.6e-1, 2.3e-2, 1.9e-2, 1.3e-3, 6.7e-3, 1.6e-2},
        { 8.5e-2, 3.1e-2, 4.2e-2, 3.8e-2, 1.3e-2, 2.6e-2, 4.0e-2, 5.1e-2, 1.5e-2, 2.3e-2, 3.2e-2, 4.2e-2, 1.2e-2, 1.6e-2, 2.3e-2, 1.7e-1, 6.3e-2, 4.0e-3, 1.3e-2, 3.2e-2},
        { 5.0e-2, 2.4e-2, 3.0e-2, 2.6e-2, 1.2e-2, 1.9e-2, 2.7e-2, 3.0e-2, 9.4e-3, 3.6e-2, 4.5e-2, 3.1e-2, 1.3e-2, 1.6e-2, 1.9e-2, 6.3e-2, 1.7e-1, 4.0e-3, 1.2e-2, 4.9e-2},
        { 5.4e-3, 4.0e-3, 2.7e-3, 2.7e-3, 1.3e-3, 2.7e-3, 4.0e-3, 5.4e-3, 2.7e-3, 5.4e-3, 9.4e-3, 4.0e-3, 2.7e-3, 1.1e-2, 1.3e-3, 4.0e-3, 4.0e-3, 8.8e-2, 1.2e-2, 5.4e-3},
        { 1.8e-2, 1.2e-2, 9.4e-3, 8.1e-3, 4.0e-3, 9.4e-3, 1.2e-2, 1.1e-2, 2.0e-2, 1.9e-2, 3.0e-2, 1.3e-2, 8.1e-3, 5.7e-2, 6.7e-3, 1.3e-2, 1.2e-2, 1.2e-2, 1.4e-1, 2.0e-2},
        { 6.9e-2, 2.2e-2, 1.6e-2, 1.8e-2, 1.9e-2, 1.6e-2, 2.3e-2, 2.4e-2, 8.1e-3, 1.6e-1, 1.3e-1, 2.6e-2, 3.1e-2, 3.5e-2, 1.6e-2, 3.2e-2, 4.9e-2, 5.4e-3, 2.0e-2, 2.6e-1},
        { 8.3e-2, 5.5e-2, 4.1e-2, 5.5e-2, 1.4e-2, 3.9e-2, 6.8e-2, 7.1e-2, 2.3e-2, 6.0e-2, 9.7e-2, 5.8e-2, 2.4e-2, 3.9e-2, 4.7e-2, 6.6e-2, 5.3e-2, 1.1e-2, 2.9e-2, 6.9e-2}
      };

      // blosum62 matrix without any log/entropic scaling
      linal::VectorConstReference< double> blosum_row
      (
        size_t( biol::AATypes::s_NumberStandardAATypes),
        s_BLOSUM62_f[ std::min( native_type.GetIndex(), size_t( biol::AATypes::s_NumberStandardAATypes))]
      );

      linal::Vector< double> pseudocounts( size_t( biol::AATypes::s_NumberStandardAATypes), blosum_row.Begin());
      const double pseudocount_weight( m_Pseudocount + m_FractionalPseudocount * n_seq + m_EffFractionalPseudocount * n_unique_types);
      if( m_UseMSABackgroundFlag->GetFlag())
      {
        pseudocounts = BACKGROUND_AA_FREQ;
      }
      else if( m_LogOnlyFlag->GetFlag())
      {
        pseudocounts = double( 0.05);
      }
      pseudocounts *= pseudocount_weight;
      probs += pseudocounts;

      const double n_seqm1( n_seq - 1);
      const double alignment_weight( n_seqm1 / ( n_seqm1 + pseudocount_weight));
      if( !m_LogOnlyFlag->GetFlag())
      {
        probs /= probs.Sum();

        for( size_t i( 0), n_aa( biol::AATypes::s_NumberStandardAATypes); i < n_aa; ++i)
        {
          probs( i) /= s_BLOSUM62_f[ 20][ i]; // background frequency
        }
      }

      // obtain a value from the BLOSUM matrix (this is what psiblast does as well)

      // log-odds, base doesn't matter so long as we're consistent
      double( *log_ptr)( double) = std::log;
      std::transform( probs.Begin(), probs.End(), probs.Begin(), log_ptr);

      return biol::BlastProfile( probs, alignment_weight, 0.0);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    const ApplicationType MSA2PSSM::MSA2PSSM_Instance
    (
      GetAppGroups().AddAppToGroup( new MSA2PSSM(), GetAppGroups().e_BioInfo)
    );

  } // namespace app
} // namespace bcl
