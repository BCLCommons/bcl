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
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_aligner_dp.h"
#include "align/bcl_align_handler_fasta.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "function/bcl_function_binary_sum.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "score/bcl_score_aa_assignment_blast_profile.h"
#include "score/bcl_score_aa_assignment_blosum.h"
#include "score/bcl_score_aa_assignment_identity.h"
#include "score/bcl_score_aa_assignment_pam.h"
#include "score/bcl_score_aa_assignment_phat.h"
#include "score/bcl_score_aa_assignment_property.h"
#include "score/bcl_score_aa_assignment_ss_prediction.h"

// external includes

namespace bcl
{
  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignToMSA
    //! @brief Aligns one target fasta to a given MSA sequence and outputs a new MSA with the target sequence in the
    //!        first positio
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_app_AlignToMSA.cpp @endlink
    //! @author teixeipl
    //! @date June 2, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AlignToMSA :
      public Interface
    {
    private:
    //////////
    // data //
    //////////

      // pair scores
      enum s_PairScores
      {
        IDENTITY,
        PAM100,
        PAM120,
        PAM160,
        PAM250,
        BLOSUM90,
        BLOSUM80,
        BLOSUM62,
        BLOSUM45,
        PHAT85,
        PHAT80,
        PHAT75,
        PHAT70,
        BLAST,
        PSIPRED,
        JUFO,
        SAM,
        TMHMM,
        TMMOD,
        B2TMPRED,
        PROFTMB,
        CONPRED,
        STERICAL_PARAMETER,
        POLARIZABILITY,
        VOLUME,
        HYDROPHOBICITY,
        ISOELECTRIC_POINT,
        TFE_WHITE,
        TFE_ENGELMAN
      };

      static const size_t s_NumberPairScores = 29; //!< number pair scores
      static const storage::Pair< double, double> s_PairScoresToZScoreConversion[]; //!< parameters to convert PairScores to ZScores

      //! path for input files
      util::ShPtr< command::FlagInterface> m_PathFlag;

      //! flag to define the input target fasta filename
      util::ShPtr< command::FlagInterface> m_InputFastaPrefixFlag;

      //! flag to define the input MSA filename
      util::ShPtr< command::FlagInterface> m_InputMSAPrefixFlag;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AlignToMSA();

      //! @brief Clone function
      //! @return pointer to new AlignToMSA
      AlignToMSA *Clone() const
      {
        return new AlignToMSA( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      //! @return initialized command object
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief Main function
      //! @return return value of the application
      int Main() const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief adds the given pair score with the given weight to the scoring function
      //! @param SCORE the pair score to add
      //! @param SCORE_WEIGHT the weight for the pair score
      //! @param SCORE_FCT the scoring function (will be modified)
      //! @return if score was added successfully
      bool AddScore
      (
        const s_PairScores SCORE,
        double SCORE_WEIGHT,
        function::BinarySum< const biol::AABase, const biol::AABase, double> &SCORE_FCT
      ) const;

      //! static instance of this class
      static const ApplicationType AlignToMSA_Instance;
    }; // class AlignToMSA

    //! @brief default constructor
    AlignToMSA::AlignToMSA() :
      m_PathFlag
      (
        new command::FlagStatic
        (
          "path",
          "flag for setting path where input files can be found",
          command::Parameter( "path_param", "path where input files can be found", ".")
        )
      ),
      m_InputFastaPrefixFlag
      (
        new command::FlagStatic
        (
          "fasta", "Flag to input fasta filename",
          command::Parameter
          (
            "fasta_prefix", "\tInput filename", "target_file"
          )
        )
      ),
      m_InputMSAPrefixFlag
      (
        new command::FlagStatic
        (
          "msa", "Flag to input msa filename",
          command::Parameter
          (
            "fasta_prefix", "\tInput filename", "msa_file"
          )
        )
      )
    {
    }

    //! @brief initializes the command object for that executable
    //! @return initialized command object
    util::ShPtr< command::Command> AlignToMSA::InitializeCommand() const
    {
      // initialize a new command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // path for files
      sp_cmd->AddFlag( m_PathFlag);
      sp_cmd->AddFlag( m_InputFastaPrefixFlag);
      sp_cmd->AddFlag( m_InputMSAPrefixFlag);
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief Main function
    //! @return return value of the application
    int AlignToMSA::Main() const
    {
    /////////////////////
    // initializations //
    /////////////////////

      // Initialize with appropriate directory path
      const std::string path_initializer( m_PathFlag->GetFirstParameter()->GetValue() + PATH_SEPARATOR);

      // Initialize directory
      const io::Directory directory( path_initializer);

      const std::string input_fasta( m_InputFastaPrefixFlag->GetFirstParameter()->GetValue() + ".fasta");
      const std::string msa_prefix( m_InputMSAPrefixFlag->GetFirstParameter()->GetValue());
      const std::string input_msa( msa_prefix + ".mfasta");

      BCL_MessageDbg( "Given path:" + path_initializer);
      BCL_MessageDbg( "Given target fasta file with extension:" + input_fasta);
      BCL_MessageDbg( "Given msa file with extension:" + input_msa);

      // Create Score for aligner
      // create assignment score, fold recognition score, see Align-Paper, not user-modifiable
      function::BinarySum< const biol::AABase, const biol::AABase, double> assign_score;
      AddScore( PAM250, 2.5, assign_score);
      AddScore( BLOSUM45, 2.5, assign_score);

//      AddScore( BLAST, 40.0, assign_score);   // SAM method is not supported by BCL anymore; redistribute the weight of
//      AddScore( PSIPRED, 40, assign_score); // 0.07 for SAM equally over the other ss-prediction scores
//      AddScore( JUFO, 13.0, assign_score);    // BLAST and JUFO: 0.18+(0.07/3); PSIPRED: 0.2+(0.07/3)
//      AddScore( STERICAL_PARAMETER, 0.034, assign_score); // the total weight of 0.24 for all chemical properties
//      AddScore( POLARIZABILITY, 0.034, assign_score);     // is equally distributed over them
//      AddScore( VOLUME, 0.034, assign_score);
//      AddScore( HYDROPHOBICITY, 0.034, assign_score);
//      AddScore( ISOELECTRIC_POINT, 0.034, assign_score);
//      AddScore( TFE_WHITE, 0.034, assign_score);
//      AddScore( TFE_ENGELMAN, 0.034, assign_score);

      // Based on Liz's Weight tables -0.672 -0.907 -0.635  0.00

      // fold recognition scores ENCLOSED_SINGLE_GAP, ENCLOSED_MULTIPLE_GAP, BOUNDARY_SINGLE_GAP, BOUNDARY_MULTIPLE_GAP
      score::AssignmentWithGap< biol::AABase> assign_gap_score
      (
        util::CloneToShPtr( assign_score),
        -0.672, // ENCLOSED_SINGLE_GAP
        -0.907, // ENCLOSED_MULTIPLE_GAP
        -0.635, // BOUNDARY_SINGLE_GAP
        0.0 // BOUNDARY_MULTIPLE_GAP
      );

      // Add pdb sequence as first member of alignment since it is the target protein
      align::AlignerDP< biol::AABase> aligner;
      aligner.SetScoringFunction( assign_gap_score);

      // Open file, create object, and read in
      io::IFStream read;
      io::File::MustOpenIFStream( read, path_initializer + input_fasta);

      BCL_MessageDbg( "Reading fasta sequence " + input_fasta);

      // Read in fasta sequence
      util::ShPtr< biol::AASequence> sp_target_seq
      (
        new biol::AASequence( biol::AASequenceFactory::BuildSequenceFromFASTA( read))
      );
      util::ShPtr< align::AlignmentInterface< biol::AABase> > sp_seq_to_align( new align::AlignmentLeaf< biol::AABase>( sp_target_seq));

      // DEBUG to make sure sequence is read in correctly!
      sp_target_seq->WriteFasta( util::GetLogger());
      BCL_MessageDbg( "Reading fasta sequence " + input_fasta + " done");

      BCL_MessageDbg( "Reading alignment portion " + input_msa);
      std::string alignment_file( directory.AppendFilename( input_msa));
      util::ShPtr< align::AlignmentNode< biol::AABase> > alignment_node_ptr;

      // Create necessary PIR handler
      align::HandlerFasta< biol::AABase> mfasta_handler;
      // Check for mfasta file first
      if( io::DirectoryEntry( alignment_file).DoesExist())
      {
        // Create necessary stream and open it
        io::IFStream read_stream;
        io::File::MustOpenIFStream( read_stream, alignment_file);
        alignment_node_ptr = mfasta_handler.ReadAlignment( read_stream, biol::AASequence());
      }

      // Read in MSA
      BCL_MessageDbg( "Reading alignment portion " + input_msa + " done");

      util::ShPtr< align::AlignmentInterface< biol::AABase> > sp_alignment( alignment_node_ptr);

      BCL_MessageDbg( "Generate Alignment");

      // Align target fasta to MSA
      storage::Pair< align::AlignmentNode< biol::AABase>, double> result_pair( aligner.AlignPair( sp_seq_to_align, sp_alignment));
      const align::AlignmentNode< biol::AABase> target_seq_and_alignment
      (
        result_pair.First()
      );
      BCL_MessageDbg( "score test is: " + util::Format()( result_pair.Second()));
      BCL_MessageDbg( "Alignment size afterwards is: "+ util::Format()( target_seq_and_alignment.GetSequences().FirstElement()->GetSize()));
      BCL_MessageDbg( "Generate Alignment done");

    ////////////
    // output //
    ////////////

      // Open output file
      const std::string output_filename( directory.AppendFilename( msa_prefix + ".aligned_mfasta"));
      // Write output new alignment

      // OFstream for writing into bcl file
      io::OFStream write;

      // Write out contact cases
      io::File::MustOpenOFStream( write, output_filename);

      mfasta_handler.WriteAlignment( write, target_seq_and_alignment);
      io::File::CloseClearFStream( write);

      // Return successful int value
      return 0;
    } // Main

    //! @brief adds the given pair score with the given weight to the scoring function
    //! @param SCORE the pair score to add
    //! @param SCORE_WEIGHT the weight for the pair score
    //! @param SCORE_FCT the scoring function (will be modified)
    //! @return if score was added successfully
    bool AlignToMSA::AddScore
    (
      const s_PairScores SCORE,
      double SCORE_WEIGHT,
      function::BinarySum< const biol::AABase, const biol::AABase, double> &SCORE_FCT
    ) const
    {
      double offset( s_PairScoresToZScoreConversion[ SCORE].First());
      double stddev( s_PairScoresToZScoreConversion[ SCORE].Second());

      // correct weight to be a Z-score
      double Zscore_weight( SCORE_WEIGHT / stddev);

      // add different scoring functions
      switch( SCORE)
      {
        case IDENTITY:           SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentIdentity() - offset)); break;
        case PAM250:             SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_250) - offset)); break;
        case PAM160:             SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_160) - offset)); break;
        case PAM120:             SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_120) - offset)); break;
        case PAM100:             SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentPAM( score::AAAssignmentPAM::e_PAM_100) - offset)); break;
        case BLOSUM90:           SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_90) - offset)); break;
        case BLOSUM80:           SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_80) - offset)); break;
        case BLOSUM62:           SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_62) - offset)); break;
        case BLOSUM45:           SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentBLOSUM( score::AAAssignmentBLOSUM::e_BLOSUM_45) - offset)); break;
        case PHAT85:             SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentPHAT( score::AAAssignmentPHAT::e_PHAT_85) - offset)); break;
        case PHAT80:             SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentPHAT( score::AAAssignmentPHAT::e_PHAT_80) - offset)); break;
        case PHAT75:             SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentPHAT( score::AAAssignmentPHAT::e_PHAT_75) - offset)); break;
        case PHAT70:             SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentPHAT( score::AAAssignmentPHAT::e_PHAT_70) - offset)); break;
        case BLAST:              SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentBlastProfile() - offset)); break;
        case PSIPRED:            SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentSSPrediction( sspred::GetMethods().e_PSIPRED) - offset)); break;
        case JUFO:               SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentSSPrediction( sspred::GetMethods().e_JUFO) - offset)); break;
        case SAM:                SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentSSPrediction( sspred::GetMethods().e_SAM) - offset)); break;
        case TMHMM:              SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentSSPrediction( sspred::GetMethods().e_TMHMM) - offset)); break;
        case TMMOD:              SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentSSPrediction( sspred::GetMethods().e_TMMOD) - offset)); break;
        case B2TMPRED:           SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentSSPrediction( sspred::GetMethods().e_B2TMPRED) - offset)); break;
        case PROFTMB:            SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentSSPrediction( sspred::GetMethods().e_PROFTMB) - offset)); break;
        case CONPRED:            SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentSSPrediction( sspred::GetMethods().e_CONPRED) - offset)); break;
        case STERICAL_PARAMETER: SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentProperty( biol::AATypeData::e_StericalParameter) - offset)); break;
        case POLARIZABILITY:     SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentProperty( biol::AATypeData::e_Polarizability) - offset)); break;
        case VOLUME:             SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentProperty( biol::AATypeData::e_Volume) - offset)); break;
        case HYDROPHOBICITY:     SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentProperty( biol::AATypeData::e_Hydrophobicity) - offset)); break;
        case ISOELECTRIC_POINT:  SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentProperty( biol::AATypeData::e_IsoelectricPoint) - offset)); break;
        case TFE_WHITE:          SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentProperty( biol::AATypeData::e_TransferFreeEnergyWhimleyWhite) - offset)); break;
        case TFE_ENGELMAN:       SCORE_FCT += ( Zscore_weight * ( score::AAAssignmentProperty( biol::AATypeData::e_TransferFreeEnergyEngelmanSeitzGoldman) - offset)); break;
        default: return false;
      }

      return true;
    }

    //! parameters to convert PairScores to ZScores
    const storage::Pair< double, double> AlignToMSA::s_PairScoresToZScoreConversion[ s_NumberPairScores] =
    {
      storage::Pair< double, double>(   0.0000,   1.0000), // -identity
      storage::Pair< double, double>(  -0.1951,   0.3109), // -pam100
      storage::Pair< double, double>(  -0.1581,   0.2803), // -pam120
      storage::Pair< double, double>(  -0.1108,   0.2340), // -pam160
      storage::Pair< double, double>(  -0.0824,   0.2498), // -pam250
      storage::Pair< double, double>(  -0.1698,   0.2586), // -blosum90
      storage::Pair< double, double>(  -0.2111,   0.3580), // -blosum80
      storage::Pair< double, double>(  -0.0967,   0.2088), // -blosum62
      storage::Pair< double, double>(  -0.0821,   0.2273), // -blosum45
      storage::Pair< double, double>(  -0.1431,   0.3114), // -phat85
      storage::Pair< double, double>(  -0.1676,   0.4248), // -phat80
      storage::Pair< double, double>(  -0.1469,   0.3990), // -phat75
      storage::Pair< double, double>(  -0.0895,   0.3552), // -phat70
      storage::Pair< double, double>(  -0.0072,   0.0881), // -blast
      storage::Pair< double, double>(  -0.1431,   0.4728), // -psipred
      storage::Pair< double, double>(  -0.0388,   0.2451), // -jufo
      storage::Pair< double, double>(  -0.0056,   0.2076), // -sam
      storage::Pair< double, double>(   0.0000,   0.5000), // -tmhmm
      storage::Pair< double, double>(   0.0000,   0.5000), // -tmmod
      storage::Pair< double, double>(   0.0000,   0.5000), // -b2tmpred
      storage::Pair< double, double>(   0.0000,   0.5000), // -proftmb
      storage::Pair< double, double>(   0.0000,   0.5000), // -conpred
      storage::Pair< double, double>(  -1.1514,   0.8981), // -steric
      storage::Pair< double, double>(  -0.1061,   0.0814), // -polarizability
      storage::Pair< double, double>(  -1.9938,   1.5660), // -volume
      storage::Pair< double, double>(  -1.0737,   0.7871), // -hydrophobicity
      storage::Pair< double, double>(  -1.6180,   1.8058), // -isoelectric
      storage::Pair< double, double>(  -1.8252,   1.4175), // -tfe_white
      storage::Pair< double, double>(  -5.3870,   4.7486)  // -tfe_engelman
    };

    const ApplicationType AlignToMSA::AlignToMSA_Instance
    (
      GetAppGroups().AddAppToGroup( new AlignToMSA(), GetAppGroups().e_InternalBiol)
    );
  } // namespace app
} // namespace bcl
