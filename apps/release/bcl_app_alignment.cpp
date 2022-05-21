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
#include "bcl_app_alignment.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_handler_classes.h"
#include "align/bcl_align_multiple_aligner_classes.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_parameter_check_extension.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "function/bcl_function_binary_adapter.h"
#include "function/bcl_function_binary_sum.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_z_score.h"
#include "score/bcl_score_alignment_assignment.h"
#include "sspred/bcl_sspred_jufo.h"
#include "sspred/bcl_sspred_method_handler.h"

namespace bcl
{
  namespace app
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Alignment::Alignment() :
      m_FastaListFlag
      (
        new command::FlagDynamic
        (
          "fastas",
          "this is a list of fasta sequence files to be aligned",
          command::Parameter( "fasta file", "any fasta file", command::ParameterCheckExtension( ".fasta")),
          0
        )
      ),
      m_OutputprefixFlag
      (
        new command::FlagDynamic
        (
          "outputprefix",
          "prefix for filenames for the written alignments - if not given, output is written to console",
          command::Parameter( "filename", "filename to be used for output \"{filename}.{outputformat}\""),
          0,
          1
        )
      ),
      m_FlagAligner
      (
        new command::FlagStatic
        (
          "aligner",
          "select the aligner used for alignment calculations",
          command::Parameter
          (
            "aligner_class",
            "choice of aligner class",
            command::ParameterCheckEnumerate< align::MultipleAlignerClasses< biol::AABase> >(),
            align::GetMultipleAlignerClasses< biol::AABase>().e_AlignerProgressiveDP.GetName()
          )
        )
      ),
      m_ScoreAlignmentFlag
      (
        new command::FlagStatic
        (
          "score_alignment",
          "score the given alignment",
          command::Parameter
          (
            "alignment file",
            "any sequence alignment in a known format, which is determined by the file extension",
            command::ParameterCheckFileExistence(),
            ""
          )
        )
      )
    {
      // flags with parameters: weights and gap scores
      for
      (
        score::AAAssignments::const_iterator
          itr( score::GetAAAssignments().Begin()), itr_end( score::GetAAAssignments().End());
        itr != itr_end;
        ++itr
      )
      {
        // create description for score parameter
        std::stringstream str_description;
        str_description << "weight for " << itr->GetName() << " score";
        if( !score::GetAAAssignments().GetFileExtension( *itr).empty())
        {
          str_description << ", <sequence>" << score::GetAAAssignments().GetFileExtension( *itr) << " must exist";
        }
        // create flag with static parameter list and parameter with allowed range
        util::ShPtr< command::FlagStatic> pair_score_flag
        (
          new command::FlagStatic
          (
            itr->GetName(),
            "weight for " + itr->GetName() + " score",
            command::Parameter
            (
              "score_weight",
              str_description.str(),
              command::ParameterCheckRanged< double>( double( 0.0), double( 1.0)),
              "1.0"
            )
          )
        );
        // push back parameter list and parameter in util::ShPtrVector
        m_PairScoreFlags[ *itr] = pair_score_flag;
      }

      // flags with parameters: gap scores
      for( size_t i( 0); i < s_NumberGapPenalties; ++i)
      {
        // create flag with static parameter list and parameter with allowed range
        util::ShPtr< command::FlagStatic> gap_penalty_flag
        (
          new command::FlagStatic
          (
            s_GapPenaltiesDescription[ i], s_GapPenaltiesDescription[ i] + " penalty",
            command::Parameter
            (
              "gap_penalty",
              s_GapPenaltiesDescription[ i] + " penalty",
              command::ParameterCheckRanged< double>( double( -10.0), double( 0.0)),
              s_DefaultGapPenalties[ i]
            )
          )
        );
        // push back parameter list and parameter in util::ShPtrVector
        m_GapPenaltyFlags.PushBack( gap_penalty_flag);
      }
    }

    //! @brief Clone function
    //! @return pointer to new Alignment
    Alignment *Alignment::Clone() const
    {
      return new Alignment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &Alignment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the bcl::commons name
    //! @return string for the bcl::commons name of that application
    std::string Alignment::GetBCLScopedName() const
    {
      return "BCL::Align";
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string Alignment::GetDescription() const
    {
      return "Sequence alignment with choice of alignment algorithm, scores and weights";
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> Alignment::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // input ascii and fasta files
      sp_cmd->AddFlag( m_FastaListFlag);

      // output format flag
      sp_cmd->AddFlag( align::HandlerClasses< biol::AABase>::GetFlagOutputFormats());

      // prefix to the output alignment file
      sp_cmd->AddFlag( m_OutputprefixFlag);

      sp_cmd->AddFlag( m_ScoreAlignmentFlag);

      // flags with parameters: weights and gap scores
      for
      (
        storage::Map< score::AAAssignment, util::ShPtr< command::FlagInterface> >::const_iterator
          itr( m_PairScoreFlags.Begin()), itr_end( m_PairScoreFlags.End());
        itr != itr_end;
        ++itr
      )
      {
        sp_cmd->AddFlag( itr->second);
      }

      // flags with parameters: gap scores
      for
      (
        util::ShPtrVector< command::FlagInterface>::const_iterator
          itr( m_GapPenaltyFlags.Begin()), itr_end( m_GapPenaltyFlags.End());
        itr != itr_end;
        ++itr
      )
      {
        sp_cmd->AddFlag( *itr);
      }

      // aligner to calculate alignments
      sp_cmd->AddFlag( m_FlagAligner);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &Alignment::GetReadMe() const
    {
      static const std::string readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::Align, terms of use, appropriate citation, installation "
        "procedures, BCL::Align execution, technical support, and future research directions.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::Align?\n"
        "BCL::Align is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Commons.  BCL::Align is a sequence alignment utility which "
        "allows the scoring function to be flexibly specified.\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::Align.\n"
        "When using BCL::Align in a publication, please cite the following publication describing the application's "
        "development:\n"
        "E. Dong, J. Smith, S. Heinze, N. Alexander, J. Meiler. BCL::Align-Sequence alignment and fold recognition "
        "with a custom scoring function online. Gene 422, 41 (Oct, 2008).\n"
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL::Align.\n"
        "Consists of four steps :\n"
        "1) Gather your fasta files for the sequences you desire to align.\n"
        "2) Run any needed prediction algorithms (primarily secondary structure prediction methods such as JUFO or PSIPRED) on the sequences. This is dependent on the scoring functions you choose to use.\n"
        "3) Set weights for each of the scoring function parameters you would like to include while aligning your sequences.\n"
        "4) Use BCL::Align to perform sequence alignment by specifying command line arguments.\n"
        "\n"
        "ADDITIONAL INFORMATION\n"
        "\n"
        "Information concerning syntax and flags can be obtained by typing : <bcl.exe> Alignment -help\n"
        "\n"
        "For further explanation, examples of the flags, example command lines, input and output format information,\n"
        "and example output please see the documentation file.\n"
        "\n"
        "For more general information about the product, type : <bcl.exe> Alignment -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::Align.\n"
        "BCL::Align is under ongoing further development. A word-based alignment algorithm is being developed, which "
        "will improve BCL::Align's capability of detecting locally homologous sequences and can speed global alignment. "
        "In addition, a neural network based assignment score that takes the local sequence composition into account, "
        "is currently in preparation.\n"
        "\n"
        + DefaultSectionSeparator() +
        "IX. EXAMPLE COMMANDLINE"
        "<bcl.exe> Alignment -blast 1.0 -blosum62 1.0 -extend_boundary_gap 0.0 -extend_gap -0.5 -fastas 1IE9A.fasta "
        "1IE9A.fasta -hydrophobicity 1.0 -isoelectric 1.0 -open_boundary_gap 0.0 -open_gap -1.0 -outputformat pir "
        "-outputprefix 1IE9A -pam250 1.0 -polarizability 1.0 -psipred 1.0 -tfe_white 1.0 -volume 1.0 "
        + DefaultSectionSeparator()
      );
      return readme;
    }

    //! Main
    int Alignment::Main() const
    {
      // read options
      BCL_MessageStd( "Detect scores to be used:\n");
      function::BinarySum< const biol::AABase, const biol::AABase, double> score;
      size_t score_count( 0);
      double sum_score_weight( 0);
      storage::Set< score::AAAssignment> selected_scores;

      // read pair scores and weights
      for
      (
        storage::Map< score::AAAssignment, util::ShPtr< command::FlagInterface> >::const_iterator
          itr( m_PairScoreFlags.Begin()), itr_end( m_PairScoreFlags.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->second->GetFlag())
        {
          double score_weight( itr->second->GetFirstParameter()->GetNumericalValue< double>());
          if( score_weight != 0.0)
          {
            AddScore
            (
              itr->first,
              score_weight,
              score,
              util::GetLogger()
            );
            selected_scores.Insert( itr->first);
            sum_score_weight += score_weight;
            score_count++;
          }
        }
      }

      // check scores
      BCL_Assert( score_count != 0, "At least one score needs to be input");
      // normalize scores
      if( sum_score_weight != 1)
      {
        BCL_MessageStd
        (
          "Normalize score weights from " + util::Format().W( 6).FFP( 2)( sum_score_weight) + " to 1.00"
        );
        score /= sum_score_weight;
      }

      // output gap penalties
      BCL_MessageStd( "Current gap penalties:");

      for( size_t i( 0); i < s_NumberGapPenalties; ++i)
      {
        BCL_MessageStd
        (
          s_GapPenaltiesDescription[ i] + " penalty : " +
          util::Format().W( 6).FFP( 2)( m_GapPenaltyFlags( i)->GetFirstParameter()->GetNumericalValue< double>())
        );
      }

      // build scoring function
      score::AssignmentWithGap< biol::AABase> assignment_score
      (
        util::CloneToShPtr( score),
        m_GapPenaltyFlags( 0)->GetFirstParameter()->GetNumericalValue< double>(),
        m_GapPenaltyFlags( 1)->GetFirstParameter()->GetNumericalValue< double>(),
        m_GapPenaltyFlags( 2)->GetFirstParameter()->GetNumericalValue< double>(),
        m_GapPenaltyFlags( 3)->GetFirstParameter()->GetNumericalValue< double>()
      );

      // create Pair alignment_and_score to hold the alignment and score and calculate the alignment
      storage::Pair< align::AlignmentNode< biol::AABase>, double>
        alignment_and_score;

      if( m_ScoreAlignmentFlag->GetFirstParameter()->GetWasSetInCommandLine())
      {
        const std::string filename( m_ScoreAlignmentFlag->GetFirstParameter()->GetValue());
        const std::string format_ext( io::File::GetLastExtension( filename));
        const align::HandlerClasses< biol::AABase>::HandlerClass handler( format_ext);
        io::IFStream read;
        io::File::MustOpenIFStream( read, filename);
        alignment_and_score.First() = *( *handler)->ReadAlignment( read, biol::AASequence());
        io::File::CloseClearFStream( read);
        score::AlignmentAssignment< biol::AABase> score_alignment( util::CloneToShPtr( assignment_score));
        alignment_and_score.Second() = score_alignment( alignment_and_score.First());
      }
      else
      {
        // move filenames from the paramlist to sequence_names
        storage::Vector< std::string> sequence_names( m_FastaListFlag->GetStringList());
        for
        (
          storage::Vector< std::string>::iterator name_itr( sequence_names.Begin()), name_itr_end( sequence_names.End());
          name_itr != name_itr_end;
          ++name_itr
        )
        {
          // calculate position of last point and last slash in the filename
          const size_t pos_point( name_itr->rfind(".", name_itr->length()));
          const size_t pos_slash( name_itr->rfind("/", name_itr->length()));

          // remove the extension from the filenames
          // if a dot exists and either no slash exists or the slash comes before the dot
          if( pos_point != std::string::npos && ( pos_slash == std::string::npos || pos_point > pos_slash))
          {
            ( *name_itr) = name_itr->substr( 0, pos_point);
          }
        }

        // read AASequences
        BCL_MessageStd( "Read sequence data:");
        storage::Vector< biol::AASequence> sequences;
        sequences.AllocateMemory( sequence_names.GetSize());
        {
          for
          (
            std::vector< std::string>::const_iterator
              name_itr( sequence_names.Begin()),
              name_itr_end( sequence_names.End());
            name_itr != name_itr_end;
            ++name_itr
          )
          {
            sequences.PushBack( ReadAllSequenceData( *name_itr, selected_scores, util::GetLogger()));
            sequences.LastElement().WriteFasta( util::GetLogger(), 100);
          }
        }

        // test if Z-score parameters are to be computed
        if( s_ZScoreParameter)
        {
          ComputeParametersToConvertPairScoresToZScores( sequences, selected_scores, util::GetLogger());
          return 0;
        }

        // convert sequences to single-sequence alignments
        util::ShPtrList< align::AlignmentInterface< biol::AABase> > alignments;
        for
        (
          std::vector< biol::AASequence>::const_iterator
            sequence_itr( sequences.Begin()),
            sequence_itr_end( sequences.End());
          sequence_itr != sequence_itr_end;
          ++sequence_itr
        )
        {
          util::ShPtr< align::SequenceInterface< biol::AABase> > sp_sequence_interface( sequence_itr->Clone());
          align::AlignmentLeaf< biol::AABase> alignment_leaf( sp_sequence_interface);
          util::ShPtr< align::AlignmentInterface< biol::AABase> > sp_alignment_interface( alignment_leaf.Clone());
          alignments.Append( sp_alignment_interface);
        }

        // create aligner to align the sequences and set scoring function
        util::ShPtr< align::MultipleAlignerInterface< biol::AABase> >
          aligner( *align::MultipleAlignerClasses< biol::AABase>::MultipleAlignerClass( m_FlagAligner->GetFirstParameter()->GetValue()));
        aligner->SetScoringFunction( assignment_score);

        alignment_and_score = aligner->AlignMultiple( alignments);

        // output
        size_t shortest( util::GetUndefined< size_t>());
        for
        (
          std::vector< biol::AASequence>::const_iterator
             sequence_itr( sequences.Begin()),
             sequence_itr_end( sequences.End());
           sequence_itr != sequence_itr_end;
           ++sequence_itr
        )
        {
          if( sequence_itr->GetSize() < shortest)
          {
            shortest = sequence_itr->GetSize();
          }
        }

        BCL_MessageStd( "Shortest Sequence:\t" + util::Format()( shortest));
        BCL_MessageStd
        (
          "Normalized Score:\t" + util::Format()( alignment_and_score.Second() / shortest)
        );
      }
      BCL_MessageStd( "Score:\t" + util::Format()( alignment_and_score.Second()));

      // no output format given -> write pir
      if( align::HandlerClasses< biol::AABase>::GetFlagOutputFormats()->GetParameterList().IsEmpty())
      {
        // write the Alignment in pir format
        align::HandlerPIR< biol::AABase> handler;
        handler.WriteAlignment( util::GetLogger(), alignment_and_score.First());
      }
      else
      {
        //iterate over all given output formats
        for
        (
          util::ShPtrVector< command::ParameterInterface>::const_iterator
            format_itr( align::HandlerClasses< biol::AABase>::GetFlagOutputFormats()->GetParameterList().Begin()),
            format_itr_end( align::HandlerClasses< biol::AABase>::GetFlagOutputFormats()->GetParameterList().End());
          format_itr != format_itr_end;
          ++format_itr
        )
        {
          // get align::Handler from command line parameter
          util::ShPtr< align::HandlerInterface< biol::AABase> >
            handler( align::HandlerClasses< biol::AABase>::HandlerClass( ( *format_itr)->GetValue())->HardCopy());

          BCL_MessageStd( "output in format: " + ( *format_itr)->GetValue());

          //stream to be used for output of the alignment
          std::ostream *filestream;
          io::OFStream write;
          bool writing_to_std_cout( false);

          //if there was no prefix for output given, use the output stream
          if
          (
            m_OutputprefixFlag->GetParameterList().IsEmpty() ||
            m_OutputprefixFlag->GetFirstParameter()->GetValue().empty()
          )
          {
            filestream = &util::GetLogger();
            writing_to_std_cout = true;
          }
          //if there was a prefix given, open a write stream {filename}{outputformat_file_extension}
          else
          {
            const std::string file_name( m_OutputprefixFlag->GetFirstParameter()->GetValue() + handler->GetFileExtension());

            io::File::MustOpenOFStream( write, file_name);
            filestream = &write;
          }

          // if output is not written to a file, we add the identifier after options and before the actual alignment results
          // so that scripts can locate the exact position where the alignments start
          if( writing_to_std_cout)
          {
            alignment_and_score.First().WriteIdentifier( *filestream, 0);
          }

          // score to be used, if the handler supports it
          handler->SetAssignmentScore( util::CloneToShPtr( assignment_score));

          // write alignment
          handler->WriteAlignment( *filestream, alignment_and_score.First());

          //clean up stream
          io::File::CloseClearFStream( write);
        }
      }

      // end
      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Alignment::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &Alignment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

    //! read sequence profile data
    biol::AASequence Alignment::ReadAllSequenceData
    (
      const std::string &CODE,
      const storage::Set< score::AAAssignment> &READ,
      std::ostream &OUTSTREAM
    ) const
    {
      io::IFStream read;

      // read seq1 from fasta
      io::File::MustOpenIFStream( read, CODE + ".fasta");

      //TODO: Find a solution for OUTSTREAM problem
      biol::AASequence seq( biol::AASequenceFactory::BuildSequenceFromFASTA( read));
      io::File::CloseClearFStream( read);

      // read blast profile
      if( READ.Contains( score::GetAAAssignments().e_BLAST))
      {
        io::File::MustOpenIFStream( read, CODE + score::GetAAAssignments().GetFileExtension( score::GetAAAssignments().e_BLAST));
        biol::BlastProfileHandler::ReadProfileForAASequence( read, seq);
        io::File::CloseClearFStream( read);
      }

      // read secondary structure prediction methods
      if( READ.Contains( score::GetAAAssignments().e_PSIPRED))
      {
        ReadSSPrediction( seq, CODE, sspred::GetMethods().e_PSIPRED); // read psipred ss prediction
      }
      if( READ.Contains( score::GetAAAssignments().e_JUFO))
      {
        ReadSSPrediction( seq, CODE, sspred::GetMethods().e_JUFO9D); // read jufo9d ss prediction
      }
      // if jufo is not provided calculate it on the go for correct one letter SS display in the output
      else
      {
        if( READ.Contains( score::GetAAAssignments().e_BLAST))
        {
          BCL_MessageStd
          (
            "jufo will be calculated since blast was given for that sequence: " + CODE
          );
          sspred::JUFO::Calculate( seq);
        }
      }
      if( READ.Contains( score::GetAAAssignments().e_SAM))
      {
        ReadSSPrediction( seq, CODE, sspred::GetMethods().e_SAM      ); // read sam ss prediction
      }
      if( READ.Contains( score::GetAAAssignments().e_TMHMM))
      {
        ReadSSPrediction( seq, CODE, sspred::GetMethods().e_TMHMM    ); // read tmhmm ss prediction
      }
      if( READ.Contains( score::GetAAAssignments().e_TMMOD))
      {
        ReadSSPrediction( seq, CODE, sspred::GetMethods().e_TMMOD    ); // read tmmod ss prediction
      }
      if( READ.Contains( score::GetAAAssignments().e_B2TMPRED))
      {
        ReadSSPrediction( seq, CODE, sspred::GetMethods().e_B2TMPRED); // read b2tmpred ss prediction
      }
      if( READ.Contains( score::GetAAAssignments().e_PROFTMB))
      {
        ReadSSPrediction( seq, CODE, sspred::GetMethods().e_PROFTMB); // read proftmb ss prediction
      }
      if( READ.Contains( score::GetAAAssignments().e_CONPRED))
      {
        ReadSSPrediction( seq, CODE, sspred::GetMethods().e_CONPRED); // read compred ss prediction
      }

      return seq;
    }

    //! AddScore
    bool Alignment::AddScore
    (
      const score::AAAssignment &SCORE,
      const double SCORE_WEIGHT,
      function::BinarySum< const biol::AABase, const biol::AABase, double> &SCORE_FCT,
      std::ostream &OUTSTREAM
    ) const
    {
      util::ShPtr< math::ZScore> z_score( score::GetAAAssignments().GetZScore( SCORE));

      if( !z_score.IsDefined())
      {
        return false;
      }
      else
      {
        z_score = z_score.HardCopy();
      }

      // correct Z-score by the weight
      *z_score *= SCORE_WEIGHT;

      // add to scoring function
      SCORE_FCT += function::BinaryAdapter< const biol::AABase, const biol::AABase, double, const double, double>( *SCORE, z_score);

      if( !s_ZScoreParameter)
      {
        OUTSTREAM << "Weight for option " << SCORE->GetName() << " is set to "
          << util::Format().W( 6).FFP( 2)( SCORE_WEIGHT) << '\n';
      }
      return true;
    }

    //! estimate parameters to convert pair scores to Z-scores
    void Alignment::ComputeParametersToConvertPairScoresToZScores
    (
      const storage::Vector< biol::AASequence> &SEQUENCES,
      const storage::Set< score::AAAssignment> &SELECTED_SCORES,
      std::ostream &OUTSTREAM
    ) const
    {
      for
      (
        storage::Set< score::AAAssignment>::const_iterator
          itr( SELECTED_SCORES.Begin()), itr_end( SELECTED_SCORES.End());
        itr != itr_end;
        ++itr
      )
      {
        // score function
        function::BinarySum< const biol::AABase, const biol::AABase, double> score;
        AddScore( *itr, double( 1), score, OUTSTREAM);
        score::AssignmentWithGap< biol::AABase> assign( util::CloneToShPtr( score), 0, 0, 0, 0);

        // accumulate data for average and stddev computation
        math::RunningAverageSD< double> mean_sd_score;

        // iterate over every sequence
        for
        (
          storage::Vector< biol::AASequence>::const_iterator
            seq_itr_a( SEQUENCES.Begin()),
            seq_itr_end( SEQUENCES.End());
          seq_itr_a != seq_itr_end;
          ++seq_itr_a
        )
        {
          // iterate over every sequence
          for
          (
            storage::Vector< biol::AASequence>::const_iterator seq_itr_b( seq_itr_a);
            seq_itr_b != seq_itr_end;
            ++seq_itr_b
          )
          {
            // iterate over every residue
            for
            (
              biol::AASequence::const_iterator
                aa_itr_a( seq_itr_a->Begin()),
                aa_itr_end_a( seq_itr_a->End());
              aa_itr_a != aa_itr_end_a;
              ++aa_itr_a
            )
            {
              // iterate over every residue
              for
              (
                biol::AASequence::const_iterator
                  aa_itr_b( seq_itr_b->Begin()),
                  aa_itr_end_b( seq_itr_b->End());
                aa_itr_b != aa_itr_end_b;
                ++aa_itr_b
              )
              {
                const double curr_score
                (
                  assign
                  (
                    util::ToSiPtr( **aa_itr_a),
                    util::ToSiPtr( **aa_itr_b)
                  )
                );
                mean_sd_score += curr_score;
              } // end aa_itr_b
            } // end aa_itr_a
          } // end seq_itr_b
        } // end seq_itr_a

        // output
        BCL_MessageStd
        (
          "  storage::Pair< double, double>( " + util::Format().W( 8).FFP( 4)( mean_sd_score.GetAverage()) + ").F( " +
          util::Format().W( 8).FFP( 4)( mean_sd_score.GetStandardDeviation()) + ")).F( // " + itr->GetName()
        );
      }

      // end
      return;
    }

    //! read secondary structure prediction
    biol::AASequence &Alignment::ReadSSPrediction
    (
      biol::AASequence &SEQUENCE,
      const std::string &CODE,
      const sspred::Method &METHOD
    ) const
    {
      io::IFStream read;

      io::File::MustOpenIFStream( read, CODE + ( *METHOD)->GetFileExtension());
      sspred::MethodHandler::ReadPredictionsForAASequence( read, SEQUENCE, METHOD);
      io::File::CloseClearFStream( read);

      return SEQUENCE;
    }

    // option strings gap penalties
    const std::string Alignment::s_GapPenaltiesDescription[ s_NumberGapPenalties] =
    {
      "open_gap", "extend_gap", "open_boundary_gap", "extend_boundary_gap"
    };

    // default gap penalties
    const std::string Alignment::s_DefaultGapPenalties[ 4] =
    {
      "-1", "-1", "0", "0"
    };

    const ApplicationType Alignment::Alignment_Instance
    (
      GetAppGroups().AddAppToGroup( new Alignment(), GetAppGroups().e_Sequence)
    );

  } // namespace app
} // namespace bcl
