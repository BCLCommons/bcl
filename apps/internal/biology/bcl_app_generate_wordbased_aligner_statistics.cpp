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
#include "align/bcl_align_aligner_progressive.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "biol/bcl_biol_blast_profile_handler.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_extension.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "function/bcl_function_binary_sum.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_roc_curve.h"
#include "score/bcl_score_aa_assignment_blast_profile.h"
#include "score/bcl_score_aa_assignment_blosum.h"
#include "score/bcl_score_aa_assignment_identity.h"
#include "score/bcl_score_aa_assignment_pam.h"
#include "score/bcl_score_aa_assignment_phat.h"
#include "score/bcl_score_aa_assignment_property.h"
#include "score/bcl_score_aa_assignment_ss_prediction.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes

namespace bcl
{
  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GenerateWordbasedAlignerStatistics
    //! @brief generate statistics for wordbased aligner
    //!
    //! @author heinzes1
    //! @date Mar 25, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GenerateWordbasedAlignerStatistics :
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
      static const std::string s_PairScoresFileExtensions[]; //!< option strings pair scores

      util::ShPtr< command::FlagInterface> m_FlagFastaList; //!< input fasta files
      util::ShPtr< command::FlagInterface> m_FlagGoldStandardAlignment; //!< input gold standard alignment
      util::ShPtr< command::FlagInterface> m_FlagCalculateExtendedHits; //!< calculate and print extended hits

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      GenerateWordbasedAlignerStatistics();

      //! @brief Clone function
      //! @return pointer to new GenerateWordbasedAlignerStatistics
      GenerateWordbasedAlignerStatistics *Clone() const
      {
        return new GenerateWordbasedAlignerStatistics( *this);
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

      //! @brief the Main function
      //! @return error code - 0 for success
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

      //! @brief initializes the command object for that executable
      //! @return ShPtr to a command object
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief calculates word statistics from alignments (and optional gold standard alignment) given on commandline
      //! @param ALIGNMENTS non-const reference to alignments to calculate statistics from
      //! @param SCORE scoring function
      void CalculateExtendedHits
      (
        util::ShPtrList< align::AlignmentInterface< biol::AABase> > &ALIGNMENTS,
        const score::AssignmentWithGap< biol::AABase> &SCORE
      ) const;

      //! @brief calculates word statistics from alignments (and optional gold standard alignment) given on commandline
      //! @param FILENAMES filenames for the sequence files
      //! @param SEQUENCES sequences (must be same than in alignments)
      //! @param ALIGNMENTS non-const reference to alignments to calculate statistics from
      //! @param SCORE scoring function
      void CalculateWordStatistics
      (
        const storage::Vector< std::string> &FILENAMES,
        const util::ShPtrVector< biol::AASequence> &SEQUENCES,
        util::ShPtrList< align::AlignmentInterface< biol::AABase> > &ALIGNMENTS,
        const score::AssignmentWithGap< biol::AABase> &SCORE
      ) const;

      //! @brief test the given list of AlignmentHits against a gold standard alignment
      //! @param ALIGNMENT_GOLD_STANDARD the gold standard alignment to compare the hits with
      //! @param HITS the list of hit alignments
      //! @param ASSIGN_SCORE the assignment scoring function
      //! @return a list of double+bool encoding the score of a hit and if its true positive
      storage::List< storage::Pair< double, bool> > TestAlignmentHits
      (
        const align::AlignmentInterface< biol::AABase> &ALIGNMENT_GOLD_STANDARD,
        const storage::List< align::AlignmentHit< biol::AABase> > &HITS,
        const score::AssignmentWithGap< biol::AABase> &ASSIGN_SCORE
      ) const;

      //! @brief print information about the read sequences to the current logger
      //! @param FILENAMES the files read
      //! @param SEQUENCES the sequences read from the files
      void PrintSequenceData
      (
        const storage::Vector< std::string> &FILENAMES,
        const util::ShPtrVector< biol::AASequence> &SEQUENCES
      ) const;

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

      //! @brief read secondary structure prediction
      //! @param SEQUENCE the AASequence to read/modify
      //! @param CODE the filename
      //! @param SCORE the pair score to determine the file extension
      //! @param METHOD the ss prediction method to use
      //! @return the AASequence passed as SEQUENCE
      biol::AASequence &ReadSSPrediction
      (
        biol::AASequence &SEQUENCE,
        const std::string &CODE,
        const s_PairScores SCORE,
        const sspred::Method METHOD
      ) const;

      static const ApplicationType GenerateWordbasedAlignerStatistics_Instance; //!< static instance

    private:

      //! @brief operator taking two pairs and comparing them
      //! @param LHS left pair
      //! @param RHS right pair
      //! @return if the score of LHS > RHS
      static bool LessThanPairFirst
      (
        const storage::Pair< double, bool> &LHS,
        const storage::Pair< double, bool> &RHS
      )
      {
        return LHS.First() < RHS.First();
      }

      //! @brief operator taking two pairs and comparing them
      //! @param LHS left pair
      //! @param RHS right pair
      //! @return if the score of LHS > RHS
      static bool LessThanPairSecond
      (
        const storage::Pair< align::AlignmentNode< biol::AABase>, double> &LHS,
        const storage::Pair< align::AlignmentNode< biol::AABase>, double> &RHS
      )
      {
        return LHS.Second() < RHS.Second();
      }

    }; // class GenerateWordStatistics

    //! @brief default constructor
    GenerateWordbasedAlignerStatistics::GenerateWordbasedAlignerStatistics() :
      m_FlagFastaList
      (
        new command::FlagDynamic
        (
          "fastas",
          "these are fasta files to be optimally aligned and then compared used to calculate the words",
          command::Parameter( "fasta file", "any fasta file", command::ParameterCheckExtension( ".fasta")),
          2 // minimum of 2 parameters
        )
      ),
      m_FlagGoldStandardAlignment
      (
        new command::FlagStatic
        (
          "gold_standard_alignment",
          "gold standard alignment in pir format",
          command::Parameter( "pir file", "any pir alignment file", command::ParameterCheckFileExistence(), "")
        )
      ),
      m_FlagCalculateExtendedHits
      (
        new command::FlagStatic( "calculate_extended_hits", "calculate and print extended hits")
      )
    {
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int GenerateWordbasedAlignerStatistics::Main() const
    {
      BCL_Assert
      (
        !( m_FlagGoldStandardAlignment->GetFirstParameter()->GetWasSetInCommandLine()
          && m_FlagCalculateExtendedHits->GetFlag()),
        "Error: -gold_standard_alignment and -calculate_extended_hits CANNOT be used at the same time!"
      );

      // read and print fasta files
      BCL_MessageStd( "Read fasta data: ");
      const storage::Vector< std::string> fasta_filenames( m_FlagFastaList->GetStringList());
      util::ShPtrVector< biol::AASequence>
        fasta_sequences( biol::AASequenceFactory::ReadSequenceData( m_FlagFastaList->GetStringList()));
      PrintSequenceData( fasta_filenames, fasta_sequences);
      BCL_MessageStd( "Found " + util::Format()( fasta_sequences.GetSize()) + " fasta sequences.");

      // do per-sequence work
      util::ShPtrList< align::AlignmentInterface< biol::AABase> > alignments;
      storage::Vector< std::string>::const_iterator
        itr_filename( fasta_filenames.Begin()),
        itr_filename_end( fasta_filenames.End());
      for
      (
        util::ShPtrVector< biol::AASequence>::iterator itr( fasta_sequences.Begin()), itr_end( fasta_sequences.End());
        itr != itr_end && itr_filename != itr_filename_end;
        ++itr, ++itr_filename
      )
      {
        const std::string filename_without_extension( io::File::RemoveLastExtension( *itr_filename));

        // read ss predictions
        io::IFStream read;
        io::File::MustOpenIFStream( read, filename_without_extension + s_PairScoresFileExtensions[ BLAST]);
        biol::BlastProfileHandler::ReadProfileForAASequence( read, **itr);
        io::File::CloseClearFStream( read);
        ReadSSPrediction( **itr, filename_without_extension, PSIPRED, sspred::GetMethods().e_PSIPRED);
        ReadSSPrediction( **itr, filename_without_extension, JUFO, sspred::GetMethods().e_JUFO);

        // create single sequence alignments from the sequences
        util::ShPtr< align::AlignmentInterface< biol::AABase> >
          sp_alignment( new align::AlignmentLeaf< biol::AABase>( *itr));
        alignments.PushBack( sp_alignment);
      }

      // create assignment score
      function::BinarySum< const biol::AABase, const biol::AABase, double> assign_score;
      AddScore( BLOSUM45, 0.13, assign_score);
      AddScore( BLAST, 0.203, assign_score);   // SAM method is not supported by BCL anymore; redistribute the weight of
      AddScore( PSIPRED, 0.223, assign_score); // 0.07 for SAM equally over the other ss-prediction scores
      AddScore( JUFO, 0.203, assign_score);    // BLAST and JUFO: 0.18+(0.07/3); PSIPRED: 0.2+(0.07/3)
      AddScore( STERICAL_PARAMETER, 0.034, assign_score); // the total weight of 0.24 for all chemical properties
      AddScore( POLARIZABILITY, 0.034, assign_score);     // is equally distributed over them
      AddScore( VOLUME, 0.034, assign_score);
      AddScore( HYDROPHOBICITY, 0.034, assign_score);
      AddScore( ISOELECTRIC_POINT, 0.034, assign_score);
      AddScore( TFE_WHITE, 0.034, assign_score);
      AddScore( TFE_ENGELMAN, 0.034, assign_score);

      score::AssignmentWithGap< biol::AABase> assign_gap_score
      (
        util::CloneToShPtr( assign_score),
        -0.8, // ENCLOSED_SINGLE_GAP
        -1.4, // ENCLOSED_MULTIPLE_GAP
        -1.7, // BOUNDARY_SINGLE_GAP
        -0.1  // BOUNDARY_MULTIPLE_GAP
      );

      if( m_FlagCalculateExtendedHits->GetFlag())
      {
        CalculateExtendedHits( alignments, assign_gap_score);
      }
      else
      {
        CalculateWordStatistics( fasta_filenames, fasta_sequences, alignments, assign_gap_score);
      }

      return 0;
    }

    //! @brief initializes the command object for that executable
    //! @return ShPtr to a command object
    util::ShPtr< command::Command> GenerateWordbasedAlignerStatistics::InitializeCommand() const
    {
      // initialize ShPtr to a new command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add all flags
      sp_cmd->AddFlag( m_FlagFastaList);
      sp_cmd->AddFlag( m_FlagGoldStandardAlignment);
      sp_cmd->AddFlag( m_FlagCalculateExtendedHits);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      return sp_cmd;
    }

    //! @brief calculates word statistics from alignments (and optional gold standard alignment) given on commandline
    //! @param ALIGNMENTS non-const reference to alignments to calculate statistics from
    //! @param SCORE scoring function
    void GenerateWordbasedAlignerStatistics::CalculateExtendedHits
    (
      util::ShPtrList< align::AlignmentInterface< biol::AABase> > &ALIGNMENTS,
      const score::AssignmentWithGap< biol::AABase> &SCORE
    ) const
    {
      align::AlignerWordbased< biol::AABase> aligner( SCORE, 2); // construct with assign score

      // create first alignment and alignment_list for the wordbased aligner
      util::ShPtr< align::AlignmentInterface< biol::AABase> > first_alignment( ALIGNMENTS.FirstElement());
      util::ShPtrList< align::AlignmentInterface< biol::AABase> > alignments_except_first( ++ALIGNMENTS.Begin(), ALIGNMENTS.End());

      // align alignment_a to all alignments in the list
      storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, double> >
        result_alignments( aligner.AlignPairwise( first_alignment, alignments_except_first));

      // sort by score
      result_alignments.Sort( &LessThanPairSecond);

      // remove all except top 100
      storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, double> >::iterator
        itr_begin( result_alignments.Begin()),
        itr_end( itr_begin);
      storage::AdvanceIterator( itr_end, result_alignments.End(), result_alignments.GetSize() - 200);
      result_alignments.Remove( itr_begin, itr_end);

      // print alignments
      align::HandlerPIR< biol::AABase> handler_pir;
      for
      (
        storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, double> >::const_iterator
          itr( result_alignments.Begin()),
          itr_end( result_alignments.End());
        itr != itr_end;
        ++itr
      )
      {
        util::GetLogger() << "Alignment: score=" << itr->Second();
//        util::GetLogger() << "; child_alignments_size=" << itr->First().GetChildAlignments().GetSize();
        util::GetLogger() << '\n';
        handler_pir.WriteAlignment( util::GetLogger(), itr->First());
//        util::GetLogger() << "child_alignment_a:\n";
//        handler_pir.WriteAlignment( util::GetLogger(), *itr->First().GetChildAlignments().FirstElement());
//        util::GetLogger() << "child_alignment_b:\n";
//        handler_pir.WriteAlignment( util::GetLogger(), *itr->First().GetChildAlignments().LastElement());
      }
    }

    //! @brief calculates word statistics from alignments (and optional gold standard alignment) given on commandline
    //! @param FILENAMES filenames for the sequence files
    //! @param SEQUENCES sequences (must be same than in alignments)
    //! @param ALIGNMENTS non-const reference to alignments to calculate statistics from
    //! @param SCORE scoring function
    void GenerateWordbasedAlignerStatistics::CalculateWordStatistics
    (
      const storage::Vector< std::string> &FILENAMES,
      const util::ShPtrVector< biol::AASequence> &SEQUENCES,
      util::ShPtrList< align::AlignmentInterface< biol::AABase> > &ALIGNMENTS,
      const score::AssignmentWithGap< biol::AABase> &SCORE
    ) const
    {
      std::string filename_prefix;
      for
      (
        storage::Vector< std::string>::const_iterator itr( FILENAMES.Begin()), itr_end( FILENAMES.End());
        itr != itr_end;
        ++itr
      )
      {
        filename_prefix.append( io::File::RemovePath( io::File::RemoveLastExtension( *itr)) + "-");
      }

      // start writing the gnuplot file
      const std::string plot_filename( filename_prefix + "generate_wordbased_aligner_statistics");
      io::OFStream write_plot;
      io::File::MustOpenOFStream( write_plot, plot_filename + ".plot");

      // write header to plot ROC curve with gnuplot
      write_plot << "set terminal png enhanced # transparent size 2160,640\n";
      write_plot << "set output \"" << plot_filename << ".png\"\n";
      write_plot << "set encoding iso\n";
      write_plot << "set key right bottom\n";
      write_plot << "set title \"roc curve\"\n";
      write_plot << "set xlabel \"false positive rate\"\n";
      write_plot << "set xrange [0:1]\n";
      write_plot << "set ylabel \"true positive rate\"\n";
      write_plot << "set yrange [0:1]\n";
      write_plot << "set multiplot\n";
      write_plot << "plot x with lines linetype 1 linewidth 1 linecolor rgb 'black'\n";
      write_plot << "set style line 1 linetype 1 linewidth 0 pointtype 1 linecolor rgb 'red'\n";
      write_plot << "set style line 2 linetype 1 linewidth 0 pointtype 1 linecolor rgb 'green'\n";
      write_plot << "set style line 3 linetype 1 linewidth 0 pointtype 1 linecolor rgb 'blue'\n";
      write_plot << "set style line 4 linetype 1 linewidth 0 pointtype 1 linecolor rgb 'magenta'\n";
      write_plot << "set style line 5 linetype 1 linewidth 0 pointtype 1 linecolor rgb 'cyan'\n";
      write_plot << "set style line 6 linetype 1 linewidth 0 pointtype 1 linecolor rgb 'black'\n";
      write_plot << "set style line 7 linetype 1 linewidth 0 pointtype 1 linecolor rgb 'orange'\n";
      write_plot << "set style line 8 linetype 1 linewidth 0 pointtype 1 linecolor rgb 'brown'\n";

      // create dynamic programming/progressive aligner to get an optimal alignment, assumed to be the ideal alignment
      util::ShPtr< align::AlignerDynamicProgramming< biol::AABase> >
        sp_aligner_dyn_prog( new align::AlignerDynamicProgramming< biol::AABase>());
      align::AlignerProgressive< biol::AABase> aligner_progressive( sp_aligner_dyn_prog, SCORE);

      // calculate "ideal" alignment and print result alignment to user
      const storage::Pair< align::AlignmentNode< biol::AABase>, double>
        calculated_gold_standard_alignment_score_pair( aligner_progressive.AlignMultiple( ALIGNMENTS));
      align::HandlerPIR< biol::AABase> pir_handler;
      pir_handler.WriteAlignment( util::GetLogger(), calculated_gold_standard_alignment_score_pair.First());

      // create first alignment and alignment_list for the wordbased aligner
      util::ShPtr< align::AlignmentInterface< biol::AABase> > first_alignment( ALIGNMENTS.FirstElement());
      util::ShPtrList< align::AlignmentInterface< biol::AABase> > alignments_except_first( ++ALIGNMENTS.Begin(), ALIGNMENTS.End());

      // calculate roc curve from gold standard alignment calculated from given sequences
      align::AlignerWordbased< biol::AABase> aligner_wordbased( SCORE, 2);
      storage::List< align::AlignmentHit< biol::AABase> >
        hits( aligner_wordbased.GenerateHits( first_alignment, alignments_except_first));
      math::ROCCurve roc_curve( TestAlignmentHits( calculated_gold_standard_alignment_score_pair.First(), hits, SCORE));
      double integral( roc_curve.Integral());

      const std::string dat_filename( plot_filename + "_calc_gold");
      io::OFStream write_dat;
      io::File::MustOpenOFStream( write_dat, dat_filename + ".dat");
      roc_curve.WriteRatePlottingTable( write_dat);
      io::File::CloseClearFStream( write_dat);

      write_plot << "plot \\\n";
      write_plot << "\"" << dat_filename << ".dat\" using 1:2 with line linestyle 1 "
                 << "title \"Hits from DynamicProgramming alignment: rocint=" << integral << "\"";

      // calculate roc curve from read gold standard alignment
      if( m_FlagGoldStandardAlignment->GetFirstParameter()->GetWasSetInCommandLine())
      {
        const std::string gold_standard_alignment_filename( m_FlagGoldStandardAlignment->GetFirstParameter()->GetValue());
        BCL_MessageStd( "Read gold standard alignment: " + gold_standard_alignment_filename);

        // read in the file and compare
        io::IFStream readstream;
        io::File::MustOpenIFStream( readstream, gold_standard_alignment_filename);
        align::AlignmentNode< biol::AABase> read_gold_standard_alignment( *pir_handler.ReadAlignment( readstream, SEQUENCES));
        io::File::CloseClearFStream( readstream);

        pir_handler.WriteAlignment( util::GetLogger(), read_gold_standard_alignment);

        // test hits, calculate roc curve
        math::ROCCurve roc_curve( TestAlignmentHits( read_gold_standard_alignment, hits, SCORE));
        double integral( roc_curve.Integral());

        const std::string dat_filename( plot_filename + "_read_gold");
        io::OFStream write_dat;
        io::File::MustOpenOFStream( write_dat, dat_filename + ".dat");
        roc_curve.WriteRatePlottingTable( write_dat);
        io::File::CloseClearFStream( write_dat);

        write_plot << ", \\\n"
                   << "\"" << dat_filename << ".dat\" using 1:2 with line linestyle 2 "
                   << "title \"Hits from SABmark alignment: rocint=" << integral << "\"\n";
      }

      write_plot << '\n';
      io::File::CloseClearFStream( write_plot);
    }

    //! @brief test the given list of AlignmentHits against a gold standard alignment
    //! @param ALIGNMENT_GOLD_STANDARD the gold standard alignment to compare the hits with
    //! @param HITS the list of hit alignments
    //! @param ASSIGN_SCORE the assignment scoring ffunction
    //! @return a list of double+bool encoding the score of a hit and if its true positive
    storage::List< storage::Pair< double, bool> >
    GenerateWordbasedAlignerStatistics::TestAlignmentHits
    (
      const align::AlignmentInterface< biol::AABase> &ALIGNMENT_GOLD_STANDARD,
      const storage::List< align::AlignmentHit< biol::AABase> > &HITS,
      const score::AssignmentWithGap< biol::AABase> &ASSIGN_SCORE
    ) const
    {
      // storage for the test and results
      storage::List< storage::Pair< double, bool> > test_results_classified;
      size_t number_true( 0);
      bool is_sub_alignment;
      align::HandlerPIR< biol::AABase> handler_pir;

      // iterate over all given hit alignments
      for
      (
        storage::List< align::AlignmentHit< biol::AABase> >::const_iterator
          itr( HITS.Begin()),
          itr_end( HITS.End());
        itr != itr_end;
        ++itr
      )
      {
        double score( itr->Score( ASSIGN_SCORE));
        is_sub_alignment = ALIGNMENT_GOLD_STANDARD.IsSubAlignment( *itr);
        if( is_sub_alignment)
        {
          number_true++;
        }

        // print general information; print alignment only when message level at least verbose
        if( util::GetMessenger().GetCurrentMessageLevel() > util::Message::e_Standard)
        {
          util::GetLogger() << "Alignment"
            << ": score=" << score
            << "; is_sub_alignment=" << is_sub_alignment
            << "; true=" << number_true << "/" << ( test_results_classified.GetSize() + 1) << '\n';
          handler_pir.WriteAlignment( util::GetLogger(), *itr);
        }

        // store negative of score
        test_results_classified.PushBack( storage::Pair< double, bool>( -score, is_sub_alignment));
      }

      test_results_classified.Sort( &LessThanPairFirst); // sort

      return test_results_classified;
    }

    //! @brief print information about the read sequences to the current logger
    //! @param FILENAMES the files read
    //! @param SEQUENCES the sequences read from the files
    void GenerateWordbasedAlignerStatistics::PrintSequenceData
    (
      const storage::Vector< std::string> &FILENAMES,
      const util::ShPtrVector< biol::AASequence> &SEQUENCES
    ) const
    {
      // just write the filenames, not the sequences if MessageLevel is e_Standard or lower
      if( util::GetMessenger().GetCurrentMessageLevel() <= util::Message::e_Standard)
      {
        std::copy( FILENAMES.Begin(), FILENAMES.End(), std::ostream_iterator< std::string>( util::GetLogger(), "\n"));
        return;
      }

      util::ShPtrVector< biol::AASequence>::const_iterator itr_seq( SEQUENCES.Begin()), itr_seq_end( SEQUENCES.End());
      storage::Vector< std::string>::const_iterator itr_name( FILENAMES.Begin()), itr_name_end( FILENAMES.End());
      for( ; itr_name != itr_name_end && itr_seq != itr_seq_end; ++itr_name, ++itr_seq)
      {
        util::GetLogger() << *itr_name << '\n';
        ( **itr_seq).WriteFasta( util::GetLogger(), 100);
      }
    }

    //! @brief adds the given pair score with the given weight to the scoring function
    //! @param SCORE the pair score to add
    //! @param SCORE_WEIGHT the weight for the pair score
    //! @param SCORE_FCT the scoring function (will be modified)
    //! @return if score was added successfully
    bool GenerateWordbasedAlignerStatistics::AddScore
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

    //! @brief read secondary structure prediction
    //! @param SEQUENCE the AASequence to read/modify
    //! @param CODE the filename
    //! @param SCORE the pair score to determine the file extension
    //! @param METHOD the ss prediction method to use
    //! @return the AASequence passed as SEQUENCE
    biol::AASequence &GenerateWordbasedAlignerStatistics::ReadSSPrediction
    (
      biol::AASequence &SEQUENCE,
      const std::string &CODE,
      const s_PairScores SCORE,
      const sspred::Method METHOD
    ) const
    {
      io::IFStream read;

      io::File::MustOpenIFStream( read, CODE + s_PairScoresFileExtensions[ SCORE]);
      sspred::MethodHandler::ReadPredictionsForAASequence( read, SEQUENCE, METHOD);
      io::File::CloseClearFStream( read);

      return SEQUENCE;
    }

    //! parameters to convert PairScores to ZScores
    const storage::Pair< double, double> GenerateWordbasedAlignerStatistics::s_PairScoresToZScoreConversion[ s_NumberPairScores] =
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
      storage::Pair< double, double>(  -5.3870,   4.7486), // -tfe_engelman
    };

    //! option strings pair scores
    const std::string GenerateWordbasedAlignerStatistics::s_PairScoresFileExtensions[ s_NumberPairScores] =
    {
      "", "", "", "", "", "", "", "", "", "", "", "", "",
      ".ascii", ".psipred_ss2", ".jufo", ".rdb6Prof", ".tmhmm", ".tmmod", ".tmpdb", ".proftmb", ".conpred",
      "", "", "", "", "", "", ""
    };

    //! static instance
    const ApplicationType GenerateWordbasedAlignerStatistics::GenerateWordbasedAlignerStatistics_Instance
    (
      GetAppGroups().AddAppToGroup( new GenerateWordbasedAlignerStatistics(), GetAppGroups().e_InternalBiol)
    );

  } // namespace app
} // namespace bcl
