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
#include "align/bcl_align_aligner_shift.h"
#include "align/bcl_align_handler_classes.h"
#include "assemble/bcl_assemble_biomolecule.h"
#include "assemble/bcl_assemble_domain.h"
#include "assemble/bcl_assemble_locator_aa.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_quality.h"
#include "biol/bcl_biol_dssp.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "fold/bcl_fold_protocol_loop_close.h"
#include "fold/bcl_fold_protocol_loop_coordinate_add.h"
#include "fold/bcl_fold_score_weight_set.h"
#include "function/bcl_function_binary_adapter.h"
#include "function/bcl_function_binary_sum.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_histogram_2d.h"
#include "math/bcl_math_mutate_result.h"
#include "math/bcl_math_running_min_max.h"
#include "math/bcl_math_z_score.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "pdb/bcl_pdb_site.h"
#include "sched/bcl_sched_binary_function_job_with_data.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "score/bcl_score_aa_assignments.h"
#include "score/bcl_score_alignment_assignment.h"
#include "score/bcl_score_alignment_quality.h"
#include "score/bcl_score_assignment.h"
#include "util/bcl_util_wrapper.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FusionProtein
    //! @brief this application seeks to provide a tool for a rational design of chimeric proteins through the fusion of
    //!        proteins by superimposing/aligning SSEs in a scaffold and donor protein and finding the optimal cut point
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_app_fusion_protein.cpp @endlink
    //! @author woetzen
    //! @date Mar 13, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FusionProtein :
      public InterfaceRelease
    {

    private:

    //////////
    // data //
    //////////

      //! input pdb for scaffold
      util::ShPtr< command::FlagStatic> m_ScaffoldProteinFlag;

      //! input pdb for donor
      util::ShPtr< command::FlagDynamic> m_DonorProteinFlag;

      //! select fragments within scaffold as anchor points
      util::ShPtr< command::FlagDynamic> m_ScaffoldFragmentFlag;

      //! select fragment within donor protein to fuse
      util::ShPtr< command::FlagDynamic> m_DonorSitesFlag;

      //! quality measure used for the alignment and the cutoff to write out the file
      util::ShPtr< command::FlagStatic> m_QualityMeasureFlag;
      util::ShPtr< command::Parameter>  m_QualityCutoffParam;

      //! flag to pass all assignment scores that are added and reported for each assignment in the alignment output
      util::ShPtr< command::FlagDynamic> m_AssignmentScoreFlag;

      //! flag to trigger the output of the alignment used for each fusion
      util::ShPtr< command::FlagStatic> m_WriteAlignmentFlag;

      //! flag to trigger the output of the rosetta resfile needed for design
      util::ShPtr< command::FlagStatic> m_WriteRosettaResfileFlag;

      //! flag for the minimal number of aligned residues before the fused protein is written
      util::ShPtr< command::FlagStatic> m_MinNumberAlignedResiduesFlag;

      //! flag to write only selected models
      util::ShPtr< command::FlagDynamic> m_WriteSelectedModelsFlag;

      //! flag to write cross over heatmaps
      util::ShPtr< command::FlagStatic> m_WriteCrossOverHeatmapFlag;

      //! prefix for all output file names
      util::ShPtr< command::FlagStatic> m_OutputPrefixFlag;

      //! flag to enable the appending of the donor instead of the replacement of the region between the anchor
      util::ShPtr< command::FlagStatic> m_AppendFlag;

      //! flag to replace terms on scaffold
      util::ShPtr< command::FlagStatic> m_ReplaceTermFlag;
      util::ShPtr< command::Parameter> m_NTermDonorParam;
      util::ShPtr< command::Parameter> m_NTermScaffoldLocatorParam;
      util::ShPtr< command::Parameter> m_NTermDonorLocatorParam;
      util::ShPtr< command::Parameter> m_CTermDonorParam;
      util::ShPtr< command::Parameter> m_CTermScaffoldLocatorParam;
      util::ShPtr< command::Parameter> m_CTermDonorLocatorParam;

      //! flag to pass file containing math::Vector with weight set to combine score
      util::ShPtr< command::FlagInterface> m_WeightSetFlag;

      //! flag to define wether cutpoints are generated systematically (default) or at the most optimal peptide bond
      util::ShPtr< command::FlagStatic> m_OptimalPeptideBondCutPointFlag;

      //! flag to create models for the combination of all replacement sites
      util::ShPtr< command::FlagStatic> m_ReplaceCombinatorialFlag;

      // output directory
      mutable io::Directory m_OutputDirectory;

      //! weight set
      mutable storage::Table< double> m_Weights;

      //! scores
      mutable score::ProteinModelScoreSum m_ScoreFunction;

      //! results
      mutable storage::Table< double> m_Results;

      //! pdb factory
      pdb::Factory            m_Factory;

      //! pymol script header
      mutable std::string m_PymolScriptHeader;

      //! scaffold protein
      mutable assemble::ProteinModel m_Scaffold;

      //! atom types used for superimposition
      mutable storage::Set< biol::AtomType> m_AtomTypes;

      //! measure used for superimposition
      mutable util::ShPtr< quality::SuperimposeInterface> m_SPMeasure;

      mutable util::ShPtr< function::BinarySum< const biol::AABase, const biol::AABase, double> > m_SPAAAssignmentScore;
      mutable util::ShPtr< score::Assignment< biol::AABase> >                                     m_SPAssignmentScore;
      mutable util::ShPtr< function::UnaryInterface< const align::AlignmentInterface< biol::AABase>, double> >
        m_SPScoreAlignmentAssignment;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FusionProtein();

    public:

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class Result
      //! @brief class that represents the result of a fusion of two or more sequences
      //! @details this class contains the final fused sequence, a pymol script coloring the pdb file adequaetly
      //!          the assignment scores for the aligned resgion and the number of residues that were fused between the
      //!          scaffold and donor anchors
      //!
      //! @author woetzen
      //! @date Feb 01, 2013
      //!
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      class Result :
        public util::ObjectInterface
      {
      public:
        biol::AASequence m_Sequence;
        std::string      m_PymolScriptTop;
        std::string      m_PymolScript;
        std::string      m_PymolTailScript;
        std::string      m_RosettaResFile; //!< rosetta resfile according to https://www.rosettacommons.org/manuals/archive/rosetta3.4_user_guide/d1/d97/resfiles.html
        std::string      m_ReplaceInformation; //!< information about the res ids at the crossover points
        int              m_LastScaffoldSeqID;
        int              m_LastDonorSeqID;
        double           m_AlignFuseScore;
        size_t           m_NrFusedResiduesNTerm;
        size_t           m_NrFusedResiduesCTerm;
        size_t           m_NrReplacedResidues;
        std::string      m_FilenameBase;

        //! @brief default constructor
        Result() :
          m_Sequence(),
          m_PymolScriptTop(),
          m_PymolScript(),
          m_PymolTailScript(),
          m_RosettaResFile(),
          m_ReplaceInformation(),
          m_LastScaffoldSeqID( -util::GetUndefined< int>()),
          m_LastDonorSeqID( -util::GetUndefined< int>()),
          m_AlignFuseScore( 0.0),
          m_NrFusedResiduesNTerm( 0),
          m_NrFusedResiduesCTerm( 0),
          m_NrReplacedResidues( 0),
          m_FilenameBase()
        {
        }

        Result *Clone() const
        {
          return new Result( *this);
        }

        const std::string &GetClassIdentifier() const
        {
          return GetStaticClassName< Result>();
        }

        void Append( const Result &RESULT)
        {
          m_Sequence.AppendSequence( RESULT.m_Sequence);
          m_PymolScriptTop       += RESULT.m_PymolScriptTop;
          m_PymolScript          += RESULT.m_PymolScript;
          m_PymolTailScript      += RESULT.m_PymolTailScript;
          m_RosettaResFile       += RESULT.m_RosettaResFile;
          m_ReplaceInformation   += RESULT.m_ReplaceInformation;
          m_LastScaffoldSeqID     = std::max( m_LastScaffoldSeqID, RESULT.m_LastScaffoldSeqID);
          m_LastDonorSeqID        = std::max( m_LastDonorSeqID   , RESULT.m_LastDonorSeqID);
          m_AlignFuseScore       += RESULT.m_AlignFuseScore;
          m_NrFusedResiduesNTerm += RESULT.m_NrFusedResiduesNTerm;
          m_NrFusedResiduesCTerm += RESULT.m_NrFusedResiduesCTerm;
          m_NrReplacedResidues   += RESULT.m_NrReplacedResidues;
        }

        //! @brief read from std::istream
        //! @param ISTREAM input stream
        //! @return istream which was read from
        std::istream &Read( std::istream &ISTREAM)
        {
          return ISTREAM;
        }

        //! @brief write to std::ostream
        //! @param OSTREAM output stream to write to
        //! @param INDENT - number of indentations
        //! @return output stream which was written to
        std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
        {
          return OSTREAM;
        }

      }; // class Result

      //! @brief Clone function
      //! @return pointer to new FusionProtein
      FusionProtein *Clone() const
      {
        return new FusionProtein( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief return the bcl::commons name
      //! @return string for the bcl::commons name of that application
      std::string GetBCLScopedName() const
      {
        return "BCL::FusionProtein";
      }

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const
      {
        // create a static string to hold readme information
        static const std::string s_readme_text
        (
          DefaultSectionSeparator() +
          "I. OVERVIEW.\n"
          "This document provides a description of BCL::FusionProtein, terms of use, appropriate citation,"
          " installation procedures, BCL::FusionProtein execution, technical support, and future research directions.\n"
          "\n"
          + DefaultSectionSeparator() +
          "II. WHAT IS BCL::FusionProtein?\n"
          "BCL::FusionProtein is a C++ based application, created by Max-Planck-Institute of Developmental Biology's "
          "Hoecker group, which is part of a larger library of applications called BCL::Commons.\n"
          "BCL::FusionProtein fuses fragments of donor proteins to a scaffold protein to create a chimeric protein. "
          "It superimposes the pair of secondary structure elements flanking the donor's fragment to identically typed "
          "secondary structure elements in the scaffold protein and detects an optimal cross over site between the two "
          "protein's sequences.\n"
          "\n"
          + DefaultSectionSeparator() +
          "III. TERMS OF USE.\n"
          + DefaultTermsOfUseString() +
          "\n"
          + DefaultSectionSeparator() +
          "IV. APPROPRIATE CITATIONS FOR USING BCL::FusionProtein.\n"
          "When using BCL::FusionProtein in a publication, please cite the publications describing the application's "
          "development, which is currently under review. Any news will be published at "
          "http://www.eb.tuebingen.mpg.de/research/research-groups/birte-hoecker.html\n"
          "\n"
          + DefaultSectionSeparator() +
          "V. INSTALLATION PROCEDURES.\n"
          + DefaultInstallationProcedure() +
          "\n"
          + DefaultSectionSeparator() +
          "VI. RUNNING BCL::FusionProtein.\n"
          "Running BCL::FusionProtein consists of three main steps.\n"
          "\n"
          "1) Identify a scaffold protein and at least one donor protein that is to be fused to the scaffold protein.\n"
          "These are advised to be single chain files. Use PDB::Convert to extract single chains and if desired to "
          "correct secondary structure assignments (with the \"-dssp\" flag).\n"
          "\n"
          "2) Define the \"active site\" withthin the donor protein.\n"
          "A site definition file needs to be provided that looks like this at defines all residues that need to be "
          "present in the fused chimeric protein:\n"
          "file: 1thf_loop1.site:\n"
          "bcl::pdb::Site\n"
          "AC1\n"
          "SOFTWARE\n"
          "BINDING loop1\n"
          "bcl::storage::List<bcl::pdb::ResidueSimple>\n"
          "  2\n"
          "  bcl::pdb::ResidueSimple\n"
          "    \"LYS\"     D       13      ' '\n"
          "  bcl::pdb::ResidueSimple\n"
          "    \"ASP\"     D       31      ' '\n"
          "bcl::storage::List<bcl::pdb::ResidueSimple>\n"
          "  0\n"
          "bcl::util::ShPtr<bcl::pdb::Ligand>\n"
          "  0\n"
          "  NULL\n\n"
          "The algorithm will automatically use the closest secondary structure elements (helix or strand) which do "
          "not contain any of the site residues.\n"
          "\n"
          "3) Define scaffold anchor point locators.\n"
          "Locators define the anchor points within the scaffold protein the \"active site\" is fused to. If these are "
          "not given, all possible pairs of secondary structure elements in the scaffold protein will be tested as "
          "potential fusion sites. The locator file is defined like this:\n"
          "file: 1qo2_loop1.locators:\n"
          "bcl::assemble::LocatorAA\n"
          "  LocatorAA(locator_chain=A,seq_id=10,use_pdb_id=1)\n"
          "bcl::assemble::LocatorAA\n"
          "  LocatorAA(locator_chain=A,seq_id=31,use_pdb_id=1)\n\n"
          "This file defines any residue that is contained in the left and right anhor secondary structure element "
          "that is to be used as cross over region.\n"
          "\n"
          "3) Run BCL::FusionProtein:\n"
          "Please read the publication on the different parameters and there meaning; also call BCL::FusionProtein on "
          "the command line with \"-help\" which documents the meaning of each flag. A sample command line could be:\n"
          "\n"
          "protein:FusionProtein "
          "-scaffold 1qo2.pdb "
          "-scaffold_fragment 1qo2_loop1.locators "
          "-donor 1thf.pdb "
          "-donor_sites 1thf_loop1.site "
          "-helix_classes 1 5 "
          "-merge_overlapping_sses "
          "-assignment_score steric hydrophobicity "
          "-cutpoint_optimal_peptide "
          "-min_number_aligned_residues 3 "
          "-quality RMSD 4 "
          "-output_prefix 1thf_1qo2/ "
          "-histogram_path histogram/ "
          "\n"
          "FLAGS:\n"
          "\n"
          "-scaffold 1qo2.pdb                      -> the scaffold protein\n"
          "-scaffold_fragment 1qo2_loop1.locators  -> the locator defining the anchor SSEs in the scaffold protein\n"
          "-donor 1thf.pdb                         -> the donor protein\n"
          "-donor_sites 1thf_loop1.site            -> the residues defining the active site in the donor\n"
          "-helix_classes 1 5                      -> pdb helix classes considered"
          "-merge_overlapping_sses                 -> merge secondary structure elements in given pdbs if their definition overlaps\n"
          "-assignment_score steric hydrophobicity -> assignment scores considered\n"
          "-cutpoint_optimal_peptide               -> choose the sequence cross over so that the resulting peptide bond is optimal\n"
          "-min_number_aligned_residues 3          -> after superimposing donor and scaffold SSE at least 3 residues have to align\n"
          "-quality RMSD 4                         -> the RMSD of the superimposed SSEs has to be at least 4 Angstrom\n"
          "-output_prefix 1thf_1qo2/               -> path to write results to\n"
          "-histogram_path histogram/              -> bcl histograms necessary to initialize the scoring terms\n"
          "\n"
          "INPUT AND OUTPUT.\n"
          "\n"
          "BCL::FusionProtein usually uses at least 4 input files, two pdb files for scaffold and donor, and a site "
          "definition file as well as a fragment locator. Theses are documented above.\n"
          "The output files are:\n"
          "donor1D.fasta       -> fasta sequence of first donor\n"
          "donor1_fusedA.fasta -> fasta sequence of the protein after fusing the first donor\n"
          "donor1_fused.pdb    -> model file of the protein after fusing the first donor\n"
          "donor1_fused.pml    -> pymol session to visualize the fused protein, cross over points and other details\n"
          "donor1.pdb          -> pdb file of the donor (necessary for the pymol session)\n"
          "results.table       -> table with detailed information and scores of all structures\n"
          "scaffoldA.fasta     -> fasta sequence of the scaffold\n"
          "scaffold.pdb        -> pdb file of the scaffold (necessary for the pymol session)\n"
          "\n"
          + DefaultSectionSeparator() +
          "VII. TECHNICAL SUPPORT.\n"
          + DefaultTechnicalSupportString() +
          "\n"
          + DefaultSectionSeparator() +
          "VIII. FUTURE DEVELOPMENT OF BCL::FusionProtein.\n"
          "BCL::FusionProtein is under ongoing further development. For current research please refer to "
          "http://www.eb.tuebingen.mpg.de/research/research-groups/birte-hoecker.html\n"
          "http://www.meilerlab.org\n"
          + DefaultSectionSeparator()
        );

        // return readme information
        return s_readme_text;
      }

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const
      {
        return "designs chimeric proteins from a scaffold by fusing mutliple donor protein fragments";
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const
      {
        util::ShPtr< command::Command> sp_command( new command::Command());

        sp_command->AddFlag( m_ScaffoldProteinFlag);
        sp_command->AddFlag( m_DonorProteinFlag);
        sp_command->AddFlag( m_ScaffoldFragmentFlag);
        sp_command->AddFlag( m_DonorSitesFlag);

        sp_command->AddFlag( m_QualityMeasureFlag);
        sp_command->AddFlag( m_AssignmentScoreFlag);
        sp_command->AddFlag( m_WriteAlignmentFlag);
        sp_command->AddFlag( m_WriteRosettaResfileFlag);
        sp_command->AddFlag( m_MinNumberAlignedResiduesFlag);
        sp_command->AddFlag( m_OutputPrefixFlag);
        sp_command->AddFlag( m_WriteSelectedModelsFlag);
        sp_command->AddFlag( m_WriteCrossOverHeatmapFlag);
        sp_command->AddFlag( m_AppendFlag);
        sp_command->AddFlag( m_ReplaceTermFlag);
        sp_command->AddFlag( m_WeightSetFlag);
        sp_command->AddFlag( m_OptimalPeptideBondCutPointFlag);
        sp_command->AddFlag( m_ReplaceCombinatorialFlag);

        sp_command->AddFlag( pdb::Handler::GetFlagHelixClasses());
        sp_command->AddFlag( pdb::Handler::GetFlagMergeOverlappingSSEs());

        sp_command->AddFlag( pdb::Factory::GetFlagBiomolecule());

        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_command);

        return sp_command;
      }

      //! Main
      int Main() const
      {
        fold::ProtocolLoopCoordinateAdd::GetInstance().InitializeScores();

        // initialize weight set
        fold::ScoreWeightSet score_weight_set;

        // set weights from file
        if( m_WeightSetFlag->GetFirstParameter()->GetWasSetInCommandLine())
        {
          // initialize read
          BCL_MessageStd( "reading score weightset from file");
          io::IFStream read;
          io::File::MustOpenIFStream( read, m_WeightSetFlag->GetFirstParameter()->GetValue());
          m_Weights.ReadFormatted( read);
          // close stream
          io::File::CloseClearFStream( read);
          score_weight_set = fold::ScoreWeightSet( m_Weights);
        }
        // from the names in the scoring functions map
        else
        {
          for
          (
            util::Enumerate< util::ShPtr< score::ProteinModel>, fold::Scores>::const_iterator
              itr( fold::GetScores().Begin()), itr_end( fold::GetScores().End());
            itr != itr_end; ++itr
          )
          {
            score_weight_set.SetWeight( *itr, 1.0);
          }

          // create weights tables with weight one for each of the scores
          m_Weights = score_weight_set.CreateTable();
        }

        // construct the scoring function
        m_ScoreFunction = *( score_weight_set.ConstructScoreSum());

        // table for all results for each model generated
        storage::Vector< std::string> table_header;
        table_header.PushBack( "quality");
        table_header.PushBack( "align");
        table_header.PushBack( "align_fuse");
        table_header.PushBack( "align_norm");
        table_header.PushBack( "align_fuse_norm");
        table_header.PushBack( "nr_coord_pairs");
        table_header.PushBack( "nr_fused_pairs_first");
        table_header.PushBack( "nr_fused_pairs_last");
        table_header.PushBack( "nr_fused_pairs");
        table_header.PushBack( "nr_replaced");
        table_header.PushBack( "crossover_scaffold");
        table_header.PushBack( "crossover_donor");
        table_header.PushBack( "length");
        table_header.PushBack( "file");
        table_header.Append( m_ScoreFunction.GetFunctionSchemes());
        m_Results = storage::Table< double>( ( storage::TableHeader( table_header)));

        const storage::Map< biol::SSType, size_t> min_sse_size( pdb::Factory::GetCommandlineSSETypeMinSizes());

        // initialize the assignment score
        m_SPAAAssignmentScore = util::ShPtr< function::BinarySum< const biol::AABase, const biol::AABase, double> >
        (
          AssignmentScore()
        );

        m_SPAssignmentScore = util::ShPtr< score::Assignment< biol::AABase> >
        (
          new score::Assignment< biol::AABase>( m_SPAAAssignmentScore)
        );

        // score for alignment
        m_SPScoreAlignmentAssignment = util::ShPtr< function::UnaryInterface< const align::AlignmentInterface< biol::AABase>, double> >
          (
            new score::AlignmentAssignment< biol::AABase>( m_SPAssignmentScore)
          );

        m_AtomTypes = biol::GetAtomTypes().CA;
        m_SPMeasure = quality::Measure( m_QualityMeasureFlag->GetFirstParameter()->GetValue())->HardCopy();

        // initialize the pymol script header
        m_PymolScriptHeader = PymolScriptHeader();

        // initialize the scaffold protein
        const std::string scaffold_pdb_id( InitializeScaffoldProtein());

        m_OutputDirectory = io::Directory
        (
          io::File::SplitToPathAndFileName( m_OutputPrefixFlag->GetFirstParameter()->GetValue()).First()
        );

        // check if c and n term should be replaced
        if( m_ReplaceTermFlag->GetFlag())
        {
          return ReplaceTerm( scaffold_pdb_id);
        }

        // read in all donors
        storage::List< storage::Pair< assemble::ProteinModel, pdb::Site> > donors( Donors());
        BCL_MessageStd( "using that many donors: " + util::Format()( donors.GetSize()));

        // add sites to donors
        if( m_DonorSitesFlag->GetFlag())
        {
          AddSitesToDonors( donors);
        }

        // get locators for anchor, if given in commandline
        storage::Vector< storage::VectorND< 2, util::ShPtr< assemble::LocatorAA> > > scaffold_anchor_locators;
        if( m_ScaffoldFragmentFlag->GetFlag())
        {
          const storage::Vector< std::string> locator_file_names( m_ScaffoldFragmentFlag->GetStringList());
          BCL_Assert( donors.GetSize() == locator_file_names.GetSize(), "need as many locator definitions as donors");

          for
          (
            storage::Vector< std::string>::const_iterator
              loc_itr( locator_file_names.Begin()), loc_itr_end( locator_file_names.End());
            loc_itr != loc_itr_end;
            ++loc_itr
          )
          {
            io::IFStream read;
            io::File::MustOpenIFStream( read, *loc_itr);
            util::ShPtr< assemble::LocatorAA> sp_locator_anchor_first( new assemble::LocatorAA());
            util::ShPtr< assemble::LocatorAA> sp_locator_anchor_last( new assemble::LocatorAA());
            read >> *sp_locator_anchor_first;
            read >> *sp_locator_anchor_last;
            BCL_Assert
            (
                 sp_locator_anchor_first->GetLocatorChain().GetChainID()
              == sp_locator_anchor_last->GetLocatorChain().GetChainID(),
              "given locators for scaffold are from different chains"
            );
            io::File::CloseClearFStream( read);
            BCL_MessageStd
            (
              "selecting fragment to replace:\n" + util::Format()( *sp_locator_anchor_first) +
              "\nand\n" + util::Format()( *sp_locator_anchor_last)
            );
            scaffold_anchor_locators.PushBack
            (
              storage::VectorND< 2, util::ShPtr< assemble::LocatorAA> >
              (
                sp_locator_anchor_first, sp_locator_anchor_last
              )
            );
          }
        }

        storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, util::SiPtr< const assemble::ProteinModel> > > alignments
        (
          CreateAlignments( donors, scaffold_pdb_id)
        );

        BCL_MessageStd
        (
          "creating fused proteins from alignments: " + util::Format()( alignments.GetSize())
        );

        storage::List< std::string> donor_filenames;
        size_t count( 1);
        double align_score_sum( 0);
        size_t nr_coord_pairs_sum( 0);

        storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > > heatmap_boxes;
        std::string heatmap_title( scaffold_pdb_id);
        std::string heatmap_xaxis( "donor crossover");

        // iterate through alignments
        for
        (
          storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, util::SiPtr< const assemble::ProteinModel> > >::iterator
            align_itr( alignments.Begin()), align_itr_end( alignments.End());
          align_itr != align_itr_end;
          ++align_itr, ++count
        )
        {
          const std::string donor_name( "donor" + util::Format()( count));
          const align::AlignmentNode< biol::AABase> &alignment( align_itr->First());
          assemble::ProteinModel donor_protein( *util::ShPtr< assemble::ProteinModel>( align_itr->Second()->HardCopy()));
          const util::ShPtr< util::Wrapper< std::string> > sp_donor_id
          (
            donor_protein.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Identification)
          );
          BCL_Assert( sp_donor_id.IsDefined(), "cannot find pdb id for donor");
          const std::string donor_pdb_id( sp_donor_id->GetData());
          heatmap_title += "->" + donor_pdb_id + "->" + scaffold_pdb_id;
          heatmap_xaxis += ' ' + donor_pdb_id;

          // alignment has scaffold on top and donor on bottom
          const char scaffold_chain_id( alignment.GetSequences().FirstElement()->GetFirstMember()->GetChainID());
          const char donor_chain_id( ( *( ++alignment.GetSequences().Begin()))->GetFirstMember()->GetChainID());

          // write donor fasta and pdb
          const std::string file_name_base( m_OutputPrefixFlag->GetFirstParameter()->GetValue() + donor_name);

          io::OFStream write;
          // write donor and scaffold information
          const std::string donor_fasta_filename( file_name_base + donor_chain_id + ".fasta");
          if( io::File::TryOpenOFStream( write, donor_fasta_filename))
          {
            donor_protein.GetChain( donor_chain_id)->GetSequence()->WriteFasta( write);
            io::File::CloseClearFStream( write);
          }
          else
          {
            BCL_MessageCrt( "unable to open file for writing: " + donor_fasta_filename);
          }

          const io::DirectoryEntry scaffold_fasta_filename
          (
            m_OutputPrefixFlag->GetFirstParameter()->GetValue() + "scaffold" + scaffold_chain_id + ".fasta"
          );
          if( !scaffold_fasta_filename.DoesExist() && io::File::TryOpenOFStream( write, scaffold_fasta_filename.GetFullName()))
          {
            m_Scaffold.GetChain( scaffold_chain_id)->GetSequence()->WriteFasta( write);
            io::File::CloseClearFStream( write);
          }
          else
          {
            BCL_MessageCrt( "unable to open file for writing: " + scaffold_fasta_filename.GetFullName());
          }

          // score of the alignment
          const double align_score( m_SPScoreAlignmentAssignment->operator ()( alignment));
          align_score_sum += align_score;

          const storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > coord_pair
          (
            assemble::Quality::CoordinatesFromAlignment( alignment, m_AtomTypes)
          );
          const size_t nr_coord_pairs( coord_pair.First().GetSize());
          nr_coord_pairs_sum += nr_coord_pairs;
          BCL_MessageStd( "number coord pairs: " + util::Format()( nr_coord_pairs));

          const double quality( m_SPMeasure->CalculateMeasure( coord_pair.First(), coord_pair.Second()));
          const math::TransformationMatrix3D transformation
          (
            m_SPMeasure->CalculateSuperimposition( coord_pair.Second(), coord_pair.First())
          );

          const std::string donor_filename( file_name_base + ".pdb");
          if( io::File::TryOpenOFStream( write, donor_filename))
          {
            donor_filenames.PushBack( donor_filename);
            // transform sequences within donor as well
            donor_protein.Transform( transformation);
            m_Factory.WriteModelToPDB( donor_protein, write);
            io::File::CloseClearFStream( write);
          }
          else
          {
            BCL_MessageCrt( "unable to open file for writing: " + donor_filename);
          }

          // add row for donor protein to score table
          // reassign secondary structure for scoring
          {
            // actual scores for scaffold
            storage::Row< double> &row( m_Results.InsertRow( donor_name));
            biol::DSSP dssp;
            math::MutateResult< assemble::ProteinModel> result( dssp( donor_protein));
            if( result.GetArgument().IsDefined())
            {
              // actual scores
              AddScoresToRow( *result.GetArgument(), Result(), row);
            }
            else
            {
              // actual scores
              AddScoresToRow( donor_protein, Result(), row);
            }
          }

          // write alignment
          // get align::Handler from command line parameter
          util::ShPtr< align::HandlerInterface< biol::AABase> > handler
          (
            align::HandlerClasses< biol::AABase>::HandlerClass
            (
              m_WriteAlignmentFlag->GetFirstParameter()->GetValue()
            )->HardCopy()
          );
          handler->SetAssignmentScore( util::CloneToShPtr( score::Assignment< biol::AABase>( m_SPAAAssignmentScore)));
          if( m_WriteAlignmentFlag->GetFlag())
          {
            const std::string file_name_align( file_name_base + "alignment" + handler->GetFileExtension());
            if( io::File::TryOpenOFStream( write, file_name_align))
            {
              handler->WriteAlignment( write, alignment);
              io::File::CloseClearFStream( write);
              BCL_MessageVrb( "wrote alignment to file: " + file_name_align);
            }
            else
            {
              BCL_MessageCrt( "unable to open file for writing: " + file_name_align);
            }
          }
          else if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Verbose))
          {
            handler->WriteAlignment( util::GetLogger(), alignment);
          }

          storage::List< Result> fused_sequences;

          if( m_AppendFlag->GetFlag())
          {
            fused_sequences.PushBack
            (
              Fuse
              (
                *m_Scaffold.GetChain( scaffold_chain_id)->GetSequence(),
                *donor_protein.GetChain( donor_chain_id)->GetSequence(),
                transformation,
                alignment
              )
            );
          }
          else
          {
            storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, util::SiPtr< const assemble::ProteinModel> > > single_alignment;
            single_alignment.PushBack( *align_itr);

            if( m_OptimalPeptideBondCutPointFlag->GetFlag())
            {
              fused_sequences.PushBack
              (
                Replace
                (
                  *m_Scaffold.GetChain( scaffold_chain_id)->GetSequence(),
                  single_alignment
                )
              );
            }
            else
            {
              fused_sequences = ReplaceSystematic
                                (
                                  *m_Scaffold.GetChain( scaffold_chain_id)->GetSequence(),
                                  single_alignment
                                );
            }
          }

          const util::Format count_format( util::Format().W( size_t( std::log10( fused_sequences.GetSize()) + 1)).ForceW().Fill( '0').R());
          const bool write_selected_models( m_WriteSelectedModelsFlag->GetFlag());
          storage::Set< size_t> selected_models;
          {
            const storage::Vector< size_t> temp( m_WriteSelectedModelsFlag->GetNumericalList< size_t>());
            selected_models.InsertElements( temp.Begin(), temp.End());
          }

          util::ShPtrList< sched::JobInterface> jobs;

          size_t count( 0);
          for
          (
            storage::List< Result>::iterator itr( fused_sequences.Begin()), itr_end( fused_sequences.End());
            itr != itr_end;
            ++itr, ++count
          )
          {
            Result &fused_sequence( *itr);
            fused_sequence.m_Sequence.SetFastaHeader( fused_sequence.m_ReplaceInformation);
            fused_sequence.m_FilenameBase = file_name_base + ( m_OptimalPeptideBondCutPointFlag->GetFlag() ? "" : '_' + count_format( count));
            fused_sequence.m_PymolScriptTop += "load " + donor_filename + ", " + "donor1_model\n";
            fused_sequence.m_PymolScriptTop += "hide all\n";
            fused_sequence.m_PymolScriptTop += "cd " + m_OutputDirectory.GetPath() + '\n';
            fused_sequence.m_PymolScriptTop += "\n#END HEADER\n";

            // write result to table
            storage::Row< double> &row( m_Results.InsertRow( fused_sequence.m_FilenameBase));
            row[ "quality"] = quality;
            row[ "align"] = align_score;
            row[ "align_norm"] = align_score / nr_coord_pairs;
            row[ "nr_coord_pairs"] = nr_coord_pairs;
            row[ "crossover_scaffold"] = fused_sequence.m_LastScaffoldSeqID;
            row[ "crossover_donor"] = fused_sequence.m_LastDonorSeqID;
            row[ "length"] = fused_sequence.m_Sequence.GetSize();

            // write a file for the fused model
            row[ "file"] = size_t( !write_selected_models || selected_models.Erase( count) != 0);

            // ProcessFusedSequence( fused_sequence, row);
            util::ShPtr< sched::JobInterface> sp_job
            (
              new sched::BinaryFunctionJobWithData< const Result, storage::Row< double>, void, FusionProtein>
              (
                0,
                *this,
                &FusionProtein::ProcessFusedSequence,
                fused_sequence,
                row,
                sched::JobInterface::e_READY,
                NULL
              )
            );

            sched::GetScheduler().SubmitJob( sp_job);
            jobs.PushBack( sp_job);
          }

          // iterate through all results and join
          while( jobs.GetSize() != 0)
          {
            util::ShPtr< sched::JobInterface> &sp_job( jobs.FirstElement());
            sched::GetScheduler().Join( sp_job);

            util::ShPtr< sched::BinaryFunctionJobWithData< const Result, storage::Row< double>, void, FusionProtein> > sp_job_with_data( sp_job);
            if( !sp_job.IsDefined())
            {
              continue;
            }
            jobs.PopFront();

//            const storage::Row< double> &row( sp_job_with_data->GetSecondArgument());
//            const Result &fused_sequence( sp_job_with_data->GetFirstArgument());
//            const int last_scaffold_seqid( row[ "crossover_scaffold"]);
//            const int last_donor_seqid( row[ "crossover_donor"]);
//            BCL_MessageStd
//            (
//              "heatmap field: " +
//              util::Format()( last_scaffold_seqid) + ", " +
//              util::Format()( last_donor_seqid) + ", " +
//              fused_sequence.m_ReplaceInformation + ", " +
//              util::Format()( row[ "sum"])
//            );
//            score_histogram_crossover.PushBack
//            (
//              storage::VectorND< 2, double>( last_scaffold_seqid, last_donor_seqid),
//              row[ "sum"] / fused_sequence.m_Sequence.GetSize()
//            );
          }

          // write histogram to heatmap
          if( !m_OptimalPeptideBondCutPointFlag->GetFlag())
          {
            // find the optimal cutpoint
            storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, util::SiPtr< const assemble::ProteinModel> > > single_alignment;
            single_alignment.PushBack( *align_itr);

            const Result optimal_fused_sequence
            (
              Replace
              (
                *m_Scaffold.GetChain( scaffold_chain_id)->GetSequence(),
                single_alignment
              )
            );

            heatmap_boxes.PushBack
            (
              storage::VectorND< 2, storage::VectorND< 2, double> >
              (
                storage::VectorND< 2, double>( optimal_fused_sequence.m_LastScaffoldSeqID - 1.5, optimal_fused_sequence.m_LastDonorSeqID - 1.5),
                storage::VectorND< 2, double>( optimal_fused_sequence.m_LastScaffoldSeqID - 0.5, optimal_fused_sequence.m_LastDonorSeqID - 0.5)
              )
            );
          }
        }

        // write heatmaps
        if( m_WriteCrossOverHeatmapFlag->GetFlag())
        {
          const math::Histogram2D histogram_template( CreateCrossoverHistogram());
          storage::Map< std::string, math::Histogram2D> scores_histogram_crossover;
          storage::Vector< std::string> schemes( m_ScoreFunction.GetFunctionSchemes());
          schemes.PushBack( "quality");
          schemes.PushBack( "align_fuse_norm");
          schemes.PushBack( "sum");
          for
          (
            storage::Vector< std::string>::const_iterator scheme_itr( schemes.Begin()), scheme_itr_end( schemes.End());
            scheme_itr != scheme_itr_end;
            ++scheme_itr
          )
          {
            scores_histogram_crossover[ *scheme_itr] = histogram_template;
            scores_histogram_crossover[ *scheme_itr + "_seqnorm"] = histogram_template;
          }

          storage::Map< std::string, math::RunningMinMax< double> > scores_min_max;
          for
          (
            storage::Table< double>::const_iterator itr( m_Results.Begin()), itr_end( m_Results.End());
            itr != itr_end;
            ++itr
          )
          {
            const storage::Row< double> &row( itr->Second());
            const int last_scaffold_seqid( row[ "crossover_scaffold"]);
            const int last_donor_seqid( row[ "crossover_donor"]);
            const size_t length( row[ "length"]);

            if( last_scaffold_seqid < 0 || last_donor_seqid < 0 || length == 0)
            {
              continue;
            }

            storage::VectorND< 2, double> bin_coord( last_scaffold_seqid, last_donor_seqid);

            for
            (
              storage::Vector< std::string>::const_iterator scheme_itr( schemes.Begin()), scheme_itr_end( schemes.End());
              scheme_itr != scheme_itr_end;
              ++scheme_itr
            )
            {
              const std::string &current_scheme( *scheme_itr);
              const double current_score( row[ current_scheme]);
              scores_histogram_crossover[ current_scheme].PushBack( bin_coord, current_score);
              scores_min_max[ current_scheme] += current_score;

              // normalized score
              const std::string norm_scheme( current_scheme + "_seqnorm");
              const double score_norm( current_score / length);
              scores_histogram_crossover[ norm_scheme].PushBack( bin_coord, score_norm);
              scores_min_max[ norm_scheme] += score_norm;
            }
          }

          // heatmap of all crossover points
          for
          (
            storage::Map< std::string, math::RunningMinMax< double> >::const_iterator
              min_max_itr( scores_min_max.Begin()), min_max_end( scores_min_max.End());
            min_max_itr != min_max_end;
            ++min_max_itr
          )
          {
            const std::string &current_scheme( min_max_itr->first);
            const math::RunningMinMax< double> &min_max( min_max_itr->second);
            const math::Histogram2D &current_histogram( scores_histogram_crossover.Find( current_scheme)->second);
            io::OFStream write;
            math::GnuplotHeatmap heatmap;

            // if the min max range of the score is empty then there is no meaningful heatmap
            if( min_max.GetRange() < 0.000001)
            {
              continue;
            }

            heatmap.SetFromHistogram( current_histogram, true, true);
            heatmap.SetTitleAndLabel
            (
              heatmap_title + " score: " + current_scheme,
              heatmap_xaxis,
              "scaffold crossover " + scaffold_pdb_id, "score: " + current_scheme
            );

            heatmap.SetBoxes
            (
              heatmap_boxes,
              current_histogram.GetNumberOfBinsX() * current_histogram.GetBinSizeXY().First(),
              current_histogram.GetNumberOfBinsY() * current_histogram.GetBinSizeXY().Second(),
              current_histogram.GetBoundariesX().First(),
              current_histogram.GetBoundariesY().First()
            );
            heatmap.SetMinMaxZ( min_max.GetMin(), min_max.GetMax());
            heatmap.SetPalette( ( min_max.GetMax() <= 0.0) ? math::GnuplotHeatmap::e_BlueGreenRedWhite : math::GnuplotHeatmap::e_WhiteBlueGreenRed);

            heatmap.SetFilename( m_OutputPrefixFlag->GetFirstParameter()->GetValue() + "score_heatmap_" + current_scheme);
            const std::string filename_heatmap( m_OutputPrefixFlag->GetFirstParameter()->GetValue() + "score_heatmap_" + current_scheme + ".gnuplot");
            if( io::File::TryOpenOFStream( write, filename_heatmap))
            {
              heatmap.WriteScript( write);
              io::File::CloseClearFStream( write);
            }
            else
            {
              BCL_MessageCrt( "unable to open file for writing: " + filename_heatmap);
            }
          }
        }

        // write one fused model for all donors, sites and locators
        if( m_ReplaceCombinatorialFlag->GetFlag() && !m_AppendFlag->GetFlag() && donors.GetSize() > 1)
        {
          const char scaffold_chain_id( alignments.FirstElement().First().GetSequences().FirstElement()->GetFirstMember()->GetChainID());

          storage::List< Result> fused_sequences;
          if( m_OptimalPeptideBondCutPointFlag->GetFlag())
          {
            fused_sequences.PushBack( Replace( *m_Scaffold.GetChain( scaffold_chain_id)->GetSequence(), alignments));
          }
          else
          {
            fused_sequences = ReplaceSystematic
                              (
                                *m_Scaffold.GetChain( scaffold_chain_id)->GetSequence(),
                                alignments
                              );
          }

          const util::Format count_format( util::Format().W( size_t( std::log10( fused_sequences.GetSize()) + 1)).ForceW().Fill( '0').R());
          const bool write_selected_models( m_WriteSelectedModelsFlag->GetFlag());
          storage::Set< size_t> selected_models;
          {
            const storage::Vector< size_t> temp( m_WriteSelectedModelsFlag->GetNumericalList< size_t>());
            selected_models.InsertElements( temp.Begin(), temp.End());
          }

          util::ShPtrList< sched::JobInterface> jobs;

          size_t count( 0);
          for
          (
            storage::List< Result>::iterator itr( fused_sequences.Begin()), itr_end( fused_sequences.End());
            itr != itr_end;
            ++itr, ++count
          )
          {
            Result &fused_sequence( *itr);
            fused_sequence.m_Sequence.SetFastaHeader( fused_sequence.m_ReplaceInformation);
            fused_sequence.m_FilenameBase = m_OutputPrefixFlag->GetFirstParameter()->GetValue() + "combined_" + count_format( count);
            size_t current( 1);
            for
            (
              storage::List< std::string>::const_iterator
                don_itr( donor_filenames.Begin()), don_itr_end( donor_filenames.End());
              don_itr != don_itr_end;
              ++don_itr, ++current
            )
            {
              fused_sequence.m_PymolScriptTop += "load " + *don_itr + ", " + "donor" + util::Format()( current) + "_model\n";
            }
            fused_sequence.m_PymolScriptTop += "hide all\n";
            fused_sequence.m_PymolScriptTop += "cd " + m_OutputDirectory.GetPath() + '\n';
            fused_sequence.m_PymolScriptTop += "\n#END HEADER\n";

            // write result to table
            storage::Row< double> &row( m_Results.InsertRow( fused_sequence.m_FilenameBase));
//          row[ "quality"] = quality;
            row[ "align"] = align_score_sum;
            row[ "align_norm"] = align_score_sum / nr_coord_pairs_sum;
            row[ "nr_coord_pairs"] = nr_coord_pairs_sum;
            row[ "crossover_scaffold"] = -1;
            row[ "crossover_donor"] = -1;

            // write a file for the fused model
            row[ "file"] = size_t( !write_selected_models || selected_models.Erase( count) != 0);

            // ProcessFusedSequence( fused_sequence, row);
            util::ShPtr< sched::JobInterface> sp_job
            (
              new sched::BinaryFunctionJobWithData< const Result, storage::Row< double>, void, FusionProtein>
              (
                0,
                *this,
                &FusionProtein::ProcessFusedSequence,
                fused_sequence,
                row,
                sched::JobInterface::e_READY,
                NULL
              )
            );

            sched::GetScheduler().SubmitJob( sp_job);
            jobs.PushBack( sp_job);
          }

          // iterate through all results and join
          while( jobs.GetSize() != 0)
          {
            util::ShPtr< sched::JobInterface> &sp_job( jobs.FirstElement());
            sched::GetScheduler().Join( sp_job);
            jobs.PopFront();
          }
        }

        // write the results table - sorted by if the output model was written to a file, to a results file
        m_Results.SortByColumn( "file");
        const std::string file_name_results( m_OutputPrefixFlag->GetFirstParameter()->GetValue() + "results.table");
        io::OFStream write;
        if( io::File::TryOpenOFStream( write, file_name_results))
        {
          m_Results.WriteFormatted( write);
          io::File::CloseClearFStream( write);
        }
        else
        {
          BCL_MessageCrt( "unable to open file for writing: " + file_name_results);
        }

        // end
        return 0;
      }

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

      storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, util::SiPtr< const assemble::ProteinModel> > >
      CreateAlignments
      (
        const storage::List< storage::Pair< assemble::ProteinModel, pdb::Site> > &DONORS,
        const std::string &SCAFFOLD_PDB_ID
      ) const;

      void AddSitesToDonors( storage::List< storage::Pair< assemble::ProteinModel, pdb::Site> > &DONORS) const;

      storage::List< storage::Pair< assemble::ProteinModel, pdb::Site> > Donors() const;

      int ReplaceTerm( const std::string &SCAFFOLD_PDB_ID) const;

      std::string InitializeScaffoldProtein() const;

      static const std::string &PymolScriptHeader();

      void AddScoresToRow
      (
        const assemble::ProteinModel &MODEL,
        const Result &FUSED_SEQUENCE,
        storage::Row< double> &ROW
      ) const;

      assemble::Domain DomainFromSite
      (
        const assemble::ProteinModel &PROTEIN,
        const pdb::Site &SITE
      ) const;

      util::ShPtr< function::BinarySum< const biol::AABase, const biol::AABase, double> > AssignmentScore() const;

      align::AlignmentNode< biol::AABase> AlignSingleAnchor
      (
        const std::string &PDB_ID_SITE, const std::string &PDB_ID_ANCHOR,
        const assemble::SSE &SSE_SITE, const assemble::SSE &SSE_ANCHOR,
        const storage::Set< biol::AtomType> &ATOM_TYPES,
        const util::ShPtr< quality::MeasureInterface> &SP_MEASURE
      ) const;

      //! @brief
      static storage::Map< char, storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > > >
      CollectAnchorSSEs
      (
        const assemble::ProteinModel &PROTEIN,
        const assemble::SSE &FIRST_SSE_SITE,
        const assemble::SSE &LAST_SSE_SITE,
        const util::ShPtr< assemble::LocatorAA> &SP_LOCATOR_ANCHOR_FIRST,
        const util::ShPtr< assemble::LocatorAA> &SP_LOCATOR_ANCHOR_LAST
      );

      align::AlignmentNode< biol::AABase> AlignTwoAnchors
      (
        const std::string &PDB_ID_SITE, const std::string &PDB_ID_ANCHOR,
        const assemble::SSE &FIRST_SSE_SITE, const assemble::SSE &FIRST_SSE_ANCHOR,
        const assemble::SSE &LAST_SSE_SITE, const assemble::SSE &LAST_SSE_ANCHOR
      ) const;

      Result Fuse
      (
        const biol::AASequence &SCAFFOLD,
        const biol::AASequence &DONOR,
        const math::TransformationMatrix3D &DONOR_TRANSFORMATION,
        const align::AlignmentInterface< biol::AABase> &ALIGNMENT
      ) const;

      //! @brief
      //! @param RESULT the result contains the fused sequence so far and will be elongated
      //! @param ALIGN_ITR alignment iterator at the starting position
      //! @param ALIGN_ITR_SWITCH itr where the fused sequence switches from donor to scaffold
      //! @param ALIGN_ITR_END end of alignment
      //! @param AA_SCAF_ITR current pos in scaffold sequence
      //! @param AA_SCAF_ITR_END end of scaffold sequence
      //! @param AA_DON_ITR current pos in donor sequence
      //! @param AA_DON_ITR_END current pos in scaffold sequence
      //! @param DONOR_NAME name of donor in pymol script
      //! @param DONOR_SITE_NAME name of site within fused sequence
      //! @param TRANSFORMATION transformation to apply to donor amino acids
      //! @param COUNT the cutpoint number
      //! @param SEQ_POS sequence position, if RESULT is empty, but the sequence does not start at 0
      //! @return a new result with elongated sequence and additional inforation
      Result HandleCutpointDonorScaffold
      (
        const FusionProtein::Result &RESULT,
        const align::AlignmentInterface< biol::AABase>::const_iterator &ALIGN_ITR,
        const align::AlignmentInterface< biol::AABase>::const_iterator &ALIGN_ITR_SWITCH,
        const align::AlignmentInterface< biol::AABase>::const_iterator &ALIGN_ITR_END,
        const biol::AASequence::const_iterator &AA_SCAF_ITR,
        const biol::AASequence::const_iterator &AA_SCAF_ITR_END,
        const biol::AASequence::const_iterator &AA_DON_ITR,
        const biol::AASequence::const_iterator &AA_DON_ITR_END,
        const std::string &DONOR_NAME,
        const std::string &DONOR_SITE_NAME,
        const math::TransformationMatrix3D &TRANSFORMATION,
        const size_t &COUNT,
        const size_t SEQ_POS
      ) const;

      //! @brief
      //! @param RESULT the result contains the fused sequence so far and will be elongated
      //! @param ALIGN_ITR alignment iterator at the starting position
      //! @param ALIGN_ITR_SWITCH itr where the fused sequence switches from donor to scaffold
      //! @param ALIGN_ITR_END end of alignment
      //! @param AA_SCAF_ITR current pos in scaffold sequence
      //! @param AA_SCAF_ITR_END end of scaffold sequence
      //! @param AA_DON_ITR current pos in donor sequence
      //! @param AA_DON_ITR_END current pos in scaffold sequence
      //! @param DONOR_NAME name of donor in pymol script
      //! @param DONOR_SITE_NAME name of site within fused sequence
      //! @param TRANSFORMATION transformation to apply to donor amino acids
      //! @param COUNT the cutpoint number
      //! @param SEQ_POS sequence position, if RESULT is empty, but the sequence does not start at 0
      //! @return a new result with elongated sequence and additional inforation
      Result HandleCutpointScaffoldDonor
      (
        const FusionProtein::Result &RESULT,
        const align::AlignmentInterface< biol::AABase>::const_iterator &ALIGN_ITR,
        const align::AlignmentInterface< biol::AABase>::const_iterator &ALIGN_ITR_SWITCH,
        const align::AlignmentInterface< biol::AABase>::const_iterator &ALIGN_ITR_END,
        const biol::AASequence::const_iterator &AA_SCAF_ITR,
        const biol::AASequence::const_iterator &AA_SCAF_ITR_END,
        const biol::AASequence::const_iterator &AA_DON_ITR,
        const biol::AASequence::const_iterator &AA_DON_ITR_END,
        const std::string &DONOR_NAME,
        const std::string &DONOR_SITE_NAME,
        const math::TransformationMatrix3D &TRANSFORMATION,
        const size_t &COUNT,
        const size_t SEQ_POS
      ) const;

      storage::List< Result>
      ReplaceSystematic
      (
        const biol::AASequence &SCAFFOLD,
        const storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, util::SiPtr< const assemble::ProteinModel> > > &ALIGNMENTS
      ) const;

      Result Replace
      (
        const biol::AASequence &SCAFFOLD,
        const storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, util::SiPtr< const assemble::ProteinModel> > > &ALIGNMENTS
      ) const;

      storage::List< Result>
      ReplaceNandCTermSystematic
      (
        const biol::AASequence &SCAFFOLD,
        const align::AlignmentNode< biol::AABase> &NTERM_ALIGN,
        const assemble::ProteinModel &NTERM_DONOR,
        const align::AlignmentNode< biol::AABase> &CTERM_ALIGN,
        const assemble::ProteinModel &CTERM_DONOR
      ) const;

      Result ReplaceNandCTerm
      (
        const biol::AASequence &SCAFFOLD,
        const align::AlignmentNode< biol::AABase> &NTERM_ALIGN,
        const assemble::ProteinModel &NTERM_DONOR,
        const align::AlignmentNode< biol::AABase> &CTERM_ALIGN,
        const assemble::ProteinModel &CTERM_DONOR
      ) const;

      static std::string PymolClosure
      (
        const biol::AABase &AA_LEFT,
        const biol::AABase &AA_RIGHT,
        const char SIDE,
        const size_t COUNT
      );

      static std::string &PymolImages
      (
        std::string &PYMOL_SCRIPT,
        const std::string &DONOR_NAME,
        const std::string &DONOR_SITE_NAME
      );

      static util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator FindBestPeptideBond
      (
        const align::AlignmentInterface< biol::AABase>::const_iterator &ITR_FIRST,
        const align::AlignmentInterface< biol::AABase>::const_iterator &ITR_END,
        const math::TransformationMatrix3D &TRANSFORMATION_BOTTOM,
        const bool TOP
      );

      //! @brief elongate a sequence by a given amino acid applying a transformation to it
      //! @param SEQUENCE the sequence to elongate
      //! @param AMINO_ACID the amino acid that is added to the sequence
      //! @param TRANSFORMATION the transformation to apply to the amino acid that is added
      //! @param SEQ_ID seq_id to use if given sequence is empty
      //! @return reference to the added amino acid
      static biol::AABase &ElongateSequence
      (
        biol::AASequence &SEQUENCE, const biol::AABase &AMINO_ACID,
        const math::TransformationMatrix3D &TRANSFORMATION, const int SEQ_ID
      );

      //! @brief elongate a sequence by a given amino acid
      //! @param SEQUENCE the sequence to elongate
      //! @param AMINO_ACID the amino acid that is added to the sequence
      //! @param SEQ_ID seq_id to use if given sequence is empty
      //! @return reference to the added amino acid
      static biol::AABase &ElongateSequence
      (
        biol::AASequence &SEQUENCE, const biol::AABase &AMINO_ACID, const int SEQ_ID
      );

      //! @brief create fused model from fused sequence adding all remaining scaffold chains
      //! @param FUSED_SEQUENCE the fused sequence
      //! @param SCAFFOLD scaffold protein
      //! @return protein with fused chain and remaining chains
      static assemble::ProteinModel CreateFusedModel
      (
        const biol::AASequence &FUSED_SEQUENCE,
        const assemble::ProteinModel &SCAFFOLD
      );

      //! @brief assemble multiple segments into a list of all combinations of sequence segments
      //! @param SEGMENTS a list of list of segments with varying cut points
      //! @return a list of complete sequences containing all combinations of segments
      static storage::List< Result> CombineSegments( const storage::List< storage::List< Result> > &SEGMENTS);

      //! @brief add a line for the rosetta resfile for the pair of amino acids for the given sequence position in the fused protein
      //! @param AA_FUSED the amino acid in the fused sequence
      //! @param AA_SCAFFOLD the scaffold amino acid
      //! @param AA_DONOR the donor amino acid
      //! @return the string containg the information required to define the specifics for that residue in the rosetta res file
      static std::string RosettaResFileEntry
      (
        const biol::AAData &AA_FUSED,
        const biol::AAData &AA_SCAFFOLD,
        const biol::AAData &AA_DONOR
      );

      //! @brief create a histogram that can hold scores for two cross overs (scaf->don and don->scaf)
      //! @param FUSED_SEQUENCE the sequences representing the cross over results
      //! @return Gnuplotheatmap that covers the seqids for the crossovers
      math::Histogram2D CreateCrossoverHistogram() const;

      void ProcessFusedSequence( const Result &FUSED_SEQUENCE, storage::Row< double> &ROW) const
      {
        const std::string file_name_pdb(   FUSED_SEQUENCE.m_FilenameBase + "_fused.pdb");
        const std::string file_name_fasta( FUSED_SEQUENCE.m_FilenameBase + "_fused" + FUSED_SEQUENCE.m_Sequence.GetChainID() + ".fasta");
        const std::string file_name_pml(   FUSED_SEQUENCE.m_FilenameBase + "_fused.pml");

        const assemble::ProteinModel fused_model( CreateFusedModel( FUSED_SEQUENCE.m_Sequence, m_Scaffold));

        // actual scores
        AddScoresToRow( fused_model, FUSED_SEQUENCE, ROW);

        const double quality( ROW[ "quality"]);
        const size_t nr_coord_pairs( ROW[ "nr_coord_pairs"]);

        io::OFStream write;

        // only good superimposition
        if( !m_SPMeasure->GetComparisonFunction()( quality, m_QualityCutoffParam->GetNumericalValue< double>()))
        {
          return;
        }

        // sufficiently high number of assignments
        if( nr_coord_pairs < m_MinNumberAlignedResiduesFlag->GetFirstParameter()->GetNumericalValue< size_t>())
        {
          return;
        }

        // write file for that model
        if( !bool( ROW[ "file"]))
        {
          return;
        }

        if( m_WriteRosettaResfileFlag->GetFlag())
        {
          const std::string file_name_resfile( FUSED_SEQUENCE.m_FilenameBase + "_rosetta.res");
          if( io::File::TryOpenOFStream( write, file_name_resfile))
          {
            write << "NATAA\nstart\n";
            write << FUSED_SEQUENCE.m_RosettaResFile;
            io::File::CloseClearFStream( write);
            BCL_MessageVrb( "wrote rosetta resfile to: " + file_name_resfile);
          }
          else
          {
            BCL_MessageCrt( "unable to open file for writing: " + file_name_resfile);
          }
        }

        // fasta
        if( io::File::TryOpenOFStream( write, file_name_fasta))
        {
          FUSED_SEQUENCE.m_Sequence.WriteFasta( write);
          io::File::CloseClearFStream( write);
        }
        else
        {
          BCL_MessageCrt( "unable to open file for writing: " + file_name_fasta);
        }

        // pdb
        if( io::File::TryOpenOFStream( write, file_name_pdb))
        {
          m_Factory.WriteModelToPDB( fused_model, write);
          io::File::CloseClearFStream( write);
          BCL_MessageStd( "wrote fused model to: " + file_name_pdb);
        }
        else
        {
          BCL_MessageCrt( "unable to open file for writing: " + file_name_pdb);
        }

        // pymol session file
        if( io::File::TryOpenOFStream( write, file_name_pml))
        {
          write << m_PymolScriptHeader;
          write << "load " + file_name_pdb + ", fused_model\n";
          write << FUSED_SEQUENCE.m_PymolScriptTop;
          write << FUSED_SEQUENCE.m_PymolScript;
          write << FUSED_SEQUENCE.m_PymolTailScript;
          io::File::CloseClearFStream( write);
          BCL_MessageStd( "wrote pymol_script to: " + file_name_pml);
        }
        else
        {
          BCL_MessageCrt( "unable to open file for writing: " + file_name_pml);
        }
      }

      static const ApplicationType FusionProtein_Instance;

    }; // class FusionProtein

    //! @brief default constructor
    FusionProtein::FusionProtein() :
      m_ScaffoldProteinFlag
      (
        new command::FlagStatic
        (
          "scaffold", "scaffold protein to fuse a part of the donor on",
          command::Parameter
          (
            "scaffold_pdb", "pdb file name of the scaffold pdb", command::ParameterCheckFileExistence()
          )
        )
      ),
      m_DonorProteinFlag
      (
        new command::FlagDynamic
        (
          "donor", "donor protein to select part from to add to scaffold",
          command::Parameter
          (
            "donor_pdb", "pdb file name of the donor pdb", command::ParameterCheckFileExistence()
          ),
          1, 10
        )
      ),
      m_ScaffoldFragmentFlag
      (
        new command::FlagDynamic
        (
          "scaffold_fragment", "file containing aa locators defining substitution fragment",
          command::Parameter
          (
            "file_locators", "file with two aa locators for beginning and end AA of the fragment to be replaced",
            command::ParameterCheckFileExistence(), ""
          ),
          0, 10
        )
      ),
      m_DonorSitesFlag
      (
        new command::FlagDynamic
        (
          "donor_sites", "file containing bcl site definition restricting the fragment that should be selected for fusion",
          command::Parameter
          (
            "file_sites", "file with a single pdb::Site for beginning and end AA of the fragment to be replaced",
            command::ParameterCheckFileExistence(), ""
          ),
          0, 10
        )
      ),
      m_QualityMeasureFlag
      (
        new command::FlagStatic
        (
          "quality", "quality for anchor superimposition",
          command::Parameter
          (
            "quality_measure", "superimposition measure",
            command::ParameterCheckEnumerate< quality::Measures>(), quality::GetMeasures().e_RMSD.GetName()
          )
        )
      ),
      m_QualityCutoffParam
      (
        new command::Parameter
        (
          "quality_cutoff", "cutoff value above/below which the anchor superimposition is skipped", "2.5"
        )
      ),
      m_AssignmentScoreFlag
      (
        new command::FlagDynamic
        (
          "assignment_score", "scores to use for the scoring of sequence alignments",
          command::Parameter( "score", "assignment_score with weight 1.0", command::ParameterCheckEnumerate< score::AAAssignments>()),
          0, score::GetAAAssignments().GetEnumCount()
        )
      ),
      m_WriteAlignmentFlag
      (
        new command::FlagStatic
        (
          "write_alignment", "write the alignment",
          command::Parameter
          (
            "format", "output format of the alignment",
            command::ParameterCheckEnumerate< align::HandlerClasses< biol::AABase> >(),
            align::GetHandlerClasses< biol::AABase>().e_Standard.GetName()
          )
        )
      ),
      m_WriteRosettaResfileFlag
      (
        new command::FlagStatic
        (
          "write_rosetta_res_file", "write a rosetta resfile for designing the aligned regions"
        )
      ),
      m_MinNumberAlignedResiduesFlag
      (
        new command::FlagStatic
        (
          "min_number_aligned_residues", "minimum number of aligned residues required",
          command::Parameter
          (
            "min_num_aligned_res", "minimum number of aligned residues for anchor",
            command::ParameterCheckRanged< size_t>( 0, 100),
            "2"
          )
        )
      ),
      m_WriteSelectedModelsFlag
      (
        new command::FlagDynamic
        (
          "write_selected_models",
          "list of model numbers to be written",
          command::Parameter
          (
            "model_number",
            "number of the model to be written",
            command::ParameterCheckRanged< size_t>( 0, 999999),
            "0"
          ),
          0,
          100
        )
      ),
      m_WriteCrossOverHeatmapFlag
      (
        new command::FlagStatic
        (
          "cross_over_heatmap",
          "write heatmaps for all scores and cross over points"
        )
      ),
      m_OutputPrefixFlag
      (
        new command::FlagStatic
        (
          "output_prefix", "prefix for output",
          command::Parameter( "prefix", "prefix to use for all file names", "")
        )
      ),
      m_AppendFlag
      (
        new command::FlagStatic
        (
          "append", "append the site to the scaffold instead of replacing it"
        )
      ),
      m_ReplaceTermFlag
      (
        new command::FlagStatic
        (
          "replace_term", "replace C and N term of the scaffold"
        )
      ),
      m_NTermDonorParam
      (
        new command::Parameter
        (
          "nterm_donor", "pdb file for N-term donor",
          command::ParameterCheckFileExistence(),
          ""
        )
      ),
      m_NTermScaffoldLocatorParam
      (
        new command::Parameter
        (
          "nterm_scaffold_locator", "SSE locator for N-term fusion in scaffold",
          command::ParameterCheckFileExistence(),
          ""
        )
      ),
      m_NTermDonorLocatorParam
      (
        new command::Parameter
        (
          "nterm_donor_locator", "SSE locator for N-term fusion in donor",
          command::ParameterCheckFileExistence(),
          ""
        )
      ),
      m_CTermDonorParam
      (
        new command::Parameter
        (
          "cterm_donor", "pdb file for C-term donor",
          command::ParameterCheckFileExistence(),
          ""
        )
      ),
      m_CTermScaffoldLocatorParam
      (
        new command::Parameter
        (
          "cterm_scaffold_locator", "SSE locator for C-term fusion in scaffold",
          command::ParameterCheckFileExistence(),
          ""
        )
      ),
      m_CTermDonorLocatorParam
      (
        new command::Parameter
        (
          "cterm_donor_locator", "sse locator for C-term fusion in donor",
          command::ParameterCheckFileExistence(),
          ""
        )
      ),
      m_WeightSetFlag
      (
        new command::FlagStatic
        (
          "weight_set",
          "file containing storage::Table< double> with weight set",
          command::Parameter( "weight_set_file_name", "filename of file containing weight set", "")
        )
      ),
      m_OptimalPeptideBondCutPointFlag
      (
        new command::FlagStatic
        (
          "cutpoint_optimal_peptide",
          "choose the cutpoint between two sequences in the alignment by the geometrically best peptide bond instead "
          "of systematically generating all possible fusion proteins"
        )
      ),
      m_ReplaceCombinatorialFlag
      (
        new command::FlagStatic
        (
          "replace_combinatorial",
          "create the combination of all crossover points over all replacement sites as models and score them. Use with"
          " \"write_selected_models\" flag so that not every possible model is written!!"
        )
      ),
      m_Factory( biol::GetAAClasses().e_AAComplete)
    {
      m_QualityMeasureFlag->PushBack( m_QualityCutoffParam);

      m_ReplaceTermFlag->PushBack( m_NTermDonorParam);
      m_ReplaceTermFlag->PushBack( m_NTermScaffoldLocatorParam);
      m_ReplaceTermFlag->PushBack( m_NTermDonorLocatorParam);
      m_ReplaceTermFlag->PushBack( m_CTermDonorParam);
      m_ReplaceTermFlag->PushBack( m_CTermScaffoldLocatorParam);
      m_ReplaceTermFlag->PushBack( m_CTermDonorLocatorParam);
    }

    storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, util::SiPtr< const assemble::ProteinModel> > >
    FusionProtein::CreateAlignments
    (
      const storage::List< storage::Pair< assemble::ProteinModel, pdb::Site> > &DONORS,
      const std::string &SCAFFOLD_PDB_ID
    ) const
    {
      // get locators for anchor, if given in commandline
      storage::Vector< storage::VectorND< 2, util::ShPtr< assemble::LocatorAA> > > scaffold_anchor_locators;
      if( m_ScaffoldFragmentFlag->GetFlag())
      {
        const storage::Vector< std::string> locator_file_names( m_ScaffoldFragmentFlag->GetStringList());
        BCL_Assert( DONORS.GetSize() == locator_file_names.GetSize(), "need as many locator definitions as donors");

        for
        (
          storage::Vector< std::string>::const_iterator
            loc_itr( locator_file_names.Begin()), loc_itr_end( locator_file_names.End());
          loc_itr != loc_itr_end;
          ++loc_itr
        )
        {
          io::IFStream read;
          io::File::MustOpenIFStream( read, *loc_itr);
          util::ShPtr< assemble::LocatorAA> sp_locator_anchor_first( new assemble::LocatorAA());
          util::ShPtr< assemble::LocatorAA> sp_locator_anchor_last( new assemble::LocatorAA());
          read >> *sp_locator_anchor_first;
          read >> *sp_locator_anchor_last;
          BCL_Assert
          (
               sp_locator_anchor_first->GetLocatorChain().GetChainID()
            == sp_locator_anchor_last->GetLocatorChain().GetChainID(),
            "given locators for scaffold are from different chains"
          );
          io::File::CloseClearFStream( read);
          BCL_MessageStd
          (
            "selecting fragment to replace:\n" + util::Format()( *sp_locator_anchor_first) +
            "\nand\n" + util::Format()( *sp_locator_anchor_last)
          );
          scaffold_anchor_locators.PushBack
          (
            storage::VectorND< 2, util::ShPtr< assemble::LocatorAA> >
            (
              sp_locator_anchor_first, sp_locator_anchor_last
            )
          );
        }
      }

      size_t donor_nr( 0);
      storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, util::SiPtr< const assemble::ProteinModel> > > alignments;

      for
      (
        storage::List< storage::Pair< assemble::ProteinModel, pdb::Site> >::const_iterator
          don_itr( DONORS.Begin()), don_itr_end( DONORS.End());
        don_itr != don_itr_end;
        ++don_itr, ++donor_nr
      )
      {
        util::ShPtr< assemble::LocatorAA> sp_locator_anchor_first;
        util::ShPtr< assemble::LocatorAA> sp_locator_anchor_last;
        if( !scaffold_anchor_locators.IsEmpty())
        {
          sp_locator_anchor_first = scaffold_anchor_locators( donor_nr).First();
          sp_locator_anchor_last = scaffold_anchor_locators( donor_nr).Second();
        }
        // donor protein
        const assemble::ProteinModel &donor_protein( don_itr->First());

        // site
        const pdb::Site &current_site( don_itr->Second());

        // skipping sites without ligand
        if( !current_site.GetLigand().IsDefined() && current_site.GetDescription().find( "BINDING") == std::string::npos)
        {
          BCL_MessageVrb( "skipping site without ligand: " + current_site.GetName());
          continue;
        }

        BCL_MessageVrb( "trying to fuse site: " + current_site.GetName());

//          // ligand
//          const pdb::Ligand &current_ligand( *current_site.GetLigand());
//
        const assemble::Domain site_domain( DomainFromSite( donor_protein, current_site));

        // search for an sse type equal to the domains first sse
        const util::SiPtr< const assemble::SSE> first_sse_site( site_domain.GetSSEs().FirstElement());
        const util::SiPtr< const assemble::SSE> last_sse_site( site_domain.GetSSEs().LastElement());

        const storage::Map< char, storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > > > anchor_sses
        (
          CollectAnchorSSEs( m_Scaffold, *first_sse_site, *last_sse_site, sp_locator_anchor_first, sp_locator_anchor_last)
        );

        BCL_MessageDbg( "number of chains in anchor SSEs: " + util::Format()( anchor_sses.GetSize()));

        // iterate through anchors
        for
        (
          storage::Map< char, storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > > >::const_iterator
            chain_itr( anchor_sses.Begin()), chain_itr_end( anchor_sses.End());
          chain_itr != chain_itr_end;
          ++chain_itr
        )
        {
          // iterate through anchors
          for
          (
            storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > >::const_iterator
              anchor_itr( chain_itr->second.Begin()), anchor_itr_end( chain_itr->second.End());
            anchor_itr != anchor_itr_end;
            ++anchor_itr
          )
          {
            const assemble::SSE &first_sse( *anchor_itr->First());
            const assemble::SSE &last_sse( *anchor_itr->Second());

            BCL_MessageStd
            (
              "aligning:\n" +
              first_sse_site->GetIdentification() + '\t' + first_sse.GetIdentification() + '\t' +
              last_sse_site->GetIdentification() + '\t' + last_sse.GetIdentification()
            );

            const util::ShPtr< util::Wrapper< std::string> > sp_id
            (
              donor_protein.GetProteinModelData()->GetData( assemble::ProteinModelData::e_Identification)
            );
            BCL_Assert( sp_id.IsDefined(), "cannot find pdb id for donor");

            // found two sses, that could be potential anchor sses
            const align::AlignmentNode< biol::AABase> alignment
            (
              AlignTwoAnchors
              (
                sp_id->GetData(), SCAFFOLD_PDB_ID,
                first_sse, *first_sse_site,
                last_sse, *last_sse_site
              )
            );

            if( alignment.IsEmpty())
            {
              BCL_MessageVrb( "Superimpose returned empty alignment");
              continue;
            }

            alignments.PushBack
            (
              storage::Pair< align::AlignmentNode< biol::AABase>, util::SiPtr< const assemble::ProteinModel> >
              (
                alignment, util::SiPtr< const assemble::ProteinModel>( donor_protein)
              )
            );
          } // iterate through anchors
        } // iterate through anchors chain wise
      } // iterate through donors

      return alignments;
    }

    //! @brief add site definition to the donors
    //! @param DONORS list all donors
    void FusionProtein::AddSitesToDonors( storage::List< storage::Pair< assemble::ProteinModel, pdb::Site> > &DONORS) const
    {
      const storage::Vector< std::string> site_file_names( m_DonorSitesFlag->GetStringList());
      BCL_Assert( DONORS.GetSize() == site_file_names.GetSize(), "need as many site definitions as donors");

      storage::List< storage::Pair< assemble::ProteinModel, pdb::Site> >::iterator don_itr( DONORS.Begin());

      for
      (
        storage::Vector< std::string>::const_iterator
          site_itr( site_file_names.Begin()), site_itr_end( site_file_names.End());
        site_itr != site_itr_end;
        ++site_itr, ++don_itr
      )
      {
        io::IFStream read;
        io::File::MustOpenIFStream( read, *site_itr);
        read >> don_itr->Second();
        io::File::CloseClearFStream( read);
        BCL_MessageStd( "selecting site from donor:\n" + util::Format()( don_itr->Second()));
      }
    }

    //! @brief read all donors from the commandline
    //! @return list of donors
    storage::List< storage::Pair< assemble::ProteinModel, pdb::Site> > FusionProtein::Donors() const
    {
      io::IFStream read;

      const storage::Map< biol::SSType, size_t> min_sse_size( pdb::Factory::GetCommandlineSSETypeMinSizes());

      const storage::Vector< std::string> donor_file_names( m_DonorProteinFlag->GetStringList());
      storage::List< storage::Pair< assemble::ProteinModel, pdb::Site> > donors;

      for
      (
        storage::Vector< std::string>::const_iterator
          don_itr( donor_file_names.Begin()), don_itr_end( donor_file_names.End());
        don_itr != don_itr_end;
        ++don_itr
      )
      {
        io::File::MustOpenIFStream( read, *don_itr);
        pdb::Handler donor_handler( read);
        io::File::CloseClearFStream( read);

        util::ShPtrList< pdb::Site> sites;
        sites = donor_handler.GetSites();
        const pdb::Site site( sites.IsEmpty() ? pdb::Site() : *sites.FirstElement());
        // donor protein
        assemble::ProteinModel donor_protein( m_Factory.ProteinModelFromPDB( donor_handler, min_sse_size));
        {
          util::ShPtr< util::Wrapper< std::string> > sp_id
          (
            new util::Wrapper< std::string>( donor_handler.GetHead().GetPDBID())
          );
          if( sp_id->GetData() == pdb::Head::GetBCLPdbID())
          {
            sp_id->GetData() = io::File::RemoveFullExtension( io::File::RemovePath( *don_itr));
          }

          util::ShPtr< assemble::ProteinModelData> sp_model_data( donor_protein.GetProteinModelData());
          sp_model_data->Insert( assemble::ProteinModelData::e_Identification, sp_id);
        }

        donors.PushBack( storage::Pair< assemble::ProteinModel, pdb::Site>( donor_protein, site));
      }

      return donors;
    }

    //! @brief replace just the terminal ends of the scaffold
    int FusionProtein::ReplaceTerm( const std::string &SCAFFOLD_PDB_ID) const
    {
      io::IFStream read;

      const storage::Map< biol::SSType, size_t> min_sse_size( pdb::Factory::GetCommandlineSSETypeMinSizes());

      // read locators and donor proteins
      const io::DirectoryEntry nterm_donor( m_NTermDonorParam->GetValue());
      const io::DirectoryEntry cterm_donor( m_CTermDonorParam->GetValue());

      io::File::MustOpenIFStream( read, nterm_donor.GetFullName());
      pdb::Handler nterm_donor_handler( read);
      io::File::CloseClearFStream( read);
      io::File::MustOpenIFStream( read, cterm_donor.GetFullName());
      pdb::Handler cterm_donor_handler( read);
      io::File::CloseClearFStream( read);

      // nterm donor protein
      assemble::ProteinModel nterm_donor_protein( m_Factory.ProteinModelFromPDB( nterm_donor_handler, min_sse_size));
      util::ShPtr< util::Wrapper< std::string> > sp_nterm_donor_id
      (
        new util::Wrapper< std::string>( nterm_donor_handler.GetHead().GetPDBID())
      );
      if( sp_nterm_donor_id->GetData() == pdb::Head::GetBCLPdbID())
      {
        sp_nterm_donor_id->GetData() = io::File::RemoveFullExtension( nterm_donor.GetName());
      }
      {
        util::ShPtr< assemble::ProteinModelData> sp_model_data( nterm_donor_protein.GetProteinModelData());
        sp_model_data->Insert( assemble::ProteinModelData::e_Identification, sp_nterm_donor_id);
      }

      // cterm donor protein
      assemble::ProteinModel cterm_donor_protein( m_Factory.ProteinModelFromPDB( cterm_donor_handler, min_sse_size));
      util::ShPtr< util::Wrapper< std::string> > sp_cterm_donor_id
      (
        new util::Wrapper< std::string>( cterm_donor_handler.GetHead().GetPDBID())
      );
      if( sp_cterm_donor_id->GetData() == pdb::Head::GetBCLPdbID())
      {
        sp_cterm_donor_id->GetData() = io::File::RemoveFullExtension( cterm_donor.GetName());
      }
      {
        util::ShPtr< assemble::ProteinModelData> sp_model_data( cterm_donor_protein.GetProteinModelData());
        sp_model_data->Insert( assemble::ProteinModelData::e_Identification, sp_cterm_donor_id);
      }

      assemble::LocatorSSE nterm_donor_locator;
      io::File::MustOpenIFStream( read, m_NTermDonorLocatorParam->GetValue());
      read >> nterm_donor_locator;
      io::File::CloseClearFStream( read);

      assemble::LocatorSSE nterm_scaffold_locator;
      io::File::MustOpenIFStream( read, m_NTermScaffoldLocatorParam->GetValue());
      read >> nterm_scaffold_locator;
      io::File::CloseClearFStream( read);

      assemble::LocatorSSE cterm_donor_locator;
      io::File::MustOpenIFStream( read, m_CTermDonorLocatorParam->GetValue());
      read >> cterm_donor_locator;
      io::File::CloseClearFStream( read);

      assemble::LocatorSSE cterm_scaffold_locator;
      io::File::MustOpenIFStream( read, m_CTermScaffoldLocatorParam->GetValue());
      read >> cterm_scaffold_locator;
      io::File::CloseClearFStream( read);

      // locate nterm sses and create alignment
      const util::SiPtr< const assemble::SSE> sp_sse_nterm_scaffold( nterm_scaffold_locator.Locate( m_Scaffold));
      const util::SiPtr< const assemble::SSE> sp_sse_nterm_donor( nterm_donor_locator.Locate( nterm_donor_protein));
      BCL_Assert( sp_sse_nterm_scaffold.IsDefined() && sp_sse_nterm_donor.IsDefined(), "cannot locate N-term SSEs");
      const align::AlignmentNode< biol::AABase> nterm_alignment
      (
        AlignSingleAnchor
        (
          SCAFFOLD_PDB_ID, sp_nterm_donor_id->GetData(),
          *sp_sse_nterm_scaffold, *sp_sse_nterm_donor,
          m_AtomTypes, m_SPMeasure
        )
      );

      // locate nterm sses and create alignment
      const util::SiPtr< const assemble::SSE> sp_sse_cterm_scaffold( cterm_scaffold_locator.Locate( m_Scaffold));
      const util::SiPtr< const assemble::SSE> sp_sse_cterm_donor( cterm_donor_locator.Locate( cterm_donor_protein));
      BCL_Assert( sp_sse_cterm_scaffold.IsDefined() && sp_sse_cterm_donor.IsDefined(), "cannot locate C-term SSEs");
      const align::AlignmentNode< biol::AABase> cterm_alignment
      (
        AlignSingleAnchor
        (
          SCAFFOLD_PDB_ID, sp_cterm_donor_id->GetData(),
          *sp_sse_cterm_scaffold, *sp_sse_cterm_donor,
          m_AtomTypes, m_SPMeasure
        )
      );

      // sequence from the scaffold with associated with the nterm SSE
      const biol::AASequence &scaffold_sequence
      (
        *m_Scaffold.GetChain( sp_sse_nterm_scaffold->GetChainID())->GetSequence()
      );

      size_t nr_coord_pairs_sum( 0);
      double align_fuse_score( 0.0);

      // nterm
      const std::string nterm_donor_filename
      (
        m_OutputPrefixFlag->GetFirstParameter()->GetValue() + "nterm_donor.pdb"
      );

      {
        io::OFStream write;
        // score of the alignment
        const double align_score( m_SPScoreAlignmentAssignment->operator ()( nterm_alignment));
        align_fuse_score += align_score;

        const storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > coord_pair
        (
          assemble::Quality::CoordinatesFromAlignment( nterm_alignment, m_AtomTypes)
        );
        const size_t nr_coord_pairs( coord_pair.First().GetSize());
        nr_coord_pairs_sum += nr_coord_pairs;
        BCL_MessageStd( "number coord pairs: " + util::Format()( nr_coord_pairs));

//            const double quality( measure->CalculateMeasure( coord_pair.First(), coord_pair.Second()));
        const math::TransformationMatrix3D transformation
        (
          m_SPMeasure->CalculateSuperimposition( coord_pair.Second(), coord_pair.First())
        );

        // write donors
        if( io::File::TryOpenOFStream( write, nterm_donor_filename))
        {
          nterm_donor_protein.Transform( transformation);
          m_Factory.WriteModelToPDB( nterm_donor_protein, write);
          io::File::CloseClearFStream( write);
        }
        else
        {
          BCL_MessageCrt( "unable to open file for writing: " + nterm_donor_filename);
        }
      }

      // cterm
      const std::string cterm_donor_filename( m_OutputPrefixFlag->GetFirstParameter()->GetValue() + "cterm_donor.pdb");
      {
        io::OFStream write;
        // score of the alignment
        const double align_score( m_SPScoreAlignmentAssignment->operator ()( cterm_alignment));
        align_fuse_score += align_score;

        const storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > coord_pair
        (
          assemble::Quality::CoordinatesFromAlignment( cterm_alignment, m_AtomTypes)
        );
        const size_t nr_coord_pairs( coord_pair.First().GetSize());
        nr_coord_pairs_sum += nr_coord_pairs;
        BCL_MessageStd( "number coord pairs: " + util::Format()( nr_coord_pairs));

//            const double quality( measure->CalculateMeasure( coord_pair.First(), coord_pair.Second()));
        const math::TransformationMatrix3D transformation
        (
          m_SPMeasure->CalculateSuperimposition( coord_pair.Second(), coord_pair.First())
        );

        // write donors
        if( io::File::TryOpenOFStream( write, cterm_donor_filename))
        {
          cterm_donor_protein.Transform( transformation);
          m_Factory.WriteModelToPDB( cterm_donor_protein, write);
          io::File::CloseClearFStream( write);
        }
        else
        {
          BCL_MessageCrt( "unable to open file for writing: " + cterm_donor_filename);
        }
      }

      storage::List< Result> fused_sequences;

      if( m_OptimalPeptideBondCutPointFlag->GetFlag())
      {
        fused_sequences.PushBack
        (
          ReplaceNandCTerm
          (
            scaffold_sequence,
            nterm_alignment, nterm_donor_protein,
            cterm_alignment, cterm_donor_protein
          )
        );
      }
      else
      {
        fused_sequences = ReplaceNandCTermSystematic
                          (
                            scaffold_sequence,
                            nterm_alignment, nterm_donor_protein,
                            cterm_alignment, cterm_donor_protein
                          );
      }

      util::Format number_format;
      number_format.Fill( '0');
      number_format.W( std::log10( fused_sequences.GetSize()) + 1);

      size_t count( 0);
      // iterate through all fused sequences
      for
      (
        storage::List< Result>::const_iterator result_itr( fused_sequences.Begin()), result_itr_end( fused_sequences.End());
        result_itr != result_itr_end;
        ++result_itr, ++count
      )
      {
        const Result fused_sequence( *result_itr);

        io::OFStream write;
        const std::string file_name_base( m_OutputPrefixFlag->GetFirstParameter()->GetValue() + "ncterm_" + number_format( count));
        const std::string file_name_pdb( file_name_base + "fused.pdb");
        const std::string file_name_fasta( file_name_base + "sequence.fasta");
        const std::string file_name_pml( file_name_base + "pymol.pml");

        // actual scores
        const assemble::ProteinModel fused_model( CreateFusedModel( fused_sequence.m_Sequence, m_Scaffold));

        // write result to table
        storage::Row< double> &row( m_Results.InsertRow( file_name_base));
        AddScoresToRow( fused_model, fused_sequence, row);
        // row[ "quality"] = quality;
        row[ "align"] = align_fuse_score;
        row[ "align_norm"] = align_fuse_score / nr_coord_pairs_sum;
        row[ "nr_coord_pairs"] = nr_coord_pairs_sum;
        // write to file
        row[ "file"] = 1;

        // fasta
        if( io::File::TryOpenOFStream( write, file_name_fasta))
        {
          fused_sequence.m_Sequence.WriteFasta( write);
          io::File::CloseClearFStream( write);
        }
        else
        {
          BCL_MessageCrt( "unable to open file for writing: " + file_name_fasta);
        }

        // pdb and pymol session file
        if( io::File::TryOpenOFStream( write, file_name_pdb))
        {
          m_Factory.WriteModelToPDB( fused_model, write);
          io::File::CloseClearFStream( write);
          BCL_MessageStd( "wrote fused model to: " + file_name_pdb);
        }
        else
        {
          BCL_MessageCrt( "unable to open file for writing: " + file_name_pdb);
        }

        if( io::File::TryOpenOFStream( write, file_name_pml))
        {
          const io::Directory dir
          (
            io::File::SplitToPathAndFileName( m_OutputPrefixFlag->GetFirstParameter()->GetValue()).First()
          );

          write << m_PymolScriptHeader << "\n#END HEADER\n";
          write << "load " + file_name_pdb + ", fused_model\n";
          write << "load " + nterm_donor_filename + ", donor_nterm_model\n";
          write << "load " + cterm_donor_filename + ", donor_cterm_model\n";
          write << "hide all\n";
          write << "cd " + dir.GetPath() + '\n';
          write << fused_sequence.m_PymolScript;
          write << fused_sequence.m_PymolTailScript;
          write << "enable all\n";
          write << "disable peptide_good\n";
          write << "disable peptide_bad\n";
          write << "disable peptide_atoms\n";

          io::File::CloseClearFStream( write);
          BCL_MessageStd( "wrote pymol_script to: " + file_name_pml);
        }
        else
        {
          BCL_MessageCrt( "unable to open file for writing: " + file_name_pml);
        }
      }

      m_Results.SortByColumn( "file");
      m_Results.WriteFormatted( util::GetLogger());

      return 0;
    }

    //! @brief initialize the scaffold protein object from filename that was passed
    //! @return the pdb_id of the scaffold
    std::string FusionProtein::InitializeScaffoldProtein() const
    {
      const storage::Map< biol::SSType, size_t> min_sse_size( pdb::Factory::GetCommandlineSSETypeMinSizes());

      io::IFStream read;
      io::File::MustOpenIFStream( read, m_ScaffoldProteinFlag->GetFirstParameter()->GetValue());
      pdb::Handler scaffold_handler( read);
      io::File::CloseClearFStream( read);

      m_Scaffold = m_Factory.ProteinModelFromPDB( scaffold_handler, min_sse_size);
      {
        // reduce scaffold to non symmetrical chains:
        const assemble::Biomolecule biomolecule
        (
          biol::GetAtomTypes().GetBackBoneAtomTypes(),
          *quality::SuperimposeMeasure( m_QualityMeasureFlag->GetFirstParameter()->GetValue()),
          m_QualityCutoffParam->GetNumericalValue< double>()
        );
        math::MutateResult< assemble::ProteinModel> result( biomolecule( m_Scaffold));
        if( result.GetArgument().IsDefined())
        {
          BCL_MessageStd
          (
            "reduced the size of the scaffold by chain superimposition from " +
            util::Format()( m_Scaffold.GetNumberOfChains()) + " chains to " +
            util::Format()( result.GetArgument()->GetNumberOfChains()) + " chains!"
          );
          m_Scaffold = *result.GetArgument();
        }
      }

      std::string scaffold_pdb_id( scaffold_handler.GetHead().GetPDBID());
      if( scaffold_pdb_id == pdb::Head::GetBCLPdbID())
      {
        scaffold_pdb_id = io::File::RemoveFullExtension
            (
              io::File::RemovePath( m_ScaffoldProteinFlag->GetFirstParameter()->GetValue())
            );
      }

      // reassign secondary structure for scoring
      {
        // actual scores for scaffold
        storage::Row< double> &row( m_Results.InsertRow( "scaffold"));
        biol::DSSP dssp;
        math::MutateResult< assemble::ProteinModel> result( dssp( m_Scaffold));
        if( result.GetArgument().IsDefined())
        {
          // actual scores
          AddScoresToRow( *result.GetArgument(), Result(), row);
        }
        else
        {
          // actual scores
          AddScoresToRow( m_Scaffold, Result(), row);
        }
      }

      // write the scaffold
      io::OFStream write;
      const std::string scaffold_filename( m_OutputPrefixFlag->GetFirstParameter()->GetValue() + "scaffold.pdb");
      if( io::File::TryOpenOFStream( write, scaffold_filename))
      {
        m_Factory.WriteModelToPDB( m_Scaffold, write);
        io::File::CloseClearFStream( write);
      }
      else
      {
        BCL_MessageCrt( "unable to open file for writing: " + scaffold_filename);
      }

      m_PymolScriptHeader += "load " + scaffold_filename + ", scaffold_model\n";

      return scaffold_pdb_id;
    }

    const std::string &FusionProtein::PymolScriptHeader()
    {
      static const std::string s_pymol_script_header
      (
        "#BEGIN HEADER\n"
        "bg_color white\n"
        "set_color scaffold_color= [0.95 , 0.78 , 0.00]\n"
        "set_color overlap_color=  [0.00 , 0.53 , 0.22]\n"
        "set_color donor_color=    [0.02 , 0.50 , 0.72]\n"
        "set_color ignore_color=   [1.00 , 0.00 , 0.00]\n"
      );

      return s_pymol_script_header;
    }

    //! @brief helper function to score a model, and add all scores to argument ROW
    //! @param MODEL the model to score
    //! @param ROW the row to set the fields according to the scoring function
    void FusionProtein::AddScoresToRow
    (
      const assemble::ProteinModel &MODEL,
      const Result &FUSED_SEQUENCE,
      storage::Row< double> &ROW
    ) const
    {
      // actual scores
      const storage::Table< double> scores( m_ScoreFunction.CreateValueTableHorizontal( MODEL));

      // insert all the scores into the ROW
      for
      (
        storage::TableHeader::const_iterator itr( scores.GetHeader().Begin()), itr_end( scores.GetHeader().End());
        itr != itr_end;
        ++itr
      )
      {
        // if their is no such col in the scoring table, just skip
        if( !ROW.GetHeader().HasColumn( *itr))
        {
          continue;
        }

        // copy value from scores to the scoring table
        ROW[ *itr] = scores[ "weighted_value"][ *itr];
      }

      // copy weighted sum
      ROW[ "sum"] = scores[ "weighted_value"][ "sum"];

      // copy the fused seuqence stats to the row as well
      const size_t nr_fused_res( FUSED_SEQUENCE.m_NrFusedResiduesNTerm + FUSED_SEQUENCE.m_NrFusedResiduesCTerm);
      ROW[ "align_fuse"] = FUSED_SEQUENCE.m_AlignFuseScore;
      ROW[ "align_fuse_norm"] = FUSED_SEQUENCE.m_AlignFuseScore / nr_fused_res;
      ROW[ "nr_fused_pairs"] = nr_fused_res;
      ROW[ "nr_fused_pairs_first"] = FUSED_SEQUENCE.m_NrFusedResiduesNTerm;
      ROW[ "nr_fused_pairs_last"] = FUSED_SEQUENCE.m_NrFusedResiduesCTerm;
      ROW[ "nr_replaced"] = FUSED_SEQUENCE.m_NrReplacedResidues;
    }

    //! @brief collect all sses that belong to a SITE within the MODEL
    //! @param PROTEIN the protein the site is a part of
    //! @param SITE the site definition containing residue information
    //! @return the domain containing all SSES between the first and the last residue within a site
    assemble::Domain FusionProtein::DomainFromSite
    (
      const assemble::ProteinModel &PROTEIN,
      const pdb::Site &SITE
    ) const
    {
      const storage::List< assemble::LocatorAA> aa_locators( SITE.AALocators());

      BCL_MessageVrb
      (
        "number of residues associate with site: " + SITE.GetName() + " " +
        util::Format()( aa_locators.GetSize())
      );

      assemble::Domain site_domain;

      // iterate through all sses in the model, and check if the sse contains one of the aas
      const util::SiPtrVector< const assemble::SSE> model_sses( PROTEIN.GetSSEs());

      for
      (
        util::SiPtrVector< const assemble::SSE>::const_iterator
          sse_itr( model_sses.Begin()), sse_itr_end( model_sses.End());
        sse_itr != sse_itr_end;
        ++sse_itr
      )
      {
        // iterate through all aa locators
        for
        (
          storage::List< assemble::LocatorAA>::const_iterator
            loc_itr( aa_locators.Begin()), loc_itr_end( aa_locators.End());
          loc_itr != loc_itr_end;
          ++loc_itr
        )
        {
          // TODO: check for the first residue, how many residues are to the left; for the last how many residues are to the right;
          // if they are within SSEs (helix, strand) and there are not enough residues space to the border, take the next helix or strand, as attachment points for the fragment
          // if the amino acid can be found in this sse
          if( loc_itr->Locate( **sse_itr).IsDefined())
          {
            // find this sse in the model and insert it into the domain
            const util::ShPtr< assemble::SSE> sp_sse( PROTEIN.FindSSE( **sse_itr));
            if( sp_sse.IsDefined())
            {
              site_domain.Insert( sp_sse);
            }

            // move on to next sse
            break;
          }
        } // locator
      } // sses

      // domain for that site
      // check that is only a single chain
      if( site_domain.GetChainIds().GetSize() > 1)
      {
        BCL_MessageVrb( "skipping site with residues in multiple chains: " + SITE.GetName());
        return assemble::Domain();
      }

      // also add the sse before the first and after the last SSE in that fragment, and add the missing SSEs in between
      if( site_domain.IsEmpty())
      {
        BCL_MessageVrb( "unable to identify domain for site: " + SITE.GetName());
        return assemble::Domain();
      }

      const util::ShPtr< assemble::SSE> &first_sse_domain( *site_domain.GetData().Begin());
      const util::ShPtr< assemble::SSE> &last_sse_domain( *site_domain.GetData().ReverseBegin());

      // iterate over model sses
      util::SiPtrVector< const assemble::SSE>::const_iterator
        sse_itr( model_sses.Begin()), sse_itr_end( model_sses.End());

      util::SiPtr< const assemble::SSE> preceeding_sse( **sse_itr);

      // there is no preceeding sse, if the first sse in the model is the first sse in the site domain
      if( first_sse_domain != preceeding_sse)
      {
        while( sse_itr != sse_itr_end)
        {
          if( ( *sse_itr)->DoesPrecede( *first_sse_domain))
          {
            preceeding_sse = *sse_itr;
            break;
          }
          ++sse_itr;
        }
      }

      // fill the site with all sses in between
      while( sse_itr != sse_itr_end && *sse_itr != last_sse_domain)
      {
        // find this sse in the model and insert it into the domain
        const util::ShPtr< assemble::SSE> sp_sse( PROTEIN.FindSSE( **sse_itr));
        if( sp_sse.IsDefined())
        {
          site_domain.Insert( sp_sse);
        }
        ++sse_itr;
      }
      // goto proceeding sse
      util::SiPtr< const assemble::SSE> proceeding_sse( last_sse_domain);
      if( sse_itr != sse_itr_end && ++sse_itr != sse_itr_end)
      {
        proceeding_sse = *sse_itr;
      }

      BCL_MessageStd
      (
        "site:\npreceeding_sse: " + preceeding_sse->GetIdentification() +
        "\nfirst_sse:      " + first_sse_domain->GetIdentification() +
        "\nlast_sse:       " + last_sse_domain->GetIdentification() +
        "\nproceeding_sse: " + proceeding_sse->GetIdentification()
      );

      if( preceeding_sse->GetType() != biol::GetSSTypes().COIL)
      {
        // insert into domain
        // find this sse in the model and insert it into the domain
        const util::ShPtr< assemble::SSE> sp_sse( PROTEIN.FindSSE( *preceeding_sse));
        if( sp_sse.IsDefined())
        {
          site_domain.Insert( sp_sse);
        }
      }
      if( proceeding_sse->GetType() != biol::GetSSTypes().COIL)
      {
        // insert into domain
        // find this sse in the model and insert it into the domain
        const util::ShPtr< assemble::SSE> sp_sse( PROTEIN.FindSSE( *proceeding_sse));
        if( sp_sse.IsDefined())
        {
          site_domain.Insert( sp_sse);
        }
      }

      // end
      return site_domain;
    }

    //! @brief generate the assignment score function according to the commandline arguments
    //! @return ShPtr to the assignment scoring function
    util::ShPtr< function::BinarySum< const biol::AABase, const biol::AABase, double> > FusionProtein::AssignmentScore() const
    {
      // initialize the assignment score
      util::ShPtr< function::BinarySum< const biol::AABase, const biol::AABase, double> > sp_assignment_score
      (
        new function::BinarySum< const biol::AABase, const biol::AABase, double>()
      );

      // get set of assignment scores
      storage::Set< score::AAAssignment> assign_score_enums( m_AssignmentScoreFlag->GetObjectSet< score::AAAssignment>());

      size_t number_scores( 0);
      // iterate through all scores
      for
      (
        storage::Set< score::AAAssignment>::const_iterator
          itr( assign_score_enums.Begin()), itr_end( assign_score_enums.End());
        itr != itr_end;
        ++itr
      )
      {
        const util::ShPtr< math::ZScore> z_score( score::GetAAAssignments().GetZScore( *itr));
        if( !z_score.IsDefined())
        {
          BCL_MessageStd( "skipping assignment score, due to missing z-score: " + itr->GetName());
          continue;
        }

        // add to scoring function
        *sp_assignment_score += function::BinaryAdapter< const biol::AABase, const biol::AABase, double, const double, double>( **itr, z_score);
        ++number_scores;
      }

      *sp_assignment_score /= double( number_scores);

      return sp_assignment_score;
    }

    //! @brief given two sses (one from the scaffold structure and one from the site from another structure) an
    //!        alignment with optimal superimposition is created
    //! @param PDB_ID_SITE descriptor for the structure of the site protein
    //! @param PDB_ID_ANCHOR descriptor for the structure of the scaffold protein
    //! @param SSE_SITE the SSE at the left/right side of the site
    //! @param SSE_ANCHOR the SSE at the right/left side of the scaffold
    //! @param ATOM_TYPES the atoms to be used for the superimpostion
    //! @param SP_MEASURE the measure to be used for the superimposition
    //! @return alignment of the sequences of the site SSE (top) and scaffold/anchor sse (bottom)
    align::AlignmentNode< biol::AABase> FusionProtein::AlignSingleAnchor
    (
      const std::string &PDB_ID_SITE, const std::string &PDB_ID_ANCHOR,
      const assemble::SSE &SSE_SITE, const assemble::SSE &SSE_ANCHOR,
      const storage::Set< biol::AtomType> &ATOM_TYPES,
      const util::ShPtr< quality::MeasureInterface> &SP_MEASURE
    ) const
    {
      util::ShPtr< biol::AASequence> sp_sse_site_seq( util::CloneToShPtr( SSE_SITE));
      sp_sse_site_seq->SetFastaHeader( PDB_ID_SITE);

      util::ShPtr< biol::AASequence> sp_sse_anchor_seq( util::CloneToShPtr( SSE_ANCHOR));
      sp_sse_anchor_seq->SetFastaHeader( PDB_ID_ANCHOR);

      util::ShPtr< align::AlignmentInterface< biol::AABase> > sp_sse_site_align(   new align::AlignmentLeaf< biol::AABase>( sp_sse_site_seq));
      util::ShPtr< align::AlignmentInterface< biol::AABase> > sp_sse_anchor_align( new align::AlignmentLeaf< biol::AABase>( sp_sse_anchor_seq));

      // align first sses
      align::AlignerShift< biol::AABase> aligner;
      aligner.SetMinNumberAssignmentsWithoutGap( SSE_SITE.GetType()->GetFragmentLength());
      util::ShPtr< score::AlignmentQuality> sp_score( new score::AlignmentQuality());
      sp_score->SetAtoms( ATOM_TYPES);
      sp_score->SetMeasure( SP_MEASURE);
      aligner.SetScoreComparisonFunction( SP_MEASURE->GetComparisonFunction());
      aligner.SetScoringFunction( sp_score);

      const storage::Pair< align::AlignmentNode< biol::AABase>, double> alignment
      (
        aligner.AlignPair
        (
          sp_sse_site_align, // left side top
          sp_sse_anchor_align   // right side top
        )
      );

      return alignment.First();
    }

    //! @brief collect all possible pairs of anchor SSEs in a given protein, that has the same sse types as the site
    //!        sses and store the in a list per chain
    //! @param PROTEIN the model to search for SSEs of the same type as the site sses
    //! @param FIRST_SSE_SITE first sse of the site
    //! @param LAST_SSE_SITE second sse of the site
    //! @param SP_LOCATOR_ANCHOR_FIRST if the locator is defined, use only sses that match that locator
    //! @param SP_LOCATOR_ANCHOR_LAST if the locator is defined, use only sses that match that locator
    //! @return a map of chain ids with list of pairs of matching sses (according to type)
    storage::Map< char, storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > > >
    FusionProtein::CollectAnchorSSEs
    (
      const assemble::ProteinModel &PROTEIN,
      const assemble::SSE &FIRST_SSE_SITE,
      const assemble::SSE &LAST_SSE_SITE,
      const util::ShPtr< assemble::LocatorAA> &SP_LOCATOR_ANCHOR_FIRST,
      const util::ShPtr< assemble::LocatorAA> &SP_LOCATOR_ANCHOR_LAST
    )
    {
      storage::Map< char, storage::List< storage::VectorND< 2, util::SiPtr< const assemble::SSE> > > > anchor_sses;

      for
      (
        assemble::ProteinModel::const_iterator
          itr_chain( PROTEIN.GetChains().Begin()), itr_chain_end( PROTEIN.GetChains().End());
        itr_chain != itr_chain_end;
        ++itr_chain
      )
      {
        if
        (
             SP_LOCATOR_ANCHOR_FIRST.IsDefined()
          && SP_LOCATOR_ANCHOR_FIRST->GetLocatorChain().GetChainID() != ( *itr_chain)->GetChainID()
        )
        {
          continue;
        }

        // iterate over sses within the chains
        for
        (
          assemble::Chain::const_ierator
            itr_sse1( ( *itr_chain)->GetData().Begin()), itr_sse_end( ( *itr_chain)->GetData().End());
          itr_sse1 != itr_sse_end;
          ++itr_sse1
        )
        {
          if( !biol::GetSSTypes().AreSimilar( ( *itr_sse1)->GetType(), FIRST_SSE_SITE.GetType()))
          {
            continue;
          }

          if( SP_LOCATOR_ANCHOR_FIRST.IsDefined())
          {
            if( !SP_LOCATOR_ANCHOR_FIRST->Locate( **itr_sse1).IsDefined())
            {
              BCL_MessageDbg
              (
                "missing: " + ( *itr_sse1)->GetIdentification() + "\nby: " +
                SP_LOCATOR_ANCHOR_FIRST->GetIdentification()
              );
              continue;
            }
            else
            {
              BCL_MessageVrb
              (
                "considering: " + ( *itr_sse1)->GetIdentification() + "\nby: " +
                SP_LOCATOR_ANCHOR_FIRST->GetIdentification()
              );
            }
          }

          //              // for debugging
          //              if( ( *itr_sse1)->GetFirstAA()->GetSeqID() != first_sse_site->GetFirstAA()->GetSeqID()) continue;

          // iterate over sses within the chains
          assemble::Chain::const_ierator itr_sse2( itr_sse1);
          ++itr_sse2;
          for( ; itr_sse2 != itr_sse_end; ++itr_sse2)
          {
            if( !biol::GetSSTypes().AreSimilar( ( *itr_sse2)->GetType(), LAST_SSE_SITE.GetType()))
            {
              continue;
            }

            if( SP_LOCATOR_ANCHOR_LAST.IsDefined())
            {
              if( !SP_LOCATOR_ANCHOR_LAST->Locate( **itr_sse2).IsDefined())
              {
                BCL_MessageDbg
                (
                  "missing: " + ( *itr_sse2)->GetIdentification() + "\nby: " +
                  SP_LOCATOR_ANCHOR_LAST->GetIdentification()
                );
                continue;
              }
              else
              {
                BCL_MessageVrb
                (
                  "and considering: " + ( *itr_sse2)->GetIdentification() + "\nby: " +
                  SP_LOCATOR_ANCHOR_LAST->GetIdentification()
                );
              }
            }

            anchor_sses[ ( *itr_chain)->GetChainID()].PushBack
            (
              storage::VectorND< 2, util::SiPtr< const assemble::SSE> >( **itr_sse1, **itr_sse2)
            );

          } // sse 2
        } // sse 1
      } // chains

      return anchor_sses;
    }

    align::AlignmentNode< biol::AABase> FusionProtein::AlignTwoAnchors
    (
      const std::string &PDB_ID_SITE, const std::string &PDB_ID_ANCHOR,
      const assemble::SSE &FIRST_SSE_SITE, const assemble::SSE &FIRST_SSE_ANCHOR,
      const assemble::SSE &LAST_SSE_SITE, const assemble::SSE &LAST_SSE_ANCHOR
    ) const
    {
      util::ShPtr< biol::AASequence> sp_first_sse_site_seq( util::CloneToShPtr( FIRST_SSE_SITE));
      sp_first_sse_site_seq->SetFastaHeader( PDB_ID_SITE);
      util::ShPtr< biol::AASequence> sp_last_sse_site_seq( util::CloneToShPtr( LAST_SSE_SITE));
      sp_last_sse_site_seq->SetFastaHeader( PDB_ID_SITE);

      util::ShPtr< biol::AASequence> sp_first_sse_anchor_seq( util::CloneToShPtr( FIRST_SSE_ANCHOR));
      sp_first_sse_anchor_seq->SetFastaHeader( PDB_ID_ANCHOR);
      util::ShPtr< biol::AASequence> sp_last_sse_anchor_seq( util::CloneToShPtr( LAST_SSE_ANCHOR));
      sp_last_sse_anchor_seq->SetFastaHeader( PDB_ID_ANCHOR);

      util::ShPtr< align::AlignmentInterface< biol::AABase> > sp_first_sse_site_align
      (
        new align::AlignmentLeaf< biol::AABase>( sp_first_sse_site_seq)
      );
      util::ShPtr< align::AlignmentInterface< biol::AABase> > sp_last_sse_site_align
      (
        new align::AlignmentLeaf< biol::AABase>( sp_last_sse_site_seq)
      );
      util::ShPtr< align::AlignmentInterface< biol::AABase> > sp_first_sse_anchor_align
      (
        new align::AlignmentLeaf< biol::AABase>( sp_first_sse_anchor_seq)
      );
      util::ShPtr< align::AlignmentInterface< biol::AABase> > sp_last_sse_anchor_align
      (
        new align::AlignmentLeaf< biol::AABase>( sp_last_sse_anchor_seq)
      );

      // align first sses
      align::AlignerShift< biol::AABase> aligner;
      aligner.SetMinNumberAssignmentsWithoutGap
      (
        std::min( FIRST_SSE_SITE.GetType()->GetFragmentLength(), LAST_SSE_SITE.GetType()->GetFragmentLength())
      );
      util::ShPtr< score::AlignmentQuality> sp_score( new score::AlignmentQuality());
      sp_score->SetAtoms( m_AtomTypes);
      sp_score->SetMeasure( m_SPMeasure);
      aligner.SetScoreComparisonFunction( m_SPMeasure->GetComparisonFunction());
      aligner.SetScoringFunction( sp_score);

      const storage::Pair< align::AlignmentNode< biol::AABase>, double> alignment
      (
        aligner.AlignPairPair
        (
          sp_first_sse_site_align,   // left side top
          sp_first_sse_anchor_align, // left side bottom
          sp_last_sse_site_align,    // right side top
          sp_last_sse_anchor_align,  // right side bottom
          FIRST_SSE_SITE.GetType()->GetFragmentLength(), LAST_SSE_SITE.GetType()->GetFragmentLength(),
          true
        )
      );

      return alignment.First();
    }

    FusionProtein::Result FusionProtein::Fuse
    (
      const biol::AASequence &SCAFFOLD,
      const biol::AASequence &DONOR,
      const math::TransformationMatrix3D &DONOR_TRANSFORMATION,
      const align::AlignmentInterface< biol::AABase> &ALIGNMENT
    ) const
    {
      // residues at fusion point
      util::SiPtr< const biol::AABase> sp_left_scaffold;
      util::SiPtr< const biol::AABase> sp_left_donor;
      util::SiPtr< const biol::AABase> sp_right_donor;
      util::SiPtr< const biol::AABase> sp_right_scaffold;
      util::SiPtr< const biol::AABase> sp_left_focus;
      util::SiPtr< const biol::AABase> sp_right_focus;
      util::SiPtr< const biol::AABase> sp_left_insertion;
      util::SiPtr< const biol::AABase> sp_right_insertion;

      std::string replaced_scaffold_seq;
      std::string inserted_donor_seq;

      Result new_sequence;

      if( ALIGNMENT.IsEmpty() || ALIGNMENT.GetDepth() != 2)
      {
        return new_sequence;
      }

      align::AlignmentInterface< biol::AABase>::const_iterator
        align_itr( ALIGNMENT.GetAssignments().Begin()), align_itr_end( ALIGNMENT.GetAssignments().End());

      const char scaf_chain_id( SCAFFOLD.GetChainID());
      biol::AASequence::const_iterator aa_scaf_itr( SCAFFOLD.Begin()), aa_scaf_itr_end( SCAFFOLD.End());

      // find first assignment
      for( ; align_itr != align_itr_end && !( *align_itr)->GetMembers().IsDefined(); ++align_itr);

      if( align_itr == align_itr_end)
      {
        return new_sequence;
      }

      util::SiPtr< const biol::AABase> current_scaffold_align_aa( ( *align_itr)->GetMembers().FirstElement());
      const char new_chain_id( current_scaffold_align_aa->GetChainID());
      new_sequence.m_Sequence.SetChainID( new_chain_id);

      new_sequence.m_PymolScript += std::string( "show cartoon, /fused_model//") + new_chain_id + '\n';
      new_sequence.m_PymolScript += std::string( "color scaffold_color, /fused_model//")    + new_chain_id + '\n';
      new_sequence.m_PymolScript += std::string( "show cartoon, /scaffold_model//") + new_chain_id + '\n';
      new_sequence.m_PymolScript += std::string( "color scaffold_color, /scaffold_model//") + new_chain_id + '\n';

      // add all amino acids until start of the alignment
      for( ; aa_scaf_itr != aa_scaf_itr_end && ( *aa_scaf_itr)->GetData() != current_scaffold_align_aa->GetData(); ++aa_scaf_itr)
      {
        new_sequence.m_Sequence.PushBack( **aa_scaf_itr);
      }

      // ignore all amio acids from donor until the start of the alignment
      const char donor_chain_id( DONOR.GetChainID());
      biol::AASequence::const_iterator aa_don_itr( DONOR.Begin()), aa_don_itr_end( DONOR.End());
      util::SiPtr< const biol::AABase> current_donor_align_aa( *( ++( *align_itr)->GetMembers().Begin()));
      new_sequence.m_PymolScript += std::string( "show cartoon, /donor1_model//") + donor_chain_id + '/' + '\n';
      new_sequence.m_PymolScript += std::string( "color ignore_color, /donor1_model//") + donor_chain_id + '/' + '\n';

      while( aa_don_itr != aa_don_itr_end && ( *aa_don_itr)->GetData() != current_donor_align_aa->GetData())
      {
        ++aa_don_itr;
      }

      // color overlapping region before the replacement region
      {
        util::SiPtr< const biol::AABase> sp_last_scaf;
        util::SiPtr< const biol::AABase> sp_last_donor;

        std::string color_fused_overlap( std::string( "color overlap_color, /fused_model//")    + scaf_chain_id  + '/' + util::Format()( ( *aa_scaf_itr)->GetSeqID()) + '-');
        std::string color_scaf_overlap(  std::string( "color overlap_color, /scaffold_model//") + scaf_chain_id  + '/' + util::Format()( ( *aa_scaf_itr)->GetSeqID()) + '-');
        std::string color_donor_overlap( std::string( "color overlap_color, /donor1_model//")   + donor_chain_id + '/' + util::Format()( ( *aa_don_itr)->GetSeqID()) + '-');

        // add all amino acids from scaffold for aligned residues at first SSE
        while( align_itr != align_itr_end && aa_scaf_itr != aa_scaf_itr_end && aa_don_itr != aa_don_itr_end && ( *align_itr)->GetMembers().IsDefined())
        {
          // are the residues in alignment matching the scaffold and donor residues
          current_scaffold_align_aa = ( *align_itr)->GetMembers().FirstElement();
          current_donor_align_aa = *( ++( *align_itr)->GetMembers().Begin());

          // as long the scaffold and donor amino acid in the alignment agrees with the amino acid in the sequence
          if( ( *aa_scaf_itr)->GetData() != current_scaffold_align_aa->GetData() || ( *aa_don_itr)->GetData() != current_donor_align_aa->GetData())
          {
            break;
          }

          // update the first aa of the alignment for the pymol focus
          if( !sp_left_focus.IsDefined())
          {
            sp_left_focus = current_scaffold_align_aa;
          }

          //          BCL_MessageStd( "adding scaf align: " + ( *aa_scaf_itr)->GetIdentification());
          new_sequence.m_Sequence.PushBack( **aa_scaf_itr);
          sp_last_scaf = util::SiPtr< const biol::AABase>( *aa_scaf_itr);
          sp_last_donor = util::SiPtr< const biol::AABase>( *aa_don_itr);

          ++aa_scaf_itr;
          ++aa_don_itr;
          ++align_itr;
        }

        color_fused_overlap += util::Format()( sp_last_scaf->GetSeqID())  + '\n';
        color_scaf_overlap  += util::Format()( sp_last_scaf->GetSeqID())  + '\n';
        color_donor_overlap += util::Format()( sp_last_donor->GetSeqID()) + '\n';

        new_sequence.m_PymolScript += color_fused_overlap;
        new_sequence.m_PymolScript += color_scaf_overlap;
        new_sequence.m_PymolScript += color_donor_overlap;
      }

      // move align itr to next aligned scaffold aa
      while( align_itr != align_itr_end && !( *align_itr)->GetMembers().IsDefined())
      {
        ++align_itr;
      }
      BCL_Assert( align_itr != align_itr_end, "unexpected end of alignment");

      current_donor_align_aa = *( ++( *align_itr)->GetMembers().Begin());

      // transform and add AAs from donor
      if( aa_don_itr != aa_don_itr_end)
      {
        sp_left_donor = util::SiPtr< const biol::AABase>( **aa_don_itr);
      }
      for( ; aa_don_itr != aa_don_itr_end && ( *aa_don_itr)->GetData() != current_donor_align_aa->GetData(); ++aa_don_itr)
      {
//          BCL_MessageStd( "adding donor: " + ( *aa_don_itr)->GetIdentification());
        const biol::AABase &inserted_aa( ElongateSequence( new_sequence.m_Sequence, **aa_don_itr, DONOR_TRANSFORMATION, 0));
        if( !sp_left_insertion.IsDefined())
        {
          sp_left_insertion = &inserted_aa;
        }
        inserted_donor_seq.push_back( inserted_aa.GetType()->GetOneLetterCode());
        sp_right_donor = util::SiPtr< const biol::AABase>( **aa_don_itr);
        sp_right_insertion = &inserted_aa;
      }
      new_sequence.m_PymolScript += std::string( "color donor_color, /fused_model//") + new_chain_id + '/' +
                      util::Format()( sp_left_insertion->GetSeqID()) + '-' +
                      util::Format()( sp_right_insertion->GetSeqID()) + '\n';
      new_sequence.m_PymolScript += std::string( "color donor_color, /donor1_model//") + new_chain_id + '/' +
                      util::Format()( sp_left_donor->GetSeqID())      + '-' +
                      util::Format()( sp_right_donor->GetSeqID())     + '\n';

      // discard AAs from the scaffold
      current_scaffold_align_aa = ( *align_itr)->GetMembers().FirstElement();
      if( aa_scaf_itr != aa_scaf_itr_end)
      {
        sp_left_scaffold = util::SiPtr< const biol::AABase>( **aa_scaf_itr);
      }
      for( ; aa_scaf_itr != aa_scaf_itr_end && ( *aa_scaf_itr)->GetData() != current_scaffold_align_aa->GetData(); ++aa_scaf_itr)
      {
        new_sequence.m_PymolScript += std::string( "color ignore_color, /scaffold_model//") + scaf_chain_id + '/' +
                        util::Format()( ( *aa_scaf_itr)->GetSeqID()) + '\n';
        replaced_scaffold_seq.push_back( ( *aa_scaf_itr)->GetType()->GetOneLetterCode());
        sp_right_scaffold = util::SiPtr< const biol::AABase>( **aa_scaf_itr);
      }

      // add all amino acids from donor for aligned residues at last SSE
      while
      (
           align_itr != align_itr_end && aa_scaf_itr != aa_scaf_itr_end && aa_don_itr != aa_don_itr_end
        && ( *align_itr)->GetMembers().IsDefined()
      )
      {
        // are the residues in alignment matching the scaffold and donor residues
        current_scaffold_align_aa = ( *align_itr)->GetMembers().FirstElement();
        current_donor_align_aa = *( ++( *align_itr)->GetMembers().Begin());

        // as long the scaffold and donor amino acid in the alignment agrees with the amino acid in the sequence
        if
        (
             ( *aa_scaf_itr)->GetData() != current_scaffold_align_aa->GetData()
          || ( *aa_don_itr)->GetData() != current_donor_align_aa->GetData()
        )
        {
          break;
        }
        const biol::AABase &inserted_aa( ElongateSequence( new_sequence.m_Sequence, **aa_don_itr, DONOR_TRANSFORMATION, 0));

        new_sequence.m_PymolScript += std::string( "color overlap_color, /fused_model//")    + new_chain_id   + '/' +
                        util::Format()( inserted_aa.GetSeqID()) + "\n";
        new_sequence.m_PymolScript += std::string( "color overlap_color, /scaffold_model//") + scaf_chain_id  + '/' +
                        util::Format()( ( *aa_scaf_itr)->GetSeqID()) + "\n";
        new_sequence.m_PymolScript += std::string( "color overlap_color, /donor1_model//")   + donor_chain_id + '/' +
                        util::Format()( ( *aa_don_itr)->GetSeqID()) + "\n";

        // update the last aa of the alignment for pymol focus
        sp_right_focus = util::SiPtr< const biol::AABase>( inserted_aa);

        ++aa_scaf_itr;
        ++aa_don_itr;
        ++align_itr;
      }

      // color rest of scaffold
      for( ; aa_scaf_itr != aa_scaf_itr_end; ++aa_scaf_itr)
      {
        new_sequence.m_PymolScript += std::string( "color ignore_color, /scaffold_model//") + scaf_chain_id + '/' +
                        util::Format()( ( *aa_scaf_itr)->GetSeqID()) + '\n';
      }

      // add the rest from the donor
      while( aa_don_itr != aa_don_itr_end)
      {
        const biol::AABase &inserted_aa( ElongateSequence( new_sequence.m_Sequence, **aa_don_itr, DONOR_TRANSFORMATION, 0));

        new_sequence.m_PymolScript += std::string( "color donor_color, /donor1_model//") + donor_chain_id + '/' +
                        util::Format()( ( *aa_don_itr)->GetSeqID()) + '\n';
        new_sequence.m_PymolScript += std::string( "color donor_color, /fused_model//")  + new_chain_id   + '/' +
                        util::Format()( inserted_aa.GetSeqID()) + '\n';
        ++aa_don_itr;
      }

      new_sequence.m_PymolScript += "bg_color white\n";
      new_sequence.m_PymolScript += "disable all\n";

      // fused model
      new_sequence.m_PymolScript += "enable fused_model\n";
      if( sp_left_focus.IsDefined() && sp_right_focus.IsDefined())
      {
        // zoom onto replacement
        new_sequence.m_PymolScript += std::string( "orient /fused_model//") + new_chain_id + '/' +
                        util::Format()( sp_left_focus->GetSeqID()) + '-' +
                        util::Format()( sp_right_focus->GetSeqID()) + '\n';
        new_sequence.m_PymolScript += "rotate y, 90\n";
        new_sequence.m_PymolScript += "ray\n";
        new_sequence.m_PymolScript += "png fused_site.png, dpi=300\n";
      }

      new_sequence.m_PymolScript += "orient fused_model\n";
      new_sequence.m_PymolScript += "ray\n";
      new_sequence.m_PymolScript += "png fused_model.png, dpi=300\n";

      new_sequence.m_PymolScript += "disable all\n";
      new_sequence.m_PymolScript += "enable donor1_model\n";
      new_sequence.m_PymolScript += "ray\n";
      new_sequence.m_PymolScript += "png donor1_model.png, dpi=300\n";

      new_sequence.m_PymolScript += "disable all\n";
      new_sequence.m_PymolScript += "enable scaffold_model\n";
      new_sequence.m_PymolScript += "ray\n";
      new_sequence.m_PymolScript += "png scaffold_model.png, dpi=300\n";

      new_sequence.m_PymolScript += "system convert -pointsize 32 -delay 0 scaffold_model.png -delay 100 label:'scaffold' -delay 0 "
                      "donor1_model.png -delay 100 label:'donor' -delay 0 fused_model.png -delay 100 label:'fused' "
                      "-delay 0 fused_site.png -delay 100 label:'site' animation.gif\n";

      new_sequence.m_PymolScript += "enable all\n";

//        size_t closure_count( 1);
//        PYMOL_SCRIPT += PymolClosure( sp_left_insertion->GetChainID(), sp_left_insertion->GetSeqID() - 1, sp_left_insertion->GetSeqID(), 'l', closure_count);
//        PYMOL_SCRIPT += PymolClosure( sp_right_insertion->GetChainID(), sp_right_insertion->GetSeqID(), sp_right_insertion->GetSeqID() + 1, 'r', closure_count);

      BCL_MessageStd( "sequence replaced in scaffold: " + replaced_scaffold_seq);
      BCL_MessageStd( "sequence inserted from donor:  " + inserted_donor_seq);

      std::string replace_info;

      // scaffold
      if( sp_left_scaffold.IsDefined())
      {
        BCL_MessageStd
        (
          "replaced in scaffold first: " + sp_left_scaffold->GetIdentification() + ' ' +
          util::Format()( sp_left_scaffold->GetPdbID())
        );
        replace_info += sp_left_scaffold->GetType()->GetThreeLetterCode() + '_' + scaf_chain_id +
                        '_' + util::Format()( sp_left_scaffold->GetSeqID()) + '_' +
                        util::Format()( sp_left_scaffold->GetPdbID()) + '_';
      }
      if( sp_right_scaffold.IsDefined())
      {
        BCL_MessageStd
        (
          "replaced in scaffold last:  " + sp_right_scaffold->GetIdentification() + ' ' +
          util::Format()( sp_right_scaffold->GetPdbID())
        );
        replace_info += sp_right_scaffold->GetType()->GetThreeLetterCode() + '_' + scaf_chain_id +
                        '_' + util::Format()( sp_right_scaffold->GetSeqID()) + '_' +
                        util::Format()( sp_right_scaffold->GetPdbID()) + '_';
      }
      replace_info += replaced_scaffold_seq + '_';

      // donor
      if( sp_left_donor.IsDefined())
      {
        BCL_MessageStd
        (
          "took from donor start: " + sp_left_donor->GetIdentification() + ' ' +
          util::Format()( sp_left_donor->GetPdbID())
        );
        replace_info += sp_left_donor->GetType()->GetThreeLetterCode() + '_' + donor_chain_id + '_' +
                        util::Format()( sp_left_donor->GetSeqID()) + '_' +
                        util::Format()( sp_left_donor->GetPdbID()) + '_';
      }
      if( sp_right_donor.IsDefined())
      {
        BCL_MessageStd
        (
          "took form donor last:  " + sp_right_donor->GetIdentification() + ' ' +
          util::Format()( sp_right_donor->GetPdbID())
        );
        replace_info += sp_right_donor->GetType()->GetThreeLetterCode() + '_' + donor_chain_id + '_' +
                        util::Format()( sp_right_donor->GetSeqID()) + '_' +
                        util::Format()( sp_right_donor->GetPdbID()) + '_';
      }
      replace_info += inserted_donor_seq;
      replace_info += '_';

      new_sequence.m_Sequence.SetFastaHeader( replace_info);

      return new_sequence;
    }

    //! @brief
    //! @param RESULT the result contains the fused sequence so far and will be elongated
    //! @param ALIGN_ITR alignment iterator at the starting position
    //! @param ALIGN_ITR_SWITCH itr where the fused sequence switches from donor to scaffold
    //! @param ALIGN_ITR_END end of alignment
    //! @param AA_SCAF_ITR current pos in scaffold sequence
    //! @param AA_SCAF_ITR_END end of scaffold sequence
    //! @param AA_DON_ITR current pos in donor sequence
    //! @param AA_DON_ITR_END current pos in scaffold sequence
    //! @param DONOR_NAME name of donor in pymol script
    //! @param DONOR_SITE_NAME name of site within fused sequence
    //! @param TRANSFORMATION transformation to apply to donor amino acids
    //! @param COUNT the cutpoint number
    //! @param SEQ_POS sequence position, if RESULT is empty, but the sequence does not start at 0
    //! @return a new result with elongated sequence and additional information
    FusionProtein::Result FusionProtein::HandleCutpointDonorScaffold
    (
      const FusionProtein::Result &RESULT,
      const align::AlignmentInterface< biol::AABase>::const_iterator &ALIGN_ITR,
      const align::AlignmentInterface< biol::AABase>::const_iterator &ALIGN_ITR_SWITCH,
      const align::AlignmentInterface< biol::AABase>::const_iterator &ALIGN_ITR_END,
      const biol::AASequence::const_iterator &AA_SCAF_ITR,
      const biol::AASequence::const_iterator &AA_SCAF_ITR_END,
      const biol::AASequence::const_iterator &AA_DON_ITR,
      const biol::AASequence::const_iterator &AA_DON_ITR_END,
      const std::string &DONOR_NAME,
      const std::string &DONOR_SITE_NAME,
      const math::TransformationMatrix3D &TRANSFORMATION,
      const size_t &COUNT,
      const size_t SEQ_POS
    ) const
    {
      Result fused_sequence( RESULT);

      align::AlignmentInterface< biol::AABase>::const_iterator align_itr( ALIGN_ITR);

      biol::AASequence::const_iterator aa_scaf_itr( AA_SCAF_ITR), aa_don_itr( AA_DON_ITR);

      util::SiPtr< const biol::AABase> current_scaffold_align_aa( ( *align_itr)->GetMembers().FirstElement());
      util::SiPtr< const biol::AABase> current_donor_align_aa( *( ++( *align_itr)->GetMembers().Begin()));

      const char new_chain_id( fused_sequence.m_Sequence.GetChainID());
      const char scaffold_chain_id( current_scaffold_align_aa->GetChainID());
      const char donor_chain_id( current_donor_align_aa->GetChainID());

      std::string replaced_scaffold_seq;
      std::string inserted_donor_seq;

      util::SiPtr< const biol::AABase> sp_scaffold_ignored_begin;
      util::SiPtr< const biol::AABase> sp_scaffold_ignored_last;
      util::SiPtr< const biol::AABase> sp_donor_used_begin;
      util::SiPtr< const biol::AABase> sp_donor_used_last;

      util::SiPtr< const biol::AABase> sp_scaffold_used_begin;
      util::SiPtr< const biol::AABase> sp_scaffold_used_last;
      util::SiPtr< const biol::AABase> sp_donor_ignored_begin;
      util::SiPtr< const biol::AABase> sp_donor_ignored_last;

      util::SiPtr< const biol::AABase> sp_fused_begin;
      util::SiPtr< const biol::AABase> sp_fused_last;
      util::SiPtr< const biol::AABase> sp_fused_right;

      // add all amino acids from donor for aligned residues at first SSE
      while
      (
           align_itr != ALIGN_ITR_SWITCH && aa_scaf_itr != AA_SCAF_ITR_END && aa_don_itr != AA_DON_ITR_END
        && ( *align_itr)->GetMembers().IsDefined()
      )
      {
        // are the residues in alignment matching the scaffold and donor residues
        current_scaffold_align_aa = ( *align_itr)->GetMembers().FirstElement();
        current_donor_align_aa = *( ++( *align_itr)->GetMembers().Begin());

        // as long the scaffold and donor amino acid in the alignment agrees with the amino acid in the sequence
        if
        (
             ( *aa_scaf_itr)->GetData() != current_scaffold_align_aa->GetData()
          || ( *aa_don_itr)->GetData() != current_donor_align_aa->GetData()
        )
        {
          break;
        }

        fused_sequence.m_AlignFuseScore += m_SPAssignmentScore->operator()( **align_itr);
        ++fused_sequence.m_NrFusedResiduesCTerm;

        sp_fused_last = ElongateSequence( fused_sequence.m_Sequence, **aa_don_itr, TRANSFORMATION, SEQ_POS);

        // add entry to rosetta res file
        fused_sequence.m_RosettaResFile +=
            RosettaResFileEntry( *sp_fused_last->GetData(), *( *aa_scaf_itr)->GetData(), *( *aa_don_itr)->GetData());

        if( !sp_fused_begin.IsDefined())
        {
          sp_fused_begin = sp_fused_last;
        }

        if( !sp_donor_used_begin.IsDefined())
        {
          sp_donor_used_begin = *aa_don_itr;
        }
        sp_donor_used_last = *aa_don_itr;

        if( !sp_scaffold_ignored_begin.IsDefined())
        {
          sp_scaffold_ignored_begin = *aa_scaf_itr;
        }
        sp_scaffold_ignored_last = *aa_scaf_itr;

        inserted_donor_seq.push_back( sp_fused_last->GetType()->GetOneLetterCode());
        replaced_scaffold_seq.push_back( ( *aa_scaf_itr)->GetType()->GetOneLetterCode());

        ++aa_scaf_itr;
        ++aa_don_itr;
        ++align_itr;
      }

      fused_sequence.m_PymolScript += "#cutpoint donor->scaffold BEGIN\n";
      if( sp_donor_used_begin.IsDefined())
      {
        // replace information
        fused_sequence.m_LastDonorSeqID = sp_donor_used_last->GetSeqID();
        fused_sequence.m_ReplaceInformation += '-' + util::Format()( fused_sequence.m_LastDonorSeqID);
        fused_sequence.m_ReplaceInformation += "/scaf:" + util::Format()( sp_scaffold_ignored_last->GetSeqID() + 1);

        fused_sequence.m_PymolScript += "color donor_color, /" + DONOR_NAME + "//" + donor_chain_id + '/'
            + util::Format()( sp_donor_used_begin->GetSeqID()) + '-'
            + util::Format()( sp_donor_used_last->GetSeqID()) + '\n';

        fused_sequence.m_PymolScript += std::string( "color ignore_color, /scaffold_model//") + scaffold_chain_id + '/'
            + util::Format()( sp_scaffold_ignored_begin->GetSeqID()) + '-'
            + util::Format()( sp_scaffold_ignored_last->GetSeqID()) + '\n';
        // selection for site in fused model
        fused_sequence.m_PymolScript += "select " + DONOR_SITE_NAME + ", " + DONOR_SITE_NAME + " | /fused_model//"
            + new_chain_id + '/'
            + util::Format()( sp_fused_begin->GetSeqID()) + '-'
            + util::Format()( sp_fused_last->GetSeqID()) + '\n';

        fused_sequence.m_PymolScript += "color donor_color, " + DONOR_SITE_NAME + '\n';
      }
      // no donor aa was added anymore
      else if( aa_don_itr != AA_DON_ITR_END)
      {
        fused_sequence.m_LastDonorSeqID = ( *aa_don_itr)->GetSeqID() - 1;
        fused_sequence.m_ReplaceInformation += '-' + util::Format()( fused_sequence.m_LastDonorSeqID);
        if( aa_scaf_itr != AA_SCAF_ITR_END)
        {
          fused_sequence.m_ReplaceInformation += "/scaf:" + util::Format()( ( *aa_scaf_itr)->GetSeqID());
        }
      }

      // add all amino acids from scaffold for aligned residues at first SSE
      fused_sequence.m_RosettaResFile += "# donor->scaffold " + DONOR_SITE_NAME + '\n';
      while
      (
           align_itr != ALIGN_ITR_END && aa_scaf_itr != AA_SCAF_ITR_END && aa_don_itr != AA_DON_ITR_END
        && ( *align_itr)->GetMembers().IsDefined()
      )
      {
        // are the residues in alignment matching the scaffold and donor residues
        current_scaffold_align_aa = ( *align_itr)->GetMembers().FirstElement();
        current_donor_align_aa    = *( ++( *align_itr)->GetMembers().Begin());

        // as long the scaffold and donor amino acid in the alignment agrees with the amino acid in the sequence
        if
        (
             ( *aa_scaf_itr)->GetData() != current_scaffold_align_aa->GetData()
          || ( *aa_don_itr)->GetData() != current_donor_align_aa->GetData()
        )
        {
          break;
        }

        const biol::AABase &inserted_aa( ElongateSequence( fused_sequence.m_Sequence, **aa_scaf_itr, SEQ_POS));

        // add entry to rosetta res file
        fused_sequence.m_RosettaResFile +=
            RosettaResFileEntry( *inserted_aa.GetData(), *( *aa_scaf_itr)->GetData(), *( *aa_don_itr)->GetData());

        if( !sp_fused_right.IsDefined())
        {
          sp_fused_right = &inserted_aa;
        }

        if( !sp_scaffold_used_begin.IsDefined())
        {
          sp_scaffold_used_begin = *aa_scaf_itr;
        }
        sp_scaffold_used_last = *aa_scaf_itr;

        if( !sp_donor_ignored_begin.IsDefined())
        {
          sp_donor_ignored_begin = *aa_don_itr;
        }
        sp_donor_ignored_last = *aa_don_itr;

        ++aa_scaf_itr;
        ++aa_don_itr;
        ++align_itr;
      }

      if( sp_scaffold_used_begin.IsDefined())
      {
        fused_sequence.m_PymolScript += std::string( "color overlap_color, /scaffold_model//") + scaffold_chain_id + '/'
           + util::Format()( sp_scaffold_used_begin->GetSeqID()) + '-'
           + util::Format()( sp_scaffold_used_last->GetSeqID()) + '\n';
        fused_sequence.m_PymolScript += std::string( "color overlap_color, /") + DONOR_NAME + "//" + donor_chain_id
            + '/' + util::Format()( sp_donor_ignored_begin->GetSeqID()) + '-'
            + util::Format()( sp_donor_ignored_last->GetSeqID()) + '\n';
      }

      if( sp_fused_last.IsDefined() && sp_fused_right.IsDefined())
      {
        fused_sequence.m_PymolTailScript += PymolClosure( *sp_fused_last, *sp_fused_right, 'r', COUNT);
      }

      BCL_MessageStd( "sequence replaced in scaffold: " + replaced_scaffold_seq);
      BCL_MessageStd( "sequence inserted from donor:  " + inserted_donor_seq);
      fused_sequence.m_NrReplacedResidues += replaced_scaffold_seq.size();

      fused_sequence.m_PymolScript += "#cutpoint donor->scaffold END\n";

      return fused_sequence;
    }

    //! @brief
    //! @param RESULT the result contains the fused sequence so far and will be elongated
    //! @param ALIGN_ITR alignment iterator at the starting position
    //! @param ALIGN_ITR_SWITCH itr where the fused sequence switches from donor to scaffold
    //! @param ALIGN_ITR_END end of alignment
    //! @param AA_SCAF_ITR current pos in scaffold sequence
    //! @param AA_SCAF_ITR_END end of scaffold sequence
    //! @param AA_DON_ITR current pos in donor sequence
    //! @param AA_DON_ITR_END current pos in scaffold sequence
    //! @param DONOR_NAME name of donor in pymol script
    //! @param DONOR_SITE_NAME name of site within fused sequence
    //! @param TRANSFORMATION transformation to apply to donor amino acids
    //! @param COUNT the cutpoint number
    //! @param SEQ_POS sequence position, if RESULT is empty, but the sequence does not start at 0
    //! @return a new result with elongated sequence and additional information
    FusionProtein::Result FusionProtein::HandleCutpointScaffoldDonor
    (
      const FusionProtein::Result &RESULT,
      const align::AlignmentInterface< biol::AABase>::const_iterator &ALIGN_ITR,
      const align::AlignmentInterface< biol::AABase>::const_iterator &ALIGN_ITR_SWITCH,
      const align::AlignmentInterface< biol::AABase>::const_iterator &ALIGN_ITR_END,
      const biol::AASequence::const_iterator &AA_SCAF_ITR,
      const biol::AASequence::const_iterator &AA_SCAF_ITR_END,
      const biol::AASequence::const_iterator &AA_DON_ITR,
      const biol::AASequence::const_iterator &AA_DON_ITR_END,
      const std::string &DONOR_NAME,
      const std::string &DONOR_SITE_NAME,
      const math::TransformationMatrix3D &TRANSFORMATION,
      const size_t &COUNT,
      const size_t SEQ_POS
    ) const
    {
      Result fused_sequence( RESULT);

      align::AlignmentInterface< biol::AABase>::const_iterator align_itr( ALIGN_ITR);

      biol::AASequence::const_iterator aa_scaf_itr( AA_SCAF_ITR), aa_don_itr( AA_DON_ITR);

      util::SiPtr< const biol::AABase> current_scaffold_align_aa( ( *align_itr)->GetMembers().FirstElement());
      util::SiPtr< const biol::AABase> current_donor_align_aa( *( ++( *align_itr)->GetMembers().Begin()));

      const char new_chain_id( fused_sequence.m_Sequence.GetChainID());
      const char scaffold_chain_id( current_scaffold_align_aa->GetChainID());
      const char donor_chain_id( current_donor_align_aa->GetChainID());

      std::string replaced_scaffold_seq;
      std::string inserted_donor_seq;

      util::SiPtr< const biol::AABase> sp_scaffold_used_begin;
      util::SiPtr< const biol::AABase> sp_scaffold_used_last;
      util::SiPtr< const biol::AABase> sp_donor_ignored_begin;
      util::SiPtr< const biol::AABase> sp_donor_ignored_last;

      util::SiPtr< const biol::AABase> sp_scaffold_ignored_begin;
      util::SiPtr< const biol::AABase> sp_scaffold_ignored_last;
      util::SiPtr< const biol::AABase> sp_donor_used_begin;
      util::SiPtr< const biol::AABase> sp_donor_used_last;

      util::SiPtr< const biol::AABase> sp_fused_left;
      util::SiPtr< const biol::AABase> sp_fused_begin;
      util::SiPtr< const biol::AABase> sp_fused_last;

      fused_sequence.m_RosettaResFile += "# scaffold->donor " + DONOR_SITE_NAME + '\n';
      while
      (
           align_itr != ALIGN_ITR_SWITCH && aa_scaf_itr != AA_SCAF_ITR_END && aa_don_itr != AA_DON_ITR_END
        && ( *align_itr)->GetMembers().IsDefined()
      )
      {
        // are the residues in alignment matching the scaffold and donor residues
        current_scaffold_align_aa = ( *align_itr)->GetMembers().FirstElement();
        current_donor_align_aa = *( ++( *align_itr)->GetMembers().Begin());

        // as long the scaffold and donor amino acid in the alignment agrees with the amino acid in the sequence
        if
        (
             ( *aa_scaf_itr)->GetData() != current_scaffold_align_aa->GetData()
          || ( *aa_don_itr)->GetData() != current_donor_align_aa->GetData()
        )
        {
          break;
        }

        sp_fused_left = ElongateSequence( fused_sequence.m_Sequence, **aa_scaf_itr, SEQ_POS);

        // add entry to rosetta res file
        fused_sequence.m_RosettaResFile +=
            RosettaResFileEntry( *sp_fused_left->GetData(), *( *aa_scaf_itr)->GetData(), *( *aa_don_itr)->GetData());

        if( !sp_scaffold_used_begin.IsDefined())
        {
          sp_scaffold_used_begin = *aa_scaf_itr;
        }
        sp_scaffold_used_last = *aa_scaf_itr;

        if( !sp_donor_ignored_begin.IsDefined())
        {
          sp_donor_ignored_begin = *aa_don_itr;
        }
        sp_donor_ignored_last = *aa_don_itr;

        ++aa_scaf_itr;
        ++aa_don_itr;
        ++align_itr;
      }

      fused_sequence.m_PymolScript += "#cutpoint scaffold->donor BEGIN\n";
      if( sp_scaffold_used_begin.IsDefined())
      {
        // replace information
        fused_sequence.m_LastScaffoldSeqID = sp_scaffold_used_last->GetSeqID();
        fused_sequence.m_ReplaceInformation += '-' + util::Format()( fused_sequence.m_LastScaffoldSeqID);
        fused_sequence.m_ReplaceInformation += "\\don:" + util::Format()( sp_donor_ignored_last->GetSeqID() + 1);

        fused_sequence.m_PymolScript += std::string( "color overlap_color, /scaffold_model//") + scaffold_chain_id + '/'
           + util::Format()( sp_scaffold_used_begin->GetSeqID()) + '-'
           + util::Format()( sp_scaffold_used_last->GetSeqID()) + '\n';
        fused_sequence.m_PymolScript += std::string( "color overlap_color, /") + DONOR_NAME + "//" + donor_chain_id
            + '/' + util::Format()( sp_donor_ignored_begin->GetSeqID()) + '-'
            + util::Format()( sp_donor_ignored_last->GetSeqID()) + '\n';
      }
      // no scaffold aa was added anymore
      else if( aa_scaf_itr != AA_SCAF_ITR_END)
      {
        fused_sequence.m_LastScaffoldSeqID = ( *aa_scaf_itr)->GetSeqID() - 1;
        fused_sequence.m_ReplaceInformation += '-' + util::Format()( fused_sequence.m_LastScaffoldSeqID);
        if( aa_don_itr != AA_DON_ITR_END)
        {
          fused_sequence.m_ReplaceInformation += "\\don:" + util::Format()( ( *aa_don_itr)->GetSeqID());
        }
      }

      // add all amino acids from donor for aligned residues at first SSE
      while
      (
           align_itr != ALIGN_ITR_END && aa_scaf_itr != AA_SCAF_ITR_END && aa_don_itr != AA_DON_ITR_END
        && ( *align_itr)->GetMembers().IsDefined()
      )
      {
        // are the residues in alignment matching the scaffold and donor residues
        current_scaffold_align_aa = ( *align_itr)->GetMembers().FirstElement();
        current_donor_align_aa    = *( ++( *align_itr)->GetMembers().Begin());

        // as long the scaffold and donor amino acid in the alignment agrees with the amino acid in the sequence
        if
        (
             ( *aa_scaf_itr)->GetData() != current_scaffold_align_aa->GetData()
          || ( *aa_don_itr)->GetData() != current_donor_align_aa->GetData()
        )
        {
          break;
        }

        if( !sp_scaffold_ignored_begin.IsDefined())
        {
          sp_scaffold_ignored_begin = *aa_scaf_itr;
        }
        sp_scaffold_ignored_last = *aa_scaf_itr;

        if( !sp_donor_used_begin.IsDefined())
        {
          sp_donor_used_begin = *aa_don_itr;
        }
        sp_donor_used_last = *aa_don_itr;

        fused_sequence.m_AlignFuseScore += m_SPAssignmentScore->operator()( **align_itr);
        ++fused_sequence.m_NrFusedResiduesNTerm;

        sp_fused_last = ElongateSequence( fused_sequence.m_Sequence, **aa_don_itr, TRANSFORMATION, SEQ_POS);

        // add entry to rosetta res file
        fused_sequence.m_RosettaResFile +=
            RosettaResFileEntry( *sp_fused_last->GetData(), *( *aa_scaf_itr)->GetData(), *( *aa_don_itr)->GetData());

        // update the first aa of the alignment for the pymol focus
        if( !sp_fused_begin.IsDefined())
        {
          sp_fused_begin = sp_fused_last;
        }

        inserted_donor_seq.push_back( sp_fused_last->GetType()->GetOneLetterCode());
        replaced_scaffold_seq.push_back( ( *aa_scaf_itr)->GetType()->GetOneLetterCode());

        ++aa_scaf_itr;
        ++aa_don_itr;
        ++align_itr;
      }

      // selection onto replacement
      if( sp_donor_used_begin.IsDefined())
      {
        fused_sequence.m_PymolScript += "color donor_color, /" + DONOR_NAME + "//" + donor_chain_id + '/'
            + util::Format()( sp_donor_used_begin->GetSeqID()) + '-'
            + util::Format()( sp_donor_used_last->GetSeqID()) + '\n';

        fused_sequence.m_PymolScript += std::string( "color ignore_color, /scaffold_model//") + scaffold_chain_id + '/'
            + util::Format()( sp_scaffold_ignored_begin->GetSeqID()) + '-'
            + util::Format()( sp_scaffold_ignored_last->GetSeqID()) + '\n';
        // selection for site in fused model
        fused_sequence.m_PymolScript += "select " + DONOR_SITE_NAME + ", " + DONOR_SITE_NAME + " | /fused_model//"
            + new_chain_id + '/'
            + util::Format()( sp_fused_begin->GetSeqID()) + '-'
            + util::Format()( sp_fused_last->GetSeqID()) + '\n';

        fused_sequence.m_PymolScript += "color donor_color, " + DONOR_SITE_NAME + '\n';
      }

      if( sp_fused_left.IsDefined() && sp_fused_begin.IsDefined())
      {
        fused_sequence.m_PymolTailScript += PymolClosure( *sp_fused_left, *sp_fused_begin, 'l', COUNT);
      }

      BCL_MessageStd( "sequence replaced in scaffold: " + replaced_scaffold_seq);
      BCL_MessageStd( "sequence inserted from donor:  " + inserted_donor_seq);
      fused_sequence.m_NrReplacedResidues += replaced_scaffold_seq.size();

      fused_sequence.m_PymolScript += "#cutpoint scaffold->donor END\n";

      return fused_sequence;
    }

    storage::List< FusionProtein::Result> FusionProtein::ReplaceSystematic
    (
      const biol::AASequence &SCAFFOLD,
      const storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, util::SiPtr< const assemble::ProteinModel> > > &ALIGNMENTS
    ) const
    {
      const char scaf_chain_id( SCAFFOLD.GetChainID());
      const char new_chain_id( scaf_chain_id);

      biol::AASequence::const_iterator aa_scaf_itr( SCAFFOLD.Begin()), aa_scaf_itr_end( SCAFFOLD.End());
      size_t count( 1);

      // collect all fused sequences
      storage::List< storage::List< Result> > fused_sequence_segments;
      {
        Result start_sequence;

        start_sequence.m_PymolScript += std::string( "show cartoon, /fused_model//") + new_chain_id + '\n';
        start_sequence.m_PymolScript += std::string( "show cartoon, /scaffold_model//") + new_chain_id + '\n';
        start_sequence.m_PymolScript += std::string( "color scaffold_color, /fused_model//") + new_chain_id + '\n';
        start_sequence.m_PymolScript += std::string( "color scaffold_color, /scaffold_model//") + new_chain_id + '\n';

        if( aa_scaf_itr != aa_scaf_itr_end)
        {
          start_sequence.m_ReplaceInformation += "scaf:" + util::Format()( ( *aa_scaf_itr)->GetSeqID());
        }

        fused_sequence_segments.PushBack( storage::List< Result>( 1, start_sequence));
      }

      // keep track of the total sequence length (which is the same for all constructs) so that identical subsequences
      // that are generate as intermediate steps start with the proper seq id
      size_t seq_length( 0);

      // iterate through all alignments and donors
      for
      (
        storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, util::SiPtr< const assemble::ProteinModel> > >::const_iterator
          itr( ALIGNMENTS.Begin()), itr_end( ALIGNMENTS.End());
        itr != itr_end;
        ++itr, ++count
      )
      {
        const align::AlignmentNode< biol::AABase> &current_alignment( itr->First());
        if( current_alignment.IsEmpty() || current_alignment.GetDepth() != 2)
        {
          continue;
        }

        const std::string donor_name( "donor" + util::Format()( count) + "_model");
        const std::string donor_site_name( "donor" + util::Format()( count) + "_site");

        const storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > coord_pair
        (
          assemble::Quality::CoordinatesFromAlignment( current_alignment, m_AtomTypes)
        );
        const math::TransformationMatrix3D transformation( m_SPMeasure->CalculateSuperimposition( coord_pair.Second(), coord_pair.First()));

        // residues at fusion point
        util::SiPtr< const biol::AABase> sp_scaffold_ignore_begin;
        util::SiPtr< const biol::AABase> sp_scaffold_ignore_last;

        util::SiPtr< const biol::AABase> sp_donor_used_begin;
        util::SiPtr< const biol::AABase> sp_donor_used_last;

        util::SiPtr< const biol::AABase> sp_fused_begin;
        util::SiPtr< const biol::AABase> sp_fused_last;

        std::string replaced_scaffold_seq;
        std::string inserted_donor_seq;

        align::AlignmentInterface< biol::AABase>::const_iterator
          align_itr( current_alignment.GetAssignments().Begin()),
          align_itr_end( current_alignment.GetAssignments().End());

        // find first assignment
        for( ; align_itr != align_itr_end && !( *align_itr)->GetMembers().IsDefined(); ++align_itr);

        if( align_itr == align_itr_end)
        {
          continue;
        }

        const char donor_chain_id( ( *( ++current_alignment.GetSequences().Begin()))->GetFirstMember()->GetChainID());
        const biol::AASequence &donor_sequence( *itr->Second()->GetChain( donor_chain_id)->GetSequence());
        util::SiPtr< const biol::AABase> current_donor_align_aa( *( ++( *align_itr)->GetMembers().Begin()));

        // sequence before the alignment
        {
          Result new_sequence;
          new_sequence.m_Sequence.SetChainID( new_chain_id);
          new_sequence.m_PymolScript += "select " + donor_site_name + ", " + donor_name + " & scaffold_model\n";
          new_sequence.m_PymolScript += std::string( "show cartoon, /")       + donor_name + "//" + donor_chain_id + '/' + '\n';
          new_sequence.m_PymolScript += std::string( "color ignore_color, /") + donor_name + "//" + donor_chain_id + '/' + '\n';

          util::SiPtr< const biol::AABase> current_scaffold_align_aa( ( *align_itr)->GetMembers().FirstElement());

          // add all amino acids from scaffold to fused sequence until start of the alignment
          for
          (
            ; aa_scaf_itr != aa_scaf_itr_end && ( *aa_scaf_itr)->GetData() != current_scaffold_align_aa->GetData();
            ++aa_scaf_itr, ++seq_length
          )
          {
            ElongateSequence( new_sequence.m_Sequence, **aa_scaf_itr, seq_length);
          }
          fused_sequence_segments.PushBack( storage::List< Result>( 1, new_sequence));
          BCL_MessageVrb
          (
            "starting sequence before fusion " + util::Format()( count) + ':' +
            new_sequence.m_Sequence.GetSequenceIdentification()
          );
        }

        biol::AASequence::const_iterator aa_don_itr( donor_sequence.Begin()), aa_don_itr_end( donor_sequence.End());

        // discard donor amino acids until the start of the alignment
        while( aa_don_itr != aa_don_itr_end && ( *aa_don_itr)->GetData() != current_donor_align_aa->GetData())
        {
          ++aa_don_itr;
        }

        // for possible cut points from scaffold to donor, create a new sequence
        {
          storage::List< Result> fused_sequences;
          for
          (
            align::AlignmentInterface< biol::AABase>::const_iterator align_itr_switch( align_itr);
              align_itr_switch != align_itr_end && ( *align_itr_switch)->GetMembers().IsDefined();
            ++align_itr_switch
          )
          {
            Result fused_sequence;
            fused_sequence.m_Sequence.SetChainID( new_chain_id);

            fused_sequence = HandleCutpointScaffoldDonor
                (
                  fused_sequence,
                  align_itr, align_itr_switch, align_itr_end,
                  aa_scaf_itr, aa_scaf_itr_end,
                  aa_don_itr, aa_don_itr_end,
                  donor_name, donor_site_name,
                  transformation, count, seq_length
                );
            fused_sequences.PushBack( fused_sequence);
            BCL_MessageVrb( "sequence scaffold donor: " + util::Format()( count) + ':' + fused_sequence.m_Sequence.GetSequenceIdentification());
          }

          // add all fused sequences for possible cut points for this segment
          fused_sequence_segments.PushBack( fused_sequences);
          BCL_MessageVrb( "finished cutpoint scaffold donor " + util::Format()( count));

          // move align, donor and scaffold itr to the end of the alignment
          for
          (
            ; align_itr != align_itr_end && ( *align_itr)->GetMembers().IsDefined() &&
              aa_don_itr != aa_don_itr_end && aa_scaf_itr != aa_scaf_itr_end;
            ++align_itr, ++aa_scaf_itr, ++aa_don_itr, ++seq_length
          );
        }

        // sequence between the two secondary structure elements that are the cut a different points
        {
          Result new_sequence;
          new_sequence.m_Sequence.SetChainID( new_chain_id);

          // move align itr to next aligned scaffold aa, since the two anchor point alignments are separated by an empty assignment
          for( ; align_itr != align_itr_end && !( *align_itr)->GetMembers().IsDefined(); ++align_itr);

          // for replacement mode - not append
          current_donor_align_aa = *( ++( *align_itr)->GetMembers().Begin());

          // transform and add AAs from donor
          if( aa_don_itr != aa_don_itr_end)
          {
            sp_donor_used_begin = util::SiPtr< const biol::AABase>( **aa_don_itr);
          }

          for
          (
            ; aa_don_itr != aa_don_itr_end && ( *aa_don_itr)->GetData() != current_donor_align_aa->GetData();
            ++aa_don_itr, ++seq_length
          )
          {
            const biol::AABase &inserted_aa( ElongateSequence( new_sequence.m_Sequence, **aa_don_itr, transformation, seq_length));

            if( !sp_fused_begin.IsDefined())
            {
              sp_fused_begin = &inserted_aa;
            }

            sp_fused_last = &inserted_aa;
            sp_donor_used_last = util::SiPtr< const biol::AABase>( **aa_don_itr);

            inserted_donor_seq.push_back( inserted_aa.GetType()->GetOneLetterCode());
          }

          // selection onto replacement
          if( sp_donor_used_begin.IsDefined() && sp_donor_used_last.IsDefined())
          {
            new_sequence.m_PymolScript += std::string( "color donor_color, /") + donor_name + "//" + donor_chain_id + '/' +
                            util::Format()( sp_donor_used_begin->GetSeqID()) + '-' +
                            util::Format()( sp_donor_used_last->GetSeqID()) + '\n';
          }

          if( sp_fused_begin.IsDefined() && sp_fused_last.IsDefined())
          {
            new_sequence.m_PymolScript += "select " + donor_site_name + ", " + donor_site_name + " | /fused_model//"
              + new_chain_id + '/'
              + util::Format()( sp_fused_begin->GetSeqID()) + '-'
              + util::Format()( sp_fused_last->GetSeqID()) + '\n';
            new_sequence.m_PymolScript += "color donor_color, " + donor_site_name + '\n';
          }

          util::SiPtr< const biol::AABase> current_scaffold_align_aa( ( *align_itr)->GetMembers().FirstElement());
          if( !sp_scaffold_ignore_begin.IsDefined() && aa_scaf_itr != aa_scaf_itr_end)
          {
            sp_scaffold_ignore_begin = util::SiPtr< const biol::AABase>( **aa_scaf_itr);
            sp_scaffold_ignore_last  = util::SiPtr< const biol::AABase>( **aa_scaf_itr);
          }

          for
          (
            ; aa_scaf_itr != aa_scaf_itr_end && ( *aa_scaf_itr)->GetData() != current_scaffold_align_aa->GetData();
            ++aa_scaf_itr
          )
          {
            replaced_scaffold_seq.push_back( ( *aa_scaf_itr)->GetType()->GetOneLetterCode());
            sp_scaffold_ignore_last = util::SiPtr< const biol::AABase>( **aa_scaf_itr);
          }

          if( sp_scaffold_ignore_begin.IsDefined() && sp_scaffold_ignore_last.IsDefined())
          {
            new_sequence.m_PymolScript += std::string( "color ignore_color, /scaffold_model//") + scaf_chain_id + '/'
                + util::Format()( sp_scaffold_ignore_begin->GetSeqID()) + '-'
                + util::Format()( sp_scaffold_ignore_last->GetSeqID()) + '\n';
          }
          fused_sequence_segments.PushBack( storage::List< Result>( 1, new_sequence));
          BCL_MessageVrb( "sequence from donor " + util::Format()( count) + ':' + new_sequence.m_Sequence.GetSequenceIdentification());
        }

        {
          storage::List< Result> fused_sequences;
          for
          (
            align::AlignmentInterface< biol::AABase>::const_iterator align_itr_switch( align_itr);
              align_itr_switch != align_itr_end && ( *align_itr_switch)->GetMembers().IsDefined();
            ++align_itr_switch
          )
          {
            Result fused_sequence;
            fused_sequence.m_Sequence.SetChainID( new_chain_id);

            fused_sequence = HandleCutpointDonorScaffold
            (
              fused_sequence,
              align_itr, align_itr_switch, align_itr_end,
              aa_scaf_itr, aa_scaf_itr_end,
              aa_don_itr, aa_don_itr_end,
              donor_name, donor_site_name,
              transformation, count, seq_length
            );

            // create images
            PymolImages( fused_sequence.m_PymolScript, donor_name, donor_site_name);
            fused_sequences.PushBack( fused_sequence);
            BCL_MessageVrb( "sequence donor scaffold: " + util::Format()( count) + ':' + fused_sequence.m_Sequence.GetSequenceIdentification());
          }

          // add all fused sequences for possible cut points for this segment
          fused_sequence_segments.PushBack( fused_sequences);

          // move align, donor and scaffold itr to the end of the alignment
          for
          (
            ; align_itr != align_itr_end && ( *align_itr)->GetMembers().IsDefined() &&
              aa_don_itr != aa_don_itr_end && aa_scaf_itr != aa_scaf_itr_end;
            ++align_itr, ++aa_scaf_itr, ++aa_don_itr, ++seq_length
          );
        }
      }

      // add the rest from the scaffold
      if( aa_scaf_itr != aa_scaf_itr_end)
      {
        Result new_sequence;
        new_sequence.m_Sequence.SetChainID( new_chain_id);

        util::SiPtr< const biol::AABase> sp_last_fused;
        util::SiPtr< const biol::AABase> sp_last_scaf;

        std::string color_scaf( std::string( "color scaffold_color, /scaffold_model//") + scaf_chain_id + '/' + util::Format()( ( *aa_scaf_itr)->GetSeqID()) + '-');

        while( aa_scaf_itr != aa_scaf_itr_end)
        {
          //          BCL_MessageStd( "adding scaf end: " + ( *aa_scaf_itr)->GetIdentification());
          const biol::AABase &inserted_aa( ElongateSequence( new_sequence.m_Sequence, **aa_scaf_itr, seq_length));

          sp_last_fused = &inserted_aa;
          sp_last_scaf = util::SiPtr< const biol::AABase>( *aa_scaf_itr);
          ++aa_scaf_itr;
          ++seq_length;
        }

        if( sp_last_scaf.IsDefined())
        {
          new_sequence.m_ReplaceInformation += '-' + util::Format()( sp_last_scaf->GetSeqID());

          color_scaf += util::Format()( sp_last_scaf->GetSeqID()) + '\n';

          new_sequence.m_PymolScript += color_scaf;
        }

        new_sequence.m_PymolScript += "group sites, donor*_site\n";
        new_sequence.m_PymolScript += "disable sites\n";
        new_sequence.m_PymolScript += "group donor_models, donor*_model\n";
        new_sequence.m_PymolScript += "disable donor_models\n";

        new_sequence.m_PymolScript += "#system convert -pointsize 32 -delay 0 scaffold_model.png -delay 100 label:'scaffold' -delay 0 "
                        "donor1_model.png -delay 100 label:'donor' -delay 0 fused_model.png -delay 100 label:'fused' "
                        "-delay 0 fused_site.png -delay 100 label:'site' animation.gif\n";

        new_sequence.m_PymolScript += "enable all\n";
        new_sequence.m_PymolTailScript += "disable peptide_good\n";
        new_sequence.m_PymolTailScript += "disable peptide_bad\n";
        new_sequence.m_PymolTailScript += "disable peptide_atoms\n";

        fused_sequence_segments.PushBack( storage::List< Result>( 1, new_sequence));
      }

      return CombineSegments( fused_sequence_segments);
    }

    FusionProtein::Result FusionProtein::Replace
    (
      const biol::AASequence &SCAFFOLD,
      const storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, util::SiPtr< const assemble::ProteinModel> > > &ALIGNMENTS
    ) const
    {
      Result new_sequence;
      const char new_chain_id( SCAFFOLD.GetChainID());
      new_sequence.m_Sequence.SetChainID( new_chain_id);

      std::string replace_info;

      new_sequence.m_PymolScript += std::string( "show cartoon, /fused_model//") + new_chain_id + '\n';
      new_sequence.m_PymolScript += std::string( "show cartoon, /scaffold_model//") + new_chain_id + '\n';
      new_sequence.m_PymolScript += std::string( "color scaffold_color, /fused_model//") + new_chain_id + '\n';
      new_sequence.m_PymolScript += std::string( "color scaffold_color, /scaffold_model//") + new_chain_id + '\n';

      const char scaf_chain_id( SCAFFOLD.GetChainID());
      biol::AASequence::const_iterator aa_scaf_itr( SCAFFOLD.Begin()), aa_scaf_itr_end( SCAFFOLD.End());
      size_t count( 1);

      if( aa_scaf_itr != aa_scaf_itr_end)
      {
        new_sequence.m_ReplaceInformation += "scaf:" + util::Format()( ( *aa_scaf_itr)->GetSeqID());
      }

      // iterate through all alignments and donors
      for
      (
        storage::List< storage::Pair< align::AlignmentNode< biol::AABase>, util::SiPtr< const assemble::ProteinModel> > >::const_iterator
          itr( ALIGNMENTS.Begin()), itr_end( ALIGNMENTS.End());
        itr != itr_end;
        ++itr, ++count
      )
      {
        const std::string donor_name( "donor" + util::Format()( count) + "_model");
        const std::string donor_site_name( "donor" + util::Format()( count) + "_site");
        new_sequence.m_PymolScript += "select " + donor_site_name + ", " + donor_name + " & scaffold_model\n";
        const align::AlignmentNode< biol::AABase> &current_alignment( itr->First());

        const storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > coord_pair
        (
          assemble::Quality::CoordinatesFromAlignment( current_alignment, m_AtomTypes)
        );
        const math::TransformationMatrix3D transformation( m_SPMeasure->CalculateSuperimposition( coord_pair.Second(), coord_pair.First()));

        // residues at fusion point
        util::SiPtr< const biol::AABase> sp_scaffold_ignore_begin;
        util::SiPtr< const biol::AABase> sp_scaffold_ignore_last;

        util::SiPtr< const biol::AABase> sp_donor_used_begin;
        util::SiPtr< const biol::AABase> sp_donor_used_last;

        util::SiPtr< const biol::AABase> sp_fused_begin;
        util::SiPtr< const biol::AABase> sp_fused_last;

        std::string replaced_scaffold_seq;
        std::string inserted_donor_seq;

        if( current_alignment.IsEmpty() || current_alignment.GetDepth() != 2)
        {
          continue;
        }

        align::AlignmentInterface< biol::AABase>::const_iterator
          align_itr( current_alignment.GetAssignments().Begin()),
          align_itr_end( current_alignment.GetAssignments().End());

        // find first assignment
        for( ; align_itr != align_itr_end && !( *align_itr)->GetMembers().IsDefined(); ++align_itr);

        if( align_itr == align_itr_end)
        {
          continue;
        }

        util::SiPtr< const biol::AABase> current_scaffold_align_aa( ( *align_itr)->GetMembers().FirstElement());

        // add all amino acids from scaffold to fused sequence until start of the alignment
        for
        (
          ; aa_scaf_itr != aa_scaf_itr_end && ( *aa_scaf_itr)->GetData() != current_scaffold_align_aa->GetData();
          ++aa_scaf_itr
        )
        {
          ElongateSequence( new_sequence.m_Sequence, **aa_scaf_itr, 0);
        }

        const char donor_chain_id( ( *( ++current_alignment.GetSequences().Begin()))->GetFirstMember()->GetChainID());
        const biol::AASequence &donor_sequence( *itr->Second()->GetChain( donor_chain_id)->GetSequence());
        util::SiPtr< const biol::AABase> current_donor_align_aa( *( ++( *align_itr)->GetMembers().Begin()));

        biol::AASequence::const_iterator aa_don_itr( donor_sequence.Begin()), aa_don_itr_end( donor_sequence.End());
        new_sequence.m_PymolScript += std::string( "show cartoon, /")       + donor_name + "//" + donor_chain_id + '/' + '\n';
        new_sequence.m_PymolScript += std::string( "color ignore_color, /") + donor_name + "//" + donor_chain_id + '/' + '\n';

        // discard donor amino acids until the start of the alignment
        while( aa_don_itr != aa_don_itr_end && ( *aa_don_itr)->GetData() != current_donor_align_aa->GetData())
        {
          ++aa_don_itr;
        }

        {
          new_sequence = HandleCutpointScaffoldDonor
              (
                new_sequence,
                align_itr, FindBestPeptideBond( align_itr, align_itr_end, transformation, true), align_itr_end,
                aa_scaf_itr, aa_scaf_itr_end,
                aa_don_itr, aa_don_itr_end,
                donor_name, donor_site_name,
                transformation, count, 0
              );
        }

        // move align, donor and scaffold itr to the end of the alignment
        for
        (
          ; align_itr != align_itr_end && ( *align_itr)->GetMembers().IsDefined() &&
            aa_don_itr != aa_don_itr_end && aa_scaf_itr != aa_scaf_itr_end;
          ++align_itr, ++aa_scaf_itr, ++aa_don_itr
        );

        // move align itr to next aligned scaffold aa, since the two anchor point alignments are separated by an empty assignment
        for( ; align_itr != align_itr_end && !( *align_itr)->GetMembers().IsDefined(); ++align_itr);

        // for replacement mode - not append
        current_donor_align_aa = *( ++( *align_itr)->GetMembers().Begin());

        // transform and add AAs from donor
        if( aa_don_itr != aa_don_itr_end)
        {
          sp_donor_used_begin = util::SiPtr< const biol::AABase>( **aa_don_itr);
        }

        for
        (
          ; aa_don_itr != aa_don_itr_end && ( *aa_don_itr)->GetData() != current_donor_align_aa->GetData();
          ++aa_don_itr
        )
        {
          const biol::AABase &inserted_aa( ElongateSequence( new_sequence.m_Sequence, **aa_don_itr, transformation, 0));

          if( !sp_fused_begin.IsDefined())
          {
            sp_fused_begin = &inserted_aa;
          }

          sp_fused_last = &inserted_aa;
          sp_donor_used_last = util::SiPtr< const biol::AABase>( **aa_don_itr);

          inserted_donor_seq.push_back( inserted_aa.GetType()->GetOneLetterCode());
        }

        // selection onto replacement
        if( sp_donor_used_begin.IsDefined() && sp_donor_used_last.IsDefined())
        {
          new_sequence.m_PymolScript += std::string( "color donor_color, /") + donor_name + "//" + donor_chain_id + '/' +
                          util::Format()( sp_donor_used_begin->GetSeqID()) + '-' +
                          util::Format()( sp_donor_used_last->GetSeqID()) + '\n';
        }

        if( sp_fused_begin.IsDefined() && sp_fused_last.IsDefined())
        {
          new_sequence.m_PymolScript += "select " + donor_site_name + ", " + donor_site_name + " | /fused_model//"
            + new_chain_id + '/'
            + util::Format()( sp_fused_begin->GetSeqID()) + '-'
            + util::Format()( sp_fused_last->GetSeqID()) + '\n';
          new_sequence.m_PymolScript += "color donor_color, " + donor_site_name + '\n';
        }

        current_scaffold_align_aa = ( *align_itr)->GetMembers().FirstElement();
        if( !sp_scaffold_ignore_begin.IsDefined() && aa_scaf_itr != aa_scaf_itr_end)
        {
          sp_scaffold_ignore_begin = util::SiPtr< const biol::AABase>( **aa_scaf_itr);
          sp_scaffold_ignore_last  = util::SiPtr< const biol::AABase>( **aa_scaf_itr);
        }

        for
        (
          ; aa_scaf_itr != aa_scaf_itr_end && ( *aa_scaf_itr)->GetData() != current_scaffold_align_aa->GetData();
          ++aa_scaf_itr
        )
        {
          replaced_scaffold_seq.push_back( ( *aa_scaf_itr)->GetType()->GetOneLetterCode());
          sp_scaffold_ignore_last = util::SiPtr< const biol::AABase>( **aa_scaf_itr);
        }

        if( sp_scaffold_ignore_begin.IsDefined() && sp_scaffold_ignore_last.IsDefined())
        {
          new_sequence.m_PymolScript += std::string( "color ignore_color, /scaffold_model//") + scaf_chain_id + '/'
              + util::Format()( sp_scaffold_ignore_begin->GetSeqID()) + '-'
              + util::Format()( sp_scaffold_ignore_last->GetSeqID()) + '\n';
        }

        new_sequence = HandleCutpointDonorScaffold
          (
            new_sequence,
            align_itr, FindBestPeptideBond( align_itr, align_itr_end, transformation, false), align_itr_end,
            aa_scaf_itr, aa_scaf_itr_end,
            aa_don_itr, aa_don_itr_end,
            donor_name, donor_site_name,
            transformation, count, 0
          );

        // create images
        PymolImages( new_sequence.m_PymolScript, donor_name, donor_site_name);

        // move align, donor and scaffold itr to the end of the alignment
        for
        (
          ; align_itr != align_itr_end && ( *align_itr)->GetMembers().IsDefined() &&
            aa_don_itr != aa_don_itr_end && aa_scaf_itr != aa_scaf_itr_end;
          ++align_itr, ++aa_scaf_itr, ++aa_don_itr
        );

        replace_info += inserted_donor_seq;
        replace_info += '_';

        // create images
        PymolImages( new_sequence.m_PymolScript, donor_name, donor_site_name);

        BCL_MessageStd( "sequence replaced in scaffold: " + replaced_scaffold_seq);
        BCL_MessageStd( "sequence inserted from donor:  " + inserted_donor_seq);
        new_sequence.m_NrReplacedResidues += replaced_scaffold_seq.size();
      }

      // add the rest from the scaffold
      if( aa_scaf_itr != aa_scaf_itr_end)
      {
        util::SiPtr< const biol::AABase> sp_last_fused;
        util::SiPtr< const biol::AABase> sp_last_scaf;

        std::string color_fused( std::string( "color scaffold_color, /fused_model//")    + new_chain_id  + '/' + util::Format()(  new_sequence.m_Sequence.GetLastAA()->GetSeqID() + 1) + '-');
        std::string color_scaf( std::string( "color scaffold_color, /scaffold_model//") + scaf_chain_id + '/' + util::Format()( ( *aa_scaf_itr)->GetSeqID()) + '-');

        while( aa_scaf_itr != aa_scaf_itr_end)
        {
          //          BCL_MessageStd( "adding scaf end: " + ( *aa_scaf_itr)->GetIdentification());
          const biol::AABase &inserted_aa( ElongateSequence( new_sequence.m_Sequence, **aa_scaf_itr, 0));

          sp_last_fused = &inserted_aa;
          sp_last_scaf = util::SiPtr< const biol::AABase>( *aa_scaf_itr);
          ++aa_scaf_itr;
        }

        if( sp_last_scaf.IsDefined())
        {
          new_sequence.m_ReplaceInformation += '-' + util::Format()( sp_last_scaf->GetSeqID());

          color_fused += util::Format()( sp_last_fused->GetSeqID()) + '\n';
          color_scaf += util::Format()( sp_last_scaf->GetSeqID()) + '\n';

          new_sequence.m_PymolScript += color_fused;
          new_sequence.m_PymolScript += color_scaf;
        }
      }

      new_sequence.m_PymolScript += "group sites, donor*_site\n";
      new_sequence.m_PymolScript += "disable sites\n";
      new_sequence.m_PymolScript += "group donor_models, donor*_model\n";
      new_sequence.m_PymolScript += "disable donor_models\n";

      new_sequence.m_PymolScript += "#system convert -pointsize 32 -delay 0 scaffold_model.png -delay 100 label:'scaffold' -delay 0 "
                      "donor1_model.png -delay 100 label:'donor' -delay 0 fused_model.png -delay 100 label:'fused' "
                      "-delay 0 fused_site.png -delay 100 label:'site' animation.gif\n";

      new_sequence.m_PymolScript += "enable all\n";
      new_sequence.m_PymolTailScript += "disable peptide_good\n";
      new_sequence.m_PymolTailScript += "disable peptide_bad\n";
      new_sequence.m_PymolTailScript += "disable peptide_atoms\n";

      new_sequence.m_Sequence.SetFastaHeader( replace_info);

      return new_sequence;
    }

    storage::List< FusionProtein::Result>
    FusionProtein::ReplaceNandCTermSystematic
    (
      const biol::AASequence &SCAFFOLD,
      const align::AlignmentNode< biol::AABase> &NTERM_ALIGN,
      const assemble::ProteinModel &NTERM_DONOR,
      const align::AlignmentNode< biol::AABase> &CTERM_ALIGN,
      const assemble::ProteinModel &CTERM_DONOR
    ) const
    {
      const char scaf_chain_id( SCAFFOLD.GetChainID());
      const char new_chain_id( scaf_chain_id);

      std::string script_header;

      script_header += std::string( "show cartoon, /fused_model//") + new_chain_id + '\n';
      script_header += std::string( "show cartoon, /scaffold_model//") + new_chain_id + '\n';
      script_header += std::string( "color scaffold_color, /fused_model//") + new_chain_id + '\n';
      script_header += std::string( "color scaffold_color, /scaffold_model//") + new_chain_id + '\n';

      std::string label_distances_script;
      std::string replace_info;

      biol::AASequence new_sequence;
      new_sequence.SetChainID( new_chain_id);

      biol::AASequence::const_iterator aa_scaf_itr( SCAFFOLD.Begin()), aa_scaf_itr_end( SCAFFOLD.End());

      // collect all fused sequences
      storage::List< FusionProtein::Result> fused_sequences_nterm;

      // keep track of the total sequence length (which is the same for all constructs) so that identical subsequences
      // that are generate as intermediate steps start with the proper seq id
      size_t seq_length( 0);

      // nterm
      if( !NTERM_ALIGN.IsEmpty() && NTERM_ALIGN.GetDepth() == 2)
      {
        const char donor_chain_id( ( *( ++NTERM_ALIGN.GetSequences().Begin()))->GetFirstMember()->GetChainID());
        const biol::AASequence &donor_sequence( *NTERM_DONOR.GetChain( donor_chain_id)->GetSequence());
        biol::AASequence::const_iterator aa_don_itr( donor_sequence.Begin()), aa_don_itr_end( donor_sequence.End());

        const std::string donor_name( "donor_nterm_model");
        const std::string donor_site_name( "donor_nterm_site");

        script_header += std::string( "show cartoon, /")       + donor_name + "//" + donor_chain_id + '/' + '\n';
        script_header += std::string( "color ignore_color, /") + donor_name + "//" + donor_chain_id + '/' + '\n';
        script_header += "select " + donor_site_name + ", " + donor_name + " & scaffold_model\n";

        const storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > coord_pair
        (
          assemble::Quality::CoordinatesFromAlignment( NTERM_ALIGN, m_AtomTypes)
        );
        const math::TransformationMatrix3D transformation( m_SPMeasure->CalculateSuperimposition( coord_pair.Second(), coord_pair.First()));

        // residues at fusion point
        util::SiPtr< const biol::AABase> sp_left_scaffold;
        util::SiPtr< const biol::AABase> sp_right_scaffold;
        util::SiPtr< const biol::AABase> sp_right_donor;

        std::string replaced_scaffold_seq;
        std::string inserted_donor_seq;

        align::AlignmentInterface< biol::AABase>::const_iterator
          align_itr( NTERM_ALIGN.GetAssignments().Begin()),
          align_itr_end( NTERM_ALIGN.GetAssignments().End());

        // find first assignment
        for( ; align_itr != align_itr_end && !( *align_itr)->GetMembers().IsDefined(); ++align_itr);

        if( align_itr == align_itr_end)
        {
          return fused_sequences_nterm;
        }

        // donor amino acid is second within assignment
        util::SiPtr< const biol::AABase> current_donor_align_aa( *( ++( *align_itr)->GetMembers().Begin()));

        for
        (
          ; aa_don_itr != aa_don_itr_end && ( *aa_don_itr)->GetData() != current_donor_align_aa->GetData();
          ++aa_don_itr, ++seq_length
        )
        {
          ElongateSequence( new_sequence, **aa_don_itr, transformation, seq_length);
          sp_right_donor = util::SiPtr< const biol::AABase>( **aa_don_itr);
        }

        // selection onto replacement
        script_header += std::string( "color donor_color, /") + donor_name + "//" + donor_chain_id + '/' +
                         util::Format()( donor_sequence.GetFirstAA()->GetSeqID()) + '-' +
                         util::Format()( sp_right_donor->GetSeqID()) + '\n';

        // discard scaffold amino acids
        util::SiPtr< const biol::AABase> current_scaffold_align_aa( ( *align_itr)->GetMembers().FirstElement());
        if( !sp_left_scaffold.IsDefined() && aa_scaf_itr != aa_scaf_itr_end)
        {
          sp_left_scaffold = util::SiPtr< const biol::AABase>( **aa_scaf_itr);
          sp_right_scaffold = util::SiPtr< const biol::AABase>( **aa_scaf_itr);
        }

        for
        (
          ; aa_scaf_itr != aa_scaf_itr_end && ( *aa_scaf_itr)->GetData() != current_scaffold_align_aa->GetData();
          ++aa_scaf_itr
        )
        {
          replaced_scaffold_seq.push_back( ( *aa_scaf_itr)->GetType()->GetOneLetterCode());
          sp_right_scaffold = util::SiPtr< const biol::AABase>( **aa_scaf_itr);
        }
        script_header += std::string( "color ignore_color, /scaffold_model//") + scaf_chain_id + '/' +
                         util::Format()( sp_left_scaffold->GetSeqID()) + '-' +
                         util::Format()( sp_right_scaffold->GetSeqID()) + '\n';

        // color overlapping region before the replacement region
        util::SiPtr< const biol::AABase> sp_last_scaf;
        util::SiPtr< const biol::AABase> sp_last_donor;

        // for possible cut points, create a new sequence
        for
        (
          align::AlignmentInterface< biol::AABase>::const_iterator align_itr_switch( align_itr);
            align_itr_switch != align_itr_end && ( *align_itr_switch)->GetMembers().IsDefined();
          ++align_itr_switch
        )
        {
          Result fused_sequence;
          fused_sequence.m_Sequence = new_sequence;
          fused_sequence.m_PymolScript = script_header;

          fused_sequence = HandleCutpointDonorScaffold
              (
                fused_sequence,
                align_itr, align_itr_switch, align_itr_end,
                aa_scaf_itr, aa_scaf_itr_end,
                aa_don_itr, aa_don_itr_end,
                donor_name, donor_site_name,
                transformation, 0, seq_length
              );
          fused_sequences_nterm.PushBack( fused_sequence);
        }

        // move align, donor and scaffold itr to the end of the alignment
        for
        (
          ; align_itr != align_itr_end && ( *align_itr)->GetMembers().IsDefined() &&
            aa_don_itr != aa_don_itr_end && aa_scaf_itr != aa_scaf_itr_end;
          ++align_itr, ++aa_scaf_itr, ++aa_don_itr, ++seq_length
        );
      }
      else // did not replace N-term
      {
        Result fused_sequence;
        fused_sequence.m_Sequence = new_sequence;
        fused_sequence.m_PymolScript = script_header;
        fused_sequences_nterm.PushBack( fused_sequence);
      }

      // collect all fused sequences
      storage::List< FusionProtein::Result> fused_sequences;

      // cterm
      if( !CTERM_ALIGN.IsEmpty() && CTERM_ALIGN.GetDepth() == 2)
      {
        std::string pymol_script_cterm;

        const char donor_chain_id( ( *( ++CTERM_ALIGN.GetSequences().Begin()))->GetFirstMember()->GetChainID());
        const biol::AASequence &donor_sequence( *CTERM_DONOR.GetChain( donor_chain_id)->GetSequence());
        biol::AASequence::const_iterator aa_don_itr( donor_sequence.Begin()), aa_don_itr_end( donor_sequence.End());

        const std::string donor_name( "donor_cterm_model");
        const std::string donor_site_name( "donor_cterm_site");

        pymol_script_cterm += std::string( "show cartoon, /")       + donor_name + "//" + donor_chain_id + '/' + '\n';
        pymol_script_cterm += std::string( "color ignore_color, /") + donor_name + "//" + donor_chain_id + '/' + '\n';
        pymol_script_cterm += "select " + donor_site_name + ", " + donor_name + " & scaffold_model\n";

        const storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > coord_pair
        (
          assemble::Quality::CoordinatesFromAlignment( CTERM_ALIGN, m_AtomTypes)
        );
        const math::TransformationMatrix3D transformation( m_SPMeasure->CalculateSuperimposition( coord_pair.Second(), coord_pair.First()));

        // residues at fusion point
        util::SiPtr< const biol::AABase> sp_left_scaffold;
        util::SiPtr< const biol::AABase> sp_right_scaffold;
        util::SiPtr< const biol::AABase> sp_left_donor;

        std::string replaced_scaffold_seq;
        std::string inserted_donor_seq;

        align::AlignmentInterface< biol::AABase>::const_iterator
          align_itr( CTERM_ALIGN.GetAssignments().Begin()),
          align_itr_end( CTERM_ALIGN.GetAssignments().End());

        // find first assignment
        for( ; align_itr != align_itr_end && !( *align_itr)->GetMembers().IsDefined(); ++align_itr);

        if( align_itr == align_itr_end)
        {
          return fused_sequences_nterm;
        }

        biol::AASequence middle_sequence;
        middle_sequence.SetChainID( new_chain_id);

        util::SiPtr< const biol::AABase> current_scaffold_align_aa( ( *align_itr)->GetMembers().FirstElement());

        // add all scaffold amino acids until start of the alignment
        for
        (
          ; aa_scaf_itr != aa_scaf_itr_end && ( *aa_scaf_itr)->GetData() != current_scaffold_align_aa->GetData();
          ++aa_scaf_itr, ++seq_length
        )
        {
          ElongateSequence( middle_sequence, **aa_scaf_itr, seq_length);
        }

        // discard donor amino acids
        util::SiPtr< const biol::AABase> current_donor_align_aa( *( ++( *align_itr)->GetMembers().Begin()));
        while( aa_don_itr != aa_don_itr_end && ( *aa_don_itr)->GetData() != current_donor_align_aa->GetData())
        {
          ++aa_don_itr;
        }

        // iterate over all previously generated results
        for
        (
          storage::List< Result>::const_iterator result_itr( fused_sequences_nterm.Begin()), result_itr_end( fused_sequences_nterm.End());
          result_itr != result_itr_end;
          ++result_itr
        )
        {
          Result middle_fused_sequence( *result_itr);
          middle_fused_sequence.m_Sequence.AppendSequence( middle_sequence);
          middle_fused_sequence.m_PymolScript += pymol_script_cterm;

          // for possible cut points, create a new sequence
          for
          (
            align::AlignmentInterface< biol::AABase>::const_iterator align_itr_switch( align_itr);
            align_itr_switch != align_itr_end && ( *align_itr_switch)->GetMembers().IsDefined();
            ++align_itr_switch
          )
          {
            Result fused_sequence( middle_fused_sequence);

            fused_sequence = HandleCutpointScaffoldDonor
                (
                  fused_sequence,
                  align_itr, align_itr_switch, align_itr_end,
                  aa_scaf_itr, aa_scaf_itr_end,
                  aa_don_itr, aa_don_itr_end,
                  donor_name, donor_site_name,
                  transformation, 0, seq_length
                );

            fused_sequences.PushBack( fused_sequence);
          }
        }

        // move align, donor and scaffold itr to the end of the alignment
        for
        (
          ; align_itr != align_itr_end && ( *align_itr)->GetMembers().IsDefined() &&
            aa_don_itr != aa_don_itr_end && aa_scaf_itr != aa_scaf_itr_end;
          ++align_itr, ++aa_scaf_itr, ++aa_don_itr, ++seq_length
        );

        // add the rest of the donor sequence
        biol::AASequence term_sequence;
        term_sequence.SetChainID( new_chain_id);

        // add all scaffold amino acids until start of the alignment
        for( ; aa_don_itr != aa_don_itr_end; ++aa_don_itr, ++seq_length)
        {
          ElongateSequence( term_sequence, **aa_don_itr, seq_length);
        }

        // iterate over all previously generated results
        for
        (
          storage::List< Result>::iterator result_itr( fused_sequences.Begin()), result_itr_end( fused_sequences.End());
          result_itr != result_itr_end;
          ++result_itr
        )
        {
          result_itr->m_Sequence.AppendSequence( term_sequence);
        }
      }
      else // add the rest of the scaffold
      {
        // terminal sequence
        biol::AASequence term_sequence;
        term_sequence.SetChainID( new_chain_id);

        // add all scaffold amino acids until start of the alignment
        for( ; aa_scaf_itr != aa_scaf_itr_end; ++aa_scaf_itr, ++seq_length)
        {
          ElongateSequence( term_sequence, **aa_scaf_itr, seq_length);
        }

        // iterate over all previously generated results
        for
        (
          storage::List< Result>::iterator result_itr( fused_sequences.Begin()), result_itr_end( fused_sequences.End());
          result_itr != result_itr_end;
          ++result_itr
        )
        {
          result_itr->m_Sequence.AppendSequence( term_sequence);
        }
      }

      return fused_sequences;
    }

    FusionProtein::Result FusionProtein::ReplaceNandCTerm
    (
      const biol::AASequence &SCAFFOLD,
      const align::AlignmentNode< biol::AABase> &NTERM_ALIGN,
      const assemble::ProteinModel &NTERM_DONOR,
      const align::AlignmentNode< biol::AABase> &CTERM_ALIGN,
      const assemble::ProteinModel &CTERM_DONOR
    ) const
    {
      Result fused_sequence;

      const char new_chain_id( SCAFFOLD.GetChainID());
      fused_sequence.m_Sequence.SetChainID( new_chain_id);

      std::string replace_info;

      fused_sequence.m_PymolScript += std::string( "show cartoon, /fused_model//") + new_chain_id + '\n';
      fused_sequence.m_PymolScript += std::string( "show cartoon, /scaffold_model//") + new_chain_id + '\n';
      fused_sequence.m_PymolScript += std::string( "color scaffold_color, /fused_model//") + new_chain_id + '\n';
      fused_sequence.m_PymolScript += std::string( "color scaffold_color, /scaffold_model//") + new_chain_id + '\n';

      const char scaf_chain_id( SCAFFOLD.GetChainID());
      biol::AASequence::const_iterator aa_scaf_itr( SCAFFOLD.Begin()), aa_scaf_itr_end( SCAFFOLD.End());

      // nterm
      if( !NTERM_ALIGN.IsEmpty() && NTERM_ALIGN.GetDepth() == 2)
      {
        const char donor_chain_id( ( *( ++NTERM_ALIGN.GetSequences().Begin()))->GetFirstMember()->GetChainID());
        const biol::AASequence &donor_sequence( *NTERM_DONOR.GetChain( donor_chain_id)->GetSequence());
        biol::AASequence::const_iterator aa_don_itr( donor_sequence.Begin()), aa_don_itr_end( donor_sequence.End());

        const std::string donor_name( "donor_nterm_model");
        const std::string donor_site_name( "donor_nterm_site");

        fused_sequence.m_PymolScript += "select " + donor_site_name + ", " + donor_name + " & scaffold_model\n";

        fused_sequence.m_PymolScript += std::string( "show cartoon, /")       + donor_name + "//" + donor_chain_id + '/' + '\n';
        fused_sequence.m_PymolScript += std::string( "color ignore_color, /") + donor_name + "//" + donor_chain_id + '/' + '\n';

        const storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > coord_pair
        (
          assemble::Quality::CoordinatesFromAlignment( NTERM_ALIGN, m_AtomTypes)
        );
        const math::TransformationMatrix3D transformation( m_SPMeasure->CalculateSuperimposition( coord_pair.Second(), coord_pair.First()));

        // residues at fusion point
        util::SiPtr< const biol::AABase> sp_left_scaffold;
        util::SiPtr< const biol::AABase> sp_right_scaffold;
        util::SiPtr< const biol::AABase> sp_right_donor;

        std::string replaced_scaffold_seq;
        std::string inserted_donor_seq;

        align::AlignmentInterface< biol::AABase>::const_iterator
          align_itr( NTERM_ALIGN.GetAssignments().Begin()),
          align_itr_end( NTERM_ALIGN.GetAssignments().End());

        // find first assignment
        for( ; align_itr != align_itr_end && !( *align_itr)->GetMembers().IsDefined(); ++align_itr);

        if( align_itr == align_itr_end)
        {
          return fused_sequence;
        }

        util::SiPtr< const biol::AABase> current_donor_align_aa( *( ++( *align_itr)->GetMembers().Begin()));

        for
        (
          ; aa_don_itr != aa_don_itr_end && ( *aa_don_itr)->GetData() != current_donor_align_aa->GetData();
          ++aa_don_itr
        )
        {
          ElongateSequence( fused_sequence.m_Sequence, **aa_don_itr, transformation, 0);
          sp_right_donor = util::SiPtr< const biol::AABase>( **aa_don_itr);
        }

        // selection onto replacement
        fused_sequence.m_PymolScript += std::string( "color donor_color, /") + donor_name + "//" + donor_chain_id + '/' +
                        util::Format()( donor_sequence.GetFirstAA()->GetSeqID()) + '-' +
                        util::Format()( sp_right_donor->GetSeqID()) + '\n';

        // discard scaffold amino acids
        util::SiPtr< const biol::AABase> current_scaffold_align_aa( ( *align_itr)->GetMembers().FirstElement());
        if( !sp_left_scaffold.IsDefined() && aa_scaf_itr != aa_scaf_itr_end)
        {
          sp_left_scaffold = util::SiPtr< const biol::AABase>( **aa_scaf_itr);
          sp_right_scaffold = util::SiPtr< const biol::AABase>( **aa_scaf_itr);
        }

        for
        (
          ; aa_scaf_itr != aa_scaf_itr_end && ( *aa_scaf_itr)->GetData() != current_scaffold_align_aa->GetData();
          ++aa_scaf_itr
        )
        {
          replaced_scaffold_seq.push_back( ( *aa_scaf_itr)->GetType()->GetOneLetterCode());
          sp_right_scaffold = util::SiPtr< const biol::AABase>( **aa_scaf_itr);
        }

        fused_sequence.m_PymolScript += std::string( "color ignore_color, /scaffold_model//") + scaf_chain_id + '/' +
                        util::Format()( sp_left_scaffold->GetSeqID()) + '-' +
                        util::Format()( sp_right_scaffold->GetSeqID()) + '\n';

        fused_sequence = HandleCutpointDonorScaffold
          (
            fused_sequence,
            align_itr, FindBestPeptideBond( align_itr, align_itr_end, transformation, false), align_itr_end,
            aa_scaf_itr, aa_scaf_itr_end,
            aa_don_itr, aa_don_itr_end,
            donor_name, donor_site_name,
            transformation, 0, 0
          );

        // images for current site
        PymolImages( fused_sequence.m_PymolScript, donor_name, donor_site_name);

        // move align, donor and scaffold itr to the end of the alignment
        for
        (
          ; align_itr != align_itr_end && ( *align_itr)->GetMembers().IsDefined() &&
            aa_don_itr != aa_don_itr_end && aa_scaf_itr != aa_scaf_itr_end;
          ++align_itr, ++aa_scaf_itr, ++aa_don_itr
        );
      }

      // cterm
      if( !CTERM_ALIGN.IsEmpty() && CTERM_ALIGN.GetDepth() == 2)
      {
        const char donor_chain_id( ( *( ++CTERM_ALIGN.GetSequences().Begin()))->GetFirstMember()->GetChainID());
        const biol::AASequence &donor_sequence( *CTERM_DONOR.GetChain( donor_chain_id)->GetSequence());
        biol::AASequence::const_iterator aa_don_itr( donor_sequence.Begin()), aa_don_itr_end( donor_sequence.End());

        const std::string donor_name( "donor_cterm_model");
        const std::string donor_site_name( "donor_cterm_site");

        fused_sequence.m_PymolScript += "select " + donor_site_name + ", " + donor_name + " & scaffold_model\n";

        fused_sequence.m_PymolScript += std::string( "show cartoon, /")       + donor_name + "//" + donor_chain_id + '/' + '\n';
        fused_sequence.m_PymolScript += std::string( "color ignore_color, /") + donor_name + "//" + donor_chain_id + '/' + '\n';

        const storage::VectorND< 2, util::SiPtrVector< const linal::Vector3D> > coord_pair
        (
          assemble::Quality::CoordinatesFromAlignment( CTERM_ALIGN, m_AtomTypes)
        );
        const math::TransformationMatrix3D transformation( m_SPMeasure->CalculateSuperimposition( coord_pair.Second(), coord_pair.First()));

        // residues at fusion point
        util::SiPtr< const biol::AABase> sp_left_scaffold;
        util::SiPtr< const biol::AABase> sp_right_scaffold;
        util::SiPtr< const biol::AABase> sp_left_donor;
        util::SiPtr< const biol::AABase> sp_right_donor;
        util::SiPtr< const biol::AABase> sp_fused_before;
        util::SiPtr< const biol::AABase> sp_fused_after;

        std::string replaced_scaffold_seq;
        std::string inserted_donor_seq;

        align::AlignmentInterface< biol::AABase>::const_iterator
          align_itr( CTERM_ALIGN.GetAssignments().Begin()),
          align_itr_end( CTERM_ALIGN.GetAssignments().End());

        // find first assignment
        for( ; align_itr != align_itr_end && !( *align_itr)->GetMembers().IsDefined(); ++align_itr);

        if( align_itr == align_itr_end)
        {
          return fused_sequence;
        }

        // add all scaffold amino acids until start of the alignment
        for
        (
          util::SiPtr< const biol::AABase> current_scaffold_align_aa( ( *align_itr)->GetMembers().FirstElement());
          aa_scaf_itr != aa_scaf_itr_end && ( *aa_scaf_itr)->GetData() != current_scaffold_align_aa->GetData();
          ++aa_scaf_itr
        )
        {
          ElongateSequence( fused_sequence.m_Sequence, **aa_scaf_itr, 0);
        }

        // discard donor amino acids
        util::SiPtr< const biol::AABase> current_donor_align_aa( *( ++( *align_itr)->GetMembers().Begin()));
        while( aa_don_itr != aa_don_itr_end && ( *aa_don_itr)->GetData() != current_donor_align_aa->GetData())
        {
          ++aa_don_itr;
        }

        // color overlapping region before the replacement region
        {
          util::SiPtr< const biol::AABase> sp_last_scaf;
          util::SiPtr< const biol::AABase> sp_last_donor;

          fused_sequence = HandleCutpointScaffoldDonor
              (
                fused_sequence,
                align_itr, FindBestPeptideBond( align_itr, align_itr_end, transformation, true), align_itr_end,
                aa_scaf_itr, aa_scaf_itr_end,
                aa_don_itr, aa_don_itr_end,
                donor_name, donor_site_name,
                transformation, 0, 0
              );

          // images for current site
          PymolImages( fused_sequence.m_PymolScript, donor_name, donor_site_name);

          // move align, donor and scaffold itr to the end of the alignment
          for
          (
            ; align_itr != align_itr_end && ( *align_itr)->GetMembers().IsDefined() &&
              aa_don_itr != aa_don_itr_end && aa_scaf_itr != aa_scaf_itr_end;
            ++align_itr, ++aa_scaf_itr, ++aa_don_itr
          );

          // add all remaining residues from the donor
          for( ; aa_don_itr != aa_don_itr_end; ++aa_don_itr)
          {
            const biol::AABase &inserted_aa( ElongateSequence( fused_sequence.m_Sequence, **aa_don_itr, transformation, 0));

            // update the first aa of the alignment for the pymol focus
            if( !sp_fused_after.IsDefined())
            {
              sp_fused_after = &inserted_aa;
            }

            inserted_donor_seq.push_back( inserted_aa.GetType()->GetOneLetterCode());
          }
        }
      }
      else // add the rest of the scaffold
      {
        // add all scaffold amino acids until start of the alignment
        for( ; aa_scaf_itr != aa_scaf_itr_end; ++aa_scaf_itr)
        {
          ElongateSequence( fused_sequence.m_Sequence, **aa_scaf_itr, 0);
        }
      }

      fused_sequence.m_Sequence.SetFastaHeader( replace_info);

      return fused_sequence;
    }

    //! @brief generate a pymol script that highlights the fusion point by desiplaying the petide bond and grouping them
    //!        into good or bad peptide bonds
    //! @param AA_LEFT the left aa of the peptide bond (CA + C)
    //! @param AA_RIGHT the right aa of the peptide bond (N + CA)
    //! @param SIDE define a side uniquely as 'r', 'l' or 'n' right of site, left of site or fusion with single anchor
    //! @param COUNT numbering along the sequence
    //! @return string that can be appended to the pymol script to highlight funsion bonds
    std::string FusionProtein::PymolClosure
    (
      const biol::AABase &AA_LEFT,
      const biol::AABase &AA_RIGHT,
      const char SIDE,
      const size_t COUNT
    )
    {
      const bool proper_bond( biol::AABase::AreAminoAcidsPeptideBonded( AA_LEFT, AA_RIGHT, false));
      const char chain_id( AA_LEFT.GetChainID());
      std::string script;

      const std::string identifier( util::Format()( COUNT) + SIDE);
      const std::string select_cal( std::string( "/fused_model//") + chain_id + '/' + util::Format()( AA_LEFT.GetSeqID())  + "/CA");
      const std::string select_c(   std::string( "/fused_model//") + chain_id + '/' + util::Format()( AA_LEFT.GetSeqID())  + "/C");
      const std::string select_n(   std::string( "/fused_model//") + chain_id + '/' + util::Format()( AA_RIGHT.GetSeqID()) + "/N");
      const std::string select_car( std::string( "/fused_model//") + chain_id + '/' + util::Format()( AA_RIGHT.GetSeqID()) + "/CA");

      const std::string name_cal( "peptide_atom_cal" + identifier);
      const std::string name_c(   "peptide_atom_c"   + identifier);
      const std::string name_n(   "peptide_atom_n"   + identifier);
      const std::string name_car( "peptide_atom_car" + identifier);

      script += "pseudoatom " + name_cal + ", selection=" + select_cal + '\n';
      script += "pseudoatom " + name_c   + ", selection=" + select_c + '\n';
      script += "pseudoatom " + name_n   + ", selection=" + select_n + '\n';
      script += "pseudoatom " + name_car + ", selection=" + select_car + '\n';
      script += "group atoms" + identifier + ", " + name_cal + ' ' + name_c + ' ' + name_n + ' ' + name_car + '\n';
      script += "group peptide_atoms, atoms" + identifier + '\n';

      script += "show sticks, " + select_cal + ' ' + select_c + '\n';
      script += "show sticks, " + select_n + ' ' + select_car + '\n';
      script += "color carbon, "   + select_c + '\n';
      script += "color nitrogen, " + select_n + '\n';

      script += std::string( "dihedral angle") + identifier + ", " + name_cal + ", " + name_c + ", " + name_n + ", " + name_car + '\n';
      script += std::string( "distance dist") + identifier + ", " + name_c + ", " + name_n + '\n';

      if( proper_bond)
      {
        script += "group peptide_good, angle" + identifier + '\n';
        script += "group peptide_good, dist" + identifier + '\n';
      }
      else
      {
        script += "group peptide_bad, angle" + identifier + '\n';
        script += "group peptide_bad, dist" + identifier + '\n';
      }

      return script;
    }

    std::string &FusionProtein::PymolImages
    (
      std::string &PYMOL_SCRIPT,
      const std::string &DONOR_NAME,
      const std::string &DONOR_SITE_NAME
    )
    {
      PYMOL_SCRIPT += "disable all\n";

      // fused model
      PYMOL_SCRIPT += "enable fused_model\n";

      // zoom onto replacement
      PYMOL_SCRIPT += std::string( "orient ") + DONOR_SITE_NAME + '\n';
      PYMOL_SCRIPT += "rotate y, 90\n";
      PYMOL_SCRIPT += "#ray\n";
      PYMOL_SCRIPT += "png fused_site_" + DONOR_NAME + ".png, dpi=300\n";

      PYMOL_SCRIPT += "orient fused_model\n";
      PYMOL_SCRIPT += "#ray\n";
      PYMOL_SCRIPT += "png fused_model_" + DONOR_NAME + ".png, dpi=300\n";

      PYMOL_SCRIPT += "disable all\n";
      PYMOL_SCRIPT += "enable " + DONOR_NAME + '\n';
      PYMOL_SCRIPT += "#ray\n";
      PYMOL_SCRIPT += "png " + DONOR_NAME + ".png, dpi=300\n";

      PYMOL_SCRIPT += "disable all\n";
      PYMOL_SCRIPT += "enable scaffold_model\n";
      PYMOL_SCRIPT += "#ray\n";
      PYMOL_SCRIPT += "png scaffold_model_" + DONOR_NAME + ".png, dpi=300\n";

      // end
      return PYMOL_SCRIPT;
    }

    //! @brief for an alignment find the two neighboring assignments where the
    //!        TOP=true:  first amino acid in the assignment has a most optimal connection to the next assignments second amino acid
    //!        TOP=false: second amino acid in the assignment has a most optimal connection to the next assignments first amino acid
    //! @param ITR_FIRST first position within the alignment
    //! @param ITR_END end position within the alignment
    //! @param TRANSFORMATION_BOTTOM transformation to apply to the second amino acid's coordinates before calcuating
    //!        the peptide bond parameters
    //! @param TOP calcuate from top first to following second assignments amino acid, or from second to following
    //!        assignments first amino acid
    //! @return iterator on the assignment for the following amino acid in the peptide bond
    util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator FusionProtein::FindBestPeptideBond
    (
      const align::AlignmentInterface< biol::AABase>::const_iterator &ITR_FIRST,
      const align::AlignmentInterface< biol::AABase>::const_iterator &ITR_END,
      const math::TransformationMatrix3D &TRANSFORMATION_BOTTOM,
      const bool TOP
    )
    {
      static const double s_peptide_bond_length( 0.9 * biol::GetAtomTypes().C->GetBondLength( biol::GetAtomTypes().N));

      util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator itr1( ITR_FIRST);
      if( itr1 == ITR_END)
      {
        return ITR_END;
      }

      util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator itr2( ITR_FIRST);
      ++itr2;
      if( itr2 == ITR_END)
      {
        return ITR_END;
      }

      if( ( *itr1)->GetSize() != 2)
      {
        return ITR_END;
      }

      util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator best_itr( ITR_END);

      storage::VectorND< 2, double> best_bond_length_angle( 3.0, 0.0);

      for( ; itr2 != ITR_END && ( *itr1)->GetMembers().IsDefined() && ( *itr2)->GetMembers().IsDefined(); ++itr1, ++itr2)
      {
        util::ShPtr< biol::AABase> sp_transformed_aa;
        // are the residues in alignment matching the scaffold and donor residues
        util::SiPtr< const biol::AABase> sp_left_aa, sp_right_aa;
        if( TOP)
        {
          sp_left_aa = ( *itr1)->GetMembers().FirstElement();
          sp_transformed_aa = util::ShPtr< biol::AABase>( ( *( ++( *itr2)->GetMembers().Begin()))->Clone());
          sp_transformed_aa->Transform( TRANSFORMATION_BOTTOM);
          sp_right_aa = sp_transformed_aa;
        }
        else
        {
          sp_transformed_aa = util::ShPtr< biol::AABase>( ( *( ++( *itr1)->GetMembers().Begin()))->Clone());
          sp_transformed_aa->Transform( TRANSFORMATION_BOTTOM);
          sp_left_aa = sp_transformed_aa;
          sp_right_aa = ( *itr2)->GetMembers().FirstElement();
        }

        const storage::VectorND< 2, double> current_bond_length_angle( biol::AABase::PeptideBondLengthAndAngle( *sp_left_aa, *sp_right_aa));

        if( current_bond_length_angle.First() > s_peptide_bond_length && current_bond_length_angle.First() < best_bond_length_angle.First())
        {
          best_bond_length_angle = current_bond_length_angle;
          best_itr = itr2;
        }
      }
      BCL_MessageStd( "best bond length and angle: " + util::Format()( best_bond_length_angle.First()));
      return best_itr;
    }

    //! @brief elongate a sequence by a given amino acid applying a transformation to it
    //! @param SEQUENCE the sequence to elongate
    //! @param AMINO_ACID the amino acid that is added to the sequence
    //! @param TRANSFORMATION the transformation to apply to the amino acid that is added
    //! @param SEQ_ID seq_id to use if given sequence is empty
    //! @return reference to the added amino acid
    biol::AABase &FusionProtein::ElongateSequence
    (
      biol::AASequence &SEQUENCE, const biol::AABase &AMINO_ACID,
      const math::TransformationMatrix3D &TRANSFORMATION, const int SEQ_ID
    )
    {
      int current_seq_id( SEQ_ID);
      if( SEQUENCE.GetSize() > 0)
      {
        current_seq_id = SEQUENCE.GetLastAA()->GetSeqID();
      }
      util::ShPtr< biol::AABase> copy( AMINO_ACID.Clone());
      util::ShPtr< biol::AAData> new_aadata
      (
        new biol::AAData( AMINO_ACID.GetType(), ++current_seq_id, AMINO_ACID.GetPdbID(), ' ', SEQUENCE.GetChainID())
      );
      copy->SetData( new_aadata);
      copy->Transform( TRANSFORMATION);
      SEQUENCE.PushBack( copy);

      // end
      return *copy;
    }

    //! @brief elongate a sequence by a given amino acid
    //! @param SEQUENCE the sequence to elongate
    //! @param AMINO_ACID the amino acid that is added to the sequence
    //! @param SEQ_ID seq_id to use if given sequence is empty
    //! @return reference to the added amino acid
    biol::AABase &FusionProtein::ElongateSequence( biol::AASequence &SEQUENCE, const biol::AABase &AMINO_ACID, const int SEQ_ID)
    {
      int current_seq_id( SEQ_ID);
      if( SEQUENCE.GetSize() > 0)
      {
        current_seq_id = SEQUENCE.GetLastAA()->GetSeqID();
      }
      util::ShPtr< biol::AABase> copy( AMINO_ACID.Clone());
      util::ShPtr< biol::AAData> new_aadata
      (
        new biol::AAData( AMINO_ACID.GetType(), ++current_seq_id, AMINO_ACID.GetPdbID(), ' ', SEQUENCE.GetChainID())
      );
      copy->SetData( new_aadata);
      SEQUENCE.PushBack( copy);

      // end
      return *copy;
    }

    //! @brief create fused model from fused sequence adding all remaining scaffold chains
    //! @param FUSED_SEQUENCE the fused sequence
    //! @param SCAFFOLD scaffold protein
    //! @return protein with fused chain and remaining chains
    assemble::ProteinModel FusionProtein::CreateFusedModel
    (
      const biol::AASequence &FUSED_SEQUENCE,
      const assemble::ProteinModel &SCAFFOLD
    )
    {
      const char scaffold_chain_id( FUSED_SEQUENCE.GetChainID());
      util::ShPtr< assemble::Chain> fused_chain( new assemble::Chain( util::CloneToShPtr( FUSED_SEQUENCE)));
      fused_chain->Insert( util::ShPtr< assemble::SSE>( new assemble::SSE( FUSED_SEQUENCE, biol::GetSSTypes().COIL)));
      assemble::ProteinModel fused_model( fused_chain);

      const std::string scaffold_chain_ids( SCAFFOLD.GetChainIDs());

      // add unused chains from the scaffold
      for( assemble::ProteinModel::const_iterator itr( SCAFFOLD.GetChains().Begin()), itr_end( SCAFFOLD.GetChains().End()); itr != itr_end; ++itr)
      {
        if( ( *itr)->GetChainID() == scaffold_chain_id)
        {
          continue;
        }
        fused_model.Insert( *itr);
      }

      fold::ProtocolLoopClose::SplitCoilsAtNonPetideBond( fused_model);
      // assign secondary structure
      {
        biol::DSSP dssp;
        math::MutateResult< assemble::ProteinModel> result( dssp( fused_model));
        if( result.GetArgument().IsDefined())
        {
          fused_model = *result.GetArgument();
        }
      }

      // end
      return fused_model;
    }

    //! @brief assemble multiple segments into a list of all combinations of sequence segments
    //! @param SEGMENTS a list of list of segments with varying cut points
    //! @return a list of complete sequences containing all combinations of segments
    storage::List< FusionProtein::Result>
    FusionProtein::CombineSegments( const storage::List< storage::List< FusionProtein::Result> > &SEGMENTS)
    {
      size_t result_count( 1);
      for
      (
        storage::List< storage::List< Result> >::const_iterator itr( SEGMENTS.Begin()), itr_end( SEGMENTS.End());
        itr != itr_end;
        ++itr
      )
      {
        result_count *= itr->GetSize();
      }

      BCL_MessageVrb
      (
        "combining " + util::Format()( SEGMENTS.GetSize()) + " different segments along the sequence into " +
        util::Format()( result_count) + " results!"
      );

      storage::List< Result> combined_sequences;
      if( SEGMENTS.IsEmpty())
      {
        return combined_sequences;
      }

      storage::List< storage::List< Result> >::const_iterator
        itr( SEGMENTS.Begin()), itr_end( SEGMENTS.End());

      combined_sequences = *itr;

      while( ++itr != itr_end)
      {
        if( itr->IsEmpty())
        {
          continue;
        }
        storage::List< Result> elongated_sequences;
        for
        (
          storage::List< Result>::const_iterator itr1( combined_sequences.Begin()), itr1_end( combined_sequences.End());
          itr1 != itr1_end;
          ++itr1
        )
        {
          std::string previous_sequence;

          for
          (
            storage::List< Result>::const_iterator itr2( itr->Begin()), itr2_end( itr->End());
            itr2 != itr2_end;
            previous_sequence = itr2->m_Sequence.Sequence(), ++itr2
          )
          {
            std::string current_sequence( itr2->m_Sequence.Sequence());
//            if( current_sequence == previous_sequence)
//            {
//              continue;
//            }
            Result new_sequence( *itr1);
            new_sequence.Append( *itr2);
            BCL_MessageVrb( "new sequ comb: " + new_sequence.m_ReplaceInformation);
            elongated_sequences.PushBack( new_sequence);
          }
        }
        combined_sequences.InternalData().swap( elongated_sequences.InternalData());
      }

      BCL_MessageVrb( "combined into " + util::Format()( combined_sequences.GetSize()) + " sequences!");

      return combined_sequences;
    }

    //! @brief add a line for the rosetta resfile for the pair of amino acids for the given sequence position in the fused protein
    //! @param AA_FUSED the amino acid in the fused sequence
    //! @param AA_SCAFFOLD the scaffold amino acid
    //! @param AA_DONOR the donor amino acid
    //! @return the string containg the information required to define the specifics for that residue in the rosetta res file
    std::string FusionProtein::RosettaResFileEntry
    (
      const biol::AAData &AA_FUSED,
      const biol::AAData &AA_SCAFFOLD,
      const biol::AAData &AA_DONOR
    )
    {
      std::string entry;

      // for identical AAs comment the line, but still write it
      if( AA_SCAFFOLD.GetType() == AA_DONOR.GetType())
      {
        entry = "# no need to design identical AA\n#";
      }
      // position in the fused protein with the possible amino acid identies at this position
      entry += util::Format()( AA_FUSED.GetSeqID()) + ' ' + AA_FUSED.GetChainID() +
               " PIKAA " + AA_DONOR.GetType()->GetOneLetterCode() + AA_SCAFFOLD.GetType()->GetOneLetterCode() + '\n';

      return entry;
    }

    //! @brief create a histogram that can hold scores for two cross overs (scaf->don and don->scaf)
    //! @param FUSED_SEQUENCE the sequences representing the cross over results
    //! @return Gnuplotheatmap that covers the seqids for the crossovers
    math::Histogram2D FusionProtein::CreateCrossoverHistogram() const
    {
      math::RunningMinMax< double> seqid_range_scaf;
      math::RunningMinMax< double> seqid_range_don;

      // iterate through the fused sequences
      for
      (
        storage::Table< double>::const_iterator itr( m_Results.Begin()), itr_end( m_Results.End());
        itr != itr_end;
        ++itr
      )
      {
        const int last_scaffold_seqid( itr->Second()[ "crossover_scaffold"]);
        const int last_donor_seqid( itr->Second()[ "crossover_donor"]);

        if( last_scaffold_seqid < 1 || last_donor_seqid < 1)
        {
          continue;
        }

        seqid_range_scaf += double( last_scaffold_seqid);
        seqid_range_don += double( last_donor_seqid);
      }

      math::Histogram2D heatmap
      (
        storage::VectorND< 2, double>( seqid_range_scaf.GetMin() - 1.5, seqid_range_don.GetMin() - 1.5),
        storage::VectorND< 2, double>( 1, 1),
        storage::VectorND< 2, size_t>
        (
          seqid_range_scaf.GetMax() - seqid_range_scaf.GetMin() + 3,
          seqid_range_don.GetMax() - seqid_range_don.GetMin() + 3
        )
      );

      return heatmap;
    }

    const ApplicationType FusionProtein::FusionProtein_Instance
    (
      GetAppGroups().AddAppToGroup( new FusionProtein(), GetAppGroups().e_Protein)
    );

  } // namespace app
} // namespace bcl
