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
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "fold/bcl_fold_default_flags.h"
#include "fold/bcl_fold_stage_factory.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Fold
    //! @brief Protein Structure Prediction program
    //! @detail Main application for protein structure prediction. The Monte Carlo Minimizer is the main feature.
    //! Multiple chains are input as fasta files over the command line. The pool of SSEs is given in a file. Moves and
    //! scoring function are defaulted. Scoring weights can be given in a weights file, which is a bcl table. Scoring
    //! and Move objects can also be given as bcl objects in a file.
    //! Command line options are available to fold membrane proteins or to use experimentally derived restraints like
    //! NMR, EPR and Cryo-EM.
    //!
    //! @author karakam, woetzen, fischea, pinojc
    //! @date Feb 25, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Fold :
      public InterfaceRelease
    {

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    private:

      //! @brief default constructor
      Fold();

    public:

      //! @brief Clone function
      //! @return pointer to new Fold
      Fold *Clone() const
      {
        return new Fold( *this);
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

      //! @brief get a description for the app
      //! @return a brief (no more than 3 line) description for the application
      std::string GetDescription() const
      {
        return "de-novo protein structure prediction algorithm that assembles a pool of secondary structure elements. "
               "Allows for usage of experimental restraints";
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      //! @return shared pointer to the command object
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize new command
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // pushback all flags and parameters
        sp_cmd->PushBack( fold::Setup::GetAllFlags());

        // add default bcl parameters
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        return sp_cmd;
      }

      //! @brief conducts the de novo protein structure prediction
      //! @detail assembles the tertiary protein structure from the provided SSE pool according to the given stage file
      //! or command line arguments
      //! @return error code - 0 for success
      int Main() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief returns web text information
      //! @return text (html allowed but not required) that will be displayed on the website
      //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
      const std::string &GetWebText() const;

      // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
      // but which should not be displayed, e.g. for help
      storage::Vector< std::string> GetDeprecatedAppNames() const
      {
        return storage::Vector< std::string>( size_t( 1), "EMFold");
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
        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }

    private:

      //! static instance of this class
      static const ApplicationType Fold_Instance;

    }; // class Fold

      //! @brief conducts the de novo protein structure prediction
      //! @detail assembles the tertiary protein structure from the provided SSE pool according to the given stage file
      //! or command line arguments
      //! @return error code - 0 for success
    int Fold::Main() const
    {
      // initialize the fold setup which contains the information regarding the chosen options for the prediction process
      fold::Setup::InitializeStaticInstance();
      const fold::Setup &setup( fold::GetSetup());
      const fold::DefaultFlags flags;

      // output information about the SSE pool used
      BCL_MessageCrt( "The pool used is");

      util::ShPtr< assemble::SSEPool> sp_pool
      (
        setup.GetEmptyModel()->GetProteinModelData()->GetData( assemble::ProteinModelData::e_Pool)
      );

      BCL_Assert( sp_pool.IsDefined(), "the pool is not initialized!");
      sp_pool->WriteSSEPool( util::GetLogger());

      // construct the stages from the given stage file or command line options
      util::ShPtrVector< fold::StageInterface> stages( fold::StageFactory::CreateStages());

      // initialize counters that keep track of model and stage number
      const size_t total_model_number( setup.GetNumberRounds());
      const size_t total_stage_number( stages.GetSize());

      // shared pointer to the result of the previous stage
      util::ShPtr< assemble::ProteinModel> sp_previous_result;

      // loop to create the specified number of models
      for( size_t current_model_number( 0); current_model_number != total_model_number; ++current_model_number)
      {
        BCL_MessageCrt
        (
          "Folding model #" + util::Format()( current_model_number + 1) + "/" + util::Format()( total_model_number)
        );

        // iterate over the specified stages to generate the model
        for( size_t current_stage_number( 0); current_stage_number != total_stage_number; ++current_stage_number)
        {
          BCL_MessageCrt
          (
            "stage #" + util::Format()( current_stage_number + 1) + "/" + util::Format()( total_stage_number) +
            " of round #" + util::Format()( current_model_number + 1) + "/" + util::Format()( total_model_number)
          );

          // initialize the starting model
          BCL_MessageStd( "initializing starting model");
          bool is_first_stage( current_stage_number == 0);
          util::ShPtr< assemble::ProteinModel> sp_starting_model
          (
            is_first_stage ? // in the first stage use a specified starting model if provided
              (
                fold::DefaultFlags::GetFlagStartModel()->GetFlag() ?
                util::ShPtr< assemble::ProteinModel>( setup.GetStartModel()->HardCopy()) :
                util::ShPtr< assemble::ProteinModel>( setup.GetEmptyModel()->HardCopy())
              ) : // in later stages use the results from previous stages as starting model
              (
                util::ShPtr< assemble::ProteinModel>( sp_previous_result)
              )
          );
          BCL_MessageStd( "done initializing starting model");

          // start the approximation with the algorithm encapsulated in the stage
          util::ShPtr< fold::StageInterface> sp_stage( stages( current_stage_number));
          const bool is_last_stage( current_stage_number == total_stage_number - 1);
          sp_stage->Initialize( current_model_number, current_stage_number, is_last_stage);
          util::ShPtr< assemble::ProteinModel> sp_approximation_result( sp_stage->Approximate( sp_starting_model));

          // update the previous final model
          sp_previous_result = sp_approximation_result;
        }
      }

      // end
      return 0;
     }

    //! @brief returns web text information
    //! @return text (html allowed but not required) that will be displayed on the website
    //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
    const std::string &Fold::GetWebText() const
    {
      // create a static string to hold readme information
      static const std::string s_web_text
      (
        "BCL::Fold: De novo prediction of protein structures\n\n"
        "BCL::Fold uses a Monte Carlo sampling approach with simulated annealing to place secondary structure elements "
        "(SSEs) in the three-dimension space. After each Monte Carlo step the structures are scored with knowledge-based "
        "scoring functions to determine the accuracy of the model. Those scoring functions include among others the "
        "radius of gyration, SSE pairing and packing, amino acid environment and amino acid pair distances. The usage "
        "of precomputed SSEs favors the sampling of non-local contacts and thereby allows the prediction of complex "
        "protein topologies.\n"
        "!fold_workflow.jpeg!\n"
        "Fig. 1:\n"
        "\n\nIn a first step a pool of SSEs in predicted from  the primary structure of the protein (A). Random moves "
        " are applied to the model which include adding SSEs from the pool and transformations of the SSEs in the model "
        "(B). After each move the model is scored with the knowledge-based scoring functions in order to determine the "
        " accuracy of the prediction. Depending on the score and a Metropolis criterion the move is either in accepted "
        " or rejected (C). This process is repeated until the energy minimum is reached (D).\n\n"
        "BCL::Fold was benchmarked on soluble proteins on 66 proteins with lengths between 83 and 293 amino acids. For "
        "61 out of these proteins the best SSE-only models obtained have an RMSD100 below 8.0 Å and recover more than "
        "20% of the native contacts.\n\n"
        "BCL:Fold is able to incorporate experimentally determined data from EPR, NMR and Cryo-EM experiments into the "
        "folding process. Adding protein-specific information significantly improves the accuracy of the sampling and "
        "model selection. This was tested on a benchmark set of 25 monomeric and 5 multimeric membrane proteins with "
        "sequence lengths ranging from 69 to 568 residues which were folded with and without EPR distance and accessibility "
        "restraints. For 28 out of 30 cases BCL::Fold was able to sample structures with an RMSD100 below 8.0 Å when "
        "compared to the native structure."
      );
      return s_web_text;
    }

    //! @brief default constructor
    Fold::Fold()
    {
    }

    //! static instance of this class
    const ApplicationType Fold::Fold_Instance( GetAppGroups().AddAppToGroup( new Fold(), GetAppGroups().e_Protein));

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &Fold::GetReadMe() const
    {
      // create a static string to hold readme information
      static const std::string s_readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::Fold, terms of use, "
        "appropriate citation, installation procedures, BCL::Fold execution, "
        "technical support, and future research directions.\n\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::Fold?\n"
        "BCL::Fold is a C++ based application, created by Vanderbilt University's "
        "Meiler Laboratory, which is part of a larger library of applications "
        "called BCL::Commons.  BCL::Fold is a de novo protein structure prediction "
        "algorithm. A pool of secondary structure elements is assembled using a Monte Carlo "
        "assembly process and evaluated using knowledge-based energy functions.\n\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::Fold.\n"
        "When using BCL::Fold in a publication, please cite the publication describing the application's development:\n"
        "Karakas, M., Woetzel, N., Staritzbichler, R., Alexander, N., Weiner B.E., and Meiler, J. (2012) "
        "BCL::Fold--de novo prediction of complex and large protein topologies by assembly of secondary structure "
        "elements. PLoS One, 7(11).\n"
        "N. Woetzel, M. Karakaş, R. Staritzbichler, R. Müller, B. E. Weiner, and J. Meiler, BCL::Score--knowledge "
        "based energy potentials for ranking protein models represented by idealized secondary structure elements., "
        "PLoS One, vol. 7, no. 11, p. e49242, Jan. 2012.\n"

        "\nFor membrane protein structure prediciton, please cite in addition:\n"
        "B. E. Weiner, N. Woetzel, M. Karakaş, N. Alexander, and J. Meiler, BCL::MP-Fold: Folding Membrane Proteins "
        "through Assembly of Transmembrane Helices., Structure, May 2013.\n"

        "\nFor EM-Fold please cite also:\n"
        "S. Lindert, N. Alexander, N. Wötzel, M. Karakaş, P. L. Stewart, and J. Meiler, EM-Fold: De Novo Atomic-Detail "
        "Protein Structure Determination from Medium-Resolution Density Maps, Structure, vol. 20, no. 3, pp. 464–478, "
        "Mar. 2012.\n"
        "S. Lindert, R. Staritzbichler, N. Wötzel, M. Karakaş, P. L. Stewart, and J. Meiler, EM-fold: De novo folding "
        "of alpha-helical proteins guided by intermediate-resolution electron microscopy density maps., Structure, vol. "
        "17, no. 7, pp. 990–1003, Jul. 2009.\n"

        "\nFor using BCL::Fold with NMR restraints please cite also:\n"
        "B. E. Weiner, N. Alexander, L. R. Akin, N. Woetzel, M. Karakaş, J. Meiler, BCL::Fold—Protein topology "
        "determination from limited NMR restraints., Proteins, pp. 1-9, Oct. 2013."
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL::Fold.\n"
        "\n"
        "Running BCL::Fold consists of five main steps.\n"
        "\n"
        "1) Create the fasta sequence file for the protein to be studied: "
        "You will need the protein sequence in fasta format and it should "
        "be stored in a <.fasta> file.  An example is given "
        "below.  For more information about fasta formats, please visit "
        "http://www.ncbi.nlm.nih.gov/Blast/fasta.shtml.\n"
        "\n"
        "2) Create the pool file. This could be done by hand or using the CreatePool "
        "application in the BCL. An example is given below.\n"
        "\n"
        "3) Create the optional stage file. This is only required if multiple stages, "
        "i.e. assembly and refinement are required.  See below for details.\n"
        "\n"
        "4) Create the optional score file. This file is a table listing scores and weights. "
        "Omit to use the default scores and weights.  See below for details.\n"
        "\n"
        "5) Run BCL::Fold:\n"
        "At a command prompt, navigate to the location of your BCL::Fold executable "
        "program.  The syntax for running the application looks like the following:\n"
        "\n"
        "bcl.exe Fold -native input/????A.pdb -pool_separate 1 -pool_min_sse_lengths 5 3 "
        "-quality RMSD GDT_TS -superimpose RMSD -message_level Critical "
        "-sspred JUFO PSIPRED -sspred_path_prefix input/ ?????A -stages_read input/stages.txt "
        "-pool input/????A.pool -loop_closure_threshold 0.1 -loop_rama_mutate_prob 0.0  "
        "-ccd_fraction \'[0.5,1.0]\' -nmodels 5 -protein_storage out/ -random_seed 1\n\n"
        "FLAGS:\n\n"
        "-native input/????A.pdb -> native PDB used for quality calculations\n"
        "-pool_min_sse_lengths 5 3 -> helices must be at least 5 residues and strands at least 3\n"
        "-quality RMSD GDT_TS -> use RMSD and GDT quality measures\n"
        "-superimpose RMSD -> superimpose models to native using RMSD\n"
        "-message_level Critical -> only critical message output\n"
        "-sspred JUFO PSIPRED -> use secondary structure predictions for extending/shrinking SSEs\n"
        "-sequence_data input/ ?????A -> directory and prefix for secondary structure predictions\n"
        "-stages_read input/stages.txt -> read stage file\n"
        "-pool input/????A.pool -> the pool file containing the secondary structure prediction of the protein. \n"
        "-loop_closure_threshold 0.1\n"
        "-loop_rama_mutate_prob 0.0\n"
        "-ccd_fraction \'[0.5,1.0]\'\n"
        "-nmodels 5 -> number of models that a single run creates (5 in that case)\n"
        "-protein_storage out/ -> directory that PDBs will be written to\n"
        "-random_seed 1 -> specify random seed\n\n"
        "INPUT AND OUTPUT.\n"
        "\n"
        "BCL::Fold requires two inputs, a) a fasta file, and b) a pool file.\n"
        "Additional optional files, such as c) stage file, d) score file, e) native PDB, or f) SSE predictions "
        "may also be included.\n"
        "\n"
        "a) The fasta file uses one letter codes for protein sequence and looks \n"
        "like the following:\n"
        "\n"
        ">1UBI:A|PDBID|CHAIN|SEQUENCE\n"
        "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG\n"
        "\n"
        "b) The pool file specifies the predicted secondary structure elements and looks \n"
        "like the following:\n"
        "\n"
        "bcl::assemble::SSEPool\n"
        "HELIX    1   1 TRP A    4  ALA A   19                                     16\n"
        "HELIX    2   2 TRP A    4  LEU A   36                                     33\n"
        "HELIX    3   3 LYS A    5  GLU A   18                                     14\n"
        "HELIX    4   4 LYS A    5  SER A   33                                     29\n"
        "HELIX    5   5 LYS A    7  THR A   22                                     16\n"
        "HELIX    6   6 GLU A   20  SER A   33                                     14\n"
        "HELIX    7   7 LEU A   21  LEU A   36                                     16\n"
        "HELIX    8   8 LEU A   41  TRP A   51                                     11\n"
        "\n"
        "c) The stage file specifies which protocols, scores, mutates, etc. to use, and has the following form:\n"
        "NUMBER_CYCLES 1                                          \n"
        "STAGE Stage_assembly                                     \n"
        "  TYPE MCM                                               \n"
        "  FOLD_PROTOCOLS Default Assembly                        \n"
        "  SCORE_PROTOCOLS Default                                \n"
        "  SCORE_WEIGHTSET_FILE assembly.scoreweights             \n"
        "  MUTATE_PROTOCOLS Default Assembly                      \n"
        "  NUMBER_ITERATIONS 5000 1000                            \n"
        "STAGE_END                                                \n"
        "STAGE Stage_refinement                                   \n"
        "  TYPE MCM                                               \n"
        "  FOLD_PROTOCOLS Default Refinement                      \n"
        "  SCORE_PROTOCOLS Default                                \n"
        "  SCORE_WEIGHTSET_FILE assembly.scoreweights             \n"
        "  MUTATE_PROTOCOLS Default Refinement                    \n"
        "  NUMBER_ITERATIONS 2000 400                             \n"
        "  PRINT_END_MODEL true                                   \n"
        "STAGE_END                                                \n"
        "STAGE Stage_loop_grow                                    \n"
        "  TYPE MCM                                               \n"
        "  FOLD_PROTOCOLS Default LoopCoordinateAdd               \n"
        "  SCORE_WEIGHTSET_FILE loop_add_coordinates.scoreweights \n"
        "  NUMBER_ITERATIONS 3500 600                             \n"
        "  MODIFY_START_MODEL true                                \n"
        "STAGE_END                                                \n"
        "STAGE Stage_close                                        \n"
        "  TYPE MCM                                               \n"
        "  FOLD_PROTOCOLS Default LoopClose                       \n"
        "  SCORE_WEIGHTSET_FILE loop_close.scoreweights           \n"
        "  NUMBER_ITERATIONS   5000     250                       \n"
        "  MODIFY_START_MODEL true                                \n"
        "STAGE_END                                                \n"
        "STAGE Stage_force_close                                  \n"
        "  TYPE MCM                                               \n"
        "  FOLD_PROTOCOLS Default LoopClose                       \n"
        "  SCORE_WEIGHTSET_FILE loop_close_force.scoreweights     \n"
        "  NUMBER_ITERATIONS   2500     125                       \n"
        "  MODIFY_START_MODEL true                                \n"
        "STAGE_END                                                \n"
        "\n"
        "d) The score file specifies used scores and should look like the following:\n"
        "\n"
        "bcl::storage::Table<double>    aaclash     aadist    aaneigh       aaneigh_ent  loop loop_closure       rgyr"
        "   sseclash ssepack_fr  strand_fr co_score ss_PSIPRED ss_PSIPRED_ent ss_JUFO ss_JUFO_ent\n"
        "weights                            500       0.35         50              50.0  10.0          500        5.0"
        "        500        8.0         20      0.5       20.0           20.0     5.0         5.0\n"
        "\n"
        "ADDITIONAL INFORMATION\n"
        "\n"
        "Helpful information concerning syntax and flags can be obtained by typing bcl.exe -help\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        "\n"
        + DefaultTechnicalSupportString() +
        "\n"
        "VIII. FUTURE DEVELOPMENT OF BCL::Fold.\n"
        + DefaultSectionSeparator() +
        "\n"
        "BCL::Fold is under ongoing further development.  For current research please refer to www.meilerlab.org "
        "and navigate to research.\n"
        + DefaultSectionSeparator() +
        "\n"
        "IX.  PROTOCOL-SPECIFIC README\n"
        + fold::Setup::GetReadme()
      );

      // end
      return s_readme;
    }

  } // namespace app
} // namespace bcl
