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
#include "fold/bcl_fold_default_flags.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_storage_file.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_extension.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "fold/bcl_fold_protocols.h"
#include "mc/bcl_mc_movie_printers.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "pdb/bcl_pdb_factory.h"
#include "quality/bcl_quality_measures.h"
#include "quality/bcl_quality_superimpose_measures.h"
#include "score/bcl_score_aa_neighborhood_exposure_prediction.h"
#include "score/bcl_score_symmetry.h"
#include "util/bcl_util_runtime_environments.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace fold
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DefaultFlags::s_Instance
    (
      GetObjectInstances().AddInstance( new DefaultFlags())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DefaultFlags::DefaultFlags()
    {
    }

    //! @brief Clone function
    //! @return pointer to new DefaultFlags
    DefaultFlags *DefaultFlags::Clone() const
    {
      return new DefaultFlags( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DefaultFlags::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////
  // flags //
  ///////////

    //! @brief return all command line flags defined by this class
    //! @return all command line flags defined by this class
    const util::ShPtrVector< command::FlagInterface> &DefaultFlags::GetAllFlags()
    {
      // initialize static ShPtrVector of FlagInterfaces to form the comment line
      static util::ShPtrVector< command::FlagInterface> s_all_flags_vector;

      // if the flag vector is initialize for the first time
      if( s_all_flags_vector.IsEmpty())
      {
        // flags from this class
        s_all_flags_vector.PushBack( GetProtocols().GetFlagProtocols());
        s_all_flags_vector.PushBack( GetProtocols().GetFlagMutateProtocols());
        s_all_flags_vector.PushBack( GetProtocols().GetFlagScoreProtocols());
        s_all_flags_vector.PushBack( GetFlagFastaRead());
        s_all_flags_vector.PushBack( GetFlagChainIdRead());
        s_all_flags_vector.PushBack( GetFlagPrefix());
        s_all_flags_vector.PushBack( GetFlagNumberModels());
        s_all_flags_vector.PushBack( GetFlagNativeModel());
        s_all_flags_vector.PushBack( GetFlagStartModel());
        s_all_flags_vector.PushBack( GetFlagUseNativeSSEsAsPool());
        s_all_flags_vector.PushBack( GetFlagEnableSSEResize());
        s_all_flags_vector.PushBack( GetFlagPoolSeparate());
        s_all_flags_vector.PushBack( GetFlagPrintMinimization());
        s_all_flags_vector.PushBack( GetFlagPrintTrackerHistory());
        s_all_flags_vector.PushBack( GetFlagMCNumberIterations());
        s_all_flags_vector.PushBack( GetFlagPrintStartModel());
        s_all_flags_vector.PushBack( GetFlagReadSequenceDataPath());
        s_all_flags_vector.PushBack( GetFlagScoreRead());
        s_all_flags_vector.PushBack( GetFlagScoreWrite());
        s_all_flags_vector.PushBack( GetFlagScoreWeightSetRead());
        s_all_flags_vector.PushBack( GetFlagScoreWeightSetWrite());
        s_all_flags_vector.PushBack( GetFlagMutateRead());
        s_all_flags_vector.PushBack( GetFlagMutateWrite());
        s_all_flags_vector.PushBack( GetFlagMutateWeightSetRead());
        s_all_flags_vector.PushBack( GetFlagMutateWeightSetWrite());
        s_all_flags_vector.PushBack( GetFlagStagesFileRead());
        s_all_flags_vector.PushBack( GetFlagStagesFileWrite());
        s_all_flags_vector.PushBack( GetFlagStagesNumberCycles());
        s_all_flags_vector.PushBack( GetFlagFitSwappedSSEs());
        s_all_flags_vector.PushBack( GetFlagPDBIDNumbering());

        // flags from other classes
        s_all_flags_vector.PushBack( pdb::Factory::GetFlagMinSSESize());
        s_all_flags_vector.PushBack( pdb::Factory::GetFlagAAClass());
        s_all_flags_vector.PushBack( assemble::SSEPool::GetFlagPoolRead());
        s_all_flags_vector.PushBack( assemble::SSEPool::GetFlagMinSSELengths());
        s_all_flags_vector.PushBack( mc::MoviePrinterInterface::GetFlagMoviePrinter());
        s_all_flags_vector.PushBack( mc::TemperatureAccepted::GetFlagTemperature());
        s_all_flags_vector.PushBack( quality::Measures::GetFlagQualityMeasures());
        s_all_flags_vector.PushBack( quality::SuperimposeMeasures::GetFlagSuperimposeMeasure());
        s_all_flags_vector.PushBack( score::Symmetry< assemble::ProteinModel>::GetFlagScoreSymmetry());
        s_all_flags_vector.PushBack( sspred::Methods::GetFlagReadSSPredictions());
        s_all_flags_vector.PushBack( assemble::ProteinStorageFile::GetDefaultStorageFlag());
        s_all_flags_vector.PushBack( assemble::SSEPool::GetFlagPoolPrefix());
        s_all_flags_vector.PushBack( score::AANeighborhoodExposurePrediction::GetFlagScoreExposure());
      }

      // return
      return s_all_flags_vector;
    }

  ///////////////////
  // general flags //
  ///////////////////

    //! @brief return command line flag for writing the minimization to files - only of given step statuses
    //! @return command line flag for writing the minimization to files - only of given step statuses
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagPrintMinimization()
    {
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "write_minimization", "write the minimization to files - only of given step statuses",
          command::Parameter
          (
            "stepstatuses",
            "any step status from the list",
            command::ParameterCheckAllowed( opti::StepStatusEnum::GetStringVector()),
            opti::GetStepStatusName( opti::e_Improved)
          ),
          0,
          opti::s_NumberStepStatus
        )
      );

      return s_flag;
    }

    //! @brief return command line parameter for specifying the path for where tracker files are generated
    //! @return command line parameter for specifying the path for where tracker files are generated
    util::ShPtr< command::ParameterInterface> &DefaultFlags::GetParameterPrintTrackerHistoryPath()
    {
      static util::ShPtr< command::ParameterInterface> s_parameter
      (
        new command::Parameter
        (
          "print_tracker_history_path",
          "path for where the tracker history is saved",
          ""
        )
      );

      return s_parameter;
    }

    //! @brief return command line flag for a detailed analysis of each step in a tabulated format
    //! @return command line flag for a detailed analysis of each step in a tabulated format
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagPrintTrackerHistory()
    {
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "print_tracker_history",
          "\tprint <path>" + util::GetRuntimeEnvironment().GetPathSeperator()
            + "<prefix>.tracker_history that has a detailed analysis of each step in a tabulated format"
        )
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameter into flag
        flag->PushBack( GetParameterPrintTrackerHistoryPath());
      }

      return s_flag;
    }

    //! @brief return command line flag for monte carlo minimization, max number of rejected steps and max iterations
    //! @return command line flag for monte carlo minimization, max number of rejected steps and max iterations
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagMCNumberIterations()
    {
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "mc_number_iterations",
          "\tmodify the number of iterations and max steps without improvement"
        )
      );
      static util::ShPtr< command::ParameterInterface> s_total
      (
        new command::Parameter
        (
          "mc_number_total_iterations", "\ttotal number of iterations for minimization", "100"
        )
      );
      static util::ShPtr< command::ParameterInterface> s_max_unimproved
      (
        new command::Parameter
        (
          "mc_max_unimproved_steps", "\tnumber of steps without improvement before terminating", "50"
        )
      );

      // if the flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);
        // insert parameters into flag
        flag->PushBack( s_total);
        flag->PushBack( s_max_unimproved);
      }

      // end
      return s_flag;
    }

    //! @brief return command line flag for specifying whether or not to print the start model
    //! @return command line flag for specifying  whether or not to print the start model
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagPrintStartModel()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic( "print_start_model", "\t\tflag to determine whether to print start model or not")
      );

      // end
      return s_flag;
    }

    //! @brief return command line flag for specifying a prefix to be used for writing files
    //! @return command line flag for specifying a prefix to be used for writing files
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagPrefix()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "prefix", "\tfile prefix to be used for writing files",
          command::Parameter( "file_prefix", "\twill be prepended to each file - can contain a path", "")
        )
      );

      // end
      return s_flag;
    }

    //! @brief return command line flag for specifying the number of models to build
    //! @return command line flag for specifying the number of models to build
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagNumberModels()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "nmodels", "the number of models that will be created",
          command::Parameter
          (
            "number_of_models", "\tthe number of models that will be created", "1"
          )
        )
      );

      // end
      return s_flag;
    }

    //! @brief return command line flag for providing one or more fasta files for complete de novo folding
    //! @return command line flag for providing one or more fasta files for complete de novo folding
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagFastaRead()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "fasta",
          "\t filenames of one or more fasta files to be used, one file per chain in alphabetical order",
          command::Parameter
          (
            "fasta_filename",
            "\tfilename for input amino acid sequence of the form {name}{chainid}.fasta, or {name}.fasta and chain id given by chain id flag",
            command::ParameterCheckFileExistence()
          ),
          0, 20
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for providing one or more chain ids for fasta files for complete de novo folding
    //! @return command line flag for providing one or more chain ids for fasta files for complete de novo folding
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagChainIdRead()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "fasta_chain_id",
          "\t chain ids for each amino acid sequence given by fasta files",
          command::Parameter
          (
            "fasta_filename",
            "\tfilename single character A to Z for each amino acid sequence"
          ),
          0, 20
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for providing a native model for comparison
    //! @return command line flag for providing a native model for comparison
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagNativeModel()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "native", "\t\tpdb file of native structure or template structure for comparison",
          command::Parameter
          (
            "native_pdb_filename", "\tfilename for native or template pdb to which rmsd will be calculated",
             command::ParameterCheckExtension( ".pdb"), ""
          )
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for providing a starting model
    //! @return command line flag for providing a starting model
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagStartModel()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "start_model", "\t\tpdb file of starting model structure for refinement or folding",
          command::Parameter
          (
            "start_model_filename", "\tfilename for pdb file of starting model structure",
             command::ParameterCheckExtension( ".pdb"), ""
          )
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for using native SSE definitions as the SSE pool
    //! @return command line flag for using native SSE definitions as the SSE pool
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagUseNativeSSEsAsPool()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "use_native_pool",
          "\tflag to use the SSE definitions in the native structure as the pool",
          command::Parameter
          (
            "idealize",
            "idealize the SSEs or use the native conformations",
            command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "ideal", "native")),
            "ideal"
          )
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for enable resizing of SSEs
    //! @return command line flag for enable resizing of SSEs
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagEnableSSEResize()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "enable_sse_resize", "\tflag to allow resizing and split moves for SSEs. "
          "Requires a high (>20) weight for sspred_PSIPRED and similar"
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for separating adjoining SSE pools with specified number of loop residues
    //! @return command line flag for separating adjoining SSE pools with specified number of loop residues
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagPoolSeparate()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "pool_separate", "\t\tseparate the adjoining SSEs in the pool with specified number of loop residues",
          command::Parameter
          (
            "pool_separate_nr_res", "\tnumber of residues each adjoining SSE will be shortened by",
             command::ParameterCheckRanged< size_t>( 1, 3),
             "1"
          )
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for specifying the path where the ss predictions should be read from
    //! @return command line flag for specifying the path where the ss predictions should be read from
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagReadSequenceDataPath()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "sequence_data",
          "\tflag for specifying the path and prefix for sequence data files, such as secondary structure prediction"
          ": {data_path} {data_prefix};"
          "to be used in conjunction with sspred flag. The {data_prefix} should be identifier for your protein ( "
          " most likely the *FOUR* letter code)."
          "The chain id for secondary structure predictions will be determined from the chains that are in your protein"
          " model. A data file of the specified types will be searched for for each chain as "
          "\"{data_path}{data_prefix}{chainid}.{extension}\"."
        )
      );

      // initialize parameters
      static util::ShPtr< command::ParameterInterface> s_path_param
      (
        new command::Parameter
        (
          "data_path",
          "\tpath of sequence data file",
          "."
        )
      );
      static util::ShPtr< command::ParameterInterface> s_prefix_param
      (
        new command::Parameter
        (
          "data_prefix",
          "\tprefix of sequence data file",
          "XXXX"
        )
      );

      // if this flag is initialized for the first time
      if( s_flag->GetParameterList().IsEmpty())
      {
        util::ShPtr< command::FlagStatic> flag( s_flag);

        // insert parameters
        flag->PushBack( s_path_param);
        flag->PushBack( s_prefix_param);
      }

      // end
      return s_flag;
    }

    //! @brief return command line flag for reading scores from a file
    //! @return command line flag for reading scores from a file
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagScoreRead()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "score_read", "\t\tread scoring function from a given file",
          command::Parameter( "score_read_filename", "\tfilename for the input scoring function", "")
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for writing scores to a file
    //! @return command line flag for writing scores to a file
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagScoreWrite()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "score_write", "\t\twrite scoring function to a file",
          command::Parameter( "score_write_filename", "\tfilename for writing scoring function", "")
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for reading score weight set from a file
    //! @return command line flag for reading score weight set from a file
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagScoreWeightSetRead()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "score_weightset_read", "\t\tfile with weightset for the scores to be used",
          command::Parameter( "score_weightset_read_filename", "\tfilename for score weightset", "")
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for writing score weight set to a file
    //! @return command line flag for writing score weight set to a file
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagScoreWeightSetWrite()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "score_weightset_write", "\t\twrite score weightset to a file",
          command::Parameter( "score_weightset_filename_prefix", "\tfilename prefix for score weightset", "")
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for reading mutates from a file
    //! @return command line flag for reading mutates from a file
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagMutateRead()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "mutate_read", "\t\tread mutate objects from a given file",
          command::Parameter( "mutate_read_filename", "\tfilename for reading mutate objects", "")
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for writing mutates to a file
    //! @return command line flag for writing mutates to a file
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagMutateWrite()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "mutate_write", "\t\twrite mutate objects to a file",
          command::Parameter( "mutate_write_filename", "\tfilename for writing mutate objects", "")
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for reading mutate weight set from a file
    //! @return command line flag for reading mutate weight set from a file
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagMutateWeightSetRead()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "mutate_weightset_read", "\t\tfile with weightset for the mutates to be used",
          command::Parameter( "mutate_weightset_filename", "\tfilename for mutate weightset", "")
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for writing mutate weight set to a file
    //! @return command line flag for writing mutate weight set to a file
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagMutateWeightSetWrite()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "mutate_weightset_write", "\t\twrite mutate weightset to a file",
          command::Parameter( "mutate_weightset_filename_prefix", "\tfilename prefix for mutate weightset", "")
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for defining the filename to read the stages from
    //! @return command line flag for defining the filename to read the stages from
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagStagesFileRead()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "stages_read", "\t\tfile with stages",
          command::Parameter( "stages_filename", "\tfilename for stages", "")
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for defining the filename to write the stages to
    //! @return command line flag for defining the filename to write the stages to
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagStagesFileWrite()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "stages_write", "\t\tfile with stages to be used",
          command::Parameter( "stage_filename", "\tfilename for stages", "")
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for defining the number of cycles stages should go through
    //! @return command line flag for defining the number of cycles stages should go through
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagStagesNumberCycles()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "stages_cycle", "\t\t the number of cycles stages should go through",
          command::Parameter
          (
            "stage_number_cycles", "\tthe number of cycles",
            command::ParameterCheckRanged< size_t>( 1, 100), "1"
          )
        )
      );
      // end
      return s_flag;
    }

    //! @brief return command line flag for defining whether to fit swapped sses
    //! @return command line flag for defining whether to fit swapped sses
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagFitSwappedSSEs()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic( "fit_swapped", "\t\t whether to fit swapped SSEs to the original")
      );
      // end
      return s_flag;
    }

    //! @brief command line flag indicating the input file numbering is PDB numbering instead of sequence id numbering
    //! @return flag indicating the input file numbering is PDB numbering instead of sequence id numbering
    util::ShPtr< command::FlagInterface> &DefaultFlags::GetFlagPDBIDNumbering()
    {
      // initialize static flag
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "use_pdbid_numbering",
          "\tindicates the input file (e.g. restraint or domain files) numbering is PDB numbering (from the atom section) instead of sequence id numbering"
        )
      );

      return s_flag;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DefaultFlags::Read( std::istream &ISTREAM)
    {
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DefaultFlags::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // return the stream
      return OSTREAM;
    }

  } // namespace fold

} // namespace bcl
