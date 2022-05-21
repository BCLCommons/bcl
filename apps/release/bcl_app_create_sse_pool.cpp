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
#include "bcl_app_create_sse_pool.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_factories.h"
#include "assemble/bcl_assemble_sse_pool_agreement.h"
#include "biol/bcl_biol_aa_sequence_factory.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "sspred/bcl_sspred_sse_factory_threshold.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! default constructor
    CreateSSEPool::CreateSSEPool() :
      m_PrefixFlag
      (
        new command::FlagStatic
        (
          "prefix",
          "\tflag for providing a prefix for fasta  and sspredictions from which the pool is going to be generated",
          command::Parameter
          (
            "prefix_param",
            "\tfprefix for fasta  and sspredictions from which the pool is going to be generated"
          )
        )
      ),
      m_OutputPrefixFlag
      (
        new command::FlagStatic
        (
          "output_prefix",
          "\tflag for providing an output prefix for files to be written such as pool",
          command::Parameter
          (
            "output_prefix_param",
            "\tfprefix for output files to be written such as pool",
            ""
          )
        )
      ),
      m_ChainIdFlag
      (
        new command::FlagStatic
        (
          "chain_id",
          "\tflag to indicate, that the prefix does end with the chain id - otherwise, assume chain id A"
        )
      ),
      m_PdbFlag
      (
        new command::FlagStatic
        (
          pdb::GetDefaultFileExtension(),
          "\tflag for providing a real pdb filename so the real sse definitions are also outputted",
          command::Parameter
          (
            "pdb_filename",
            "\tfilename for the pdb to be read so the real sse definitions are also outputted",
            "default.pdb"
          )
        )
      ),
      m_PdbPoolFileOutParam
      (
        new command::Parameter
        (
          "pdb_pool_filename",
          "filename for pool generated from the given pdb file",
          ""
        )
      ),
      m_EvaluatePoolFlag
      (
        new command::FlagStatic
        (
          "evaluate_pool",
          "flag for evaluating a given pool in comparison to the SSE definitions from the given pdb",
          command::Parameter
          (
            "evaluate_pool_file", "filename of the pool to be evaluated", "test.pool"
          )
        )
      ),
      m_SsMethods
      (
        new command::FlagDynamic
        (
          "ssmethods",
          "\tone or more ssmethods to be used in generation of the pool",
          command::Parameter
          (
            "ssmethod",
            "any ssmethod from the list",
            command::ParameterCheckEnumerate< sspred::Methods>()
          ),
          1,
          sspred::GetMethods().GetEnumCount()
        )
      ),
      m_ChopSsesFlag
      (
        new command::FlagStatic
        (
          "chop_sses",
          "\tgive flag for chopping long sses"
        )
      ),
      m_SystematicExtendFlag
      (
        new command::FlagStatic
        (
          "extend_sses",
          "\tgive flag for adding systematically extended copies of sses to pool"
        )
      ),
      m_SseThresholdFlag
      (
        new command::FlagStatic
        (
          "sse_threshold",
          "\tthreshold to be used to identify whether a residue in SSE or not, between 0 and 1"
        )
      ),
      m_HelixThresholdParam
      (
        new command::Parameter
        (
          "helix_threshold_value",
          "\tthreshold to be used to identify whether a residue in helix or not, between 0 and 1",
          command::ParameterCheckRanged< double>( 0.0, 1.0),
          "0.5"
        )
      ),
      m_StrandThresholdParam
      (
        new command::Parameter
        (
          "strand_threshold_value",
          "\tthreshold to be used to identify whether a residue in strand or not, between 0 and 1",
          command::ParameterCheckRanged< double>( 0.0, 1.0),
          "0.5"
        )
      ),
      m_CoilThresholdParam
      (
        new command::Parameter
        (
          "coil_threshold_value",
          "\tthreshold to be used to identify whether a residue in coil or not, between 0 and 1",
          command::ParameterCheckRanged< double>( 0.0, 1.0),
          "0.5"
        )
      ),
      m_FactoryFlag
      (
        new command::FlagStatic
        (
          "factory",
          "use SSE factory of choice"
        )
      ),
      m_FactoryParam
      (
        new command::Parameter
        (
          "factory",
          "the factory of choice",
          command::ParameterCheckEnumerate< assemble::SSEFactories>(),
          assemble::GetSSEFactories().e_PredictionThreshold.GetName()
        )
      ),
      m_FactoryStreamParam
      (
        new command::Parameter
        (
          "factory_file",
          "bcl file to read factory from; parameters given on commandline will overwrite members that were read eventually",
          command::ParameterCheckFileExistence(),
          ""
        )
      ),
      m_JoinSeparateFlag
      (
        new command::FlagStatic
        (
          "join_separate",
          "joins sses below given lengths and separates directly adjacent sses by inserting a two residue coil the\n"
          "given; joining is only done, if the sses are shorter than the min sizes given; separation is only done, if\n"
          "the resulting sses are still long enough"
        )
      )
    {
      m_PdbFlag->PushBack( m_PdbPoolFileOutParam);
      m_SseThresholdFlag->PushBack( m_HelixThresholdParam);
      m_SseThresholdFlag->PushBack( m_StrandThresholdParam);
      m_SseThresholdFlag->PushBack( m_CoilThresholdParam);
      m_FactoryFlag->PushBack( m_FactoryParam);
      m_FactoryFlag->PushBack( m_FactoryStreamParam);
    }

    // add this app to the protein group only; as its input is tailored only for further use in the protein:Fold workflow
    const ApplicationType CreateSSEPool::CreateSSEPool_Instance
    (
      GetAppGroups().AddAppToGroup( new CreateSSEPool(), GetAppGroups().e_Protein)
    );

    //! @brief Clone function
    //! @return pointer to new CreateSSEPool
    CreateSSEPool *CreateSSEPool::Clone() const
    {
      return new CreateSSEPool( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &CreateSSEPool::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string CreateSSEPool::GetDescription() const
    {
      return "Generates a pool of secondary structure elements (SSEs) used by other potein applications including Fold";
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief the Main function
    //! @return error code - 0 for success
    int CreateSSEPool::Main() const
    {

    ////////////////
    // initialize //
    ////////////////

      // create storage map with the sse thresholds for helix and strand
      m_ThresholdsMap[ biol::GetSSTypes().HELIX] = m_HelixThresholdParam->GetNumericalValue< double>();
      m_ThresholdsMap[ biol::GetSSTypes().STRAND] = m_StrandThresholdParam->GetNumericalValue< double>();
      m_ThresholdsMap[ biol::GetSSTypes().COIL] = m_CoilThresholdParam->GetNumericalValue< double>();

      // min size
      m_SSEMinSizeMap = assemble::SSEPool::GetCommandLineMinSSELengths();

      // initialize storage::Set of SSMethods to be used
      const storage::Set< sspred::Method> ss_methods( m_SsMethods->GetObjectSet< sspred::Method>());

      //initialize write and read stream objects
      io::OFStream write;
      io::IFStream read;

    ////////////////
    // read fasta //
    ////////////////

      // get fasta filename
      std::string prefix( m_PrefixFlag->GetFirstParameter()->GetValue());

      // parse the store the fasta path and tag
      const storage::VectorND< 2, std::string> fasta_path_and_tag
      (
        pdb::Handler::ExtractPathAndPDBTag( prefix + ".fasta")
      );
      std::string fasta_tag( fasta_path_and_tag.Second());
      const std::string fasta_path( fasta_path_and_tag.First());
      char chain_id( 'A');

      // if the chain id was provided
      if( m_ChainIdFlag->GetFlag())
      {
        // set the chain id from the prefix
        chain_id = fasta_tag[ fasta_tag.length() - 1];
        fasta_tag = fasta_tag.substr( 0, fasta_tag.length() - 1);
      }

      const std::string fasta_filename( prefix + ".fasta");

      // open fasta file for reading
      BCL_MessageStd( "read fasta file from " + fasta_filename);
      io::File::MustOpenIFStream( read, fasta_filename);

      // initialize sequence w/ fasta
      biol::AASequence sequence
      (
        biol::AASequenceFactory::BuildSequenceFromFASTA( read, biol::GetAAClasses().e_AA, chain_id)
      );

      // reset stream
      io::File::CloseClearFStream( read);

      // write the fasta
      BCL_MessageStd( "sequence read: ");
      sequence.WriteFasta( util::GetLogger());

    ////////////////////////
    // read SSPredictions //
    ////////////////////////

      // read the secondary structure predictions
      BCL_MessageStd( "Reading the ss methods");
      BCL_Assert
      (
        sspred::MethodHandler::ReadPredictionsForAASequence( ss_methods, sequence, fasta_tag, fasta_path),
        "reading ss predictions has failed!"
      );

      // sse pool from pdb
      assemble::SSEPool sse_pool_from_pdb;

      BCL_UserAssert
      (
        m_FactoryFlag->GetFlag() || m_EvaluatePoolFlag->GetFlag()
          || ( m_PdbFlag->GetFlag() && m_PdbPoolFileOutParam->GetWasSetInCommandLine()),
        "Please select one flag of -" + m_FactoryFlag->GetName() + " or -" + m_EvaluatePoolFlag->GetName()
          + " or -" + m_PdbFlag->GetName() + " " + m_PdbPoolFileOutParam->GetName()
      );
      
      if( m_PdbFlag->GetFirstParameter()->GetWasSetInCommandLine())
      {
        // issue message
        BCL_MessageStd( "Reading the provided pdb file to get real pool");

        const std::string pdb_filename( m_PdbFlag->GetFirstParameter()->GetValue());

        // initialize the pdb handler and factory
        io::File::MustOpenIFStream( read, pdb_filename);
        pdb::Handler pdb_handler( read);
        io::File::CloseClearFStream( read);

        pdb::Factory pdb_factory( biol::GetAAClasses().e_AA);

        // read in the sequences from pdb file
        assemble::ProteinModel protein_model( pdb_factory.ProteinModelFromPDB( pdb_handler));

        // create pool from sse definitions
        sse_pool_from_pdb = assemble::SSEPool( protein_model.GetSSEs());
        ProcessPool( sse_pool_from_pdb);

        // print out the real pool
        if( m_PdbPoolFileOutParam->GetWasSetInCommandLine())
        {
          io::File::MustOpenOFStream( write, m_PdbPoolFileOutParam->GetValue());
          BCL_MessageStd( "writing native pool to file: " + m_PdbPoolFileOutParam->GetValue());
          sse_pool_from_pdb.WriteSSEPool( write);
          io::File::CloseClearFStream( write);
        }
      }
      
      // sse pool generated from factory
      assemble::SSEPool sse_pool_factory;

      // if factory flag is given
      if( m_FactoryFlag->GetFlag())
      {

      /////////////////////
      // Create the pool //
      /////////////////////

        BCL_MessageStd( "Creating the sse pool from provided secondary structure information");

        std::string file_extension( ".");

        const assemble::SSEFactory sse_factory_enum( m_FactoryParam->GetValue());
        file_extension += sse_factory_enum.GetName();

        util::ShPtr< assemble::SSEFactoryInterface> sse_factory; // create before entering the loop
        
        for
        (
          storage::Set< sspred::Method>::const_iterator met_itr( ss_methods.Begin()), met_itr_end( ss_methods.End());
          met_itr != met_itr_end;
          ++met_itr
        )
        {
          // factory from file
          if( m_FactoryStreamParam->GetWasSetInCommandLine())
          {
            io::IFStream read;
            io::File::MustOpenIFStream( read, m_FactoryStreamParam->GetValue());
            sse_factory = util::ShPtr< assemble::SSEFactoryInterface>::GetShPtrToNewObjectFromStream( read);
            io::File::CloseClearFStream( read);
          }
          // factory from enum given in commandline
          else
          {
            sse_factory = assemble::GetSSEFactories().Create( sse_factory_enum, *met_itr, m_ThresholdsMap);
          }

          // systematic extand for threshold factory
          if( sse_factory_enum == assemble::GetSSEFactories().e_PredictionThreshold)
          {
            util::ShPtr< sspred::SSEFactoryThreshold>( sse_factory)->SetExtendSSEs( m_SystematicExtendFlag->GetFlag());
          }

          assemble::SSEPool current_pool( sse_factory->operator()( sequence));
          sse_pool_factory.InsertElements( current_pool);

          file_extension += '_' + met_itr->GetName(); // add each method use to the pool filename
        }

        ProcessPool( sse_pool_factory);

        // output pool to file
        std::string pool_filename
        (
          m_OutputPrefixFlag->GetFirstParameter()->GetValue() + fasta_tag + chain_id + file_extension + ".pool"
        );
        BCL_MessageStd( "Writing the pool to file: " + pool_filename);
        io::File::MustOpenOFStream( write, pool_filename);
        sse_pool_factory.WriteSSEPool( write);
        io::File::CloseClearFStream( write);
      }

      if( m_EvaluatePoolFlag->GetFlag())
      {
        // calculate agreement
        assemble::SSEPoolAgreement agreement_function;

        if( m_EvaluatePoolFlag->GetFirstParameter()->GetWasSetInCommandLine())
        {
          // read the pool from given file
          io::File::MustOpenIFStream( read, m_EvaluatePoolFlag->GetFirstParameter()->GetValue());
          assemble::SSEPool sse_pool_evaluate;
          sse_pool_evaluate.ReadSSEPool( read, sequence, 0, 0);
          io::File::CloseClearFStream( read);
          ProcessPool( sse_pool_evaluate);

          if( !sse_pool_from_pdb.IsEmpty())
          {
            const double agreement_pdb_evaluate( agreement_function( sse_pool_from_pdb, sse_pool_evaluate));
            BCL_MessageStd( "agreement between pdb and evaluate pool: " + util::Format()( agreement_pdb_evaluate));
            const double q3_pdb_evaluate( agreement_function.Q3Score( sse_pool_from_pdb, sse_pool_evaluate));
            BCL_MessageStd( "Q3 between pdb and evaluate pool: " + util::Format()( q3_pdb_evaluate));
            BCL_MessageStd( "pool analysis compared to given native model");
            sse_pool_evaluate.CalculateStatistics( sse_pool_from_pdb).GetTransposedTable().WriteFormatted( util::GetLogger());
          }

          if( !sse_pool_factory.IsEmpty())
          {
            const double agreement_factory_evaluate( agreement_function( sse_pool_factory, sse_pool_evaluate));
            BCL_MessageStd( "agreement between factory and evaluate pool: " + util::Format()( agreement_factory_evaluate));
            const double q3_factory_evaluate( agreement_function.Q3Score( sse_pool_factory, sse_pool_evaluate));
            BCL_MessageStd( "Q3 between factory and evaluate pool: " + util::Format()( q3_factory_evaluate));
          }
        }
        const double agreement_factory_pdb( agreement_function( sse_pool_factory, sse_pool_from_pdb));
        BCL_MessageStd( "agreement between factory and pdb pool: " + util::Format()( agreement_factory_pdb));
        const double q3_factory_pdb( agreement_function.Q3Score( sse_pool_factory, sse_pool_from_pdb));
        BCL_MessageStd( "Q3 between factory and pdb pool: " + util::Format()( q3_factory_pdb));
        BCL_MessageStd( "pool analysis of the factory generated pool to compared to given native model");
        sse_pool_factory.CalculateStatistics( sse_pool_from_pdb).GetTransposedTable().WriteFormatted( util::GetLogger());
      }

      return 0;
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &CreateSSEPool::GetReadMe() const
    {
      // create a static string to hold readme information
      static const std::string s_readme_text
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::CreateSSEPool, terms of use, appropriate citation, installation "
        "procedures, BCL::CreateSSEPool execution, technical support, and future research directions.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::CreateSSEPool?\n"
        "\n"
        "BCL::CreateSSEPool is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is"
        "part of a larger library of applications called BCL::Commons.  BCL::CreateSSEPool generates the secondary"
        "structure pool for given proteins using individual and consensus secondary structure predictions.  The "
        "application first reads a fasta file as input, then reads secondary structure predication parameters entered"
        "in the commandline, and finally creates the secondary structure pool.  The Output is an ASCII files that shows"
        "the predicted SSEs based on input paramaters."
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        "\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::Jufo.\n"
        "\n"
        "When using BCL::CreateSSEPool in a publication, please cite the following publications describing the "
        "application's development, which is currently under review.  Any news will be published at "
        "http://www.meilerlab.org/\n"
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        "1) BCL::CreateSSEPool:\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL::CreateSSEPool.\n"
        "\n"
        "Running BCL::CreateSSEPool consists of three main steps.\n"
        "\n"
        "1) Create the fasta sequence file for the protein to be studied:\n"
        "You will need the protein sequence in fasta format, and it should be stored "
        "in a <.fasta> file.  An example is given below.  For more information about fasta formats, please visit "
        "http://www.ncbi.nlm.nih.gov/Blast/fasta.shtml.\n"
        "\n"
        "2) Create the Secondary Structure Prediction Files or use PDB format for protein to be studied:\n"
        "\n"
        "3) Run BCL::CreateSSEPool:\n"
        "At a command prompt, navigate to the location of your BCL::CreateSSEPool executable program.  The syntax for "
        "running the application looks like the following:\n"
        "\n"
        "bcl.exe CreateSSEPool -prefix <protein> -ssmethods <sse prediction type> -factory\n"
        "\n"
        "BCL::CreateSSEPool needs a fasta sequence with the extension <.fasta> and a secondary structure element (sse)"
        "file type to exist in the same directory.\n"
        "The allowed file types and parameters -help flag. \n"
        "INPUT AND OUTPUT.\n"
        "\n"
        "BCL::CreateSSEPool requires two inputs, a fasta file and a corresponding sse prediction type. "
        "The fasta file uses one letter codes for protein sequence and looks like the following:\n"
        "\n"
        ">1UBI:A|PDBID|CHAIN|SEQUENCE\n"
        "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG\n"
        "\n"
        "The allowed sse prection type files are listed in the -help flag. \n"
        "\n"
        "The output file is formatted as following:\n"
        "\n"
        " HELIX    1   1 GLY A   98  ALA A  108  1                                  11\n"
        " HELIX    2   2 GLY A   98  LYS A  120  1                                  23\n"
        " SHEET    3     GLY A  100  SER A  106                                       \n"
        " HELIX    4   4 ASP A  110  LYS A  120  1                                  11\n"
        " HELIX    5   5 ASP A  113  ARG A  123  1                                  11\n"
        " HELIX    6   6 SER A  127  ASN A  138  1                                  12\n"
        " HELIX    7   7 SER A  127  LEU A  150  1                                  24\n"
        " SHEET    8     SER A  128  ASN A  131                                       \n"
        " HELIX    9   9 SER A  133  ILE A  142  1                                  10\n"
        "\n"
        "The individual columns represent following:\n"
        "Column  1: Secondary structure type\n"
        "Column  2: Secondary structure number\n"
        "Column  3: Helix number\n"
        "Column  4: Three letter amino acid code at first position of SSE\n"
        "Column  5: Chain type\n"
        "Column  6: Starting Amino Acid position of SSE\n"
        "Column  7: Three letter amino acid code at last position of SSE\n"
        "Column  8: Chain type\n"
        "Column  9: Ending amino acid position of SSE\n"
        "Column 10: Helix type\n"
        "Column 11: Length of SSE\n"
        "\n"
        "ADDITIONAL INFORMATION\n"
        "\n"
        "Helpful information concerning syntax and flags can be obtained by typing bcl.exe CreateSSEPool -help\n"
        "\n"
        "For more general information about the product, type bcl.exe CreateSSEPool -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        "\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::CreateSSEPool.\n"
        "\n"
        "BCL::CreateSSEPool is under ongoing further development.  For current research please refer to "
        "www.meilerlab.org and navigate to research.\n"
        + DefaultSectionSeparator()
      );

      // return readme information
      return s_readme_text;
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> CreateSSEPool::InitializeCommand() const
    {
      // initialize ShPtr to a new command object
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // fasta flag
      sp_cmd->AddFlag( m_PrefixFlag);

      // output flag
      sp_cmd->AddFlag( m_OutputPrefixFlag);

      // chain id flag
      sp_cmd->AddFlag( m_ChainIdFlag);

      sp_cmd->AddFlag( m_PdbFlag);

      // evaluate_pool flag and param
      sp_cmd->AddFlag( m_EvaluatePoolFlag);

      // minimum pool sse lengths flag
      sp_cmd->AddFlag( assemble::SSEPool::GetFlagMinSSELengths());

      // ssmethods flag and parameters
      sp_cmd->AddFlag( m_SsMethods);

      // flag for chopping sses
      sp_cmd->AddFlag( m_ChopSsesFlag);

      // flag for extending sses
      sp_cmd->AddFlag( m_SystematicExtendFlag);

      // adjust minimal sse lengths
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());

      // threshold
      sp_cmd->AddFlag( m_SseThresholdFlag);

      // mc factory
      sp_cmd->AddFlag( m_FactoryFlag);

      // join separate
      sp_cmd->AddFlag( m_JoinSeparateFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief process pool by joining, separating, chopping ...
    //! @param POOL the pool to process
    void CreateSSEPool::ProcessPool( assemble::SSEPool &POOL) const
    {
      if( m_JoinSeparateFlag->GetFlag())
      {
        // join adjacent that are too short
        while( POOL.Join( m_SSEMinSizeMap));

        // separate adjacent, that are long enough, but should be separated by a coil
        POOL.Separate( m_SSEMinSizeMap, 1);
      }

      // prune of short sses
      POOL.Prune( m_SSEMinSizeMap);

      // chop
      if( m_ChopSsesFlag->GetFlag())
      {
        POOL.ChopSSEs( m_SSEMinSizeMap);
      }
    }

  } // namespace app
} // namespace bcl
