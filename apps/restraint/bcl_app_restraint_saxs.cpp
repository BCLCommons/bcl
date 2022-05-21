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
#include "bcl_app_restraint_saxs.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_sasa_data.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "fold/bcl_fold_add_parabolic_loops.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "restraint/bcl_restraint_sas_analysis.h"
#include "restraint/bcl_restraint_sas_debye.h"
#include "restraint/bcl_restraint_sas_debye_interface.h"
#include "restraint/bcl_restraint_sas_density_data.h"
#include "restraint/bcl_restraint_sas_optimization.h"
#include "restraint/bcl_restraint_saxs_data_reduction.h"
#include "restraint/bcl_restraint_saxs_opti_result.h"
#include "score/bcl_score_sas_type.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RestraintSaxs::RestraintSaxs() :
      m_PDBFilenameFlag
      (
        new command::FlagStatic
        (
          "pdb",
          "\tpdb file for which restraints will be created",
          command::Parameter
          (
            "pdb_filename",
            "\tfilename for which restraints will be calculated",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_ExpRestraintFileFlag
      (
        new command::FlagStatic
        (
          "exp_data",
          "\texperimental data file for given pdb",
          command::Parameter
          (
            "data_filename",
            "\tfilename containing experimental data",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_OutputFileFlag
      (
        new command::FlagStatic
        (
          "output_file",
          "path and name of the output file which will hold the simulated restraints.",
          command::Parameter
          (
            "string",
            "\tpath and name of the outputfile",
            "saxs.data"
          )
        )
      ),
      m_SasaFileFlag
      (
        new command::FlagStatic
        (
          "sasa_data",
          "\tSurface Accessible Solvent Area data from MSMS for given pdb",
          command::Parameter
          (
            "sasa_filename",
            "\tfilename containing sasa data",
            ""
          )
        )
      ),
      m_PrintRMSD
      (
        new command::FlagStatic
        (
          "rmsd",
          "pass this flag to print rmsd scores"
        )
      ),
      m_PrintStovgaard
      (
        new command::FlagStatic
        (
          "stovgaard",
          "pass this flag to print stovgaard scores"
        )
      ),
      m_PrintProteinModelWithLoopsFileFlag
      (
        new command::FlagStatic
        (
          "print_loops",
          "pass this flag to print out the protein model with loops to a filename",
          command::Parameter
          (
            "filename",
            "\tpath and name of the pdb file to output with loops",
            "loop_approximation.pdb"
          )
        )
      ),
      m_PrintDerivative
      (
        new command::FlagStatic
        (
          "print_derivative",
          "pass this flag to print out derivative to a filename",
          command::Parameter
          (
            "filename",
            "\tpath and name of the text file to output derivative",
            "derivative.data"
          )
        )
      ),
      m_ReduceData
      (
        new command::FlagStatic
        (
          "reduce_data",
          "pass this flag if you want to use Shannon Sampling to reduce experimental dataset"
        )
      ),
      m_UseErrors
      (
        new command::FlagStatic
        (
          "use_errors",
          "pass this flag use Experimental/Simulated Errors in Chi Calculation"
        )
      ),
      m_ApproximateSideChains
      (
        new command::FlagStatic
        (
          "apx_sc",
          "pass this flag to approximate side chains"
        )
      ),
      m_ApproximateLoops
      (
        new command::FlagStatic
        (
          "apx_loops",
          "pass this flag to approximate loop regions"
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new Quality
    RestraintSaxs *RestraintSaxs::Clone() const
    {
      return new RestraintSaxs( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RestraintSaxs::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
    // but which should not be displayed, e.g. for help
    storage::Vector< std::string> RestraintSaxs::GetDeprecatedAppNames() const
    {
      return storage::Vector< std::string>( size_t( 1), "SimulateSaxsData");
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string RestraintSaxs::GetDescription() const
    {
      return "Create SAXS profiles from input pdb file and compare the profile generated with an experimental profile";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &RestraintSaxs::GetReadMe() const
    {
      static const std::string readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::SAXS, terms of use, appropriate citation, installation "
        "procedures, BCL::SAXS execution, technical support, and future research directions.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::SAXS?\n"
        "BCL::SAXS is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Fold.  BCL::SAXS computes SAXS profiles from an input pdb "
        "file and compares the profile generated with the experimental SAXS profiles.  There are three levels of "
        "approximation: complete models, approximated side chains, approximated side chains and loop regions \n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::SAXS.\n"
        "When using BCL::SAXS in a publication, please cite the following publication describing the application's "
        "development:  BCL::SAXS:GPU accelerated Debye method for computation of small angle X-ray scattering profiles.\n"
        "Putnam DK, Weiner BE, Woetzel N, Lowe ED Jr, Meiler J, Proteins. 2015 Aug 83(8):1500-12"
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL::SAXS.\n"
        "Consists of four steps :\n"
        "1) Generate Clean BCL Protein Model from PDB\n"
        " bcl.exe protein:PDBConvert [target].pdb -bcl_pdb -output_prefix [target]_"
        "2) Generate Solvent Accessible Surface Area file \n"
        " pdb_to_xyzr [target}_bcl.pdb > [target]_bcl.xyzr"
        " msms -if [target]_bcl.xyzr -af [target]_bcl.area -probe_radius 1.399"
        "3) Remove the header information on the experimental file"
        " The format should be: q_value intensity error"
        "4) Run restraint: SimulateSaxsData ( see example commandlines below). \n"
        "\n"
        "ADDITIONAL INFORMATION\n"
        "\n"
        "Information concerning syntax and flags can be obtained by typing : <bcl.exe> SimulateSaxsData -help\n"
        "\n"
        "For further explanation, examples of the flags, example command lines, input and output format information,\n"
        "and example output please see the documentation file.\n"
        "\n"
        "For more general information about the product, type : <bcl.exe> SimulateSaxsData -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::SAXS.\n"
        "BCL::SAXS is under ongoing further development.\n"
        "\n"
        + DefaultSectionSeparator() +
        "IX. EXAMPLE COMMANDLINE Complete Models"
        "<bcl.exe> -pdb input_protein.pdb "
        "-exp_data experimental_data.dat "
        "-aaclass AAComplete -sasa_data sasa_data.area -use_errors"
        "\n"
        "IX. EXAMPLE COMMANDLINE Approximated Side Chains"
        "<bcl.exe> -pdb input_protein.pdb "
         "-exp_data experimental_data.dat "
         "-aaclass AABackBone -apx_sc"
        "\n"
        "IX. EXAMPLE COMMANDLINE Approximated Side Chains and Loop regions"
        "<bcl.exe> -pdb input_protein.pdb - min_sse_size 5 3 999 "
         "-exp_data experimental_data.dat "
         "-aaclass AABackBone -apx_sc -apx_loops"

        "\n"
        + DefaultSectionSeparator()
      );
      return readme;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> RestraintSaxs::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add flag to give PDB File
      sp_cmd->AddFlag( m_PDBFilenameFlag);

      // add flag to give experimental restraint file
      sp_cmd->AddFlag( m_ExpRestraintFileFlag);
      sp_cmd->AddFlag( m_SasaFileFlag);

      // add output file flag
      sp_cmd->AddFlag( m_OutputFileFlag);

      // pdb factory flags
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());
      sp_cmd->AddFlag( pdb::Factory::GetFlagSSEsFromBackBone());
      sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());

      // add printing flags
      sp_cmd->AddFlag( m_PrintRMSD);
      sp_cmd->AddFlag( m_PrintStovgaard);
      sp_cmd->AddFlag( m_PrintProteinModelWithLoopsFileFlag);
      sp_cmd->AddFlag( m_PrintDerivative);

      // add flag for using Simulated/Experimental Errors and approximations
      sp_cmd->AddFlag( m_UseErrors);
      sp_cmd->AddFlag( m_ApproximateSideChains);
      sp_cmd->AddFlag( m_ApproximateLoops);

      // add flag for reducing data set using Shannon Sampling
      sp_cmd->AddFlag( m_ReduceData);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int RestraintSaxs::Main() const
    {

      // Read in Protein Model data file from application Parameter
      assemble::ProteinModel protein_model
      (
        pdb::Factory().ProteinModelFromPDBFilename( m_PDBFilenameFlag->GetFirstParameter()->GetValue())
      );

      // Read and store Unprocessed Experimental Data can be raw or reduced
      util::ShPtr< restraint::SasScatteringData> sp_experimental_data( new restraint::SasScatteringData());
      util::ShPtr< restraint::SasDensityData> sp_density_data( new restraint::SasDensityData());

      io::IFStream read;
      io::File::MustOpenIFStream( read, m_ExpRestraintFileFlag->GetFirstParameter()->GetValue());
      sp_experimental_data->ReadFromDataFile( read);
      io::File::CloseClearFStream( read);

      // Reduce the Experimental data to the Shannon Set
      if( m_ReduceData->GetFlag())
      {
        // Create Data object to hold reduced_data and populate the object with the reduced data set
        util::ShPtr< restraint::SasScatteringData> sp_reduced_experimental_data
        (
          restraint::SaxsDataReduction::SasSignalRecovery( protein_model, 1001, sp_experimental_data, 65)
        );

        sp_experimental_data = sp_reduced_experimental_data;
      }

      // Sasa stands for Solvent Accesible surface area and is used for hydration shell optimization.  This step is for
      // complete protein models to fine tune the fit to experimental data

      if( m_SasaFileFlag->GetFlag())
      {
        // Tell the protein model how to access SAXS Data and store the experimental data
        util::ShPtr< assemble::ProteinModelData> sp_protein_model_data( protein_model.GetProteinModelData());

        // Read in sasa data file from application Parameter
        // Sasa File must be preprocessed from PDBConvert and MSMS
        util::ShPtr< biol::SasaData> sp_sasa_data( new biol::SasaData());

        // establish input/output stream
        io::IFStream read_sasa;
        io::File::MustOpenIFStream( read_sasa, m_SasaFileFlag->GetFirstParameter()->GetValue());

        // Read Sasa Data
        sp_sasa_data->ReadFromDataFile( read_sasa);

        // Close Stream
        io::File::CloseClearFStream( read_sasa);

        //  Store the sasa data
        sp_protein_model_data->Insert( assemble::ProteinModelData::e_Sasa, sp_sasa_data);

        // Insert the Data into the protein model
        protein_model.SetProteinModelData( sp_protein_model_data);
      }

      util::ShPtr< restraint::SasExperimentalAndCalculatedData> sp_data_sets;

      float C1( 1.0);
      float C2( 0.0);

      // set default transformations
      storage::Vector< restraint::SasTransformation::TransformationTypeEnum> default_transformations;
      default_transformations.PushBack( restraint::SasTransformation::e_NormalizeData);
      default_transformations.PushBack( restraint::SasTransformation::e_Log10Profiles);
      default_transformations.PushBack( restraint::SasTransformation::e_ScaleCalculatedProfile);
      restraint::SasTransformation transform( default_transformations, false, m_UseErrors->GetFlag(), 1.0);

      // Fit the Experimental and Calculated scattering profile
      if( m_SasaFileFlag->GetFlag())
      {
        // use default parameters in SaxsOptimization object
        restraint::SasOptimization find_parameters;

        // set stage of find_parameters object
        find_parameters.SetScoreFunction( score::SasType::e_chi);
        find_parameters.SetErrorFlag( true);
        find_parameters.SetTransformationTypes( transform);
        find_parameters.SetExperimentalData( sp_experimental_data);
        find_parameters.SetApproximateSideChainsFlag( false);
        find_parameters.SetHardwareType( false);
        find_parameters.SetSasType( false);
        find_parameters.SetDeuteriumExchangeParameter( 0.0);

        restraint::SaxsOptiResult optimal_values( find_parameters( protein_model));

        C1 = optimal_values.GetC1();
        C2 = optimal_values.GetC2();

        BCL_MessageStd( "Optimal c1: " + util::Format()( C1));
        BCL_MessageStd( "Optimal c2: " + util::Format()( C2));
      }

      util::Implementation< restraint::SasDebyeInterface> saxs
      (
        restraint::SasAnalysis::SetDebyeImplementation
        (
          m_ApproximateLoops->GetFlag(),
          m_ApproximateSideChains->GetFlag(),
          C1,
          C2,
          false,
          false,
          0.0
        )
      );

      // Set the Experimental Data
      saxs->SetExperimentalData( sp_experimental_data);

      //Create container to hold experimental and computed SAXS profiles and store results in container
      util::ShPtr< restraint::SasExperimentalAndCalculatedData> sp_experimental_and_calculated_data
      (
        new restraint::SasExperimentalAndCalculatedData( saxs->operator()( protein_model))
      );

      // Transform the data
      restraint::SasExperimentalAndCalculatedData transformed_data( transform( *sp_experimental_and_calculated_data));

      // Score the transformed data
      double result
      (
        score::SasType( transform.GetUseErrors(), score::SasType::e_chi)( transformed_data)
      );

      BCL_MessageStd( "Score: " + util::Format()( result));

      // open stream to write the data
      io::OFStream write;
      io::File::MustOpenOFStream( write, m_OutputFileFlag->GetFirstParameter()->GetValue());
      transformed_data.WriteToGnuplot( write);

      // close the ofstream
      io::File::CloseClearFStream( write);

      // write out protein model if desired
      if( m_PrintProteinModelWithLoopsFileFlag->GetFlag())
      {
        pdb::Factory factory( biol::GetAAClasses().e_AAComplete);
        io::OFStream output;
        io::File::MustOpenOFStream( output, m_PrintProteinModelWithLoopsFileFlag->GetFirstParameter()->GetValue());

        if( saxs->GetLabel().GetArgument( 1).GetValue() == "0")
        {
          factory.WriteModelToPDB( fold::AddParabolicLoops( false)( protein_model), output);
        }
        else
        {
          factory.WriteModelToPDB( fold::AddParabolicLoops( true)( protein_model), output);
        }
        io::File::CloseClearFStream( output);
      }

      if( m_PrintDerivative->GetFlag())
      {
        io::OFStream output;
        io::File::MustOpenOFStream( output, m_PrintDerivative->GetFirstParameter()->GetValue());

        // Normalize Data and take log base 10
        restraint::SasExperimentalAndCalculatedData scaled_data( *sp_data_sets);

        // normalize the data to the maximum value in the experimental saxs curve
        scaled_data.SetYScale( 1.0);

        // log10 the data
        scaled_data.Log10();

        // Calculate Scaling Factor based on sigma
        const double scaling_factor( restraint::SasAnalysis::CalculateScalingWeight( scaled_data, false));

        // Scale the Calculated Data Intensity values by the scaling factor
        scaled_data.ScaleCalculatedData( scaling_factor);

        // Derivative Function
        scaled_data.Derivative();

        // write
        scaled_data.WriteToGnuplot( output);

        io::File::CloseClearFStream( output);
      }

      //successful end
      return 0;
    }

    //! @brief read from std::istream  find_parameters.SetHardwareType( false);
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RestraintSaxs::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &RestraintSaxs::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

    const ApplicationType RestraintSaxs::RestraintSaxs_Instance
    (
      GetAppGroups().AddAppToGroup( new RestraintSaxs(), GetAppGroups().e_Restraint)
    );

  } // namespace app
} // namespace bcl
