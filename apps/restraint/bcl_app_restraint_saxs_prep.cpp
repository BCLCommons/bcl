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
#include "bcl_app_restraint_saxs_prep.h"

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
#include "restraint/bcl_restraint_sas_experimental_and_calculated_data.h"
#include "restraint/bcl_restraint_sas_optimization.h"
#include "restraint/bcl_restraint_saxs_data_reduction.h"
#include "restraint/bcl_restraint_saxs_opti_result.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RestraintSaxsPrep::RestraintSaxsPrep() :
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
            ""
          )
        )
      ),
      m_PDBFilenameFlag
      (
        new command::FlagStatic
        (
          "pdb_file",
          "\tpdb file for which Shannon Sampling will be evaluated on",
          command::Parameter
          (
            "pdb_filename",
            "\tfilename for which Shannon Sampling will be evaluated on",
            ""
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
            "reduced_saxs.data"
          )
        )
      ),
      m_OutputProteinModelFlag
      (
        new command::FlagStatic
        (
          "output_model",
          "path and name of the output file which will hold protein model with approximated loops.",
          command::Parameter
          (
            "string",
            "\tpath and name of the outputfile",
            "approximated_model.pdb"
          )
        )
      ),
      m_UseAnalyticNormFactorFlag
      (
        new command::FlagStatic
         (
           "use_analytic",
           "pass this flag to use regula falsi for normalization factor, default is triangualar approximation"
         )
      ),
      m_DmaxFlag
      (
        new command::FlagStatic
        (
          "dmax", "This flag is the maximum dimension of the protein",
          command::Parameter
          (
            "maximum_particle_dimension",
            "The value of the largest distance in the protein - default undefined",
            util::Format()( util::GetUndefined< double>())
          )
        )
      ),
      m_GnomeFlag
      (
        new command::FlagStatic
        (
           "gnome_format",
           "pass this flag to output data in format gnome uses for input"
        )
      ),
      m_SamplingRoundsFlag
      (
        new command::FlagStatic
        (
          "sampling_rounds", "This flag is the number of iterations for shannon sampling",
          command::Parameter
          (
            "k",
            "the number of iterations for shannon sampling - default 1001",
            util::Format()( 1001)
          )
        )
      ),
      m_SimulateErrorFlag
      (
        new command::FlagStatic
        (
          "simulate_errors",
          "pass this flag if you want to simulate experimental errors on the dataset"
        )
      ),
      m_ScoreFileFlag
      (
        new command::FlagStatic
        (
          "score_file",
          "\tfile with experimental saxs profiles and calculated saxs profiles formatted through bcl ",
          command::Parameter
          (
            "score_filename",
            "\tfilename with bcl formated data",
            ""
          )
        )
      ),
      m_ReadGnomeFitFlag
      (
        new command::FlagStatic
        (
          "gnome_fit",
          "pass this flag if you want to read the transformed p(r) curve as the experimental data representation"
        )
      ),
      m_SetErrorFlag
      (
        new command::FlagStatic
        (
          "set_error", "This is a constant value to set error values to",
          command::Parameter
          (
            "set_error_values",
            "Constant error value - default undefined",
            util::Format()( util::GetUndefined< double>())
          )
        )
      ),
      m_ReduceDataFlag
      (
        new command::FlagStatic
        (
          "reduce_data",
          "pass this flag if you want to use Shannon Sampling to reduce experimental dataset"
        )
      ),
      m_ReduceDataMinErrorFlag
      (
        new command::FlagStatic
        (
          "reduce_data_min_error",
          "pass this flag if you want to use Min Error to reduce experimental dataset"
        )
      ),
      m_UseErrorFlag
      (
        new command::FlagStatic
        (
          "use_error",
          "pass this flag if you want to use error in the scorefile command"
        )
      ),
      m_BCLFileNameFlag
      (
        new command::FlagStatic
        (
          "bcl_file",
          "\tfile with experimental saxs profiles and calculated saxs profiles through bcl ",
          command::Parameter
          (
            "bcl_filename",
            "\tfilename with bcl computed data",
            ""
          )
        )
      ),
      m_FoxsFileNameFlag
      (
        new command::FlagStatic
        (
          "foxs_file",
          "\tfile with experimental saxs profiles and calculated saxs profiles through foxs ",
          command::Parameter
          (
            "foxs_filename",
            "\tfilename with foxs computed data",
            ""
          )
        )
      ),
      m_CrysolFileNameFlag
      (
        new command::FlagStatic
        (
          "crysol_file",
          "\tfile with experimental saxs profiles and calculated saxs profiles through crysol ",
          command::Parameter
          (
            "crysol_filename",
            "\tfilename with crysol computed data",
            ""
          )
        )
      ),
      m_LogScaleFlag
      (
        new command::FlagStatic
        (
          "log_scale",
          "pass this flag if you want to convert from absolute to log scale"
        )
      ),
      m_DerivativeFlag
      (
        new command::FlagStatic
        (
          "derivative",
          "pass this flag if you want to use derivative of input"
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new Quality
    RestraintSaxsPrep *RestraintSaxsPrep::Clone() const
    {
      return new RestraintSaxsPrep( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RestraintSaxsPrep::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
    // but which should not be displayed, e.g. for help
    storage::Vector< std::string> RestraintSaxsPrep::GetDeprecatedAppNames() const
    {
      return storage::Vector< std::string>( size_t( 1), "SimulateSaxsDataPrep");
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string RestraintSaxsPrep::GetDescription() const
    {
      return "Add Simulated error and produce clean saxs signal via shannon sampling";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &RestraintSaxsPrep::GetReadMe() const
    {
      static const std::string readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::SAXSPREP, terms of use, appropriate citation, installation "
        "procedures, BCL::SAXS execution, technical support, and future research directions.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::SAXSPREP?\n"
        "BCL::SAXSPREP is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Fold.  BCL::SAXSPREP simulates experimental error of computed"
        "SAXS profiles from an input PDB file and / or reduces the number of experimental data points using Shannon"
        "sampling \n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
//        + DefaultSectionSeparator() +
//        "IV. APPROPRIATE CITATIONS FOR USING BCL::SAXS.\n"
//        "When using BCL::SAXS in a publication, please cite the following publication describing the application's "
//        "development:\n"
//        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL::SAXSPREP.\n"
        "Consists of four steps :\n"
        "1) Gather your experimental SAXS data either simulated i.e crysol or provided through experiment. \n"
        "2) Prepare Solvent Accessible Surface Area file ( using methods such as MSMS). \n"
        "3) Run Application \n"
        "\n"
        "ADDITIONAL INFORMATION\n"
        "\n"
        "Information concerning syntax and flags can be obtained by typing : <bcl.exe> Application -help\n"
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
        "IX. EXAMPLE COMMANDLINE Read Fit from Gnome and set experimental error to constant for folding"
        "<bcl.exe> restraint:SaxsPrep -exp_data input.data"
        "-output_file output.data -gnome_fit -set_error 1.0"
        "\n"
        "IX. EXAMPLE COMMANDLINE Reduce Experimental data using Shannon Sampling"
        "<bcl.exe> SimulateSaxsData 'SasDebye(consider loops=0, analytic=0)' "
        "<bcl.exe> restraint:SaxsPrep -exp_data input.data"
         "-output_file output.data -pdb_file input.pdb -dmax 65 -sampling_rounds 1001 -simulate_errors -reduce_data"
        "\n"
        + DefaultSectionSeparator()
      );
      return readme;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> RestraintSaxsPrep::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // Add flags to command object
      sp_cmd->AddFlag( m_ExpRestraintFileFlag);
      sp_cmd->AddFlag( m_PDBFilenameFlag);
      sp_cmd->AddFlag( m_DmaxFlag);
      sp_cmd->AddFlag( m_GnomeFlag);
      sp_cmd->AddFlag( m_SamplingRoundsFlag);
      sp_cmd->AddFlag( m_OutputFileFlag);
      sp_cmd->AddFlag( m_SimulateErrorFlag);
      sp_cmd->AddFlag( m_SetErrorFlag);
      sp_cmd->AddFlag( m_UseErrorFlag);
      sp_cmd->AddFlag( m_ReadGnomeFitFlag);
      sp_cmd->AddFlag( m_ReduceDataFlag);
      sp_cmd->AddFlag( m_ReduceDataMinErrorFlag);
      sp_cmd->AddFlag( m_OutputProteinModelFlag);
      sp_cmd->AddFlag( m_UseAnalyticNormFactorFlag);
      sp_cmd->AddFlag( m_BCLFileNameFlag);
      sp_cmd->AddFlag( m_FoxsFileNameFlag);
      sp_cmd->AddFlag( m_CrysolFileNameFlag);
      sp_cmd->AddFlag( m_LogScaleFlag);
      sp_cmd->AddFlag( m_DerivativeFlag);
      sp_cmd->AddFlag( m_ScoreFileFlag);

      // Add All Flags from pdb::Factory
      for
      (
        util::ShPtrVector< command::FlagInterface>::const_iterator
         itr( pdb::Factory::GetAllFlags().Begin()), itr_end( pdb::Factory::GetAllFlags().End());
        itr != itr_end;
        ++itr
      )
      {
        sp_cmd->AddFlag( *itr);
      }

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int RestraintSaxsPrep::Main() const
    {
      if( m_ScoreFileFlag->GetFlag())
      {
        // Create object to hold data
        util::ShPtr< restraint::SasExperimentalAndCalculatedData> sp_Data( new restraint::SasExperimentalAndCalculatedData());

        // Read the BCL Data
        io::IFStream read;
        io::File::MustOpenIFStream( read, m_ScoreFileFlag->GetFirstParameter()->GetValue());
        sp_Data->ReadFromDataFile( read);
        io::File::CloseClearFStream( read);

        double score( 0.0);

        if( m_UseErrorFlag->GetFlag())
        {
          score = score::SasType( true, score::SasType::e_chi)( *sp_Data);
        }
        else
        {
          score = score::SasType( false, score::SasType::e_chi)( *sp_Data);
        }

        std::string error;

        ( m_UseErrorFlag->GetFlag()) ? ( error = "true") : ( error = "false");

        std::string result
        (
          "\n ScoreFunction: Chi \n"
          " Use errors: " + error + "\n"
          " Score: " + util::Format()( score)
        );

        BCL_MessageStd( result);

        if( m_DerivativeFlag->GetFlag())
        {

          restraint::SasExperimentalAndCalculatedData original( *sp_Data);
          original.Log10();
          original.Derivative();
          original.WriteToGnuplotFileName( "scored_derivative.data");
        }

        return 0;
      }

      if(
          m_BCLFileNameFlag->GetFlag() ||
          m_CrysolFileNameFlag->GetFlag() ||
          m_FoxsFileNameFlag->GetFlag() ||
          m_ExpRestraintFileFlag->GetFlag()
        )
      {
        // Create Objects to hold data
        util::ShPtr< restraint::SasExperimentalAndCalculatedData> sp_BCLData( new restraint::SasExperimentalAndCalculatedData());
        util::ShPtr< restraint::SasExperimentalAndCalculatedData> sp_CrysolData( new restraint::SasExperimentalAndCalculatedData());
        util::ShPtr< restraint::SasExperimentalAndCalculatedData> sp_FoxsData( new restraint::SasExperimentalAndCalculatedData());

        util::ShPtr< restraint::SasScatteringData> sp_ExperimentalData( new restraint::SasScatteringData);

        double bcl_score( util::GetUndefined< double>());
        double crysol_score( util::GetUndefined< double>());
        double foxs_score( util::GetUndefined< double>());

        if( m_ExpRestraintFileFlag->GetFlag())
        {
          io::IFStream read;
          io::File::MustOpenIFStream( read, m_ExpRestraintFileFlag->GetFirstParameter()->GetValue());
          sp_ExperimentalData->ReadFromDataFile( read);
          io::File::CloseClearFStream( read);

          if( m_LogScaleFlag->GetFlag())
          {
            BCL_MessageStd( "Inside Code covert to log scale:");
            restraint::SasScatteringData log10ExperimentalData( restraint::SasAnalysis::Log10( *sp_ExperimentalData));

            BCL_MessageStd( "write the data:");
            log10ExperimentalData.WriteToDataFileName( "log10_experimental.data", false);
          }
        }

        if( m_BCLFileNameFlag->GetFlag())
        {
          // Read the BCL Data
          io::IFStream read;
          io::File::MustOpenIFStream( read, m_BCLFileNameFlag->GetFirstParameter()->GetValue());
          sp_BCLData->ReadFromDataFile( read);
          io::File::CloseClearFStream( read);

          // Write the Experimental Data to an object
          util::ShPtr< restraint::SasScatteringData> tmp_ptr( sp_BCLData->GetExperimentalData().Clone());
          sp_ExperimentalData = tmp_ptr;

          // Transform the data to log scale if specified
          if( m_LogScaleFlag->GetFlag())
          {
            sp_BCLData->Log10();
            sp_BCLData->WriteToGnuplotFileName( "log10_bcl.data");
            if( m_GnomeFlag->GetFlag())
            {
              sp_BCLData->WriteToGnomeFileName( "log10_bcl_gnom.data");
            }

            // Score the transformed data
            bcl_score = score::SasType( true, score::SasType::e_chi)( *sp_BCLData);
          }
          else
          {
            sp_BCLData->WriteToGnuplotFileName( "absolute_bcl.data");
            if( m_GnomeFlag->GetFlag())
            {
              sp_BCLData->WriteToGnomeFileName( "absolute_bcl_gnom.data");
              sp_ExperimentalData->WriteToDataFileName( "absolute_experimental_gnom.data", false);
            }

            bcl_score = score::SasType( true, score::SasType::e_chi)( *sp_BCLData);
          }
        }

        if( m_CrysolFileNameFlag->GetFlag())
        {
          // Read the Crysol Data
          io::IFStream read;
          io::File::MustOpenIFStream( read, m_CrysolFileNameFlag->GetFirstParameter()->GetValue());
          sp_CrysolData->ReadFromCrysolFitFile( read, sp_ExperimentalData->Begin()->GetQvalue());
          io::File::CloseClearFStream( read);

          // insert experimental errors into crysol data
          sp_CrysolData->SetExperimentalData( *sp_ExperimentalData);

          if( m_LogScaleFlag->GetFlag())
          {
            sp_CrysolData->Log10();
            sp_CrysolData->WriteToGnuplotFileName( "log10_crysol.data");
            if( m_GnomeFlag->GetFlag())
            {
              sp_CrysolData->WriteToGnomeFileName( "log10_crysol_gnom.data");
            }

            // Score the transformed data
            crysol_score = score::SasType( true, score::SasType::e_chi)( *sp_CrysolData);
          }
          else
          {
            sp_CrysolData->WriteToGnuplotFileName( "absolute_crysol.data");
            if( m_GnomeFlag->GetFlag())
            {
              sp_CrysolData->WriteToGnomeFileName( "absolute_crysol_gnom.data");
            }

            crysol_score = score::SasType( true, score::SasType::e_chi)( *sp_CrysolData);
          }
        }

        if( m_FoxsFileNameFlag->GetFlag())
        {
          // Read the Foxs Data
          io::IFStream read;
          io::File::MustOpenIFStream( read, m_FoxsFileNameFlag->GetFirstParameter()->GetValue());
          sp_FoxsData->ReadFromFoxsFile( read);
          io::File::CloseClearFStream( read);

          // insert experimental errors into foxs data
          sp_FoxsData->SetExperimentalData( *sp_ExperimentalData);

          if( m_LogScaleFlag->GetFlag())
          {
            sp_FoxsData->Log10();
            sp_FoxsData->WriteToGnuplotFileName( "log10_foxs.data");
            if( m_GnomeFlag->GetFlag())
            {
              sp_FoxsData->WriteToGnomeFileName( "log10_foxs_gnom.data");
            }

            // Score the transformed data
            foxs_score = score::SasType( true, score::SasType::e_chi)( *sp_FoxsData);
          }
          else
          {
            sp_FoxsData->WriteToGnuplotFileName( "absolute_foxs.data");
            if( m_GnomeFlag->GetFlag())
            {
              sp_FoxsData->WriteToGnomeFileName( "absolute_foxs_gnom.data");
            }

            foxs_score = score::SasType( true, score::SasType::e_chi)( *sp_FoxsData);
          }
        }

        std::string state;

        ( m_LogScaleFlag->GetFlag()) ? ( state = "log10") : ( state = "absolute");

        std::string summary
        (
          "\n ScoreFunction: Chi \n"
          " Use errors: true \n"
          " State: " + state + "\n"
          " BCL Score: " + util::Format()( bcl_score) + "\n"
          " Crysol Score: " + util::Format()( crysol_score) + "\n"
          " Foxs Score: " + util::Format()( foxs_score) + "\n"
        );

        BCL_MessageStd( summary);
      }

      if( m_ExpRestraintFileFlag->GetFlag())
      {
        // Read and store Experimental Data
        // First create Shared pointers to the appropriate data type
        util::ShPtr< restraint::SasScatteringData> sp_experimental_data( new restraint::SasScatteringData());

        // Read transformation of P(r) function from Gnome.  In this case, experimental error is set by commandline
        // If experimental error is not set at the commandline, it will be set to nan.
        if( m_ReadGnomeFitFlag->GetFlag())
        {
          io::IFStream read;
          io::File::MustOpenIFStream( read, m_ExpRestraintFileFlag->GetFirstParameter()->GetValue());
          sp_experimental_data->ReadFitFromGnom( read);
          io::File::CloseClearFStream( read);

          *sp_experimental_data =
            restraint::SasAnalysis::SetErrors
            (
              *sp_experimental_data,
              m_SetErrorFlag->GetFirstParameter()->GetNumericalValue< double>()
            );
        }
        else
        {
          io::IFStream read;
          io::File::MustOpenIFStream( read, m_ExpRestraintFileFlag->GetFirstParameter()->GetValue());
          sp_experimental_data->ReadFromDataFile( read);
          io::File::CloseClearFStream( read);
        }

        // Simulate error if errors are not present
        if( m_SimulateErrorFlag->GetFlag())
        {
          // Generate Errors into new SasScatteringData Object
          restraint::SasScatteringData simulated_experimental_error
          (
            restraint::SasAnalysis::AddErrors( *sp_experimental_data)
          );

          // Clone the Created object into a shared pointer
          util::ShPtr< restraint::SasScatteringData> sp_simulated_experimental_error
          (
            simulated_experimental_error.Clone()
          );

          // Point the experimental data to the new object
          sp_experimental_data = sp_simulated_experimental_error;
        }

        // Reduce the Experimental data using min error
        if( m_ReduceDataMinErrorFlag->GetFlag())
        {
          // assert that the errors are defined
          BCL_Assert
          (
            sp_experimental_data->IsErrorDefined(),
            "Dataset has undefined error values, you must simulate errors before data reduction can happen"
          );

          sp_experimental_data =
            restraint::SaxsDataReduction::SasSignalRecoveryEstimate
            (
              sp_experimental_data,
              m_DmaxFlag->GetFirstParameter()->GetNumericalValue< double>()
            );
        }

        // Reduce the Experimental data using Shannon Sampling
        if( m_ReduceDataFlag->GetFlag())
        {

          // assert that the errors are defined
          BCL_Assert
          (
            sp_experimental_data->IsErrorDefined(),
            "Dataset has undefined error values, you must simulate errors before data reduction can happen"
          );

          // Read in Protein Model specified in Command line
          assemble::ProteinModel protein_model
          (
            pdb::Factory( biol::GetAAClasses().e_AAComplete).ProteinModelFromPDBFilename
            (
              m_PDBFilenameFlag->GetFirstParameter()->GetValue()
            )
          );

          BCL_MessageDbg( "Protein Model read in ");

          sp_experimental_data =
            restraint::SaxsDataReduction::SasSignalRecovery
            (
              protein_model,
              m_SamplingRoundsFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
              sp_experimental_data,
              m_DmaxFlag->GetFirstParameter()->GetNumericalValue< double>()
            );
        }

        if( m_SetErrorFlag->GetFlag())
        {
          *sp_experimental_data =
            restraint::SasAnalysis::SetErrors
            (
              *sp_experimental_data,
              m_SetErrorFlag->GetFirstParameter()->GetNumericalValue< double>()
            );
        }

        // open stream to write the data
        io::OFStream write;
        io::File::MustOpenOFStream( write, m_OutputFileFlag->GetFirstParameter()->GetValue());
        sp_experimental_data->WriteToDataFile( write);

        // close the ofstream
        io::File::CloseClearFStream( write);
      }

      if( m_OutputProteinModelFlag->GetFlag())
      {
        // Read in Protein Model data file from application Parameter
        assemble::ProteinModel protein_model
        (
          pdb::Factory( biol::GetAAClasses().e_AABackBone).ProteinModelFromPDBFilename
          (
            m_PDBFilenameFlag->GetFirstParameter()->GetValue()
          )
        );

        pdb::Factory factory( biol::GetAAClasses().e_AAComplete);
        io::OFStream output;
        io::File::MustOpenOFStream( output, m_OutputProteinModelFlag->GetFirstParameter()->GetValue());

        if( m_UseAnalyticNormFactorFlag->GetFlag())
        {
          // Case to use Regula Falsi Solution
          factory.WriteModelToPDB( fold::AddParabolicLoops( true)( protein_model), output);
        }
        else
        {
          // Default case to use Triangular approximation
          factory.WriteModelToPDB( fold::AddParabolicLoops( false)( protein_model), output);
        }

        io::File::CloseClearFStream( output);
      }

      //successful end
      return 0;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RestraintSaxsPrep::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &RestraintSaxsPrep::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

    const ApplicationType RestraintSaxsPrep::RestraintSaxsPrep_Instance
    (
      GetAppGroups().AddAppToGroup( new RestraintSaxsPrep(), GetAppGroups().e_Restraint)
    );

  } // namespace app
} // namespace bcl
