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
#include "bcl_app_restraint_pofr.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "fold/bcl_fold_add_parabolic_loops.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "restraint/bcl_restraint_sas_analysis.h"
#include "restraint/bcl_restraint_sas_experimental_and_calculated_density.h"
#include "restraint/bcl_restraint_sas_pofr.h"
#include "score/bcl_score_pofr.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RestraintPofr::RestraintPofr() :
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
          "\tmodel to be evaluated",
          command::Parameter
          (
            "pdb_filename",
            "\tfilename of model to be evaluated",
            ""
          )
        )
      ),
      m_ApproximateLoops
      (
        new command::FlagStatic
        (
          "approximate_loops",
          "pass this flag if you want to approximate loops prior to scoring"
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new Quality
    RestraintPofr *RestraintPofr::Clone() const
    {
      return new RestraintPofr( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RestraintPofr::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
    // but which should not be displayed, e.g. for help
    storage::Vector< std::string> RestraintPofr::GetDeprecatedAppNames() const
    {
      return storage::Vector< std::string>( size_t( 1), "RestraintPofr");
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string RestraintPofr::GetDescription() const
    {
      return "Compute pairwise distance histogram of bcl_model and scale to experimental data";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &RestraintPofr::GetReadMe() const
    {
      static const std::string readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::POFR, terms of use, appropriate citation, installation "
        "procedures, BCL::POFR execution, technical support, and future research directions.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL::POFR?\n"
        "BCL::POFR is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
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
    util::ShPtr< command::Command> RestraintPofr::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // Add flags to command object
      sp_cmd->AddFlag( m_ExpRestraintFileFlag);
      sp_cmd->AddFlag( m_PDBFilenameFlag);
      sp_cmd->AddFlag( m_ApproximateLoops);

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
    int RestraintPofr::Main() const
    {
      // Create object to hold data
      util::Implementation< restraint::SasDensityData> sp_data( new restraint::SasDensityData());

      // Read Experimental Data
      io::IFStream read;
      io::File::MustOpenIFStream( read, m_ExpRestraintFileFlag->GetFirstParameter()->GetValue());
      sp_data->ReadFromDataFile( read);
      io::File::CloseClearFStream( read);

      // Read BCL Model from application Parameter
      assemble::ProteinModel protein_model
      (
        pdb::Factory( biol::GetAAClasses().e_AAComplete).ProteinModelFromPDBFilename
        (
          m_PDBFilenameFlag->GetFirstParameter()->GetValue()
        )
      );

      // Create a pointer to point either at a protein model or a protein model with approximated loop regions
      util::ShPtr< assemble::ProteinModel> sp_protein_model( protein_model.Clone());

      if( m_ApproximateLoops->GetFlag())
      {
        // Write BCL Model with approximated loops
        pdb::Factory factory( biol::GetAAClasses().e_AAComplete);
        io::OFStream output;
        io::File::MustOpenOFStream( output, "apx_loops.pdb");

        // Store the model with the approximated loops using triangular approximation
        assemble::ProteinModel loop_apx_model( fold::AddParabolicLoops( false)( protein_model));

        // Write the approximated model to a pdb file
        factory.WriteModelToPDB( loop_apx_model, output);

        io::File::CloseClearFStream( output);

        // point the pointer to the approximated model
        sp_protein_model = util::ShPtr< assemble::ProteinModel>( loop_apx_model.Clone());
      }

      // Add the restraints to the restraint::SasPofR object for storage
      restraint::SasPofR density_result;
      density_result.SetExperimentalDensity( sp_data);

      // Compute the p of r profile
      restraint::SasExperimentalAndCalculatedDensity data_sets( density_result( *sp_protein_model));

      double auc_difference( restraint::SasAnalysis::CalculatePofRIntegralScore( data_sets));

      double cal_dmax( data_sets.GetCalculatedDensity().GetDmax());
      double exp_dmax( data_sets.GetExperimentalDensity().GetDmax());

      double excess_integral( restraint::SasAnalysis::CalculatePofRExcessIntegralScore( data_sets));
      double cal_oscillation( restraint::SasAnalysis::CalculatePofROscillationScore( data_sets.GetCalculatedDensity()));
      double exp_oscillation( restraint::SasAnalysis::CalculatePofROscillationScore( data_sets.GetExperimentalDensity()));

      double oscillation_score( math::Absolute( exp_oscillation - cal_oscillation));
      double dmax_difference( exp_dmax - cal_dmax);

      if( excess_integral > 0.1)
      {
        excess_integral = 10 * excess_integral;
      }

      double total_score( ( dmax_difference) + ( auc_difference) + ( 10 * excess_integral) + 10 * oscillation_score);
//      double total_score( dmax_difference + auc_difference);

      std::string summary
      (
        "\n Experimental Dmax: " + util::Format()( exp_dmax) + "\n"
        " Dmax difference: " + util::Format()( dmax_difference) + "\n"
        " Integral difference: " + util::Format()( auc_difference) + "\n"
        " Excess: " + util::Format()( excess_integral) + "\n"
        " Oscillation score: " + util::Format()( oscillation_score) + "\n"
        " Total Score: " + util::Format()( total_score) + "\n"
      );

      BCL_MessageStd( summary);

      // Print the normalized data for visualization
      restraint::SasExperimentalAndCalculatedDensity normalized_data( data_sets);

      normalized_data.ScaleCalculatedDensity( 1 / data_sets.GetCalculatedDensity().GetHmax());
      normalized_data.ScaleExperimentalDensity( 1 / data_sets.GetExperimentalDensity().GetHmax());

      normalized_data.WriteToGnuplotFileName( "pofr.data");

      normalized_data.ShiftDensity();

      normalized_data.WriteToGnuplotFileName( "shift.data");

      //successful end
      return 0;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RestraintPofr::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &RestraintPofr::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

    const ApplicationType RestraintPofr::RestraintPofr_Instance
    (
      GetAppGroups().AddAppToGroup( new RestraintPofr(), GetAppGroups().e_Restraint)
    );

  } // namespace app
} // namespace bcl
