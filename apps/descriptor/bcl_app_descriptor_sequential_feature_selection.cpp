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
#include "bcl_app_descriptor_sequential_feature_selection.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "model/bcl_model_descriptor_selection_interface.h"
#include "model/bcl_model_meta_data_storage_interface.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {

  /////////////////
  // data access //
  /////////////////

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> DescriptorSequentialFeatureSelection::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // filename for feature code object
      sp_cmd->AddParameter( m_EntireDescriptorCodeObjectFileName);
      // filename for assembled feature code object
      sp_cmd->AddParameter( m_FinalCodeObjectOutputFilename);
      // flag for round number
      sp_cmd->AddFlag( m_FlagRoundNumber);
      // flag for writing out initital descriptors for this iteration
      sp_cmd->AddFlag( m_FlagWriteInitialCodeObjectFile);
      // flag for setting the descriptor selection type
      sp_cmd->AddFlag( m_FlagDescriptorSelectionType);
      // flag for model::MetaData storage used in descriptor selection
      sp_cmd->AddFlag( m_FlagMetaDataStorage);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
    // but which should not be displayed, e.g. for help
    storage::Vector< std::string> DescriptorSequentialFeatureSelection::GetDeprecatedAppNames() const
    {
      return storage::Vector< std::string>( size_t( 1), "TrainModelDescriptorSelection");
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string DescriptorSequentialFeatureSelection::GetDescription() const
    {
      return "takes a descriptor groups file and generates new groups of features to test,"
             " based on the best descriptor set in a prior round";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &DescriptorSequentialFeatureSelection::GetReadMe() const
    {
      static const std::string readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL::DescriptorSequentialFeatureSelection, terms of use, appropriate citation, "
        "installation procedures, BCL::DescriptorSequentialFeatureSelection execution, technical support, and future research directions.\n"
        "\n"
        + Interface::DefaultSectionSeparator() +
        "II. WHAT IS BCL::DescriptorSequentialFeatureSelection?\n"
        "BCL::DescriptorSequentialFeatureSelection is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Commons. BCL::DescriptorSequentialFeatureSelection is an application\n"
        "that takes a descriptor groups file and generates individual files with sub groups. BCL::DescriptorSequentialFeatureSelection\n"
        "is part of a software suite termed BCL::ChemInfo and is essential for descriptor selection.\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + Interface::DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL::DescriptorSequentialFeatureSelection.\n"
        "When using BCL::DescriptorSequentialFeatureSelection in a publication, please cite the following publication describing the application's "
        "development:\n"
        "\n"
        "Butkiewicz, M.; Lowe, E.W., Jr.; Mueller, R.; Mendenhall, J.L.; Teixeira, P.L.; Weaver, C.D.; Meiler, J.\n"
        "'Benchmarking Ligand-Based Virtual High-Throughput Screening with the PubChem Database'. Molecules 2013, 18, 735-756.\n"
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. Information concerning syntax and flags can be obtained by typing <bcl.exe> DescriptorSequentialFeatureSelection -help -help\n"
        "\n"
        "For further explanation, examples of the flags, example command lines, input and output format information,\n"
        "and example output please see the documentation file.\n"
        "\n"
        "For more general information about the product, type <bcl.exe> DescriptorSequentialFeatureSelection -help\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL::DescriptorSequentialFeatureSelection.\n"
        "BCL::DescriptorSequentialFeatureSelection is under ongoing further development.\n"
        "\n"
        + DefaultSectionSeparator() +
        "IX. EXAMPLE COMMANDLINE"
        "<bcl.exe>  DescriptorSequentialFeatureSelection <entire_descriptor_code_object_filename> <final_code_output_filename>"
        "[flags]"
        + DefaultSectionSeparator()
      );

      return readme;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int DescriptorSequentialFeatureSelection::Main() const
    {
      // output stream
      io::OFStream output;

      // given round number
      const size_t round_number( m_FlagRoundNumber->GetFirstParameter()->GetNumericalValue< size_t>());

      // read in descriptor set(s)
      util::ObjectDataLabel entire_descriptor_set;

      // initialize input stream
      io::IFStream read;
      io::File::MustOpenIFStream( read, m_EntireDescriptorCodeObjectFileName->GetValue());
      entire_descriptor_set = util::ObjectDataLabel( read);
      io::File::CloseClearFStream( read);

      // set up the initial descriptors, give it the same name and label as the entire descriptor set
      util::ObjectDataLabel initial_descriptors( entire_descriptor_set.GetName(), entire_descriptor_set.GetValue());

      BCL_MessageStd( "Loading initial descriptors for this round");
      BCL_MessageStd
      (
        "Entire descriptor set contains " + util::Format()( entire_descriptor_set.GetNumberArguments())
        + " feature labels"
      );

      // model interface storage
      util::Implementation< model::MetaDataStorageInterface> model_storage
      (
        util::ObjectDataLabel( m_FlagMetaDataStorage->GetFirstParameter()->GetValue())
      );

      // initialize descriptor selection scheme
      util::Implementation< model::DescriptorSelectionInterface> descriptor_selection
      (
        util::ObjectDataLabel( m_FlagDescriptorSelectionType->GetFirstParameter()->GetValue())
      );

      // initialize descriptor selection scheme
      if( round_number > entire_descriptor_set.GetNumberArguments())
      {
        // then there can be no initial descriptors, the round number is too high
        BCL_MessageStd( "there are no initial descriptors, the round number is too high!");
        initial_descriptors = descriptor_selection->GetInitialDescriptorSet( entire_descriptor_set).ToStringDefaultWidth();
      }
      else
      {
        // get best descriptor set so far in storage
        initial_descriptors = model_storage->RetrieveBestDescriptorSetByResult();
      }

      BCL_MessageStd
      (
        "Found best descriptors for this round: " + initial_descriptors.ToStringDefaultWidth()
      );

      const std::string filename_base
      (
        round_number > 0
        ? m_FinalCodeObjectOutputFilename->GetValue() + "." + util::Format()( round_number - 1) + "."
        : m_FinalCodeObjectOutputFilename->GetValue() + ".final."
      );

      if( m_FlagWriteInitialCodeObjectFile->GetFlag())
      {
        io::File::MustOpenOFStream( output, filename_base + "initial");

        // write the initial descriptor set too
        output << initial_descriptors;

        io::File::CloseClearFStream( output);
      }

      // stop if round number is zero
      if( round_number == 0)
      {
        return 0;
      }

      // create vector of descriptor combinations
      storage::Vector< util::ObjectDataLabel> descriptor_ensemble_list
      (
        descriptor_selection->operator ()( initial_descriptors, entire_descriptor_set)
      );

      BCL_MessageDbg
      (
        "entire list of QSAR code: \n" + util::Format()( entire_descriptor_set)
        + "\ninitial list of QSAR code: \n" + util::Format()( initial_descriptors)
        + "\nlist of QSAR code: \n" + util::Format()( descriptor_ensemble_list)
      );

      // write out all possible ways to add/or eliminate one more descriptor to/from the current set
      for
      (
        size_t iteration_number( 0), final_itr( descriptor_ensemble_list.GetSize());
        iteration_number < final_itr;
        ++iteration_number
      )
      {
        io::File::MustOpenOFStream( output, filename_base + util::Format()( iteration_number));
        output << descriptor_ensemble_list( iteration_number).ToStringDefaultWidth();
        io::File::CloseClearFStream( output);
      }

      // end
      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief standard constructor
    DescriptorSequentialFeatureSelection::DescriptorSequentialFeatureSelection() :
      m_EntireDescriptorCodeObjectFileName
      (
        new command::Parameter
        (
          "entire_descriptor_code_object_filename", "\tcode object filename with all descriptors available"
        )
      ),
      m_FinalCodeObjectOutputFilename
      (
        new command::Parameter( "final_code_output_filename", "\tfilename of final code object file")
      ),
      m_FlagRoundNumber
      (
        new command::FlagStatic
        (
          "flag_round_number",
          "flag for round number!",
          command::Parameter
          (
            "flag_round_number",
            "flag for round number!",
            "0"
          )
        )
      ),
      m_FlagWriteInitialCodeObjectFile
      (
        new command::FlagStatic
        (
          "get_initial_descriptor_set_for_this_round",
          "flag for writing out the initial descriptor set object for this round, base descriptor set used for adding"
          " descriptors! output: {filename}.{round_number}.initial"
        )
      ),
      m_FlagDescriptorSelectionType
      (
        new command::FlagStatic
        (
          "descriptor_selection_type",
          "flag for setting the descriptor selection type label",
          command::Parameter
          (
            "descriptor_selection_type",
            "flag for setting the descriptor selection type label",
            command::ParameterCheckSerializable( util::Implementation< model::DescriptorSelectionInterface>()),
            "FeatureForward"
          )
        )
      ),
      m_FlagMetaDataStorage
      (
        new command::FlagStatic
        (
          "storage_descriptor_selection",
          "choice of storage for meta data (used in descriptor selection)",
          command::Parameter
          (
            "storage_descriptor_selection",
            "choice of storage for meta data (used in descriptor selection)",
            command::ParameterCheckSerializable( util::Implementation< model::MetaDataStorageInterface>()),
            "File"
          )
        )
      )
    {
    }

    const ApplicationType DescriptorSequentialFeatureSelection::DescriptorSequentialFeatureSelection_Instance
    (
      GetAppGroups().AddAppToGroup( new DescriptorSequentialFeatureSelection(), GetAppGroups().e_Descriptor)
    );

  } // namespace app
} // namespace bcl
