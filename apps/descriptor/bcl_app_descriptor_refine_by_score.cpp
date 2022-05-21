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

// include the header of this class
#include "bcl_app_descriptor_refine_by_score.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_statistics.h"
#include "model/bcl_model_collect_features_interface.h"
#include "model/bcl_model_data_set_score.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief standard constructor
    DescriptorRefineByScore::DescriptorRefineByScore() :
      m_FlagDataSetScoreFilename
      (
        new command::FlagStatic
        (
          "score_file",
          "file containing the dataset score object",
          command::Parameter
          (
            "source",
            "",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_FlagSelectionMethod
      (
        new command::FlagStatic
        (
          "select",
          "method of choosing descriptors to keep",
          command::Parameter
          (
            "",
            "",
            command::ParameterCheckSerializable( util::Implementation< model::CollectFeaturesInterface>()),
            "Top(10)"
          )
        )
      ),
      m_FlagOutputFilename
      (
        new command::FlagStatic
        (
          "output",
          "filename for output features",
          command::Parameter
          (
            "output",
            "filename for output features",
            "descriptors.refined.object"
          )
        )
      ),
      m_FlagOutputInfo
      (
        new command::FlagStatic
        (
          "info",
          "if set, also write out a text file containing information with # feature columns removed, "
          "remaining, and which were removed"
        )
      ),
      m_FlagCompareScoreFiles
      (
        new command::FlagStatic
        (
          "compare_score_file",
          "filename for score file to compare",
          command::Parameter
          (
            "compare_score_file",
            "filename for score file to compare",
            ""
          )
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new RefineDescriptors
    DescriptorRefineByScore *DescriptorRefineByScore::Clone() const
    {
      return new DescriptorRefineByScore( *this);
    }

    // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
    // but which should not be displayed, e.g. for help
    storage::Vector< std::string> DescriptorRefineByScore::GetDeprecatedAppNames() const
    {
      return storage::Vector< std::string>( size_t( 1), "RefineDescriptors");
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string DescriptorRefineByScore::GetDescription() const
    {
      return "Create a new descriptor file using the best descriptors by score";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &DescriptorRefineByScore::GetReadMe() const
    {
      // construct static readme
      static const std::string s_readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of the BCL application descriptor:RefineByScore, its terms of use, "
        "appropriate citation, installation procedures, descriptor:RefineByScore execution, and technical support\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. What is descriptor:RefineByScore?\n"
        "descriptor:RefineByScore is a C++ based application, created by Vanderbilt University's Meiler Laboratory, "
        "which is part a larger library of applications called BCL::Commons.  descriptor:RefineByScore is used create "
        "descriptor files using score files generated by descriptor:ScoreDataset.  The descriptors may be selected by "
        "one of two (currently) methods: either the top descriptors by score or all descriptors above a certain score "
        "threshold. Secondarily, descriptor:RefineByScore can be used to compute the spearman correlation coefficient "
        "between two descriptor files."
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING descriptor:RefineByScore.\n"
        "When using descriptor:RefineByScore in a publication, please cite the following publications describing the "
        "application's development:\n"
        "\n"
        "Butkiewicz M, Lowe EW, Mueller R, Mendenhall JL, Teixeira PL, Weaver CD, Meiler J. "
        "Benchmarking Ligand-Based Virtual High-Throughput Screening with the PubChem Database "
        "Molecules, 18, (1), 735-756. ; 2013\n"
        "Link:  www.http://meilerlab.org/index.php/publications/show/2013\nJournal link: http://www.mdpi.com/1420-3049/18/1/735\n"
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING descriptor:RefineByScore.\n"
        "Running descriptor:RefineByScore is very straightforward; the three basic command lines that can be used are:\n"
        "1. bcl.exe descriptor:RefineByScore -output new_descriptors.object -score_file infogain.scores -select 'Top(50)'\n"
        "  This selects the top 50 descriptors from infogain.scores and places them in new_descriptors.object\n"
        "2. bcl.exe descriptor:RefineByScore -output new_descriptors.object -score_file infogain.scores -select 'Above(0.05)'\n"
        "  This selects the descriptors with scores above 0.05 from infogain.scores and places them in new_descriptors.object\n"
        "3. bcl.exe descriptor:RefineByScore -score_file infogain.scores -compare_score_file fscore.scores\n"
        "  This command the spearman correlation coefficient between infogain.scores and fscore.scores\n"
        "\nADDITIONAL INFORMATION\n"
        "\n"
        "Helpful information concerning syntax and flags can be obtained by typing\n"
        "bcl.exe descriptor:RefineByScore  -help\n"
        "\n"
        "For this read me, type\n"
        "bcl.exe descriptor:RefineByScore -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF descriptor:RefineByScore.\n"
        "descriptor:RefineByScore is under ongoing further development. For current research please refer to "
        "www.meilerlab.org and navigate to research\n"
        + DefaultSectionSeparator()
      );
      return s_readme;
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> DescriptorRefineByScore::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      sp_cmd->AddFlag( m_FlagDataSetScoreFilename);
      sp_cmd->AddFlag( m_FlagSelectionMethod);
      sp_cmd->AddFlag( m_FlagOutputFilename);
      sp_cmd->AddFlag( m_FlagOutputInfo);
      sp_cmd->AddFlag( m_FlagCompareScoreFiles);

    ///////////////////
    // default flags //
    ///////////////////

      // default flags are unnecessary for this application, but message level and read-me are useful
      command::GetAppDefaultFlags().AddRequiredCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int DescriptorRefineByScore::Main() const
    {
      model::DataSetScore score;

      // get the score file
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_FlagDataSetScoreFilename->GetFirstParameter()->GetValue());
      io::Serialize::Read( score, input);
      io::File::CloseClearFStream( input);

      if( m_FlagCompareScoreFiles->GetFlag())
      {
        model::DataSetScore score_compare;

        // get the score file for comparison
        io::File::MustOpenIFStream( input, m_FlagCompareScoreFiles->GetFirstParameter()->GetValue());
        io::Serialize::Read( score_compare, input);
        io::File::CloseClearFStream( input);

        const float spearman_corr_coeff
        (
          math::Statistics::CorrelationSpearman
          (
            score.GetScores().Begin(),
            score.GetScores().End(),
            score_compare.GetScores().Begin(),
            score_compare.GetScores().End()
          )
        );

        BCL_MessageStd
        (
          "Ranking comparison of highest weighted columns for score file " + m_FlagCompareScoreFiles->GetFirstParameter()->GetValue()
          + "\n with score file " + m_FlagCompareScoreFiles->GetFirstParameter()->GetValue()
          + "\n\n"
          + "Spearman Correlation Coefficient: " + util::Format()( spearman_corr_coeff)
          + "\n"
        );

        return 0;
      }

      // get the labels
      const model::FeatureLabelSet labels( score.GetFeatures());

      // create an implementation to select the best columns
      util::Implementation< model::CollectFeaturesInterface>
        selector( m_FlagSelectionMethod->GetFirstParameter()->GetValue());

      // select the scores
      const storage::Vector< size_t> selected_columns( selector->Collect( score.GetScores()));

      // update the labels
      const model::FeatureLabelSet refined_labels( labels.CreateSubFeatureLabelSet( selected_columns));

      // write out the labels
      io::OFStream output;
      io::File::MustOpenOFStream( output, m_FlagOutputFilename->GetFirstParameter()->GetValue());
      output << util::ObjectDataLabel( refined_labels.GetString()).ToStringDefaultWidth();
      io::File::CloseClearFStream( output);

      if( m_FlagOutputInfo->GetFlag())
      {
        io::OFStream output;
        const std::string info_filename
        (
          io::File::RemoveLastExtension( m_FlagOutputFilename->GetFirstParameter()->GetValue()) + ".info"
        );
        io::File::MustOpenOFStream( output, info_filename);
        const size_t n_features_original( score.GetScores().GetSize());
        output << n_features_original - selected_columns.GetSize() << " columns removed\n";
        output << selected_columns.GetSize() << " columns remain\n";
        output << "Removed features:\n";
        storage::Vector< size_t> removed_cols;

        size_t feature( 0);
        for
        (
          storage::Vector< size_t>::const_iterator
            itr_selected( selected_columns.Begin()), itr_selected_end( selected_columns.End());
          feature < n_features_original && itr_selected != itr_selected_end;
          ++feature
        )
        {
          if( feature != *itr_selected)
          {
            removed_cols.PushBack( feature);
          }
          else
          {
            ++itr_selected;
          }
        }
        while( feature < n_features_original)
        {
          removed_cols.PushBack( feature);
          ++feature;
        }
        const model::FeatureLabelSet removed_labels( labels.CreateSubFeatureLabelSet( removed_cols));

        // write out the labels
        output << util::ObjectDataLabel( removed_labels.GetString()).ToStringDefaultWidth();
        io::File::CloseClearFStream( output);
      }

      // end
      return 0;
    }

    const ApplicationType DescriptorRefineByScore::s_Instance
    (
      GetAppGroups().AddAppToGroup( new DescriptorRefineByScore(), GetAppGroups().e_Descriptor)
    );

  } // namespace app
} // namespace bcl
