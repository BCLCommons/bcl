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
#include "bcl_app_descriptor_score_dataset.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_limits.h"
#include "model/bcl_model_data_set_score.h"
#include "model/bcl_model_retrieve_data_set_base.h"
#include "model/bcl_model_score_dataset_interface.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief standard constructor
    DescriptorScoreDataset::DescriptorScoreDataset() :
      m_FlagDataSet
      (
        new command::FlagStatic
        (
          "source",
          "method to retrieve the dataset",
          command::Parameter
          (
            "source",
            "",
            command::ParameterCheckSerializable( util::Implementation< model::RetrieveDataSetBase>())
          )
        )
      ),
      m_FlagFeatureCode
      (
        new command::FlagStatic
        (
          "feature_labels",
          "label or file containing the label for the feature descriptors",
          command::Parameter
          (
            "feature_labels",
            "",
            ""
          )
        )
      ),
      m_FlagResultCode
      (
        new command::FlagStatic
        (
          "result_labels",
          "label or file containing the result label",
          command::Parameter
          (
            "result_labels",
            "",
            ""
          )
        )
      ),
      m_FlagOutputFilename
      (
        new command::FlagStatic
        (
          "output",
          "filename for output scores",
          command::Parameter
          (
            "output",
            "filename for output scores"
          )
        )
      ),
      m_FlagScore
      (
        new command::FlagStatic
        (
          "score",
          "how to score the dataset",
          command::Parameter
          (
            "score",
            "method to score the dataset",
            command::ParameterCheckSerializable( util::Implementation< model::ScoreDatasetInterface>())
          )
        )
      ),
      m_FlagTerse
      (
        new command::FlagStatic
        (
          "terse",
          "set to true to write a short score file, skipping the writing of the sorted descriptors, which are "
          "for user reference only. "
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new DescriptorScoreDataset
    DescriptorScoreDataset *DescriptorScoreDataset::Clone() const
    {
      return new DescriptorScoreDataset( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string DescriptorScoreDataset::GetDescription() const
    {
      return "Scores feature columns using various criteria. Used during descriptor development and selection";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &DescriptorScoreDataset::GetReadMe() const
    {
      // construct static readme
      static const std::string s_readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of the BCL application descriptor:ScoreDataset, its terms of use, "
        "appropriate citation, installation procedures, descriptor:ScoreDataset execution, and technical support\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. What is descriptor:ScoreDataset?\n"
        "descriptor:ScoreDataset is a C++ based application, created by Vanderbilt University's Meiler Laboratory, "
        "which is part a larger library of applications called BCL::Commons.  descriptor:ScoreDataset is used to "
        "score descriptor columns, based on e.g. information gain, fscore, input sensitity, etc., and various elementary"
        "operations of these elementary scores. "
        "\nThe output of this application is a score file, which contains a listing of all the in the dataset, "
        "their sizes, and a vector of scores for all members of the dataset.  This information is used by other bcl "
        "applications (primarily descriptor:RefineByScore) for descriptor selection.\n"
        "The score file's final section contains a list of the descriptors (column 2) sorted by ascending scores (column 1). "
        "While not directly used by the bcl, this section makes it easy to obtain the scores for each descriptor. "
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING descriptor:ScoreDataset.\n"
        "When using descriptor:ScoreDataset in a publication, please cite the following publications describing the "
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
        "VI. RUNNING descriptor:ScoreDataset.\n"
        "\n"
        "1) Specify your data source(s):\n"
        "This is what you will pass to the -source flag.  Detailed explanations of data sources and the overall framework "
        "are available via bcl.exe descriptor:GenerateDataset -readme, under Running descriptor:GenerateDataset\n"
        "The data source for this application is usually .bin files generated by descriptor:GenerateDataset, owing "
        "to their high speed of reading, though any dataset retriever can be used"
        "There are many dataset retrievers in the bcl, so type bcl.exe descriptor:ScoreDataset -help shows a full list. "
        "\n"
        "2) Choose which descriptors to construct:\n"
        "These will be passed to the -feature_labels and -result_labels flags.\n"
        "See bcl.exe descriptor:GenerateDataset -readme for details\n"
        "3) Select the score type\n"
        "A listing of all scores is available, with bcl.exe descriptor:ScoreDataset -help\n"
        "The most common methods of scoring datasets for classification type tasks "
        "(e.g. training a model to learn yes-no like outputs) are information gain, gini index, f-score, and input sensitivity.\n"
        "For regression-type tasks, input-sensitivity and correlation (coming soon) are usually the best fits, though "
        "other scores may work well enough if a cutoff value can be obtained for the result that split the data into two "
        "roughly-equal halves\n"
        "For input sensitivity, one must already have trained a model (or a set of models) on the same descriptor set.\n"
        "4) Run descriptor:ScoreDataset\n"
        "At a command prompt, navigate to the location of your BCL executable program. The syntax for running "
        "the application looks like the following:\n"
        "\n"
        "bcl.exe descriptor:ScoreDataset -source <source-specifier> -output mydataset.scores "
        "-feature_labels <label or file with labels> -result_labels <label or file with labels> -score <score-specifier>"
        "For file-based datasets (e.g. those specified with the File, or Subset retriever), the feature and result_labels "
        "flags are usually omitted, resulting in scoring of all feature columns against all result columns. "
        "For certain applications it is useful to score against a subset of the results, especially if different thresholds "
        "are used for different result columns."
        "Example command lines: "
        "bcl.exe descriptor:ScoreDataset -source \"Subset(filename=dataset.bin)\" -output dataset_ig.scores "
        "-score \"Partition( partitioner = InformationGain, cutoff = 0.5)\"\n"
        "It is necessary to adjust the cutoff value "
        " to match the threshold between classification states in the results of your dataset. "
        "Some alternative scores that are often used: "
        "1. \"FScore(cutoff=0.5)\"\n"
        "2. \"Partition(partitioner=Gini,cutoff=0.5)\"\n"
        "3. \"Multiply(FScore(cutoff=0.5), Partition( partitioner = InformationGain, cutoff = 0.5))\" Multiplication benefits scores that are rated high by both scorers"
        "4. \"InputSensitivity( storage=File(directory=./models/,prefix=model), delta = 0.1)\""
        "5. \"InputSensitivity( storage=File(directory=./models/,key=2,prefix=model), delta = 0.1)\""
        "Notes on using input sensitivity: \n"
        "  * the directory and prefix parameters of File must be the same as when TrainModel was called.\n"
        "  * The optional key parameter (in the storage), corresponds to the number in the filename, e.g. model000001.model has key=1\n"
        "  * If no key is specified, all models in the given storage are used.  One can speed up calculation by selecting a representative model\n"
        "  * the delta determines the amount that individual feature columns are perturbed\n"
        "  ** Use small values, like 0.01 - 0.05 for trained ANNs and SVMs; medium values (0.5 - 1) for kohonen networks, and large values (1 - 10) for decision trees, which have discontinuous derivatives\n"
        "  * Use the number chunks and chunk parameter of your dataset source to use only a part of the dataset for massive datasets"
        "Input sensitivity is usually quite stable after 10,000 - 100,000 features."
        "\nADDITIONAL INFORMATION\n"
        "\n"
        "Helpful information concerning syntax and flags can be obtained by typing\n"
        "bcl.exe descriptor:ScoreDataset  -help\n"
        "\n"
        "For this read me, type\n"
        "bcl.exe descriptor:ScoreDataset -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF descriptor:ScoreDataset.\n"
        "descriptor:ScoreDataset is under ongoing further development. For current research please refer to "
        "www.meilerlab.org and navigate to research\n"
        + DefaultSectionSeparator()
      );
      return s_readme;
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> DescriptorScoreDataset::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      sp_cmd->AddFlag( m_FlagDataSet);
      sp_cmd->AddFlag( m_FlagFeatureCode);
      sp_cmd->AddFlag( m_FlagResultCode);
      sp_cmd->AddFlag( m_FlagOutputFilename);
      sp_cmd->AddFlag( m_FlagScore);
      sp_cmd->AddFlag( m_FlagTerse);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int DescriptorScoreDataset::Main() const
    {
      // initialize monitor data
      util::Implementation< model::RetrieveDataSetBase> source( m_FlagDataSet->GetFirstParameter()->GetValue());

      // read in feature / result codes
      source->SelectFeaturesGivenFilenameFlag( *m_FlagFeatureCode);
      source->SelectResultsGivenFilenameFlag( *m_FlagResultCode);

      util::Implementation< model::ScoreDatasetInterface> scorer( m_FlagScore->GetFirstParameter()->GetValue());

      model::DataSetScore score;
      score.SetFeatures( source->GetFeatureLabelsWithSizes());

      // generate the features and results for this dataset
      util::ShPtr< descriptor::Dataset> dataset( source->GenerateDataSet());
      score.SetScores( scorer->Score( *dataset));

      io::OFStream output;
      io::File::MustOpenOFStream( output, m_FlagOutputFilename->GetFirstParameter()->GetValue());
      io::Serialize::Write( score, output) << '\n';

      if( !m_FlagTerse->GetFlag())
      {
        output << "\nSorted Scores:\n";

        // split the feature result dataset, to write the sorted scores as well to the file
        model::FeatureLabelSet features( source->GetFeatureLabelsWithSizes().SplitFeatureLabelSet( true));

        // create a map with the score value and the descriptor index to sort the scores
        std::map< std::pair< float, size_t>, util::ObjectDataLabel> sorted_scores;
        for
        (
          size_t descriptor_index( 0), n_descriptor_cols( features.GetMemberLabels().GetSize());
          descriptor_index < n_descriptor_cols;
          ++descriptor_index
        )
        {
          // add the feature to the map, indexed by score and object data label
          sorted_scores[ std::make_pair( score.GetScores()( descriptor_index), descriptor_index)]
            = features.GetMemberLabels()( descriptor_index);
        }

        for
        (
          std::map< std::pair< float, size_t>, util::ObjectDataLabel>::const_iterator
            itr( sorted_scores.begin()), itr_end( sorted_scores.end());
          itr != itr_end;
          ++itr
        )
        {
          output << itr->first.first << ' ' << itr->second.ToString() << '\n';
        }
      }
      io::File::CloseClearFStream( output);
      // end
      return 0;
    }

    const ApplicationType DescriptorScoreDataset::s_Instance
    (
      GetAppGroups().AddAppToGroup( new DescriptorScoreDataset(), GetAppGroups().e_Descriptor)
    );

  } // namespace app
} // namespace bcl
