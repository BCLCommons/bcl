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
#include "bcl_app_descriptor_generate_dataset.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "model/bcl_model_data_set_statistics.h"
#include "model/bcl_model_retrieve_dataset_subset.h"
#include "model/bcl_model_score_dataset_non_redundant.h"
#include "storage/bcl_storage_table.h"
#include "util/bcl_util_loggers.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  namespace app
  {

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief standard constructor
    DescriptorGenerateDataset::DescriptorGenerateDataset() :
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
            command::ParameterCheckSerializable( util::Implementation< model::RetrieveDataSetBase>()),
            ""
          )
        )
      ),
      m_FlagFeatureCode
      (
        new command::FlagStatic
        (
          "feature_labels",
          "label or file containing the label for the feature descriptors",
          command::Parameter( "feature_labels", "", "")
        )
      ),
      m_FlagResultCode
      (
        new command::FlagStatic
        (
          "result_labels",
          "label or file containing the label for the result descriptors",
          command::Parameter( "result_labels", "", "")
        )
      ),
      m_FlagIDCode
      (
        new command::FlagStatic
        (
          "id_labels",
          "label or file containing the id label",
          command::Parameter( "id_labels", "", "")
        )
      ),
      m_FlagOutputFilename
      (
        new command::FlagStatic
        (
          "output",
          "filename for output dataset with suffix .bin (binary) or .csv (comma separated values) or .arff (for use with WEKA)",
          command::Parameter
          (
            "output",
            "filename for output dataset with suffix .bin (binary) or .csv (comma separated values) or .arff (for use with WEKA)",
            ""
          )
        )
      ),
      m_FlagStatisticsFilenames
      (
        new command::FlagDynamic
        (
          "compare",
          "filenames to compare generated output bin file with",
          command::Parameter
          (
            "filename",
            "bin file for comparison",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_FlagInfo
      (
        new command::FlagStatic
        (
          "info",
          "if set, do not generate the given dataset, just retrieve # features, results, memory, etc. for the dataset"
        )
      ),
      m_FlagRedundancyRemoval
      (
        new command::FlagStatic
        (
          "nonredundant",
          "Allows dataset to be post-filtered to remove redundant descriptors; criteria for which can be set in the flag",
          command::Parameter
          (
            "NonRedundantParams",
            "criteria for deciding whether a descriptor is redundant or not, type -nonredundant help for detailed "
            "description",
            command::ParameterCheckSerializable( model::ScoreDatasetNonRedundant()),
            "(max outliers=20,tol=0.05,min span=0.001,min std=0.0001,min rsq=0.95)"
          )
        )
      ),
      m_FlagBlockSize
      (
        new command::FlagStatic
        (
          "block_size",
          "number of MB in the output dataset to generate before writing out. Use smaller values if running out of memory,"
          "which may happen when descriptors are predictions from models",
          command::Parameter
          (
            "block_size",
            "number of MB in the output dataset to generate before writing out. Use smaller values if running out of memory,"
            "which may happen when descriptors are predictions from models",
            command::ParameterCheckRanged< double>( 1.0e-6, 128.0),
            "32"
          )
        )
      ),
      m_FlagForbidIncompleteResults
      (
        new command::FlagStatic
        (
          "forbid_incomplete_records",
          "whether to forbid records that have incomplete results (some, but not all, results nan). "
          "By default, records with all result nan are omitted"
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new GenerateDataset
    DescriptorGenerateDataset *DescriptorGenerateDataset::Clone() const
    {
      return new DescriptorGenerateDataset( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string DescriptorGenerateDataset::GetDescription() const
    {
      return "Generates bin files (fast, small, easier to use with other applications), csv files (human readable), or "
              " .arff files (for use with WEKA) from any combination of dataset sources."
             "Can compute dataset statistics and compare bin files";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &DescriptorGenerateDataset::GetReadMe() const
    {
      // construct static readme
      static const std::string s_readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of the BCL application descriptor:GenerateDataset, its terms of use, "
        "appropriate citation, installation procedures, descriptor:GenerateDataset execution, and technical support\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. What is descriptor:GenerateDataset?\n"
        "descriptor:GenerateDataset is a C++ based application, created by Vanderbilt University's Meiler Laboratory, "
        "which is part a larger library of applications called BCL::Commons.  descriptor:GenerateDataset is used to "
        "generate files containing numerical descriptions of e.g. molecules, atoms, proteins, and amino acid sequences."
        "  Three types of dataset output formats are supported.  The first is .csv files (comma-separated value), which "
        "are easy to read, and intended solely for use external to the bcl. Related to csv format is the .arff format "
        "which is the native format for WEKA (http://www.cs.waikato.ac.nz/ml/weka/), another machine learning package. \n"
        "The third and most common format written by this application is .bin files, which are ideal for use with other "
        "applications in the descriptor: and model: groups, owing to its small size and binary"
        ", fixed width format that make it several times faster for the bcl to read than ASCII files.  Other "
        "applications read .bin files using the Subset retriever, e.g. Subset(filename=x.bin).\n"
        "descriptor:GenerateDataset can also compute statistics on each descriptor in existing .bin files and pairwise "
        "statistics when the -compare flag is used with multiple .bin files\n"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING descriptor:GenerateDataset.\n"
        "When using descriptor:GenerateDataset in a publication, please cite the following publications describing the "
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
        "VI. RUNNING descriptor:GenerateDataset.\n"
        "Running descriptor:GenerateDataset consists of three main steps.\n"
        "\n"
        "1) Specify your data source(s):\n"
        "This is what you will pass to the -source flag. "
        "There are many dataset retrievers in the bcl, so type bcl.exe descriptor:GenerateDataset -help shows a full list. "
        "Typically, the primary data source comes from files containing one of the following:\n"
        "1. molecules (in sdf format) using SdfFile\n"
        "2. sequences (in fasta format) using SequenceDirectory\n"
        "3. proteins (in pdb format) using ProteinDirectory retriever\n"
        "4. another bcl application that writes out native bcl vector formats using the File retriever.\n"
        "5. bcl-native dataset format known as the .bin file (using the Subset retriever).\n"
        "Other dataset retrievers perform operations on the dataset before storing it, such as Combine, which directly\n"
        "concatenates datasets, and KMeans, which performs clustering. Read the help for a complete list of available sources\n"
        "Whichever dataset retrievers you decide to use, read the help for them using "
        "e.g. bcl.exe descriptor:GenerateDataset -source \"SdfFile(help)\" and use it to specify the internal parameters of the retriever"
        "\n"
        "2) Choose which descriptors to construct:\n"
        "\nA. Background\n"
        "Descriptors are broadly divided into three categories in the bcl:\n"
        "\"Feature\" descriptors provide numerical values specific to the given objects (e.g. molecules or proteins) "
        "that are commonly used as the input to trained machine learning algorithms.\n"
        "\"Result\" descriptors are just like Feature descriptors except that they are the desired output values "
        "(commonly the experimental values) for the machine learning algorithm\n"
        "\"ID\" descriptors are used primarily by you, the end user, to identify what each feature / result row describes"
        "  They may also be used by certain descriptors, see their descriptions for more info\n"
        "Descriptors of all types are given with the -feature_labels flag (for features), which takes either a filename or"
        "the descriptor directly.  Results and ID descriptors can be given using separate flags.\n"
        "For bcl-native vector format, the only descriptors available are the column indices (0-indexed)\n"
        "For bin files, the provided descriptors must be a subset of those used when the file was originally created.\n"
        "If no descriptors are explicitly given for vector or bin files, all descriptors from them will be used\n\n"
        "For all other primary sources, feature and results descriptors are required (ids are always optional).\n"
        "To see a list of all available descriptors that can be created for your primary source (other than .bin files "
        "and vector files), type\n"
        "bcl.exe descriptor:GenerateDataset -source _PrimarySource_ -feature_labels help\n"
        "e.g.\n"
        "bcl.exe descriptor:GenerateDataset -source \"ProteinDirectory(.)\" -feature_labels help\n"
        "Similar commands can be used to obtain all available descriptors for results (identical to features) and ids\n"
        "A descriptor file that computes a few descriptors for molecules is as simple as:\n"
        "Combine( NAtoms, Polarizability, 2DA( property = Atom_Identity, steps=12, normalized = 0))\n"
        "Some of the descriptor files used in the lab during research routinely contain several thousand individual\n"
        "descriptors though; and since descriptors can be combined and numerous operations like addition, multiplication,"
        "statistics, etc., the possibilities are infinite\n"
        "Notes on syntax for descriptor labels:\n"
        "  a. Spacing is irrelevant except that if the label is split onto multiple lines, the first line must contain\n"
        "at least the opening parenthesis, e.g.: Combine\n( NAtoms)\nwill not work because the first line has balanced "
        "() so the parser does not continue.\n"
        "  b. Where a string is expected, characters such as (), or spaces require double quotes\n"
        "3) Run descriptor:GenerateDataset:\n"
        "At a command prompt, navigate to the location of your BCL executable program. The syntax for running "
        "the application looks like the following:\n"
        "\n"
        "bcl.exe descriptor:GenerateDataset -source <source-specifier> -output dataset.(bin|csv|arff) "
        "-feature_labels <label or file with labels> -result_labels <label or file with labels> -id_labels <label or filename with labels> "
        "\n"
        "The -compare flag may also be given, without any arguments, to compute statistics for the generated file."
        "e.g."
        "bcl.exe descriptor:GenerateDataset -source \"Subset(filename=dataset.bin)\" -compare\n"
        "or -compare may be given separately in the absence of a source argument to compare multiple bin files:\n"
        "including statistics on the average difference between descriptor values, e.g.\n"
        "bcl.exe descriptor:GenerateDataset -compare MyDataset.bin MyOtherDataset.bin\n"
        "\nif you have a bin file and want to be able to read its values using other tools, write it out as:\n"
        "bcl.exe descriptor:GenerateDataset -source \"Subset(filename=dataset.bin)\" -output dataset.csv\n"
        "Note that the extension of the file given after -output controls the output type: .bin for bin files, .csv for comma separated values\n"
        "\nADDITIONAL INFORMATION\n"
        "\n"
        "Helpful information concerning syntax and flags can be obtained by typing\n"
        "bcl.exe descriptor:GenerateDataset  -help\n"
        "\n"
        "For this read me, type\n"
        "bcl.exe descriptor:GenerateDataset -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF descriptor:GenerateDataset.\n"
        "descriptor:GenerateDataset is under ongoing further development. For current research please refer to "
        "www.meilerlab.org and navigate to research\n"
        + DefaultSectionSeparator()
      );
      return s_readme;
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> DescriptorGenerateDataset::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      sp_cmd->AddFlag( m_FlagDataSet);
      sp_cmd->AddFlag( m_FlagFeatureCode);
      sp_cmd->AddFlag( m_FlagResultCode);
      sp_cmd->AddFlag( m_FlagIDCode);
      sp_cmd->AddFlag( m_FlagOutputFilename);
      sp_cmd->AddFlag( m_FlagStatisticsFilenames);
      sp_cmd->AddFlag( m_FlagInfo);
      sp_cmd->AddFlag( m_FlagRedundancyRemoval);
      sp_cmd->AddFlag( m_FlagBlockSize);
      sp_cmd->AddFlag( m_FlagForbidIncompleteResults);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int DescriptorGenerateDataset::Main() const
    {
      std::string output_filename( m_FlagOutputFilename->GetFirstParameter()->GetValue());

      const bool allow_incomplete( !m_FlagForbidIncompleteResults->GetFlag());
      // check whether the user wants a new bin file created
      if( m_FlagDataSet->GetFlag())
      {
        // initialize the source dataset retriever
        util::Implementation< model::RetrieveDataSetBase> source( m_FlagDataSet->GetFirstParameter()->GetValue());

        // set up the feature result codes
        source->SelectFeaturesGivenFilenameFlag( *m_FlagFeatureCode);
        source->SelectResultsGivenFilenameFlag( *m_FlagResultCode);
        source->SelectIdsGivenFilenameFlag( *m_FlagIDCode);

        if( m_FlagOutputFilename->GetFlag() && !m_FlagInfo->GetFlag())
        {
          // for non-redundant removal, ensure that a bin file is being created
          if( m_FlagRedundancyRemoval->GetFlag())
          {
            BCL_Assert
            (
              util::EndsWith( m_FlagOutputFilename->GetFirstParameter()->GetValue(), ".bin"),
              "Redundancy removal requires the output to be in a .bin file"
            );
            // write to a temporary bin file first
            output_filename = output_filename.substr( 0, output_filename.size() - 4) + ".tmp.bin";
          }
          const double block_size( m_FlagBlockSize->GetNumericalList< double>()( 0));
          // store the dataset
          model::RetrieveDatasetSubset::StoreMasterDataset( output_filename, *source, block_size, allow_incomplete);
        }
        else if( m_FlagInfo->GetFlag())
        {
          // get the desired information from the source
          const size_t size( source->GetNominalSize());

          const size_t number_result_cols( source->GetResultCodeWithSizes().GetSize());
          const size_t number_feature_cols( source->GetFeatureLabelsWithSizes().GetSize());
          const size_t number_id_cols( source->GetIdCodeWithSizes().GetSize());
          const size_t bytes_per_row
          (
            sizeof( float) * ( number_result_cols + number_feature_cols) + sizeof( char) * number_id_cols
          );
          const size_t total_megabytes( ( ( bytes_per_row * size) >> 20) + 1);

          // output the info to the given file, or the screen, if desired
          io::OFStream output;
          if( m_FlagOutputFilename->GetFlag())
          {
            io::File::MustOpenOFStream( output, m_FlagOutputFilename->GetFirstParameter()->GetValue());
          }
          std::ostream &output_info_stream
          (
            m_FlagOutputFilename->GetFlag()
            ? (std::ostream &)output
            : (std::ostream &)( util::GetLogger())
          );
          output_info_stream << "Information for: " << source->GetLabel().ToString() << '\n';
          output_info_stream << total_megabytes << "MB\n";
          output_info_stream << size << " rows\n";
          output_info_stream << number_feature_cols << " feature columns\n";
          output_info_stream << number_result_cols << " result columns\n";
          output_info_stream << number_id_cols << " id columns" << std::endl;

          if( m_FlagOutputFilename->GetFlag())
          {
            io::File::CloseClearFStream( output);
          }
        }
        else
        {
          BCL_Exit( "-source and -output must be given together", -1);
        }
      }

      // check whether the user wants to take statistics
      if
      (
        m_FlagStatisticsFilenames->GetFlag()
        &&
        (
          m_FlagStatisticsFilenames->GetParameterList().GetSize()
          || ( m_FlagOutputFilename->GetFlag() && m_FlagDataSet->GetFlag() && !m_FlagInfo->GetFlag())
        )
      )
      {
        WriteStatistics( util::GetLogger());
      }

      // check whether the user requested non-redundant descriptor removal
      if( m_FlagRedundancyRemoval->GetFlag())
      {
        const std::string final_outputf( m_FlagOutputFilename->GetFirstParameter()->GetValue());
        const std::string object_filename( final_outputf.substr( 0, final_outputf.size() - 4) + ".nonredundant.obj");
        const std::string object_filelog( final_outputf.substr( 0, final_outputf.size() - 4) + ".nonredundant.log");

        BCL_MessageStd( "Initial dataset generated. Removing redundant descriptors");
        model::ScoreDatasetNonRedundant nr;
        nr.AssertRead( m_FlagRedundancyRemoval->GetFirstParameter()->GetValue());
        linal::Vector< float> scores;
        model::FeatureLabelSet feature_labels;
        {
          model::RetrieveDatasetSubset subset_retriever( output_filename);
          util::ShPtr< descriptor::Dataset> dataset( subset_retriever.GenerateDataSet());
          feature_labels = *dataset->GetFeatures().GetFeatureLabelSet();
          bool set_logger( false);
          if( !util::GetLoggers().GetFlagLogger()->GetFlag())
          {
            BCL_MessageStd( "See " + object_filelog + " for details of which descriptors were removed and why");
            set_logger = true;
            std::stringstream er;
            remove( object_filelog.c_str());
            util::GetLoggers().GetFlagLogger()->ReadFromList
            (
              storage::Vector< std::string>::Create( "File", object_filelog),
              er
            );
          }
          scores = nr.Score( *dataset);
          if( set_logger)
          {
            util::GetLoggers().GetFlagLogger()->ResetFlag();
          }
        }
        storage::Vector< size_t> non_redundant_indices( size_t( scores.Sum()));
        // select the non-redundant indices
        const size_t number_non_redundant_desc( non_redundant_indices.GetSize());
        for
        (
          size_t nonredundant_col_n( 0), col_n( 0);
          nonredundant_col_n < number_non_redundant_desc;
          ++nonredundant_col_n, ++col_n
        )
        {
          while( scores( col_n) < 0.5)
          {
            ++col_n;
          }
          non_redundant_indices( nonredundant_col_n) = col_n;
        }
        const double nr_percent( float( number_non_redundant_desc) / float( scores.GetSize()) * 100.0);
        BCL_MessageStd
        (
          util::Format()( number_non_redundant_desc)
          + " (" + util::Format().FFP( 3)( nr_percent) + "%) "
          "of the original descriptors were non-redundant"
        );
        BCL_MessageStd
        (
          util::Format()( scores.GetSize() - number_non_redundant_desc)
          + " (" + util::Format().FFP( 3)( 100.0 - nr_percent) + "%) "
          "of the original descriptors were redundant (criteria: " + nr.GetLabel().ArgumentsToString( false) + ")"
        );
        model::RetrieveDatasetSubset subset_retriever( output_filename);

        // update the labels
        const model::FeatureLabelSet refined_labels( feature_labels.CreateSubFeatureLabelSet( non_redundant_indices));

        subset_retriever.SelectFeatures( refined_labels.GetLabel());
        subset_retriever.SelectResults( subset_retriever.GetResultCodeWithSizes().GetLabel());

        // store the dataset
        model::RetrieveDatasetSubset::StoreMasterDataset( final_outputf, subset_retriever, 32.0, allow_incomplete);
        remove( output_filename.c_str());

        // write out the updated descriptor file
        io::OFStream os;
        io::File::MustOpenOFStream( os, object_filename);
        os << refined_labels.GetLabel().ToNamedString( 100, 0, 1);
        io::File::CloseClearFStream( os);
        BCL_MessageStd( object_filename + " was written with the non-redundant descriptors");
      }

      // end
      return 0;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief compute statistics
    //! @param OSTREAM output stream for statistics
    void DescriptorGenerateDataset::WriteStatistics( std::ostream &OSTREAM) const
    {
      // generate stats for each file
      // also, if the # of rows is equal across all files, also take the standard deviation
      storage::Vector< std::string> filenames;
      if( m_FlagDataSet->GetFlag() && m_FlagOutputFilename->GetFlag())
      {
        filenames.PushBack( m_FlagOutputFilename->GetFirstParameter()->GetValue());
      }
      filenames.Append( m_FlagStatisticsFilenames->GetStringList());

      // default headings; always compute these statistics
      const storage::Vector< std::string> default_headings
      (
        storage::Vector< std::string>::Create( "Ave", "Std", "Min", "Max")
      );

      // storage for all the datasets
      util::ShPtrList< descriptor::Dataset> datasets;

      // track whether all datasets have the same # of rows and columns
      bool all_rows_cols_equal( true);

      // set these members after loading the first dataset
      size_t number_rows( util::GetUndefined< size_t>());
      size_t number_cols( util::GetUndefined< size_t>());
      model::FeatureLabelSet features_with_sizes;
      model::FeatureLabelSet results_with_sizes;

      // iterate over all filenames to load datasets
      for
      (
        storage::Vector< std::string>::const_iterator itr_a( filenames.Begin()), itr_end( filenames.End());
        itr_a != itr_end;
        ++itr_a
      )
      {
        // determine whether the file is txt or binary; use the corresponding dataset retriever
        const std::string alias_a
        (
          itr_a->size() < size_t( 4) || itr_a->substr( itr_a->size() - 4) == ".txt"
          ? "File(filename="
          : "Subset(filename="
        );
        util::Implementation< model::RetrieveDataSetBase> impl_a( alias_a + *itr_a + ")");
        impl_a->SelectFeaturesGivenFilenameFlag( *m_FlagFeatureCode);
        impl_a->SelectResultsGivenFilenameFlag( *m_FlagResultCode);
        impl_a->SelectIdsGivenFilenameFlag( *m_FlagIDCode);

        BCL_Assert( impl_a.IsDefined(), "undefined implementation of dataset retriever");

        datasets.PushBack( impl_a->GenerateDataSet());
        if( !util::IsDefined( number_rows))
        {
          number_rows = datasets.LastElement()->GetSize();
          number_cols = datasets.LastElement()->GetFeatureSize();
          features_with_sizes = impl_a->GetFeatureLabelsWithSizes().SplitFeatureLabelSet( true);
          results_with_sizes = impl_a->GetResultCodeWithSizes().SplitFeatureLabelSet( true);
        }
        else if( all_rows_cols_equal)
        {
          if( number_rows != datasets.LastElement()->GetSize())
          {
            all_rows_cols_equal = false;
            BCL_MessageStd
            (
              "Datasets have different numbers of rows; difference statistics will not be computed"
            );
          }
          if( number_cols != datasets.LastElement()->GetFeatureSize())
          {
            BCL_Exit( "Datasets have different numbers of columns; cannot place all statistics on one table", -1);
          }
        }
      }

      // get labels for all the descriptors
      storage::Vector< util::ObjectDataLabel> labels( features_with_sizes.GetMemberLabels());
      labels.Append( results_with_sizes.GetMemberLabels());
      storage::Vector< size_t> sizes( features_with_sizes.GetPropertySizes());
      sizes.Append( results_with_sizes.GetPropertySizes());
      storage::Vector< std::string> rows;
      rows.AllocateMemory( labels.GetSize());
      storage::Vector< size_t>::const_iterator itr_sizes( sizes.Begin());
      for
      (
        util::ObjectDataLabel::const_iterator itr( labels.Begin()), itr_end( labels.End());
        itr != itr_end;
        ++itr, ++itr_sizes
      )
      {
        rows.PushBack( itr->ToString());
      }

      // create a table to hold all the statistics
      storage::Table< float> statistics( rows);

      // initially, construct the table as
      // Table Descriptor1 Descriptor2 Descriptor3 Descriptor4 ...
      // File1Ave
      // File2Ave
      // File1Std
      // File2Std
      // File1Max
      // File2Max
      // File1Min
      // File2Min
      // Std File1-File2, etc.
      // These will be transposed in the end because there are usually far more descriptors than files,
      // and it is easier to read/process excess rows than columns

      // track the rows that will be used for each file
      storage::Vector< util::SiPtrVector< storage::Row< float> > > rows_for_file( filenames.GetSize());
      // if there are multiple filenames, add each filename to each column heading
      for
      (
        storage::Vector< std::string>::const_iterator
          itr_col( default_headings.Begin()), itr_col_end( default_headings.End());
        itr_col != itr_col_end;
        ++itr_col
      )
      {
        size_t filenumber( 0);
        for
        (
          storage::Vector< std::string>::const_iterator itr( filenames.Begin()), itr_end( filenames.End());
          itr != itr_end;
          ++itr, ++filenumber
        )
        {
          const std::string descriptor_suffix( filenames.GetSize() > 1 ? "_" + *itr : std::string());
          const std::string col_name( *itr_col + descriptor_suffix);
          rows_for_file( filenumber).PushBack( statistics.InsertRow( col_name));
        }
      }

      // if all the datasets had the same # of rows and columns, compute standard deviations between each column
      if( all_rows_cols_equal)
      {
        // add columns for comparison of different features
        size_t filenumber( 0);
        for
        (
          storage::Vector< std::string>::const_iterator itr( filenames.Begin()), itr_end( filenames.End());
          itr != itr_end;
          ++itr, ++filenumber
        )
        {
          for
          (
            storage::Vector< std::string>::const_iterator itr_cmp( filenames.Begin());
            itr_cmp != itr;
            ++itr_cmp
          )
          {
            const std::string col_name_cmp( "Std_" + *itr + "-" + *itr_cmp);
            rows_for_file( filenumber).PushBack( statistics.InsertRow( col_name_cmp));
          }
        }
      }

      // for each file, iterate over datasets and rows in the statistics table to populate the statistics table
      storage::Vector< util::SiPtrVector< storage::Row< float> > >::iterator itr_rows( rows_for_file.Begin());
      size_t file_count( 0);

      // split all the features into 1 feature per row
      for
      (
        util::ShPtrList< descriptor::Dataset>::const_iterator
          itr_dataset( datasets.Begin()), itr_dataset_end( datasets.End());
        itr_dataset != itr_dataset_end;
        ++itr_dataset, ++itr_rows, ++file_count
      )
      {
        // get the dataset
        const descriptor::Dataset &dataset_a( **itr_dataset);

        // compute primary statistics on the dataset
        model::DataSetStatistics stats_a;
        stats_a.AddDataSet( dataset_a);

        // make a vector to hold references to all the statistical vectors
        storage::Vector< linal::Vector< float> > row_vects;
        row_vects.PushBack( stats_a.GetAveStdFeatures().GetAverage());
        row_vects.PushBack( stats_a.GetAveStdResults().GetAverage());
        row_vects.PushBack( stats_a.GetAveStdFeatures().GetStandardDeviation());
        row_vects.PushBack( stats_a.GetAveStdResults().GetStandardDeviation());
        row_vects.PushBack( stats_a.GetMinMaxFeatures().GetMin());
        row_vects.PushBack( stats_a.GetMinMaxResults().GetMin());
        row_vects.PushBack( stats_a.GetMinMaxFeatures().GetMax());
        row_vects.PushBack( stats_a.GetMinMaxResults().GetMax());

        // create statistics for all pairs of datasets if they all have the same # of rows and columns
        storage::List< model::DataSetStatistics> cmp_stats_list;

        if( all_rows_cols_equal && file_count)
        {
          for
          (
            util::ShPtrList< descriptor::Dataset>::const_iterator
              itr_dataset_b( datasets.Begin());
            itr_dataset_b != itr_dataset;
            ++itr_dataset_b
          )
          {
            cmp_stats_list.PushBack();
            cmp_stats_list.LastElement().AddDifference( **itr_dataset, **itr_dataset_b);
            row_vects.PushBack( cmp_stats_list.LastElement().GetAveStdFeatures().GetStandardDeviation());
            row_vects.PushBack( cmp_stats_list.LastElement().GetAveStdResults().GetStandardDeviation());
          }
        }
        BCL_Assert
        (
          2 * itr_rows->GetSize() == row_vects.GetSize(),
          "Unequal ref and index size " + util::Format()( itr_rows->GetSize())
          + " " + util::Format()( row_vects.GetSize())
        );

        // populate the statistics table with the statistical vectors
        storage::Vector< linal::Vector< float> >::const_iterator itr_stats( row_vects.Begin());
        for
        (
          util::SiPtrVector< storage::Row< float> >::iterator itr_row_id( itr_rows->Begin()), itr_row_id_end( itr_rows->End());
          itr_row_id != itr_row_id_end;
          ++itr_row_id, ++itr_stats
        )
        {
          storage::Vector< linal::Vector< float> >::const_iterator itr_stats_features( itr_stats);
          storage::Vector< linal::Vector< float> >::const_iterator itr_stats_results( ++itr_stats);
          // copy the features, store the iterator to the next element (which is where the result goes)
          storage::Vector< float>::iterator itr_result
          (
            std::copy( itr_stats_features->Begin(), itr_stats_features->End(), ( *itr_row_id)->Begin())
          );
          std::copy( itr_stats_results->Begin(), itr_stats_results->End(), itr_result);
        }
      }

      // transpose the table so that descriptors are in each row
      storage::Table< float> statistics_transposed( statistics.GetTransposedTable());
      statistics_transposed.WriteFormatted( OSTREAM, util::Format().FFP( 3), "Statistics");
    }

    const ApplicationType DescriptorGenerateDataset::GenerateDataset_Instance
    (
      GetAppGroups().AddAppToGroup( new DescriptorGenerateDataset(), GetAppGroups().e_Descriptor)
    );

  } // namespace app
} // namespace bcl
