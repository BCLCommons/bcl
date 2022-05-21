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
#include "bcl_app_descriptor_analyze.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_range.h"
#include "model/bcl_model_retrieve_data_set_base.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {
    //! @brief Clone function
    //! @return pointer to new FoldProtein
    DescriptorAnalyze *DescriptorAnalyze::Clone() const
    {
      return new DescriptorAnalyze( *this);
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string DescriptorAnalyze::GetDescription() const
    {
      return "Analyze the tanimoto overlap between descriptor sets";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &DescriptorAnalyze::GetReadMe() const
    {
      static const std::string s_readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of the BCL application descriptor:Analyze, its terms of use, "
        "appropriate citation, installation procedures, descriptor:Analyze execution, and technical support\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. What is descriptor:Analyze?\n"
        "descriptor:Analyze is a C++ based application, created by Vanderbilt University's Meiler Laboratory, "
        "which is part a larger library of applications called BCL::Commons.  descriptor:Analyze is used to "
        "analyze files containing descriptor objects.  Currently, it is limited to analyzing the % of features various\n"
        "descriptor files have in common"
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING descriptor:Analyze.\n"
        "When using descriptor:Analyze in a publication, please cite the following publications describing the "
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
        "VI. RUNNING descriptor:Analyze.\n"
        "Running descriptor:Analyze consists of three main steps.\n"
        "\n"
        "1) Specify your data source(s):\n"
        "This is what you will pass to the -retriever flag. "
        "There are many dataset retrievers in the bcl, so type bcl.exe descriptor:Analyze -help shows a full list. "
        "Typically, the primary data source comes from files containing one of the following:\n"
        "1. molecules (in sdf format) using SdfFile\n"
        "2. sequences (in fasta format) using SequenceDirectory\n"
        "3. proteins (in pdb format) using ProteinDirectory retriever\n"
        "4. another bcl application that writes out native bcl vector formats using the File retriever.\n"
        "5. bcl-native dataset format known as the .bin file (using the Subset retriever).\n"
        "Whichever dataset retrievers you use, read the help for them using "
        "e.g. bcl.exe descriptor:Analyze -retriever \"SdfFile(help)\" and use it to specify the internal parameters of the retriever"
        "Note that a retriever is not strictly necessary; if used without a retriever, this application will look at the % of "
        "overlapping descriptor groups, instead of the percentage of overlapping descriptor columns.  This is appropriate only "
        "when the descriptor file contains no groups of size > 1 or otherwise contains no groups that were split"
        " differently (using Partial) across different files.  Generally, it is best to use the retriever flag unless you "
        "really understand what you're doing"
        "\n"
        "2) Choose which descriptors to construct:\n"
        "\nA. Background\n"
        "For bcl-native vector format, the only descriptors available are the column indices (0-indexed)\n"
        "For bin files, the provided descriptors must be a subset of those used when the file was originally created.\n"
        "To see a list of all available descriptors that can be created for your primary source (other than .bin files "
        "and vector files), type\n"
        "bcl.exe descriptor:Analyze -retriever _PrimarySource_ -code_object_files_horizontal help\n"
        "e.g.\n"
        "bcl.exe descriptor:Analyze -retriever \"ProteinDirectory(.)\" -code_object_files_horizontal help\n"
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
        "3) Run descriptor:Analyze:\n"
        "At a command prompt, navigate to the location of your BCL executable program. The syntax for running "
        "the application looks like the following:\n"
        "\n"
        "bcl.exe descriptor:Analyze -retriever <source-specifier> -code_object_files_horizontal file1 file2 ... "
        "-code_object_files_vertical file3 file4 ... -output_matrix_filename overlap.txt "
        "\nADDITIONAL INFORMATION\n"
        "\n"
        "Helpful information concerning syntax and flags can be obtained by typing\n"
        "bcl.exe descriptor:Analyze  -help\n"
        "\n"
        "For this read me, type\n"
        "bcl.exe descriptor:Analyze -readme\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF descriptor:Analyze.\n"
        "descriptor:Analyze is under ongoing further development. For current research please refer to "
        "www.meilerlab.org and navigate to research\n"
        + DefaultSectionSeparator()
      );
      return s_readme;
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> DescriptorAnalyze::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // code object file of interest
      sp_cmd->AddFlag( m_CodeObjectFilesRow);

      // code object file of interest
      sp_cmd->AddFlag( m_CodeObjectFilesColumn);

      // filename of the output matrix
      sp_cmd->AddFlag( m_OutputMatrixFileName);

      // dataset retriever; used to determine the size of each feature
      sp_cmd->AddFlag( m_DatasetRetriever);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    } // InitializeCommand

    //! @brief the Main function
    //! @return error code - 0 for success
    int DescriptorAnalyze::Main() const
    {
      // container with all code object files
      storage::Vector< model::FeatureLabelSet> feature_label_sets_rows, feature_label_sets_cols;

      // container with sizes of each label
      storage::Vector< size_t> code_object_sizes;

      // create a retriever for its SelectFeatures() method
      util::Implementation< model::RetrieveDataSetBase> retriever;
      if( m_DatasetRetriever->GetFlag())
      {
        retriever =
          util::Implementation< model::RetrieveDataSetBase>( m_DatasetRetriever->GetFirstParameter()->GetValue());
      }
      else
      {
        BCL_MessageCrt
        (
          "Warning: No descriptor retriever specified; will return percentage of overlapping descriptor groups, assuming "
          "that the descriptor groups are indivisible.  If the files contains Partial descriptors, "
          "then results will not be reliable without the -retriever flag"
        );
      }

      // iterate over all code object files per row
      for
      (
        util::ShPtrVector< command::ParameterInterface>::const_iterator
          itr( m_CodeObjectFilesRow->GetParameterList().Begin()),
          itr_end( m_CodeObjectFilesRow->GetParameterList().End());
        itr != itr_end;
        ++itr
      )
      {
        io::IFStream input;
        io::File::MustOpenIFStream( input, ( *itr)->GetValue());

        // convert code object file into its Partial components of the retriever was given
        if( retriever.IsDefined())
        {
          retriever->SelectFeatures( util::ObjectDataLabel( input));
          feature_label_sets_rows.PushBack( retriever->GetFeatureLabelsWithSizes().SplitFeatureLabelSet());
        }
        else
        {
          // retriever omitted, compare features directly, assuming they cannot be split
          util::ObjectDataLabel complete_label( input);
          model::FeatureLabelSet label_set
          (
            complete_label.GetValue(),
            complete_label.GetArguments(),
            storage::Vector< size_t>( complete_label.GetNumberArguments(), size_t( 1))
          );
          feature_label_sets_rows.PushBack( label_set);
        }
        io::File::CloseClearFStream( input);
      }

      // iterate over all code object files per row
      for
      (
        util::ShPtrVector< command::ParameterInterface>::const_iterator
          itr( m_CodeObjectFilesColumn->GetParameterList().Begin()),
          itr_end( m_CodeObjectFilesColumn->GetParameterList().End());
        itr != itr_end;
        ++itr
      )
      {
        io::IFStream input;
        io::File::MustOpenIFStream( input, ( *itr)->GetValue());

        // convert code object file into its Partial components of the retriever was given
        if( retriever.IsDefined())
        {
          retriever->SelectFeatures( util::ObjectDataLabel( input));
          feature_label_sets_cols.PushBack( retriever->GetFeatureLabelsWithSizes().SplitFeatureLabelSet());
        }
        else
        {
          // retriever omitted, compare features directly, assuming they cannot be split
          util::ObjectDataLabel complete_label( input);
          model::FeatureLabelSet label_set
          (
            complete_label.GetName(),
            complete_label.GetArguments(),
            storage::Vector< size_t>( complete_label.GetNumberArguments(), size_t( 1))
          );
          feature_label_sets_cols.PushBack( label_set);
        }

        io::File::CloseClearFStream( input);
      }

      linal::Matrix< float> overlap_matrix
      (
        m_CodeObjectFilesRow->GetStringList().GetSize(),
        m_CodeObjectFilesColumn->GetStringList().GetSize(),
        float( 0)
      );

      // iterate over all code object files per column
      size_t col_id( 0);
      for
      (
          storage::Vector< model::FeatureLabelSet>::const_iterator
          itr( feature_label_sets_cols.Begin()), itr_end( feature_label_sets_cols.End());
        itr != itr_end;
        ++itr, ++col_id
      )
      {
        size_t row_id( 0);
        for
        (
            storage::Vector< model::FeatureLabelSet>::const_iterator
            itr_compare( feature_label_sets_rows.Begin()), itr_end_compare( feature_label_sets_rows.End());
          itr_compare != itr_end_compare;
          ++itr_compare, ++row_id
        )
        {

          const model::FeatureLabelSet compare_feature_label( itr->GetCommonFeatures( *itr_compare));

          // compute overlap of Partial descriptors in each code object file
          const float percent( compare_feature_label.GetSize() / float( itr_compare->GetSize()));
          overlap_matrix( row_id, col_id) = percent;

          BCL_MessageDbg
          (
            "--feat label compare: " + util::Format()( itr_compare->GetSize())
            + " overlap: " + util::Format()( compare_feature_label.GetSize())
            + " (" + util::Format()( compare_feature_label.GetSize() / float( itr_compare->GetSize()))
            + ")"
          );
        }
      }

      // write overlap overlap matrix if flag is set
      if( m_OutputMatrixFileName->GetFlag())
      {
        io::OFStream output;
        io::File::MustOpenOFStream( output, m_OutputMatrixFileName->GetFirstParameter()->GetValue());
        output << overlap_matrix;
        io::File::CloseClearFStream( output);

        BCL_MessageStd
        (
          "overlap matrix written out to: " + m_OutputMatrixFileName->GetFirstParameter()->GetValue()
        );
      }
      else
      {
        BCL_MessageStd( "overlap matrix: " + util::Format()( overlap_matrix));
      }

      return 0;
    } // Main

    //! @brief standard constructor
    DescriptorAnalyze::DescriptorAnalyze() :
      m_CodeObjectFilesRow
      (
        new command::FlagDynamic
        (
          "code_object_files_horizontal",
          "code object files along the x-axis of the overlap matrix",
          storage::Vector< command::Parameter>
          (
            1,
            command::Parameter
            (
              "code_object_file_vertical",
              "code object files along the x-axis of the overlap matrix"
            )
          )
        )
      ),
      m_CodeObjectFilesColumn
      (
        new command::FlagDynamic
        (
          "code_object_files_vertical",
          "code object files along the y-axis of the overlap matrix",
          storage::Vector< command::Parameter>
          (
            1,
            command::Parameter
            (
              "code_object_file_vertical",
              "code object files along the y-axis of the overlap matrix"
            )
          )
        )
      ),
      m_OutputMatrixFileName
      (
        new command::FlagStatic
        (
          "output_matrix_filename",
          "",
          command::Parameter
          (
            "output_matrix_filename",
            "filename of output matrix",
            "overlap.txt"
          )
        )
      ),
      m_DatasetRetriever
      (
        new command::FlagStatic
        (
          "retriever",
          "dataset retriever; used to determine the size of the features. "
          "If omitted, features will not be split and the overlap values returned will be based on the overall"
          "# of overlapping descriptor groups, rather than the % of overlapping columns",
          command::Parameter
          (
            "retriever",
            "",
            command::ParameterCheckSerializable( util::Implementation< model::RetrieveDataSetBase>()),
            ""
          )
        )
      )
    {
    }

    const ApplicationType DescriptorAnalyze::DescriptorAnalyze_Instance
    (
      GetAppGroups().AddAppToGroup( new DescriptorAnalyze(), GetAppGroups().e_Descriptor)
    );

  } // namespace app
} // namespace bcl
