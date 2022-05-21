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
#include "bcl_app_descriptor_convert_code_object_file.h"

// includes from the bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_range.h"
#include "model/bcl_model_retrieve_data_set_base.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string DescriptorConvertCodeObjectFile::GetDescription() const
    {
      return "converts a code object file "
          "from an CSV-style formating version into the new version given a complete "
          "set of CSV-formatted descriptors and standard BCL code-object file descriptors.";
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &DescriptorConvertCodeObjectFile::GetReadMe() const
    {
      static std::string s_read_me =
          "DescriptorConvertCodeObjectFile is an application that converts a code object file "
          "from an CSV-style formating version into the new version given a complete "
          "set of CSV-formatted descriptors and standard BCL code-object file descriptors.";
      return s_read_me;
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> DescriptorConvertCodeObjectFile::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // code object file of interest
      sp_cmd->AddParameter( m_CodeObjectFile);
      // code object file of new converted code object file
      sp_cmd->AddParameter( m_ConvertedCodeObjectFile);
      // filenames for mapping one code object file onto another
      sp_cmd->AddFlag( m_CodeMapping);
      // flag for source to use for splitting a code object file
      sp_cmd->AddFlag( m_Split);
      // flag for merging code object files
      sp_cmd->AddFlag( m_Merge);
      // flag for updating code object files
      sp_cmd->AddFlag( m_Update);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    } // InitializeCommand

    //! @brief the Main function
    //! @return error code - 0 for success
    int DescriptorConvertCodeObjectFile::Main() const
    {
      BCL_Assert
      (
        int( m_CodeMapping->GetFlag()) + int( m_Split->GetFlag()) + int( m_Merge->GetFlag()) + int( m_Update->GetFlag())
        == int( 1),
        "Exactly one of -split, -map, -merge, or -update must be given!"
      );

      io::IFStream input;
      io::File::MustOpenIFStream( input, m_CodeObjectFile->GetValue());
      util::ObjectDataLabel code_object( input);
      io::File::CloseClearFStream( input);
      io::OFStream out;
      io::File::MustOpenOFStream( out, m_ConvertedCodeObjectFile->GetValue());

      if( m_CodeMapping->GetFlag())
      {
        const storage::Vector< std::string> filenames( m_CodeMapping->GetStringList());
        BCL_Assert( filenames.GetSize() == size_t( 2), "Exactly two filenames must be given for map flag");

        io::File::MustOpenIFStream( input, filenames( 0));
        util::ObjectDataLabel complete_old( input);
        io::File::CloseClearFStream( input);

        io::File::MustOpenIFStream( input, filenames( 1));
        util::ObjectDataLabel complete_new( input);
        io::File::CloseClearFStream( input);

        BCL_Assert
        (
          complete_old.GetNumberArguments() == complete_old.GetNumberArguments(),
          "Number of properties in complete code object files differs!"
          " old: " + util::Format()( complete_old.GetNumberArguments()) +
          " new: " + util::Format()( complete_new.GetNumberArguments())
        );

        BCL_MessageStd
        (
          "Number of properties in complete code object files \n"
          "old: " + util::Format()( complete_old.GetNumberArguments()) +
          " new: " + util::Format()( complete_new.GetNumberArguments())
        );

        storage::Map< util::ObjectDataLabel, util::ObjectDataLabel> mapping;

        util::ObjectDataLabel::const_iterator itr_old( complete_old.Begin());

        for
        (
          util::ObjectDataLabel::const_iterator
            itr_new( complete_new.Begin()), itr_new_end( complete_new.End());
          itr_new != itr_new_end;
          ++itr_new, ++itr_old
        )
        {
          mapping[ *itr_old] = *itr_new;
        }

        storage::Vector< util::ObjectDataLabel> new_labels;
        new_labels.AllocateMemory( code_object.GetNumberArguments());

        for
        (
          util::ObjectDataLabel::const_iterator
            itr( code_object.Begin()), itr_end( code_object.End());
          itr != itr_end;
          ++itr
        )
        {
          new_labels.PushBack( mapping[ *itr]);
        }

        util::ObjectDataLabel out_file( code_object.GetName(), code_object.GetValue(), new_labels);
        out << out_file.ToNamedString( 100, 0, 1) << "\n";
      }
      else if( m_Split->GetFlag())
      {
        util::Implementation< model::RetrieveDataSetBase> retriever( m_Split->GetFirstParameter()->GetValue());
        retriever->SelectFeatures( code_object);
        model::FeatureLabelSet features( retriever->GetFeatureLabelsWithSizes().SplitFeatureLabelSet());

        util::ObjectDataLabel out_file( code_object.GetName(), code_object.GetValue(), features.GetMemberLabels());
        io::OFStream out;

        io::File::MustOpenOFStream( out, m_ConvertedCodeObjectFile->GetValue());
        out << out_file.ToNamedString( 100, 0, 1) << "\n";
      }
      else if( m_Merge->GetFlag())
      {
        io::File::MustOpenIFStream( input, m_Merge->GetFirstParameter()->GetValue());
        util::ObjectDataLabel code_object_to_merge( input);
        io::File::CloseClearFStream( input);

        out << model::FeatureLabelSet::MergeConsideringPartials( code_object, code_object_to_merge).ToNamedString( 360)
            << "\n";
      }
      else if( m_Update->GetFlag())
      {
        util::Implementation< model::RetrieveDataSetBase> retriever( m_Update->GetFirstParameter()->GetValue());
        retriever->SelectFeatures( code_object);
        model::FeatureLabelSet features( retriever->GetFeatureLabelsWithSizes());

        util::ObjectDataLabel out_file( code_object.GetName(), code_object.GetValue(), features.GetMemberLabels());
        out << out_file.ToNamedString( 100, 0, 1) << "\n";
      }
      io::File::CloseClearFStream( out);
      return 0;
    } // Main

    //! @brief standard constructor
    DescriptorConvertCodeObjectFile::DescriptorConvertCodeObjectFile() :
      m_CodeObjectFile
      (
        new command::Parameter
        (
          "code_object_file",
          "file containing the feature labels you want to convert",
          command::ParameterCheckFileExistence()
        )
      ),
      m_ConvertedCodeObjectFile
      (
        new command::Parameter
        (
          "output_file",
          "output file name of convert code object file",
          "converted.txt"
        )
      ),
      m_CodeMapping
      (
        new command::FlagStatic
        (
          "map",
          "map object labels from an old code object file to corresponding labels in a new file",
          util::ShPtrVector< command::ParameterInterface>::Create
          (
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "old_code_object_file",
                "code object file containing all available descriptors in old format",
                command::ParameterCheckFileExistence(),
                ""
              )
            ),
            util::ShPtr< command::ParameterInterface>
            (
              new command::Parameter
              (
                "new_code_object_file",
                "code object file containing all available descriptors in new format",
                command::ParameterCheckFileExistence(),
                ""
              )
            )
          )
        )
      ),
      m_Split
      (
        new command::FlagStatic
        (
          "split",
          "split descriptors for multiple columns in each code object file into individual descriptors, one per line",
          command::Parameter
          (
            "dataset_retriever",
            "dataset retriever for obtaining the sizes of each label",
            command::ParameterCheckSerializable( util::Implementation< model::RetrieveDataSetBase>()),
            ""
          )
        )
      ),
      m_Merge
      (
        new command::FlagStatic
        (
          "merge",
          "merge descriptors, including partials",
          command::Parameter
          (
            "file_to_merge",
            "file name of file to merge",
            command::ParameterCheckFileExistence(),
            ""
          )
        )
      ),
      m_Update
      (
        new command::FlagStatic
        (
          "update",
          "update descriptors using the given dataset retriever; allows descriptor files that were not used in dataset "
          "generation to be used for dataset retrieval",
          command::Parameter
          (
            "dataset_retriever",
            "dataset retriever for obtaining the sizes of each label",
            command::ParameterCheckSerializable( util::Implementation< model::RetrieveDataSetBase>()),
            ""
          )
        )
      )
    {
    }

    const ApplicationType DescriptorConvertCodeObjectFile::DescriptorConvertCodeObjectFile_Instance
    (
      GetAppGroups().AddAppToGroup( new DescriptorConvertCodeObjectFile(), GetAppGroups().e_Descriptor)
    );

  } // namespace app
} // namespace bcl
