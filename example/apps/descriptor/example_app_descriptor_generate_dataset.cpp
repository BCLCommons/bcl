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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "descriptor/bcl_app_descriptor_generate_dataset.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "model/bcl_model_retrieve_data_set_base.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_descriptor_generate_dataset.cpp
  //!
  //! @author mendenjl
  //! @date May 24, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppDescriptorGenerateDataset :
    public ExampleInterface
  {
  public:

    ExampleAppDescriptorGenerateDataset *Clone() const
    {
      return new ExampleAppDescriptorGenerateDataset( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      ApplicationExampleHelper generate_dataset_helper( app::DescriptorGenerateDataset::GenerateDataset_Instance);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      const std::string generated_891_descriptors_file_filename
      (
        AddExampleOutputPathToFilename( model::GetNamespaceIdentifier(), "generated_code_891_file.bin")
      );

      // check that the example input path exists.  If it does not, register an example check failure
      if( !io::DirectoryEntry( generate_dataset_helper.GetThisApplicationExampleInputPath()).DoesExist())
      {
        BCL_ExampleIndirectAssert
        (
          io::DirectoryEntry( generate_dataset_helper.GetThisApplicationExampleInputPath()).DoesExist(),
          true,
          "application example directory for this example, controlled by -application_example_path ("
          + generate_dataset_helper.GetThisApplicationExampleInputPath()
          + ") should exist to allow testing of this example"
        );
      }
      // file containing the correctly generated dataset in binary format
      const std::string correct_891_descriptors_filename_file
      (
        generate_dataset_helper.GetThisApplicationExampleInputPath() + "generated_code_from_file.bin.correct"
      );

      // file containing the correctly generated dataset in binary format for windows targets, which round somewhat differently
      const std::string correct_891_descriptors_filename_file_win
      (
        generate_dataset_helper.GetThisApplicationExampleInputPath() + "generated_code_from_file.bin.correct.win"
      );
      const std::string correct_891_descriptors_filename_file_mac
      (
        generate_dataset_helper.GetThisApplicationExampleInputPath() + "generated_code_from_file.bin.correct.mac"
      );

      // file containing descriptors to generate for each molecule
      const std::string filename_descriptor_labels( AddExampleInputPathToFilename( e_Model, "code_features.object"));

      // file containing descriptors to generate for each molecule
      const std::string filename_result_labels_db( AddExampleInputPathToFilename( e_Model, "code_result_db.object"));
      const std::string filename_result_labels_file( AddExampleInputPathToFilename( e_Model, "code_result_file.object"));

      // create a command line that will generate descriptors from a stored sdf
      generate_dataset_helper.ResetFlagsAndParameters();

      generate_dataset_helper.SetFlag( "feature_labels", filename_descriptor_labels);
      generate_dataset_helper.SetFlag( "result_labels", filename_result_labels_db);

      // Test the command with an sdf file
      const std::string aid_filename_actives
      (
        generate_dataset_helper.GetThisApplicationExampleInputPath() + "AID891_actives.sdf"
      );
      generate_dataset_helper.SetFlag( "source", "SdfFile( filename=" + aid_filename_actives + ")");
      generate_dataset_helper.SetFlag( "result_labels", filename_result_labels_file);

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( generate_dataset_helper.CheckCommandString( true), true);

      generate_dataset_helper.SetFlag( "output", generated_891_descriptors_file_filename);

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( generate_dataset_helper.RunCommand(), 0))
      {
        // if the application ran successfully, check that binary files match
        // check that the output is correct
        // open a binary stream to the correct file
        if
        (
          BCL_ExampleIndirectCheck
          (
            io::File::BinaryFilesMatch( generated_891_descriptors_file_filename, correct_891_descriptors_filename_file)
            || io::File::BinaryFilesMatch( generated_891_descriptors_file_filename, correct_891_descriptors_filename_file_win)
            || io::File::BinaryFilesMatch( generated_891_descriptors_file_filename, correct_891_descriptors_filename_file_mac),
            true,
            "Code generation produces correct binary file when molecules are loaded from an sdf file"
          )
        )
        {
          // remove the generated file since it is too large to commit to the repository and it is no longer needed
          remove( generated_891_descriptors_file_filename.c_str());
        }
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end Exampleexample_name

  const ExampleClass::EnumType ExampleAppDescriptorGenerateDataset::s_Instance
  (
    GetExamples().AddEnum( ExampleAppDescriptorGenerateDataset())
  );

} // namespace bcl
