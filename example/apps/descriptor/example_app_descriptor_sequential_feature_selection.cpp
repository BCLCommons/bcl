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
#include "descriptor/bcl_app_descriptor_sequential_feature_selection.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "io/bcl_io_directory_entry.h"
#include "model/bcl_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_descriptor_sequential_feature_selection.cpp
  //! @brief application example for application app::DescriptorSequentialFeatureSelection.
  //!
  //! @author butkiem1
  //! @date Feb 23, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppDescriptorSequentialFeatureSelection :
    public ExampleInterface
  {
  public:

    ExampleAppDescriptorSequentialFeatureSelection *Clone() const
    {
      return new ExampleAppDescriptorSequentialFeatureSelection( *this);
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

    ////////////////
    // unit tests //
    ////////////////

      // unit tests of application class functions, if applicable, go here

    ////////////////////
    // initialization //
    ////////////////////

      // create a helper to run this application
      ApplicationExampleHelper train_model_desc_sel
      (
        app::DescriptorSequentialFeatureSelection::DescriptorSequentialFeatureSelection_Instance
      );

      // specify output directory
      std::string output_directory_name( AddExampleOutputPathToFilename( model::GetNamespaceIdentifier(), "desc_select_ffs/"));

      // create a command line
      train_model_desc_sel.ResetFlagsAndParameters();

      // parameters
      train_model_desc_sel.AddParameter( AddExampleInputPathToFilename( e_Model, "code_features.object"));
      train_model_desc_sel.AddParameter( AddExampleOutputPathToFilename( model::GetNamespaceIdentifier(), "desc_select_ffs/round"));

      // flags
      train_model_desc_sel.SetFlag( "descriptor_selection_type", "FeatureForwardSelection");
      train_model_desc_sel.SetFlag( "storage_descriptor_selection", "File(directory=" + output_directory_name + ")");
      train_model_desc_sel.SetFlag( "get_initial_descriptor_set_for_this_round");
      train_model_desc_sel.SetFlag( "flag_round_number", "27");

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( train_model_desc_sel.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( train_model_desc_sel.RunCommand(), 0))
      {
        BCL_MessageStd( "command line executed successfully!");
      }

      // check if directory was created
      BCL_ExampleCheck( io::DirectoryEntry( output_directory_name).DoesExist(), true);

      // check number of generated files in output directory
      BCL_ExampleCheck( io::Directory( output_directory_name).ListEntries().GetSize(), size_t( 28));

      // clean up generated output directory
      BCL_ExampleCheck( io::Directory( output_directory_name).Remove( true), true);

      // check if directory was created
      BCL_ExampleCheck( io::DirectoryEntry( output_directory_name).DoesExist(), false);

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // testing use case where no initial descriptor file with all descriptor groups left for this round is generated
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // create a command line
      train_model_desc_sel.ResetFlagsAndParameters();

      // parameters
      train_model_desc_sel.AddParameter( AddExampleInputPathToFilename( e_Model, "code_features.object"));
      train_model_desc_sel.AddParameter( AddExampleOutputPathToFilename( model::GetNamespaceIdentifier(), "desc_select_ffs/round"));

      // flags
      train_model_desc_sel.SetFlag( "descriptor_selection_type", "FeatureForwardSelection");
      train_model_desc_sel.SetFlag( "storage_descriptor_selection", "File(directory=" + output_directory_name + ")");
      train_model_desc_sel.SetFlag( "flag_round_number", "27");

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( train_model_desc_sel.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( train_model_desc_sel.RunCommand(), 0))
      {
        BCL_MessageStd( "command line executed successfully!");
      }

      // check if directory was created
      BCL_ExampleCheck( io::DirectoryEntry( output_directory_name).DoesExist(), true);

      // check number of generated files in output directory
      BCL_ExampleCheck( io::Directory( output_directory_name).ListEntries().GetSize(), size_t( 27));

      // clean up generated output directory
      BCL_ExampleCheck( io::Directory( output_directory_name).Remove( true), true);

      // check if directory was created
      BCL_ExampleCheck( io::DirectoryEntry( output_directory_name).DoesExist(), false);

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // testing use case when last round is reached
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // create a command line
      train_model_desc_sel.ResetFlagsAndParameters();

      // parameters
      train_model_desc_sel.AddParameter( AddExampleInputPathToFilename( e_Model, "code_features.object"));
      train_model_desc_sel.AddParameter( AddExampleOutputPathToFilename( model::GetNamespaceIdentifier(), "desc_select_ffs/round"));

      // flags
      train_model_desc_sel.SetFlag( "descriptor_selection_type", "FeatureForwardSelection");
      train_model_desc_sel.SetFlag( "storage_descriptor_selection", "File(directory=" + output_directory_name + ")");
      train_model_desc_sel.SetFlag( "get_initial_descriptor_set_for_this_round");
      train_model_desc_sel.SetFlag( "flag_round_number", "0");

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( train_model_desc_sel.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( train_model_desc_sel.RunCommand(), 0))
      {
        BCL_MessageStd( "command line executed successfully!");
      }

      // check if directory was created
      BCL_ExampleCheck( io::DirectoryEntry( output_directory_name).DoesExist(), true);

      // check number of generated files in output directory
      BCL_ExampleCheck( io::Directory( output_directory_name).ListEntries().GetSize(), size_t( 1));

      // clean up generated output directory
      BCL_ExampleCheck( io::Directory( output_directory_name).Remove( true), true);

      // check if directory was created
      BCL_ExampleCheck( io::DirectoryEntry( output_directory_name).DoesExist(), false);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAppDescriptorSequentialFeatureSelection

  const ExampleClass::EnumType ExampleAppDescriptorSequentialFeatureSelection::s_Instance
  (
    GetExamples().AddEnum( ExampleAppDescriptorSequentialFeatureSelection())
  );

} // namespace bcl
