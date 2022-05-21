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
#include "descriptor/bcl_app_descriptor_score_dataset.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "descriptor/bcl_descriptor.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_descriptor_refine_by_score.cpp
  //!
  //! @author mendenjl
  //! @date Mar 19, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppDescriptorScoreDataset :
    public ExampleInterface
  {
  public:

    ExampleAppDescriptorScoreDataset *Clone() const
    {
      return new ExampleAppDescriptorScoreDataset( *this);
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

      ApplicationExampleHelper descriptor_score_helper( app::DescriptorScoreDataset::s_Instance);

    /////////////////
    // data access //
    /////////////////

      // try running the command with no options, which should fail
      BCL_ExampleCheck( descriptor_score_helper.CheckCommandString( false), false);

      const std::string input_bin_file
      (
        descriptor_score_helper.GetThisApplicationExampleInputPath() + "aid891.bin"
      );

      // check that the example input path exists.  If it does not, register an example check failure
      if( !io::DirectoryEntry( descriptor_score_helper.GetThisApplicationExampleInputPath()).DoesExist())
      {
        BCL_ExampleIndirectAssert
        (
          io::DirectoryEntry( descriptor_score_helper.GetThisApplicationExampleInputPath()).DoesExist(),
          true,
          "application example directory for this example, controlled by -application_example_path ("
          + descriptor_score_helper.GetThisApplicationExampleInputPath()
          + ") should exist to allow testing of this example"
        );
      }

      // add the terse flag. this is used to prevent the output of the sorted scores, the precise order of which
      // varies across platforms
      descriptor_score_helper.SetFlag( "terse");

      // get the path to the bin file
      const std::string output_ig_scores_file
      (
        AddExampleOutputPathToFilename( descriptor::GetNamespaceIdentifier(), "aid891_infogain.scores")
      );
      const std::string output_f_scores_file
      (
        AddExampleOutputPathToFilename( descriptor::GetNamespaceIdentifier(), "aid891_fscore.scores")
      );

      // add the flag to access this score file
      descriptor_score_helper.SetFlag( "score", "Partition(partitioner=InformationGain,cutoff=3.5)");

      // add the source flag
      descriptor_score_helper.SetFlag( "source", "Subset(filename=" + input_bin_file + ")");

      // add the output flag
      descriptor_score_helper.SetFlag( "output", output_ig_scores_file);

      // test that the command can be run
      if( BCL_ExampleCheck( descriptor_score_helper.CheckCommandString( true), true))
      {
        // run the command; expect 0 errors
        if( BCL_ExampleCheck( descriptor_score_helper.RunCommand(), 0))
        {
          // command ran, check files
          if
          (
            BCL_ExampleIndirectCheck
            (
              io::File::FilesMatchWithinAbsoluteTolerance( output_ig_scores_file, output_ig_scores_file + ".correct", 0.1),
              true,
              descriptor_score_helper.GetCurrentCommandLine()
            )
          )
          {
            // file was fine, remove the newly-generated output file
            remove( output_ig_scores_file.c_str());
          }
        }
      }

      // add the output flag
      descriptor_score_helper.SetFlag( "output", output_f_scores_file);

      // add the flag to set the score type
      descriptor_score_helper.SetFlag( "score", "FScore(cutoff=3.5)");

      // test that the command can be run
      if( BCL_ExampleCheck( descriptor_score_helper.CheckCommandString( true), true))
      {
        // run the command; expect 0 errors
        if( BCL_ExampleCheck( descriptor_score_helper.RunCommand(), 0))
        {
          // command ran, check files
          if
          (
            BCL_ExampleIndirectCheck
            (
              io::File::FilesMatchWithinAbsoluteTolerance( output_f_scores_file, output_f_scores_file + ".correct", 0.1),
              true,
              descriptor_score_helper.GetCurrentCommandLine()
            )
          )
          {
            // file was fine, remove the newly-generated output file
            remove( output_f_scores_file.c_str());
          }
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

  }; //end ExampleAppDescriptorScoreDataset

  const ExampleClass::EnumType ExampleAppDescriptorScoreDataset::s_Instance
  (
    GetExamples().AddEnum( ExampleAppDescriptorScoreDataset())
  );

} // namespace bcl
