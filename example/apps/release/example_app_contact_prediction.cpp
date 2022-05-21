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
#include "release/bcl_app_contact_prediction.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "contact/bcl_contact.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_contact_prediction.cpp
  //!
  //! @author teixeipl, weinerbe
  //! @date Oct 15, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppContactPrediction :
    public ExampleInterface
  {
  public:

    ExampleAppContactPrediction *Clone() const
    {
      return new ExampleAppContactPrediction( *this);
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

      ApplicationExampleHelper contact_helper( app::ContactPrediction::ContactPrediction_Instance);

    ////////////////
    // operations //
    ////////////////

      // check that flags are needed
      BCL_ExampleCheck( contact_helper.CheckCommandString( false), false);

    ///////////////
    // operators //
    ///////////////

      // set input and output filenames/dirs
      const std::string input_fasta( AddExampleInputPathToFilename( e_Biology, "1ubiA.fasta"));
      const std::string input_pdb_list( AddExampleInputPathToFilename( e_Biology, "contact_pred_test_1ubi.ls"));
      const std::string output_dir( AddExampleOutputPathToFilename( contact::GetNamespaceIdentifier(), ""));

      // set fasta file as parameter
      contact_helper.AddParameter( input_fasta);

      // set output flag
      contact_helper.SetFlag( "output_path", output_dir);

      // check command line
      BCL_ExampleAssert( contact_helper.CheckCommandString( true), true);

      // contact files
      const std::string generated_output( output_dir + "1ubiA.contact");
      const std::string generated_pdb_output( output_dir + "1ubiA.contacts");
      const std::string ref_output( output_dir + "1ubiA.contact.ref");

      // run command using fasta file
      if( BCL_ExampleCheck( contact_helper.RunCommand(), 0))
      {
        // make sure the output file agrees with the correct version
        if
        (
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatch( generated_output, ref_output),
            true,
            "Comparison of output files generated from FASTA"
          )
        )
        {
          // remove the generated file
          remove( generated_output.c_str());
        }
      }

      // reset parameters
      contact_helper.ResetParameters();

      // set new PDB list parameter
      contact_helper.AddParameter( input_pdb_list);

      // set flags
      contact_helper.SetFlag( "pdb_list");
      contact_helper.SetFlag( "real_contacts");

      // check command line
      BCL_ExampleAssert( contact_helper.CheckCommandString( true), true);

      // run command using pdb list
      if( BCL_ExampleCheck( contact_helper.RunCommand(), 0))
      {
        // make sure the output file agrees with the correct version
        BCL_ExampleIndirectCheck
        (
          io::File::FilesMatch( generated_output, ref_output),
          true,
          "Comparison of output files generated from PDB"
        );

        // remove the generated files
        remove( generated_output.c_str());
        remove( generated_pdb_output.c_str());
      }

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAppContactPrediction

  const ExampleClass::EnumType ExampleAppContactPrediction::s_Instance
  (
    GetExamples().AddEnum( ExampleAppContactPrediction())
  );

} // namespace bcl
