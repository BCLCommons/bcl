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
#include "molecule/bcl_app_molecule_split.h"

// other bcl headers - sorted alphabetically
#include "example_application_example_helper.h"
#include "chemistry/bcl_chemistry.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example bcl_app_molecule_split.cpp
  //! @brief this example tests the implementation of the application MoleculeSplit, which fragments a given
  //! ensemble of molecules
  //!
  //! @see @link example_app_molecule_split.cpp @endlink
  //! @author geanesar
  //! @date Apr 18, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleAppMoleculeSplit :
    public ExampleInterface
  {

  //////////
  // data //
  //////////

  public:

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief clone function
    //! @return pointer to a new ExampleAppLoopTemplate
    ExampleAppMoleculeSplit *Clone() const
    {
      return new ExampleAppMoleculeSplit( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! @detail this is performing the execution of the example
    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      ApplicationExampleHelper molecule_split_helper( app::MoleculeSplit::MoleculeSplit_Instance);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // IO files
      const std::string input_ensemble( AddExampleInputPathToFilename( e_Chemistry, "mGluR5_five_actives.sdf"));
      const std::string output_basename( "mGluR5_scaffolds.sdf");
      const std::string output_filename( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), output_basename));

      // Set flags
      molecule_split_helper.ResetFlagsAndParameters();
      molecule_split_helper.SetFlag( "input_filenames", input_ensemble);
      molecule_split_helper.SetFlag( "output", output_filename);
      molecule_split_helper.SetFlag( "implementation", "Scaffolds");

      BCL_ExampleCheck( molecule_split_helper.CheckCommandString( true), true);

      // Run the command that generates scaffolds for molecules
      if( BCL_ExampleCheck( molecule_split_helper.RunCommand(), 0))
      {
        std::string correct_output( output_filename + ".correct");

        // check if the generated file is right; if so remove it
        if
        (
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatch( output_filename, correct_output),
            true,
            this->GetClassIdentifier() + " correct output"
          )
        )
        {
          remove( output_filename.c_str());
        }
        else
        {
          BCL_MessageStd( "Test file \"" + output_filename + "\" did not match expected output in \"" + correct_output + "\"");
        }
      }

      return 0;
    }

  }; // class ExampleAppMoleculeSplit

  //! single instance of this class
  const ExampleClass::EnumType ExampleAppMoleculeSplit::s_Instance
  (
     GetExamples().AddEnum( ExampleAppMoleculeSplit())
  );

} // namespace bcl
