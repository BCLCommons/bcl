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
#include "internal/chemistry/bcl_app_conformer_generator.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "chemistry/bcl_chemistry.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_rotamer_library_interface.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example bcl_app_conformer_generator.cpp
  //! @brief this example tests the implementation of the application ConformerGenerator, which fragments a given
  //! ensemble of molecules
  //!
  //! @author kothiwsk
  //! @date Jan 28, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleAppConformerGenerator :
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
    ExampleAppConformerGenerator *Clone() const
    {
      return new ExampleAppConformerGenerator( *this);
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

      ApplicationExampleHelper conformer_generator_helper( app::ConformerGenerator::ConformerGenerator_Instance);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      std::string extensions( ".sdf");
      // determine whether file compression can be used
      if( io::GetStreamBufferClasses().GetCompressionFromExtension( "bz2").IsDefined())
      {
        extensions += ".bz2";
      }

      io::IFStream input;

      const std::string input_ensemble( AddExampleInputPathToFilename( e_Chemistry, "mGluR5_five_actives.sdf"));

      const std::string output_filename( AddExampleOutputPathToFilename( chemistry::GetNamespaceIdentifier(), "sample_conformations.sdf"));

      // create a command line that will generate descriptors from a stored sdf
      conformer_generator_helper.ResetFlagsAndParameters();

      conformer_generator_helper.SetFlag( "add_h");
      conformer_generator_helper.SetFlag( "ensemble_filenames", input_ensemble);
      conformer_generator_helper.SetFlag( "conformers_single_file", output_filename);
      conformer_generator_helper.SetFlag( "rnd_dihedral_mutate_weight", "0.01");
      conformer_generator_helper.SetFlag( "top_models", "5");
      conformer_generator_helper.SetFlag( "max_iterations", "10");
      conformer_generator_helper.SetFlag( "cluster");
      conformer_generator_helper.SetFlag
      (
        "conformation_comparer",
        storage::Vector< std::string>::Create( "SymmetryRMSD", "0.25")
      );

      // check a valid set of flags.  Since all later commands depend on this command succeeding, make it an assert
      BCL_ExampleCheck( conformer_generator_helper.CheckCommandString( true), true);

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( conformer_generator_helper.RunCommand(), 0))
      {
        // sample conformer output differs on mac and linux almost entirely due to differences in how std::random_shuffle
        // is implemented across libraries. Likewise we have to check two equally valid files
        std::string correct_filename( output_filename + ".correct");
        std::string correct_filename_mac( correct_filename + ".mac");
        if( io::File::FilesMatchWithinAbsoluteTolerance( output_filename, correct_filename, 0.1))
        {
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatchWithinAbsoluteTolerance( output_filename, correct_filename, 0.1),
            true,
            "AppConformerGenerator Linux-like std::random_shuffle"
          );
          remove( output_filename.c_str());
        }
        else if( io::File::FilesMatchWithinAbsoluteTolerance( output_filename, correct_filename_mac, 0.1))
        {
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatchWithinAbsoluteTolerance( output_filename, correct_filename_mac, 0.1),
            true,
            "AppConformerGenerator Mac/clang-like std::random_shuffle"
          );
          remove( output_filename.c_str());
        }
        else
        {
          BCL_ExampleIndirectCheck
          (
            io::File::FilesMatchWithinAbsoluteTolerance( output_filename, correct_filename, 0.1),
            true,
            "AppConformerGenerator Mac/Linux"
          );
        }
      }

      // check if the generated file is right

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
      }

  }; // class ExampleAppConformerGenerator

  //! single instance of this class
  const ExampleClass::EnumType ExampleAppConformerGenerator::s_Instance
  (
     GetExamples().AddEnum( ExampleAppConformerGenerator())
  );

} // namespace bcl
