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
#include "release/bcl_app_create_sse_pool.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_create_sse_pool.cpp
  //!
  //! @author mendenjl, putnamdk
  //! @date Jan 14, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAppCreateSSEPool :
    public ExampleInterface
  {
  public:

    ExampleAppCreateSSEPool *Clone() const
    {
      return new ExampleAppCreateSSEPool( *this);
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

    ////////////////////
    // initialization //
    ////////////////////

      // create a helper to run this application
      ApplicationExampleHelper create_sse_pool( app::ApplicationType( "CreateSSEPool"));

      // check that flags are needed
      BCL_ExampleCheck( create_sse_pool.CheckCommandString( false), false);

    ///////////
    // files //
    ///////////

      // create filenames for any input/output files used/generated by the test
      const std::string output_file_highest_three_method
      (
        AddExampleOutputPathToFilename
        (
          sspred::GetNamespaceIdentifier(),
          "1IE9A.SSPredHighest_PSIPRED_JUFO_PROFphd.pool"
        )
      );
      const std::string output_file_highest_three_method_correct( output_file_highest_three_method + ".correct");

      // try to create a mini-pool from jufo on 1IE9
      create_sse_pool.SetFlag( "prefix", AddExampleInputPathToFilename( e_Biology, "1IE9"));
      create_sse_pool.SetFlag( "ssmethods", sspred::GetMethods().e_JUFO);
      create_sse_pool.AddParameterToFlag( "ssmethods", sspred::GetMethods().e_PSIPRED);
      create_sse_pool.AddParameterToFlag( "ssmethods", sspred::GetMethods().e_PROFphd);
      create_sse_pool.SetFlag( "chop_sses");
      create_sse_pool.SetFlag( "output_prefix" , AddExampleOutputPathToFilename( sspred::GetNamespaceIdentifier(), ""));
      create_sse_pool.SetFlag( "pool_min_sse_lengths", "10"); // min helix length
      create_sse_pool.AddParameterToFlag( "pool_min_sse_lengths", "4"); // min coil length
      create_sse_pool.SetFlag( "sse_threshold", "0.4"); // helix threshold value
      create_sse_pool.AddParameterToFlag( "sse_threshold", "0.2"); // coil threshold value
      create_sse_pool.SetFlag( "factory", "SSPredHighest"); // arbitrary type of factory

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( create_sse_pool.RunCommand(), 0))
      {
        // check that the file is identical to the correctly generated one
        BCL_ExampleIndirectCheck
        (
          io::File::FilesMatch
          (
            output_file_highest_three_method,
            output_file_highest_three_method_correct
          ),

          true,
          "Validation of last command"
        );
      }

      // try the same command with just jufo, and changing various other flags
      create_sse_pool.SetFlag( "ssmethods", sspred::GetMethods().e_JUFO);
      create_sse_pool.SetFlag( "factory", "SSPredThreshold");
      create_sse_pool.UnsetFlag( "chop_sses");

      const std::string output_file_threshold_jufo
      (
        AddExampleOutputPathToFilename
        (
          sspred::GetNamespaceIdentifier(),
          "1IE9A.SSPredThreshold_JUFO.pool"
        )
      );
      const std::string output_file_threshold_jufo_correct( output_file_threshold_jufo + ".correct");

      // run a valid set of flags, check that the return status is 0
      if( BCL_ExampleCheck( create_sse_pool.RunCommand(), 0))
      {
        // check that the file is identical to the correctly generated one
        BCL_ExampleIndirectCheck
        (
          io::File::FilesMatch( output_file_threshold_jufo, output_file_threshold_jufo_correct),
          true,
          "Validation of last command"
        );
      }

      const std::string filename_1IE9_pool
      (
        AddExampleOutputPathToFilename( sspred::GetNamespaceIdentifier(), "1IE9.pool")
      );
      const std::string filename_1IE9_pool_correct( filename_1IE9_pool + ".correct");

      BCL_MessageStd( "filename is: " + filename_1IE9_pool_correct);

      // now try creating the pool from the pdb file
      create_sse_pool.ResetFlagsAndParameters();
      create_sse_pool.SetFlag( "pdb", AddExampleInputPathToFilename( e_Biology, "1IE9.pdb"));
      create_sse_pool.AddParameterToFlag( "pdb", filename_1IE9_pool);

      // Added to remove nightly build crashes
      create_sse_pool.SetFlag( "prefix", AddExampleInputPathToFilename( e_Biology, "1IE9"));
      create_sse_pool.SetFlag( "ssmethods", sspred::GetMethods().e_JUFO);

      if( BCL_ExampleCheck( create_sse_pool.RunCommand(), 0))
      {
        // check that the file is identical to the correctly generated one
        BCL_ExampleIndirectCheck
        (
          io::File::FilesMatch( filename_1IE9_pool, filename_1IE9_pool_correct),
          true,
          "creating the pool from the pdb file"
        );
      }

      BCL_MessageStd( "generated filename is: " + filename_1IE9_pool);

    ///////////////////////
    // integration tests //
    ///////////////////////

      // set flags; check valid and invalid command lines; check results

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAppCreateSSEPool

  const ExampleClass::EnumType ExampleAppCreateSSEPool::s_Instance
  (
    GetExamples().AddEnum( ExampleAppCreateSSEPool())
  );

} // namespace bcl
