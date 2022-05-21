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
#include "fold/bcl_fold_locator_missing_coordinates.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"

// external includes - sorted alphabetically

namespace bcl
{

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_locator_loop.cpp
  //! @brief this example tests the implementation of fold::LocatorMissingCoordinates
  //!
  //! @author fischea
  //! @date Mar 2, 2017
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleFoldLocatorMissingCoordinates :
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
    //! @return pointer to a new ExampleFoldLocatorMissingCoordinates
    ExampleFoldLocatorMissingCoordinates *Clone() const
    {
      return new ExampleFoldLocatorMissingCoordinates( *this);
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
    //! this is performing the execution of the example
    int Run() const
    {
      // read in protein model without loop coordinates for testing
      const std::string pdb_file_name( AddExampleInputPathToFilename( e_Fold, "1x91_no_loops.pdb"));
      io::IFStream pdb_file;
      BCL_ExampleMustOpenInputFile( pdb_file, pdb_file_name);
      const pdb::Handler pdb_handler( pdb_file, true);
      const pdb::Factory factory( biol::GetAAClasses().e_AABackBone);
      const assemble::ProteinModel model_no_loops( factory.ProteinModelFromPDB( pdb_handler));
      io::File::CloseClearFStream( pdb_file);

      // read in protein model with some loop coordinates for testing
      const std::string pdb_file_name_2( AddExampleInputPathToFilename( e_Fold, "1x91_some_loop_coords.pdb"));
      io::IFStream pdb_file_2;
      BCL_ExampleMustOpenInputFile( pdb_file_2, pdb_file_name_2);
      const pdb::Handler pdb_handler_2( pdb_file_2, true);
      const assemble::ProteinModel model_some_loops( factory.ProteinModelFromPDB( pdb_handler_2));
      io::File::CloseClearFStream( pdb_file_2);

      // read in protein model with all loop coordinates for testing
      const std::string pdb_file_name_3( AddExampleInputPathToFilename( e_Fold, "1x91.pdb"));
      io::IFStream pdb_file_3;
      BCL_ExampleMustOpenInputFile( pdb_file_3, pdb_file_name_3);
      const pdb::Handler pdb_handler_3( pdb_file_3, true);
      const assemble::ProteinModel model_all_loops( factory.ProteinModelFromPDB( pdb_handler_3));
      io::File::CloseClearFStream( pdb_file_3);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // ignore terminal loop regions and find undefined loops
      fold::LocatorMissingCoordinates locator_default;

      // don't ignore terminal loop regions and find undefined loops
      fold::LocatorMissingCoordinates locator_termini( false);

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( locator_default.GetClassIdentifier(), ( GetStaticClassName< fold::LocatorMissingCoordinates>()));

    ////////////////
    // operations //
    ////////////////

      // test locator for a model with no undefined loops besides missing terminal loops
      typedef fold::LocatorMissingCoordinates::Span Span;
      const storage::Vector< Span> spans_1( locator_default.Locate( model_all_loops));
      BCL_ExampleCheck( spans_1.IsEmpty(), true);

      // test locator_termini for a model with no undefined loops besides missing terminal loops
      const storage::Vector< Span> spans_2( locator_termini.Locate( model_all_loops));
      BCL_ExampleCheck( spans_2.GetSize(), 1);
      const Span span_2_1( spans_2( 0));
      BCL_ExampleCheck( span_2_1.First() == 1 && span_2_1.Second() == 4 && span_2_1.Third() == 'A', true);

      // test locator for a model with some undefined loops besides missing terminal loops
      const storage::Vector< Span> spans_3( locator_default.Locate( model_some_loops));
      BCL_ExampleCheck( spans_3.GetSize(), 3);

      // test locator for a model with all loops undefined
      const storage::Vector< Span> spans_4( locator_default.Locate( model_no_loops));
      BCL_ExampleCheck( spans_4.GetSize(), 5);

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

  }; // class ExampleFoldLocatorMissingCoordinates

  //! single instance of this class
  const ExampleClass::EnumType ExampleFoldLocatorMissingCoordinates::s_Instance
  (
    GetExamples().AddEnum( ExampleFoldLocatorMissingCoordinates())
  );

} // namespace bcl
