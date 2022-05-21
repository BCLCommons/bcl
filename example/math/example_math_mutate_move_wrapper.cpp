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
#include "math/bcl_math_mutate_move_wrapper.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_atom.h"
#include "coord/bcl_coord_move_translate_defined.h"
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_mutate_move_wrapper.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathMutateMoveWrapper :
    public ExampleInterface
  {
  public:

    //! @brief clone function
    //! @return pointer to a new ExampleMathMutateMoveWrapper
    ExampleMathMutateMoveWrapper *Clone() const
    {
      return new ExampleMathMutateMoveWrapper( *this);
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
      // initialize pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_helix.pdb"));

      // construct the model
      assemble::ProteinModel model( Proteins::GetModel( pdb_filename));

      // construct translation
      linal::Vector3D translation( 1.0, 2.0, 3.0);

      // construct move interface
      util::ShPtr< coord::MoveInterface> move( new coord::MoveTranslateDefined( translation, false));

      // expect static class name for this move
      const std::string move_class_name( GetStaticClassName( *move));

      // initialize locator
      assemble::LocatorAtom atom_locator( 'A', 23, biol::GetAtomTypes().N);

      // store the expected location
      linal::Vector3D expected_coords( 32.230, 22.705, 18.738);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      BCL_MessageStd( "test default constructor");
      math::MutateMoveWrapper< assemble::ProteinModel> mutate_def;

      // constructor from a MoveInterface
      BCL_MessageStd( "test constructor from a MoveInterface");
      math::MutateMoveWrapper< assemble::ProteinModel> mutate_a( *move);

      // default constructor
      BCL_MessageStd( "test constructor from a ShPtr to a MoveInterface");
      math::MutateMoveWrapper< assemble::ProteinModel> mutate_b( move);

      // default constructor
      BCL_MessageStd( "test clone constructor");
      util::ShPtr< math::MutateMoveWrapper< assemble::ProteinModel> > sp_mutate( mutate_a.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck
      (
        GetStaticClassName< math::MutateMoveWrapper< assemble::ProteinModel> >(),
        math::MutateMoveWrapper< assemble::ProteinModel>().GetClassIdentifier()
      );

    ///////////////
    // operators //
    ///////////////

      BCL_MessageStd( "test operator function");

      // call the mutate and store the result
      math::MutateResult< assemble::ProteinModel> mutate_result( mutate_b( model));

      // make sure mutate returned a success
      BCL_Example_Check
      (
        mutate_result.GetArgument().IsDefined(),
        "The mutate failed!"
      );

      // get the coordinates
      linal::Vector3D coords( atom_locator.Locate( *mutate_result.GetArgument()));

      BCL_Example_Check
      (
        math::EqualWithinTolerance( coords, expected_coords),
        "The coordinates of the N of residue 23 should be\n" + util::Format()( expected_coords) +
        "\n but instead it is\n" + util::Format()( coords)
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write object
      WriteBCLObject( mutate_a);

      // read object
      math::MutateMoveWrapper< assemble::ProteinModel> mutate_read;
      ReadBCLObject( mutate_read);
      math::MutateResult< assemble::ProteinModel> mutate_read_result( mutate_read( model));
      // make sure mutate returned a success
      BCL_Example_Check
      (
        mutate_read_result.GetArgument().IsDefined(),
        "The mutate failed!"
      );

      // get the coordinates
      linal::Vector3D coords_read( atom_locator.Locate( *mutate_read_result.GetArgument()));

      BCL_Example_Check
      (
        math::EqualWithinTolerance( coords_read, expected_coords),
        "The atom coordinates from mutate_read should be\n" + util::Format()( expected_coords) +
        "\n but instead it is\n" + util::Format()( coords_read)
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathMutateMoveWrapper

  const ExampleClass::EnumType ExampleMathMutateMoveWrapper::s_Instance
  (
    GetExamples().AddEnum( ExampleMathMutateMoveWrapper())
  );

} // namespace bcl
