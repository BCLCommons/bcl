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
#include "coord/bcl_coord_move_combine.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "coord/bcl_coord_move_translate_defined.h"
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_coord_move_combine.cpp
  //!
  //! @author karakam
  //! @date Feb 22, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCoordMoveCombine :
    public ExampleInterface
  {
  public:

    ExampleCoordMoveCombine *Clone() const
    {
      return new ExampleCoordMoveCombine( *this);
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

      // construct list of translates
      const linal::Vector3D vector_a( 0.0, 5.0, 0.0);
      const linal::Vector3D vector_b( -5.0, 5.0, 2.0);
      util::ShPtr< coord::MoveTranslateDefined> sp_translate_a( new coord::MoveTranslateDefined( vector_a));
      util::ShPtr< coord::MoveTranslateDefined> sp_translate_b( new coord::MoveTranslateDefined( vector_b));
      const linal::Vector3D vector_c( -5.0, 10.0, 2.0);

      // form list
      util::ShPtrList< coord::MoveInterface> translate_list;
      translate_list.PushBack( sp_translate_a);
      translate_list.PushBack( sp_translate_b);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      const coord::MoveCombine move_default;

      // construct from a list of moves
      const coord::MoveCombine move( translate_list);

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

      // construct atom
      const linal::Vector3D coordinates_a( 0.0, 0.0, 0.0);
      biol::Atom atom_a;
      atom_a.SetCoordinates( coordinates_a);

      // move the atom and check the end coordinates
      move.Move( atom_a);
      BCL_ExampleCheck( atom_a.GetCoordinates(), vector_c);

    //////////////////////
    // input and output //
    //////////////////////

      // write object out
      WriteBCLObject( move);
      coord::MoveCombine move_read;
      ReadBCLObject( move_read);

      // create atom
      biol::Atom atom_b;
      atom_b.SetCoordinates( coordinates_a);
      // apply move and check the coordinates
      move_read.Move( atom_b);
      BCL_ExampleCheck( atom_b.GetCoordinates(), vector_c);
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleCoordMoveCombine

  const ExampleClass::EnumType ExampleCoordMoveCombine::s_Instance
  (
    GetExamples().AddEnum( ExampleCoordMoveCombine())
  );
  
} // namespace bcl
