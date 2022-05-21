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
#include "restraint/bcl_restraint_atom_distance_assignment.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_atom.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_atom_distance_assignment.cpp
  //!
  //! @author weinerbe
  //! @date Mar 16, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintAtomDistanceAssignment :
    public ExampleInterface
  {
  public:

    ExampleRestraintAtomDistanceAssignment *Clone() const
    {
      return new ExampleRestraintAtomDistanceAssignment( *this);
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
      // create locators
      const assemble::LocatorAtom locator_43_cb( 'A', 43, biol::GetAtomTypes().CB);
      const assemble::LocatorAtom locator_30_cb( 'A', 30, biol::GetAtomTypes().CB);

      // create distance
      const util::ShPtr< restraint::Distance> sp_distance( new restraint::Distance( 3.0, 3.75, 1.8));

      // get a protein model
      const assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"))
      );

      // locate the atoms
      const biol::Atom atom_43_cb( *( locator_43_cb.LocateAtom( protein_model)));
      const biol::Atom atom_30_cb( *( locator_30_cb.LocateAtom( protein_model)));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      restraint::AtomDistanceAssignment def_construct;
      BCL_ExampleIndirectCheck( util::IsDefined( def_construct.CalculateAtomDistance()), false, "default constructor");

      // test constructor from atoms and distance
      const restraint::AtomDistanceAssignment atom_construct( atom_43_cb, atom_30_cb, sp_distance);

    /////////////////
    // data access //
    /////////////////

      // test GetAtomA
      BCL_ExampleCheck( atom_construct.GetAtomA().GetType(), atom_43_cb.GetType());

      // test GetAtomB
      BCL_ExampleCheck( atom_construct.GetAtomB().GetType(), atom_30_cb.GetType());

      // test GetDistance
      BCL_ExampleCheck( atom_construct.GetDistance(), 3.0);

      // test GetUpperBound
      BCL_ExampleCheck( atom_construct.GetUpperBound(), 3.75);

      // test GetLowerBound
      BCL_ExampleCheck( atom_construct.GetLowerBound(), 1.8);

      // test CalculateAtomDistance
      const double correct_distance( 8.52192);
      BCL_ExampleCheckWithinTolerance( atom_construct.CalculateAtomDistance(), correct_distance, 0.001);

    //////////////////////
    // input and output //
    //////////////////////

      // test read and write
      WriteBCLObject( atom_construct);
      ReadBCLObject( def_construct);
      BCL_ExampleIndirectCheck( atom_construct.GetDistance(), def_construct.GetDistance(), "read and write");

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintAtomDistanceAssignment

  const ExampleClass::EnumType ExampleRestraintAtomDistanceAssignment::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintAtomDistanceAssignment())
  );
  
} // namespace bcl
