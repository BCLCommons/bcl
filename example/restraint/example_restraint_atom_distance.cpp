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
#include "restraint/bcl_restraint_atom_distance.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_atom.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "restraint/bcl_restraint_atom_distance_assignment.h"
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_restraint_atom_distance.cpp
  //!
  //! @author weinerbe
  //! @date Mar 16, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleRestraintAtomDistance :
    public ExampleInterface
  {
  public:

    ExampleRestraintAtomDistance *Clone() const
    {
      return new ExampleRestraintAtomDistance( *this);
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
      const restraint::LocatorCoordinatesHydrogen locator_43_hd( 'A', 43, biol::GetAtomTypes().HD21);
      const restraint::LocatorCoordinatesHydrogen locator_30_hd( 'A', 30, biol::GetAtomTypes().HD13);
      const assemble::LocatorAtom locator_43_cb( 'A', 43, biol::GetAtomTypes().CB);
      const assemble::LocatorAtom locator_30_cb( 'A', 30, biol::GetAtomTypes().CB);

      // create distance
      const util::ShPtr< restraint::Distance> sp_distance( new restraint::Distance( 3.0, 3.75, 1.8));

      // get a protein model
      const assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"))
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      restraint::AtomDistance def_construct;
      BCL_ExampleIndirectCheck( def_construct.IsDefined(), false, "default constructor");

      // test constructor with cb locators
      const restraint::AtomDistance cb_construct( locator_43_cb, locator_30_cb, sp_distance);

      // test constructor with shptrs
      const restraint::AtomDistance shptr_h_construct
      (
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( locator_43_hd.Clone()),
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>( locator_30_hd.Clone()),
        sp_distance
      );
      BCL_ExampleIndirectCheck( shptr_h_construct.IsDefined(), true, "constructor from shptrs");

    /////////////////
    // data access //
    /////////////////

      // write out the identification
      BCL_MessageStd( "proton restraint: " + shptr_h_construct.GetIdentification());

      // test GetData
      const restraint::DataPairwise &data_pair( shptr_h_construct.GetData());
      BCL_ExampleCheck( data_pair.IsSet(), true);

      // test GetDistance
      BCL_ExampleCheck( shptr_h_construct.GetDistance()->GetDistance(), sp_distance->GetDistance());

      // test is defined
      BCL_ExampleCheck( shptr_h_construct.IsDefined(), true);

    ////////////////
    // operations //
    ////////////////

      // test GenerateAssignment
      const restraint::AtomDistanceAssignment h_assignment( shptr_h_construct.GenerateAssignment( protein_model));
      BCL_ExampleIndirectCheck
      (
        util::IsDefined( h_assignment.CalculateAtomDistance()), true, "GenerateAssignment with protons"
      );
      const restraint::AtomDistanceAssignment cb_assignmnet( cb_construct.GenerateAssignment( protein_model));
      BCL_ExampleIndirectCheck
      (
        util::IsDefined( cb_assignmnet.CalculateAtomDistance()), true, "GenerateAssignment with CBs"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test read and write
      WriteBCLObject( shptr_h_construct);
      ReadBCLObject( def_construct);
      BCL_ExampleIndirectCheck
      (
        shptr_h_construct.GetIdentification(),
        def_construct.GetIdentification(),
        "read and write"
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleRestraintAtomDistance

  const ExampleClass::EnumType ExampleRestraintAtomDistance::s_Instance
  (
    GetExamples().AddEnum( ExampleRestraintAtomDistance())
  );

} // namespace bcl
