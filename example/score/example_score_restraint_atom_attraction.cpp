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
#include "score/bcl_score_restraint_atom_attraction.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_const_function.h"
#include "restraint/bcl_restraint_atom_distance.h"
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_restraint_atom_attraction.cpp
  //!
  //! @author akinlr, weinerbe
  //! @date Jul 06, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreRestraintAtomAttraction :
    public ExampleInterface
  {
  public:

    ExampleScoreRestraintAtomAttraction *Clone() const
    {
      return new ExampleScoreRestraintAtomAttraction( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // ensure that s_Instance is defined for math::ConstFunction
      math::ConstFunction< double, double>::s_Instance.IsDefined();

      // get a protein model
      assemble::ProteinModel protein_model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb")));

      // create locators to find atoms to be used in restraints
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_43_hd
      (
        new restraint::LocatorCoordinatesHydrogen( 'A', 43, biol::GetAtomTypes().HD21)
      );
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_30_hd
      (
        new restraint::LocatorCoordinatesHydrogen( 'A', 30, biol::GetAtomTypes().HD13)
      );
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_17_hg
      (
        new restraint::LocatorCoordinatesHydrogen( 'A', 17, biol::GetAtomTypes().HG11)
      );
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_13_hg
      (
        new restraint::LocatorCoordinatesHydrogen( 'A', 13, biol::GetAtomTypes().HG13)
      );
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_1_hb
      (
        new restraint::LocatorCoordinatesHydrogen( 'A', 1, biol::GetAtomTypes().HB3)
      );
      util::ShPtr< assemble::LocatorAtomCoordinatesInterface> locator_71_hb
      (
        new restraint::LocatorCoordinatesHydrogen( 'A', 71, biol::GetAtomTypes().HB3)
      );

      // create restraints
      restraint::AtomDistance restraint_a
      (
        locator_43_hd,
        locator_30_hd,
        util::ShPtr< restraint::Distance>( new restraint::Distance( 3.0, 3.75, 1.8))
      );
      restraint::AtomDistance restraint_b
      (
        locator_17_hg,
        locator_13_hg,
        util::ShPtr< restraint::Distance>( new restraint::Distance( 6.0, 7.0, 1.8))
      );
      restraint::AtomDistance restraint_c
      (
        locator_1_hb,
        locator_71_hb,
        util::ShPtr< restraint::Distance>( new restraint::Distance( 5.0, 6.5, 1.8))
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test Default constructor
      BCL_MessageStd( "test Default def_constructor");
      score::RestraintAtomAttraction def_constr;
      BCL_Example_Check
      (
        score::RestraintAtomAttraction::GetDefaultScheme() == def_constr.GetScheme(),
        "Default def_constructor should return " + util::Format()( def_constr.GetScheme()) +
        " but instead returns " + util::Format()( score::RestraintAtomAttraction::GetDefaultScheme())
      );

      // test Clone constructor
      BCL_MessageStd( "test Clone def_constructor");
      const util::ShPtr< score::RestraintAtomAttraction> clone_constr( def_constr.Clone());
      BCL_Example_Check
      (
        def_constr.GetScheme() == clone_constr->GetScheme(),
        "Clone def_constructor should return " + util::Format()( def_constr.GetScheme()) +
        " but instead returned " + util::Format()( clone_constr->GetScheme())
      );

      // create a normal constructor
      score::RestraintAtomAttraction orig_constr( math::Range< double>( -1.0, 0.0), 2.75, 15.0, true);

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_MessageStd( "test GetClassIdentifier function");
      const std::string correct_static_class_name( "bcl::score::RestraintAtomAttraction");
      BCL_Example_Check
      (
        GetStaticClassName< score::RestraintAtomAttraction>() == clone_constr->GetClassIdentifier() &&
        GetStaticClassName< score::RestraintAtomAttraction>() == correct_static_class_name,
        "GetClassIdentifier gives " + clone_constr->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

      // test GetScheme
      const std::string correct_scheme( "atom_attraction");
      BCL_ExampleCheck( def_constr.GetScheme(), correct_scheme);

    ///////////////
    // operators //
    ///////////////

      // test () operator for good restraint
      const double correct_score_a( -0.419517);
      const double calc_score_a( orig_constr( restraint_a.GenerateAssignment( protein_model)));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( calc_score_a, correct_score_a),
        true,
        "operator () with good restraint should be " + util::Format()( correct_score_a) +
        " but is " + util::Format()( calc_score_a)
      );

      // test () operator for ok restraint
      const double correct_score_b( -0.381599);
      const double calc_score_b( orig_constr( restraint_b.GenerateAssignment( protein_model)));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( calc_score_b, correct_score_b),
        true,
        "operator () with ok restraint should be " + util::Format()( correct_score_b) +
        " but is " + util::Format()( calc_score_b)
      );

      // test () operator for bad restraint
      const double correct_score_c( 0.0);
      const double calc_score_c( orig_constr( restraint_c.GenerateAssignment( protein_model)));
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( calc_score_c, correct_score_c),
        true,
        "operator () with bad restraint should be " + util::Format()( correct_score_c) +
        " but is " + util::Format()( calc_score_c)
      );
      // if message level is debug
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        // iterate over distances
        for( double i( 45); i >= 0; i -= 0.1)
        {
          // create an assignment
          restraint::AtomDistanceAssignment test_assignment
          (
            biol::Atom( linal::Vector3D( 0, 0, 0), biol::GetAtomTypes().H),
            biol::Atom( linal::Vector3D( 0, 0, i), biol::GetAtomTypes().HD1),
            util::ShPtr< restraint::Distance>( new restraint::Distance( 15.0, 20.0, 1.8))
          );

          // write out sc-cb and score
          BCL_MessageDbg
          (
            util::Format()( 15 - i) + "\t" + util::Format()( orig_constr( test_assignment))
          );
        }
      }
    //////////////////////
    // input and output //
    //////////////////////

      // test read and write
      WriteBCLObject( orig_constr);
      score::RestraintAtomAttraction read_constr;
      ReadBCLObject( read_constr);
      BCL_ExampleIndirectCheck
      (
        orig_constr.GetScheme() == read_constr.GetScheme(),
        true,
        "read and write"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreRestraintAtomAttraction

  const ExampleClass::EnumType ExampleScoreRestraintAtomAttraction::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreRestraintAtomAttraction())
  );

} // namespace bcl
