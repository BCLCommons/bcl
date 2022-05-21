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
#include "score/bcl_score_fuzzy_logic_filter.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "restraint/bcl_restraint_atom_distance.h"
#include "restraint/bcl_restraint_locator_coordinates_hydrogen.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_fuzzy_logic_filter.cpp
  //!
  //! @author alexanns, weinerbe
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreFuzzyLogicFilter :
    public ExampleInterface
  {
  public:

    ExampleScoreFuzzyLogicFilter *Clone() const
    {
      return new ExampleScoreFuzzyLogicFilter( *this);
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
      score::FuzzyLogicFilter def_constr;
      BCL_Example_Check
      (
        score::FuzzyLogicFilter::GetDefaultScheme() == def_constr.GetScheme(),
        "Default def_constructor should return " + util::Format()( def_constr.GetScheme()) +
        " but instead returns " + util::Format()( score::FuzzyLogicFilter::GetDefaultScheme())
      );

    /////////////////
    // data access //
    /////////////////

      // test GetScheme
      const std::string correct_scheme( "fuzzy_logic_filter");
      BCL_ExampleCheck( def_constr.GetScheme(), correct_scheme);

    ///////////////
    // operators //
    ///////////////

      // test () operator for good restraint
      const double correct_score_a( -0.816141);
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( def_constr( restraint_a.GenerateAssignment( protein_model)), correct_score_a),
        true,
        "operator () with good restraint"
      );

      // test () operator for ok restraint
      const double correct_score_b( -0.167149);
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( def_constr( restraint_b.GenerateAssignment( protein_model)), correct_score_b),
        true,
        "operator () with ok restraint"
      );

      // test () operator for bad restraint
      const double correct_score_c( 0.0);
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance( def_constr( restraint_c.GenerateAssignment( protein_model)), correct_score_c),
        true,
        "operator () with bad restraint"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test read and write
      score::FuzzyLogicFilter orig_constr( "read_write_test");
      WriteBCLObject( orig_constr);
      score::FuzzyLogicFilter read_constr;
      ReadBCLObject( read_constr);
      BCL_ExampleIndirectCheck( orig_constr.GetScheme() == read_constr.GetScheme(), true, "read and write");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

   }; //end example_score_fuzzy_logic_filter

  const ExampleClass::EnumType ExampleScoreFuzzyLogicFilter::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreFuzzyLogicFilter())
  );

} // namespace bcl
