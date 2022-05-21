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
#include "score/bcl_score_aa_pair_contact_energy.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa.h"
#include "biol/bcl_biol_aa_back_bone.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_gnuplot_heatmap.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_aa_pair_contact_energy.cpp
  //!
  //! @author mendenjl
  //! @date Dec 16, 2016
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAAPairContactEnergy :
    public ExampleInterface
  {
  public:

    ExampleScoreAAPairContactEnergy *Clone() const
    {
      return new ExampleScoreAAPairContactEnergy( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // testing default constructor
      score::AAPairContactEnergy score;

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::score::AAPairContactEnergy");
      const std::string given_static_class_name( GetStaticClassName< score::AAPairContactEnergy>());
      BCL_Example_Check
      (
        given_static_class_name == correct_static_class_name,
        "GetStaticClassName gives " + given_static_class_name + " but should give " + correct_static_class_name
      );

    ///////////////
    // operators //
    ///////////////

      // test the scoring operator

      // create amino acid a and fill it with coordinates, type, etc.
      biol::AA simple_aa_a( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().LYS, 1, 8)));
      biol::AACaCb ca_cb_aa_a
      (
        simple_aa_a,
        biol::Atom( linal::Vector3D( 0, 0, 0), biol::GetAtomTypes().CA),
        biol::Atom( linal::Vector3D( 0, 0, 1), biol::GetAtomTypes().CB)
      );

      // create amino acid b and fill it with coordinates, type, etc.
      biol::AA simple_aa_b( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ARG, 2, 9)));
      biol::AACaCb ca_cb_aa_b
      (
        simple_aa_b,
        biol::Atom( linal::Vector3D( 0, 0, 2), biol::GetAtomTypes().CA),
        biol::Atom( linal::Vector3D( 0, 0, 1), biol::GetAtomTypes().CB)
      );

      BCL_MessageStd
      (
        "moving one aminoacid in 0.9 Angstrom steps in z direction and score the distance"
      );

      BCL_MessageStd
      (
        "distance of: " + util::Format()( ca_cb_aa_a.GetType()->GetThreeLetterCode())
        + " and " + util::Format()( ca_cb_aa_b.GetType()->GetThreeLetterCode()) + " and its score:"
      );

      // do iterative translation of "ca_cb_aa_b" in the z-direction and score its interaction with "ca_cb_aa_a"
      for( size_t i( 0); i < 14; ++i)
      {
        // move "ca_cb_aa_b"
        ca_cb_aa_b.Translate( linal::Vector3D( 0, 0, 0.9));

        // calculate the score between "ca_cb_aa_a" and "ca_cb_aa_b"
        const double this_score
        (
          score
          (
            ca_cb_aa_a,
            ca_cb_aa_b,
            biol::FirstSidechainAtomDistance( ca_cb_aa_a, ca_cb_aa_b),
            biol::GetSSTypes().HELIX,
            biol::GetSSTypes().HELIX
          )
        );

        // message distance between "ca_cb_aa_a" and "ca_cb_aa_b"
        BCL_MessageStd
        (
          util::Format()
          (
            biol::FirstSidechainAtomDistance( ca_cb_aa_a, ca_cb_aa_b)
          ) + " " + util::Format()( this_score)
        );

        if( i == 3)
        {
          // check the operator results
          BCL_ExampleCheckWithinTolerance
          (
            score
            (
              ca_cb_aa_a,
              ca_cb_aa_b,
              biol::FirstSidechainAtomDistance( ca_cb_aa_a, ca_cb_aa_b),
              biol::GetSSTypes().HELIX,
              biol::GetSSTypes().HELIX
            ),
            0.217,
            0.00001
          );
        }
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities for score::AAPairContactEnergy");
      // open output file

      // test WriteDetailedSchemeAndValues

      // create stringstream "str_stream" to hold the output from WriteDetailedSchemeAndValues
      std::stringstream str_stream;

      // get the WriteDetailedSchemeAndValues
      score.WriteDetailedSchemeAndValues( ca_cb_aa_a, ca_cb_aa_b, biol::GetSSTypes().HELIX, biol::GetSSTypes().HELIX, str_stream);

      // create string "expected_detailed_output" and initialize with what the detailed scheme and values should be

      // test WriteDetailedSchemeAndValues
      BCL_ExampleIndirectCheck
      (
        str_stream.str(),
        "1\tLYS\t2\tARG\t12.6\t1\t1\t0.326\n",
        "WriteDetailedSchemeAndValues"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAAPairContactEnergy

  const ExampleClass::EnumType ExampleScoreAAPairContactEnergy::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAAPairContactEnergy())
  );

} // namespace bcl

