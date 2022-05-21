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
#include "score/bcl_score_aa_pair_clash.h"

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
  //! @example example_score_aa_pair_clash.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAAPairClash :
    public ExampleInterface
  {
  public:

    ExampleScoreAAPairClash *Clone() const
    {
      return new ExampleScoreAAPairClash( *this);
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
      score::AAPairClash score;

      // make sure the histogram filename is correct
      BCL_Example_Check
      (
        score.GetHistogramFilename() == score::AAPairClash::GetDefaultHistogramFilename(),
        "The histogram file name for the default constructor should be " + score::AAPairClash::GetDefaultHistogramFilename()
        + " but is " + score.GetHistogramFilename()
      );

      // make sure that the scheme is as is expected
      BCL_Example_Check
      (
        score.GetScheme() == score::AAPairClash::GetDefaultScheme(),
        "The scheme for the default constructor should be " + score::AAPairClash::GetDefaultScheme()
        + " but is " + score.GetDefaultScheme()
      );

      // test clone function
      util::ShPtr< score::AAPairClash> clone_constr( score.Clone());

      // make sure the histogram filename is correct
      BCL_Example_Check
      (
        clone_constr->GetHistogramFilename() == score::AAPairClash::GetDefaultHistogramFilename(),
        "The histogram file name for the default constructor should be " + score::AAPairClash::GetDefaultHistogramFilename()
        + " but is " + clone_constr->GetHistogramFilename()
      );

      // make sure that the scheme is as is expected
      BCL_Example_Check
      (
        clone_constr->GetScheme() == score::AAPairClash::GetDefaultScheme(),
        "The scheme for the default constructor should be " + score::AAPairClash::GetDefaultScheme()
        + " but is " + clone_constr->GetDefaultScheme()
      );

    /////////////////
    // data access //
    /////////////////

      // check that the default histogram filename is not empty
      BCL_Example_Check
      (
        !score::AAPairClash::GetDefaultHistogramFilename().empty(),
        "Warning, the default histogram file name is empty"
      );
      // check that the default scheme is not empty
      BCL_Example_Check
      (
        !score::AAPairClash::GetDefaultScheme().empty(),
        "Warning, the default scheme is not defined"
      );

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::score::AAPairClash");
      const std::string given_static_class_name( GetStaticClassName< score::AAPairClash>());
      BCL_Example_Check
      (
        given_static_class_name == correct_static_class_name,
        "GetStaticClassName gives " + given_static_class_name + " but should give " + correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_Example_Check
      (
        correct_static_class_name == clone_constr->GetClassIdentifier(),
        "GetClassIdentifier gives " + clone_constr->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

      // test GetHistogramFilename
      BCL_Example_Check
      (
        clone_constr->GetHistogramFilename() == score::AAPairClash::GetDefaultHistogramFilename(),
        "The histogram file name for the default constructor should be " + score::AAPairClash::GetDefaultHistogramFilename()
        + " but is " + clone_constr->GetHistogramFilename()
      );

      // test GetScheme
      BCL_Example_Check
      (
        clone_constr->GetScheme() == score::AAPairClash::GetDefaultScheme(),
        "The scheme for the default constructor should be " + score::AAPairClash::GetDefaultScheme()
        + " but is " + clone_constr->GetDefaultScheme()
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
        biol::Atom( linal::Vector3D( 5, 5, -6.7), biol::GetAtomTypes().CA),
        biol::Atom( linal::Vector3D( 4, 4, -4.7), biol::GetAtomTypes().CB)
      );

      // create amino acid b and fill it with coordinates, type, etc.
      biol::AA simple_aa_b( util::ShPtr< biol::AAData>( new biol::AAData( biol::GetAATypes().ARG, 2, 9)));
      biol::AACaCb ca_cb_aa_b
      (
        simple_aa_b,
        biol::Atom( linal::Vector3D( 5, 5, -10), biol::GetAtomTypes().CA),
        biol::Atom( linal::Vector3D( 4, 4, -7), biol::GetAtomTypes().CB)
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
        ca_cb_aa_b.Translate( linal::Vector3D( 0, 0, -0.1));

        const double distance( biol::FirstSidechainAtomDistance( ca_cb_aa_a, ca_cb_aa_b));
        const double clash_score( score( ca_cb_aa_a, ca_cb_aa_b));

        // message score between "ca_cb_aa_a" and "ca_cb_aa_b"
        BCL_MessageStd
        (
          "dist:\t" + util::Format()( distance) + "\tscore:\t" + util::Format()( clash_score)
        );

        if( i == 6)
        {
          // check the operator results
          const double expected_value( 0.578217);
          BCL_Example_Check
          (
            math::EqualWithinTolerance
            (
              clash_score,
              expected_value
            ),
            "The score at this AAPairClash should be " + util::Format()( expected_value) + " but is " +
            util::Format()( clash_score)
          );
        }
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities for score::AAPairClash");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( *clone_constr);
      score::AAPairClash read_obj;
      BCL_MessageVrb( "read object");
      ReadBCLObject( read_obj);
      BCL_Example_Check
      (
        read_obj.GetHistogramFilename() == score::AAPairClash::GetDefaultHistogramFilename(),
        "The histogram file name for the ReadBCLObject should be " + score::AAPairClash::GetDefaultHistogramFilename()
        + " but is " + read_obj.GetHistogramFilename()
      );
      BCL_Example_Check
      (
        read_obj.GetScheme() == score::AAPairClash::GetDefaultScheme(),
        "The scheme for the ReadBCLObject should be " + score::AAPairClash::GetDefaultScheme()
        + " but is " + read_obj.GetDefaultScheme()
      );

      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        storage::Vector< std::string> tics;
        tics.AllocateMemory( biol::AATypes::s_NumberStandardAATypes);
        for( biol::AATypes::const_iterator itr( biol::GetAATypes().ALA.GetIterator()), itr_end( biol::GetAATypes().VAL.GetIterator()); itr <= itr_end; ++itr)
        {
          tics.PushBack( ( *itr)->GetThreeLetterCode());
        }
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, "aa_pair_clash_distances.gnuplot");
        math::GnuplotHeatmap heatmap;
        heatmap.SetFromMatrix( score.GetShortestObservedDistanceMatrix());
        heatmap.SetTitleAndLabel( "", "", "", "min distance [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]");
        heatmap.SetPixelAndRatio( 1080, 800, -1);
        heatmap.SetTicsY( tics, true, 1);
        heatmap.SetTicsX( tics, true, 1);
        heatmap.SetRotationXTics( 90);
        heatmap.SetFont( "arialbd", 20);
        heatmap.SetMinMaxZ( 2.0, 4.0);
        heatmap.SetFilename( "aa_pair_clash_distances");
        heatmap.WriteScript( write);
        io::File::CloseClearFStream( write);
      }

      // test WriteDetailedSchemeAndValues
      BCL_Example_Check
      (
        read_obj.GetHistogramFilename() == score::AAPairClash::GetDefaultHistogramFilename(),
        "The histogram file name for the ReadBCLObject should be " + score::AAPairClash::GetDefaultHistogramFilename()
        + " but is " + read_obj.GetHistogramFilename()
      );

      // create stringstream "str_stream" to hold the output from WriteDetailedSchemeAndValues
      std::stringstream str_stream;

      // get the WriteDetailedSchemeAndValues
      score.WriteDetailedSchemeAndValues( ca_cb_aa_a, ca_cb_aa_b, str_stream);

      BCL_MessageDbg( "detailed scheme and value is " + str_stream.str());

      // create string "expected_detailed_output" and initialize with what the detailed scheme and values should be
      const std::string expected_detailed_output( "1\tLYS\t2\tARG\t3.7\t0\n");

      // test WriteDetailedSchemeAndValues
      BCL_Example_Check
      (
        str_stream.str() == expected_detailed_output,
        "The detailed scheme and values should be " + expected_detailed_output
        + " but is " + str_stream.str()
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAAPairClash

  const ExampleClass::EnumType ExampleScoreAAPairClash::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAAPairClash())
  );

} // namespace bcl

