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
#include "score/bcl_score_aa_pair_distance.h"

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
  //! @example example_score_aa_pair_distance.cpp
  //!
  //! @author rouvelgh
  //! @date Nov 6, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAAPairDistance :
    public ExampleInterface
  {
  public:

    ExampleScoreAAPairDistance *Clone() const
    {
      return new ExampleScoreAAPairDistance( *this);
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
      score::AAPairDistance score;

      // make sure the histogram filename is correct
      BCL_Example_Check
      (
        score.GetHistogramFilename() == score::AAPairDistance::GetDefaultHistogramFilename(),
        "The histogram file name for the default constructor should be " + score::AAPairDistance::GetDefaultHistogramFilename()
        + " but is " + score.GetHistogramFilename()
      );

      // make sure that the scheme is as is expected
      BCL_Example_Check
      (
        score.GetScheme() == score::AAPairDistance::GetDefaultScheme(),
        "The scheme for the default constructor should be " + score::AAPairDistance::GetDefaultScheme()
        + " but is " + score.GetDefaultScheme()
      );

      // make sure that the energy function map has the expected size
      // create size_t "expected_size" and initialize with the size the energy map should have
      const size_t expected_size( 400);
      BCL_Example_Check
      (
        score.GetEnergyFunctionMap().GetSize() == expected_size,
        "The size of the energy function map should be " + util::Format()( expected_size)
        + " but is size " + util::Format()( score.GetEnergyFunctionMap().GetSize())
      );

      // test clone function
      util::ShPtr< score::AAPairDistance> clone_constr( score.Clone());

      // make sure the histogram filename is correct
      BCL_Example_Check
      (
        clone_constr->GetHistogramFilename() == score::AAPairDistance::GetDefaultHistogramFilename(),
        "The histogram file name for the default constructor should be " + score::AAPairDistance::GetDefaultHistogramFilename()
        + " but is " + clone_constr->GetHistogramFilename()
      );

      // make sure that the scheme is as is expected
      BCL_Example_Check
      (
        clone_constr->GetScheme() == score::AAPairDistance::GetDefaultScheme(),
        "The scheme for the default constructor should be " + score::AAPairDistance::GetDefaultScheme()
        + " but is " + clone_constr->GetDefaultScheme()
      );

      // make sure that the energy map has the correct size
      BCL_Example_Check
      (
        clone_constr->GetEnergyFunctionMap().GetSize() == expected_size,
        "The size of the energy function map should be " + util::Format()( expected_size)
        + " but is size " + util::Format()( clone_constr->GetEnergyFunctionMap().GetSize())
      );

    /////////////////
    // data access //
    /////////////////

      // check that the default historam filename is not empty
      BCL_Example_Check
      (
        !score::AAPairDistance::GetDefaultHistogramFilename().empty(),
        "Warning, the default histogram file name is empty"
      );
      // check that the default scheme is not empty
      BCL_Example_Check
      (
        !score::AAPairDistance::GetDefaultScheme().empty(),
        "Warning, the default scheme is not defined"
      );

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::score::AAPairDistance");
      const std::string given_static_class_name( GetStaticClassName< score::AAPairDistance>());
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
        clone_constr->GetHistogramFilename() == score::AAPairDistance::GetDefaultHistogramFilename(),
        "The histogram file name for the default constructor should be " + score::AAPairDistance::GetDefaultHistogramFilename()
        + " but is " + clone_constr->GetHistogramFilename()
      );

      // test GetScheme
      BCL_Example_Check
      (
        clone_constr->GetScheme() == score::AAPairDistance::GetDefaultScheme(),
        "The scheme for the default constructor should be " + score::AAPairDistance::GetDefaultScheme()
        + " but is " + clone_constr->GetDefaultScheme()
      );

      // test GetEnergyFunctionMap
      BCL_Example_Check
      (
        clone_constr->GetEnergyFunctionMap().GetSize() == expected_size,
        "The size of the energy function map should be " + util::Format()( expected_size)
        + " but is size " + util::Format()( clone_constr->GetEnergyFunctionMap().GetSize())
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
        biol::Atom( linal::Vector3D( 5, 5, -9), biol::GetAtomTypes().CA),
        biol::Atom( linal::Vector3D( 4, 4, -7), biol::GetAtomTypes().CB)
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
        ca_cb_aa_b.Translate( linal::Vector3D( 0, 0, -0.9));

        // calculate the score between "ca_cb_aa_a" and "ca_cb_aa_b"
        const double this_score( score( ca_cb_aa_a, ca_cb_aa_b));

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
          BCL_ExampleCheckWithinTolerance( score( ca_cb_aa_a, ca_cb_aa_b), 0.890225, 0.00001);
        }
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      BCL_MessageStd( "testing read and write functionalities for score::AAPairDistance");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( *clone_constr);
      score::AAPairDistance read_obj;
      BCL_MessageVrb( "read object");
      ReadBCLObject( read_obj);

      // make sure read and write worked
      BCL_Example_Check
      (
        read_obj.GetEnergyFunctionMap().GetSize() == expected_size,
        "The size of the energy function map for the ReadBCLObject should be " + util::Format()( expected_size)
        + " but is size " + util::Format()( read_obj.GetEnergyFunctionMap().GetSize())
      );
      BCL_Example_Check
      (
        read_obj.GetHistogramFilename() == score::AAPairDistance::GetDefaultHistogramFilename(),
        "The histogram file name for the ReadBCLObject should be " + score::AAPairDistance::GetDefaultHistogramFilename()
        + " but is " + read_obj.GetHistogramFilename()
      );
      BCL_Example_Check
      (
        read_obj.GetScheme() == score::AAPairDistance::GetDefaultScheme(),
        "The scheme for the ReadBCLObject should be " + score::AAPairDistance::GetDefaultScheme()
        + " but is " + read_obj.GetDefaultScheme()
      );

      // print gnuplot heat map
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        util::SiPtrVector< const math::CubicSplineDamped> splines;
        storage::Vector< std::string> spline_descriptors;

        // iterate over all aa type pairs
        for
        (
          storage::Map
          <
            storage::Pair< biol::AAType, biol::AAType>,
            util::ShPtr< math::CubicSplineDamped>
          >::const_iterator itr( score.GetEnergyFunctionMap().Begin()), itr_end( score.GetEnergyFunctionMap().End());
          itr != itr_end;
          ++itr
        )
        {
          // skip redundant pairs
          if( itr->first.First() > itr->first.Second())
          {
            continue;
          }
          splines.PushBack( itr->second);
          spline_descriptors.PushBack( itr->first.First()->GetThreeLetterCode() + " " + itr->first.Second()->GetThreeLetterCode());
        }

        // write splines to gnuplot file
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, "aa_pair_distance_potentials.gnuplot");
        math::GnuplotHeatmap heatmap;
        heatmap.SetFromCubicSplines( splines, false, false);
        heatmap.SetTicsY( spline_descriptors, true, 1);
        heatmap.SetPixelAndRatio( 1080, 8000, -1);
        heatmap.SetTitleAndLabel( "amino acid pair distance potentials", "distance [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]", "", "-log(p/b)");
        heatmap.SetFont( "arialbd", 20);
        heatmap.SetMinMaxZ( -1.0, 1.0);
        heatmap.SetFilename( "aa_pair_distance_potentials");
        heatmap.WriteScript( write);
        io::File::CloseClearFStream( write);
      }

      // print gnuplot heat map of selected interactions
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        util::SiPtrVector< const math::CubicSplineDamped> splines;
        storage::Vector< std::string> spline_descriptors;

        // select a set of aa pair types - so that each aa type is contained at least once
        storage::Set< storage::Pair< biol::AAType, biol::AAType> > selected_aa_pairs;
        selected_aa_pairs.Insert( storage::Pair< biol::AAType, biol::AAType>( biol::GetAATypes().ARG, biol::GetAATypes().LYS));
        selected_aa_pairs.Insert( storage::Pair< biol::AAType, biol::AAType>( biol::GetAATypes().ARG, biol::GetAATypes().GLU));
        selected_aa_pairs.Insert( storage::Pair< biol::AAType, biol::AAType>( biol::GetAATypes().ILE, biol::GetAATypes().LEU));
        selected_aa_pairs.Insert( storage::Pair< biol::AAType, biol::AAType>( biol::GetAATypes().TRP, biol::GetAATypes().TRP));
        selected_aa_pairs.Insert( storage::Pair< biol::AAType, biol::AAType>( biol::GetAATypes().ALA, biol::GetAATypes().GLY));
        selected_aa_pairs.Insert( storage::Pair< biol::AAType, biol::AAType>( biol::GetAATypes().ASN, biol::GetAATypes().VAL));
        selected_aa_pairs.Insert( storage::Pair< biol::AAType, biol::AAType>( biol::GetAATypes().ASP, biol::GetAATypes().PHE));
        selected_aa_pairs.Insert( storage::Pair< biol::AAType, biol::AAType>( biol::GetAATypes().CYS, biol::GetAATypes().PRO));
        selected_aa_pairs.Insert( storage::Pair< biol::AAType, biol::AAType>( biol::GetAATypes().HIS, biol::GetAATypes().SER));
        selected_aa_pairs.Insert( storage::Pair< biol::AAType, biol::AAType>( biol::GetAATypes().THR, biol::GetAATypes().TYR));
        selected_aa_pairs.Insert( storage::Pair< biol::AAType, biol::AAType>( biol::GetAATypes().GLN, biol::GetAATypes().MET));

        for
        (
          storage::Set< storage::Pair< biol::AAType, biol::AAType> >::iterator
            sel_itr( selected_aa_pairs.Begin()), sel_itr_end( selected_aa_pairs.End());
          sel_itr != sel_itr_end;
          ++sel_itr
        )
        {
          const storage::Map
          <
            storage::Pair< biol::AAType, biol::AAType>,
            util::ShPtr< math::CubicSplineDamped>
          >::const_iterator
            itr( score.GetEnergyFunctionMap().Find( *sel_itr));
          splines.PushBack( itr->second);
          spline_descriptors.PushBack( itr->first.First()->GetThreeLetterCode() + " " + itr->first.Second()->GetThreeLetterCode());
        }

        // write splines to gnuplot file
        io::OFStream write;
        BCL_ExampleMustOpenOutputFile( write, "aa_pair_distance_potentials_selected.gnuplot");
        math::GnuplotHeatmap heatmap;
        heatmap.SetFromCubicSplines( splines, false, false);
        heatmap.SetTicsY( spline_descriptors, true, 1);
        heatmap.SetPixelAndRatio( 1080, 800, -1.0);
        heatmap.SetTitleAndLabel( "amino acid pair distance potentials", "distance [" + math::GnuplotHeatmap::s_AngstromSymbolGnuplot + "]", "", "-log(p/b)");
        heatmap.SetFont( "arialbd", 20);
        heatmap.SetMinMaxZ( -0.5, 0.5);
        heatmap.SetFilename( "aa_pair_distance_potentials_selected");
        heatmap.WriteScript( write);
        io::File::CloseClearFStream( write);
      }
      // test WriteDetailedSchemeAndValues
      BCL_Example_Check
      (
        read_obj.GetHistogramFilename() == score::AAPairDistance::GetDefaultHistogramFilename(),
        "The histogram file name for the ReadBCLObject should be " + score::AAPairDistance::GetDefaultHistogramFilename()
        + " but is " + read_obj.GetHistogramFilename()
      );

      // create stringstream "str_stream" to hold the output from WriteDetailedSchemeAndValues
      std::stringstream str_stream;

      // get the WriteDetailedSchemeAndValues
      score.WriteDetailedSchemeAndValues( ca_cb_aa_a, ca_cb_aa_b, str_stream);

      // create string "expected_detailed_output" and initialize with what the detailed scheme and values should be

      // test WriteDetailedSchemeAndValues
      BCL_ExampleIndirectCheck
      (
        str_stream.str(),
        "1\tLYS\t2\tARG\t12.6\t0.00157516\n",
        "WriteDetailedSchemeAndValues"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAAPairDistance

  const ExampleClass::EnumType ExampleScoreAAPairDistance::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAAPairDistance())
  );

} // namespace bcl

