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
#include "score/bcl_score_data_set_pairwise_sequence_separation.h"

// includes from bcl - sorted alphabetically
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "restraint/bcl_restraint_locator_coordinates_first_side_chain_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_data_set_pairwise_sequence_separation.cpp
  //!
  //! @author alexanns
  //! @date May 26, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreDataSetPairwiseSequenceSeparation :
    public ExampleInterface
  {
  public:

    ExampleScoreDataSetPairwiseSequenceSeparation *Clone() const
    {
      return new ExampleScoreDataSetPairwiseSequenceSeparation( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      score::DataSetPairwiseSequenceSeparation score( 100);

      // test data set with some data
      {
        restraint::DataSetPairwise data_set;

        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> a( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 32));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> b( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 67));
        data_set.Insert( a, b);

        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> c( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 10));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> d( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 15));
        data_set.Insert( c, d);

        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> e( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 1));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> f( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 99));
        data_set.Insert( e, f);

        const double score_of_dataset( score( data_set));
        const double expected_score( -0.705711);

        BCL_MessageDbg( "the score of the dataset is " + util::Format()( score_of_dataset));
        BCL_MessageDbg( "the expected score of the dataset is " + util::Format()( expected_score));
        BCL_ExampleCheckWithinTolerance( score_of_dataset, expected_score, 0.001);
      }

      // test data set with some data
      {
        restraint::DataSetPairwise data_set;

        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> a( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 30));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> b( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 60));
        data_set.Insert( a, b);

        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> c( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 10));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> d( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 40));
        data_set.Insert( c, d);

        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> e( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 20));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> f( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 80));
        data_set.Insert( e, f);

        const double score_of_dataset( score( data_set));
        const double expected_score( -0.788732);

        BCL_MessageDbg( "the score of the dataset is " + util::Format()( score_of_dataset));
        BCL_MessageDbg( "the expected score of the dataset is " + util::Format()( expected_score));
        BCL_ExampleCheckWithinTolerance( score_of_dataset, expected_score, 0.001);
      }

      // test empty data set
      {
        BCL_MessageDbg( "test empty dataset");
        restraint::DataSetPairwise empty_data_set;
        const double score_of_dataset( score( empty_data_set));
        const double expected_score( 0);

        BCL_MessageDbg( "the score of the empty dataset is " + util::Format()( score_of_dataset));

        BCL_ExampleCheck( score_of_dataset, expected_score);
      }

      // test data set with data point on same position
      {
        BCL_MessageDbg( "test dataset using the same residue at least twice");
        restraint::DataSetPairwise data_set_one;

        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> a( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 10));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> b( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 30));
        data_set_one.Insert( a, b);

        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> c( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 10));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> d( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 50));
        data_set_one.Insert( c, d);

        const double score_of_dataset( score( data_set_one));
        const double expected_score( -0.725772);

        BCL_MessageDbg( "the score of the dataset is with one" + util::Format()( score_of_dataset));

        BCL_ExampleCheckWithinTolerance( score_of_dataset, expected_score, 0.001);
      }

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreDataSetPairwiseSequenceSeparation

  const ExampleClass::EnumType ExampleScoreDataSetPairwiseSequenceSeparation::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreDataSetPairwiseSequenceSeparation())
  );

} // namespace bcl
