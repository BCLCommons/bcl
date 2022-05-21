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
#include "score/bcl_score_data_set_pairwise_bipolar.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "restraint/bcl_restraint_data_set_pairwise.h"
#include "restraint/bcl_restraint_locator_coordinates_first_side_chain_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_data_set_pairwise_bipolar.cpp
  //!
  //! @author alexanns
  //! @date Jun 29, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreDataSetPairwiseBipolar :
    public ExampleInterface
  {
  public:

    ExampleScoreDataSetPairwiseBipolar *Clone() const
    {
      return new ExampleScoreDataSetPairwiseBipolar( *this);
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

      const assemble::ProteinModel model( Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2LZM.pdb")));

      util::SiPtrVector< const assemble::SSE> model_sses( model.GetSSEs());

      util::ShPtr< assemble::SSEPool> sse_pool( new assemble::SSEPool( model_sses));

      sse_pool->WriteSSEPool( std::cout);

      // the sses from the pool
      storage::Set< util::SiPtr< const assemble::SSE>, assemble::SSELessThanNoOverlap> sses
      (
        sse_pool->GetRandomNonOverlappingSet()
      );

      // test calculate weight function
      {
        const assemble::SSE &sse( **sses.Begin());
        const biol::AABase &resi( *sse.GetFirstAA());
        const bool nterminus( false);
        const double weight( score::DataSetPairwiseBipolar::CalculateWeight( sse, resi, nterminus));

        BCL_MessageDbg
        (
          "weight of residue " + resi.GetIdentification() + " in sse " + sse.GetIdentification()
          + " with nterminus " + util::Format()( nterminus) + " is " + util::Format()( weight)
        );
        const double expected_weight( 0);
        BCL_ExampleCheck( weight, expected_weight);
      }
      {
        const assemble::SSE &sse( **sses.Begin());
        const biol::AABase &resi( *sse.GetLastAA());
        const bool nterminus( false);
        const double weight( score::DataSetPairwiseBipolar::CalculateWeight( sse, resi, nterminus));

        BCL_MessageDbg
        (
          "weight of residue " + resi.GetIdentification() + " in sse " + sse.GetIdentification()
          + " with nterminus " + util::Format()( nterminus) + " is " + util::Format()( weight)
        );
        const double expected_weight( 1);
        BCL_ExampleCheck( weight, expected_weight);
      }
      // test calculate weight function
      {
        const assemble::SSE &sse( **sses.Begin());
        const biol::AABase &resi( *sse.GetFirstAA());
        const bool nterminus( true);
        const double weight( score::DataSetPairwiseBipolar::CalculateWeight( sse, resi, nterminus));

        BCL_MessageDbg
        (
          "weight of residue " + resi.GetIdentification() + " in sse " + sse.GetIdentification()
          + " with nterminus " + util::Format()( nterminus) + " is " + util::Format()( weight)
        );
        const double expected_weight( 1);
        BCL_ExampleCheck( weight, expected_weight);
      }
      {
        const assemble::SSE &sse( **sses.Begin());
        const biol::AABase &resi( *sse.GetLastAA());
        const bool nterminus( true);
        const double weight( score::DataSetPairwiseBipolar::CalculateWeight( sse, resi, nterminus));

        BCL_MessageDbg
        (
          "weight of residue " + resi.GetIdentification() + " in sse " + sse.GetIdentification()
          + " with nterminus " + util::Format()( nterminus) + " is " + util::Format()( weight)
        );
        const double expected_weight( 0);
        BCL_ExampleCheck( weight, expected_weight);
      }

      // test CalculateResiduePositionWeight
      BCL_MessageDbg( "test CalculateResiduePositionWeight");
      storage::Map
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, storage::Pair< double, util::ShPtr< assemble::SSE> >,
        assemble::LocatorAtomCoordinatesInterface::PtrResidueLessThan
      > resis_weights( score::DataSetPairwiseBipolar::CalculateResiduePositionWeight( *sse_pool));

      BCL_MessageDbg( "test resis_weights size");
      BCL_ExampleCheck( resis_weights.GetSize(), 101);

      // check weight of first resi in first sse
      {
        BCL_MessageDbg( "check weight of first resi in first sse");
        BCL_ExampleCheck( resis_weights.Begin()->second.First(), 0);
      }
      // check weight of last resi in first sse
      {
        BCL_MessageDbg( "check weight of last resi in first sse");
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> resi
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 11)
        );
        storage::Map
        <
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, storage::Pair< double, util::ShPtr< assemble::SSE> >,
          assemble::LocatorAtomCoordinatesInterface::PtrLessThan
        >::const_iterator resi_itr( resis_weights.Find( resi));
        BCL_ExampleAssert( resi_itr != resis_weights.End(), true);
        BCL_ExampleCheck( resi_itr->second.First(), 1);
      }

      // check weight of middle resi in first sse
      {
        BCL_MessageDbg( "check weight of middle resi in first sse");
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> resi
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 5)
        );
        storage::Map
        <
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, storage::Pair< double, util::ShPtr< assemble::SSE> >,
          assemble::LocatorAtomCoordinatesInterface::PtrLessThan
        >::const_iterator resi_itr( resis_weights.Find( resi));
        BCL_ExampleAssert( resi_itr != resis_weights.End(), true);
        BCL_ExampleCheck( resi_itr->second.First(), 0.25);
      }

      // check weight of first resi in second sse
      {
        BCL_MessageDbg( "check weight of first resi in second sse");
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> resi
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 39)
        );
        storage::Map
        <
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, storage::Pair< double, util::ShPtr< assemble::SSE> >,
          assemble::LocatorAtomCoordinatesInterface::PtrLessThan
        >::const_iterator resi_itr( resis_weights.Find( resi));
        BCL_ExampleAssert( resi_itr != resis_weights.End(), true);
        BCL_ExampleCheck( resi_itr->second.First(), 1);
      }

      // check weight of last resi in second sse
      {
        BCL_MessageDbg( "check weight of last resi in second sse");
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> resi
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 50)
        );
        storage::Map
        <
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, storage::Pair< double, util::ShPtr< assemble::SSE> >,
          assemble::LocatorAtomCoordinatesInterface::PtrLessThan
        >::const_iterator resi_itr( resis_weights.Find( resi));
        BCL_ExampleAssert( resi_itr != resis_weights.End(), true);
        BCL_ExampleCheck( resi_itr->second.First(), 0);
      }

      // check weight of middle resi in second sse
      {
        BCL_MessageDbg( "check weight of middle resi in second sse");
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> resi
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 45)
        );
        storage::Map
        <
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, storage::Pair< double, util::ShPtr< assemble::SSE> >,
          assemble::LocatorAtomCoordinatesInterface::PtrLessThan
        >::const_iterator resi_itr( resis_weights.Find( resi));
        BCL_ExampleAssert( resi_itr != resis_weights.End(), true);
        BCL_ExampleCheckWithinTolerance( resi_itr->second.First(), double( 5.0) / double( 11.0), 0.001);
      }

      // check weight of first resi in third sse
      {
        BCL_MessageDbg( "check weight of first resi in third sse");
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> resi
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 60)
        );
        storage::Map
        <
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, storage::Pair< double, util::ShPtr< assemble::SSE> >,
          assemble::LocatorAtomCoordinatesInterface::PtrLessThan
        >::const_iterator resi_itr( resis_weights.Find( resi));
        BCL_ExampleAssert( resi_itr != resis_weights.End(), true);
        BCL_ExampleCheck( resi_itr->second.First(), 0);
      }

      // check weight of last resi in third sse
      {
        BCL_MessageDbg( "check weight of last resi in third sse");
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> resi
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 80)
        );
        storage::Map
        <
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, storage::Pair< double, util::ShPtr< assemble::SSE> >,
          assemble::LocatorAtomCoordinatesInterface::PtrLessThan
        >::const_iterator resi_itr( resis_weights.Find( resi));
        BCL_ExampleAssert( resi_itr != resis_weights.End(), true);
        BCL_ExampleCheck( resi_itr->second.First(), 1);
      }

      // check weight of middle resi in third sse
      {
        BCL_MessageDbg( "check weight of middle resi in third sse");
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> resi
        (
          new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 70)
        );
        storage::Map
        <
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, storage::Pair< double, util::ShPtr< assemble::SSE> >,
          assemble::LocatorAtomCoordinatesInterface::PtrLessThan
        >::const_iterator resi_itr( resis_weights.Find( resi));
        BCL_ExampleAssert( resi_itr != resis_weights.End(), true);
        BCL_ExampleCheck( resi_itr->second.First(), 0.5);
      }

      score::DataSetPairwiseBipolar score_function( *sse_pool);

      // test score operator
      {
        restraint::DataSetPairwise data_set;
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> a( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 3));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> b( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 39));
        data_set.Insert( a, b);

        const double score( score_function( data_set));
        const double expected_score( -0.5);
        BCL_ExampleCheck( score, expected_score);
      }
      // test score operator
      {
        restraint::DataSetPairwise data_set;
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> a( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 3));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> b( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 39));
        data_set.Insert( a, b);

        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> c( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 3));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> d( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 11));
        data_set.Insert( c, d);

        const double score( score_function( data_set));
        const double expected_score( -0.55555);
        BCL_ExampleCheckWithinTolerance( score, expected_score, 0.00001);
      }
      // test score operator
      {
        restraint::DataSetPairwise data_set;
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> a( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 3));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> b( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 39));
        data_set.Insert( a, b);

        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> c( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 3));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> d( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 11));
        data_set.Insert( c, d);

        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> e( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 3));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> f( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 50));
        data_set.Insert( e, f);

        const double score( score_function( data_set));
        const double expected_score( -0.680556);
        BCL_MessageDbg
        (
          "score is " + util::Format()( score) + " but expected score is " + util::Format()( expected_score)
        );
        BCL_ExampleCheckWithinTolerance( score, expected_score, 0.000001);
      }
      // test score operator
      {
        restraint::DataSetPairwise data_set;
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> a( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 3));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> b( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 39));
        data_set.Insert( a, b);

        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> c( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 3));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> d( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 11));
        data_set.Insert( c, d);

        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> e( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 3));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> f( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 50));
        data_set.Insert( e, f);

        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> g( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 37));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> h( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 50));
        data_set.Insert( g, h);

        const double score( score_function( data_set));
        const double expected_score( -0.569444);
        BCL_ExampleCheckWithinTolerance( score, expected_score, 0.000001);
      }

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

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

  }; //end ExampleScoreDataSetPairwiseBipolar

  const ExampleClass::EnumType ExampleScoreDataSetPairwiseBipolar::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreDataSetPairwiseBipolar())
  );

} // namespace bcl
