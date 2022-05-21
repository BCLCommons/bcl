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
#include "score/bcl_score_data_set_pairwise_sse_center.h"

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
  //! @example example_score_data_set_pairwise_sse_center.cpp
  //!
  //! @author alexanns
  //! @date Jun 29, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreDataSetPairwiseSSECenter :
    public ExampleInterface
  {
  public:

    ExampleScoreDataSetPairwiseSSECenter *Clone() const
    {
      return new ExampleScoreDataSetPairwiseSSECenter( *this);
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
        const double weight( score::DataSetPairwiseSSECenter::CalculateWeight( sse, resi));

        BCL_MessageDbg
        (
          "weight of residue " + resi.GetIdentification() + " in sse " + sse.GetIdentification()
          + " is " + util::Format()( weight)
        );
        const double expected_weight( 4.0 / 9.0);
        BCL_ExampleCheckWithinTolerance( weight, expected_weight, 0.001);
      }
      {
        const assemble::SSE &sse( **sses.Begin());
        const biol::AABase &resi( *sse.GetLastAA());
        const double weight( score::DataSetPairwiseSSECenter::CalculateWeight( sse, resi));

        BCL_MessageDbg
        (
          "weight of residue " + resi.GetIdentification() + " in sse " + sse.GetIdentification()
          + " is " + util::Format()( weight)
        );
        const double expected_weight( 4.0 / 9.0);
        BCL_ExampleCheckWithinTolerance( weight, expected_weight, 0.001);
      }
      // test calculate weight function with second to last residue of first sse
      {
        const assemble::SSE &sse( **sses.Begin());
        const biol::AABase &resi( **----sse.End()); //< seconcd to last residue of first sse
        const double weight( score::DataSetPairwiseSSECenter::CalculateWeight( sse, resi));

        BCL_MessageDbg
        (
          "weight of residue " + resi.GetIdentification() + " in sse " + sse.GetIdentification()
          + " is " + util::Format()( weight)
        );
        const double expected_weight( 3.0 / 9.0);
        BCL_ExampleCheckWithinTolerance( weight, expected_weight, 0.001);
      }

      // test calculate weight function with middle residue of second sse
      {
        const assemble::SSE &sse( **++++++++sses.Begin());
        const biol::AABase &resi( **++++++++++++sse.Begin()); //< residue 45 in the second sse
        const double weight( score::DataSetPairwiseSSECenter::CalculateWeight( sse, resi));

        BCL_MessageDbg
        (
          "weight of residue " + resi.GetIdentification() + " in sse " + sse.GetIdentification()
          + " is " + util::Format()( weight)
        );
        const double expected_weight( 0.0416667);
        BCL_ExampleCheckWithinTolerance( weight, expected_weight, 0.001);
      }

      // test CalculateResiduePositionWeight
      BCL_MessageDbg( "test CalculateResiduePositionWeight");
      storage::Map
      <
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, double,
        assemble::LocatorAtomCoordinatesInterface::PtrResidueLessThan
      > resis_weights( score::DataSetPairwiseSSECenter::CalculateResiduePositionWeight( *sse_pool));

      BCL_MessageDbg( "test resis_weights size");
      BCL_ExampleCheck( resis_weights.GetSize(), 101);

      // check weight of first resi in first sse
      {
        BCL_MessageDbg( "check weight of first resi in first sse");
        BCL_ExampleCheck( resis_weights.Begin()->second, 4.0 / 9.0);
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
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, double,
          assemble::LocatorAtomCoordinatesInterface::PtrLessThan
        >::const_iterator resi_itr( resis_weights.Find( resi));
        BCL_ExampleAssert( resi_itr != resis_weights.End(), true);
        BCL_ExampleCheck( resi_itr->second, 4.0 / 9.0);
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
          util::ShPtr< assemble::LocatorAtomCoordinatesInterface>, double,
          assemble::LocatorAtomCoordinatesInterface::PtrLessThan
        >::const_iterator resi_itr( resis_weights.Find( resi));
        BCL_ExampleAssert( resi_itr != resis_weights.End(), true);
        BCL_ExampleCheckWithinTolerance( resi_itr->second, 0.2222, 0.001);
      }

      score::DataSetPairwiseSSECenter score_function( *sse_pool);

      // test score operator
      {
        restraint::DataSetPairwise data_set;
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> a( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 3));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> b( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 39));
        data_set.Insert( a, b);

        const double score( score_function( data_set));
        const double expected_score( -0.0972222);
        BCL_MessageDbg
        (
          "score is " + util::Format()( score) + " expected score is "
          + util::Format()( expected_score)
        );
        BCL_ExampleCheckWithinTolerance( score, expected_score, 0.001);
      }
      // test score operator
      {
        restraint::DataSetPairwise data_set;
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> a( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 3));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> b( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 39));
        data_set.Insert( a, b);

        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> c( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 44));
        util::ShPtr< assemble::LocatorAtomCoordinatesInterface> d( new restraint::LocatorCoordinatesFirstSideChainAtom( 'A', 70));
        data_set.Insert( c, d);

        const double score( score_function( data_set));
        const double expected_score( -0.527778);
        BCL_MessageDbg
        (
          "score is " + util::Format()( score) + " expected score is "
          + util::Format()( expected_score)
        );
        BCL_ExampleCheckWithinTolerance( score, expected_score, 0.001);
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

  }; //end ExampleScoreDataSetPairwiseSSECenter

  const ExampleClass::EnumType ExampleScoreDataSetPairwiseSSECenter::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreDataSetPairwiseSSECenter())
  );

} // namespace bcl
