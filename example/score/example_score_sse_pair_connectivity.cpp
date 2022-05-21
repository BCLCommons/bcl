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
#include "score/bcl_score_sse_pair_connectivity.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "score/bcl_score_protein_model_sse_neighbors.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_sse_pair_connectivity.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreSSEPairConnectivity :
    public ExampleInterface
  {
  public:

    ExampleScoreSSEPairConnectivity *Clone() const
    {
      return new ExampleScoreSSEPairConnectivity( *this);
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
      // create const string which will be used as a scheme for the scoring function
      const std::string scheme( "dummy_scheme");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      score::SSEPairConnectivity def_constr;

      // make sure the scheme is correct
      BCL_ExampleCheck( def_constr.GetScheme(), score::SSEPairConnectivity::GetDefaultScheme());

      // test constructor taking scheme
      score::SSEPairConnectivity param_constr( scheme);

      // make sure the scheme is correct
      BCL_ExampleCheck( param_constr.GetScheme(), scheme);

      // test clone constructor
      util::ShPtr< score::SSEPairConnectivity> clone_constr( param_constr.Clone());

      // make sure the scheme is correct
      BCL_ExampleCheck( clone_constr->GetScheme(), scheme);

      // create score for protein model neighbors
      const score::ProteinModelSSENeighbors score_neighbors( util::CloneToShPtr( def_constr), false);

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::score::ProteinModelSSEConnectivity");

      // check GetClassIdentifier
      BCL_Example_Check
      (
        GetStaticClassName< score::SSEPairConnectivity>() == clone_constr->GetClassIdentifier(),
        "GetClassIdentifier gives " + clone_constr->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

      // test get scheme function
      BCL_ExampleCheck( clone_constr->GetScheme(), scheme);

    ////////////////
    // operations //
    ////////////////

      // create protein model whose sse connectivity can be scored
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "2LZM.pdb"));
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( pdb_filename)
      );

      // test get scoring function
      // get two sses adjacent in sequence
      const util::SiPtr< const assemble::SSE> sse_1_2( assemble::LocatorSSE( 'A', 1, 2).Locate( protein_model));
      const util::SiPtr< const assemble::SSE> sse_3_11( assemble::LocatorSSE( 'A', 3, 11).Locate( protein_model));

      // check the CalculateConnectivityScore function
      const double connectivity_score( def_constr( *sse_1_2, *sse_3_11));
      BCL_ExampleCheckWithinTolerance( connectivity_score, 0.0, 0.001);

      // check the CalculateConnectivityScore function with unconnected sses
      // get protein model
      assemble::ProteinModel protein_model_b
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2LZM_rotated_35_38.pdb"))
      );
      // get two sses adjacent in sequence
      const util::SiPtr< const assemble::SSE> sse_35_38( assemble::LocatorSSE( 'A', 35, 38).Locate( protein_model_b));
      const util::SiPtr< const assemble::SSE> sse_39_50( assemble::LocatorSSE( 'A', 39, 50).Locate( protein_model_b));

      // calculate the connectivity
      const double connectivity_score_b( def_constr( *sse_35_38, *sse_39_50));
      BCL_MessageDbg
      (
        "the connectivity score for two unconnected sses is " + util::Format()( connectivity_score_b)
      );
      BCL_ExampleCheckWithinTolerance( connectivity_score_b, 193.72, 0.001);

      // check CalculateChainSSEConnectivityScore
      BCL_MessageDbg( "check CalculateChainSSEConnectivityScore");
      // calculate the connectivity
      const double connectivity_score_c( score_neighbors( *protein_model.GetChain( 'A')));
      BCL_MessageDbg
      (
        "the connectivity score for connected chain is " + util::Format()( connectivity_score_c)
      );
      // expected connectivity score
      const double expected_connectivity_score( 1.60713E-05);
      BCL_ExampleCheckWithinTolerance( connectivity_score_c, expected_connectivity_score, 0.001);

      // check CalculateChainSSEConnectivityScore with a non continuous chain
      BCL_MessageDbg( "check CalculateChainSSEConnectivityScore with non continous chain");
      const double expected_disconnected_score( 193.72);
      // calculate the connectivity
      const double connectivity_score_d
      (
        score_neighbors( *protein_model_b.GetChain( 'A'))
      );
      BCL_MessageDbg
      (
        "the connectivity score for disconnected chain is " + util::Format()( connectivity_score_d)
      );
      BCL_ExampleCheckWithinTolerance( connectivity_score_d, expected_disconnected_score, 0.001);

      // check CalculateChainSSEConnectivityScore with a chain without loops with defined coordinates
      BCL_MessageDbg
      (
        "check CalculateChainSSEConnectivityScore with a chain without loops with defined coordinates"
      );
      // get protein model
      assemble::ProteinModel protein_model_c
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"))
      );
      // calculate the connectivity
      const double connectivity_score_e
      (
        score_neighbors( *protein_model_c.GetChain( 'A'))
      );
      BCL_MessageDbg
      (
        "the connectivity score for a chain without loops with defined coordinates is " +
        util::Format()( connectivity_score_e)
      );
      BCL_ExampleCheckWithinTolerance( connectivity_score_e, 0.0, 0.001);

      // check CalculateChainSSEConnectivityScore with a chain without loops
      BCL_MessageDbg
      (
        "check CalculateChainSSEConnectivityScore with a chain without loops"
      );
      // define min sse sizes
      storage::Map< biol::SSType, size_t> min_sse_sizes
      (
        storage::Map< biol::SSType, size_t>::Create
        (
          std::pair< biol::SSType, size_t>( biol::GetSSTypes().HELIX, 0),
          std::pair< biol::SSType, size_t>( biol::GetSSTypes().STRAND, 0),
          std::pair< biol::SSType, size_t>( biol::GetSSTypes().COIL, 999)
        )
      );

      // get protein model
      assemble::ProteinModel protein_model_d
      (
        Proteins::GetModel
        (
          AddExampleInputPathToFilename( e_Biology, "1ubi_ideal_model.pdb"),
          biol::GetAAClasses().e_AABackBone,
          min_sse_sizes
        )
      );
      // calculate the connectivity
      const double connectivity_score_f
      (
        score_neighbors( *protein_model_d.GetChain( 'A'))
      );
      BCL_MessageDbg
      (
        "the connectivity score for a chain without loops is " +
        util::Format()( connectivity_score_f)
      );
      BCL_ExampleCheckWithinTolerance( connectivity_score_f, 0.0, 0.001);

    ///////////////
    // operators //
    ///////////////

      // test operator
      BCL_MessageDbg( "check operator");
      // calculate the connectivity
      const double connectivity_score_g
      (
        score_neighbors( protein_model_b)
      );
      BCL_MessageDbg
      (
        "the connectivity score for disconnected chain is " + util::Format()( connectivity_score_g)
      );
      BCL_ExampleCheckWithinTolerance( connectivity_score_g, expected_disconnected_score, 0.001);

    //////////////////////
    // input and output //
    //////////////////////

      // test Write and read
      WriteBCLObject( param_constr);
      score::SSEPairConnectivity read_object;
      ReadBCLObject( read_object);

      // make sure the scheme is correct
      BCL_ExampleCheck( param_constr.GetScheme(), scheme);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreSSEPairConnectivity

  const ExampleClass::EnumType ExampleScoreSSEPairConnectivity::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreSSEPairConnectivity())
  );

} // namespace bcl
