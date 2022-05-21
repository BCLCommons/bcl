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
#include "score/bcl_score_sse_predictions.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_locator_sse.h"
#include "biol/bcl_biol_membrane.h"
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_protein_model_sse.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_sse_predictions.cpp
  //! @brief TODO: add an detailed description for this example
  //!
  //! @author woetzen
  //! @date Jul 1, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreSSEPredictions :
    public ExampleInterface
  {
  public:

    ExampleScoreSSEPredictions *Clone() const
    {
      return new ExampleScoreSSEPredictions( *this);
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

    /////////////////
    // preparation //
    /////////////////

      // minimum sse sizes
      storage::Map< biol::SSType, size_t> min_sse_sizes;
      min_sse_sizes[ biol::GetSSTypes().HELIX] = 9;
      min_sse_sizes[ biol::GetSSTypes().STRAND] = 5;

      // initialize pdb filename to be read
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1IE9_shifted.pdb"));

      // create factory and create the model
      pdb::Factory factory( biol::GetAAClasses().e_AABackBone);
      assemble::ProteinModel native_model( factory.ProteinModelFromPDBFilename( pdb_filename, min_sse_sizes));

      // make a copy of the ShPtr to the sequence of chain A
      util::ShPtr< biol::AASequence> sp_sequence( native_model.GetChain( 'A')->GetSequence());

      //create set of SS methods to use for score
      const sspred::Method ss_method( sspred::GetMethods().e_JUFO);
      //storage::Set< sspred::Method> ss_methods( sspred::GetMethods().e_JUFO, sspred::GetMethods().e_PSIPRED, sspred::GetMethods().e_SAM);

      sspred::MethodHandler::ReadPredictionsForProteinModel
      (
        storage::Set< sspred::Method>( ss_method), native_model, "1IE9", AddExampleInputPathToFilename( e_Biology, "")
      );

      // create the locator that finds helix 58-79
      assemble::LocatorSSE helix_58_79_locator( 'A', 58, 79);
      // now locate the SSEs from the protein model
      util::SiPtr< const assemble::SSE> sp_helix_58_79( helix_58_79_locator.Locate( native_model));

      // make sure the sse was correctly located
      BCL_ExampleIndirectAssert
      (
        sp_helix_58_79.IsDefined(),
        true,
        "Helix 58-79 should be located in 1IE9 pdb " + util::Format()( helix_58_79_locator)
      );

      const biol::Membrane sp_membrane;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // check default constructor
      score::SSEPredictions default_score;

      // construct SSEPredictionAverage object from set of ss methods
      score::SSEPredictions confidence_score( ss_method);

      score::ProteinModelSSE score_model( util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > >( confidence_score.Clone()), false);

      BCL_MessageStd( "scoring helix 58-79");
      const double total_score( -19.0543);
      //const double total_score( score_model( native_model));
      confidence_score.WriteDetailedSchemeAndValues( *sp_helix_58_79, sp_membrane, std::cout);
      // print out scores
      BCL_ExampleCheckWithinTolerance( confidence_score( *sp_helix_58_79, sp_membrane).First(), -19.0543, 1.0e-4);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      // create subsequences to prepend and append
      util::ShPtr< biol::AABase> helix_aa_to_prepend( sp_sequence->GetAA( 56));
      util::ShPtr< biol::AABase> helix_aa_to_append( sp_sequence->GetAA( 79));
      biol::AASequence helix_seq_to_prepend( sp_sequence->SubSequence( 50, 6));

      // prepend amino acid
      BCL_MessageStd( "Prepending amino acid 57 to SSEs");
      // make a copy of the helix
      util::ShPtr< assemble::SSE> helix_57_79( sp_helix_58_79->Clone());
      helix_57_79->Prepend( *helix_aa_to_prepend);
      // insert this fragment into protein model
      native_model.ReplaceWithOverlapping( helix_57_79);

      // make sure that prepending was done correctly
      //BCL_MessageStd( "helix_58_79 " + util::Format()( helix_58_79->GetIdentification()));
      //BCL_MessageStd( "helix_57_79 " + util::Format()( helix_57_79->GetIdentification()));

      BCL_MessageStd( "scoring helix 57-79");
      //const double total_score_after_prepend( score_model( native_model));
      confidence_score.WriteDetailedSchemeAndValues( *helix_57_79, sp_membrane, std::cout);
      // print out scores
      BCL_ExampleIndirectCheckWithinTolerance
      (
        confidence_score( *helix_57_79, sp_membrane).First(),
        -19.8284,
        1.0e-4,
        "Prepend modification of score"
      );

      // prepend a bunch of amino acids
      BCL_MessageStd( "Prepending amino acids 51-56 to SSEs");
      // make a copy of the helix
      util::ShPtr< assemble::SSE> helix_51_79( helix_57_79->Clone());
      helix_51_79->PrependSequence( helix_seq_to_prepend);
      // insert this fragment into protein model
      native_model.ReplaceWithOverlapping( helix_51_79);

      // make sure that prepending was done correctly
      BCL_MessageStd( "helix_51_79 " + util::Format()( helix_51_79->GetIdentification()));
      BCL_MessageStd( "scoring helix 51-79");
      //const double total_score_after_prepend_sequence( score_model( native_model));
      confidence_score.WriteDetailedSchemeAndValues( *helix_51_79, sp_membrane, std::cout);
      // print out scores
      BCL_ExampleIndirectCheckWithinTolerance
      (
        confidence_score( *helix_51_79, sp_membrane).First(),
        -21.6103,
        1.0e-4,
        " total_score_sse_after_prepend_sequence"
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write object
      WriteBCLObject( confidence_score);
      score::SSEPredictions confidence_score_read( sspred::GetMethods().e_Undefined);
      ReadBCLObject( confidence_score_read);

      // read and written object should calculate the same score
      BCL_ExampleIndirectCheck
      (
        confidence_score( *helix_51_79, sp_membrane), confidence_score_read( *helix_51_79, sp_membrane),
        "score for helix of written and read object"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreSSEPredictions

  const ExampleClass::EnumType ExampleScoreSSEPredictions::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreSSEPredictions())
  );

} // namespace bcl
