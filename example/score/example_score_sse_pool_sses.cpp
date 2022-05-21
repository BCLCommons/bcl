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
#include "score/bcl_score_sse_pool_sses.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "biol/bcl_biol_membrane.h"
#include "score/bcl_score_phi_psi.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_sse_pool_sses.cpp
  //! @brief test the score::SSPoolSSEs with the PhiPsi score of sses on a protein model's SSEs
  //!
  //! @author woetzen
  //! @date Jun 16, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreSSEPoolSSEs :
    public ExampleInterface
  {
  public:

    ExampleScoreSSEPoolSSEs *Clone() const
    {
      return new ExampleScoreSSEPoolSSEs( *this);
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
      min_sse_sizes[ biol::GetSSTypes().HELIX] = 0;
      min_sse_sizes[ biol::GetSSTypes().STRAND] = 0;
      min_sse_sizes[ biol::GetSSTypes().COIL] = 0;

      // initialize pdb filename to be read
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"));

      // create the model and the pool
      assemble::ProteinModel native_model( Proteins::GetModel( pdb_filename, biol::GetAAClasses().e_AABackBone, min_sse_sizes));
      sspred::MethodHandler::ReadPredictionsForProteinModel( storage::Set< sspred::Method>( sspred::GetMethods().e_JUFO), native_model, "1ubi", AddExampleInputPathToFilename( e_Biology, ""));

      const assemble::SSEPool pool( native_model.GetSSEs(), false);

      // score for individual sse in the pool
      const util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > > sp_sse_score_phispi
      (
        new score::PhiPsi( score::PhiPsi::GetDefaultScheme())
      );

      const biol::Membrane sp_membrane;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct from score and normalization
      score::SSEPoolSSEs score_pool_no_normalization_phipsi( sp_sse_score_phispi, false, false);

    /////////////////
    // data access //
    /////////////////

      // class identifier
      BCL_ExampleCheck( score_pool_no_normalization_phipsi.GetClassIdentifier(), GetStaticClassName< score::SSEPoolSSEs>());

      // scheme
      BCL_ExampleCheck( score_pool_no_normalization_phipsi.GetScheme(), sp_sse_score_phispi->GetScheme());

    ///////////////
    // operators //
    ///////////////

      // score the pool
      const double pool_score_phipsi( score_pool_no_normalization_phipsi( pool, sp_membrane));
      const double pool_score_expected_phipsi( -187.737);

      BCL_ExampleCheckWithinTolerance( pool_score_expected_phipsi, pool_score_phipsi, 1.0e-4);

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write detailed scheme and values
      std::stringstream sstream;

      // test read write
      WriteBCLObject( score_pool_no_normalization_phipsi);
      score::SSEPoolSSEs score_pool_read( util::ShPtr< math::BinaryFunctionInterface< assemble::SSE, biol::Membrane, storage::Pair< double, size_t> > >(), true, true);
      ReadBCLObject( score_pool_read);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreSSEPoolSSEs

  const ExampleClass::EnumType ExampleScoreSSEPoolSSEs::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreSSEPoolSSEs())
  );

} // namespace bcl
