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
#include "score/bcl_score_environment_predictions.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_locator_sse.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_membrane.h"
#include "sspred/bcl_sspred_method_handler.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_environment_predictions.cpp
  //!
  //! @author weinerbe
  //! @date Jun 20, 2012
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreEnvironmentPredictions :
    public ExampleInterface
  {
  public:

    ExampleScoreEnvironmentPredictions *Clone() const
    {
      return new ExampleScoreEnvironmentPredictions( *this);
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
      assemble::ProteinModel protein_model
      (
        Proteins::GetModel( AddExampleInputPathToFilename( e_Biology, "2K73A.pdb"))
      );

      // read in JUFO9D
      const sspred::Method ss_method( sspred::GetMethods().e_JUFO9D);
      sspred::MethodHandler::ReadPredictionsForProteinModel
      (
        storage::Set< sspred::Method>( ss_method), protein_model, "2K73", AddExampleInputPathToFilename( e_Biology, "")
      );

      // get a tm helix
      const assemble::LocatorSSE tmh_locator( 'A', 41, 62);
      const assemble::SSE tm_helix( *( tmh_locator.Locate( protein_model)));

      // create membrane
      const biol::Membrane membrane;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      score::EnvironmentPredictions def_construct;

      // test JUFO9D constructor
      const score::EnvironmentPredictions jufo_construct( ss_method);

    ////////////////
    // operations //
    ////////////////

      // test () operator
      BCL_ExampleIndirectCheck( def_construct( tm_helix, membrane).First(), 0.0, "() operator with no method specified");
      const double expected_score( -5.38411);
      const double calculated_score( jufo_construct( tm_helix, membrane).First());
      BCL_MessageStd( "Score: " + util::Format()( calculated_score));
      BCL_ExampleIndirectCheckWithinTolerance
      (
        calculated_score,
        expected_score,
        0.001,
        "() operator with JUFO9D method specified"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test WriteDetailedSchemeAndValues
      jufo_construct.WriteDetailedSchemeAndValues( tm_helix, membrane, util::GetLogger());

      // write object
      WriteBCLObject( jufo_construct);
      ReadBCLObject( def_construct);

      // read and written object should calculate the same score
      BCL_ExampleIndirectCheckWithinTolerance
      (
        def_construct( tm_helix, membrane).First(),
        expected_score,
        0.001,
        "score for read-in object"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreEnvironmentPredictions

  const ExampleClass::EnumType ExampleScoreEnvironmentPredictions::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreEnvironmentPredictions())
  );

} // namespace bcl
