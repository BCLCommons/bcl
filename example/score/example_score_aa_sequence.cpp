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
#include "score/bcl_score_aa_sequence.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_sse.h"
#include "biol/bcl_biol_membrane.h"
#include "score/bcl_score_aa_pair_clash.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_aa_sequence.cpp
  //!
  //! @author karakam, alexanns
  //! @date Jun 14, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreAASequence :
    public ExampleInterface
  {
  public:

    ExampleScoreAASequence *Clone() const
    {
      return new ExampleScoreAASequence( *this);
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

      // initialize pdb filename
      const std::string pdb_filename( AddExampleInputPathToFilename( e_Biology, "1ubi_internal_clash.pdb"));
      // get the protein model
      util::ShPtr< assemble::SSE> sp_loop( Proteins::GetSSE( pdb_filename, 'A', 46, 63));
      util::ShPtr< assemble::SSE> sp_helix( Proteins::GetSSE( pdb_filename, 'A', 23, 34));
      util::ShPtr< assemble::SSE> sp_strand( Proteins::GetSSE( pdb_filename, 'A', 40, 45));

      const biol::Membrane membrane;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct the AAPairClash score first
      util::ShPtr< score::AAPairDistanceInterface> sp_clash_score( new score::AAPairClash( 1.0, 1));

      // construct from a potential and a normalize boolean
      score::AASequence score( sp_clash_score, false);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      BCL_MessageStd( "Scoring loop " + sp_loop->GetIdentification());

      // expected values
      const double expected_score_loop( 1.4658);
      const double score_loop( score( *sp_loop, membrane).First());
      BCL_MessageStd( "score: " + util::Format()( score_loop));
      // compare the expected and calculated score
      BCL_ExampleCheckWithinTolerance( expected_score_loop, score_loop, 0.001);

      BCL_MessageStd( "Scoring helix " + sp_helix->GetIdentification());
      // expected values
      const double expected_score_helix( 0.0);
      const double score_helix( score( *sp_helix, membrane).First());
      BCL_MessageStd( "score: " + util::Format()( score_helix));
      // compare the expected and calculated score
      BCL_ExampleCheckWithinTolerance( expected_score_helix, score_helix, 0.001);

      BCL_MessageStd( "Scoring strand " + sp_strand->GetIdentification());
      // expected values
      const double expected_score_strand( 0.0);
      const double score_strand( score( *sp_strand, membrane).First());
      BCL_MessageStd( "score: " + util::Format()( score_strand));
      // compare the expected and calculated score
      BCL_ExampleCheckWithinTolerance( expected_score_strand, score_strand, 0.001);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreAASequence

  const ExampleClass::EnumType ExampleScoreAASequence::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreAASequence())
  );

} // namespace bcl
