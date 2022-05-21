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
#include "score/bcl_score_protein_model_sse_completeness.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_protein_model_sse_completeness.cpp
  //! @brief this example tests the implementation of the class scoring the completeness of protein models
  //!
  //! @author mendenjl
  //! @date Jan 23, 2017
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleScoreProteinModelSSECompleteness :
    public ExampleInterface
  {

  //////////
  // data //
  //////////

  public:

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief returns a pointer to a new ExampleScoreProteinModelSSECompleteness
    //! @return pointer to a new ExampleScoreProteinModelSSECompleteness
    ExampleScoreProteinModelSSECompleteness *Clone() const
    {
      return new ExampleScoreProteinModelSSECompleteness( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief performs the execution of the example
    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct a scoring object
      score::ProteinModelSSECompleteness default_comp_score;

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck
      (
        default_comp_score.GetClassIdentifier(), GetStaticClassName< score::ProteinModelSSECompleteness>()
      );

    ////////////////
    // operations //
    ////////////////

      // create protein models for testing
      const std::string pdb_filename_1( AddExampleInputPathToFilename( e_Fold, "1x91.pdb"));
      const assemble::ProteinModel model_1( pdb::Factory().ProteinModelFromPDBFilename( pdb_filename_1));
      const std::string pdb_filename_2( AddExampleInputPathToFilename( e_Fold, "1x91_no_loop_coords.pdb"));
      const assemble::ProteinModel model_2( pdb::Factory().ProteinModelFromPDBFilename( pdb_filename_2));

      // score the completeness of the protein models
      const double score_model_1_full( default_comp_score( model_1));
      const double score_model_2_full( default_comp_score( model_2));

      // compare the calculated scores to the expected values
      BCL_ExampleCheckWithinTolerance( score_model_1_full, 0, 0.01);
      BCL_ExampleCheckWithinTolerance( score_model_2_full, 0, 0.01);

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

  }; // class ExampleScoreProteinModelSSECompleteness

  //! single instance of this class
  const ExampleClass::EnumType ExampleScoreProteinModelSSECompleteness::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreProteinModelSSECompleteness())
  );

} // namespace bcl
