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
#include "score/bcl_score_protein_model_gap.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "pdb/bcl_pdb_factory.h"

// external includes - sorted alphabetically

namespace bcl
{

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_protein_model_gap.cpp
  //! @brief this example tests the implementation of the class scoring spans with undefined coordinates.
  //!
  //! @author fischea
  //! @date Feb 27, 2017
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleScoreProteinModelGap :
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

    //! @brief returns a pointer to a new ExampleScoreProteinModelGap
    //! @return pointer to a new ExampleScoreProteinModelGap
    ExampleScoreProteinModelGap *Clone() const
    {
      return new ExampleScoreProteinModelGap( *this);
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
      const std::string &scheme( "comp_test");
      const score::ProteinModelGap default_comp_score;
      const score::ProteinModelGap spec_comp_score( scheme);

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_ExampleCheck( default_comp_score.GetClassIdentifier(), GetStaticClassName< score::ProteinModelGap>());

      // check the schemes of the score instances
      BCL_ExampleCheck( score::ProteinModelGap::GetDefaultScheme().compare( default_comp_score.GetScheme()), 0);
      BCL_ExampleCheck( scheme.compare( spec_comp_score.GetScheme()), 0);

    ////////////////
    // operations //
    ////////////////

      // create protein models for testing
      const std::string pdb_filename_1( AddExampleInputPathToFilename( e_Fold, "1x91.pdb"));
      const assemble::ProteinModel model_1( pdb::Factory().ProteinModelFromPDBFilename( pdb_filename_1));
      const std::string pdb_filename_gap_1( AddExampleInputPathToFilename( e_Fold, "1x91_gap_1.pdb"));
      const assemble::ProteinModel model_gap_1( pdb::Factory().ProteinModelFromPDBFilename( pdb_filename_gap_1));
      const std::string pdb_filename_gap_2( AddExampleInputPathToFilename( e_Fold, "1x91_gap_2.pdb"));
      const assemble::ProteinModel model_gap_2( pdb::Factory().ProteinModelFromPDBFilename( pdb_filename_gap_2));

      // compute the gap scores of the models
      const double score_full( default_comp_score( model_1));
      const double score_gap_1( default_comp_score( model_gap_1));
      const double score_gap_2( default_comp_score( model_gap_2));

      // compare the calculated scores to the expected values
      BCL_ExampleCheckWithinTolerance( score_full, 0.0, 0.01);
      BCL_ExampleCheck( score_gap_1 > score_full, true);
      BCL_ExampleCheck( score_gap_2 > score_gap_1, true);

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

  }; // class ExampleScoreProteinModelGap

  //! single instance of this class
  const ExampleClass::EnumType ExampleScoreProteinModelGap::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreProteinModelGap())
  );

} // namespace bcl
