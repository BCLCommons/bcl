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
#include "mc/bcl_mc_mutate_loop_replace.h"

// includes from bcl - sorted alphabetically
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_protein_model_completeness.h"

// external includes - sorted alphabetically

namespace bcl
{

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_mc_mutate_loop_replace.cpp
  //! @brief this example tests the implementation of the mutate that adds a loop to a protein model
  //!
  //! @author fischea
  //! @date Mar 9, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleMcMutateLoopReplace :
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

    //! @brief clone function
    //! @return pointer to a new ExampleMcMutateLoopReplace
    ExampleMcMutateLoopReplace *Clone() const
    {
      return new ExampleMcMutateLoopReplace( *this);
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

    //! @brief run routine
    //! @detail this is performing the execution of the example
    int Run() const
    {

    //////////////////////
    // data preparation //
    //////////////////////

      // read in the loop template library for testing
      const std::string library_filename( AddExampleInputPathToFilename( e_Fold, "loop_library.ls"));

      // read in a protein model with defined loops for testing
      const std::string model_filename( AddExampleInputPathToFilename( e_Fold, "1x91.pdb"));
      const pdb::Factory factory( biol::GetAAClasses().e_AABackBone);
      assemble::ProteinModel model( factory.ProteinModelFromPDBFilename( model_filename));

      // create a scoring function to quantify model completeness
      score::ProteinModelCompleteness compl_score;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create the mutate
      const mc::MutateLoopReplace mutate( library_filename);

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( mutate.GetClassIdentifier(), ( GetStaticClassName< mc::MutateLoopReplace>()));

    ///////////////
    // operators //
    ///////////////

      // compute the initial completeness
      double last_completeness( compl_score( model));

      // construct all five missing, non-terminal loops
      for( size_t i( 0); i < 10; ++i)
      {
        // add a loop to the protein model
        const math::MutateResult< assemble::ProteinModel> mutate_result( mutate( model));
        const assemble::ProteinModel &result_model( *mutate_result.GetArgument());

        // compute the completeness of the resulting model
        const double completeness( compl_score( result_model));

        // check if completeness of the model has increased (if coordinates were added)
        BCL_Example_Check( completeness == last_completeness, "Loop was not inserted into the model.");

        // store the new results
        model = result_model;
        last_completeness = completeness;
      }

      return 0;
    }

  }; // class ExampleMcMutateLoopReplace

  //! single instance of this class
  const ExampleClass::EnumType ExampleMcMutateLoopReplace::s_Instance
  (
     GetExamples().AddEnum( ExampleMcMutateLoopReplace())
  );

} // namespace bcl
