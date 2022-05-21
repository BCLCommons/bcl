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
#include "mc/bcl_mc_mutate_loop_fragment_add.h"

// includes from bcl - sorted alphabetically
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_protein_model_completeness.h"

// external includes - sorted alphabetically

namespace bcl
{

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_mc_mutate_loop_fragment_add.cpp
  //! @brief this example tests the implementation of the mutate that adds a loop fragment to a protein model
  //!
  //! @author fischea
  //! @date Mar 1, 2017
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleMcMutateLoopFragmentAdd :
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
    //! @return pointer to a new ExampleMcMutateLoopFragmentAdd
    ExampleMcMutateLoopFragmentAdd *Clone() const
    {
      return new ExampleMcMutateLoopFragmentAdd( *this);
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

      // define which template library to use
      const std::string library_filename( AddExampleInputPathToFilename( e_Fold, "loop_library.ls"));

      // read in a protein model without loops for testing
      const std::string model_filename( AddExampleInputPathToFilename( e_Fold, "1x91_long_loop.pdb"));
      const pdb::Factory factory( biol::GetAAClasses().e_AABackBone);
      assemble::ProteinModel model( factory.ProteinModelFromPDBFilename( model_filename));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create the mutate
      const mc::MutateLoopFragmentAdd mutate( library_filename);

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( mutate.GetClassIdentifier(), ( GetStaticClassName< mc::MutateLoopFragmentAdd>()));

    ///////////////
    // operators //
    ///////////////

      // scoring function to evaluate if coordinates have been added to the model
      score::ProteinModelCompleteness score_function( true);

      // add two fragments to the protein model
      const double score_0( score_function( model));
      math::MutateResult< assemble::ProteinModel> mutate_result( mutate( model));
      const assemble::ProteinModel &model_1( *mutate_result.GetArgument());
      mutate_result = mutate( *mutate_result.GetArgument());
      const assemble::ProteinModel &model_2( *mutate_result.GetArgument());
      const double score_2( score_function( model_2));
      BCL_ExampleCheck( score_2 < score_0, true);

      return 0;
    }

  }; // class ExampleMcMutateLoopFragmentAdd

  //! single instance of this class
  const ExampleClass::EnumType ExampleMcMutateLoopFragmentAdd::s_Instance
  (
    GetExamples().AddEnum( ExampleMcMutateLoopFragmentAdd())
  );

} // namespace bcl
