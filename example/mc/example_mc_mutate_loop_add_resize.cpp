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
#include "mc/bcl_mc_mutate_loop_add_resize.h"

// includes from bcl - sorted alphabetically
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_protein_model_defined_loops.h"

// external includes - sorted alphabetically

namespace bcl
{

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_mc_mutate_loop_add_resize.cpp
  //! @brief this example tests the implementation of the mutate that adds a loop to a protein model with the
  //! possibility of resizing the anchor SSEs.
  //!
  //! @author fischea
  //! @date Mar 9, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleMcMutateLoopAddResize :
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
    //! @return pointer to a new ExampleMcMutateLoopAddResize
    ExampleMcMutateLoopAddResize *Clone() const
    {
      return new ExampleMcMutateLoopAddResize( *this);
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

      // read in a protein model without loops for testing
      const std::string model_filename( AddExampleInputPathToFilename( e_Fold, "1x91_no_loops_resize.pdb"));
      const pdb::Factory factory( biol::GetAAClasses().e_AABackBone);
      assemble::ProteinModel model( factory.ProteinModelFromPDBFilename( model_filename));

      // create a scoring function to quantify model completeness
      score::ProteinModelDefinedLoops compl_score( true);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create the mutate
      const mc::MutateLoopAddResize mutate( library_filename, mc::MutateLoopAddResize::GetDefaultMinSizes());

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( mutate.GetClassIdentifier(), ( GetStaticClassName< mc::MutateLoopAddResize>()));

    ///////////////
    // operators //
    ///////////////

      // compute the initial completeness
      double initial_completeness( compl_score( model));

      // check if loop regions in the proteins are missing
      BCL_Example_Check( initial_completeness == 5.0, "Model does not have missing loop regions.");

      // construct all five missing, non-terminal loops
      for( size_t i( 0); i < 300; ++i)
      {
        // add a loop to the protein model
        const math::MutateResult< assemble::ProteinModel> mutate_result( mutate( model));
        const assemble::ProteinModel &result_model( *mutate_result.GetArgument());
        model = result_model;
      }

      // compute the completeness of the resulting model
      const double completeness( compl_score( model));

      // check if all loops have been added to the protein model
      BCL_Example_Check( completeness < 5.0, "Model has too many missing loop regions.");

      return 0;
    }

  }; // class ExampleMcMutateLoopAddResize

  //! single instance of this class
  const ExampleClass::EnumType ExampleMcMutateLoopAddResize::s_Instance
  (
     GetExamples().AddEnum( ExampleMcMutateLoopAddResize())
  );

} // namespace bcl
