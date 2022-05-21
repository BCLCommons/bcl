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
#include "mc/bcl_mc_optimization_mcm.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_pick_sse_random.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "find/bcl_find_pick_criteria_wrapper.h"
#include "fold/bcl_fold_mutate_protein_model_sse_add.h"
#include "fold/bcl_fold_mutate_tree.h"
#include "fold/bcl_fold_placement_strand_next_to_sheet.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_mutate_decision_node.h"
#include "math/bcl_math_template_instantiations.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_criterion_unimproved.h"
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_protein_model_completeness.h"

// external includes - sorted alphabetically

namespace bcl
{

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_mc_optimization_mcm.cpp
  //! @brief this example tests the implementation of the optimizer using MCM to optimize protein models.
  //!
  //! @author fischea
  //! @date Aug 10, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleMcOptimizationMCM :
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
    //! @return pointer to a new ExampleMcOptimizationMCM
    ExampleMcOptimizationMCM *Clone() const
    {
      return new ExampleMcOptimizationMCM( *this);
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

      // read in a protein model without loops for testing
      const std::string model_without_loops_filename( AddExampleInputPathToFilename( e_Fold, "1x91_no_loops.pdb"));
      const pdb::Factory factory( biol::GetAAClasses().e_AABackBone);
      assemble::ProteinModel model_without_loops( factory.ProteinModelFromPDBFilename( model_without_loops_filename));

      // evaluating loop completeness
      const score::ProteinModelCompleteness completeness_estimator;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      mc::OptimizationMCM default_optimizer;

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( default_optimizer.GetClassIdentifier(), ( GetStaticClassName< mc::OptimizationMCM>()));

    //////////////////
    // loop hashing //
    //////////////////

      // create the optimizer
      // mc::OptimizationMCM optimizer_hash;
      // const std::string optimizer_hash_file_name( AddExampleInputPathToFilename( e_Mc, "optimizer_loop_hash"));
      // io::IFStream optimizer_hash_file;
      // BCL_ExampleMustOpenInputFile( optimizer_hash_file, optimizer_hash_file_name);
      // optimizer_hash_file >> optimizer_hash;
      // io::File::CloseClearFStream( optimizer_hash_file);

      // // construct missing loop regions in the model and check for correctness
      // const double completeness_hash_initial( completeness_estimator( model_without_loops));
      // optimizer_hash( model_without_loops);
      // const double completeness_hash_final( completeness_estimator( model_without_loops));
      // BCL_Example_Check( completeness_hash_final < completeness_hash_initial, "No loop regions were constructed.");

    /////////////////////
    // de novo folding //
    /////////////////////

      // create a start model and the reference structure
      const std::string pdb_path_native( AddExampleInputPathToFilename( e_Fold, "1x91.pdb"));
      const assemble::ProteinModel model_native( factory.ProteinModelFromPDBFilename( pdb_path_native));
      const util::ShPtr< assemble::SSEPool> sp_sse_pool( new assemble::SSEPool( model_native.GetSSEs()));
      const std::string pdb_path( AddExampleInputPathToFilename( e_Fold, "1x91_no_coordinates.pdb"));
      assemble::ProteinModel model( factory.ProteinModelFromPDBFilename( pdb_path));
      util::ShPtr< assemble::ProteinModelData> sp_model_data( new assemble::ProteinModelData( *model.GetProteinModelData()));
      sp_model_data->Insert( assemble::ProteinModelData::e_Pool, sp_sse_pool);
      model = assemble::ProteinModel( model.GetEmptyChains());
      model.SetProteinModelData( sp_model_data);

      // create the optimizer
      mc::OptimizationMCM optimizer_fold_assembly;
      const std::string optimizer_assembly_file_name( AddExampleInputPathToFilename( e_Mc, "optimizer_fold_assembly"));
      io::IFStream optimizer_assembly_file;
      BCL_ExampleMustOpenInputFile( optimizer_assembly_file, optimizer_assembly_file_name);
      optimizer_assembly_file >> optimizer_fold_assembly;
      io::File::CloseClearFStream( optimizer_assembly_file);

      // de novo predict the protein's tertiary structure
      optimizer_fold_assembly( model);

      return 0;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  private:

    template< typename t_Interface>
    static t_Interface *AddInstance( t_Interface *INSTANCE, const std::string &ALIAS = std::string())
    {
      const std::string &alias( ALIAS.empty() ? INSTANCE->GetAlias() : ALIAS);

      // validate GetSerializer
      BCL_Assert( INSTANCE->GetCompleteSerializer().IsFinalized(), "Instance named " + alias + " has invalid serializer");

      // Try to add the object to object interface
      GetObjectInstances().TryAddInstance( INSTANCE);

      return INSTANCE;
    }

  }; // class ExampleMcOptimizationMCM

  //! single instance of this class
  const ExampleClass::EnumType ExampleMcOptimizationMCM::s_Instance
  (
    GetExamples().AddEnum( ExampleMcOptimizationMCM())
  );

} // namespace bcl
