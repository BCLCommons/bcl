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
#include "release/bcl_app_optimize.h"

// includes from bcl - sorted alphabetically
#include "example_application_example_helper.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_pool.h"
#include "fold/bcl_fold.h"
#include "io/bcl_io_file.h"
#include "mc/bcl_mc_optimization_mcm.h"
#include "opti/bcl_opti_ensemble_filter.h"
#include "opti/bcl_opti_ensemble_node.h"
#include "opti/bcl_opti_pipeline.h"
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_protein_model_completeness.h"

// external includes - sorted alphabetically

namespace bcl
{

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_app_optimize.cpp
  //! @brief this example tests the implementation of the application Optimize, which optimizes protein ensembles
  //! using a user-defined optimization procedure.
  //!
  //! @author fischea
  //! @date Oct 31, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleAppOptimize :
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
    //! @return pointer to a new ExampleAppOptimize
    ExampleAppOptimize *Clone() const
    {
      return new ExampleAppOptimize( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create an application instance for testing
      const app::ApplicationType app_optimize( "Optimize");
      BCL_ExampleAssert( app_optimize.IsDefined(), true);

    ////////////////
    // operations //
    ////////////////

      // application helper for calling the object with command line arguments
      ApplicationExampleHelper optimize_helper( app_optimize);

    ///////////////////////
    // loop construction //
    ///////////////////////

      // set arguments for this application test
      // const std::string loop_input_pdb_list_filename( AddExampleInputPathToFilename( e_Mc, "loop_hash_test_set.ls"));
      // const std::string loop_optimizer_filename( AddExampleInputPathToFilename( e_Opti, "loop_hash_ccd_optimizer"));
      // const std::string loop_output_prefix
      // (
      //   AddExampleOutputPathToFilename( fold::GetNamespaceIdentifier(), "app_optimizer_loop_test_")
      // );
      // const size_t number_models( 1);

      // // test sampling five models for the test PDB
      // optimize_helper.SetFlag( "pdb_list", loop_input_pdb_list_filename);
      // optimize_helper.SetFlag( "output_prefix", loop_output_prefix);
      // optimize_helper.SetFlag( "optimizer", loop_optimizer_filename);

      // // check the command line and run the application
      // BCL_ExampleAssert( optimize_helper.CheckCommandString( true), true);
      // BCL_ExampleAssert( optimize_helper.RunCommand(), 0);

      // // check correctness of the created models
      // const pdb::Factory factory( biol::GetAAClasses().e_AABackBone);
      // score::ProteinModelCompleteness compl_score( true);
      // for( size_t i( 0); i < number_models; ++i)
      // {
      //   // check the current model for completeness
      //   const std::string model_filename( loop_output_prefix + util::Format()( i) + ".pdb");
      //   const assemble::ProteinModel model( factory.ProteinModelFromPDBFilename( model_filename));
      //   const double completeness( compl_score( model));
      //   BCL_ExampleCheckWithinTolerance( completeness, -1.0, 0.01);

      //   // remove output file
      //   remove( model_filename.c_str());
      // }

    /////////////////////////////////
    // protein ensemble prediction //
    /////////////////////////////////

      // create a start model and the reference structure
      const pdb::Factory factory( biol::GetAAClasses().e_AABackBone);
      const std::string pdb_path_native( AddExampleInputPathToFilename( e_Fold, "1x91.pdb"));
      const assemble::ProteinModel model_native( factory.ProteinModelFromPDBFilename( pdb_path_native));
      const util::ShPtr< assemble::SSEPool> sp_sse_pool( new assemble::SSEPool( model_native.GetSSEs()));
      const std::string pdb_path( AddExampleInputPathToFilename( e_Fold, "1x91_no_coordinates.pdb"));
      assemble::ProteinModel model( factory.ProteinModelFromPDBFilename( pdb_path));
      util::ShPtr< assemble::ProteinModelData> sp_model_data( new assemble::ProteinModelData( *model.GetProteinModelData()));
      sp_model_data->Insert( assemble::ProteinModelData::e_Pool, sp_sse_pool);
      model = assemble::ProteinModel( model.GetEmptyChains());
      model.SetProteinModelData( sp_model_data);
      assemble::Ensemble< assemble::ProteinModel> ensemble;
      ensemble.AddElement( model);

      // read in the assembly optimizer
      mc::OptimizationMCM optimizer_assembly_1;
      const std::string optimizer_assembly_1_filename( AddExampleInputPathToFilename( e_Mc, "optimizer_fold_assembly"));
      io::IFStream optimizer_assembly_1_file;
      BCL_ExampleMustOpenInputFile( optimizer_assembly_1_file, optimizer_assembly_1_filename);
      optimizer_assembly_1_file >> optimizer_assembly_1;
      io::File::CloseClearFStream( optimizer_assembly_1_file);
      opti::EnsembleNode< assemble::ProteinModel> assembly_1( optimizer_assembly_1, 5);

      // read in the completeness filter
      opti::EnsembleFilter ensemble_filter_completeness;
      const std::string ensemble_filter_completeness_filename( AddExampleInputPathToFilename( e_Opti, "ensemble_filter_completeness"));
      io::IFStream ensemble_filter_completeness_file;
      BCL_ExampleMustOpenInputFile( ensemble_filter_completeness_file, ensemble_filter_completeness_filename);
      ensemble_filter_completeness_file >> ensemble_filter_completeness;
      io::File::CloseClearFStream( ensemble_filter_completeness_file);

      // create the optimization pipeline
      opti::Pipeline< assemble::Ensemble< assemble::ProteinModel> > pipeline;
      pipeline.AppendModule( assembly_1);
      pipeline.AppendModule( ensemble_filter_completeness);

      // conduct the optimization
      pipeline( ensemble);

      return 0;
    }

  }; // class ExampleAppOptimize

  //! single instance of this class
  const ExampleClass::EnumType ExampleAppOptimize::s_Instance
  (
     GetExamples().AddEnum( ExampleAppOptimize())
  );

} // namespace bcl
