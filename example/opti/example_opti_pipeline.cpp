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
#include "opti/bcl_opti_pipeline.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_file.h"
#include "mc/bcl_mc_optimization_ccd.h"
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_protein_model_completeness.h"

// external includes - sorted alphabetically

namespace bcl
{

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_pipeline.cpp
  //! @brief tests the implementation of Pipeline
  //!
  //! @author fischea
  //! @date Oct 24, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleOptiPipeline :
    public ExampleInterface
  {

  //////////
  // data //
  //////////

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleOptiPipeline
    ExampleOptiPipeline *Clone() const
    {
      return new ExampleOptiPipeline( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      opti::Pipeline< assemble::ProteinModel> default_pipeline;

    /////////////////
    // data access //
    /////////////////

      // test getter for class name identifier
      BCL_ExampleCheck
      (
        default_pipeline.GetClassIdentifier(), ( GetStaticClassName< opti::Pipeline< assemble::ProteinModel> >())
      );

    ///////////////////////
    // loop construction //
    ///////////////////////

      // create a test model that is lacking loop regions
      const std::string model_without_loops_filename( AddExampleInputPathToFilename( e_Fold, "1x91_no_loops.pdb"));
      const pdb::Factory factory( biol::GetAAClasses().e_AABackBone);
      assemble::ProteinModel model_without_loops( factory.ProteinModelFromPDBFilename( model_without_loops_filename));

      // create the pipeline module for loop hashing
      mc::OptimizationMCM optimizer_hash;
      const std::string optimizer_hash_file_name( AddExampleInputPathToFilename( e_Mc, "optimizer_loop_hash"));
      io::IFStream optimizer_hash_file;
      BCL_ExampleMustOpenInputFile( optimizer_hash_file, optimizer_hash_file_name);
      optimizer_hash_file >> optimizer_hash;
      io::File::CloseClearFStream( optimizer_hash_file);

      // create the pipeline module for cyclic coordinate descent
      mc::OptimizationCCD optimizer_ccd;
      const std::string optimizer_ccd_file_name( AddExampleInputPathToFilename( e_Mc, "optimizer_loop_ccd"));
      io::IFStream optimizer_ccd_file;
      BCL_ExampleMustOpenInputFile( optimizer_ccd_file, optimizer_ccd_file_name);
      optimizer_ccd_file >> optimizer_ccd;
      io::File::CloseClearFStream( optimizer_ccd_file);

      // create a pipeline for loop construction consisting of loop hashing and cyclic coordinate descent
      opti::Pipeline< assemble::ProteinModel> loop_pipeline;
      loop_pipeline.AppendModule( optimizer_hash);
      loop_pipeline.AppendModule( optimizer_ccd);

      // apply the pipeline to the test model and check of the missing loop regions have been constructed
      const score::ProteinModelCompleteness completeness_estimator;
      const double loop_completeness_initial( completeness_estimator( model_without_loops));
      loop_pipeline( model_without_loops);
      const double loop_completeness( completeness_estimator( model_without_loops));
      BCL_Example_Check( loop_completeness_initial > loop_completeness, "No loop regions were constructed.");

      return 0;
    }

  }; // class ExampleOptiPipeline

  //! single instance of this class
  const ExampleClass::EnumType ExampleOptiPipeline::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiPipeline())
  );

} // namespace bcl
