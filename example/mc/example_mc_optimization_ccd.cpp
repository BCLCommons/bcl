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
#include "mc/bcl_mc_optimization_ccd.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "fold/bcl_fold_protocol_loop_close.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_template_instantiations.h"
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_protein_model_completeness.h"

// external includes - sorted alphabetically

namespace bcl
{

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_mc_optimization_ccd.cpp
  //! @brief this example tests the implementation of the optimizer using cyclic coordinate descent to construct
  //! missing loop regions in protein models.
  //!
  //! @author fischea
  //! @date Oct 27, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleMcOptimizationCCD :
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
    //! @return pointer to a new ExampleMcOptimizationCCD
    ExampleMcOptimizationCCD *Clone() const
    {
      return new ExampleMcOptimizationCCD( *this);
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
      mc::OptimizationCCD default_optimizer;

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( default_optimizer.GetClassIdentifier(), ( GetStaticClassName< mc::OptimizationCCD>()));

    ///////////////////////////
    // CCD loop construction //
    ///////////////////////////

      // create the optimizer
      mc::OptimizationCCD optimizer_ccd;
      const std::string optimizer_ccd_file_name( AddExampleInputPathToFilename( e_Mc, "optimizer_loop_ccd"));
      io::IFStream optimizer_ccd_file;
      BCL_ExampleMustOpenInputFile( optimizer_ccd_file, optimizer_ccd_file_name);
      optimizer_ccd_file >> optimizer_ccd;
      io::File::CloseClearFStream( optimizer_ccd_file);

      // construct missing loop regions in the model and check for correctness
      optimizer_ccd( model_without_loops);
      const double completeness_ccd_final( completeness_estimator( model_without_loops));
      BCL_ExampleCheckWithinTolerance( completeness_ccd_final, -1.0, 0.01);

      return 0;
    }

  }; // class ExampleMcOptimizationCCD

  //! single instance of this class
  const ExampleClass::EnumType ExampleMcOptimizationCCD::s_Instance
  (
    GetExamples().AddEnum( ExampleMcOptimizationCCD())
  );

} // namespace bcl
