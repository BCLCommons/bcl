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
#include "opti/bcl_opti_ensemble_filter.h"

// includes from bcl - sorted alphabetically
#include "pdb/bcl_pdb_factory.h"
#include "score/bcl_score_protein_model_completeness.h"

// external includes - sorted alphabetically

namespace bcl
{

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_ensemble_filter.cpp
  //! @brief this example tests the implementation of the ensemble filter.
  //!
  //! @author fischea
  //! @date Apr 18, 2017
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleOptiEnsembleFilter :
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
    //! @return pointer to a new ExampleOptiEnsembleFilter
    ExampleOptiEnsembleFilter *Clone() const
    {
      return new ExampleOptiEnsembleFilter( *this);
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

      // read in a few protein models with different completeness values and create an ensemble
      const pdb::Factory factory( biol::GetAAClasses().e_AABackBone);
      const std::string model_1_filename( AddExampleInputPathToFilename( e_Fold, "1x91_no_loops.pdb"));
      assemble::ProteinModel model_1( factory.ProteinModelFromPDBFilename( model_1_filename));
      const std::string model_2_filename( AddExampleInputPathToFilename( e_Fold, "1x91_partial_loops.pdb"));
      assemble::ProteinModel model_2( factory.ProteinModelFromPDBFilename( model_2_filename));
      const std::string model_3_filename( AddExampleInputPathToFilename( e_Fold, "1x91.pdb"));
      assemble::ProteinModel model_3( factory.ProteinModelFromPDBFilename( model_3_filename));
      assemble::Ensemble< assemble::ProteinModel> ensemble;
      ensemble.AddElement( model_1);
      ensemble.AddElement( model_2);
      ensemble.AddElement( model_3);

      // create a scoring function
      score::ProteinModelCompleteness score_function;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      const double keep_percentage( 0.67);
      opti::EnsembleFilter filter( score_function, keep_percentage);

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( filter.GetClassIdentifier(), ( GetStaticClassName< opti::EnsembleFilter>()));

    ////////////////
    // operations //
    ////////////////

      // filter the ensemble based on completeness
      const size_t expected_size( ensemble.GetSize() * keep_percentage);
      filter( ensemble);

      // make sure the correct number of models has been filtered output
      BCL_ExampleCheck( ensemble.GetSize(), expected_size);

      return 0;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  }; // class ExampleOptiEnsembleFilter

  //! single instance of this class
  const ExampleClass::EnumType ExampleOptiEnsembleFilter::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiEnsembleFilter())
  );

} // namespace bcl
