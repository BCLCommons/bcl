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
#include "opencl/bcl_opencl_approximator_sequential_minimial_optimization.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_criterion_number_iterations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_sequential_minimal_optimization.cpp
  //!
  //! @author mendenjl
  //! @date Sep 04, 2013
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclApproximatorSequentialMinimialOptimization :
    public ExampleInterface
  {
  public:

    ExampleOpenclApproximatorSequentialMinimialOptimization *Clone() const
    {
      return new ExampleOpenclApproximatorSequentialMinimialOptimization( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // do not try to run opencl commands if no queue was found
      if( !opencl::GetTools().HasCommandQueues())
      {
        return 1;
      }

      // dataset for training data
      util::ShPtr< storage::Vector< storage::VectorND< 2, linal::Vector< float> > > > features_results_setup
      (
        new storage::Vector< storage::VectorND< 2, linal::Vector< float> > >()
      );

      // number of training vectors
      static const size_t s_number_values( 16);

      //initializing
      storage::VectorND< 2, linal::Vector< float> > features_results[ s_number_values] =
      {
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 0.0), linal::MakeVector< float>( -1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 5.0), linal::MakeVector< float>( 2.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 5.0, 0.0), linal::MakeVector< float>( 2.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 5.0, 5.0), linal::MakeVector< float>( -1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 0.0), linal::MakeVector< float>( -1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 2.0), linal::MakeVector< float>( 2.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 2.0, 0.0), linal::MakeVector< float>( 2.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 2.0, 2.0), linal::MakeVector< float>( -1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 0.0), linal::MakeVector< float>( -1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 3.0), linal::MakeVector< float>( 2.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 3.0, 0.0), linal::MakeVector< float>( 2.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 3.0, 3.0), linal::MakeVector< float>( -1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 0.0), linal::MakeVector< float>( -1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 4.0), linal::MakeVector< float>( 2.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 4.0, 0.0), linal::MakeVector< float>( 2.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 4.0, 4.0), linal::MakeVector< float>( -1.0))
      };

      // filling training dataset
      for( size_t index( 0); index < s_number_values; index++)
      {
        features_results_setup->PushBack( features_results[ index]);
      }

      // create feature result data set
      util::ShPtr< descriptor::Dataset> features_results_frds
      (
        new descriptor::Dataset( *features_results_setup)
      );

      opencl::CommandQueue queue( opencl::GetTools().GetFirstCommandQueue());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // support vector machine model initialized with kernel
      util::ShPtr< model::Interface> svr_model( new opencl::SupportVectorMachine( queue));

      // construct iterate for svm training from properties
      opencl::ApproximatorSequentialMinimialOptimization smo
      (
        float( 1),       // parameter C
        svr_model,              // svm model
        features_results_frds,  // rescaled training dataset
        2,                       // number of internal iterations
        queue
      );

      smo.SetGamma( 0.4);

      util::ShPtr< model::ApproximatorBase> sp_iterate( smo.Clone());

      // objective function
      util::ShPtr< model::ObjectiveFunctionWrapper> sp_objective_function
      (
        new model::ObjectiveFunctionWrapper
        (
          features_results_frds,
          util::ObjectDataLabel( "Enrichment(cutoff=0.5, oversampling factor=1, FPR cutoff=0.01, parity=1)")
        )
      );

      sp_objective_function->SetData( features_results_frds, sp_iterate->GetRescaleFeatureDataSet());

      // set objective function in iterate
      sp_iterate->SetObjectiveFunction( sp_objective_function);

      // terminate criterium
      sp_iterate->SetCriterion( opti::CriterionNumberIterations< util::ShPtr< model::Interface>, float>( 3));

      BCL_MessageStd( "sp_iterate->GetCurrentModel()..");

      util::ShPtr< model::Interface> current_model( sp_iterate->GetCurrentModel());

      BCL_MessageStd( "sp_iterate->Approximate()");

      // call minimizer until criterion met
      sp_iterate->Approximate();
      const storage::Pair< util::ShPtr< model::Interface>, float> result( *sp_iterate->GetCurrentApproximation());

      BCL_ExampleCheckWithinTolerance( double( result.Second()), double( 1.17425), 0.01);

      BCL_MessageStd( "irp_approximator->Approximate.. Done");

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclApproximatorSequentialMinimialOptimization

  const ExampleClass::EnumType ExampleOpenclApproximatorSequentialMinimialOptimization::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclApproximatorSequentialMinimialOptimization())
  );

} // namespace bcl
