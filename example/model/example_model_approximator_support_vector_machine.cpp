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
#include "model/bcl_model_approximator_support_vector_machine.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_template_instantiations.h"
#include "model/bcl_model_support_vector_kernel_rbf.h"
#include "model/bcl_model_support_vector_machine.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_elapsed_time.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_criterion_result_threshold.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_sequential_minimal_optimization.cpp
  //!
  //! @author mendenjl
  //! @date
  //! @remarks status empty
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelApproximatorSupportVectorMachine :
    public ExampleInterface
  {
  public:

    ExampleModelApproximatorSupportVectorMachine *Clone() const
    {
      return new ExampleModelApproximatorSupportVectorMachine( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // kernel object for support vector machine model
      util::Implementation< model::SupportVectorKernelBase> kernel( model::SupportVectorKernelRBF( float( 0.4)));

      // support vector machine model initialized with kernel
      util::ShPtr< model::Interface> svr_model( new model::SupportVectorMachine( kernel));

      // construct iterate for svm training from properties
      util::ShPtr< model::ApproximatorBase> sp_iterate
      (
        new model::ApproximatorSupportVectorMachine
        (
          float( 1),       // parameter C
          svr_model,              // svm model
          features_results_frds,  // rescaled training dataset
          2                       // number of internal iterations
        )
      );

      // objective function
      util::ShPtr< model::ObjectiveFunctionWrapper> sp_objective_function
      (
        new model::ObjectiveFunctionWrapper
        (
          features_results_frds,
          util::ObjectDataLabel( "Enrichment(cutoff=0.5, oversampling factor=1, FPR cutoff=0.01, parity=1)")
        )
      );

      // set objective function in iterate
      sp_iterate->SetObjectiveFunction( sp_objective_function);

      // terminate criterium
      sp_iterate->SetCriterion( opti::CriterionNumberIterations< util::ShPtr< model::Interface>, float>( 2));

      BCL_MessageStd( "irp_approximator->Approximate..");

      util::ShPtr< model::Interface> current_model( sp_iterate->GetCurrentModel());

      // call minimizer until criterion met
      sp_iterate->Approximate();
      const storage::Pair< util::ShPtr< model::Interface>, float> result( *sp_iterate->GetCurrentApproximation());

      // Evaluate FINAL Model
      util::ShPtr< model::Interface> final_model( result.First());

      BCL_MessageStd( "irp_approximator->Approximate.. Done");

      // this is an actual optimization
      // create terminate Criteria
      util::ShPtr< opti::CriterionCombine< util::ShPtr< model::Interface>, float> >
        model_terminate( new opti::CriterionCombine< util::ShPtr< model::Interface>, float>());

      // number iterations, 1000 * 50 internal iterations = 50000
      model_terminate->InsertCriteria( opti::CriterionNumberIterations< util::ShPtr< model::Interface>, float>( 1000));
      model_terminate->InsertCriteria
      (
        opti::CriterionElapsedTime< util::ShPtr< model::Interface>, float>( util::Time( 60 * 60 * 5, 0))
      ); // 5 hours
      model_terminate->InsertCriteria
      (
        opti::CriterionResultThreshold< util::ShPtr< model::Interface>, float>( float( 0.1))
      ); // threshold

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelApproximatorSupportVectorMachine

  const ExampleClass::EnumType ExampleModelApproximatorSupportVectorMachine::s_Instance
  (
    GetExamples().AddEnum( ExampleModelApproximatorSupportVectorMachine())
  );

} // namespace bcl
