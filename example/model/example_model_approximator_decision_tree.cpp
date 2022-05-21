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
#include "model/bcl_model_approximator_decision_tree.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_template_instantiations.h"
#include "model/bcl_model_dtree_roc_data_partition_function.h"
#include "model/bcl_model_objective_function_accuracy.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_approximator_decision_tree.cpp
  //!
  //! @author mendenjl
  //! @date Sep 03, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelApproximatorDecisionTree :
    public ExampleInterface
  {
  public:

    ExampleModelApproximatorDecisionTree *Clone() const
    {
      return new ExampleModelApproximatorDecisionTree( *this);
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
      storage::Vector< storage::Pair< storage::VectorND< 2, linal::Vector< float> >, bool> >::s_Instance.IsDefined();

      //make parameters for next test
      const float activity_cutoff( 8);

      //make parameters for next test
      util::ShPtr< storage::Vector< storage::VectorND< 2, linal::Vector< float> > > > sp_training_data( new storage::Vector< storage::VectorND< 2, linal::Vector< float> > >());
      util::ShPtr< storage::Vector< storage::VectorND< 2, linal::Vector< float> > > > sp_monitoring_data( new storage::Vector< storage::VectorND< 2, linal::Vector< float> > >());

      util::ShPtr< descriptor::Dataset> sp_monitor_frds
      (
        new descriptor::Dataset( *sp_monitoring_data)
      );
      util::ShPtr< descriptor::Dataset> sp_training_frds
      (
        new descriptor::Dataset( *sp_training_data)
      );
      util::ShPtr< model::ObjectiveFunctionWrapper> sp_objective_function
      (
        new model::ObjectiveFunctionWrapper
        (
          sp_monitor_frds,
          util::Implementation< model::ObjectiveFunctionInterface>( model::ObjectiveFunctionAccuracy( activity_cutoff))
        )
      );

      util::Implementation< model::DtreeDataPartitionFunctionInterface> sp_dataset_partitioner
      (
        new model::DtreeRocDataPartitionFunction()
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::ApproximatorDecisionTree iterate_default;

      // constructor with parameters
      model::ApproximatorDecisionTree iterate_with_parameters
      (
        sp_dataset_partitioner,
        sp_objective_function,
        activity_cutoff,
        sp_training_frds
      );

      // check clone
      util::ShPtr< model::ApproximatorDecisionTree> iterate_clone( iterate_with_parameters.Clone());

    /////////////////
    // data access //
    /////////////////

      //default constructor
      util::ShPtr< model::ApproximatorDecisionTree> sp_roc_iterate( new model::ApproximatorDecisionTree());

      //constructor taking all necessary parameters
      util::ShPtr< model::ApproximatorDecisionTree> sp_roc_iterate_two
      (
        new model::ApproximatorDecisionTree( sp_dataset_partitioner, sp_objective_function, activity_cutoff, sp_training_frds)
      );

      util::ShPtr< model::ObjectiveFunctionWrapper> new_objective
      (
        new model::ObjectiveFunctionWrapper( sp_monitor_frds, util::ObjectDataLabel( "RMSD"))
      );

      sp_roc_iterate_two->SetObjectiveFunction( new_objective);

      BCL_ExampleAssert( sp_roc_iterate_two->GetObjectiveFunction().IsDefined(), true);

      // instantiate a new iterate
      util::ShPtr< model::ApproximatorDecisionTree> sp_iterate( BuildROCIterateToyDataOne());

      //GetActivityCutoff
      BCL_ExampleAssert( sp_iterate->GetActivityCutoff(), float( 13));

      //GetTrainingData
      BCL_ExampleAssert( sp_iterate->GetTrainingData()->GetSize(), size_t( 20));

      //SetTrainingData
      sp_iterate->SetTrainingData( sp_training_frds);
      BCL_ExampleAssert( sp_iterate->GetTrainingData()->GetSize(), size_t( 0));

      //GetObjectiveFunction
      BCL_ExampleAssert( sp_iterate->GetObjectiveFunction().IsDefined(), true);

      const util::ShPtr< model::ObjectiveFunctionWrapper> old_objective_function = sp_iterate->GetObjectiveFunction();
      sp_iterate->SetObjectiveFunction( new_objective);
      BCL_ExampleAssert
      (
        sp_iterate->GetObjectiveFunction() != old_objective_function,
        true
      );

      //GetDataPartitionFunction
      BCL_ExampleAssert( sp_iterate->GetDataPartitionFunction().IsDefined(), true);

    ///////////////
    // operators //
    ///////////////

      // instantiate a new iterate for toy data
      sp_iterate = BuildROCIterateToyDataOne();

//      util::ShPtr< model::Interface> tree( new model::DecisionTree());
//
//      util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> > prev_result
//      (
//        new storage::Pair< util::ShPtr< model::Interface>, float>
//        (
//          tree, sp_iterate->GetObjectiveFunction()->operator ()( tree)
//        )
//      );
//
//      sp_iterate->Initialize( prev_result);

      //operator()
//      util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> > result_pair = sp_iterate->operator ()();
//      BCL_Example_Check
//      (
//        ExampleClass::ExampleResult::e_Trivial,
//        result_pair->Second() >= prev_result->Second(),
//        "accuracy should have improved from "
//        + util::Format()( prev_result->Second()) + " to " + util::Format()( result_pair->Second())
//      )

    //////////////////////
    // input and output //
    //////////////////////

      //Initialize for writing
      sp_iterate = BuildROCIterateToyDataOne();

//      tree = util::ShPtr< model::Interface>( new model::DecisionTree());

//      util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> > init_result
//      (
//        new storage::Pair< util::ShPtr< model::Interface>, float>
//        (
//          tree, sp_iterate->GetObjectiveFunction()->operator ()( tree)
//        )
//      );

//      sp_iterate->Initialize( init_result);

      // write bcl object
//      WriteBCLObject( *sp_iterate);

      // create default object
//      model::DecisionTreeIterate iterate_read;

      // read bcl object
//      ReadBCLObject( iterate_read);

      // check get oversampling factor method
//      BCL_ExampleAssert( iterate_read.GetActivityCutoff(), float( 13));

    //////////////////////
    // Additional Tests //
    //////////////////////

      //Toy data set should get 100%
      sp_iterate = BuildROCIterateToyDataTwo();

//      util::ShPtr< model::Interface> sp_tree( new model::DecisionTree());

      // do the approximation
//      storage::Pair< util::ShPtr< model::Interface>, float> result( approximator->Approximate( sp_tree));

//      BCL_ExampleAssert( result.Second(), 1.0);
//
//      //Toy data set should get 100%
//      sp_iterate = BuildROCIterateToyDataOne();
//      approximator = BuildApproximator( sp_iterate);
//      sp_tree = util::ShPtr< model::Interface>( new model::DecisionTree());
//
//      // do the approximation
//      result = approximator->Approximate( sp_tree);
//
//      BCL_ExampleAssert( result.Second(), 1.0);

      return 0;
    } //end ExampleClass::ExampleDecisionTreeterate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief constructs ROCIterate object from toy data. creates a dataset with 20 items. Each item has 5 descriptors.
    //!        the first four descriptors are random, the last one is the counter value. the target output is
    //!        also the counter value.
    //! @return ShPtr to new ROCIterate object
    util::ShPtr< model::ApproximatorDecisionTree> BuildROCIterateToyDataOne() const
    {
      // cutoff for the toy example
      const float activity_cutoff( 13);

      //! create a training set
      util::ShPtr< storage::Vector< storage::VectorND< 2, linal::Vector< float> > > > sp_training_data_set
      (
        new storage::Vector< storage::VectorND< 2, linal::Vector< float> > >()
      );

      for( float counter( 0); counter < 20; ++counter)
      {
        linal::Vector< float> descriptors( 5);
        descriptors( 0) = random::GetGlobalRandom().Random( float( 100));
        descriptors( 1) = random::GetGlobalRandom().Random( float( 100));
        descriptors( 2) = random::GetGlobalRandom().Random( float( 100));
        descriptors( 3) = random::GetGlobalRandom().Random( float( 100));
        descriptors( 4) = counter;

        linal::Vector< float> target( linal::MakeVector< float>( counter));
        sp_training_data_set->PushBack
        (
          storage::VectorND< 2, linal::Vector< float> >( descriptors, target)
        );
      }

      //! create a monitoring set
      util::ShPtr< storage::Vector< storage::VectorND< 2, linal::Vector< float> > > > sp_monitoring_data_set( new storage::Vector< storage::VectorND< 2, linal::Vector< float> > >());
      for( float counter( 0); counter < 20; counter += 10)
      {
        linal::Vector< float> descriptors( 5);
        descriptors( 0) = random::GetGlobalRandom().Random( float( 100));
        descriptors( 1) = random::GetGlobalRandom().Random( float( 100));
        descriptors( 2) = random::GetGlobalRandom().Random( float( 100));
        descriptors( 3) = random::GetGlobalRandom().Random( float( 100));
        descriptors( 4) = counter;

        linal::Vector< float> target( linal::MakeVector< float>( counter));
        sp_monitoring_data_set->PushBack
        (
          storage::VectorND< 2, linal::Vector< float> >( descriptors, target)
        );
      }

      // setup parts of ROCIterate
      // objective function used by the splitter and approximator
      util::ShPtr< descriptor::Dataset> sp_monitoring_dataset
      (
        new descriptor::Dataset( *sp_monitoring_data_set)
      );
      util::ShPtr< model::ObjectiveFunctionWrapper> objective_function
      (
        new model::ObjectiveFunctionWrapper( sp_monitoring_dataset, util::ObjectDataLabel( "Accuracy( cutoff = 13)"))
      );

      util::ShPtr< descriptor::Dataset> sp_training_dataset
      (
        new descriptor::Dataset( *sp_training_data_set)
      );
      util::Implementation< model::DtreeDataPartitionFunctionInterface> data_partitioner
      (
        new model::DtreeRocDataPartitionFunction( 0.5, 1.0)
      );
      util::ShPtr< model::ApproximatorDecisionTree> sp_roc_iterate
      (
        new model::ApproximatorDecisionTree( data_partitioner, objective_function, activity_cutoff, sp_training_dataset)
      );

      return sp_roc_iterate;
    } // BuildROCIterateToyDataOne

    //! @brief constructs ROCIterate object from toy data. This dataset has the target value at index 1, and the other
    //!        descriptors are put there to distract. the target is the sin of the counter.
    //! @return ShPtr to new ROCIterate object
    util::ShPtr< model::ApproximatorDecisionTree> BuildROCIterateToyDataTwo() const
    {
      // cutoff for the toy example
      const float activity_cutoff( 0);

      // create a training dataset
      util::ShPtr< storage::Vector< storage::VectorND< 2, linal::Vector< float> > > > sp_training_data_set( new storage::Vector< storage::VectorND< 2, linal::Vector< float> > >());

      for( int counter( 0); counter < 100; counter++)
      {
        linal::Vector< float> v( 4);
        v( 0) = std::cos( float( counter));
        v( 1) = std::sin( float( counter));
        v( 2) = random::GetGlobalRandom().Random( float( 100));
        v( 3) = float( counter);

        sp_training_data_set->PushBack( storage::VectorND< 2, linal::Vector< float> >( v, linal::MakeVector< float>( std::sin( float( counter)))));
      }

      //! create a monitoring dataset for the evaluator function
      util::ShPtr< storage::Vector< storage::VectorND< 2, linal::Vector< float> > > > sp_monitoring_data_set( new storage::Vector< storage::VectorND< 2, linal::Vector< float> > >());

      // create toy data with different sampling
      for( int counter( 0); counter < 100; counter += 10)
      {
        linal::Vector< float> vector( 4);
        vector( 0) = std::cos( float( counter));
        vector( 1) = std::sin( float( counter));
        vector( 2) = random::GetGlobalRandom().Random( float( 100));
        vector( 3) = float( counter);

        linal::Vector< float> target = linal::MakeVector< float>( std::sin( float( counter)));
        sp_monitoring_data_set->PushBack( storage::VectorND< 2, linal::Vector< float> >( vector, target));
      }

      // objective function used by the splitter and approximator
      util::ShPtr< descriptor::Dataset> sp_monitoring_dataset
      (
        new descriptor::Dataset( *sp_monitoring_data_set)
      );
      util::ShPtr< model::ObjectiveFunctionWrapper> objective_function
      (
        new model::ObjectiveFunctionWrapper( sp_monitoring_dataset, util::ObjectDataLabel( "Accuracy( cutoff = 13)"))
      );
      util::ShPtr< descriptor::Dataset> sp_training_dataset
      (
        new descriptor::Dataset( *sp_training_data_set)
      );
      util::Implementation< model::DtreeDataPartitionFunctionInterface> sp_node_splitter
      (
        new model::DtreeRocDataPartitionFunction( 0.5, 1.0)
      );
      util::ShPtr< model::ApproximatorDecisionTree> sp_roc_iterate
      (
        new model::ApproximatorDecisionTree( sp_node_splitter, objective_function, activity_cutoff, sp_training_dataset)
      );
      return sp_roc_iterate;
    } // BuildROCIterateToyDataTwo

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelApproximatorDecisionTree

  const ExampleClass::EnumType ExampleModelApproximatorDecisionTree::s_Instance
  (
    GetExamples().AddEnum( ExampleModelApproximatorDecisionTree())
  );

} // namespace bcl
