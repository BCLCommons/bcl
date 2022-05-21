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
#include "model/bcl_model_approximator_neural_network.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_template_instantiations.h"
#include "model/bcl_model_neural_network.h"
#include "model/bcl_model_transfer_sigmoid.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_approximator_neural_network.cpp
  //!
  //! @author mendenjl
  //! @date Sep 04, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelApproximatorNeuralNetwork :
    public ExampleInterface
  {
  public:

    ExampleModelApproximatorNeuralNetwork *Clone() const
    {
      return new ExampleModelApproximatorNeuralNetwork( *this);
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
      const size_t number_values( 10);
      storage::Vector< size_t> architecture;
      architecture.PushBack( 16);
      architecture.PushBack( 2);
      architecture.PushBack( 1);

      storage::Vector< storage::VectorND< 2, linal::Vector< float> > > features_results_setup;

      // test the functions Y1 = sin(x1) * cos(x2) + random(-0.1,0.1), Y2 = log(x1 + 1.0) + x2 + 1.0 with x1 = index, x2 = sqrt(index)
      for( size_t index( 0); index < number_values; index++)
      {
        features_results_setup.PushBack
        (
          storage::VectorND< 2, linal::Vector< float> >
          (
            linal::Vector< float>( architecture( 0)), linal::Vector< float>( architecture( 2))
          )
        );
        linal::Vector< float> &input = features_results_setup.LastElement().First();
        linal::Vector< float> &output = features_results_setup.LastElement().Second();
        input( 0) = sin( float( index));
        input( 1) = cos( float( index));
        for( size_t a = 2; a < architecture( 0); a++)
        {
          input( a) = sin( index + math::g_Pi * cos( float( a)));
        }

        output =
          linal::MakeVector
          (
            sin( input( 0)) * cos( input( 1)) + random::GetGlobalRandom().RandomGaussian( 0.0, 0.1)
          );
      }

      // set up the data set
      util::ShPtr< storage::Vector< storage::VectorND< 2, linal::Vector< float> > > >
        features_results( new storage::Vector< storage::VectorND< 2, linal::Vector< float> > >());
      for( size_t count( 0); count < number_values; ++count)
      {
        features_results->PushBack( features_results_setup( count));
      }

      util::Implementation< model::TransferFunctionInterface> transfer_func = model::TransferSigmoid();

      util::ShPtr< descriptor::Dataset> features_results_frds
      (
        new descriptor::Dataset( *features_results)
      );

      util::ShPtr< model::RescaleFeatureDataSet> rescale_in
      (
        new model::RescaleFeatureDataSet( *features_results_frds->GetFeaturesPtr(), model::NeuralNetwork::s_DefaultInputRange)
      );

      util::ShPtr< model::RescaleFeatureDataSet> rescale_out
      (
        new model::RescaleFeatureDataSet( *features_results_frds->GetResultsPtr(), model::TransferSigmoid().GetDynamicOutputRange())
      );

      // also try training an xor network with simple propagation
      static const size_t number_values_xor( 4);

      //initializing
      storage::VectorND< 2, linal::Vector< float> > features_results_setup_xor[ number_values_xor] =
      {
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 0.0), linal::MakeVector< float>( -1.0, 2.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 5.0), linal::MakeVector< float>( 2.0, -1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 1.0, 0.0), linal::MakeVector< float>( 2.0, -1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 1.0, 5.0), linal::MakeVector< float>( -1.0, 2.0))
      };

      // set up the data set
      util::ShPtr< storage::Vector< storage::VectorND< 2, linal::Vector< float> > > >
        features_results_xor( new storage::Vector< storage::VectorND< 2, linal::Vector< float> > >());
      for( size_t count( 0); count < 4; ++count)
      {
        features_results_xor->PushBack( features_results_setup_xor[ count]);
      }

      util::ShPtr< descriptor::Dataset> features_results_xor_frds
      (
        new descriptor::Dataset( *features_results_xor)
      );

      BCL_MessageStd( "FeatureResultDataSet xor: " + util::Format()( features_results_frds));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::ApproximatorNeuralNetwork irp_default;

      const size_t update_every_nth_feature( features_results_setup.GetSize());

      util::ShPtr< model::ObjectiveFunctionWrapper> sp_objective_function
      (
        new model::ObjectiveFunctionWrapper( features_results_frds, util::ObjectDataLabel( "RMSD"))
      );

      // construct resilient propagation iterate
      util::ShPtr< model::ApproximatorBase> sp_iterate_resilient
      (
        new model::ApproximatorNeuralNetwork
        (
          features_results_frds,
          update_every_nth_feature,
          architecture,
          transfer_func,
          sp_objective_function,
          util::Implementation< model::NeuralNetworkUpdateWeightsInterface>( "Resilient")
        )
      );

      // test rprop first
      sp_objective_function->SetData
      (
        features_results_frds,
        sp_iterate_resilient->GetRescaleFeatureDataSet()
      );

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      const float initial_rmsd( ( *sp_objective_function)( sp_iterate_resilient->GetCurrentModel()));

      // check the output of the untrained network
      BCL_MessageStd
      (
        "Initially, the RMSD was \n"
        + util::Format()( initial_rmsd)
      );

      // make sure the rmsd wasn't 0
      BCL_Example_Check
      (
        initial_rmsd > 0.0,
        " initial rmsd cannot be 0, otherwise there is nothing to train!"
      );

      // iterate (train) the network one step
      sp_iterate_resilient->Next();
      const util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> > improved_pair
      (
        sp_iterate_resilient->GetCurrentApproximation()
      );
      sp_iterate_resilient->Next();
      const util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> > improved_pair2
      (
        sp_iterate_resilient->GetCurrentApproximation()
      );

      // iterate (train) the network one step
      BCL_MessageStd
      (
        "After 1 step, RMSD is \n"
        + util::Format()( improved_pair->Second())
      );

      // iterate (train) the network one step
      BCL_MessageStd
      (
        "After 2 steps, RMSD is \n"
        + util::Format()( improved_pair2->Second())
      );

      util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> > final_pair;
      const size_t total_steps( 52);
      while( sp_iterate_resilient->GetTracker().GetIteration() < total_steps)
      {
        sp_iterate_resilient->Next();
      }
      final_pair = sp_iterate_resilient->GetCurrentApproximation();

      BCL_MessageStd
      (
        "After " + util::Format()( total_steps) + " steps, RMSD is \n"
        + util::Format()( final_pair->Second())
      );

      // check that the RMSD of the trained network was better than the initial (randomly-weighted) network
      BCL_Example_Check
      (
        final_pair->Second() < improved_pair->Second(),
        " neural network did not improve RMSD after " + util::Format()( total_steps) + " steps"
      );

      const util::ShPtr< model::ObjectiveFunctionWrapper> sp_objective_function_xor
      (
        new model::ObjectiveFunctionWrapper( features_results_xor_frds, util::ObjectDataLabel( "RMSD"))
      );

      // construct simple propagation iterate
      util::ShPtr< model::ApproximatorBase> sp_iterate_simple
      (
        new model::ApproximatorNeuralNetwork
        (
          features_results_xor_frds,
          0,
          storage::Vector< size_t>( 1, size_t( 3)),
          transfer_func,
          sp_objective_function_xor,
          util::Implementation< model::NeuralNetworkUpdateWeightsInterface>( "Simple(alpha = 0.5, eta = 0.25)"),
          1000
        )
      );

      // iterate (train) the network one step
      sp_iterate_simple->Next();
      const util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> > improved_pair_xor
      (
        sp_iterate_simple->GetCurrentApproximation()
      );

      // check the output of the untrained network
      BCL_MessageStd
      (
        "NN predicts \n"
        + util::Format()( improved_pair_xor->First()->operator()( *features_results_xor_frds->GetFeaturesPtr())( 0)) +
        " for \n" + util::Format()( features_results_xor_frds->GetResultsPtr()->operator()( 0)) + " before training!"
      );

      const size_t total_steps_xor( 25001);
      for( size_t step_number( 1); step_number < total_steps_xor; ++step_number)
      {
        sp_iterate_simple->Next();
      }
      final_pair = sp_iterate_simple->GetCurrentApproximation();

      // check the output of the trained network

      linal::Vector< float> pred_result( features_results_xor_frds->GetResultSize());
      pred_result = ( final_pair->First()->operator ()( *features_results_xor_frds->GetFeaturesPtr())( 0));
      linal::Vector< float> actual_result( features_results_xor_frds->GetResultSize());
      actual_result = ( features_results_xor_frds->GetResultsPtr()->operator()( 0));

      float comp_vector[ 2] = { -1.0, 2.0};
      linal::Vector< float> comp_vec( 2, comp_vector);

      BCL_MessageStd
      (
        "NN predicts \n"
        + util::Format()( pred_result)
        + " for \n" + util::Format()( actual_result)
      );

      BCL_Example_Check
      (
        ( pred_result - comp_vec).Norm() < 0.05,
        "The deviation between the results and ( -1, 2) should be smaller than 0.05 but is: "
        + util::Format()
        (
          ( pred_result - comp_vec).Norm()
        )
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelApproximatorNeuralNetwork

  const ExampleClass::EnumType ExampleModelApproximatorNeuralNetwork::s_Instance
  (
    GetExamples().AddEnum( ExampleModelApproximatorNeuralNetwork())
  );

} // namespace bcl
