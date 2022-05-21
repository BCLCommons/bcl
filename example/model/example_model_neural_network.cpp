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
#include "model/bcl_model_neural_network.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "math/bcl_math_template_instantiations.h"
#include "model/bcl_model_transfer_sigmoid.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_neural_network.cpp
  //!
  //! @author mueller, loweew
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelNeuralNetwork :
    public ExampleInterface
  {
  public:

    ExampleModelNeuralNetwork *Clone() const
    {
      return new ExampleModelNeuralNetwork( *this);
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
      static const size_t s_number_values( 4);

      //initializing
      storage::VectorND< 2, linal::Vector< float> > features_results_setup[ s_number_values] =
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
      storage::Vector< storage::VectorND< 2, linal::Vector< float> > > features_results;
      for( size_t count( 0); count < 4; ++count)
      {
        features_results.PushBack( features_results_setup[ count]);
      }

      util::Implementation< model::TransferFunctionInterface> transfer_func( new model::TransferSigmoid());

      // create feature result data set
      util::ShPtr< descriptor::Dataset> features_results_frds
      (
        new descriptor::Dataset( features_results)
      );

      // rescale function for in an output
      util::ShPtr< model::RescaleFeatureDataSet> rescale_in
      (
        new model::RescaleFeatureDataSet( *features_results_frds->GetFeaturesPtr(), model::NeuralNetwork::s_DefaultInputRange)
      );

      util::ShPtr< model::RescaleFeatureDataSet> rescale_out
      (
        new model::RescaleFeatureDataSet( *features_results_frds->GetResultsPtr(), model::TransferSigmoid().GetDynamicOutputRange())
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::NeuralNetwork nn_default;

      storage::Vector< linal::Vector< float> > bias;
      storage::Vector< linal::Matrix< float> > weights;
      bias.PushBack( linal::Vector< float>( 3));
      bias.PushBack( linal::Vector< float>( 2));
      weights.PushBack( linal::Matrix< float>( 3, 2));
      weights.PushBack( linal::Matrix< float>( 2, 3));

      // constructor from architecture and transfer function
      model::NeuralNetwork
        nn
        (
          rescale_in,
          rescale_out,
          bias,
          weights,
          transfer_func
        );

      // clone
      util::ShPtr< model::NeuralNetwork> nn_clone( nn.Clone());

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleAssert( nn_clone->GetClassIdentifier(), GetStaticClassName< model::NeuralNetwork>());

      BCL_Example_Check
      (
        nn_clone->GetNumberLayers() == 2
        && nn_clone->GetWeight()( 0).GetNumberRows() == 3
        && nn_clone->GetWeight()( 1).GetNumberRows() == 2,
        "Neural network Copy should contain two layers (one with three hidden neurons, one with two outputs),"
        " but contains " + util::Format()( nn_clone->GetNumberLayers()) + " hidden layer(s) and "
        + util::Format()( nn_clone->GetWeight()( 0).GetNumberRows()) + " and "
        + util::Format()( nn_clone->GetWeight()( 1).GetNumberRows()) + " hidden neurons."
      );

    ///////////////
    // operators //
    ///////////////

      linal::Vector< float> pred_results( features_results_frds->GetResultSize());
      pred_results = ( nn_clone->operator()( *features_results_frds->GetFeaturesPtr())( 0));
      linal::Vector< float> actual_results( features_results_frds->GetResultSize());
      actual_results = ( features_results_frds->GetResults().DeScale()( 0));

      // check the output of the untrained network
      BCL_MessageStd
      (
        "NN predicts \n"
        + util::Format()( pred_results)
        + " for \n" + util::Format()( actual_results)
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write and read
      BCL_MessageStd( "testing read and write functionalities");
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( nn);
      BCL_MessageVrb( "read object");
      model::NeuralNetwork nn_read;
      ReadBCLObject( nn_read);

      linal::Vector< float> nn_result( features_results_frds->GetResultSize());
      nn_result = nn( *features_results_frds->GetFeaturesPtr())( 0);
      linal::Vector< float> nn_read_result( features_results_frds->GetResultSize());
      nn_read_result = nn_read( *features_results_frds->GetFeaturesPtr())( 0);

      // compare predictions
      BCL_ExampleIndirectCheckWithinTolerance( nn_result.Norm(), nn_read_result.Norm(), 0.001, "I/O");

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelNeuralNetwork

  const ExampleClass::EnumType ExampleModelNeuralNetwork::s_Instance
  (
    GetExamples().AddEnum( ExampleModelNeuralNetwork())
  );

} // namespace bcl
