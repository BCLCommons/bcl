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
#include "opencl/bcl_opencl_neural_network.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "math/bcl_math_template_instantiations.h"
#include "model/bcl_model_neural_network.h"
#include "model/bcl_model_transfer_sigmoid.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_neural_network.cpp
  //!
  //! @author loweew
  //! @date Mar 28, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclNeuralNetwork :
    public ExampleInterface
  {
  public:

    ExampleOpenclNeuralNetwork *Clone() const
    {
      return new ExampleOpenclNeuralNetwork( *this);
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

      linal::Matrix< float> features( size_t( 4), size_t( 2));
      features.ReplaceRow( 0, linal::MakeVector< float>( 0.0, 0.0));
      features.ReplaceRow( 1, linal::MakeVector< float>( 0.0, 5.0));
      features.ReplaceRow( 2, linal::MakeVector< float>( 1.0, 0.0));
      features.ReplaceRow( 3, linal::MakeVector< float>( 1.0, 5.0));

      linal::Matrix< float> results( size_t( 4), size_t( 2));
      results.ReplaceRow( 0, linal::MakeVector< float>( -1.0, 2.0));
      results.ReplaceRow( 1, linal::MakeVector< float>( 2.0, -1.0));
      results.ReplaceRow( 2, linal::MakeVector< float>( 2.0, -1.0));
      results.ReplaceRow( 3, linal::MakeVector< float>( -1.0, 2.0));

      util::Implementation< model::TransferFunctionInterface> transfer_func( new model::TransferSigmoid());

      // create feature result data set
      util::ShPtr< descriptor::Dataset> features_results_frds
      (
        new descriptor::Dataset( features, results)
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

      const opencl::CommandQueue queue( opencl::GetTools().GetFirstCommandQueue());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      opencl::NeuralNetwork nn_default( queue);

      storage::Vector< linal::Vector< float> > bias;
      storage::Vector< linal::Matrix< float> > weights;
      bias.PushBack( linal::Vector< float>( 3));
      bias.PushBack( linal::Vector< float>( 2));
      weights.PushBack( linal::Matrix< float>( 3, 2));
      weights.PushBack( linal::Matrix< float>( 2, 3));

      storage::Vector< opencl::Vector< float> > bias_buffers;
      storage::Vector< opencl::Matrix< float> > weight_buffers;

      const cl_uint block_size( 16);

      for( size_t k( 0); k < bias.GetSize(); ++k)
      {
        const cl_uint bias_pad( ( block_size - ( bias( k).GetSize() % block_size)) % block_size);
        const cl_uint weight_row_pad( ( block_size - ( weights( k).GetNumberRows() % block_size)) % block_size);
        const cl_uint weight_col_pad( ( block_size - ( weights( k).GetNumberCols() % block_size)) % block_size);
        bias_buffers.PushBack( opencl::Vector< float>( bias( k), queue, bias_pad));
        weight_buffers.PushBack( opencl::Matrix< float>( weights( k), queue, weight_row_pad, weight_col_pad));
      }

      // constructor from architecture and transfer function
      opencl::NeuralNetwork
        nn
        (
          rescale_in,
          rescale_out,
          bias,
          weights,
          transfer_func,
          queue
        );

      // gpu constructor
      opencl::NeuralNetwork gpu_nn( rescale_in, rescale_out, bias_buffers, weight_buffers, queue);

      // clone
      util::ShPtr< opencl::NeuralNetwork> nn_clone( nn.Clone());

    /////////////////
    // data access //
    /////////////////

      // class identifiers
      BCL_MessageStd( "class name: " + nn_clone->GetClassIdentifier());
      BCL_Example_Check
      (
        (
          GetStaticClassName< opencl::NeuralNetwork>( queue) == nn.GetClassIdentifier()
        ),
        "incorrect class identifier"
      );

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

      model::FeatureDataSet< float> output_results( gpu_nn( *features_results_frds->GetFeaturesPtr()));

      // check the output of the untrained network
      BCL_MessageStd
      (
        "NN predicts \n"
        + util::Format()( pred_results)
        + " for \n" + util::Format()( actual_results)
      );

      BCL_MessageStd
      (
        "NN gpu predicts \n"
        + util::Format()( output_results( 0))
        + " for \n" + util::Format()( actual_results)
      );

      // checking that cpu and gpu predict the same
      BCL_ExampleCheckWithinTolerance( double( output_results( 0)( 0)), double( pred_results( 0)), 0.001);

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write and read
      // open output file
      BCL_MessageVrb( "write object");
      WriteBCLObject( nn);
      BCL_MessageVrb( "read object");
      opencl::NeuralNetwork nn_read( queue);
      ReadBCLObject( nn_read);

      linal::Vector< float> nn_result( nn( *features_results_frds->GetFeaturesPtr())( 0));

      linal::Vector< float> nn_read_result( nn_read( *features_results_frds->GetFeaturesPtr())( 0));

      // compare predictions
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          nn_result.Norm(),
          nn_read_result.Norm()
        ),
        "written neural network predicts differently"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclNeuralNetwork

  const ExampleClass::EnumType ExampleOpenclNeuralNetwork::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclNeuralNetwork())
  );

} // namespace bcl
