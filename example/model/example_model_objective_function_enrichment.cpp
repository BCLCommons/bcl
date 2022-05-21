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
#include "model/bcl_model_objective_function_enrichment.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_template_instantiations.h"
#include "model/bcl_model_neural_network.h"
#include "model/bcl_model_objective_function_wrapper.h"
#include "model/bcl_model_transfer_sigmoid.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_objective_function_enrichment_monitoring_dataset.cpp
  //!
  //! @author butkiem1, loweew
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelObjectiveFunctionEnrichment :
    public ExampleInterface
  {
  public:

    ExampleModelObjectiveFunctionEnrichment *Clone() const
    {
      return new ExampleModelObjectiveFunctionEnrichment( *this);
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
      // monitoring dataset
      util::ShPtr< storage::Vector< storage::VectorND< 2, linal::Vector< float> > > >
      monitoring_data
      (
        new storage::Vector< storage::VectorND< 2, linal::Vector< float> > >()
      );

      // generate monitoring data for dataset
      for( size_t counter( 0); counter < 100; ++counter)
      {
        monitoring_data->PushBack
        (
          storage::VectorND< 2, linal::Vector< float> >
          (
            linal::MakeVector< float>( counter, counter + 1, counter + 2),
            linal::MakeVector< float>( counter % 2)
          )
        );
      }

      // feature result data set construction
      util::ShPtr< descriptor::Dataset> monitor_frds
      (
        new descriptor::Dataset( *monitoring_data)
      );

      // setting up architecture of neural network
      storage::Vector< size_t> architecture;
      architecture.PushBack( 3);
      architecture.PushBack( 2);
      architecture.PushBack( 1);

      storage::Vector< linal::Vector< float> > bias;
      storage::Vector< linal::Matrix< float> > weights;
      bias.PushBack( linal::Vector< float>( architecture( 1)));
      bias.PushBack( linal::Vector< float>( architecture( 2)));
      weights.PushBack( linal::Matrix< float>( architecture( 1), architecture( 0)));
      weights.PushBack( linal::Matrix< float>( architecture( 2), architecture( 1)));

      // rescale function for in an output and denormalization
      util::ShPtr< model::RescaleFeatureDataSet> rescale_in
      (
        new model::RescaleFeatureDataSet( *monitor_frds->GetFeaturesPtr(), model::NeuralNetwork::s_DefaultInputRange)
      );

      util::ShPtr< model::RescaleFeatureDataSet> rescale_out
      (
        new model::RescaleFeatureDataSet( *monitor_frds->GetResultsPtr(), model::TransferSigmoid().GetDynamicOutputRange())
      );

      // constructor from architecture and transfer function
      util::ShPtr< model::Interface> neural_network_model
      (
        new model::NeuralNetwork
        (
          rescale_in, //rescale_in,
          rescale_out, //rescale_out,
          bias,
          weights,
          util::Implementation< model::TransferFunctionInterface>( model::TransferSigmoid())
        )
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::ObjectiveFunctionWrapper objective_function_default( monitor_frds, util::ObjectDataLabel( "Enrichment"));

      // construct from properties
      model::ObjectiveFunctionWrapper objective_function
      (
        monitor_frds, util::ObjectDataLabel( "Enrichment(cutoff=0.5, FPR cutoff=0.05, parity=1)")
      );

      // construct from properties
      model::ObjectiveFunctionWrapper objective_function_40pct
      (
        monitor_frds, util::ObjectDataLabel( "Enrichment(cutoff=0.5, FPR cutoff=0.4, parity=1)")
      );

      // check clone
      const util::ShPtr< model::ObjectiveFunctionWrapper> objective_function_clone
      (
        objective_function.Clone()
      );

    /////////////////
    // data access //
    /////////////////

//      // set cutoff
//      specialization.SetCutoff( float( 0.5));
      // check GetCutoff method
      BCL_ExampleCheck( objective_function.GetThreshold(), float( 0.5));

//      // set oversampling factor
//      specialization.SetOverSamplingFactor( float( 1));
//      // check get oversampling factor method
//      BCL_ExampleAssert( specialization.GetOverSamplingFactor(), float( 1));
//      // check size of set monitoring data
//      BCL_ExampleAssert( objective_function.GetData()->GetSize(), size_t( 100));

    ///////////////
    // operators //
    ///////////////

      // check get and set methods
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          objective_function( neural_network_model),
          float( 1.0),
          0.01
        ),
        " objective_function: wrong obj function result! should be 1.0 but is "
        + util::Format()( objective_function( neural_network_model))
      );

      // operator.  Note that since the dataset had only two values, we expect the neural network to do well
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          objective_function_40pct( neural_network_model),
          float( 1.0),
          0.01
        ),
        " objective_function_40pct: wrong obj function result! should be 1.0 but is "
        + util::Format()( objective_function_40pct( neural_network_model))
      );

    //////////////////////
    // input and output //
    //////////////////////

//      // write bcl object
//      WriteBCLObject( specialization);
//      // create default object
//      model::ObjectiveFunctionEnrichment objective_function_read;
//      // read bcl object
//      ReadBCLObject( objective_function_read);
//      // check get oversampling factor method
//      BCL_ExampleAssert( objective_function_read.GetOverSamplingFactor(), float( 1));
//      // check GetCutoff method
//      BCL_ExampleAssert( objective_function_read.GetCutoff(), float( 0.5));

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelObjectiveFunctionEnrichment

  const ExampleClass::EnumType ExampleModelObjectiveFunctionEnrichment::s_Instance
  (
    GetExamples().AddEnum( ExampleModelObjectiveFunctionEnrichment())
  );

} // namespace bcl
