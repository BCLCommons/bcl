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
#include "opencl/bcl_opencl_support_vector_machine.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "model/bcl_model_support_vector_kernel_rbf.h"
#include "model/bcl_model_support_vector_machine.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_support_vector_machine.cpp
  //!
  //! @author loweew
  //! @date Mar 14, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclSupportVectorMachine :
    public ExampleInterface
  {
  public:

    ExampleOpenclSupportVectorMachine *Clone() const
    {
      return new ExampleOpenclSupportVectorMachine( *this);
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

      linal::Matrix< float> features( size_t( 5), size_t( 2));
      features.ReplaceRow( 0, linal::MakeVector< float>( 0.0, 0.0));
      features.ReplaceRow( 1, linal::MakeVector< float>( 0.0, 5.0));
      features.ReplaceRow( 2, linal::MakeVector< float>( 1.0, 0.0));
      features.ReplaceRow( 3, linal::MakeVector< float>( 1.0, 5.0));
      features.ReplaceRow( 4, linal::MakeVector< float>( 5.0, 1.0));

      linal::Matrix< float> results( size_t( 5), size_t( 1));
      results.ReplaceRow( 0, linal::MakeVector< float>( -1.0));
      results.ReplaceRow( 1, linal::MakeVector< float>( 2.0));
      results.ReplaceRow( 2, linal::MakeVector< float>( 2.0));
      results.ReplaceRow( 3, linal::MakeVector< float>( -1.0));
      results.ReplaceRow( 4, linal::MakeVector< float>( 2.0));

      // LaGrange Multipliers alpha
      linal::Vector< float> alphas( 5, float( 1.56));

      // threshold of hyperplane
      const float bias( 0.5);

      // kernel
      util::Implementation< model::SupportVectorKernelBase> kernel = model::SupportVectorKernelRBF();

      // create feature result data set
      util::ShPtr< descriptor::Dataset> features_results_frds
      (
        new descriptor::Dataset( features, results)
      );

      util::ShPtr< descriptor::Dataset> test_features_results_frds
      (
        new descriptor::Dataset( features, results)
      );

      // rescale function for in an output
      util::ShPtr< model::RescaleFeatureDataSet> rescale_feature
      (
        new model::RescaleFeatureDataSet( *features_results_frds->GetFeaturesPtr(), model::SupportVectorMachine::s_DefaultInputRange)
      );

      util::ShPtr< model::RescaleFeatureDataSet> rescale_result
      (
        new model::RescaleFeatureDataSet( *features_results_frds->GetResultsPtr(), model::SupportVectorMachine::s_DefaultInputRange)
      );

      opencl::CommandQueue queue( opencl::GetTools().GetFirstCommandQueue());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      opencl::SupportVectorMachine svm_init
      (
        bias,
        alphas,
        rescale_feature->operator()( *( features_results_frds->GetFeaturesPtr())),
        kernel,
        *rescale_feature,
        *rescale_result,
        queue
      );
      util::ShPtr< opencl::SupportVectorMachine> svm_clone( svm_init.Clone());

    /////////////////
    // data access //
    /////////////////

      // check class identifier
      BCL_Example_Check
      (
        GetStaticClassName< opencl::SupportVectorMachine>() == svm_init.GetClassIdentifier(),
        "incorrect class identifier"
      );

    ///////////////
    // operators //
    ///////////////

      util::Stopwatch timer;

      model::FeatureDataSet< float> result( svm_init( ( *test_features_results_frds->GetFeaturesPtr())));

      BCL_MessageStd( "Time to predict using gpu: " + util::Format()( timer.GetProcessDuration().GetTimeAsHourMinuteSecondMilliSeconds()));

      BCL_MessageStd( "predicted results: " + util::Format()( result));

      // check operator
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          result( 2)( 0),
          float( 12.4233),
          0.01
        ),
        "Prediction is incorrect! should be 12.4233 but is "
        + util::Format()( result( 2)( 0))
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclSupportVectorMachine

  const ExampleClass::EnumType ExampleOpenclSupportVectorMachine::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclSupportVectorMachine())
  );

} // namespace bcl
