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
#include "model/bcl_model_support_vector_machine_multi_output.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "math/bcl_math_template_instantiations.h"
#include "model/bcl_model_support_vector_kernel_rbf.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_support_vector_machine_multi_output.cpp
  //!
  //! @author butkiem1, mendenjl
  //! @date Apr 18, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelSupportVectorMachineMultiOutput :
    public ExampleInterface
  {
  public:

    ExampleModelSupportVectorMachineMultiOutput *Clone() const
    {
      return new ExampleModelSupportVectorMachineMultiOutput( *this);
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
      static const size_t s_number_values( 5);

      linal::Matrix< float> features( s_number_values, size_t( 2));
      features.ReplaceRow( 0, linal::MakeVector< float>( 0.0, 0.0));
      features.ReplaceRow( 1, linal::MakeVector< float>( 0.0, 5.0));
      features.ReplaceRow( 2, linal::MakeVector< float>( 1.0, 0.0));
      features.ReplaceRow( 3, linal::MakeVector< float>( 1.0, 5.0));
      features.ReplaceRow( 4, linal::MakeVector< float>( 5.0, 1.0));

      linal::Matrix< float> results( s_number_values, size_t( 1));
      results.ReplaceRow( 0, linal::MakeVector< float>( -1.0));
      results.ReplaceRow( 1, linal::MakeVector< float>(  2.0));
      results.ReplaceRow( 2, linal::MakeVector< float>(  2.0));
      results.ReplaceRow( 3, linal::MakeVector< float>( -1.0));
      results.ReplaceRow( 4, linal::MakeVector< float>(  2.0));

      // LaGrange Multipliers alpha
      storage::Vector< float> alphas( s_number_values, float( 1.56));

      // threshold of hyperplane
      const float bias( 0.5);

      // kernel
      util::Implementation< model::SupportVectorKernelBase> kernel = model::SupportVectorKernelRBF();

      // number of bound support vectors
      const size_t number_bound_sv( 2);

      //  number of support vectors
      const size_t number_sv( 5);

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
        new model::RescaleFeatureDataSet( *features_results_frds->GetFeaturesPtr(), model::SupportVectorMachineMultiOutput::s_DefaultInputRange)
      );

      util::ShPtr< model::RescaleFeatureDataSet> rescale_result
      (
        new model::RescaleFeatureDataSet( *features_results_frds->GetResultsPtr(), model::SupportVectorMachineMultiOutput::s_DefaultInputRange)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      model::SupportVectorMachineMultiOutput svm_default;

      model::SupportVectorMachineMultiOutput svm_init( kernel, rescale_feature, rescale_result);

      svm_init.SetAlpha( alphas);
      svm_init.SetBias( bias);
      svm_init.SetKernel( kernel);
      svm_init.SetNumberBoundSupportVectors( number_bound_sv);
      svm_init.SetNumberSupportVectors( number_sv);
      svm_init.SetSupportVectors( rescale_feature->operator ()( *( features_results_frds->GetFeaturesPtr())));
      model::SupportVectorMachineMultiOutput svm_copy( svm_init);

      util::ShPtr< model::SupportVectorMachineMultiOutput> svm_clone( svm_init.Clone());

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      model::FeatureDataSet< float> result( svm_init( *test_features_results_frds->GetFeaturesPtr()));

      BCL_MessageStd( "Prediction: " + util::Format()( result));

      // check operator
      BCL_Example_Check
      (
        math::EqualWithinTolerance
        (
          result( 2)( 0),
          float( 15.4233),
          0.01
        ),
        "Prediction is incorrect! should be 15.4233 but is "
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

  }; //end ExampleModelSupportVectorMachineMultiOutput

  const ExampleClass::EnumType ExampleModelSupportVectorMachineMultiOutput::s_Instance
  (
    GetExamples().AddEnum( ExampleModelSupportVectorMachineMultiOutput())
  );

} // namespace bcl
