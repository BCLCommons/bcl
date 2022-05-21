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
#include "model/bcl_model_feature_result_and_state.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_feature_result_and_state.cpp
  //!
  //! @author mendenjl
  //! @date May 09, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelFeatureResultAndState :
    public ExampleInterface
  {
  public:

    ExampleModelFeatureResultAndState *Clone() const
    {
      return new ExampleModelFeatureResultAndState( *this);
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
      // number elements per dimension
      const size_t number_elements( 2), results_size( 1);

      // constructing data for feature reference construction
      const linal::Matrix< float> data( number_elements, number_elements, 1.5);

      // constructing data for results construction
      const linal::Matrix< float> data_results( number_elements, results_size, 3.6);

      // constructing data for feature reference construction
      const linal::Matrix< size_t> data_results_states( number_elements, results_size, 2);

      // constructing feature references to use in construction of feature result states
      const model::FeatureReference< float> feat_ref_a( number_elements, data[ 0]);
      const model::FeatureReference< float> feat_ref_b( number_elements, data[ 1]);
      const model::FeatureReference< float> result_ref_a( results_size, data_results[ 0]);
      const model::FeatureReference< float> result_ref_b( results_size, data_results[ 1]);
      const model::FeatureReference< size_t> result_state_ref_a( results_size, data_results_states[ 0]);
      const model::FeatureReference< size_t> result_state_ref_b( results_size, data_results_states[ 1]);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // constructor
      model::FeatureResultAndState feature_result_state_a( feat_ref_a, result_ref_a, result_state_ref_a);
      model::FeatureResultAndState feature_result_state_b( feat_ref_b, result_ref_b, result_state_ref_b);

      // clone
      const util::ShPtr< model::FeatureResultAndState> feature_result_state_a_clone( feature_result_state_a.Clone());

    /////////////////
    // data access //
    /////////////////

      // get feature
      BCL_ExampleCheck( feature_result_state_a.GetFeature().Begin(), feat_ref_a.Begin());
      BCL_ExampleCheck( feature_result_state_b.GetFeature().Begin(), feat_ref_b.Begin());

      // get result
      BCL_ExampleCheck( feature_result_state_a.GetResult().Begin(), result_ref_a.Begin());
      BCL_ExampleCheck( feature_result_state_b.GetResult().Begin(), result_ref_b.Begin());

      // get result_state
      BCL_ExampleCheck( feature_result_state_a.GetResultState().Begin(), result_state_ref_a.Begin());
      BCL_ExampleCheck( feature_result_state_b.GetResultState().Begin(), result_state_ref_b.Begin());

      // end
      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelFeatureResultAndState

  const ExampleClass::EnumType ExampleModelFeatureResultAndState::s_Instance
  (
    GetExamples().AddEnum( ExampleModelFeatureResultAndState())
  );

} // namespace bcl
