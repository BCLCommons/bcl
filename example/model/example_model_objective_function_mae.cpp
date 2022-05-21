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
#include "model/bcl_model_objective_function_mae.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_objective_function_mae_monitoring_dataset.cpp
  //!
  //! @author mendenjl
  //! @date Aug 11, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelObjectiveFunctionMae :
    public ExampleInterface
  {
  public:

    ExampleModelObjectiveFunctionMae *Clone() const
    {
      return new ExampleModelObjectiveFunctionMae( *this);
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

    //! @brief compute absolute difference between two vectors

    int Run() const
    {
      // create 3 arrays with different largest columns
      float col_a[ 3] = { 6.0, 3.0, 1.0 };
      float col_b[ 3] = { 3.0, 6.0, 1.0 };
      float col_c[ 3] = { 1.0, 3.0, 6.0 };

      // make a matrix of experimental data
      linal::Matrix< float> experimental_data( 3, 3);
      experimental_data.ReplaceRow( 0, linal::Vector< float>( 3, col_a));
      experimental_data.ReplaceRow( 1, linal::Vector< float>( 3, col_b));
      experimental_data.ReplaceRow( 2, linal::Vector< float>( 3, col_c));
      model::FeatureDataSet< float> experimental_dataset( experimental_data);

      // create a matrix of predicted data
      linal::Matrix< float> predicted_data( 3, 3);
      predicted_data.ReplaceRow( 0, linal::Vector< float>( 3, col_c));
      predicted_data.ReplaceRow( 1, linal::Vector< float>( 3, col_a));
      predicted_data.ReplaceRow( 2, linal::Vector< float>( 3, col_b));
      model::FeatureDataSet< float> predicted_dataset( predicted_data);

      // sum of absolute error col_a <> col_b
      const float sae_ab
      (
          math::Absolute( col_a[ 0] - col_b[ 0])
        + math::Absolute( col_a[ 1] - col_b[ 1])
        + math::Absolute( col_a[ 2] - col_b[ 2])
      );
      // sum of absolute error col_a <> col_c
      const float sae_ac
      (
          math::Absolute( col_a[ 0] - col_c[ 0])
        + math::Absolute( col_a[ 1] - col_c[ 1])
        + math::Absolute( col_a[ 2] - col_c[ 2])
      );
      // sum of absolute error col_b <> col_b
      const float sae_bc
      (
          math::Absolute( col_b[ 0] - col_c[ 0])
        + math::Absolute( col_b[ 1] - col_c[ 1])
        + math::Absolute( col_b[ 2] - col_c[ 2])
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // objective functions are usually constructed from a label, as shown here

      util::Implementation< model::ObjectiveFunctionInterface> mae( "MAE");

      // make sure the implementation could be created, otherwise the later checks would cause a crash
      BCL_ExampleAssert( mae.IsDefined(), true);

    /////////////////
    // data access //
    /////////////////

      // check the directionality: smaller should be better
      BCL_ExampleCheck
      (
        mae->GetImprovementType(),
        opti::e_SmallerEqualIsBetter
      );

    ///////////////
    // operators //
    ///////////////

      // check the output of the function on identical datasets
      BCL_ExampleCheck
      (
        mae->operator()( experimental_dataset, experimental_dataset),
        0.0
      );

      // check the output of the function on disparate pairs of datasets
      BCL_ExampleCheck
      (
        math::EqualWithinTolerance
        (
          mae->operator()( experimental_dataset, predicted_dataset),
          float( sae_ab + sae_ac + sae_bc) / float( 3.0) / float( 3.0)
        ),
        true
      );

      // check symmetry
      BCL_ExampleCheck
      (
        math::EqualWithinTolerance
        (
          mae->operator()( predicted_dataset, experimental_dataset),
          float( sae_ab + sae_ac + sae_bc) / float( 3.0) / float( 3.0)
        ),
        true
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelObjectiveFunctionMae

  const ExampleClass::EnumType ExampleModelObjectiveFunctionMae::s_Instance
  (
    GetExamples().AddEnum( ExampleModelObjectiveFunctionMae())
  );

} // namespace bcl
