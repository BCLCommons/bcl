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
#include "model/bcl_model_objective_function_categorical_max.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_objective_function_categorical_max_monitoring_dataset.cpp
  //!
  //! @author mendenjl
  //! @date Jul 22, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelObjectiveFunctionCategoricalMax :
    public ExampleInterface
  {
  public:

    ExampleModelObjectiveFunctionCategoricalMax *Clone() const
    {
      return new ExampleModelObjectiveFunctionCategoricalMax( *this);
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
      // create 3 arrays with different largest columns
      float col_0_largest[ 3] = { 1.0, 0.4, 0.2 };
      float col_1_largest[ 3] = { 0.4, 1.0, 0.2 };
      float col_2_largest[ 3] = { 0.2, 0.4, 1.2 };

      // make a matrix of experimental data
      linal::Matrix< float> experimental_data( 6, 3);
      experimental_data.ReplaceRow( 0, linal::Vector< float>( 3, col_0_largest));
      experimental_data.ReplaceRow( 1, linal::Vector< float>( 3, col_1_largest));
      experimental_data.ReplaceRow( 2, linal::Vector< float>( 3, col_2_largest));
      experimental_data.ReplaceRow( 3, linal::Vector< float>( 3, col_0_largest));
      experimental_data.ReplaceRow( 4, linal::Vector< float>( 3, col_1_largest));
      experimental_data.ReplaceRow( 5, linal::Vector< float>( 3, col_2_largest));
      model::FeatureDataSet< float> experimental_dataset( experimental_data);

      // create a matrix of data that is never in agreement with the experimental data
      linal::Matrix< float> awfully_predicted_data( 6, 3);
      awfully_predicted_data.ReplaceRow( 0, linal::Vector< float>( 3, col_2_largest));
      awfully_predicted_data.ReplaceRow( 1, linal::Vector< float>( 3, col_0_largest));
      awfully_predicted_data.ReplaceRow( 2, linal::Vector< float>( 3, col_1_largest));
      awfully_predicted_data.ReplaceRow( 3, linal::Vector< float>( 3, col_2_largest));
      awfully_predicted_data.ReplaceRow( 4, linal::Vector< float>( 3, col_0_largest));
      awfully_predicted_data.ReplaceRow( 5, linal::Vector< float>( 3, col_1_largest));
      model::FeatureDataSet< float> awfully_predicted_dataset( awfully_predicted_data);

      // create a matrix of data that is correct 1/3 of the time because it always guesses the same value
      linal::Matrix< float> const_predicted_data( 6, 3);
      const_predicted_data.ReplaceRow( 0, linal::Vector< float>( 3, col_1_largest));
      const_predicted_data.ReplaceRow( 1, linal::Vector< float>( 3, col_1_largest));
      const_predicted_data.ReplaceRow( 2, linal::Vector< float>( 3, col_1_largest));
      const_predicted_data.ReplaceRow( 3, linal::Vector< float>( 3, col_1_largest));
      const_predicted_data.ReplaceRow( 4, linal::Vector< float>( 3, col_1_largest));
      const_predicted_data.ReplaceRow( 5, linal::Vector< float>( 3, col_1_largest));
      model::FeatureDataSet< float> const_predicted_dataset( const_predicted_data);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::ObjectiveFunctionCategoricalMax function_default;
      function_default.SetData( experimental_dataset);

    /////////////////
    // data access //
    /////////////////

      // check the directionality: larger should be better
      BCL_ExampleCheck
      (
        model::ObjectiveFunctionCategoricalMax().GetImprovementType(),
        opti::e_LargerEqualIsBetter
      );

    ///////////////
    // operators //
    ///////////////

      // check the output of the function on identical datasets
      BCL_ExampleCheckWithinTolerance
      (
        function_default( experimental_dataset, experimental_dataset),
        1.0,
        0.001
      );

      // check the output of the function on each disparate pair of datasets
      BCL_ExampleCheckWithinTolerance
      (
        function_default( experimental_dataset, awfully_predicted_dataset),
        0.0,
        0.001
      );

      // check the output of the function on each disparate pair of datasets
      BCL_ExampleCheckWithinTolerance
      (
        function_default( experimental_dataset, const_predicted_dataset),
        float( 2.0 / 6.0),
        0.001
      );

      // check the output of the function on each disparate pair of datasets
      function_default.SetData( awfully_predicted_dataset);
      BCL_ExampleCheckWithinTolerance
      (
        function_default( awfully_predicted_dataset, const_predicted_dataset),
        float( 2.0 / 6.0),
        0.001
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelObjectiveFunctionCategoricalMax

  const ExampleClass::EnumType ExampleModelObjectiveFunctionCategoricalMax::s_Instance
  (
    GetExamples().AddEnum( ExampleModelObjectiveFunctionCategoricalMax())
  );

} // namespace bcl
