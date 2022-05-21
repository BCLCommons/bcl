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
#include "model/bcl_model_objective_function_constant.h"

// includes from bcl - sorted alphabetically
#include "model/bcl_model_feature_data_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_objective_function_constant.cpp
  //!
  //! @author butkiem1, mendenjl
  //! @date May 12, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelObjectiveFunctionConstant :
    public ExampleInterface
  {
  public:

    ExampleModelObjectiveFunctionConstant *Clone() const
    {
      return new ExampleModelObjectiveFunctionConstant( *this);
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
      float col_0_largest[ 3] = { 5.0, 3.0, 1.0 };
      float col_1_largest[ 3] = { 3.0, 5.0, 1.0 };
      float col_2_largest[ 3] = { 1.0, 3.0, 5.0 };

      // make a matrix of experimental data
      linal::Matrix< float> experimental_data( 3, 3);
      experimental_data.ReplaceRow( 0, linal::Vector< float>( 3, col_0_largest));
      experimental_data.ReplaceRow( 1, linal::Vector< float>( 3, col_1_largest));
      experimental_data.ReplaceRow( 2, linal::Vector< float>( 3, col_2_largest));
      model::FeatureDataSet< float> experimental_dataset( experimental_data);

      // make a matrix of predicted data
      linal::Matrix< float> awfully_predicted_data( 3, 3);
      awfully_predicted_data.ReplaceRow( 0, linal::Vector< float>( 3, col_2_largest));
      awfully_predicted_data.ReplaceRow( 1, linal::Vector< float>( 3, col_0_largest));
      awfully_predicted_data.ReplaceRow( 2, linal::Vector< float>( 3, col_1_largest));
      model::FeatureDataSet< float> awfully_predicted_dataset( awfully_predicted_data);

      // make a matrix of predicted data
      linal::Matrix< float> const_data( 3, 3);
      const_data.ReplaceRow( 0, linal::Vector< float>( 3, col_1_largest));
      const_data.ReplaceRow( 1, linal::Vector< float>( 3, col_1_largest));
      const_data.ReplaceRow( 2, linal::Vector< float>( 3, col_1_largest));
      model::FeatureDataSet< float> const_predicted_dataset( const_data);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::ObjectiveFunctionConstant obj_const;

      // constructor with parameters
      const float value( 457.3);
      model::ObjectiveFunctionConstant obj_constant( value, opti::e_LargerEqualIsBetter);

    /////////////////
    // data access //
    /////////////////

      // check improvement type
      BCL_ExampleAssert( obj_constant.GetImprovementType(), opti::e_LargerEqualIsBetter);

    ///////////////
    // operators //
    ///////////////

      // check the output of the function on identical datasets
      BCL_ExampleCheck
      (
        obj_constant( experimental_dataset, experimental_dataset),
        value
      );

      // check the output of the function on disparate pairs of datasets
      BCL_ExampleCheck
      (
        obj_constant( experimental_dataset, const_predicted_dataset),
        value
      );
      BCL_ExampleCheck
      (
        obj_constant( awfully_predicted_dataset, const_predicted_dataset),
        value
      );

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelObjectiveFunctionConstant

  const ExampleClass::EnumType ExampleModelObjectiveFunctionConstant::s_Instance
  (
    GetExamples().AddEnum( ExampleModelObjectiveFunctionConstant())
  );

} // namespace bcl
