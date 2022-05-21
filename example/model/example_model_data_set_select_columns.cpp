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
#include "model/bcl_model_data_set_select_columns.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_data_set_select_columns.cpp
  //!
  //! @author mendenjl
  //! @date Dec 20, 2010
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelDataSetSelectColumns :
    public ExampleInterface
  {
  public:

    ExampleModelDataSetSelectColumns *Clone() const
    {
      return new ExampleModelDataSetSelectColumns( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      model::DataSetSelectColumns default_column_selector;

      // test constructors that select the 1st, second, or first and third columns
      model::DataSetSelectColumns select_column_0_of_3( 3, storage::Vector< size_t>( 1, 0));
      model::DataSetSelectColumns select_column_1_of_3( 3, storage::Vector< size_t>( 1, 1));
      model::DataSetSelectColumns select_columns_0_and_2_of_3( 3, storage::Vector< size_t>::Create( 0, 2));

    /////////////////
    // data access //
    /////////////////

      // test get input feature size
      BCL_ExampleCheck( model::DataSetSelectColumns().GetInputFeatureSize(), 0);
      BCL_ExampleCheck( model::DataSetSelectColumns( 2).GetInputFeatureSize(), 2);
      BCL_ExampleCheck( model::DataSetSelectColumns( 87).GetInputFeatureSize(), 87);

      // test get output feature size
      BCL_ExampleCheck( model::DataSetSelectColumns().GetOutputFeatureSize(), 0);
      BCL_ExampleCheck( model::DataSetSelectColumns( 2).GetOutputFeatureSize(), 0);
      BCL_ExampleCheck( select_column_0_of_3.GetOutputFeatureSize(), 1);
      BCL_ExampleCheck( select_columns_0_and_2_of_3.GetOutputFeatureSize(), 2);

      // test get column indices
      BCL_ExampleCheck( select_columns_0_and_2_of_3.GetColumnIndices()( 0), 0);
      BCL_ExampleCheck( select_columns_0_and_2_of_3.GetColumnIndices()( 1), 2);
      BCL_ExampleCheck( select_column_1_of_3.GetColumnIndices()( 0), 1);

    ////////////////
    // operations //
    ////////////////

      // make a data set of three columns, then we'll select a subset of those columns
      storage::Vector< linal::Vector< float> > features;

      // size of the data set
      const size_t data_set_size( 20);

      // constant values for the three columns
      const float first_col_value( 0.1), second_col_value( 0.01), third_col_value( 0.001);

      // initialize the data set with 20 points, all the same
      for( size_t i( 0); i < data_set_size; ++i)
      {
        features.PushBack
        (
          linal::MakeVector< float>( first_col_value, second_col_value, third_col_value)
        );
      }

      // use the select columns to convert data sets, test that the correct columns are retrieved
      linal::Matrix< float> features_results_col_0;
      linal::Matrix< float> features_results_col_1;
      linal::Matrix< float> features_results_col_0_and_2;

      ConvertDataSet( features, features_results_col_0, select_column_0_of_3);
      ConvertDataSet( features, features_results_col_1, select_column_1_of_3);
      ConvertDataSet( features, features_results_col_0_and_2, select_columns_0_and_2_of_3);

      // test that the converted data sets are the right size
      BCL_ExampleAssert( features_results_col_0.GetNumberRows(), data_set_size);
      BCL_ExampleAssert( features_results_col_1.GetNumberRows(), data_set_size);
      BCL_ExampleAssert( features_results_col_0_and_2.GetNumberRows(), data_set_size);
      BCL_ExampleAssert( features_results_col_0.GetNumberCols(), 1);
      BCL_ExampleAssert( features_results_col_1.GetNumberCols(), 1);
      BCL_ExampleAssert( features_results_col_0_and_2.GetNumberCols(), 2);

      // get a pointer to the first element in every data set
      const float *first_feature_results_col_0( features_results_col_0[ 0]);
      const float *first_feature_results_col_1( features_results_col_1[ 0]);
      const float *first_features_results_col_0_and_2( features_results_col_0_and_2[ 0]);

      // test that the right features were selected
      BCL_ExampleCheck( *first_feature_results_col_0, first_col_value);
      BCL_ExampleCheck( *first_feature_results_col_1, second_col_value);
      BCL_ExampleCheck( first_features_results_col_0_and_2[ 0], first_col_value);
      BCL_ExampleCheck( first_features_results_col_0_and_2[ 1], third_col_value);

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief convert to data sets of different type
    //! it iterates over all elements of the original data set, and pushes a converted data item in the converted data set
    //! @param ORIGINAL_DATA_SET data set with original data
    //! @param CONVERTED_DATA_SET data set, where the converted items will be added to
    //! @param COLUMN_SELECTOR  selects columnts from the original dataset
    static void ConvertDataSet
    (
      const storage::Vector< linal::Vector< float> > &ORIGINAL_DATA_SET,
      linal::Matrix< float> &CONVERTED_DATA_SET,
      const model::DataSetSelectColumns &COLUMN_SELECTOR
    )
    {
      const size_t number_elements( ORIGINAL_DATA_SET.GetSize());
      CONVERTED_DATA_SET = linal::Matrix< float>( number_elements, COLUMN_SELECTOR.GetOutputFeatureSize());
      for( size_t count( 0); count < number_elements; ++count)
      {
        COLUMN_SELECTOR( ORIGINAL_DATA_SET( count), CONVERTED_DATA_SET[ count]);
      }
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelDataSetSelectColumns

  const ExampleClass::EnumType ExampleModelDataSetSelectColumns::s_Instance
  (
    GetExamples().AddEnum( ExampleModelDataSetSelectColumns())
  );

} // namespace bcl
