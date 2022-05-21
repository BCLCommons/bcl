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
#include "model/bcl_model_retrieve_data_set_from_file.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_retrieve_data_set_from_file.cpp
  //!
  //! @author mendenjl
  //! @date Feb 14, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelRetrieveDataSetFromFile :
    public ExampleInterface
  {
  public:

    ExampleModelRetrieveDataSetFromFile *Clone() const
    {
      return new ExampleModelRetrieveDataSetFromFile( *this);
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

      // default constructor
      model::RetrieveDataSetFromFile default_retrieve_data_set_from_file;

      // file containing a data set
      const std::string data_set_filename( AddExampleInputPathToFilename( e_Model, "example_data_set.bcl"));

      // constructor from a filename
      model::RetrieveDataSetFromFile example_data_set_retriever( data_set_filename);

      // make a range set containing 0 and 2
      math::RangeSet< size_t> zero_two( math::Range< size_t>( 0, 0));
      zero_two += math::Range< size_t>( 2, 2);

      // constructor from a filename and chunk ranges
      model::RetrieveDataSetFromFile example_data_set_0_2_retriever( data_set_filename, 3, zero_two);

      // constructor from a filename and the first of two chunks
      // this demonstrates that the # of chunks need not be the same as the # of data points
      model::RetrieveDataSetFromFile example_data_set_0_1_retriever
      (
        data_set_filename,
        2,
        math::RangeSet< size_t>( math::Range< size_t>( 0, 0))
      );

      // also make sure that we can construct this type as an implementation
      util::Implementation< model::RetrieveDataSetBase> impl_example_data_set_retriever_0_2
      (
        "File(filename=" + data_set_filename + ", number chunks=3, chunks=\"[0,0]+[2,2]\")"
      );

      // load in the example data set from a file
      storage::Vector< storage::VectorND< 2, linal::Vector< double> > > example_data_set;
      {
        io::IFStream input;
        BCL_ExampleMustOpenInputFile( input, data_set_filename);
        io::Serialize::Read( example_data_set, input);
        io::File::CloseClearFStream( input);
        BCL_ExampleIndirectAssert( example_data_set.GetSize(), 3, "data set loading directly from file");
      }

      // get the actual feature and result size
      const size_t actual_data_size( example_data_set.GetSize());
      const size_t actual_feature_size( example_data_set.FirstElement().First().GetSize());
      const size_t actual_result_size( example_data_set.FirstElement().Second().GetSize());

      // convert the example_data_set into features and results stored in a linal matrix
      linal::Matrix< float> features_converted( actual_data_size, actual_feature_size);
      linal::Matrix< float> results_converted( actual_data_size, actual_result_size);

      for( size_t feature_number( 0); feature_number < actual_data_size; ++feature_number)
      {
        // feature conversion
        for( size_t col( 0); col < actual_feature_size; ++col)
        {
          features_converted( feature_number, col) = example_data_set( feature_number).First()( col);
        }
        // result conversion
        for( size_t col( 0); col < actual_result_size; ++col)
        {
          results_converted( feature_number, col) = example_data_set( feature_number).Second()( col);
        }
      }

      // make a label that indicates what the data set holds for features and results
      util::ObjectDataLabel feature_label( "Features(0,1)");
      util::ObjectDataLabel results_label( "Results(5)");

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( model::RetrieveDataSetFromFile().GetAlias(), "File");

      // set features and results labels
      model::RetrieveDataSetFromFile example_data_set_retriever_b( data_set_filename);
      example_data_set_retriever_b.SelectFeatures( feature_label);
      example_data_set_retriever_b.SelectResults( results_label);

      BCL_ExampleIndirectCheck
      (
        example_data_set_retriever_b.GetFeatureCode().ToString(),
        feature_label.ToString(),
        "SelectFeatures"
      );
      BCL_ExampleIndirectCheck
      (
        example_data_set_retriever_b.GetResultCode().ToString(),
        results_label.ToString(),
        "SelectResults"
      );

    ////////////////
    // operations //
    ////////////////

      // make sure that the data sets retrieved are correct
      if( BCL_ExampleCheck( example_data_set_retriever.GenerateDataSet().IsDefined(), true))
      {
        // make sure that the features and results are the same as expected
        BCL_ExampleCheck( example_data_set_retriever.GenerateDataSet()->GetFeaturesReference(), features_converted);
        BCL_ExampleCheck( example_data_set_retriever.GenerateDataSet()->GetResultsReference(), results_converted);
      }

      // get the 0 and 2 row from features/results_converted
      linal::Matrix< float> features_converted_0_2( 2, actual_feature_size);
      linal::Matrix< float> results_converted_0_2( 2, actual_result_size);

      for( size_t feature_number( 0); feature_number < 2; ++feature_number)
      {
        // feature conversion
        for( size_t col( 0); col < actual_feature_size; ++col)
        {
          features_converted_0_2( feature_number, col) = example_data_set( feature_number * 2).First()( col);
        }
        // result conversion
        for( size_t col( 0); col < actual_result_size; ++col)
        {
          results_converted_0_2( feature_number, col) = example_data_set( feature_number * 2).Second()( col);
        }
      }

      if( BCL_ExampleCheck( example_data_set_0_2_retriever.GenerateDataSet().IsDefined(), true))
      {
        // make sure that the features and results are the same as expected
        BCL_ExampleCheck( example_data_set_0_2_retriever.GenerateDataSet()->GetFeaturesReference(), features_converted_0_2);
        BCL_ExampleCheck( example_data_set_0_2_retriever.GenerateDataSet()->GetResultsReference(), results_converted_0_2);
      }

      // get the 0 and 2 row from features/results_converted
      linal::Matrix< float> features_converted_0_1( 2, actual_feature_size, features_converted.Begin());
      linal::Matrix< float> results_converted_0_1( 2, actual_result_size, results_converted.Begin());

      // make sure that the data sets retrieved are correct
      if( BCL_ExampleCheck( example_data_set_0_1_retriever.GenerateDataSet().IsDefined(), true))
      {
        // make sure that the features and results are the same as expected
        BCL_ExampleCheck( example_data_set_0_1_retriever.GenerateDataSet()->GetFeaturesReference(), features_converted_0_1);
        BCL_ExampleCheck( example_data_set_0_1_retriever.GenerateDataSet()->GetResultsReference(), results_converted_0_1);
      }

      // make sure that the data sets retrieved are correct from the implementation
      if
      (
        BCL_ExampleCheck( impl_example_data_set_retriever_0_2.IsDefined(), true)
        && BCL_ExampleCheck( impl_example_data_set_retriever_0_2->GenerateDataSet().IsDefined(), true)
      )
      {
        // make sure that the features and results are the same as expected
        BCL_ExampleCheck( impl_example_data_set_retriever_0_2->GenerateDataSet()->GetFeaturesReference(), features_converted_0_2);
        BCL_ExampleCheck( impl_example_data_set_retriever_0_2->GenerateDataSet()->GetResultsReference(), results_converted_0_2);
      }

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelRetrieveDataSetFromFile

  const ExampleClass::EnumType ExampleModelRetrieveDataSetFromFile::s_Instance
  (
    GetExamples().AddEnum( ExampleModelRetrieveDataSetFromFile())
  );

} // namespace bcl
