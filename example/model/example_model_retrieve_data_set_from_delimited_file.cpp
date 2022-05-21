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
#include "model/bcl_model_retrieve_data_set_from_delimited_file.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_retrieve_data_set_from_delimited_file.cpp
  //!
  //! @author butkiem1, mendenjl
  //! @date Sep 18, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelRetrieveDataSetFromDelimitedFile :
    public ExampleInterface
  {
  public:

    ExampleModelRetrieveDataSetFromDelimitedFile *Clone() const
    {
      return new ExampleModelRetrieveDataSetFromDelimitedFile( *this);
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
      model::RetrieveDataSetFromDelimitedFile default_retrieve_data_set_from_file;

      // file containing a data set
      const std::string data_set_filename( AddExampleInputPathToFilename( e_Model, "example_data_set.csv"));

      // constructor from a filename
      model::RetrieveDataSetFromDelimitedFile example_data_set_retriever( data_set_filename, 1, 2);

      // constructor from a filename and chunk ranges
      model::RetrieveDataSetFromDelimitedFile example_data_set_0_2_retriever( data_set_filename, 1, 0);

      // constructor from a filename and the first of two chunks
      // this demonstrates that the # of chunks need not be the same as the # of data points
      model::RetrieveDataSetFromDelimitedFile example_data_set_0_1_retriever
      (
        data_set_filename,
        1,
        0
      );

      // also make sure that we can construct this type as an implementation
      util::Implementation< model::RetrieveDataSetBase> impl_example_data_set_retriever_0_2
      (
        "Csv(filename=" + data_set_filename + ", number id chars = 2)"
      );

      // make a label that indicates what the data set holds for features and results
      util::ObjectDataLabel feature_label( "Features(numbers)");
      util::ObjectDataLabel results_label( "Results");

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( model::RetrieveDataSetFromDelimitedFile().GetAlias(), "Csv");

      // set features and results labels
      default_retrieve_data_set_from_file.SelectFeatures( feature_label);
      default_retrieve_data_set_from_file.SelectResults( results_label);

      BCL_ExampleCheck( default_retrieve_data_set_from_file.GetFeatureCode().ToString(), feature_label.ToString());
      BCL_ExampleCheck( default_retrieve_data_set_from_file.GetResultCode().ToString(), results_label.ToString());

    ////////////////
    // operations //
    ////////////////

      util::ShPtr< descriptor::Dataset> dataset( impl_example_data_set_retriever_0_2->GenerateDataSet());

      if( BCL_ExampleCheck( dataset.IsDefined(), true))
      {
        // make sure that the features and results are the same as expected
        BCL_ExampleCheck( dataset->GetFeaturesReference().GetNumberCols(), size_t( 3));
        BCL_ExampleCheck( dataset->GetFeaturesReference().GetNumberRows(), size_t( 5));
        BCL_ExampleCheck( dataset->GetResultSize(), size_t( 1));
        BCL_ExampleCheck( dataset->GetIdSize(), size_t( 2));
        BCL_ExampleCheck( dataset->GetIdsReference()( 1, 1), 'D');
      }

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelRetrieveDataSetFromDelimitedFile

  const ExampleClass::EnumType ExampleModelRetrieveDataSetFromDelimitedFile::s_Instance
  (
    GetExamples().AddEnum( ExampleModelRetrieveDataSetFromDelimitedFile())
  );

} // namespace bcl
