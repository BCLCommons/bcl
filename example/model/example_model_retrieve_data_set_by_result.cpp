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
#include "model/bcl_model_retrieve_data_set_by_result.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_retrieve_data_set_by_result.cpp
  //!
  //! @author mendenjl
  //! @date Jul 19, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelRetrieveDataSetByResult :
    public ExampleInterface
  {
  public:

    ExampleModelRetrieveDataSetByResult *Clone() const
    {
      return new ExampleModelRetrieveDataSetByResult( *this);
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
      // create a range set that will select everything between -1.1 - -0.8 and 0.8 to 1.1
      math::RangeSet< float> range_set_below_above_0_point_8;
      const std::string range_set_as_string( "[-1.1,-0.8]+[0.8,1.1]");
      range_set_below_above_0_point_8.FromString( range_set_as_string, util::GetLogger());

      // create an implementation of the dataset retriever to retrieve from the example data set
      const std::string model_directory( AddExampleInputPathToFilename( e_Model, ""));
      util::Implementation< model::RetrieveDataSetBase>
        dataset_retriever( "File(filename=" + model_directory + "example_data_set_score.bcl)");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test constructor from the desired ranges and implementation
      model::RetrieveDataSetByResult data_set_reduced_excluding_results_near_0
      (
        range_set_below_above_0_point_8,
        dataset_retriever
      );

    /////////////////
    // data access //
    /////////////////

      // test get ranges
      BCL_ExampleCheck( data_set_reduced_excluding_results_near_0.GetResultRanges().IsWithin( 0.79), false);
      BCL_ExampleCheck( data_set_reduced_excluding_results_near_0.GetResultRanges().IsWithin( 0.81), true);
      BCL_ExampleCheck( data_set_reduced_excluding_results_near_0.GetResultRanges().IsWithin( 0.8), true);

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      // create the reduced data set
      util::ShPtr< descriptor::Dataset>
        reduced_data_set( data_set_reduced_excluding_results_near_0.GenerateDataSet());

      // we should get 16 points within the desired ranges
      BCL_ExampleIndirectCheck( reduced_data_set->GetSize(), 16, "GenerateDataSet()");

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelRetrieveDataSetByResult

  const ExampleClass::EnumType ExampleModelRetrieveDataSetByResult::s_Instance
  (
    GetExamples().AddEnum( ExampleModelRetrieveDataSetByResult())
  );

} // namespace bcl
