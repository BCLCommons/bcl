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
#include "model/bcl_model_objective_function_segment_overlap.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "model/bcl_model_feature_data_set.h"
#include "model/bcl_model_retrieve_data_set_base.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_objective_function_segment_overlap_monitoring_dataset.cpp
  //!
  //! @author mendenjl
  //! @date Mar 18, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelObjectiveFunctionSegmentOverlap :
    public ExampleInterface
  {
  public:

    ExampleModelObjectiveFunctionSegmentOverlap *Clone() const
    {
      return new ExampleModelObjectiveFunctionSegmentOverlap( *this);
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
      // load a dataset subset that contains secondary structure analysis for several soluble proteins
      // the feature is actually the palsse prediction, the result is actually the DSSPStride prediction
      // likewise, we expect the segment overlap to be very high, since both of these methods receive the complete
      // 3d structure.  The actual values in ascii format can be seen in segment overlap.csv
      util::Implementation< model::RetrieveDataSetBase> impl_dataset_retriever
      (
        "Subset(filename="
        + AddExampleInputPathToFilename( e_Model, "segment_overlap.bin")
        + ")"
      );

      // create the dataset
      util::ShPtr< descriptor::Dataset> sp_dataset( impl_dataset_retriever->GenerateDataSet());
      BCL_ExampleIndirectAssert( impl_dataset_retriever.IsDefined(), true, "Could not create dataset retriever");
      // get the palsse-analyzed residue SSEs
      model::FeatureDataSet< float> &palsse_values( sp_dataset->GetFeatures());
      // get the DSSPStride analyzed residue SSEs
      model::FeatureDataSet< float> &dsspstride_values( sp_dataset->GetResults());
      // get the ids
      model::FeatureDataSet< char> &ids( sp_dataset->GetIds());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::ObjectiveFunctionSegmentOverlap function_default;
      function_default.WriteHelp( util::GetLogger(), 0);

      // use the TryRead function to setup this classes preferences
      BCL_ExampleCheck
      (
        function_default.TryRead
        (
          util::ObjectDataLabel
          (
            "",
            function_default.GetAlias(),
            storage::Vector< util::ObjectDataLabel>::Create
            (
              util::ObjectDataLabel( "output sequence info", "True"),
              util::ObjectDataLabel( "output subclass overlaps", "True"),
              util::ObjectDataLabel( "element weight", "True")
            )
          ),
          util::GetLogger()
        ),
        true
      );
      function_default.SetData( dsspstride_values, ids);

    /////////////////
    // data access //
    /////////////////

      // check the directionality: larger should be better
      BCL_ExampleCheck
      (
        model::ObjectiveFunctionSegmentOverlap().GetImprovementType(),
        opti::e_LargerEqualIsBetter
      );

    ///////////////
    // operators //
    ///////////////

      // check the output of the function on identical datasets
      BCL_ExampleCheckWithinTolerance
      (
        function_default( dsspstride_values, dsspstride_values),
        1.0,
        0.001
      );

      // check the output of the function on each disparate pair of datasets
      BCL_ExampleCheckWithinTolerance
      (
        function_default( dsspstride_values, palsse_values),
        0.80882,
        0.001
      );

      // call set data again; reverse the datasets to test that this objective function is asymmetric
      function_default.SetData( palsse_values, ids);
      BCL_ExampleCheckWithinTolerance
      (
        function_default( palsse_values, palsse_values),
        1.0,
        0.001
      );

      // check the output of the function on each disparate pair of datasets
      BCL_ExampleCheckWithinTolerance
      (
        function_default( palsse_values, dsspstride_values),
        0.862783,
        0.001
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelObjectiveFunctionSegmentOverlap

  const ExampleClass::EnumType ExampleModelObjectiveFunctionSegmentOverlap::s_Instance
  (
    GetExamples().AddEnum( ExampleModelObjectiveFunctionSegmentOverlap())
  );

} // namespace bcl
