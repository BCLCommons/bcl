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
#include "model/bcl_model_retrieve_data_set_balanced.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "model/bcl_model_retrieve_data_set_from_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_retrieve_data_set_balanced.cpp
  //!
  //! @author mendenjl
  //! @date Feb 14, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelRetrieveDataSetBalanced :
    public ExampleInterface
  {
  public:

    ExampleModelRetrieveDataSetBalanced *Clone() const
    {
      return new ExampleModelRetrieveDataSetBalanced( *this);
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
      model::RetrieveDataSetBalanced default_retrieve_data_set_combined;

      // file containing a data set
      const std::string data_set_filename( AddExampleInputPathToFilename( e_Model, "example_data_set.bcl"));

      // constructor from label, as on the command line
      // this label will create a chemistry::RetrieveQsarDataSet that will load diazepam
      util::ObjectDataLabel diazepam_loader
      (
        "SdfFile(filename=\"" + AddExampleInputPathToFilename( e_Chemistry, "diazepam.sdf") + "\")"
      );

      // retrieve one dataset from a file from a filename
      model::RetrieveDataSetFromFile data_set_from_file( data_set_filename);

      // create a model::RetrieveDataSetBalanced as an implementation of the retrieve data set base
      util::Implementation< model::RetrieveDataSetBase> impl_balanced_data_set
      (
        "Balanced("
        + data_set_from_file.GetString()
        + util::ObjectDataLabel::GetArgumentDelimiter()
        + diazepam_loader.ToString()
        + ")"
      );
      // make sure the implementation could be defined
      BCL_ExampleAssert( impl_balanced_data_set.IsDefined(), true);

      // load in the example data set from a file
      util::ShPtr< descriptor::Dataset> sp_example_data_set
      (
        data_set_from_file.GenerateDataSet()
      );

      // make a label that indicates what the data set holds for features and results
      util::ObjectDataLabel feature_label( "Combine(NRotBond,LogS,Dipole,Polariz,TPSA)");
      util::ObjectDataLabel results_label( "Combine(XlogP,HAcc)");

      // add these data values for diazepam
      linal::Matrix< float> diazepam_features( 1, 5);
      diazepam_features( 0, 0) = 1.0;
      diazepam_features( 0, 1) = -4.14405;
      diazepam_features( 0, 2) = 3.94236;
      diazepam_features( 0, 3) = 32.717;
      diazepam_features( 0, 4) = 32.67;

      linal::Matrix< float> diazepam_results( 1, 2);
      diazepam_results( 0, 0) = 2.92048;
      diazepam_results( 0, 1) = 3.0;

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( model::RetrieveDataSetBalanced().GetAlias(), "Balanced");

      // set features and results labels
      impl_balanced_data_set->SelectFeatures( feature_label);
      impl_balanced_data_set->SelectResults( results_label);

      BCL_ExampleCheck( impl_balanced_data_set->GetFeatureCode().ToString(), feature_label.ToString());
      BCL_ExampleCheck( impl_balanced_data_set->GetResultCode().ToString(), results_label.ToString());

    ////////////////
    // operations //
    ////////////////

      // make sure that the data sets retrieved are correct from the implementation
      if
      (
        BCL_ExampleCheck( impl_balanced_data_set->GenerateDataSet().IsDefined(), true)
        && BCL_ExampleCheck( impl_balanced_data_set->GenerateDataSet()->GetSize(), 6)
        && BCL_ExampleCheck( impl_balanced_data_set->GenerateDataSet()->GetFeatureSize(), 5)
        && BCL_ExampleCheck( impl_balanced_data_set->GenerateDataSet()->GetResultSize(), 2)
      )
      {
        util::ShPtr< descriptor::Dataset> balanced_dataset
        (
          impl_balanced_data_set->GenerateDataSet()
        );
        linal::Matrix< float> balanced_features( balanced_dataset->GetFeaturesReference());
        linal::Matrix< float> balanced_results( balanced_dataset->GetResultsReference());
        linal::Matrix< float> balanced_features_from_file( sp_example_data_set->GetFeaturesReference());
        linal::Matrix< float> balanced_results_from_file( sp_example_data_set->GetResultsReference());

        // make sure that each data element is correct
        // check each row of the features
        BCL_ExampleCheck
        (
          balanced_features.CreateSubMatrix( 1, 5),
          balanced_features_from_file.CreateSubMatrix( 1, 5)
        );
        BCL_ExampleCheck
        (
          balanced_features.CreateSubMatrix( 1, 5, 1),
          diazepam_features
        );
        BCL_ExampleCheck
        (
          balanced_features.CreateSubMatrix( 1, 5, 2),
          balanced_features_from_file.CreateSubMatrix( 1, 5, 1)
        );
        BCL_ExampleCheck
        (
          balanced_features.CreateSubMatrix( 1, 5, 3),
          diazepam_features
        );
        BCL_ExampleCheck
        (
          balanced_features.CreateSubMatrix( 1, 5, 4),
          balanced_features_from_file.CreateSubMatrix( 1, 5, 2)
        );
        BCL_ExampleCheck
        (
          balanced_features.CreateSubMatrix( 1, 5, 5),
          diazepam_features
        );

        // check results
        BCL_ExampleCheck
        (
          balanced_results.CreateSubMatrix( 1, 2),
          balanced_results_from_file.CreateSubMatrix( 1, 2)
        );
        BCL_ExampleCheck
        (
          balanced_results.CreateSubMatrix( 1, 2, 1),
          diazepam_results
        );
        BCL_ExampleCheck
        (
          balanced_results.CreateSubMatrix( 1, 2, 2),
          balanced_results_from_file.CreateSubMatrix( 1, 2, 1)
        );
        BCL_ExampleCheck
        (
          balanced_results.CreateSubMatrix( 1, 2, 3),
          diazepam_results
        );
        BCL_ExampleCheck
        (
          balanced_results.CreateSubMatrix( 1, 2, 4),
          balanced_results_from_file.CreateSubMatrix( 1, 2, 2)
        );
        BCL_ExampleCheck
        (
          balanced_results.CreateSubMatrix( 1, 2, 5),
          diazepam_results
        );
      }

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelRetrieveDataSetBalanced

  const ExampleClass::EnumType ExampleModelRetrieveDataSetBalanced::s_Instance
  (
    GetExamples().AddEnum( ExampleModelRetrieveDataSetBalanced())
  );

} // namespace bcl
