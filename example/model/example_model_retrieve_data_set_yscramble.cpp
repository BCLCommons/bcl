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
#include "model/bcl_model_retrieve_data_set_yscramble.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "linal/bcl_linal_vector_operations.h"
#include "linal/bcl_linal_vector_reference.h"
#include "model/bcl_model_retrieve_data_set_from_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_retrieve_data_set_yscramble.cpp
  //!
  //! @author mendenjl
  //! @date Apr 29, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelRetrieveDataSetYscramble :
    public ExampleInterface
  {
  public:

    ExampleModelRetrieveDataSetYscramble *Clone() const
    {
      return new ExampleModelRetrieveDataSetYscramble( *this);
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
      model::RetrieveDataSetYscramble default_retrieve_data_set_randomized;

      // file containing a data set
      const std::string data_set_filename( AddExampleInputPathToFilename( e_Model, "example_data_set.bcl"));

      // retrieve one dataset from a file from a filename
      model::RetrieveDataSetFromFile data_set_from_file( data_set_filename);

      // create a model::RetrieveDataSetYscramble as an implementation of the retrieve data set base
      util::Implementation< model::RetrieveDataSetBase> impl_randomize
      (
        "YScramble(" + data_set_from_file.GetString() + ")"
      );

      // make a label that indicates what the data set holds for features and results
      util::ObjectDataLabel feature_label( "Combine(0,1,2,3,4)");

      // make sure the implementation could be defined
      BCL_ExampleAssert( impl_randomize.IsDefined(), true);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( model::RetrieveDataSetYscramble().GetAlias(), "YScramble");

      // set features and results labels
      impl_randomize->SelectFeatures( feature_label);
      data_set_from_file.SelectFeatures( feature_label);

      BCL_ExampleCheck( impl_randomize->GetFeatureCode(), feature_label);

    ////////////////
    // operations //
    ////////////////

      // generate the dataset, random order and non-random order
      util::ShPtr< descriptor::Dataset> randomized_dataset( impl_randomize->GenerateDataSet());
      util::ShPtr< descriptor::Dataset> unrandomized_dataset( data_set_from_file.GenerateDataSet());

      if( BCL_ExampleIndirectCheck( randomized_dataset->GetSize(), unrandomized_dataset->GetSize(), "GenerateDataSet"))
      {
        // check that the generated features are the same, just ordered differently
        linal::Matrix< float> features_rndm( randomized_dataset->GetFeaturesReference());
        linal::Matrix< float> results_rndm( randomized_dataset->GetResultsReference());
        linal::Matrix< float> features_norm( unrandomized_dataset->GetFeaturesReference());
        linal::Matrix< float> results_norm( unrandomized_dataset->GetResultsReference());
        BCL_ExampleIndirectCheck( features_rndm, features_norm, "YScrambling should preserve feature order");
        BCL_ExampleIndirectCheck( results_rndm == results_norm, false, "YScrambling should not preserve result order");
        storage::Vector< linal::Vector< float> >
          results_rndm_vecs( randomized_dataset->GetSize()), results_norm_vecs( unrandomized_dataset->GetSize());
        size_t n_features( randomized_dataset->GetSize());
        for( size_t i( 0); i < n_features; ++i)
        {
          results_rndm_vecs( i) = results_rndm.GetRow( i);
          results_norm_vecs( i) = results_norm.GetRow( i);
        }
        // sort both features and results
        results_rndm_vecs.Sort( std::less< linal::Vector< float> >());
        results_norm_vecs.Sort( std::less< linal::Vector< float> >());
        BCL_ExampleIndirectCheck( results_rndm_vecs, results_norm_vecs, "YScrambling preserves rows after sorting");
      }

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelRetrieveDataSetYscramble

  const ExampleClass::EnumType ExampleModelRetrieveDataSetYscramble::s_Instance
  (
    GetExamples().AddEnum( ExampleModelRetrieveDataSetYscramble())
  );

} // namespace bcl
