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
#include "model/bcl_model_retrieve_data_set_rescaled.h"

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
  //! @example example_model_retrieve_data_set_rescaled.cpp
  //!
  //! @author mendenjl
  //! @date Jul 06, 2017
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelRetrieveDataSetRescaled :
    public ExampleInterface
  {
  public:

    ExampleModelRetrieveDataSetRescaled *Clone() const
    {
      return new ExampleModelRetrieveDataSetRescaled( *this);
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
      model::RetrieveDataSetRescaled default_retrieve_data_set_randomized;

      // file containing a data set
      const std::string data_set_filename( AddExampleInputPathToFilename( e_Model, "example_data_set.bcl"));

      // retrieve one dataset from a file from a filename
      model::RetrieveDataSetFromFile data_set_from_file( data_set_filename);

      // create a model::RetrieveDataSetRescaled as an implementation of the retrieve data set base
      util::Implementation< model::RetrieveDataSetBase> impl_rescale
      (
        "Randomize(" + data_set_from_file.GetString() + ")"
      );

      // make a label that indicates what the data set holds for features and results
      util::ObjectDataLabel feature_label( "Combine(0,1,2,3,4)");

      // make sure the implementation could be defined
      BCL_ExampleAssert( impl_rescale.IsDefined(), true);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( model::RetrieveDataSetRescaled().GetAlias(), "Rescale");

      // set features and results labels
      impl_rescale->SelectFeatures( feature_label);
      data_set_from_file.SelectFeatures( feature_label);

      BCL_ExampleCheck( impl_rescale->GetFeatureCode(), feature_label);

    ////////////////
    // operations //
    ////////////////

      // generate the dataset, random order and non-random order
      util::ShPtr< descriptor::Dataset> rescaled_dataset( impl_rescale->GenerateDataSet());

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelRetrieveDataSetRescaled

  const ExampleClass::EnumType ExampleModelRetrieveDataSetRescaled::s_Instance
  (
    GetExamples().AddEnum( ExampleModelRetrieveDataSetRescaled())
  );

} // namespace bcl
