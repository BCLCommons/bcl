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
#include "model/bcl_model_retrieve_data_set_by_id.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_retrieve_data_set_by_id.cpp
  //!
  //! @author mendenjl
  //! @date Oct 26, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelRetrieveDataSetById :
    public ExampleInterface
  {
  public:

    ExampleModelRetrieveDataSetById *Clone() const
    {
      return new ExampleModelRetrieveDataSetById( *this);
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

      // file containing a data set
      const std::string data_set_filename( AddExampleInputPathToFilename( e_Model, "example_data_set.csv"));

      // create an implementation of the dataset retriever to retrieve from the example data set
      const std::string model_directory( AddExampleInputPathToFilename( e_Model, ""));
      util::Implementation< model::RetrieveDataSetBase> dataset_retriever
      (
        "Csv(filename=" + data_set_filename + ", number id chars = 2)"
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test constructor from the desired ranges and implementation
      model::RetrieveDataSetById dataset_ab_ef
      (
        storage::Set< std::string>::Create( "AB", "EF"),
        dataset_retriever
      );
      dataset_ab_ef.SelectIds( util::ObjectDataLabel( "Combine(0,1)"));

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      // create the reduced data set
      util::ShPtr< descriptor::Dataset> reduced_data_set( dataset_ab_ef.GenerateDataSet());

      // there are 3 points with the desired ids (two repeats of AB)
      BCL_ExampleIndirectAssert( reduced_data_set->GetSize(), 3, "GenerateDataSet()");

      // check that the ids are as expected
      BCL_ExampleIndirectCheck
      (
        std::string( reduced_data_set->GetIdsReference()[ 0], size_t( 2)),
        "AB",
        "GenerateDataSet() selects the correct ids"
      );
      BCL_ExampleIndirectCheck
      (
        std::string( reduced_data_set->GetIdsReference()[ 1], size_t( 2)),
        "EF",
        "GenerateDataSet() selects the correct ids"
      );
      BCL_ExampleIndirectCheck
      (
        std::string( reduced_data_set->GetIdsReference()[ 2], size_t( 2)),
        "AB",
        "GenerateDataSet() selects the correct ids"
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

  }; //end ExampleModelRetrieveDataSetById

  const ExampleClass::EnumType ExampleModelRetrieveDataSetById::s_Instance
  (
    GetExamples().AddEnum( ExampleModelRetrieveDataSetById())
  );

} // namespace bcl
