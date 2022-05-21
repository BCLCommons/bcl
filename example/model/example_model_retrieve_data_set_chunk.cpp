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
#include "model/bcl_model_retrieve_data_set_chunk.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "model/bcl_model_retrieve_data_set_from_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_retrieve_data_set_chunk.cpp
  //!
  //! @author butkiem1
  //! @date Feb 14, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelRetrieveDataSetChunk :
    public ExampleInterface
  {
  public:

    ExampleModelRetrieveDataSetChunk *Clone() const
    {
      return new ExampleModelRetrieveDataSetChunk( *this);
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
      model::RetrieveDataSetChunk default_retrieve_data_set_chunk;

      // file containing a data set
      const std::string data_set_filename( AddExampleInputPathToFilename( e_Model, "example_data_set.bcl"));

      // constructor from label, as on the command line
      // this label will create a chemistry::RetrieveQsarDataSet that will load 5 molecules
      util::ObjectDataLabel molecule_loader
      (
        "SdfFile(filename=\"" + AddExampleInputPathToFilename( e_Chemistry, "test_set.5_structures.sdf") + "\")"
      );

      // retrieve one dataset from a file from a filename
      model::RetrieveDataSetFromFile data_set_from_file( data_set_filename);

      // create a model::RetrieveDataSetChunk as an implementation of the retrieve data set base
      util::Implementation< model::RetrieveDataSetBase> impl_chunk_data_set
      (
        "Chunks(chunks=\"[0,2)\",number chunks=8,dataset="
        " Combined("
        + data_set_from_file.GetString()
        + util::ObjectDataLabel::GetArgumentDelimiter()
        + molecule_loader.ToString()
        + ")"
        + ")"
      );

      // create a model::RetrieveDataSetChunk as an implementation of the retrieve data set base
      util::Implementation< model::RetrieveDataSetBase> impl_chunk_data_set_too_many
      (
        "Chunks(chunks=\"[2,10)\",number chunks=8,dataset="
        " Combined("
        + data_set_from_file.GetString()
        + util::ObjectDataLabel::GetArgumentDelimiter()
        + molecule_loader.ToString()
        + ")"
        + ")"
      );

      // create a model::RetrieveDataSetChunk as an implementation of the retrieve data set base
      util::Implementation< model::RetrieveDataSetBase> impl_chunk_data_set_everything
      (
        "Chunks(chunks=\"\",number chunks=8,dataset="
        " Combined("
        + data_set_from_file.GetString()
        + util::ObjectDataLabel::GetArgumentDelimiter()
        + molecule_loader.ToString()
        + ")"
        + ")"
      );

      // make a label that indicates what the data set holds for features and results
      util::ObjectDataLabel feature_label( "Combine(NRotBond,LogS,Dipole,Polariz,TPSA)");
      util::ObjectDataLabel results_label( "Combine(XlogP,HAcc)");

      // make sure the implementation could be defined
      BCL_ExampleAssert( impl_chunk_data_set.IsDefined(), true);

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( model::RetrieveDataSetChunk().GetAlias(), "Chunks");

      // set features and results labels
      impl_chunk_data_set->SelectFeatures( feature_label);
      impl_chunk_data_set->SelectResults( results_label);

      impl_chunk_data_set_too_many->SelectFeatures( feature_label);
      impl_chunk_data_set_too_many->SelectResults( results_label);

      impl_chunk_data_set_everything->SelectFeatures( feature_label);
      impl_chunk_data_set_everything->SelectResults( results_label);

      BCL_ExampleCheck( impl_chunk_data_set->GetFeatureCode().ToString(), feature_label.ToString());
      BCL_ExampleCheck( impl_chunk_data_set->GetResultCode().ToString(), results_label.ToString());

    ////////////////
    // operations //
    ////////////////

      util::ShPtr< descriptor::Dataset> chunked_dataset
      (
        impl_chunk_data_set->GenerateDataSet()
      );

      BCL_ExampleCheck( chunked_dataset->GetSize(), 2);

      util::ShPtr< descriptor::Dataset> chunked_dataset_too_many
      (
        impl_chunk_data_set_too_many->GenerateDataSet()
      );

      BCL_ExampleCheck( chunked_dataset_too_many->GetSize(), 6);

      util::ShPtr< descriptor::Dataset> chunked_dataset_everything
      (
        impl_chunk_data_set_everything->GenerateDataSet()
      );

      BCL_ExampleCheck( chunked_dataset_everything->GetSize(), 8);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelRetrieveDataSetChunk

  const ExampleClass::EnumType ExampleModelRetrieveDataSetChunk::s_Instance
  (
    GetExamples().AddEnum( ExampleModelRetrieveDataSetChunk())
  );

} // namespace bcl
