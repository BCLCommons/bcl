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
#include "model/bcl_model_retrieve_dataset_subset.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically
#include <stdio.h>

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_retrieve_dataset_subset.cpp
  //!
  //! @author butkiem1
  //! @date Nov 17, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelRetrieveDatasetSubset :
    public ExampleInterface
  {
  public:

    ExampleModelRetrieveDatasetSubset *Clone() const
    {
      return new ExampleModelRetrieveDatasetSubset( *this);
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
      model::RetrieveDatasetSubset retriever;

      // get a bin file, precreated with the first 100 actives from AID891 and some 2D descriptors
      const std::string in_filename( AddExampleInputPathToFilename( e_Model, "aid891_100_actives.bin"));

      // constructor with bin file filename as input
      model::RetrieveDatasetSubset retrieve_dataset_subset
      (
        in_filename,
        size_t( 10),
        math::RangeSet< size_t>( math::Range< size_t>( 0, 1))
      );

    /////////////////
    // data access //
    /////////////////

      // create a model::RetrieveDatasetSubset as an implementation of the retrieve data set base
      util::Implementation< model::RetrieveDataSetBase> impl_chunk_data_set( "Subset( filename=" + in_filename + ")");

      // create a model::RetrieveDatasetSubset as an implementation of the retrieve data set base
      util::Implementation< model::RetrieveDataSetBase> impl_chunk_data_set_first_chunk
      (
        "Subset( filename=" + in_filename + ", number chunks=10, chunks=[0])"
      );

    ////////////////
    // operations //
    ////////////////

      BCL_ExampleCheck( retriever.GetAlias(), "Subset");

      // retrieve all data from bin file and for number of rows
      util::ShPtr< descriptor::Dataset> data_set
      (
        impl_chunk_data_set->GenerateDataSet()
      );

      BCL_ExampleCheck( data_set->GetSize(), 100);

      // retrieve first of ten chunks from bin file and for number of rows
      util::ShPtr< descriptor::Dataset> data_set_first_chunk
      (
        impl_chunk_data_set_first_chunk->GenerateDataSet()
      );

      BCL_ExampleCheck( data_set_first_chunk->GetSize(), 10);

    ///////////////
    // operators //
    ///////////////

      // create an output filename for bin files
      const std::string out_filename_bin
      (
        AddExampleOutputPathToFilename( retriever, "aid891_actives.model_retrieve_dataset_subset.bin")
      );
      const std::string out_filename_csv
      (
        AddExampleOutputPathToFilename( retriever, "aid891_actives.model_retrieve_dataset_subset.csv")
      );

      model::RetrieveDatasetSubset::StoreMasterDataset( out_filename_bin, retrieve_dataset_subset);

      // check that the files are equal, if so, remove the generated file to avoid cluttering the svn
      if( BCL_ExampleCheck( io::File::BinaryFilesMatch( out_filename_bin, out_filename_bin + ".correct"), true))
      {
        remove( out_filename_bin.c_str());
      }

      // try storing a csv file
      model::RetrieveDatasetSubset::StoreMasterDataset( out_filename_csv, retrieve_dataset_subset);

      // check that the files are equal, if so, remove the output csv file
      if
      (
        BCL_ExampleCheck
        (
          io::File::FilesMatchWithinAbsoluteTolerance( out_filename_csv, out_filename_csv + ".correct", 0.00001)
          || io::File::FilesMatchWithinAbsoluteTolerance( out_filename_csv, out_filename_csv + ".win.correct", 0.00001),
          true
        )
      )
      {
        remove( out_filename_csv.c_str());
      }

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelRetrieveDatasetSubset

  const ExampleClass::EnumType ExampleModelRetrieveDatasetSubset::s_Instance
  (
    GetExamples().AddEnum( ExampleModelRetrieveDatasetSubset())
  );

} // namespace bcl
