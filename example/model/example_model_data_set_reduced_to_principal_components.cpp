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
#include "model/bcl_model_data_set_reduced_to_principal_components.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_data_set_reduced_to_principal_components.cpp
  //!
  //! @author loweew
  //! @date Sept 27, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelDataSetReducedToPrincipalComponents :
    public ExampleInterface
  {
  public:

    ExampleModelDataSetReducedToPrincipalComponents *Clone() const
    {
      return new ExampleModelDataSetReducedToPrincipalComponents( *this);
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
      const std::string data_set_filename( AddExampleInputPathToFilename( e_Model, "example_data_set.bcl"));
      const std::string pca_filename( AddExampleInputPathToFilename( e_Model, "pca.bcl"));

      // make a label that indicates what the data set holds for features and results
      util::ObjectDataLabel feature_label( "Combine(NRotBond,LogS,Dipole,Polariz,TPSA)");
      util::ObjectDataLabel results_label( "Combine(HAcc,HDon)");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      model::DataSetReducedToPrincipalComponents default_constructor;

      // test constructor from fraction of total sum to keep
      model::DataSetReducedToPrincipalComponents data_set_reduced_to_principal_components( 0.75, 5);

    /////////////////
    // data access //
    /////////////////

      // test get alias
      util::Implementation< model::RetrieveDataSetBase> impl_data_set_retriever
      (
        "File(filename=" + data_set_filename + ", number chunks=3, chunks=\"[0,2]\")"
      );

      data_set_reduced_to_principal_components.SetFilename( pca_filename);
      data_set_reduced_to_principal_components.SetRetriever( impl_data_set_retriever);
      BCL_ExampleCheck( data_set_reduced_to_principal_components.GetAlias(), "PCA");
      data_set_reduced_to_principal_components.SelectFeatures( feature_label);
      data_set_reduced_to_principal_components.SelectResults( results_label);

      BCL_ExampleCheck( data_set_reduced_to_principal_components.GetFeatureCode().ToString(), feature_label.ToString());
      BCL_ExampleCheck( data_set_reduced_to_principal_components.GetResultCode().ToString(), results_label.ToString());

    ////////////////
    // operations //
    ////////////////

      util::ShPtr< descriptor::Dataset> reduced_data_set( data_set_reduced_to_principal_components.GenerateDataSet());
      float expected_result_data[ 9] = { 0.0, 0.0, 0.0, 16.377, 9.31205, 11.3796, 32, 18, 23};
      linal::Matrix< float> expected_result_matrix( 3, 3, expected_result_data);

      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        reduced_data_set->GetFeaturesReference(),
        expected_result_matrix,
        0.001,
        "GenerateDataSet(), actual results: " + util::Format()( reduced_data_set)
      );

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      model::DataSetReducedToPrincipalComponents data_set_reduced_to_principal_components_read;

      std::string filename( GetExampleOutputPathForBCLObject( model::DataSetReducedToPrincipalComponents(), ".bcl"));
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, filename);
      write << data_set_reduced_to_principal_components;
      io::File::CloseClearFStream( write);

      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, filename);
      read >> data_set_reduced_to_principal_components_read;
      io::File::CloseClearFStream( read);

      util::ShPtr< descriptor::Dataset>
        reduced_data_set_read( data_set_reduced_to_principal_components_read.GenerateDataSet());

      // check read in object
      BCL_ExampleIndirectCheck( reduced_data_set_read->GetSize(), 3, "I/O");

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelDataSetReducedToPrincipalComponents

  const ExampleClass::EnumType ExampleModelDataSetReducedToPrincipalComponents::s_Instance
  (
    GetExamples().AddEnum( ExampleModelDataSetReducedToPrincipalComponents())
  );

} // namespace bcl
