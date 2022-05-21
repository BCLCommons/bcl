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
#include "model/bcl_model_data_set_reduced_to_k_means.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_data_set_reduced_to_k_means.cpp
  //!
  //! @author mendenjl
  //! @date Dec 03, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelDataSetReducedToKMeans :
    public ExampleInterface
  {
  public:

    ExampleModelDataSetReducedToKMeans *Clone() const
    {
      return new ExampleModelDataSetReducedToKMeans( *this);
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
      const size_t n( 5);

      // initialize the data set with n^2 points that run on equadistant 2D grid points (i,j), 0<=i<=1, 0<=j<=1
      descriptor::Dataset features_results( ( n + 1) * ( n + 1), size_t( 2), size_t( 1), size_t( 0));
      linal::MatrixReference< float> features( features_results.GetFeaturesReference());
      linal::MatrixReference< float> results( features_results.GetResultsReference());

      descriptor::Dataset features_results_unscaled( ( n + 1) * ( n + 1), size_t( 2), size_t( 1), size_t( 0));
      linal::MatrixReference< float> features_unscaled( features_results_unscaled.GetFeaturesReference());
      linal::MatrixReference< float> results_unscaled( features_results_unscaled.GetResultsReference());

      size_t row_number( 0);
      for( size_t i( 0); i <= n; ++i)
      {
        for( size_t j( 0); j <= n; ++j, ++row_number)
        {
          // create the result vector
          linal::Vector< float> result_vector( 1, float( j + 1) / float( i + 1));

          // create the feature vector and set its values
          linal::Vector< float> feature_vector( 2);
          feature_vector( 0) = i;
          feature_vector( 1) = j;
          features_unscaled.ReplaceRow( row_number, feature_vector);
          results_unscaled.ReplaceRow( row_number, result_vector);

          // scale the feature vector
          feature_vector /= float( n);
          features.ReplaceRow( row_number, feature_vector);

          // set the result vectors
          results.ReplaceRow( row_number, result_vector);
        }
      }

      util::ShPtr< descriptor::Dataset> sp_features_results
      (
        util::CloneToShPtr( features_results)
      );
      util::ShPtr< descriptor::Dataset> sp_features_results_unscaled
      (
        util::CloneToShPtr( features_results_unscaled)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      model::DataSetReducedToKMeans default_data_set_reduced_to_k_means;

      // test constructor from # clusters
      model::DataSetReducedToKMeans data_set_reduced_to_k_means( 15);

      // test constructor from # clusters and autoscale property
      model::DataSetReducedToKMeans data_set_reduced_to_k_means_autoscaled( 15, true);

    /////////////////
    // data access //
    /////////////////

      // test get number clusters
      BCL_ExampleCheck( default_data_set_reduced_to_k_means.GetNumberClusters(), 0);
      BCL_ExampleCheck( data_set_reduced_to_k_means.GetNumberClusters(), 15);

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      // create the reduced data set
      util::ShPtr< descriptor::Dataset>
        reduced_data_set( data_set_reduced_to_k_means( sp_features_results));

      // we should get 15 clusters as requested
      BCL_ExampleIndirectCheck( reduced_data_set->GetSize(), 15, "operator()");

      util::ShPtr< descriptor::Dataset>
        reduced_data_set_rescaled( data_set_reduced_to_k_means_autoscaled( sp_features_results_unscaled));

      // make sure the data set was actually rescaled
      const float expected_second_cluster[ 2] = { 0.0, 1.0};
      BCL_ExampleIndirectCheck
      (
        linal::Vector< float>( 2, reduced_data_set_rescaled->GetFeaturesReference()[ 1]),
        linal::Vector< float>( 2, expected_second_cluster),
        "operator() rescales the output"
      );

    //////////////////////
    // input and output //
    //////////////////////

      util::Implementation< model::RetrieveDataSetBase>
        dataset_reduced_to_k_means_impl( data_set_reduced_to_k_means),
        dataset_reduced_to_k_means_impl_read;

      std::string filename( GetExampleOutputPathForBCLObject( model::DataSetReducedToKMeans(), ".bcl"));
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, filename);
      write << dataset_reduced_to_k_means_impl;
      io::File::CloseClearFStream( write);

      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, filename);
      read >> dataset_reduced_to_k_means_impl;
      io::File::CloseClearFStream( read);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelDataSetReducedToKMeans

  const ExampleClass::EnumType ExampleModelDataSetReducedToKMeans::s_Instance
  (
    GetExamples().AddEnum( ExampleModelDataSetReducedToKMeans())
  );

} // namespace bcl
