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
#include "model/bcl_model_data_set_reduced_to_cluster_centers.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_data_set_reduced_to_cluster_centers.cpp
  //!
  //! @author mendenjl
  //! @date Dec 01, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelDataSetReducedToClusterCenters :
    public ExampleInterface
  {
  public:

    ExampleModelDataSetReducedToClusterCenters *Clone() const
    {
      return new ExampleModelDataSetReducedToClusterCenters( *this);
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
      // initialize the data set with n^2 points that run on equadistant 2D grid points (i,j), 0<=i<=1, 0<=j<=1
      const size_t n( 5);

      //initializing
      descriptor::Dataset features_results( ( n + 1) * ( n + 1), size_t( 2), size_t( 1), size_t( 0));
      linal::MatrixReference< float> features( features_results.GetFeaturesReference());
      linal::MatrixReference< float> results( features_results.GetResultsReference());

      descriptor::Dataset features_results_unscaled( ( n + 1) * ( n + 1), size_t( 2), size_t( 1), size_t( 0));
      linal::MatrixReference< float> features_unscaled( features_results_unscaled.GetFeaturesReference());
      linal::MatrixReference< float> results_unscaled( features_results_unscaled.GetResultsReference());

      size_t feature_number( 0);
      for( size_t i( 0); i <= n; ++i)
      {
        for( size_t j( 0); j <= n; ++j, ++feature_number)
        {
          features( feature_number, 0) = float( i) / float( n);
          features( feature_number, 1) = float( j) / float( n);
          results( feature_number, 0) = float( j + 1) / float( i + 1);
          features_unscaled( feature_number, 0) = float( i);
          features_unscaled( feature_number, 1) = float( j);
          results_unscaled( feature_number, 0) = float( j + 1) / float( i + 1);
        }
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      model::DataSetReducedToClusterCenters default_data_set_reduced_to_cluster_centers;

      // test constructor from the rmsd; the range is not strictly necessary here, but is given so that the object
      // that is output is consistent across 32 and 64 bit machines
      model::DataSetReducedToClusterCenters data_set_reduced_to_cluster_centers( 0.201, math::Range< size_t>( 1, 50));

      // test constructor from the rmsd and 3 clusters, giving up if the rmsd can't be found after 10 attempts
      model::DataSetReducedToClusterCenters data_set_reduced_to_three_cluster_centers
      (
        0.05,
        math::Range< size_t>( 3, 3),
        10
      );

      // test constructor from the rmsd and autoscale property
      model::DataSetReducedToClusterCenters data_set_reduced_to_any_cluster_centers_autoscaled
      (
        0.201,
        math::Range< size_t>( 0, features_results.GetSize()),
        10,
        true
      );

    /////////////////
    // data access //
    /////////////////

      // test get rmsd parameter
      BCL_ExampleCheck( default_data_set_reduced_to_cluster_centers.GetRMSDParameter(), 0.0);
      BCL_ExampleCheck( data_set_reduced_to_cluster_centers.GetRMSDParameter(), 0.201);

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      // create the reduced data set
      util::ShPtr< descriptor::Dataset>
        reduced_data_set( data_set_reduced_to_cluster_centers( features_results));

      // we should get 9 clusters, each containing the 4 nearest points
      BCL_ExampleIndirectCheck( reduced_data_set->GetSize(), 9, "operator()");

      // we won't find any clusters in the unscaled data set with the unscaled data, unless the autoscale parameter is
      // set
      util::ShPtr< descriptor::Dataset>
        reduced_data_set_unscaled( data_set_reduced_to_cluster_centers( features_results_unscaled));
      BCL_ExampleIndirectCheck
      (
        reduced_data_set_unscaled->GetSize(),
        features_results_unscaled.GetSize(),
        "operator() unscaled"
      );

      util::ShPtr< descriptor::Dataset>
        reduced_data_set_rescaled( data_set_reduced_to_any_cluster_centers_autoscaled( features_results_unscaled));
      BCL_ExampleIndirectCheck
      (
        reduced_data_set_rescaled->GetSize(),
        9,
        "operator() scaled"
      );

      // make sure the data set was actually rescaled
      BCL_ExampleIndirectCheck
      (
        reduced_data_set_rescaled->GetFeaturesReference()( 1, 1),
        2.0,
        "operator() rescales the output"
      );

    //////////////////////
    // input and output //
    //////////////////////

      model::DataSetReducedToClusterCenters data_set_reduced_to_cluster_centers_read;

      std::string filename( GetExampleOutputPathForBCLObject( model::DataSetReducedToClusterCenters(), ".bcl"));
      io::OFStream write;
      BCL_ExampleMustOpenOutputFile( write, filename);
      write << data_set_reduced_to_cluster_centers;
      io::File::CloseClearFStream( write);

      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, filename);
      read >> data_set_reduced_to_cluster_centers_read;
      io::File::CloseClearFStream( read);

      util::ShPtr< descriptor::Dataset>
        reduced_data_set_read( data_set_reduced_to_cluster_centers_read( features_results));

      // check read in object
      BCL_ExampleIndirectCheck( reduced_data_set_read->GetSize(), 9, "I/O");

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelDataSetReducedToClusterCenters

  const ExampleClass::EnumType ExampleModelDataSetReducedToClusterCenters::s_Instance
  (
    GetExamples().AddEnum( ExampleModelDataSetReducedToClusterCenters())
  );

} // namespace bcl
