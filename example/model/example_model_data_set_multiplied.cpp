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
#include "model/bcl_model_data_set_multiplied.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_data_set_multiplied.cpp
  //!
  //! @author geanesar
  //! @date Apr 28, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelDataSetMultiplied :
    public ExampleInterface
  {
  public:

    ExampleModelDataSetMultiplied *Clone() const
    {
      return new ExampleModelDataSetMultiplied( *this);
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
      descriptor::Dataset features_results( ( n + 1) * ( n + 1), size_t( 1), size_t( 1), size_t( 1));
      linal::MatrixReference< float> features( features_results.GetFeaturesReference());
      linal::MatrixReference< float> results( features_results.GetResultsReference());
      linal::MatrixReference< char> ids( features_results.GetIdsReference());

      descriptor::Dataset features_results_scaled( ( n + 1) * ( n + 1), size_t( 1), size_t( 1), size_t( 1));
      linal::MatrixReference< float> features_scaled( features_results_scaled.GetFeaturesReference());
      linal::MatrixReference< float> results_scaled( features_results_scaled.GetResultsReference());
      linal::MatrixReference< char> ids_scaled( features_results_scaled.GetIdsReference());

      size_t feature_number( 0);
      for( size_t i( 0); i <= n; ++i)
      {
        for( size_t j( 0); j <= n; ++j, ++feature_number)
        {
          features( feature_number, 0) = float( i) / float( n);
          results( feature_number, 0) = float( j + 1) / float( i + 1);
          ids( feature_number, 0) = 'O';
          features_scaled( feature_number, 0) = 2 * float( i) / float( n);
          results_scaled( feature_number, 0) = 1.5 * float( j + 1) / float( i + 1);
          ids_scaled( feature_number, 0) = 'M';
        }
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      model::DataSetMultiplied default_data_set_multiplied;

      // test constructor from the rmsd; the range is not strictly necessary here, but is given so that the object
      // that is output is consistent across 32 and 64 bit machines
      model::DataSetMultiplied data_set_multiplied( 2, 1.5, "M");

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      // create the multiplied_dataset

      BCL_ExampleIndirectCheck( data_set_multiplied( features_results)->GetSize(), features.GetNumberRows(), "operator()");

      // make sure the data set was actually rescaled
      BCL_ExampleCheckWithinAbsTolerance
      (
        data_set_multiplied( features_results)->GetFeatures().GetMatrix(),
        features_scaled,
        0.0001
      );

      BCL_ExampleCheckWithinAbsTolerance
      (
        data_set_multiplied( features_results)->GetResults().GetMatrix(),
        results_scaled,
        0.0001
      );

      BCL_ExampleCheck
      (
        data_set_multiplied( features_results)->GetIds().GetMatrix(),
        ids_scaled
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

  }; //end ExampleModelDataSetMultiplied

  const ExampleClass::EnumType ExampleModelDataSetMultiplied::s_Instance
  (
    GetExamples().AddEnum( ExampleModelDataSetMultiplied())
  );

} // namespace bcl
