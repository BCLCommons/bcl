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
#include "model/bcl_model_feature_data_reference.h"

// includes from bcl - sorted alphabetically
#include "model/bcl_model_feature_data_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_feature_data_reference.cpp
  //!
  //! @author mendenjl
  //! @date Apr 06, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelFeatureDataReference :
    public ExampleInterface
  {
  public:

    ExampleModelFeatureDataReference *Clone() const
    {
      return new ExampleModelFeatureDataReference( *this);
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
      // data set construction
      linal::Matrix< float> features( 3, 3), results( 3, 1);
      features.ReplaceRow( 0, linal::Vector< float>( size_t( 3), float( 1)));
      features.ReplaceRow( 1, linal::Vector< float>( size_t( 3), float( 2)));
      features.ReplaceRow( 2, linal::Vector< float>( size_t( 3), float( 3)));

      results.ReplaceRow( 0, linal::Vector< float>( size_t( 1), float( 3)));
      results.ReplaceRow( 1, linal::Vector< float>( size_t( 1), float( 2)));
      results.ReplaceRow( 2, linal::Vector< float>( size_t( 1), float( 1)));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // constructing feature data set of features from data set interface
      model::FeatureDataSet< float> feature_ds_own( features);

      // constructing feature data set of results from data set interface
      model::FeatureDataSet< float> result_ds_own( results);

      // constructing feature data set of features from data set interface
      model::FeatureDataReference< float> feature_ds( feature_ds_own.GetMatrix());

      // constructing feature data set of results from data set interface
      model::FeatureDataReference< float> result_ds( result_ds_own.GetMatrix());

    /////////////////
    // data access //
    /////////////////

      // GetFeatureSize()
      BCL_ExampleCheck( feature_ds.GetFeatureSize(), features.GetNumberCols());

      // GetNumberFeatures() of features FeatureDataSet
      BCL_ExampleCheck( feature_ds.GetNumberFeatures(), features.GetNumberRows());

      // GetFeatureSize() of results FeatureDataSet
      BCL_ExampleCheck( result_ds.GetFeatureSize(), results.GetNumberCols());

      // GetNumberFeatures() of results FeatureDataSet
      BCL_ExampleCheck( result_ds.GetNumberFeatures(), results.GetNumberRows());

      // GetMatrix() of features FeatureDataSet
      BCL_ExampleCheck( feature_ds.GetMatrix(), features);

      // GetMatrix() of results FeatureDataSet
      BCL_ExampleCheck( result_ds.GetMatrix(), results);

      // GetMatrix( POS, LENGTH) of features FeatureDataSet
      BCL_ExampleCheck( feature_ds.GetMatrix( 1, 2)( 0, 1), features( 1, 1));
      BCL_ExampleCheck( feature_ds.GetMatrix( 1, 2).GetNumberCols(), features.GetNumberCols());

      // GetMatrix( POS, LENGTH) of features FeatureDataSet
      BCL_ExampleCheck( result_ds.GetMatrix( 1, 2)( 1, 0), results( 2, 0));
      BCL_ExampleCheck( result_ds.GetMatrix( 1, 2).GetNumberCols(), results.GetNumberCols());

    ////////////////
    // operations //
    ////////////////

      // operator
      const math::Range< float> input_range( 0, 1);
      BCL_ExampleCheck( feature_ds.GetMatrix()( 1, 1), float( 2));

    ///////////////
    // operators //
    ///////////////

      // operator[]()
      BCL_ExampleCheck( *feature_ds[ 1], features( 1, 0));

      // operator()()
      BCL_ExampleCheck( feature_ds( 1).GetSize(), features.GetNumberCols());

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelFeatureDataReference

  const ExampleClass::EnumType ExampleModelFeatureDataReference::s_Instance
  (
    GetExamples().AddEnum( ExampleModelFeatureDataReference())
  );

} // namespace bcl
