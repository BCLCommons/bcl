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
#include "model/bcl_model_rescale_feature_data_set.h"

// includes from bcl - sorted alphabetically
#include "model/bcl_model_feature_data_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_rescale_feature_data_set.cpp
  //!
  //! @author loweew
  //! @date Nov 10, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelRescaleFeatureDataSet :
    public ExampleInterface
  {
  public:

    ExampleModelRescaleFeatureDataSet *Clone() const
    {
      return new ExampleModelRescaleFeatureDataSet( *this);
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

      model::FeatureDataSet< float> feat_ds( features);

      // setting range
      math::Range< float> range( 0.0, 1.0);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::RescaleFeatureDataSet rescale_fds_a;

      // construct taking range
      model::RescaleFeatureDataSet rescale_fds_b( feat_ds, range);

      // Clone()
      util::ShPtr< model::RescaleFeatureDataSet> rescale_clone( rescale_fds_b.Clone());

    /////////////////
    // data access //
    /////////////////

      // SetRange( RANGE)
      rescale_fds_a = model::RescaleFeatureDataSet( model::FeatureDataSet< float>(), range);

      // GetRange()
      BCL_ExampleCheck( rescale_fds_b.GetRange().GetMin(), range.GetMin());
      BCL_ExampleCheck( rescale_clone->GetRange().GetMax(), range.GetMax());
      BCL_ExampleCheck( rescale_fds_a.GetRange().GetMin(), range.GetMin());

      // GetRescaleFunction()
      util::ShPtr< model::FeatureDataSet< float> > rescaled_fds
      (
        new model::FeatureDataSet< float>( rescale_fds_b( feat_ds).GetMatrix())
      );
      BCL_ExampleCheck( rescale_fds_b.RescaleValue( 0, 1.0), float( 0));

    ///////////////
    // operators //
    ///////////////

      // check operator()( FEATURE_DATA_set)
      linal::Matrix< float> matrix( rescaled_fds->GetMatrix());
      BCL_ExampleCheck( matrix( 1, 1), 0.5);

    ////////////////
    // operations //
    ////////////////

      // DeScale( FEATURE)
      BCL_ExampleCheck( rescale_fds_b.DeScale( *rescaled_fds).GetMatrix()( 1, 2), float( 2));

    //////////////////////
    // input and output //
    //////////////////////

      // write
      WriteBCLObject( rescale_fds_b);

      model::RescaleFeatureDataSet rescale_read;
      //read
      ReadBCLObject( rescale_read);

      // comparing written and read in objects
      BCL_ExampleCheck( rescale_fds_b.GetRange().GetMin(), rescale_read.GetRange().GetMin());

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelRescaleFeatureDataSet

  const ExampleClass::EnumType ExampleModelRescaleFeatureDataSet::s_Instance
  (
    GetExamples().AddEnum( ExampleModelRescaleFeatureDataSet())
  );

} // namespace bcl
