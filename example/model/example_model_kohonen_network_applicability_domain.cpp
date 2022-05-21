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
#include "model/bcl_model_kohonen_network_applicability_domain.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_kohonen_network_applicability_domain.cpp
  //!
  //! @author mendenjl
  //! @date Sep 19, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelKohonenNetworkApplicabilityDomain :
    public ExampleInterface
  {
  public:

    ExampleModelKohonenNetworkApplicabilityDomain *Clone() const
    {
      return new ExampleModelKohonenNetworkApplicabilityDomain( *this);
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
      static const size_t s_number_values( 4);

      //initializing
      storage::VectorND< 2, linal::Vector< float> > features_results_setup[ s_number_values] =
      {
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 0.0), linal::MakeVector< float>( -1.0, 2.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 5.0), linal::MakeVector< float>( 2.0, -1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 1.0, 0.0), linal::MakeVector< float>( 2.0, -1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 1.0, 5.0), linal::MakeVector< float>( -1.0, 2.0))
      };

      // set up the data set
      storage::Vector< storage::VectorND< 2, linal::Vector< float> > > features_results;
      for( size_t count( 0); count < s_number_values; ++count)
      {
        features_results.PushBack( features_results_setup[ count]);
      }

      // create feature result data set
      util::ShPtr< descriptor::Dataset> features_results_frds
      (
        new descriptor::Dataset( features_results)
      );

      BCL_ExampleCheck( features_results_frds.IsDefined(), true);
      BCL_MessageStd( "Data set: " + util::Format()( features_results_frds));

      // rescale function for in an output
      util::ShPtr< model::RescaleFeatureDataSet> rescale_in
      (
        new model::RescaleFeatureDataSet( *features_results_frds->GetFeaturesPtr(), math::Range< float>( 0, 1))
      );

      BCL_MessageStd( "Rescale: " + util::Format()( rescale_in));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::KohonenNetworkApplicabilityDomain default_network;

      linal::Vector< double> map_dimensions( 2);
      map_dimensions( 0) = 20.0;
      map_dimensions( 1) = 50.0;

      linal::Vector< float> first_node_dimensions( 2, 0.0), last_node_dimensions( 2);
      last_node_dimensions( 0) = 19.0;
      last_node_dimensions( 1) = 49.0;

      // full constructor
      model::KohonenNetworkApplicabilityDomain network( model::KohonenNetworkAverage( map_dimensions, rescale_in));
      network.InitializeNodes( descriptor::Dataset());

      // test constructor
      BCL_ExampleIndirectCheck( network.GetCodeBook().GetSize(), 20 * 50, "Constructor");

    /////////////////
    // data access //
    /////////////////

      // test GetMapDimensions
      BCL_ExampleCheck( network.GetMapDimensions(), map_dimensions);

      // test position generation
      BCL_ExampleIndirectCheck
      (
        network.GetCodeBook().FirstElement().GetPosition(),
        first_node_dimensions,
        "position generation (first)"
      );
      BCL_ExampleIndirectCheck
      (
        network.GetCodeBook().LastElement().GetPosition(),
        last_node_dimensions,
        "position generation (last)"
      );

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write bcl object
      WriteBCLObject( network);

      // create default object
      model::KohonenNetworkApplicabilityDomain network_read;

      // read bcl object
      ReadBCLObject( network_read);

      // check read in object
      BCL_ExampleAssert( network_read.GetMapDimensions(), map_dimensions);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleClass::ExampleModelKohonenNetworkApplicabilityDomain

  const ExampleClass::EnumType ExampleModelKohonenNetworkApplicabilityDomain::s_Instance
  (
    GetExamples().AddEnum( ExampleModelKohonenNetworkApplicabilityDomain())
  );

} // namespace bcl
