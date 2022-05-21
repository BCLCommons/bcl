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
#include "model/bcl_model_kappa_nearest_neighbor.h"

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_kappa_nearest_neighbor.cpp
  //!
  //! @author loweew
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelKappaNearestNeighbor :
    public ExampleInterface
  {
  public:

    ExampleModelKappaNearestNeighbor *Clone() const
    {
      return new ExampleModelKappaNearestNeighbor( *this);
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
    ////////////////////
    // Preparing data //
    ////////////////////

      static const size_t s_number_values( 4);

      //initializing training data
      storage::VectorND< 2, linal::Vector< float> > features_results_setup_training[ s_number_values] =
      {
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 1.0), linal::MakeVector< float>( 1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 10.0), linal::MakeVector< float>( 0.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 2.0), linal::MakeVector< float>( 1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 11.0), linal::MakeVector< float>( 0.0))
      };

      //initializing query data
      storage::VectorND< 2, linal::Vector< float> > features_results_setup_query[ s_number_values] =
      {
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 2.45), linal::MakeVector< float>( 1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 11.85), linal::MakeVector< float>( 0.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 0.896), linal::MakeVector< float>( 1.0)),
        storage::VectorND< 2, linal::Vector< float> >
          ( linal::MakeVector< float>( 0.0, 8.954), linal::MakeVector< float>( 0.0))
      };

      // create training data set
      util::ShPtr< storage::Vector< storage::VectorND< 2, linal::Vector< float> > > >
        features_results_training( new storage::Vector< storage::VectorND< 2, linal::Vector< float> > >());
      for( size_t count( 0); count < s_number_values; ++count)
      {
        features_results_training->PushBack( features_results_setup_training[ count]);
      }

      // create query data set
      util::ShPtr< storage::Vector< storage::VectorND< 2, linal::Vector< float> > > >
        features_results_query( new storage::Vector< storage::VectorND< 2, linal::Vector< float> > >());
      for( size_t count( 0); count < s_number_values; ++count)
      {
        features_results_query->PushBack( features_results_setup_query[ count]);
      }

      // construct feature result data sets
      util::ShPtr< descriptor::Dataset > training( new descriptor::Dataset( *features_results_training));
      util::ShPtr< descriptor::Dataset > query   ( new descriptor::Dataset( *features_results_query   ));
      training->GetFeatures().Rescale( math::Range< float>( 0, 1));

      // print the scaled dataset
      BCL_MessageStd( "rescaled training features:\n" + util::Format()( training));

      // kappa value
      size_t kappa( 1);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::KappaNearestNeighbor knn_default;

      // constructor from properties
      util::ShPtr< model::KappaNearestNeighbor> knn( new model::KappaNearestNeighbor( training, kappa));

      // clone
      const util::ShPtr< model::KappaNearestNeighbor> knn_clone( knn_default.Clone());

    /////////////////
    // data access //
    /////////////////

      // test Set/GetKappa()
      knn_default.SetKappa( 2);
      BCL_ExampleCheck( knn_default.GetKappa(), 2);

    ///////////////
    // operators //
    ///////////////

      // running knn algorithm and return list of results (predicted and experimental)
      util::ShPtr< model::FeatureDataSet< float> > predictions
      (
        new model::FeatureDataSet< float>( knn->operator()( *query->GetFeaturesPtr()))
      );

    //////////////////////
    // checking results //
    //////////////////////

      // check predicted result
      BCL_ExampleIndirectCheck( *predictions->operator[]( 0), 1.0, "The predicted output of the query point");

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelKappaNearestNeighbor

  const ExampleClass::EnumType ExampleModelKappaNearestNeighbor::s_Instance
  (
    GetExamples().AddEnum( ExampleModelKappaNearestNeighbor())
  );

} // namespace bcl
