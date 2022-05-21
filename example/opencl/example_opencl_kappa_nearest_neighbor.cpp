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
#include "opencl/bcl_opencl_kappa_nearest_neighbor.h"

// includes from bcl - sorted alphabetically
#include "model/bcl_model_kappa_nearest_neighbor.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opencl_kappa_nearest_neighbor.cpp
  //!
  //! @author loweew
  //! @date Apr 13, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenclKappaNearestNeighbor :
    public ExampleInterface
  {
  public:

    ExampleOpenclKappaNearestNeighbor *Clone() const
    {
      return new ExampleOpenclKappaNearestNeighbor( *this);
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
      // do not try to run opencl commands if no queue was found
      if( !opencl::GetTools().HasCommandQueues())
      {
        return 1;
      }

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
      training->GetResults().Rescale( math::Range< float>( 0, 1));

      // kappa value
      size_t kappa( 1);

      opencl::CommandQueue queue( opencl::GetTools().GetFirstCommandQueue());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      opencl::KappaNearestNeighbor ocl_knn_default( queue);

      // constructor from properties
      util::ShPtr< model::KappaNearestNeighbor> knn( new model::KappaNearestNeighbor( training, kappa));

      // constructor from properties
      util::ShPtr< opencl::KappaNearestNeighbor> ocl_knn( new opencl::KappaNearestNeighbor( training, kappa, queue));

      // clone
      const util::ShPtr< opencl::KappaNearestNeighbor> knn_clone( ocl_knn_default.Clone());

    /////////////////
    // data access //
    /////////////////

      // class identifiers
      BCL_MessageStd( "class name: " + knn_clone->GetClassIdentifier());
      BCL_Example_Check
      (
        ( GetStaticClassName< opencl::KappaNearestNeighbor>() == "bcl::opencl::KappaNearestNeighbor")
        && ( knn_clone->GetClassIdentifier() == GetStaticClassName< opencl::KappaNearestNeighbor>()),
        "incorrect class name: static class name: " + ( GetStaticClassName< opencl::KappaNearestNeighbor>())
        + " or class identifier: " + knn_clone->GetClassIdentifier()
      );

      // test Set/GetKappa()
      ocl_knn_default.SetKappa( 2);
      BCL_ExampleIndirectCheck( ocl_knn_default.GetKappa(), 2, "ocl_knn_default.SetKappa( 2)");

    ///////////////
    // operators //
    ///////////////

      BCL_MessageStd( "start timing....");
      util::Stopwatch timer;

      // running knn algorithm and return list of results (predicted and experimental)
      util::ShPtr< model::FeatureDataSet< float> > predictions
      (
        new model::FeatureDataSet< float>( ocl_knn->operator()( *query->GetFeaturesPtr()))
      );

      BCL_MessageStd( "Time to predict using gpu: " + util::Format()( timer.GetProcessDuration().GetTimeAsHourMinuteSecondMilliSeconds()));

      BCL_MessageStd( "ocl knn predictions: " + util::Format()( *predictions));
      BCL_MessageStd( "cpu knn predictions: " + util::Format()( knn->operator()( *query->GetFeaturesPtr())));

    //////////////////////
    // checking results //
    //////////////////////

      // check predicted result
      BCL_ExampleCheckWithinTolerance( double( *predictions->operator[]( 0)), 1.0, 0.001);

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOpenclKappaNearestNeighbor

  const ExampleClass::EnumType ExampleOpenclKappaNearestNeighbor::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenclKappaNearestNeighbor())
  );

} // namespace bcl
