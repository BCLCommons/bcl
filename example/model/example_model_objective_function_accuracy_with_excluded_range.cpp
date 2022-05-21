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
#include "model/bcl_model_objective_function_accuracy_with_excluded_range.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_template_instantiations.h"
#include "model/bcl_model_kappa_nearest_neighbor.h"
#include "model/bcl_model_objective_function_wrapper.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_objective_function_accuracy_with_excluded_range_monitoring_dataset.cpp
  //!
  //! @author mendenjl
  //! @date Jul 18, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelObjectiveFunctionAccuracyWithExcludedRange :
    public ExampleInterface
  {
  public:

    ExampleModelObjectiveFunctionAccuracyWithExcludedRange *Clone() const
    {
      return new ExampleModelObjectiveFunctionAccuracyWithExcludedRange( *this);
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

      // create reference data set
      util::ShPtr< storage::Vector< storage::VectorND< 2, linal::Vector< float> > > >
        features_results_training( new storage::Vector< storage::VectorND< 2, linal::Vector< float> > >());
      for( size_t count( 0); count < s_number_values; ++count)
      {
        features_results_training->PushBack( features_results_setup_training[ count]);
      }

      // create monitor data set
      util::ShPtr< storage::Vector< storage::VectorND< 2, linal::Vector< float> > > >
        features_results_monitor( new storage::Vector< storage::VectorND< 2, linal::Vector< float> > >());
      for( size_t count( 0); count < s_number_values; ++count)
      {
        features_results_monitor->PushBack( features_results_setup_query[ count]);
      }

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      model::ObjectiveFunctionAccuracyWithExcludedRange function_default;

      float activity_cutoff( 8);
      util::ShPtr< descriptor::Dataset> training_dataset( new descriptor::Dataset( *features_results_training));
      util::ShPtr< descriptor::Dataset> monitor_dataset( new descriptor::Dataset( *features_results_monitor));

      //! @brief constructor taking all necessary parameters
      model::ObjectiveFunctionWrapper function
      (
        monitor_dataset,
        util::Implementation< model::ObjectiveFunctionInterface>( model::ObjectiveFunctionAccuracyWithExcludedRange( activity_cutoff))
      );

      training_dataset->GetFeatures().Rescale( math::Range< float>( 0, 1));
      util::ShPtr< model::KappaNearestNeighbor> model( new model::KappaNearestNeighbor( training_dataset, 2));

    /////////////////
    // data access //
    /////////////////

      //GetData
      BCL_ExampleAssert( function.GetData()->GetSize(), size_t( 4));

      // set monitoring and rescale function. scaling function will be ignored
      function.SetData( monitor_dataset, training_dataset->GetFeaturesPtr()->GetScaling());

    ///////////////
    // operators //
    ///////////////

      //operator()
      BCL_ExampleAssert( function( util::ShPtr< model::Interface>( model)), float( 1.0));

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelObjectiveFunctionAccuracyWithExcludedRange

  const ExampleClass::EnumType ExampleModelObjectiveFunctionAccuracyWithExcludedRange::s_Instance
  (
    GetExamples().AddEnum( ExampleModelObjectiveFunctionAccuracyWithExcludedRange())
  );

} // namespace bcl
