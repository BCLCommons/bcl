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
#include "model/bcl_model_approximator_kappa_nearest_neighbor.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_criterion_number_iterations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_approximator_kappa_nearest_neighbor.cpp
  //! @brief this example tests the implementation of model::ApproximatorKappaNearestNeighbor which trains a kappa
  //! nearest neighbor
  //!
  //! @author mendenjl, fischea
  //! @date Jan 5, 2013
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelApproximatorKappaNearestNeighbor :
    public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleModelApproximatorKappaNearestNeighbor
    ExampleModelApproximatorKappaNearestNeighbor *Clone() const
    {
      return new ExampleModelApproximatorKappaNearestNeighbor( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {
      // create a training set
      storage::Vector< storage::VectorND< 2, linal::Vector< float> > > training_data_set;

      for( size_t counter( 0); counter < 10; ++counter)
      {
        linal::Vector< float> descriptors( 3);
        descriptors( 0) = random::GetGlobalRandom().Double();
        descriptors( 1) = float( counter);
        descriptors( 2) = math::Sqrt( float( counter));
        linal::Vector< float> target = linal::MakeVector< float>( float( counter));
        training_data_set.PushBack( storage::VectorND< 2, linal::Vector< float> >( descriptors, target));
      }

      // create a monitoring set
      storage::Vector< storage::VectorND< 2, linal::Vector< float> > > monitoring_data_set;
      for( size_t counter( 0); counter < 10; counter += 2)
      {
        monitoring_data_set.PushBack( training_data_set( counter));
      }

      util::ShPtr< descriptor::Dataset> sp_monitor_frds
      (
        new descriptor::Dataset( monitoring_data_set)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::ApproximatorKappaNearestNeighbor approximator_default;
      BCL_ExampleCheck( approximator_default.GetTrainingData().IsDefined(), false);

      // complete constructor
      util::Implementation< model::ApproximatorBase> approximator
      (
        "KappaNearestNeighbor( objective function=RMSD, min kappa=1, max kappa = 10)"
      );

      opti::CriterionNumberIterations< util::ShPtr< model::Interface>, float> criterion( 10);
      approximator->SetCriterion( criterion);

      util::ShPtr< descriptor::Dataset> training_data_set_frds
      (
        new descriptor::Dataset( training_data_set)
      );

      training_data_set_frds->GetFeatures().Rescale( math::Range< float>( 0, 1));
      approximator->SetTrainingData( training_data_set_frds);

      util::ShPtr< model::ObjectiveFunctionWrapper> objective( approximator->GetObjectiveFunction());

      objective->SetData( sp_monitor_frds, approximator->GetRescaleFeatureDataSet());

      BCL_ExampleCheck( approximator->GetTrainingData()->GetSize(), training_data_set_frds->GetSize());
      BCL_MessageVrb
      (
        "Training data: " + util::Format()( training_data_set_frds)
      );

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      // start the approximation and get the result
      approximator->Approximate();
      const float &result( approximator->GetTracker().GetBest()->Second());

      BCL_MessageVrb
      (
        "Final rmsd: " + util::Format()( ( *objective)( approximator->GetTracker().GetCurrent()->First()))
      );
      BCL_MessageVrb
      (
        "Final model: " + util::Format()( approximator->GetTracker().GetCurrent()->First())
      );

      BCL_ExampleIndirectCheck
      (
        result <= float( 1.0), true,
        "the RMSD is too large, was " + util::Format()( result)
      );

      return 0;

    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleModelApproximatorKappaNearestNeighbor

  const ExampleClass::EnumType ExampleModelApproximatorKappaNearestNeighbor::s_Instance
  (
     GetExamples().AddEnum( ExampleModelApproximatorKappaNearestNeighbor())
  );

} // namespace bcl
