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
#include "model/bcl_model_approximator_applicability_domain_kohonen.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_approximator_applicability_domain_kohonen.cpp
  //!
  //! @author mendenjl
  //! @date Sep 17, 2015
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelApproximatorApplicabilityDomainKohonen :
    public ExampleInterface
  {
  public:

    ExampleModelApproximatorApplicabilityDomainKohonen *Clone() const
    {
      return new ExampleModelApproximatorApplicabilityDomainKohonen( *this);
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
      //setup parts of iterate constructor
      const size_t iterations( 20);
      const float initial_radius( 4);

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

      // setup parts of ROCIterate
      util::ShPtr< model::ObjectiveFunctionWrapper> objective
      (
        new model::ObjectiveFunctionWrapper( sp_monitor_frds, util::ObjectDataLabel( "RMSD"))
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //default constructor
      model::ApproximatorApplicabilityDomainKohonen default_constructor;
      BCL_ExampleCheck( default_constructor.GetTrainingData().IsDefined(), false);

      linal::Vector< double> dimensions( 2);
      dimensions( 0) = 5.0;
      dimensions( 1) = 2.0;

      //complete constructor
      util::ShPtr< model::ApproximatorApplicabilityDomainKohonen> sp_iterate
      (
        new model::ApproximatorApplicabilityDomainKohonen
        (
          dimensions,
          iterations,
          initial_radius,
          objective,
          10,
          model::ApproximatorApplicabilityDomainKohonen::e_Gaussian
        )
      );

      util::ShPtr< descriptor::Dataset> training_data_set_frds
      (
        new descriptor::Dataset( training_data_set)
      );

      sp_iterate->SetTrainingData( training_data_set_frds);

      objective->SetData( sp_monitor_frds, sp_iterate->GetRescaleFeatureDataSet());

      BCL_ExampleCheck( sp_iterate->GetTrainingData()->GetSize(), training_data_set_frds->GetSize());

      BCL_MessageVrb
      (
        "Initial rmsd: " + util::Format()( ( *objective)( sp_iterate->GetCurrentModel()))
      );
      BCL_MessageVrb
      (
        "Initial model: " + util::Format()( sp_iterate->GetCurrentModel())
      );
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

      util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> > result;
      for( size_t counter( 0), iterations( 10); counter < iterations; ++counter)
      {
        sp_iterate->Next();
        result = sp_iterate->GetTracker().GetCurrent();
        BCL_MessageVrb
        (
          "RMSD at iteration " + util::Format()( counter) + " = " + util::Format()( result->Second())
        );
        BCL_MessageDbg
        (
          "Model at iteration " + util::Format()( counter) + " =\n" + util::Format()( result->First())
        );
      }

      BCL_ExampleIndirectCheck
      (
        result->Second() <= float( 1.0), true,
        "the RMSD is too large, was " + util::Format()( result->Second())
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleClass::ExampleModelApproximatorApplicabilityDomainKohonen

  const ExampleClass::EnumType ExampleModelApproximatorApplicabilityDomainKohonen::s_Instance
  (
    GetExamples().AddEnum( ExampleModelApproximatorApplicabilityDomainKohonen())
  );

} // namespace bcl
