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
#include "model/bcl_model_approximator_leverage_matrix.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "model/bcl_model_leverage_matrix.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_approximator_leverage_matrix.cpp
  //!
  //! @author mendenjl
  //! @date May 06, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelApproximatorLeverageMatrix :
    public ExampleInterface
  {
  public:

    ExampleModelApproximatorLeverageMatrix *Clone() const
    {
      return new ExampleModelApproximatorLeverageMatrix( *this);
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
      // create a training set
      storage::Vector< storage::VectorND< 2, linal::Vector< float> > > training_data_set;

      // create coefficients; the result values in the training data will be these coefficients times
      // a given feature
      linal::Vector< float> coefficients_a( 3), coefficients_b( 3);
      coefficients_a( 0) = 2.2;
      coefficients_a( 1) = -1.5;
      coefficients_a( 2) = 9.2;
      coefficients_b( 0) = -14.7;
      coefficients_b( 1) = -11.1;
      coefficients_b( 2) = 0.8;

      for( size_t counter( 0); counter < 10; ++counter)
      {
        linal::Vector< float> descriptors( 3);
        descriptors( 0) = sin( float( counter));
        descriptors( 1) = float( counter);
        descriptors( 2) = math::Sqrt( float( counter));
        training_data_set.PushBack
        (
          storage::VectorND< 2, linal::Vector< float> >
          (
            descriptors,
            linal::MakeVector< float>( descriptors * coefficients_a, descriptors * coefficients_b)
          )
        );
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

      //default constructor
      model::ApproximatorLeverageMatrix default_constructor;
      BCL_ExampleCheck( default_constructor.GetTrainingData().IsDefined(), false);

      //complete constructor
      util::Implementation< model::ApproximatorBase> iterate
      (
        "Leverage( objective function=RMSD)"
      );

      util::ShPtr< descriptor::Dataset> training_data_set_frds
      (
        new descriptor::Dataset( training_data_set)
      );

      iterate->SetTrainingData( training_data_set_frds);

      util::ShPtr< model::ObjectiveFunctionWrapper> objective( iterate->GetObjectiveFunction());

      objective->SetData( sp_monitor_frds, iterate->GetRescaleFeatureDataSet());

      BCL_ExampleCheck( iterate->GetTrainingData()->GetSize(), training_data_set_frds->GetSize());

      BCL_MessageVrb
      (
        "Initial rmsd: " + util::Format()( ( *objective)( iterate->GetCurrentApproximation()->First()))
      );
      BCL_MessageVrb
      (
        "Initial model: " + util::Format()( iterate->GetCurrentApproximation()->First())
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
      iterate->Approximate();
      result = iterate->GetCurrentApproximation();
      BCL_MessageVrb
      (
        "RMSD after training: = " + util::Format()( result->Second())
      );
      BCL_MessageDbg
      (
        "Model after training: =\n" + util::Format()( result->First())
      );

      // get a shptr to the actual iterate multiple linear regression
      util::ShPtr< model::LeverageMatrix> sp_model( result->First());

      // extract the coefficients for the first result value ( ~= coefficients_a) and the second ( ~= coefficients_b)
      float expected_lev_matrix[ 3][ 3] = {
          { 2.337, -2.277, 0.385},
          { -2.277, 49.21, -46.97},
          { 0.385, -46.97, 47.06}
      };
      BCL_ExampleIndirectCheckWithinTolerance
      (
        sp_model->GetWeights(),
        linal::Matrix< float>( size_t( 3), size_t( 3), &expected_lev_matrix[ 0][ 0]),
        0.1,
        "ApproximatorLeverageMatrix trains a proper leverage matrix"
      );

      // check that the rmsd is very small
      BCL_ExampleIndirectCheck
      (
        result->Second() <= float( 0.1),
        true,
        "the RMSD is too large, was " + util::Format()( result->Second())
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleClass::ExampleModelApproximatorLeverageMatrix

  const ExampleClass::EnumType ExampleModelApproximatorLeverageMatrix::s_Instance
  (
    GetExamples().AddEnum( ExampleModelApproximatorLeverageMatrix())
  );

} // namespace bcl
