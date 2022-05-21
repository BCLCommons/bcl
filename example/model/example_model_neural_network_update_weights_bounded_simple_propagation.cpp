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
#include "model/bcl_model_neural_network_update_weights_bounded_simple_propagation.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_neural_network_update_weights_bounded_simple_propagation.cpp
  //!
  //! @author mendenjl
  //! @date Jun 27, 2017
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelNeuralNetworkUpdateWeightsBoundedSimplePropagation :
    public ExampleInterface
  {
  public:

    ExampleModelNeuralNetworkUpdateWeightsBoundedSimplePropagation *Clone() const
    {
      return new ExampleModelNeuralNetworkUpdateWeightsBoundedSimplePropagation( *this);
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
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::NeuralNetworkUpdateWeightsBoundedSimplePropagation weight_updater_default;

      // parameters
      const float eta( 0.1);
      const float alpha( 0.5);
      const float mn( -2), mx( 2);

      // constructor from label
      // create a label to make a new simple prop with different parameters
      util::ObjectDataLabel weight_updater_label
      (
        weight_updater_default.GetAlias(),
        storage::Vector< util::ObjectDataLabel>::Create
        (
          util::ObjectDataLabel( "eta", util::Format()( eta)),
          util::ObjectDataLabel( "alpha", util::Format()( alpha)),
          util::ObjectDataLabel( "min", util::Format()( mn)),
          util::ObjectDataLabel( "max", util::Format()( mx))
        )
      );
      util::Implementation< model::NeuralNetworkUpdateWeightsInterface> weight_updater_impl( weight_updater_label);

      // check that the implementation is defined
      BCL_ExampleIndirectAssert( weight_updater_impl.IsDefined(), true, "constructor from label");

      // check that the implementation returns the same data label as it was given; this ensures that the object was
      // constructed properly
      BCL_ExampleIndirectCheck
      (
        weight_updater_impl.GetLabel(),
        weight_updater_label,
        "constructor from label"
      );

    ////////////////
    // operations //
    ////////////////

      // construct a mini-data-set to test the weight update
      const size_t number_weights( 2);

      // initialize the weight updater with the correct size
      weight_updater_impl->Initialize( number_weights);

      linal::Vector< float> weights( number_weights);
      linal::Vector< float> changes( number_weights);
      linal::Vector< float> slopes( number_weights);
      linal::Vector< float> zero_vector( number_weights, float( 0.0));

      // after each call to the operator(), reinitialize weights with the initial weights vector
      const float initial_weights_arr[ number_weights] = { 0.4, 0.96};
      const linal::Vector< float> initial_weights( number_weights, initial_weights_arr);

      // setup the weights with 0.4, 0.96.  Set changes and slopes to 0 to see that the weights don't change
      weights = initial_weights;
      changes = zero_vector;
      slopes = zero_vector;
      weight_updater_impl->SetChanges( changes);
      weight_updater_impl->operator ()( weights.Begin(), slopes.Begin());
      BCL_ExampleIndirectCheck
      (
        weights,
        initial_weights,
        "updating weights with 0 changes and slopes has no effect"
      );
      BCL_ExampleIndirectCheck
      (
        slopes,
        zero_vector,
        "updating weights with 0 changes and slopes leaves changes unchanged"
      );

      // after each call to the operator(), reinitialize weights with the initial weights vector
      const float initial_changes_arr[ number_weights] = { 0.5, -0.91};
      const linal::Vector< float> initial_changes( number_weights, initial_changes_arr);

      // reset the intial weights; Set changes to 0.5, -0.91, check that changes is updated accordingly
      weights = initial_weights;
      changes = initial_changes;
      slopes = zero_vector;
      weight_updater_impl->SetChanges( changes);
      weight_updater_impl->operator ()( weights.Begin(), slopes.Begin());
      BCL_ExampleIndirectCheck
      (
        weight_updater_impl->GetChanges(),
        initial_changes * alpha,
        "change update with non-0 changes, 0 slopes"
      );
      BCL_ExampleIndirectCheck
      (
        weights,
        initial_weights + weight_updater_impl->GetChanges(),
        "change update with non-0 changes, 0 slopes"
      );

      // set the slopes up now too
      const float initial_slopes_arr[ number_weights] = { 6.1, 5.692};
      const linal::Vector< float> initial_slopes( number_weights, initial_slopes_arr);

      weights = initial_weights;
      changes = initial_changes;
      slopes = initial_slopes;
      weight_updater_impl->SetChanges( changes);
      weight_updater_impl->operator ()( weights.Begin(), slopes.Begin());

      // check all values
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        linal::Distance( weight_updater_impl->GetChanges(), initial_changes * alpha + initial_slopes * eta),
        0,
        1.0e-6,
        "change update with non-0 changes, non-0 slopes"
      );
      BCL_ExampleIndirectCheck
      (
        weights,
        initial_weights + weight_updater_impl->GetChanges(),
        "weight update with non-0 changes, non-0 slopes"
      );
      BCL_ExampleIndirectCheck
      (
        slopes,
        zero_vector,
        "updating weights should always reset slopes to zero"
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelNeuralNetworkUpdateWeightsBoundedSimplePropagation

  const ExampleClass::EnumType ExampleModelNeuralNetworkUpdateWeightsBoundedSimplePropagation::s_Instance
  (
    GetExamples().AddEnum( ExampleModelNeuralNetworkUpdateWeightsBoundedSimplePropagation())
  );

} // namespace bcl
