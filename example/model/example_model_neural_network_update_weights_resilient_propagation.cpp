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
#include "model/bcl_model_neural_network_update_weights_resilient_propagation.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_resilient_propagation.cpp
  //!
  //! @author mendenjl
  //! @date May 12, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelNeuralNetworkUpdateWeightsResilientPropagation :
    public ExampleInterface
  {
  public:

    ExampleModelNeuralNetworkUpdateWeightsResilientPropagation *Clone() const
    {
      return new ExampleModelNeuralNetworkUpdateWeightsResilientPropagation( *this);
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
      // these are the constants from the resilient propagation algorithm
      const float increase_factor( 1.2);
      const float decrease_factor( 0.5);

      // default constructor
      model::NeuralNetworkUpdateWeightsResilientPropagation weight_updater_default;

      // constructor from label
      // create a label to make a new resil prop
      util::ObjectDataLabel weight_updater_label( weight_updater_default.GetAlias());
      util::Implementation< model::NeuralNetworkUpdateWeightsInterface> weight_updater_impl( weight_updater_label);

      // check that the implementation is defined
      BCL_ExampleIndirectAssert( weight_updater_impl.IsDefined(), true, "constructor from label");

      // check that the implementation returns the same data label as it was given; this ensures that the object was
      // constructed properly
      BCL_ExampleIndirectCheck
      (
        weight_updater_impl.GetLabel(),
        weight_updater_default.GetLabel(),
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
      linal::Vector< float> slopes( number_weights);
      linal::Vector< float> zero_vector( number_weights, float( 0.0));

      // after each call to the operator(), reinitialize weights with the initial weights vector
      const float initial_weights_arr[ number_weights] = { 0.4, 0.96};
      const linal::Vector< float> initial_weights( number_weights, initial_weights_arr);

      // setup the weights with 0.4, 0.96.  Set changes and slopes to 0 to see that the weights don't change
      weights = initial_weights;
      slopes = zero_vector;
      weight_updater_impl->SetChanges( zero_vector);
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

      // after each call to the operator(), reinitialize changes with the initial changes vector
      const float initial_changes_arr[ number_weights] = { 0.5, 0.91};
      const linal::Vector< float> initial_changes( number_weights, initial_changes_arr);

      // reset the intial weights; Set changes to 0.5, 0.91, check that changes is updated accordingly
      weights = initial_weights;
      slopes = zero_vector;
      // reset the weight updater by calling initialize again
      weight_updater_impl->Initialize( number_weights);
      weight_updater_impl->SetChanges( initial_changes);
      weight_updater_impl->operator ()( weights.Begin(), slopes.Begin());
      BCL_ExampleIndirectCheck
      (
        weight_updater_impl->GetChanges(),
        initial_changes,
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

      // also set up the negative slopes so that we can check how the sign of the slopes
      // influences the result
      const linal::Vector< float> neg_initial_slopes( float( -1.0) * initial_slopes);

      weights = initial_weights;
      slopes = initial_slopes;
      // reset the weight updater by calling initialize again
      weight_updater_impl->Initialize( number_weights);
      weight_updater_impl->SetChanges( initial_changes);
      weight_updater_impl->operator ()( weights.Begin(), slopes.Begin());

      // check all values
      BCL_ExampleIndirectCheck
      (
        weight_updater_impl->GetChanges(),
        initial_changes * increase_factor,
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
        initial_slopes,
        "slopes should only be changed when the signs of the last two slopes differ"
      );

      // try setting the slopes to negative values to see the effects of having opposite signs
      weights = initial_weights;
      slopes = neg_initial_slopes;
      // reset the weight updater by calling initialize again
      weight_updater_impl->Initialize( number_weights);
      weight_updater_impl->SetChanges( initial_changes);
      weight_updater_impl->operator ()( weights.Begin(), slopes.Begin());

      // check all values
      BCL_ExampleIndirectCheck
      (
        weight_updater_impl->GetChanges(),
        initial_changes * decrease_factor,
        "change update with non-0 changes, slope of opposite sign"
      );
      BCL_ExampleIndirectCheck
      (
        weights,
        initial_weights,
        "weights should not change when previous slopes are of opposite sign"
      );
      BCL_ExampleIndirectCheck
      (
        slopes,
        zero_vector,
        "slopes should be reset to 0 when the signs of the last two slopes differ"
      );

      // call the operator again to see that the resilient prop retains a memory of the previous slope

      // reset the initial slopes to -initial_slopes again, this way the update should recognize the trend of
      // decreasing weights
      weights = initial_weights;
      slopes = neg_initial_slopes;
      weight_updater_impl->SetChanges( initial_changes);
      weight_updater_impl->operator ()( weights.Begin(), slopes.Begin());
      // check all values
      BCL_ExampleIndirectCheck
      (
        weight_updater_impl->GetChanges(),
        initial_changes,
        "memory of slope of previous iteration"
      );
      BCL_ExampleIndirectCheck
      (
        weights,
        initial_weights - weight_updater_impl->GetChanges(),
        "memory of slope of previous iteration"
      );
      BCL_ExampleIndirectCheck
      (
        slopes,
        neg_initial_slopes,
        "slopes should only be changed when the signs of the last two slopes differ"
      );

      // reset the initial slopes to -initial_slopes again, this way the update should recognize the trend of
      // decreasing weights
      weights = initial_weights;
      slopes = neg_initial_slopes;
      weight_updater_impl->SetChanges( initial_changes);
      weight_updater_impl->operator ()( weights.Begin(), slopes.Begin());
      // check all values
      BCL_ExampleIndirectCheck
      (
        weight_updater_impl->GetChanges(),
        initial_changes * increase_factor,
        "memory of slope of previous iteration"
      );
      BCL_ExampleIndirectCheckWithinTolerance
      (
        weights,
        initial_weights - weight_updater_impl->GetChanges(),
        0.001,
        "memory of slope of previous iteration"
      );
      BCL_ExampleIndirectCheck
      (
        slopes,
        neg_initial_slopes,
        "slopes should only be changed when the signs of the last two slopes differ"
      );

      // this example does not check the delta-min/delta-max values because they are arbitrary
      // indeed, the value for delta-min in the bcl (1e-3) differs from that in the paper (1e-6)
      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelNeuralNetworkUpdateWeightsResilientPropagation

  const ExampleClass::EnumType ExampleModelNeuralNetworkUpdateWeightsResilientPropagation::s_Instance
  (
    GetExamples().AddEnum( ExampleModelNeuralNetworkUpdateWeightsResilientPropagation())
  );

} // namespace bcl
