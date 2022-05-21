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
#include "model/bcl_model_neural_network_perturb_attenuate.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_resilient_propagation.cpp
  //!
  //! @author mendenjl
  //! @date Apr 30, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelNeuralNetworkPerturbAttenuate :
    public ExampleInterface
  {
  public:

    ExampleModelNeuralNetworkPerturbAttenuate *Clone() const
    {
      return new ExampleModelNeuralNetworkPerturbAttenuate( *this);
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
      // default constructor
      model::NeuralNetworkPerturbAttenuate weight_updater_default;

      // constructor from label
      // create a label to make a weight attenuator which will attenuate weights according to 0.0001 + 0.01 x^2, where
      // x is the weight
      util::ObjectDataLabel attenuate_label( weight_updater_default.GetAlias() + "(0.0001,0.1,0.01)");
      util::Implementation< model::NeuralNetworkPerturbationInterface> attenuator_impl( attenuate_label);

      // check that the implementation is defined
      BCL_ExampleIndirectAssert( attenuator_impl.IsDefined(), true, "constructor from label");

      // check that the implementation returns the same data label as it was given; this ensures that the object was
      // constructed properly
      BCL_ExampleIndirectCheck
      (
        attenuator_impl.GetLabel().GetArgument( 1),
        attenuate_label.GetArgument( 1),
        "constructor from label"
      );

    ////////////////
    // operations //
    ////////////////

      const size_t number_weights( 6);
      const float initial_weights_arr[ number_weights] = { 0.0, 1.0, 2.0, -0.5, -0.0000001, -14.0};
      linal::Matrix< float> weights( size_t( 3), size_t( 2), initial_weights_arr);
      linal::Vector< float> polynomial_coefficients( linal::MakeVector< float>( 0.0001, 0.1, 0.01));
      math::Polynomial poly( polynomial_coefficients.Begin(), polynomial_coefficients.End());
      const float expected_weights_arr[ number_weights] =
      {
        float( initial_weights_arr[ 0]),
        float( initial_weights_arr[ 1] - poly( initial_weights_arr[ 1])),
        float( initial_weights_arr[ 2] - poly( initial_weights_arr[ 2])),
        float( initial_weights_arr[ 3] + poly( math::Absolute( initial_weights_arr[ 3]))),
        float( 0.0),
        float( initial_weights_arr[ 5] + poly( math::Absolute( initial_weights_arr[ 5])))
      };
      linal::Matrix< float> expected_weights( size_t( 3), size_t( 2), expected_weights_arr);

      // setup the weights with 0.4, 0.96.  Set changes and slopes to 0 to see that the weights don't change
      attenuator_impl->operator ()( weights);
      BCL_ExampleIndirectCheck( weights, expected_weights, "operator()");

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelNeuralNetworkPerturbAttenuate

  const ExampleClass::EnumType ExampleModelNeuralNetworkPerturbAttenuate::s_Instance
  (
    GetExamples().AddEnum( ExampleModelNeuralNetworkPerturbAttenuate())
  );

} // namespace bcl
