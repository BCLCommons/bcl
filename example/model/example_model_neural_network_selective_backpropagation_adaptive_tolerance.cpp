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
#include "model/bcl_model_neural_network_selective_backpropagation_adaptive_tolerance.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_template_instantiations.h"
#include "model/bcl_model_objective_function_accuracy.h"
#include "model/bcl_model_objective_function_rmsd.h"
#include "model/bcl_model_retrieve_data_set_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_neural_network_selective_backpropagation_adaptive_tolerance.cpp
  //!
  //! @author mendenjl
  //! @date Feb 26, 2014
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelNeuralNetworkSelectiveBackpropagationAdaptiveTolerance :
    public ExampleInterface
  {
  public:

    ExampleModelNeuralNetworkSelectiveBackpropagationAdaptiveTolerance *Clone() const
    {
      return new ExampleModelNeuralNetworkSelectiveBackpropagationAdaptiveTolerance( *this);
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
      // retrieve a dataset
      util::Implementation< model::RetrieveDataSetBase>
        dataset_retriever
        (
          "File(filename=" + AddExampleInputPathToFilename( e_Model, "example_data_set_score.bcl") + ")"
        );
      descriptor::Dataset model_frds( *dataset_retriever->GenerateDataSet());
      model_frds.GetResults().Rescale( math::Range< float>( 0.0, 1.0));

      // construct the selective backpropagator
      model::NeuralNetworkSelectiveBackpropagationAdaptiveTolerance backprop_selector;

      // construct via an implementation, the more common method
      util::Implementation< model::NeuralNetworkSelectiveBackpropagationInterface>
        backprop_selector_up_to_1_off( "AdaptiveTolerant(initial tolerance=1, max tolerance=1)");

      // construct a second version that should backpropagate everything
      util::Implementation< model::NeuralNetworkSelectiveBackpropagationInterface>
        backprop_selector_all_impl( "AdaptiveTolerant(initial tolerance=0, max tolerance=1)");

      // initialize
      backprop_selector.Initialize( model_frds, model::ObjectiveFunctionAccuracy( 0.0));
      backprop_selector_up_to_1_off->Initialize( model_frds, model::ObjectiveFunctionAccuracy( 0.0));
      backprop_selector_all_impl->Initialize( model_frds, model::ObjectiveFunctionAccuracy( 0.0));

      // test alias
      BCL_ExampleCheck( backprop_selector.GetAlias(), "AdaptiveTolerant");

      // test that all features would be backpropagated and that the error vector is unchanged
      bool all_features_backpropagated( true);
      linal::Vector< float> error_vector( linal::MakeVector( float( 1)));
      for( size_t i( 0), nr_features( model_frds.GetSize()); i < nr_features; ++i)
      {
        all_features_backpropagated =
          backprop_selector_all_impl->ShouldBackpropagate
          (
            model_frds.GetResultsPtr()->GetMatrix().GetRow( 1),
            error_vector,
            i,
            0
          );
        error_vector = 1;
        if( !all_features_backpropagated)
        {
          // stop on error
          break;
        }
      }
      // call finalize
      backprop_selector_all_impl->FinalizeRound();

      BCL_ExampleIndirectCheck
      (
        all_features_backpropagated,
        true,
        "ShouldBackpropagate always returns true"
      );

      size_t number_backpropagated( 0);
      for( size_t i( 0), nr_features( model_frds.GetSize()); i < nr_features; ++i)
      {
        // test should backpropagate with a vector that is near the target
        error_vector = 0.05;
        number_backpropagated +=
          backprop_selector_up_to_1_off->ShouldBackpropagate
          (
            linal::Vector< float>( model_frds.GetResultsPtr()->GetMatrix().GetRow( i)) + float( 0.05),
            error_vector,
            i,
            0
          );
      }
      // call finalize
      backprop_selector_up_to_1_off->FinalizeRound();
      BCL_ExampleIndirectCheck
      (
        number_backpropagated,
        0,
        "ShouldBackpropagate with tolerance; all results should have been within the specified tolerance"
      );

      number_backpropagated = 0;
      for( size_t i( 0), nr_features( model_frds.GetSize()); i < nr_features; ++i)
      {
        error_vector = 0.25;
        number_backpropagated +=
          backprop_selector_up_to_1_off->ShouldBackpropagate
          (
            linal::Vector< float>( model_frds.GetResultsPtr()->GetMatrix().GetRow( i)),
            error_vector,
            i,
            0
          );
      }
      // call finalize
      backprop_selector_up_to_1_off->FinalizeRound();
      BCL_ExampleIndirectCheck
      (
        number_backpropagated,
        model_frds.GetSize(),
        "ShouldBackpropagate should notice that the error is outside the tolerance"
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelNeuralNetworkSelectiveBackpropagation

  const ExampleClass::EnumType ExampleModelNeuralNetworkSelectiveBackpropagationAdaptiveTolerance::s_Instance
  (
    GetExamples().AddEnum( ExampleModelNeuralNetworkSelectiveBackpropagationAdaptiveTolerance())
  );

} // namespace bcl
