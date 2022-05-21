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
#include "model/bcl_model_neural_network_selective_backpropagation_accuracy.h"

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
  //! @example example_model_neural_network_selective_backpropagation_accuracy.cpp
  //!
  //! @author mendenjl
  //! @date Jun 11, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelNeuralNetworkSelectiveBackpropagationAccuracy :
    public ExampleInterface
  {
  public:

    ExampleModelNeuralNetworkSelectiveBackpropagationAccuracy *Clone() const
    {
      return new ExampleModelNeuralNetworkSelectiveBackpropagationAccuracy( *this);
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
      model::NeuralNetworkSelectiveBackpropagationAccuracy backprop_selector;

      // construct via an implementation, the more common method
      util::Implementation< model::NeuralNetworkSelectiveBackpropagationInterface>
        backprop_selector_only_wrong_impl( "Accuracy(tolerance above=0.5,tolerance below=0.5)");

      // construct a second version that should backpropagate everything
      util::Implementation< model::NeuralNetworkSelectiveBackpropagationInterface>
        backprop_selector_all_impl( "Accuracy(pure classification = True,tolerance above=1,tolerance below=0)");

      // initialize
      backprop_selector.Initialize( model_frds, model::ObjectiveFunctionAccuracy( 0.0));
      backprop_selector_only_wrong_impl->Initialize( model_frds, model::ObjectiveFunctionAccuracy( 0.0));
      backprop_selector_all_impl->Initialize( model_frds, model::ObjectiveFunctionAccuracy( 0.0));

      // test alias
      BCL_ExampleCheck( backprop_selector.GetAlias(), "Accuracy");

      // test that all features would be backpropagated and that the error vector is unchanged
      bool all_features_backpropagated( true);
      linal::Vector< float> error_vector( linal::MakeVector( float( 1)));
      linal::Vector< float> too_high( linal::MakeVector( float( 0.5001)));
      linal::Vector< float> too_low( linal::MakeVector( float( 0.4999)));
      for( size_t i( 0), nr_features( model_frds.GetSize()); i < nr_features; ++i)
      {
        const float actual_value( model_frds.GetResultsPtr()->GetMatrix().GetRow( i)( 0));
        all_features_backpropagated =
          backprop_selector_all_impl->ShouldBackpropagate
          (
            actual_value > 0.0 ? too_low : too_high,
            error_vector,
            i,
            0
          );
        error_vector = float( 1.0);
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

      size_t number_backpropagated_all( 0), number_backpropagated_only_wrong( 0);
      for( size_t i( 0), nr_features( model_frds.GetSize()); i < nr_features; ++i)
      {
        // test should backpropagate with a vector just slightly closer to the cutoff than the real vector
        linal::Vector< float> test_vector
        (
          size_t( 1),
          float( 0.9 * ( model_frds.GetResultsPtr()->GetMatrix().GetRow( i)( 0) - float( 0.5)) + 0.5)
        );
        error_vector = float( 1.0);
        number_backpropagated_all += backprop_selector_all_impl->ShouldBackpropagate( test_vector, error_vector, i, 0);
        error_vector = float( 1.0);
        number_backpropagated_only_wrong
            += backprop_selector_only_wrong_impl->ShouldBackpropagate
               (
                 test_vector,
                 error_vector,
                 i,
                 0
               );
      }
      // call finalize
      backprop_selector_only_wrong_impl->FinalizeRound();
      backprop_selector_all_impl->FinalizeRound();
      BCL_ExampleIndirectCheck
      (
        number_backpropagated_only_wrong,
        0,
        "ShouldBackpropagate should not backpropagate anything if all features were on the right side of the cutoff if"
        " min fraction is set to 0"
      );
      BCL_ExampleIndirectCheck
      (
        number_backpropagated_all,
        72,
        "ShouldBackpropagate should backpropagate everything when all the predictions are closer to the cutoff than the"
        "real value and min fraction is at 1"
      );

      number_backpropagated_only_wrong = 0;
      for( size_t i( 0), nr_features( model_frds.GetSize()); i < nr_features; ++i)
      {
        number_backpropagated_only_wrong +=
          backprop_selector_only_wrong_impl->ShouldBackpropagate
          (
            model_frds.GetResultsPtr()->GetMatrix().GetRow( i),
            error_vector,
            i,
            0
          );
      }
      // call finalize
      backprop_selector_only_wrong_impl->FinalizeRound();
      BCL_ExampleIndirectCheck
      (
        number_backpropagated_only_wrong,
        0,
        "ShouldBackpropagate should not have returned true because all features were correct"
      );

      number_backpropagated_only_wrong = 0;
      for( size_t i( 0), nr_features( model_frds.GetSize()); i < nr_features; ++i)
      {
        error_vector = float( 1.0);
        number_backpropagated_only_wrong +=
          backprop_selector_only_wrong_impl->ShouldBackpropagate
          (
            linal::MakeVector< float>( 0.1),
            error_vector,
            i,
            0
          );
      }
      BCL_ExampleIndirectCheck
      (
        number_backpropagated_only_wrong,
        36,
        "ShouldBackpropagate should have returned true for almost all the dataset since the results were always wrong"
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelNeuralNetworkSelectiveBackpropagation

  const ExampleClass::EnumType ExampleModelNeuralNetworkSelectiveBackpropagationAccuracy::s_Instance
  (
    GetExamples().AddEnum( ExampleModelNeuralNetworkSelectiveBackpropagationAccuracy())
  );

} // namespace bcl
