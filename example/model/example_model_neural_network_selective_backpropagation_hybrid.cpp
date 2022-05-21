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
#include "model/bcl_model_neural_network_selective_backpropagation_hybrid.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
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
  //! @example example_model_neural_network_selective_backpropagation_hybrid.cpp
  //!
  //! @author mendenjl
  //! @date Aug 05, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelNeuralNetworkSelectiveBackpropagationHybrid :
    public ExampleInterface
  {
  public:

    ExampleModelNeuralNetworkSelectiveBackpropagationHybrid *Clone() const
    {
      return new ExampleModelNeuralNetworkSelectiveBackpropagationHybrid( *this);
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
        dataset_retriever( "File(filename=" + AddExampleInputPathToFilename( e_Model, "example_data_set_score.bcl") + ")");
      descriptor::Dataset model_frds( *dataset_retriever->GenerateDataSet());
      model_frds.GetResults().Rescale( math::Range< float>( 0.0, 1.0));

      // construct the selective backpropagator
      model::NeuralNetworkSelectiveBackpropagationHybrid backprop_selector;

      // construct via an implementation, the more common method
      util::Implementation< model::NeuralNetworkSelectiveBackpropagationInterface>
        balanced( "Hybrid(stability=0)");

      const float cutoff( 1.1); // only 15 of the 72 points are above this cutoff
      const float scaled_cutoff( model_frds.GetResultsPtr()->GetScaling()->RescaleValue( 0, cutoff));

      // initialize
      backprop_selector.Initialize( model_frds, model::ObjectiveFunctionAccuracy( cutoff));
      balanced->Initialize( model_frds, model::ObjectiveFunctionAccuracy( cutoff));

      // test alias
      BCL_ExampleCheck( backprop_selector.GetAlias(), "Hybrid");

      // get a vector of all the results, and shuffle them
      const size_t nr_features( model_frds.GetSize());
      storage::Vector< float> shuffled_results( nr_features, model_frds.GetResultsPtr()->GetMatrix().Begin());
      shuffled_results.Shuffle();

      // test that roughly equal number of features are backpropagated above and below the cutoff, given roughly equal #s
      // predicted above and below the cutoff
      // This is a realistic assumption
      int backprop_count( 0), backprop_count_above( 0), backprop_count_below( 0);
      float epsilon( float( 1.0) + std::numeric_limits< float>::epsilon());
      linal::Vector< float> prediction_high( 1, float( scaled_cutoff * ( float( 1.0) + epsilon)));
      linal::Vector< float> prediction_low( 1, float( scaled_cutoff / ( float( 1.0) + epsilon)));
      linal::MatrixConstReference< float> results( model_frds.GetResultsPtr()->GetMatrix());
      for( size_t iteration( 0); iteration < size_t( 1); ++iteration)
      {
        backprop_count = backprop_count_above = backprop_count_below = 0;
        for( size_t i( 0); i < nr_features; ++i)
        {
          linal::Vector< float> row( size_t( 1), shuffled_results( i));
          linal::Vector< float> error( float( -1.0) * results.GetRow( i));
          if( balanced->ShouldBackpropagate( row, error, i, 0))
          {
            ++backprop_count;
            if( results.GetRow( i)( 0) < scaled_cutoff)
            {
              ++backprop_count_below;
            }
            else
            {
              ++backprop_count_above;
            }
          }
        }
        // call finalize
        balanced->FinalizeRound();
      }

      // test that those above and below the cutoff were roughly equal
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        backprop_count_below,
        1,
        2,
        "Imbalance (due to better prediction of below-cutoff values)"
      );
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        backprop_count_below,
        1,
        2,
        "Imbalance (due to better prediction of below-cutoff values)"
      );

      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        backprop_count,
        16,
        4,
        "Check of absolute # backpropagated"
      );

      // now test more trivial cases: test that vanishingly little is backpropagated if the predictions are right
      for( size_t iteration( 0); iteration < size_t( 10); ++iteration)
      {
        backprop_count = backprop_count_above = backprop_count_below = 0;
        for( size_t i( 0); i < nr_features; ++i)
        {
          linal::Vector< float> row( results.GetRow( i));
          linal::Vector< float> error( float( -1.0) * results.GetRow( i));
          if( balanced->ShouldBackpropagate( row, error, i, 0))
          {
            ++backprop_count;
            if( results.GetRow( i)( 0) < scaled_cutoff)
            {
              ++backprop_count_below;
            }
            else
            {
              ++backprop_count_above;
            }
          }
        }
        // call finalize
        balanced->FinalizeRound();
      }

      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        backprop_count_below,
        backprop_count_above,
        4,
        "Balancing when results were correct"
      );
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        backprop_count_below,
        0,
        2,
        "little should be backpropagated given correct results"
      );

      // now test the case that the results are generally too high; this should result in most features below the cutoff
      // being backpropagated
      for( size_t iteration( 0); iteration < size_t( 10); ++iteration)
      {
        backprop_count = backprop_count_above = backprop_count_below = 0;
        for( size_t i( 0); i < nr_features; ++i)
        {
          linal::Vector< float> row( size_t( 1), shuffled_results( i) + float( 0.3));
          linal::Vector< float> error( float( -1.0) * results.GetRow( i));
          if( balanced->ShouldBackpropagate( row, error, i, 0))
          {
            ++backprop_count;
            if( results.GetRow( i)( 0) < scaled_cutoff)
            {
              ++backprop_count_below;
            }
            else
            {
              ++backprop_count_above;
            }
          }
        }
        // call finalize
        balanced->FinalizeRound();
      }

      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        backprop_count,
        backprop_count_below,
        4,
        "Rebalancing when average is far off should favor features on the side of the cutoff furthest from the average"
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelNeuralNetworkSelectiveBackpropagation

  const ExampleClass::EnumType ExampleModelNeuralNetworkSelectiveBackpropagationHybrid::s_Instance
  (
    GetExamples().AddEnum( ExampleModelNeuralNetworkSelectiveBackpropagationHybrid())
  );

} // namespace bcl
