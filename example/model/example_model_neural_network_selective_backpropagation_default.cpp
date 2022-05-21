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
#include "model/bcl_model_neural_network_selective_backpropagation_default.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_template_instantiations.h"
#include "model/bcl_model_objective_function_rmsd.h"
#include "model/bcl_model_retrieve_data_set_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_neural_network_selective_backpropagation_default.cpp
  //!
  //! @author mendenjl
  //! @date Jun 06, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelNeuralNetworkSelectiveBackpropagation :
    public ExampleInterface
  {
  public:

    ExampleModelNeuralNetworkSelectiveBackpropagation *Clone() const
    {
      return new ExampleModelNeuralNetworkSelectiveBackpropagation( *this);
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

      // construct the selective backpropagator
      model::NeuralNetworkSelectiveBackpropagationDefault backprop_selector;

      // initialize
      backprop_selector.Initialize( model_frds, model::ObjectiveFunctionRmsd());

      // test alias
      BCL_ExampleCheck( backprop_selector.GetAlias(), "All");

      // test that all features would be backpropagated and that the error vector is unchanged
      bool all_features_backpropagated( true), error_vector_unchanged( true);
      linal::Vector< float> error_vector( model_frds.GetResultsPtr()->GetMatrix().GetRow( 0));
      const linal::Vector< float> error_initial( error_vector);
      for( size_t i( 0), nr_features( model_frds.GetSize()); i < nr_features; ++i)
      {
        all_features_backpropagated =
          backprop_selector.ShouldBackpropagate
          (
            model_frds.GetResultsPtr()->GetMatrix().GetRow( i),
            error_vector,
            i,
            0
          );
        error_vector_unchanged = error_vector == error_initial;
        if( !all_features_backpropagated || !error_vector_unchanged)
        {
          // stop on error
          break;
        }
      }
      // call finalize; shouldn't do anything for this class
      backprop_selector.FinalizeRound();

      BCL_ExampleIndirectCheck
      (
        all_features_backpropagated,
        true,
        "ShouldBackpropagate always returns true"
      );
      BCL_ExampleIndirectCheck
      (
        error_vector_unchanged,
        true,
        "ShouldBackpropagate should not change the error vector for this class"
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelNeuralNetworkSelectiveBackpropagation

  const ExampleClass::EnumType ExampleModelNeuralNetworkSelectiveBackpropagation::s_Instance
  (
    GetExamples().AddEnum( ExampleModelNeuralNetworkSelectiveBackpropagation())
  );

} // namespace bcl
