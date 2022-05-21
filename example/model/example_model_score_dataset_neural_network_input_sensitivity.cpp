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
#include "model/bcl_model_score_dataset_neural_network_input_sensitivity.h"

// includes from bcl - sorted alphabetically
#include "model/bcl_model_retrieve_data_set_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_score_dataset_neural_network_input_sensitivity.cpp
  //!
  //! @author mendenjl
  //! @date Jul 01, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelScoreDatasetNeuralNetworkInputSensitivity :
    public ExampleInterface
  {
  public:

    ExampleModelScoreDatasetNeuralNetworkInputSensitivity *Clone() const
    {
      return new ExampleModelScoreDatasetNeuralNetworkInputSensitivity( *this);
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
      // create objects using the implementation framework, which is the only way to set this objects members
      const std::string model_directory( AddExampleInputPathToFilename( e_Model, ""));
      util::Implementation< model::ScoreDatasetInterface> input_sensitivity
      (
        "InputSensitivityNeuralNetwork(storage=File(directory=" + model_directory + ",prefix=model),key=000001,weights=(absolute=1))"
      );

      // ensure it is defined
      BCL_ExampleAssert( input_sensitivity.IsDefined(), true);

      // load in the example score dataset
      util::Implementation< model::RetrieveDataSetBase>
        dataset_retriever( "File(filename=" + model_directory + "example_data_set_score.bcl)");

      util::ShPtr< descriptor::Dataset> dataset( dataset_retriever->GenerateDataSet());

      // Calculate the input sensitivity; due to the balancing the result should be exactly 1/2 for the good descriptor (0) and
      // -1/2 for the other
      const float expected_sensitivity_col_b( 0.5);

      // the first column should get negative the second column in this case, since the first
      // column is useless, whereas the second column is good
      BCL_ExampleCheckWithinAbsTolerance
      (
        input_sensitivity->Score( *dataset)( 0),
        -expected_sensitivity_col_b,
        float( 1.0e-6)
      );

      // the second column should get 1.3125, give or take a little due to numerical error
      BCL_ExampleCheckWithinAbsTolerance
      (
        input_sensitivity->Score( *dataset)( 1),
        expected_sensitivity_col_b,
        float( 1.0e-6)
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelScoreDatasetNeuralNetworkInputSensitivity

  const ExampleClass::EnumType ExampleModelScoreDatasetNeuralNetworkInputSensitivity::s_Instance
  (
    GetExamples().AddEnum( ExampleModelScoreDatasetNeuralNetworkInputSensitivity())
  );

} // namespace bcl
