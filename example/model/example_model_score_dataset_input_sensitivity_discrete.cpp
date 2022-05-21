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
#include "model/bcl_model_score_dataset_input_sensitivity_discrete.h"

// includes from bcl - sorted alphabetically
#include "model/bcl_model_retrieve_data_set_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_score_dataset_input_sensitivity_discrete.cpp
  //!
  //! @author vuot2, mendenjl
  //! @date Feb 06, 2018
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelScoreDatasetInputSensitivityDiscrete :
    public ExampleInterface
  {
  public:

    ExampleModelScoreDatasetInputSensitivityDiscrete *Clone() const
    {
      return new ExampleModelScoreDatasetInputSensitivityDiscrete( *this);
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
      //! Create sensitity score models based on 4 different kinds of discrete derivatives
      //! the formula of the function is f(x) = x1 - 2x2
      const std::string model_directory( AddExampleInputPathToFilename( e_Model, ""));
      util::Implementation< model::ScoreDatasetInterface> model_decrement
      (
        "InputSensitivityDiscrete("
        "storage=File(directory=" + model_directory + ", prefix=discrete_linear ),"
        "key=0,"
        "weights=(average=1,balance=False,categorical=False),"
        "derivative=Decrement)"
      );

      util::Implementation< model::ScoreDatasetInterface> model_increment
      (
        "InputSensitivityDiscrete("
        "storage=File(directory=" + model_directory + ", prefix=discrete_linear ),"
        "key=0,"
        "weights=(average=1,balance=False,categorical=False),"
        "derivative=Increment)"
      );

      util::Implementation< model::ScoreDatasetInterface> model_center
      (
        "InputSensitivityDiscrete("
        "storage=File(directory=" + model_directory + ", prefix=discrete_linear ),"
        "key=0,"
        "weights=(average=1,balance=False,categorical=False),"
        "derivative=Center)"
      );

      util::Implementation< model::ScoreDatasetInterface> model_direction
      (
        "InputSensitivityDiscrete("
        "storage=File(directory=" + model_directory + ", prefix=discrete_linear ),"
        "key=0,"
        "weights=(average=1,balance=False,categorical=False),"
        "derivative=Direction)"
      );

      // Train ANN to approximate the function f(x) = x1 - x2^2
      util::Implementation< model::ScoreDatasetInterface> ann_center
      (
        "InputSensitivityDiscrete("
        "storage=File(directory=" + model_directory + ", prefix=discrete_test ),"
        "key=0,"
        "weights=(average=1,balance=False,categorical=False),"
        "derivative=Center)"
      );

      // ensure it is defined
      BCL_ExampleAssert( model_increment.IsDefined(), true);

      // load in the example score dataset with 1 pair x1, x2
      util::Implementation< model::RetrieveDataSetBase>
        dataset_retriever( "File(filename=" + model_directory + "example_data_set_1_0.bcl)");

      util::ShPtr< descriptor::Dataset> dataset( dataset_retriever->GenerateDataSet());

      // load in the example score dataset with 1 pair x1, x2
      util::Implementation< model::RetrieveDataSetBase>
        dataset_retriever_ann( "File(filename=" + model_directory + "example_data_set_1_1.bcl)");

      util::ShPtr< descriptor::Dataset> dataset_ann( dataset_retriever_ann->GenerateDataSet());

      // Check the increment derivative
      BCL_ExampleCheckWithinAbsTolerance
      (
        model_increment->Score( *dataset)( 0),
        1,
        float( 1.0e-6)
      );
      BCL_ExampleCheckWithinAbsTolerance
      (
        model_increment->Score( *dataset)( 1),
        -2,
        float( 1.0e-6)
      );

      // Check the decrement derivative
      BCL_ExampleCheckWithinAbsTolerance
      (
        model_decrement->Score( *dataset)( 0),
        -1,
        float( 1.0e-6)
      );
      BCL_ExampleCheckWithinAbsTolerance
      (
        model_decrement->Score( *dataset)( 1),
        2,
        float( 1.0e-6)
      );

      // Check the center derivative
      BCL_ExampleCheckWithinAbsTolerance
      (
        model_center->Score( *dataset)( 0),
        -1,
        float( 1.0e-6)
      );
      BCL_ExampleCheckWithinAbsTolerance
      (
        model_center->Score( *dataset)( 1),
        -2,
        float( 1.0e-6)
      );

      // Check the center derivative with ANN
      BCL_ExampleCheckWithinAbsTolerance
      (
        ann_center->Score( *dataset_ann)( 0),
        -1,
        float( 0.05)
      );
      BCL_ExampleCheckWithinAbsTolerance
      (
        ann_center->Score( *dataset_ann)( 1),
        2,
        float( 0.05)
      );

      // Check the direction derivative
      BCL_ExampleCheckWithinAbsTolerance
      (
        model_direction->Score( *dataset)( 0),
        1,
        float( 1.0e-6)
      );
      BCL_ExampleCheckWithinAbsTolerance
      (
        model_direction->Score( *dataset)( 1),
        -2,
        float( 1.0e-6)
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelScoreDatasetInputSensitivityDiscrete

  const ExampleClass::EnumType ExampleModelScoreDatasetInputSensitivityDiscrete::s_Instance
  (
    GetExamples().AddEnum( ExampleModelScoreDatasetInputSensitivityDiscrete())
  );
} // namespace bcl
