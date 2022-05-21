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
#include "model/bcl_model_pretrain_neural_network_from_file.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "model/bcl_model_objective_function_wrapper.h"
#include "model/bcl_model_retrieve_data_set_base.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_pretrain_neural_network_from_file.cpp
  //!
  //! @author mendenjl
  //! @date Aug 14, 2013
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelPretrainNeuralNetworkFromFile :
    public ExampleInterface
  {
  public:

    ExampleModelPretrainNeuralNetworkFromFile *Clone() const
    {
      return new ExampleModelPretrainNeuralNetworkFromFile( *this);
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
      model::PretrainNeuralNetworkFromFile pretrainer_default;

      // try reading from a string
      util::Implementation< model::PretrainNeuralNetworkInterface> pretrainer_from_file_impl
      (
        "File(" + AddExampleInputPathToFilename( e_Model, "model000001.model") + ")"
      );

    /////////////////
    // data access //
    /////////////////

      // check whether the implementation is available
      BCL_ExampleIndirectCheck
      (
        pretrainer_from_file_impl.IsDefined(),
        true,
        "Should be able to create implementation with name File("
        + AddExampleInputPathToFilename( e_Model, "model000001.model") + ")"
      );

      // load the associated data set
      util::Implementation< model::RetrieveDataSetBase>
        dataset_retriever( "File(filename=" + AddExampleInputPathToFilename( e_Model, "example_data_set_score.bcl)"));
      util::ShPtr< descriptor::Dataset> fds( dataset_retriever->GenerateDataSet());
      util::ShPtr< model::ObjectiveFunctionWrapper> objective( new model::ObjectiveFunctionWrapper());

      // check whether the model can be loaded
      BCL_ExampleCheck( pretrainer_from_file_impl->PretrainNetwork( fds, objective).IsDefined(), true);
      BCL_ExampleCheck( pretrainer_from_file_impl->PretrainNetwork( fds, objective)->GetNumberInputs(), 2);
      BCL_ExampleCheck( pretrainer_from_file_impl->PretrainNetwork( fds, objective)->GetNumberOutputs(), 1);

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelPretrainNeuralNetworkFromFile

  const ExampleClass::EnumType ExampleModelPretrainNeuralNetworkFromFile::s_Instance
  (
    GetExamples().AddEnum( ExampleModelPretrainNeuralNetworkFromFile())
  );

} // namespace bcl
