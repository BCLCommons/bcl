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
#include "model/bcl_model_objective_function_information_gain_ratio.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "model/bcl_model_neural_network.h"
#include "model/bcl_model_objective_function_wrapper.h"
#include "model/bcl_model_retrieve_data_set_base.h"
#include "model/bcl_model_retrieve_interface.h"
#include "model/bcl_model_store_interface.h"
#include "model/bcl_model_transfer_sigmoid.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_objective_function_information_gain_ratio.cpp
  //!
  //! @author butkiem1
  //! @date May 29, 2013
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelObjectiveFunctionInformationGainRatio :
    public ExampleInterface
  {
  public:

    ExampleModelObjectiveFunctionInformationGainRatio *Clone() const
    {
      return new ExampleModelObjectiveFunctionInformationGainRatio( *this);
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
      std::string example_path( AddExampleInputPathToFilename( e_Model, ""));

      // retrieve a given model

      //! construct a storage interface implementation
      util::Implementation< model::RetrieveInterface> model_retriever
      (
        "File( directory=" + example_path + ", prefix=model)"
      );

      model::RetrieveInterface::t_ModelPtr neural_network_model( model_retriever->Retrieve( "11"));

      // retrieve a given test data set
      util::Implementation< model::RetrieveDataSetBase> dataset_retriever
      (
        util::ObjectDataLabel( "Subset( filename=" + example_path + "dataset_aid891_100act_100inact.bin)")
      );

      util::ShPtr< descriptor::Dataset> dataset = dataset_retriever->GenerateDataSet();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      model::ObjectiveFunctionInformationGainRatio objective_function_default;

      //! constructor with parameters with specialization (here objective function rmsd)
      model::ObjectiveFunctionWrapper objective_function_default_by_label
      (
        dataset, util::ObjectDataLabel( "InformationGainRatio")
      );

      //! constructor with parameters with specialization
      model::ObjectiveFunctionWrapper objective_function
      (
        dataset, util::ObjectDataLabel( "InformationGainRatio(cutoff=4,parity=1)")
      );

      //! test clone
      util::ShPtr< model::ObjectiveFunctionWrapper> objective_function_clone( objective_function.Clone());

    /////////////////
    // data access //
    /////////////////

      // rescale function for in an output and denormalization
      util::ShPtr< model::RescaleFeatureDataSet> rescale_in
      (
        new model::RescaleFeatureDataSet( *dataset->GetFeaturesPtr(), model::NeuralNetwork::s_DefaultInputRange)
      );

      // set dataset in objective function
      objective_function_default_by_label.SetData( dataset, rescale_in);

      // check get and set methods
      BCL_ExampleAssert( objective_function_default_by_label.GetData()->GetSize(), size_t( 200));

    ///////////////
    // operators //
    ///////////////

      // check get and set methods
      BCL_ExampleIndirectCheck
      (
        math::EqualWithinTolerance
        (
          objective_function( neural_network_model),
          float( 0.76),
          0.01
        ),
        true,
        "operator(), expected the best information gain ratio at PPV 0.76, but got "
        + util::Format()( objective_function( neural_network_model))
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write bcl object
      WriteBCLObject( objective_function_default);
      // create default object
      model::ObjectiveFunctionInformationGainRatio objective_function_read;
      // read bcl object
      ReadBCLObject( objective_function_read);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelObjectiveFunctionInformationGainRatio

  const ExampleClass::EnumType ExampleModelObjectiveFunctionInformationGainRatio::s_Instance
  (
    GetExamples().AddEnum( ExampleModelObjectiveFunctionInformationGainRatio())
  );

} // namespace bcl
