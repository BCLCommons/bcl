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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "model/bcl_model_pretrain_neural_network_from_file.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_parameter_check_file_existence.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> PretrainNeuralNetworkFromFile::s_Instance
    (
      util::Enumerated< PretrainNeuralNetworkInterface>::AddInstance( new PretrainNeuralNetworkFromFile())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new PretrainNeuralNetworkFromFile
    PretrainNeuralNetworkFromFile *PretrainNeuralNetworkFromFile::Clone() const
    {
      return new PretrainNeuralNetworkFromFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PretrainNeuralNetworkFromFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &PretrainNeuralNetworkFromFile::GetAlias() const
    {
      static const std::string s_name( "File");
      return s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief the main operation, pretrains a neural network
    //! @param DATA the data for use in pretraining
    //! @param OBJECTIVE ShPtr to the objective function for the network
    util::ShPtr< NeuralNetwork> PretrainNeuralNetworkFromFile::PretrainNetwork
    (
      util::ShPtr< descriptor::Dataset> &DATA,
      const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE
    )
    {
      BCL_Assert
      (
        m_NeuralNetwork->GetNumberInputs() == DATA->GetFeatureSize(),
        "Read in network had the wrong input size!"
      );
      BCL_Assert
      (
        m_NeuralNetwork->GetNumberOutputs() == DATA->GetResultSize(),
        "Read in network had the wrong output size!"
      );
      return m_NeuralNetwork;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool PretrainNeuralNetworkFromFile::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
    )
    {
      // read in the model
      io::IFStream input_network;
      io::File::MustOpenIFStream( input_network, m_FileName);
      util::ShPtr< Interface> model;
      io::Serialize::Read( model, input_network);

      // ensure that a model was read in
      if( !model.IsDefined())
      {
        ERROR_STREAM << "Network file had an undefined model";
        return false;
      }
      io::File::CloseClearFStream( input_network);

      // ensure that the model is a neural network
      m_NeuralNetwork = util::ShPtr< NeuralNetwork>( model);
      if( !m_NeuralNetwork.IsDefined())
      {
        ERROR_STREAM << "PretrainNeuralNetworkFromFile cannot load a model of type " << model->GetClassIdentifier();
        return false;
      }
      return true;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PretrainNeuralNetworkFromFile::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "load a neural network written out as a util::ShPtr<model::Interface>"
      );

      parameters.AddInitializer
      (
        "",
        "file name of an util::ShPtr<model::Interface> to a model::NeuralNetwork",
        io::Serialization::GetAgentWithCheck
        (
          &m_FileName,
          command::ParameterCheckFileExistence()
        )
      );

      return parameters;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace model

} // namespace bcl
