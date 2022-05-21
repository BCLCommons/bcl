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

#ifndef BCL_MODEL_PRETRAIN_STACKED_AUTO_ENCODER_H_
#define BCL_MODEL_PRETRAIN_STACKED_AUTO_ENCODER_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_approximator_base.h"
#include "bcl_model_approximator_neural_network.h"
#include "bcl_model_feature_data_set.h"
#include "bcl_model_neural_network_perturbation_interface.h"
#include "bcl_model_neural_network_selective_backpropagation_interface.h"
#include "bcl_model_neural_network_update_weights_interface.h"
#include "bcl_model_objective_function_wrapper.h"
#include "bcl_model_pretrain_neural_network_interface.h"
#include "bcl_model_transfer_function_interface.h"
#include "linal/bcl_linal_vector.h"
#include "linal/bcl_linal_vector_nd.h"
#include "math/bcl_math_running_average.h"
#include "opti/bcl_opti_tracker.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PretrainStackedAutoEncoder
    //! @brief provides the iteration function for training a neural network with the
    //!        resilient propagation update algorithm, also knows how to score the resultant argument.
    //!
    //! @see @link example_model_pretrain_stacked_auto_encoder.cpp @endlink
    //! @author butkiem1, mendenjl
    //! @date Dec 05, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PretrainStackedAutoEncoder :
      public PretrainNeuralNetworkInterface
    {
    private:

    //////////
    // data //
    //////////

      ApproximatorNeuralNetwork m_Approximator;

      //! shared pointer to the training data
      util::ShPtr< descriptor::Dataset> m_TrainingData;

      //! directory name to store the auto encoder model
      std::string m_AutoEncoderStoragePath;

      //! encoded training data set
      util::ShPtr< FeatureDataSet< float> > m_EncodedFeatures;

      // track the current assembled AE network
      util::ShPtr< NeuralNetwork> m_CurrentAutoEncoderModel;

      // track the current assembled AE network
      util::ShPtr< NeuralNetwork> m_CurrentAutoDecoderModel;

      //! # of neurons in each hidden layer
      storage::Vector< size_t> m_AutoEncoderArchitecture;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_PretrainInstance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PretrainStackedAutoEncoder();

      //! @brief constructor from training data, transfer, rescale, and objective functions
      //! @param TRAINING_DATA data to train the NeuralNetwork on
      PretrainStackedAutoEncoder
      (
        util::ShPtr< descriptor::Dataset> &TRAINING_DATA
      );

      //! @brief copy constructor
      //! @return a new PretrainStackedAutoEncoder copied from this instance
      PretrainStackedAutoEncoder *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief set training data set for a specific iterate in approximater framework
      //! @param DATA training data set
      void SetTrainingData( util::ShPtr< descriptor::Dataset> &DATA);

    ////////////////
    // operations //
    ////////////////

      //! @brief the main operation, pretrains a neural network
      //! @param DATA the data for use in pretraining
      //! @param OBJECTIVE ShPtr to the objective function for the network
      util::ShPtr< NeuralNetwork> PretrainNetwork
      (
        util::ShPtr< descriptor::Dataset> &DATA,
        const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE
      );

      //! @brief construct a model from the current iterate
      //! @return shptr to the new model interface
      util::ShPtr< Interface> GetCurrentModel() const;

      //! @brief conducts the next approximation step and stores the approximation
      void TrainNextLayer( const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE);

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

    }; // class IterateResilientProagation

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_PRETRAIN_STACKED_AUTO_ENCODER_H_
