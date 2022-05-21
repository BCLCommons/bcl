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

#ifndef BCL_MODEL_APPROXIMATOR_RESTRICTED_BOLTZMANN_MACHINE_H_
#define BCL_MODEL_APPROXIMATOR_RESTRICTED_BOLTZMANN_MACHINE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_approximator_base.h"
#include "bcl_model_feature_data_set.h"
#include "bcl_model_neural_network_perturbation_interface.h"
#include "bcl_model_neural_network_update_weights_interface.h"
#include "bcl_model_objective_function_wrapper.h"
#include "bcl_model_pretrain_neural_network_interface.h"
#include "bcl_model_restricted_boltzmann_machine_layer.h"
#include "bcl_model_train_restricted_boltzmann_machine_layer.h"
#include "bcl_model_transfer_function_interface.h"
#include "linal/bcl_linal_vector.h"
#include "linal/bcl_linal_vector_nd.h"
#include "math/bcl_math_running_average.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorRestrictedBoltzmannMachine
    //! @brief provides the iteration function for training a restricted boltzmann machine
    //!
    //! @see @link example_model_approximator_restricted_boltzmann_machine.cpp @endlink
    //! @author mendenjl
    //! @date Aug 13, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ApproximatorRestrictedBoltzmannMachine :
      public ApproximatorBase,
      public PretrainNeuralNetworkInterface
    {
    private:

    //////////
    // data //
    //////////

      //! number of training features shown before updating weights
      //! 0 -> view the whole data set before updating weights
      size_t m_UpdateEveryNthFeature;

      //! number of threads
      size_t m_NumberThreads;

      //! number of batch iterations after which to report the rmsd
      size_t m_IterationsPerRMSDMessage;

      //! Network type to create
      RestrictedBoltzmannMachineLayer::TypeEnum m_Type;

      //! true to shuffle feature / result ordering after every run through the data (iteration)
      bool m_Shuffle;

      //! number of stochastic propagations per gibbs sample
      size_t m_NumberStochasticSteps;

      //! All the layers of the RBM
      storage::Vector< RestrictedBoltzmannMachineLayer> m_Layers;

      // Thread local variables
      //! Trainers, outer vector is per-thread, inner vector is per layer
      storage::Vector< storage::Vector< TrainRestrictedBoltzmannMachineLayer> > m_Trainers;

      //! vectors of data set ranges (one for each thread)_SlopesWeight;
      storage::Vector< math::Range< size_t> > m_DataSetRanges;

      //! sum of squared errors for each thread
      linal::Vector< float> m_RMSDError;
      linal::Vector< float> m_ReconstructionError; //!< Errror in first layer reconstruction

      // Thread independent variables
      //! position in m_DataSetRanges
      size_t m_DataSetRangePosition; //! Note, if the data set is written out and read back in, the data set range position is reset to 0

      //! # of neurons in each hidden layer
      storage::Vector< size_t> m_HiddenArchitecture;

      //! weight updater implementation
      util::Implementation< NeuralNetworkUpdateWeightsInterface> m_WeightUpdateType;

      //! Bias and noise update type (sometimes best to make this simple, to allow biases to evolve faster than weights)
      util::Implementation< NeuralNetworkUpdateWeightsInterface> m_BiasUpdateType;

      //! Weight decay parameter
      float m_WeightDecay;

      //! updaters for each bias
      util::ShPtrVector< NeuralNetworkUpdateWeightsInterface> m_HiddenBiasUpdaters;
      util::ShPtrVector< NeuralNetworkUpdateWeightsInterface> m_VisibleBiasUpdaters;
      util::ShPtrVector< NeuralNetworkUpdateWeightsInterface> m_HiddenNoiseUpdaters;
      util::ShPtrVector< NeuralNetworkUpdateWeightsInterface> m_VisibleNoiseUpdaters;

      //! updaters for each weight
      util::ShPtrVector< NeuralNetworkUpdateWeightsInterface> m_WeightUpdaters;

      //! order in which to visit the features (modified from consecutive if m_Shuffle is on)
      storage::Vector< size_t> m_Order;

      //! Bool: whether this iterate is just to perform pretraining
      bool m_IsPretrainer;

      //! Number of iterations, only used if this class is being treated as a pretrainer
      size_t m_NumIterations;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_IterateInstance;
      static const util::SiPtr< const util::ObjectInterface> s_PretrainInstance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param PRETRAIN whether this object is being constructed as a pretrainer or an iterate
      ApproximatorRestrictedBoltzmannMachine( const bool &PRETRAIN = false);

      //! @brief constructor from training data, transfer, rescale, and objective functions
      //! @param TRAINING_DATA data to train the NeuralNetwork on
      //! @param UPDATE_EVERY_NTH_FEATURE how often the weights get updated
      //! @param ARCHITECTURE the # of neurons in each hidden layer of the network
      //! @param TRANSFER_FUNCTION ShPtr to the transfer function between input and output of each neuron
      //! @param OBJECTIVE_FUNCTION ShPtr to objective function
      //! @param WEIGHT_UPDATE_FUNCTION method by which to update the weights
      //! @param BIAS_UPDATE_FUNCTION method by which to update the biases
      //! @param SHUFFLE set to true to automatically shuffle order of example training every iteration
      //! @param ITERATIONS_PER_RMSD_REPORT # iterations per report of the rmsd
      //! @param TYPE type of RBM to train
      ApproximatorRestrictedBoltzmannMachine
      (
        util::ShPtr< descriptor::Dataset> &TRAINING_DATA,
        const size_t UPDATE_EVERY_NTH_FEATURE,
        const storage::Vector< size_t> &ARCHITECTURE,
        const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE_FUNCTION,
        const util::Implementation< NeuralNetworkUpdateWeightsInterface> &WEIGHT_UPDATE_FUNCTION,
        const util::Implementation< NeuralNetworkUpdateWeightsInterface> &BIAS_UPDATE_FUNCTION,
        const size_t &ITERATIONS_PER_RMSD_REPORT = 1,
        const RestrictedBoltzmannMachineLayer::Type &TYPE = RestrictedBoltzmannMachineLayer::e_StochasticSigmoid
      );

      //! @brief copy constructor
      //! @return a new ApproximatorRestrictedBoltzmannMachine copied from this instance
      ApproximatorRestrictedBoltzmannMachine *Clone() const;

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

      //! @brief construct a model from the current iterate
      //! @return shptr to the new model interface
      util::ShPtr< Interface> GetCurrentModel() const;

      //! @brief the main operation, pretrains a neural network
      //! @param DATA the data for use in pretraining
      //! @param OBJECTIVE ShPtr to the objective function for the network
      util::ShPtr< NeuralNetwork> PretrainNetwork
      (
        util::ShPtr< descriptor::Dataset> &DATA,
        const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE
      );

      //! @brief returns the current approximation
      //! @return current argument result pair
      const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> > GetCurrentApproximation() const;

      //! @brief conducts the next approximation step and stores the approximation
      void Next();

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief evaluates whether the approximation can continue
      //! @return true, if the approximation can continue - otherwise false
      bool CanContinue() const
      {
        return true;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    private:

      //! train a thread with its data
      void TrainThread( const size_t &THREAD_ID);

      //! @brief update the weights and bias
      void UpdateWeights();

      //! @brief sets up the weight updaters and change weights
      //!        this has to be done unless the iterate is read from a file
      void InitializeWeightUpdaters();

      //! @brief sets the data set ranges up for the # of threads
      void SetupDataSetRanges();

    }; // class IterateResilientProagation

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_APPROXIMATOR_RESTRICTED_BOLTZMANN_MACHINE_H_
