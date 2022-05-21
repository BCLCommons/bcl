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

#ifndef BCL_MODEL_APPROXIMATOR_NEURAL_NETWORK_SELECTIVE_H_
#define BCL_MODEL_APPROXIMATOR_NEURAL_NETWORK_SELECTIVE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_approximator_base.h"
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
#include "random/bcl_random_uniform_distribution.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ApproximatorNeuralNetworkSelective
    //! @brief provides the iteration function for training a neural network with the
    //!        resilient propagation update algorithm, also knows how to score the resultant argument.
    //!
    //! @see @link example_model_approximator_neural_network_selective.cpp @endlink
    //! @author kothiwsk
    //! @date Apr 20, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ApproximatorNeuralNetworkSelective :
      public ApproximatorBase,
      public PretrainNeuralNetworkInterface
    {
    public:

      //! Methods of applying dropout to the input layer
      enum InputLayerDropoutType
      {
        e_Zero,                        //!< Set the dropped neurons to zero; fast, but biased
        e_Noise,                       //!< Set the neurons to a random gaussian value sampled according to the input mean/std
        e_CopySingleRandomFeature,     //!< Copy dropped values from a single randomly-selected feature
        e_CopyEachRandomFeature,       //!< Copy each dropped value from a randomly-selected feature (training example)
        e_CopySingleRandomPeerFeature, //!< Copy values from a peer; for classification tasks, copies values from a
                                       //!< randomly selected feature of the same class
        e_CopyEachRandomPeerFeature,   //!< Copy values from a peer; for classification tasks, copies each value from a
                                       //!< randomly selected feature of the same class
        e_CopySingleBalancedFeature,   //!< Copy values from a uniform randomly selected feature of a uniform randomly
                                       //!< selected class
        e_CopyEachBalancedFeature,     //!< Copy each value from a uniform randomly selected feature of a uniform
                                       //!< randomly selected feature of the same class
        s_NumberInputLayerDropoutTypes,
        s_InputLayerDropoutFirstClassBasedMethod = e_CopySingleRandomPeerFeature
      };

      //! @brief InputLayerDropoutType as string
      //! @param TYPE the type
      //! @return the string for the type
      static const std::string &GetInputLayerDropoutTypeName( const InputLayerDropoutType &TYPE);

      //! @brief InputLayerDropoutType enum I/O helper
      typedef util::WrapperEnum< InputLayerDropoutType, &GetInputLayerDropoutTypeName, s_NumberInputLayerDropoutTypes>
        InputLayerDropoutTypeEnum;

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

      //! Means of generating the initial network (by default, randomizes weights)
      util::Implementation< PretrainNeuralNetworkInterface> m_Pretrainer;

      //! true to shuffle feature / result ordering after every run through the data (iteration)
      bool m_Shuffle;

      //! # of features updated last round by each thread
      linal::Vector< size_t> m_NumberFeaturesUpdated;

      //! connection density; allows for sparse network formation
      float m_ConnectionDensity;

      //! whether to align cutoff to objective function, if possible
      bool m_AlignCutoff;

      //! whether the output should be rescaled to be within the dynamic output range of the transfer function
      //! This is usually best when working on regression problems. For classification, it is often preferable to
      //! rescale to the complete output range of the transfer function
      bool m_RescaleOutputDynamicRange;

      // Thread local variables
      //! last slopes of bias (one for each thread)
      storage::Vector< storage::Vector< linal::Vector< float> > > m_SlopesBias;

      //! last slopes of weights (one for each thread)
      //! Mutable to allow use as temporary storage in GetCurrentModel
      mutable storage::Vector< storage::Vector< linal::Matrix< float> > > m_SlopesWeight;

      //! last error values (one for each thread)
      storage::Vector< storage::Vector< linal::Vector< float> > > m_Errors;

      //! last hidden data (one for each thread)
      storage::Vector< storage::Vector< linal::Vector< float> > > m_Hidden;

      //! When performing dropout, last perturbed feature vector
      storage::Vector< linal::Vector< float> > m_FeatureWithDropout;

      //! last hidden input data (one for each thread)
      storage::Vector< storage::Vector< linal::Vector< float> > > m_HiddenInput;

      //! vectors of data set ranges (one for each thread)_SlopesWeight;
      storage::Vector< math::Range< size_t> > m_DataSetRanges;

      //! sum of squared errors for each thread
      linal::Vector< float> m_RMSDError;

      // Thread independent variables
      //! position in m_DataSetRanges
      size_t m_DataSetRangePosition; //! Note, if the data set is written out and read back in, the data set range position is reset to 0

      //! # of neurons in each hidden layer
      storage::Vector< size_t> m_HiddenArchitecture;

      //! Fraction of neurons to have "dropout" (set to 0) every round, for each non-output layer
      linal::Vector< float>    m_Dropout;

      //! Number of partitions in the dropout scheme.  An equal number of neurons will be dropped from each partition
      linal::Vector< size_t>    m_DropoutPartitions;

      //! Dropped columns this round, one vector for each non-output layer
      storage::Vector< storage::Vector< size_t> > m_ChosenDropped;
      storage::Vector< storage::Vector< storage::Vector< size_t> > > m_NeuronIndices; //!< Indices of all available neurons in the range
      linal::Vector< size_t>                      m_NumberToDrop;  //!< Number of neurons to drop each round per layer.  If partiions are

      //! stores the bias of the connection between different layers (one vector for connections between each consecutive layers)
      storage::Vector< linal::Vector< float> > m_Bias;

      //! stores the weight of the connection between different neurons (one matrix for connections between each consecutive layers)
      storage::Vector< linal::Matrix< float> > m_Weight;

      //! transfer function used in every neuron
      util::Implementation< TransferFunctionInterface> m_TransferFunction;

      //! rescale output from last round
      util::ShPtr< RescaleFeatureDataSet> m_RescaleOutputLastRound;

      //! weight updater implementation
      util::Implementation< NeuralNetworkUpdateWeightsInterface> m_WeightUpdateType;

      //! Bias update type (sometimes best to make this simple, to allow biases to evolve faster than weights)
      util::Implementation< NeuralNetworkUpdateWeightsInterface> m_BiasUpdateType;

      //! weight update after every complete run through the data (optional)
      util::Implementation< NeuralNetworkPerturbationInterface> m_IterationWeightUpdateType;

      //! updaters for each bias
      util::ShPtrVector< NeuralNetworkUpdateWeightsInterface> m_BiasUpdaters;

      //! updaters for each weight
      util::ShPtrVector< NeuralNetworkUpdateWeightsInterface> m_WeightUpdaters;

      //! updaters to be called after every complete run through the data
      util::ShPtrVector< NeuralNetworkPerturbationInterface> m_IterationWeightUpdaters;

      //! Selective backpropagation method
      util::Implementation< NeuralNetworkSelectiveBackpropagationInterface> m_DataSelector;

      //! order in which to visit the features (modified from consecutive if m_Shuffle is on)
      storage::Vector< size_t> m_Order;

      //! result class for each feature.  Used only with peer-based varieties of input dropout
      storage::Vector< size_t> m_ResultClass;

      //! features that have the same result class.  Used only with peer-based varieties of input dropout
      storage::Vector< storage::Vector< size_t> > m_PeerFeatures;

      //! Bool: whether this iterate is just to perform pretraining
      bool m_IsPretrainer;

      //! Number of iterations, only used if this class is being treated as a pretrainer
      size_t m_NumIterations;

      //! Bool: whether to balance the data for classification tasks
      bool m_Balance;

      //! For balancing: max number of repeated features
      size_t m_BalanceMaxRepeatedFeatures;

      //! For balancing: this is the ratio the rarer class (positives/negatives) will be oversampled to
      float m_BalanceMaxOversampling;

      //! Input feature noise (Z-score scale)
      float m_NoiseZScore;

      //! Type / method of applying dropout to the input layer
      InputLayerDropoutTypeEnum m_InputDropoutType;

      //! One uniform distribution per thread
      storage::Vector< random::UniformDistribution> m_ThreadRandomNumberGenerators;

      //! scaling type
      RescaleFeatureDataSet::TypeEnum m_RescaleType;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_IterateInstance;
      static const util::SiPtr< const util::ObjectInterface> s_PretrainInstance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param PRETRAIN whether this object is being constructed as a pretrainer or an iterate
      ApproximatorNeuralNetworkSelective( const bool &PRETRAIN = false);

      //! @brief constructor from training data, transfer, rescale, and objective functions
      //! @param TRAINING_DATA data to train the NeuralNetwork on
      //! @param UPDATE_EVERY_NTH_FEATURE how often the weights get updated
      //! @param ARCHITECTURE the # of neurons in each hidden layer of the network
      //! @param TRANSFER_FUNCTION ShPtr to the transfer function between input and output of each neuron
      //! @param OBJECTIVE_FUNCTION ShPtr to objective function
      //! @param WEIGHT_UPDATE_FUNCTION method by which to update the weights
      //! @param ITERATIONS_PER_RMSD_REPORT # iterations per report of the rmsd
      ApproximatorNeuralNetworkSelective
      (
        util::ShPtr< descriptor::Dataset> &TRAINING_DATA,
        const size_t UPDATE_EVERY_NTH_FEATURE,
        const storage::Vector< size_t> &ARCHITECTURE,
        const util::Implementation< TransferFunctionInterface> &TRANSFER_FUNCTION,
        const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE_FUNCTION,
        const util::Implementation< NeuralNetworkUpdateWeightsInterface> &WEIGHT_UPDATE_FUNCTION,
        const size_t &ITERATIONS_PER_RMSD_REPORT = 1,
        const RescaleFeatureDataSet::TypeEnum m_RescaleType = RescaleFeatureDataSet::e_AveStd
      );

      //! @brief copy constructor
      //! @return a new ApproximatorNeuralNetworkSelective copied from this instance
      ApproximatorNeuralNetworkSelective *Clone() const;

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

      //! @brief set the hidden architecture
      void SetHiddenArchitecture( const storage::Vector< size_t> &HIDDEN)
      {
        m_HiddenArchitecture = HIDDEN;
      }

      //! @brief get transfer function used in every neuron
      const storage::Vector< size_t> &GetHiddenArchitecture() const
      {
        return m_HiddenArchitecture;
      }

      //! @brief get transfer function used in every neuron
      const util::Implementation< TransferFunctionInterface> &GetTransferFunction() const
      {
        return m_TransferFunction;
      }

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

      //! @brief returns the current approximation
      //! @return current argument result pair
      const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> > GetCurrentApproximation() const;

      //! @brief conducts the next approximation step and stores the approximation
      void Next();

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

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

      //! @brief evaluates whether the approximation can continue
      //! @return true, if the approximation can continue - otherwise false
      bool CanContinue() const
      {
        return true;
      }

      //! run forward through ANN and compute err terms m_Hidden (test)
      void CalcHiddenTerms( const FeatureReference< float> &FEATURE, const size_t &THREAD_ID, const size_t &FEATURE_ID);

      //! run backward through ANN and compute err terms m_Errors (train)
      //! @return true if the error should be backpropagated
      bool CalcErrorTerms( const FeatureReference< float> &RESULT, const size_t &THREAD_ID, const size_t &FEATURE_ID);

      //! compute changes m_Changes to be applied on ANN (train)
      void CalcChangeTerms( const FeatureReference< float> &FEATURE, const size_t &THREAD_ID);

      //! train a thread with its data
      void TrainThread( const size_t &THREAD_ID);

      //! @brief update the weights and bias
      void UpdateWeights();

      //! @brief update dropped neurons
      void UpdateDroppedNeurons();

      //! @brief sets up the neural network with a particular architecture, after training data was set
      void SetupArchitecture();

      //! @brief sets up the weight updaters and change weights
      //!        this has to be done unless the iterate is read from a file
      void InitializeWeightUpdaters();

      //! @brief sets the data set ranges up for the # of threads
      void SetupDataSetRanges();

      //! @brief balance features based on the objective function
      //! @param MAX_REPEATS the maximum # of repeats allowed for any given feature; used for multi-column outputs, when
      //!        otherwise an output column with extremely few positives would see severe over-representation of certain
      //!        features
      //! @param MAX_BALANCE_RATIO 0-1; indicates the maximum ratio that the under-represented features will be balanced
      //! If underepresented features are already balanced to at least MAX_BALANCE_RATIO before calling this function,
      //! no features will be removed.  If the value of MAX_BALANCE_RATIO is greater than MAX_REPEATS would allow for the
      //! column; then each feature will just be replicated MAX_REPEATS times
      void BalanceFeatures
      (
        const size_t &MAX_REPEATS,
        const float &MAX_BALANCE_RATIO = float( 1.0)
      );

      //! @brief helper function to initialize dropout-related vectors
      void InitializeDropout();

      //! @brief determine result classes; only called if input dropout type requires peers
      void DetermineResultClasses();

    }; // class IterateResilientProagation

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_APPROXIMATOR_NEURAL_NETWORK_SELECTIVE_H_
