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

#ifndef BCL_MODEL_NEURAL_NETWORK_H_
#define BCL_MODEL_NEURAL_NETWORK_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set.h"
#include "bcl_model_interface.h"
#include "bcl_model_transfer_function_interface.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_range.h"
#include "sched/bcl_sched_mutex.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class NeuralNetwork
    //! @brief is an artificial neural network and is FunctionInterface derived
    //! @details It calculates for an argument vector (feature) of floats a prediction vector (result) with floats.
    //!
    //! @see @link example_model_neural_network.cpp @endlink
    //! @author meilerj, mueller, woetzen
    //! @date 12.04.2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API NeuralNetwork :
      public Interface
    {
    public:

    /////////////
    // friends //
    /////////////

    //////////
    // data //
    //////////

      //! the default input range for neural network transfer functions
      static const math::Range< float> s_DefaultInputRange;

    private:

    //////////
    // data //
    //////////

      //! transfer function used in every neuron
      util::Implementation< TransferFunctionInterface> m_TransferFunction;

      //! rescale input to transfer function input range
      util::ShPtr< RescaleFeatureDataSet> m_RescaleInput;

      //! rescale output from transfer function output range
      util::ShPtr< RescaleFeatureDataSet> m_RescaleOutput;

      //! bias
      storage::Vector< linal::Vector< float> > m_Bias;

      //! weight
      storage::Vector< linal::Matrix< float> > m_Weight;

      //! architecture of the neural network, stored for convenience
      storage::Vector< size_t> m_Architecture;

      //! Mutex for accessing m_Allocated and m_Available
      mutable sched::Mutex m_Mutex;

      //! array used during calculation; vectors hold hidden input and hidden vectors
      mutable storage::List< storage::Vector< linal::Vector< float> > > m_Allocated;
      mutable storage::List< storage::Vector< linal::Vector< float> > > m_Available;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      NeuralNetwork();

      //! @brief construct from all necessary parameters
      //! @param RESCALE_INPUT function to rescale input to transfer function range
      //! @param RESCALE_OUTPUT function to rescale output from transfer function range
      //! @param BIAS the bias for each neuron of a hidden layer
      //! @param WEIGHT the weight for each neuron of a hidden layer
      //! @param TRANSFER_FUNCTION function to simulate the threshold in the neuron
      NeuralNetwork
      (
        const util::ShPtr< RescaleFeatureDataSet> &RESCALE_INPUT,
        const util::ShPtr< RescaleFeatureDataSet> &RESCALE_OUTPUT,
        const storage::Vector< linal::Vector< float> > &BIAS,
        const storage::Vector< linal::Matrix< float> > &WEIGHT,
        const util::Implementation< TransferFunctionInterface> &TRANSFER_FUNCTION
      );

      //! @brief copy constructor
      //! @param NETWORK the nn that will be copied
      NeuralNetwork( const NeuralNetwork &NETWORK);

      //! copy constructor
      NeuralNetwork *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief get number layers
      size_t GetNumberLayers() const
      {
        return m_Weight.GetSize();
      }

      //! @brief returns architecture of this neural network
      //! @return first entry is the number of input neurons
      //!         middle entries are the number of hidden neurons in each layer
      //!         last number in vector is the number of output neurons
      const storage::Vector< size_t> &GetArchitecture() const;

      //! @brief set the architecture architecture of this neural network
      //! @param ARCHITECTURE first entry is the number of input neurons
      //!         middle entries are the number of hidden neurons in each layer
      //!         last number in vector is the number of output neurons
      void SetArchitecture( const storage::Vector< size_t> &ARCHITECTURE);

      //! get number inputs
      size_t GetNumberInputs() const;

      //! get number output neurons
      size_t GetNumberOutputs() const;

      //! @brief get the bias
      const storage::Vector< linal::Vector< float> > &GetBias() const
      {
        return m_Bias;
      }

      //! @brief set the bias
      void SetBias( const storage::Vector< linal::Vector< float> > &BIAS)
      {
        m_Bias = BIAS;
      }

      //! @brief get the weight
      const storage::Vector< linal::Matrix< float> > &GetWeight() const
      {
        return m_Weight;
      }

      //! @brief set the weight
      void SetWeight( const storage::Vector< linal::Matrix< float> > &WEIGHT)
      {
        m_Weight = WEIGHT;
      }

      //! @brief access the transfer function
      const util::Implementation< TransferFunctionInterface> &GetTransferFunction() const
      {
        return m_TransferFunction;
      }

      //! @brief access the rescale input function
      const util::ShPtr< RescaleFeatureDataSet> &GetRescaleInput() const
      {
        return m_RescaleInput;
      }

      //! @brief set the rescaling for the output
      //! @param RESCALING the new rescaling to use for the output
      void SetRescaleOutput( const util::ShPtr< RescaleFeatureDataSet> &RESCALING)
      {
        m_RescaleOutput = RESCALING;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Set the scaling of a feature set according to the model
      //! @param FEATURES feature set of interest
      //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
      //!       when operator() is called
      void Rescale( FeatureDataSet< float> &FEATURE) const;

      //! @brief predict result with model using a NOT rescaled feature vector
      //! @param FEATURE not rescaled feature vector
      //! @return predicted result vector using a model
      FeatureDataSet< float> PredictWithoutRescaling
      (
        const FeatureDataSetInterface< float> &FEATURE
      ) const;

      //! @brief predict result with model using a rescaled feature vector
      //! @param FEATURE normalized or rescaled feature vector
      //! @return predicted result vector using a model
      FeatureDataSet< float> operator()( const FeatureDataSetInterface< float> &FEATURE) const;

      //! @brief predict result with model and compute input sensitivity using a rescaled feature vector
      //! @param FEATURE feature of interest, MUST be rescaled
      //! @return predicted result (descaled) and input sensitivity (N_Descriptors X N_Results)
      storage::Pair< linal::Vector< float>, linal::Matrix< float> >
      ComputeResultInputSensitivity( const linal::VectorConstInterface< float> &FEATURE) const;

      //! @brief compute the result and deviation of that result using test-time dropout
      //! @param DATASET feature set of interest
      //! @param DROPOUT_RATIOS ratios of neurons to dropout from each layer (except output; excess layer ratios ignored)
      //! @param NREPEATS number of times to repeat dropout mask to get an estimate of the actual result
      storage::VectorND< 3, FeatureDataSet< float> > TestWithDropout
      (
        const FeatureDataSetInterface< float> &FEATURE,
        const storage::Vector< double> &DROPOUT_RATIOS,
        const size_t &NREPEATS
      ) const;

      //! @brief remove output layer and set last hidden layer to output layer
      void RemoveOutputLayer();

      //! @brief remove input layer and set last hidden layer to input layer
      void RemoveInputLayer();

      //! @brief append two networks if layers are compatible
      void Append( const util::ShPtr< NeuralNetwork> &NETWORK);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read NeuralNetwork from std::istream
      std::istream &Read( std::istream &ISTREAM);

      //! write NeuralNetwork into std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERROR_STREAM the stream to write errors to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

      //! @brief determines the architecture already present in the weight matrix/vectors of the neural network
      void SetImplicitArchitecture();

      //! @brief acquire a hidden input or hidden vector set (same size as m_Bias)
      //! @return iterator to the hidden input or hidden vector set
      storage::List< storage::Vector< linal::Vector< float> > >::iterator AcquireHiddenVectors() const;

      //! @brief release a given hidden vector set
      //! @param ITR iterator to the hidden input or hidden vector set
      void ReleaseHiddenVectors( const storage::List< storage::Vector< linal::Vector< float> > >::iterator &ITR) const;

      //! @brief function called by each thread to accomplish threaded forward propagation
      //! @param FEATURES features to train on
      //! @param STORAGE storage for the result
      void PredictWithThread
      (
        const FeatureDataSetInterface< float> &FEATURES,
        linal::MatrixInterface< float> &STORAGE
      ) const;

    }; // class NeuralNetwork

  } // namespace model
} // namespace bcl

#endif //BCL_MODEL_NEURAL_NETWORK_H_
