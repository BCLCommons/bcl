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

#ifndef BCL_MODEL_TRAIN_RESTRICTED_BOLTZMANN_MACHINE_LAYER_H_
#define BCL_MODEL_TRAIN_RESTRICTED_BOLTZMANN_MACHINE_LAYER_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set.h"
#include "bcl_model_interface.h"
#include "bcl_model_restricted_boltzmann_machine_layer.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_vector.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "math/bcl_math_range.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TrainRestrictedBoltzmannMachineLayer
    //! @brief A (primarily) storage class for all the gradients and data structures used for training a single
    //!        restricted boltzmann machine layer.  This class is intended to be used solely by
    //!        IterateRestrictedBoltzmannMachine.
    //!
    //! @see @link example_model_train_restricted_boltzmann_machine_layer.cpp @endlink
    //! @author mendenjl
    //! @date Aug 12, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API TrainRestrictedBoltzmannMachineLayer :
      public util::ObjectInterface
    {
    public:

    //////////
    // data //
    //////////

      //! the default input range for neural network transfer functions
      static const math::Range< float> s_DefaultInputRange;

    private:

    //////////
    // data //
    //////////

      //! Pointer to the layer that is being trained.  Constructing class owns the layers, so no need for ShPtr
      util::SiPtr< RestrictedBoltzmannMachineLayer> m_Layer;

      //! Number of stochastic propagations per feature
      //! 0 gives fast convergence but poor generalizability, 3 gives the most generality but can be slow for binomial
      //! networks. Hinton uses 1 for binomial networks to speed convergence; stochastic sigmoid should usually use 3.
      size_t m_StochasticStepCount;

      //! Gradients for weights, biases, and noise
      linal::Matrix< float> m_WeightGradient;        //!< Gradient of the weight matrix
      linal::Vector< float> m_BiasVisibleGradient;   //!< Gradient of the bias for visible neurons
      linal::Vector< float> m_BiasHiddenGradient;    //!< Gradient of the bias for hidden neurons
      linal::Vector< float> m_NoiseVisibleGradient;  //!< Gradient of the noise for visible; used for stochastic sigmoid
      linal::Vector< float> m_NoiseHiddenGradient;   //!< Gradient of the noise for hidden; used for stochastic sigmoid

      //! Arrays for gibbs sampling
      linal::Vector< float> m_Hidden;                //!< ref to vector that will hold hidden neuron activations
      linal::Vector< float> m_VisibleReconstructed;  //!< ref to vector that will hold the reconstructed visible layer
      linal::Vector< float> m_HiddenReconstructed;   //!< ref to vector that will hold the reconstructed hidden layer
      linal::Vector< float> m_VisibleSample;         //!< ref to vector that will hold visible layer sampling
      linal::Vector< float> m_HiddenSample;          //!< ref to vector that will hold hidden neuron sample

      //! Tally of the number of features seen since the last time the layer was update
      size_t                m_NumberFeaturesSeen;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      TrainRestrictedBoltzmannMachineLayer();

      //! @brief constructor from an untrained network
      //! @param LAYER the layer that this trainer will handle
      //! @param STOCHASTIC_STEP_COUNT helps control generalizability, see note above on m_StochasticStepCount
      TrainRestrictedBoltzmannMachineLayer
      (
        RestrictedBoltzmannMachineLayer &LAYER,
        const size_t &STOCHASTIC_STEP_COUNT = 3
      );

      //! copy constructor
      TrainRestrictedBoltzmannMachineLayer *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief core function : takes a feature and uses it to update the associated gradients
      //! @param INPUT the input feature of interest
      //! @return reference to the hidden layer output
      linal::VectorConstReference< float> Train( const linal::VectorConstInterface< float> &INPUT);

      //! @brief core function : takes a feature and a result and uses it to update the associated gradients
      //! @param INPUT the input feature of interest
      //! @param OUTPUT the output feature of interest
      //! @return reference to reconstructed output
      linal::VectorConstReference< float> Train
      (
        const linal::VectorConstInterface< float> &INPUT,
        const linal::VectorConstInterface< float> &OUTPUT
      );

      //! @brief Update layer to update the associated RestrictedBoltzmannMachineLayer's parameters
      //! @param WEIGHT_COST relative cost of weight
      //! @param WEIGHT_UPDATE update object for the weight
      //! @param VISIBLE_BIAS_UPDATE update object for the visible layer's bias
      //! @param HIDDEN_BIAS_UPDATE update object for the hidden layer's bias
      //! @param VISIBLE_NOISE_UPDATE update object for the visible layer's noise (for StochasticSigmoid RBMs)
      //! @param HIDDEN_NOISE_UPDATE update object for the hidden layer's noise (for StochasticSigmoid RBMs)
      void UpdateLayer
      (
        const float &WEIGHT_COST,
        NeuralNetworkUpdateWeightsInterface &WEIGHT_UPDATE,
        NeuralNetworkUpdateWeightsInterface &VISIBLE_BIAS_UPDATE,
        NeuralNetworkUpdateWeightsInterface &HIDDEN_BIAS_UPDATE,
        NeuralNetworkUpdateWeightsInterface &VISIBLE_NOISE_UPDATE,
        NeuralNetworkUpdateWeightsInterface &HIDDEN_NOISE_UPDATE
      );

      //! @brief Reset this object, to prepare for the next batch of data
      void Reset();

      //! @brief Accumulate changes from a different trainer of the same layer (used for reduction after thread-completion)
      //! @param OTHER the other trainer's whose changes should be added to this layers
      void AccumulateChangesFrom( const TrainRestrictedBoltzmannMachineLayer &OTHER);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read TrainRestrictedBoltzmannMachineLayer from std::istream
      std::istream &Read( std::istream &ISTREAM);

      //! write TrainRestrictedBoltzmannMachineLayer into std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class TrainRestrictedBoltzmannMachineLayer

  } // namespace model
} // namespace bcl

#endif //BCL_MODEL_TRAIN_RESTRICTED_BOLTZMANN_MACHINE_LAYER_H_
