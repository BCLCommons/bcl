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

#ifndef BCL_MODEL_RESTRICTED_BOLTZMANN_MACHINE_LAYER_H_
#define BCL_MODEL_RESTRICTED_BOLTZMANN_MACHINE_LAYER_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_data_set.h"
#include "bcl_model_interface.h"
#include "bcl_model_neural_network_update_weights_interface.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_range.h"
#include "sched/bcl_sched_mutex.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RestrictedBoltzmannMachineLayer
    //! @brief A single layer RBM
    //!
    //! @see @link example_model_restricted_boltzmann_machine_layer.cpp @endlink
    //! @author mendenjl
    //! @date Aug 02, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RestrictedBoltzmannMachineLayer :
      public util::ObjectInterface
    {
    public:

      friend class TrainRestrictedBoltzmannMachineLayer;

    //////////
    // enum //
    //////////

      enum Type
      {
        e_StochasticSigmoid,    //!< Noisy sigmoid approx for continous data, see http://www.ee.nthu.edu.tw/~hchen/pubs/iee2003.pdf
        e_Binomial,             //!< Appropriate for binary data or hidden layers, see http://www.cs.toronto.edu/~hinton/science.pdf
        s_NumberTypes
      };

      // The method from http://www.iro.umontreal.ca/~lisa/publications2/index.php/attachments/single/24 has been tested
      // and found to have relatively poor convergence (with eta > 0.25, the method often diverges), and so has been
      // removed.  Also, it has the unfortunate effect of changing the transfer function from a sigmoid to a truncated
      // exponential, which is more expensive to calculate and is currently not implemented for neural networks.  If you
      // are interested in the code, though check out r4607 of this class

      //! @brief Type as string
      //! @param TYPE the type
      //! @return the string for the kernel
      static const std::string &GetTypeName( const Type &TYPE);

      //! @brief TypeEnum enum I/O helper
      typedef util::WrapperEnum< Type, &GetTypeName, s_NumberTypes> TypeEnum;

    //////////
    // data //
    //////////

      //! the default input range for neural network transfer functions
      static const math::Range< float> s_DefaultInputRange;

    private:

    //////////
    // data //
    //////////

      linal::Vector< float> m_BiasHidden;  //!< bias for hidden nodes
      linal::Vector< float> m_BiasVisible; //!< bias for visible nodes
      linal::Vector< float> m_SlopeHidden; //!< slope for hidden nodes
      linal::Vector< float> m_SlopeVisible; //!< slope for visible nodes

      linal::Matrix< float> m_Weight; //!< weight for each connection hidden - visible

      size_t                m_NumberVisibleNodes; //!< Cached # of visible nodes
      size_t                m_NumberHiddenNodes;  //!< Cached # of hidden nodes

      TypeEnum              m_NetworkType; //!< Type of network to employ

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RestrictedBoltzmannMachineLayer();

      //! @brief construct an untrained network
      //! @param NUMBER_VISIBLE number of visible neurons
      //! @param NUMBER_HIDDEN number of hidden neurons
      //! @param NETWORK_TYPE type of network in terms of transfer functions
      RestrictedBoltzmannMachineLayer
      (
        const size_t &NUMBER_VISIBLE,
        const size_t &NUMBER_HIDDEN,
        const Type   &NETWORK_TYPE = e_StochasticSigmoid
      );

      //! @brief construct from all necessary parameters
      //! @param BIAS_VISIBLE the bias for each neuron of the visible layer
      //! @param BIAS_HIDDEN the bias for each neuron of the hidden layer
      //! @param WEIGHT the weight for each neuron of a hidden layer
      //! @param NETWORK_TYPE type of network in terms of transfer functions
      RestrictedBoltzmannMachineLayer
      (
        const linal::Vector< float> &BIAS_VISIBLE,
        const linal::Vector< float> &BIAS_HIDDEN,
        const linal::Matrix< float> &WEIGHT,
        const Type                  &NETWORK_TYPE = e_StochasticSigmoid
      );

      //! copy constructor
      RestrictedBoltzmannMachineLayer *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns # of input neurons to this RBM
      //! @return number of input neurons
      const size_t &GetNumberInputNeurons() const
      {
        return m_NumberVisibleNodes;
      }

      //! @brief returns # of hidden neurons to this RBM
      //! @return # of hidden neurons to this RBM
      const size_t &GetNumberHiddenNeurons() const
      {
        return m_NumberHiddenNodes;
      }

      //! @brief set the architecture of this RBM
      //! @param NUMBER_INPUTS number of input/visible neurons
      //! @param NUMBER_OUTPUTS number of output/hidden neurons
      void SetArchitecture( const size_t &NUMBER_INPUTS, const size_t &NUMBER_OUTPUTS);

      //! @brief get the bias
      const linal::Vector< float> &GetHiddenBias() const
      {
        return m_BiasHidden;
      }

      //! @brief get the bias
      const linal::Vector< float> &GetVisibleBias() const
      {
        return m_BiasVisible;
      }

      //! @brief get the weight
      const linal::Matrix< float> &GetWeight() const
      {
        return m_Weight;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Forward propagate a given feature from the input to the hidden layer
      //! @param INPUT input feature of interest
      //! @param HIDDEN reference to vector that will be overwritten by the hidden layer
      //! @param STOCHASTIC whether to apply stochastic weighting onto the hidden layer
      void ForwardPropagate
      (
        const linal::VectorConstInterface< float> &INPUT,
        linal::VectorInterface< float> &HIDDEN,
        const bool STOCHASTIC
      ) const;

      //! @brief Back propagate a given feature from the hidden to the input layer
      //! @param INPUT input feature of interest
      //! @param HIDDEN hidden vector
      //! @param STOCHASTIC whether to apply stochastic weighting onto the input layer
      void BackPropagate
      (
        linal::VectorInterface< float> &INPUT,
        const linal::VectorInterface< float> &HIDDEN,
        const bool STOCHASTIC
      ) const;

      //! @brief perform a complete gibbs sampling step
      //! @param INPUT visible input
      //! @param HIDDEN ref to vector that will hold hidden neuron activations (not sampled) at completion
      //! @param INPUT_RECONSTRUCTED ref to vector to store reconstructed input layer
      //! @param HIDDEN_RECONSTRUCTED ref to vector to store reconstructed hidden layer
      //! @param STOCHASTIC_STEP_COUNT Number of propagations that should incorporate randomness (max = 3)
      //! @note STOCHASTIC_STEP_COUNT, 0 gives fast convergence but poor generalizability, 3 gives the most generality
      //!       but is very slow for binomial networks.
      //!       Hinton uses 1 for binomial networks to speed convergence; stochastic sigmoid should usually use 3,
      //!       Truncated exponential requires slow learning rate and/or setting STOCHASTIC_STEP_COUNT to 0
      //! @param HIDDEN_SAMPLE vector to hold sampling from hidden layer, used to avoid temporary memory allocation
      //! @param VISIBLE_SAMPLE vector to hold sampling from input layer, used to avoid temporary memory allocation
      void GibbsSample
      (
        const linal::VectorConstInterface< float> &INPUT,
        linal::VectorInterface< float> &HIDDEN,
        linal::VectorInterface< float> &INPUT_RECONSTRUCTED,
        linal::VectorInterface< float> &HIDDEN_RECONSTRUCTED,
        const size_t &STOCHASTIC_STEP_COUNT,
        linal::VectorInterface< float> &HIDDEN_SAMPLE,
        linal::VectorInterface< float> &VISIBLE_SAMPLE
      ) const;

      //! @brief Invert the boltzmann machine (swaps visible with hidden neurons)
      void Invert();

      //! @brief compute the reconstruction error RMSD across a dataset
      float ComputeReconstructionError( const FeatureDataSetInterface< float> &DATA) const;

      //! @brief compute the reconstruction error RMSD across a single data point
      float ComputeReconstructionError
      (
        const linal::VectorConstInterface< float> &INPUT,
        linal::VectorInterface< float> &HIDDEN,
        linal::VectorInterface< float> &INPUT_RECONSTRUCTED
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read RestrictedBoltzmannMachineLayer from std::istream
      std::istream &Read( std::istream &ISTREAM);

      //! write RestrictedBoltzmannMachineLayer into std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class RestrictedBoltzmannMachineLayer

  } // namespace model
} // namespace bcl

#endif //BCL_MODEL_RESTRICTED_BOLTZMANN_MACHINE_LAYER_H_
