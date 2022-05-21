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
#include "model/bcl_model_restricted_boltzmann_machine_layer.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_matrix_reference.h"
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    //! the default input range for neural network transfer functions
    const math::Range< float> RestrictedBoltzmannMachineLayer::s_DefaultInputRange( 0, 1);

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> RestrictedBoltzmannMachineLayer::s_Instance
    (
      GetObjectInstances().AddInstance( new RestrictedBoltzmannMachineLayer())
    );

    //! @brief Type as string
    //! @param TYPE the type
    //! @return the string for the kernel
    const std::string &RestrictedBoltzmannMachineLayer::GetTypeName( const Type &TYPE)
    {
      static const std::string s_names[] =
      {
        "StochasticSigmoid",
        "Binomial",
        GetStaticClassName< Type>()
      };
      return s_names[ TYPE];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    RestrictedBoltzmannMachineLayer::RestrictedBoltzmannMachineLayer() :
      m_NumberVisibleNodes( 0),
      m_NumberHiddenNodes( 0),
      m_NetworkType( e_StochasticSigmoid)
    {
    }

    //! @brief construct an untrained network
    //! @param NUMBER_VISIBLE number of visible neurons
    //! @param NUMBER_HIDDEN number of hidden neurons
    //! @param NETWORK_TYPE type of network in terms of transfer functions
    RestrictedBoltzmannMachineLayer::RestrictedBoltzmannMachineLayer
    (
      const size_t &NUMBER_VISIBLE,
      const size_t &NUMBER_HIDDEN,
      const Type   &NETWORK_TYPE
    ) :
      m_NumberVisibleNodes( 0),
      m_NumberHiddenNodes( 0),
      m_NetworkType( NETWORK_TYPE)
    {
      SetArchitecture( NUMBER_VISIBLE, NUMBER_HIDDEN);
    }

    //! @brief construct from all necessary parameters
    //! @param BIAS_VISIBLE the bias for each neuron of the visible layer
    //! @param BIAS_HIDDEN the bias for each neuron of the hidden layer
    //! @param WEIGHT the weight for each neuron of a hidden layer
    //! @param TRANSFER_FUNCTION_HIDDEN transfer function for the hidden layer neurons
    //! @param TRANSFER_FUNCTION_VISIBLE transfer function for the hidden layer neurons
    //! @param HIDDEN_NOISE_VARIANCE variance of the noise to add to the hidden layer when propagating to the visible
    RestrictedBoltzmannMachineLayer::RestrictedBoltzmannMachineLayer
    (
      const linal::Vector< float> &BIAS_VISIBLE,
      const linal::Vector< float> &BIAS_HIDDEN,
      const linal::Matrix< float> &WEIGHT,
      const Type                  &NETWORK_TYPE
    ) :
      m_BiasHidden( BIAS_HIDDEN),
      m_BiasVisible( BIAS_VISIBLE),
      m_Weight( WEIGHT),
      m_NumberVisibleNodes( BIAS_VISIBLE.GetSize()),
      m_NumberHiddenNodes( BIAS_HIDDEN.GetSize()),
      m_NetworkType( NETWORK_TYPE)
    {
      m_SlopeHidden = linal::Vector< float>( m_NumberHiddenNodes, float( 1.0));
      m_SlopeVisible = linal::Vector< float>( m_NumberVisibleNodes, float( 1.0));
      BCL_Assert
      (
        m_NumberHiddenNodes == m_Weight.GetNumberCols(),
        "Incorrectly sized weight matrix or bias vector for hidden neurons"
      );
      BCL_Assert
      (
        m_NumberVisibleNodes == m_Weight.GetNumberRows(),
        "Incorrectly sized weight matrix or bias vector for hidden neurons"
      );
    }

    //! copy constructor
    RestrictedBoltzmannMachineLayer *RestrictedBoltzmannMachineLayer::Clone() const
    {
      return new RestrictedBoltzmannMachineLayer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RestrictedBoltzmannMachineLayer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the architecture of this RBM
    //! @param NUMBER_INPUTS number of input/visible neurons
    //! @param NUMBER_OUTPUTS number of output/hidden neurons
    void RestrictedBoltzmannMachineLayer::SetArchitecture( const size_t &NUMBER_INPUTS, const size_t &NUMBER_OUTPUTS)
    {
      if( m_NumberHiddenNodes == NUMBER_OUTPUTS || m_NumberVisibleNodes == NUMBER_INPUTS)
      {
        return;
      }
      m_NumberHiddenNodes = NUMBER_OUTPUTS;
      m_NumberVisibleNodes = NUMBER_INPUTS;
      m_BiasHidden = linal::Vector< float>( m_NumberHiddenNodes, float( 0.0));
      m_BiasVisible = linal::Vector< float>( m_NumberVisibleNodes, float( 0.0));
      m_SlopeHidden = linal::Vector< float>( m_NumberHiddenNodes, float( 1.0));
      m_SlopeVisible = linal::Vector< float>( m_NumberVisibleNodes, float( 1.0));
      m_Weight = linal::Matrix< float>( m_NumberVisibleNodes, m_NumberHiddenNodes, float( 0.0));
      for( float *itr( m_Weight.Begin()), *itr_end( m_Weight.End()); itr != itr_end; ++itr)
      {
        *itr = random::GetGlobalRandom().RandomGaussian( 0.0, 0.1);
      }
    }

    //! @brief stochastic activation/ sampling function for continuous distributions
    //! @see @link http://www.ee.nthu.edu.tw/~hchen/pubs/iee2003.pdf @endlink
    //! @param X the activation variable
    //! @param A the noise sensitivity parameter
    float StochasticActivation( const float &X, const float &A)
    {
      // Note: This function differs slightly from that given in http://www.ee.nthu.edu.tw/~hchen/pubs/iee2003.pdf
      // (in particular, the A also multiplies X in that paper), however, this version is both more stable and
      // mathematically justified, since in the paper's equation, A is effectively modifying all the weights and biases
      // by a multiplicative factor, which is never accounted for in the update equations. In practice, all we need is a
      // noise tuning parameter that prevents the RBM from over-fitting
      return 1.0 / ( 1.0 + exp( -random::GetGlobalRandom().RandomGaussian( X, 0.25 * A)));
    }

    //! @brief activate a hidden or visible set of neurons
    //! @param NEURONS the neurons to update
    //! @param SLOPES activation slopes for the neurons
    void Sample
    (
      linal::VectorInterface< float> &NEURONS,
      const linal::VectorConstInterface< float> &NEURON_ACTIVITY
    )
    {
      for( size_t i( 0), number_neurons( NEURONS.GetSize()); i < number_neurons; ++i)
      {
        NEURONS( i) = random::GetGlobalRandom().Double() < NEURON_ACTIVITY( i);
      }
    }

    //! @brief A simple sigmoid function
    //! @param X the activation variable
    float Sigmoid( const float &X)
    {
      return 1.0 / ( 1.0 + exp( -X));
    }

    //! @brief activate a hidden or visible set of neurons
    //! @param NEURONS the neurons to update
    //! @param SLOPES activation slopes for the neurons
    //! @param NETWORK_TYPE actual type of the network
    //! @param STOCHASTIC whether to apply stochastic weighting onto the hidden layer
    void ActivateAndSample
    (
      linal::VectorInterface< float> &NEURONS,
      const linal::Vector< float> &SLOPES,
      const RestrictedBoltzmannMachineLayer::Type TYPE,
      const bool STOCHASTIC
    )
    {
      const size_t number_neurons( NEURONS.GetSize());
      if( STOCHASTIC && RestrictedBoltzmannMachineLayer::e_StochasticSigmoid)
      {
        for( size_t i( 0); i < number_neurons; ++i)
        {
          NEURONS( i) = StochasticActivation( NEURONS( i), SLOPES( i));
        }
      }
      else
      {
        for( size_t i( 0); i < number_neurons; ++i)
        {
          NEURONS( i) = Sigmoid( NEURONS( i));
        }
      }
    }

  ////////////////
  // operations //
  ////////////////

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
    void RestrictedBoltzmannMachineLayer::GibbsSample
    (
      const linal::VectorConstInterface< float> &INPUT,
      linal::VectorInterface< float> &HIDDEN,
      linal::VectorInterface< float> &INPUT_RECONSTRUCTED,
      linal::VectorInterface< float> &HIDDEN_RECONSTRUCTED,
      const size_t &STOCHASTIC_STEP_COUNT,
      linal::VectorInterface< float> &HIDDEN_SAMPLE,
      linal::VectorInterface< float> &VISIBLE_SAMPLE
    ) const
    {
      if( !STOCHASTIC_STEP_COUNT)
      {
        // simple non-stochastic case
        ForwardPropagate( INPUT, HIDDEN, false);
        BackPropagate( INPUT_RECONSTRUCTED, HIDDEN, false);
        ForwardPropagate( INPUT_RECONSTRUCTED, HIDDEN_RECONSTRUCTED, false);
      }
      else if( m_NetworkType == e_StochasticSigmoid)
      {
        // stochastic sigmoid does not use the sampling vectors since sampling and activation are performed in a single
        // step
        ForwardPropagate( INPUT, HIDDEN, true);
        BackPropagate( INPUT_RECONSTRUCTED, HIDDEN, STOCHASTIC_STEP_COUNT > 1);
        ForwardPropagate( INPUT_RECONSTRUCTED, HIDDEN_RECONSTRUCTED, STOCHASTIC_STEP_COUNT > 2);
      }
      else
      {
        // binomial activation; sampling must occur separately from activation to compute the proper derivatives
        ForwardPropagate( INPUT, HIDDEN, false);
        Sample( HIDDEN_SAMPLE, HIDDEN);
        BackPropagate( INPUT_RECONSTRUCTED, HIDDEN_SAMPLE, false);
        if( STOCHASTIC_STEP_COUNT > 1)
        {
          Sample( VISIBLE_SAMPLE, INPUT_RECONSTRUCTED);
          ForwardPropagate( VISIBLE_SAMPLE, HIDDEN_RECONSTRUCTED, false);
        }
        else
        {
          ForwardPropagate( INPUT_RECONSTRUCTED, HIDDEN_RECONSTRUCTED, false);
        }
      }
    }

    //! @brief Forward propagate a given feature from the input to the hidden layer
    //! @param INPUT input feature of interest
    //! @param HIDDEN reference to vector that will be overwritten by the hidden layer
    //! @param STOCHASTIC whether to apply stochastic weighting onto the hidden layer
    void RestrictedBoltzmannMachineLayer::ForwardPropagate
    (
      const linal::VectorConstInterface< float> &INPUT,
      linal::VectorInterface< float> &HIDDEN,
      const bool STOCHASTIC
    ) const
    {
      // Unrolled version of
      // HIDDEN = m_Weight * INPUT + m_BiasHidden
      HIDDEN = m_BiasHidden;
      const float *itr_weight( m_Weight.Begin()), *input( INPUT.Begin()), *itr_hidden_end( HIDDEN.End());
      for( const float *itr_input_end( INPUT.End()); input != itr_input_end; ++input)
      {
        const float input_val( *input);
        for( float *itr_hidden( HIDDEN.Begin()); itr_hidden != itr_hidden_end; ++itr_hidden, ++itr_weight)
        {
          *itr_hidden += *itr_weight * input_val;
        }
      }
      ActivateAndSample( HIDDEN, m_SlopeHidden, m_NetworkType, STOCHASTIC);
    }

    //! @brief Back propagate a given feature from the hidden to the input layer
    //! @param INPUT input feature of interest
    //! @param HIDDEN hidden vector
    //! @param STOCHASTIC whether to apply stochastic weighting onto the input layer
    void RestrictedBoltzmannMachineLayer::BackPropagate
    (
      linal::VectorInterface< float> &INPUT,
      const linal::VectorInterface< float> &HIDDEN,
      const bool STOCHASTIC
    ) const
    {
      INPUT = m_BiasVisible;

      const float *itr_hidden_end( HIDDEN.End());
      const float *itr_weight( m_Weight.Begin());

      for( float *itr_input( INPUT.Begin()), *itr_end( INPUT.End()); itr_input != itr_end; ++itr_input)
      {
        float sum( 0);
        for( const float *itr_hidden( HIDDEN.Begin()); itr_hidden != itr_hidden_end; ++itr_hidden, ++itr_weight)
        {
          sum += *itr_hidden * *itr_weight;
        }
        *itr_input += sum;
      }
      ActivateAndSample( INPUT, m_SlopeVisible, m_NetworkType, STOCHASTIC);
    }

    //! @brief Invert the boltzmann machine (swaps visible with hidden neurons)
    void RestrictedBoltzmannMachineLayer::Invert()
    {
      std::swap( m_SlopeHidden, m_SlopeVisible);
      std::swap( m_BiasHidden, m_BiasVisible);
      std::swap( m_NumberVisibleNodes, m_NumberHiddenNodes);
      m_Weight = m_Weight.Transposed();
    }

    //! @brief compute the reconstruction error RMSD across a dataset
    float RestrictedBoltzmannMachineLayer::ComputeReconstructionError( const FeatureDataSetInterface< float> &DATA) const
    {
      math::RunningAverage< float> ave_deviation;
      // initialize vectors to store the reconstructed values at the visible layer
      linal::Vector< float> visible_reconstructed( m_NumberVisibleNodes), hidden( m_NumberHiddenNodes);
      const float *const itr_vis_recon_end( visible_reconstructed.End());
      for( size_t feature_id( 0), dataset_size( DATA.GetNumberFeatures()); feature_id < dataset_size; ++feature_id)
      {
        // get a reference to the actual feature
        const linal::VectorConstInterface< float> &feature( DATA( feature_id));

        // Determine non-stochastic reconstruction error
        ForwardPropagate( feature, hidden, false);
        BackPropagate( visible_reconstructed, hidden, false);

        // Update gradients
        // Update Visible bias gradient with FEATURE - visible_reconstructed
        for
        (
          const float *itr_vis_recon( visible_reconstructed.Begin()), *itr_vis( feature.Begin());
          itr_vis_recon != itr_vis_recon_end;
          ++itr_vis, ++itr_vis_recon
        )
        {
          ave_deviation += math::Sqr( *itr_vis - *itr_vis_recon);
        }
      }
      return math::Sqrt( ave_deviation.GetAverage());
    }

    //! @brief compute the reconstruction error RMSD across a single data point
    float RestrictedBoltzmannMachineLayer::ComputeReconstructionError
    (
      const linal::VectorConstInterface< float> &INPUT,
      linal::VectorInterface< float> &HIDDEN,
      linal::VectorInterface< float> &INPUT_RECONSTRUCTED
    ) const
    {
      // Determine non-stochastic reconstruction error
      ForwardPropagate( INPUT, HIDDEN, false);
      BackPropagate( INPUT_RECONSTRUCTED, HIDDEN, false);

      float deviation( 0);
      for
      (
        const float *itr_vis_recon( INPUT_RECONSTRUCTED.Begin()), *itr_vis( INPUT.Begin()), *itr_vis_end( INPUT.End());
        itr_vis != itr_vis_end;
        ++itr_vis, ++itr_vis_recon
      )
      {
        deviation += math::Sqr( *itr_vis - *itr_vis_recon);
      }
      deviation /= float( INPUT.GetSize());
      return math::Sqrt( deviation);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read RestrictedBoltzmannMachineLayer from std::istream
    std::istream &RestrictedBoltzmannMachineLayer::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_SlopeHidden, ISTREAM);
      io::Serialize::Read( m_SlopeVisible, ISTREAM);
      io::Serialize::Read( m_BiasHidden, ISTREAM);
      io::Serialize::Read( m_BiasVisible, ISTREAM);
      io::Serialize::Read( m_Weight, ISTREAM);
      io::Serialize::Read( m_NetworkType, ISTREAM);
      m_NumberVisibleNodes = m_BiasVisible.GetSize();
      m_NumberHiddenNodes = m_BiasHidden.GetSize();
      return ISTREAM;
    }

    //! write RestrictedBoltzmannMachineLayer into std::ostream
    std::ostream &RestrictedBoltzmannMachineLayer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_SlopeHidden, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SlopeVisible, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BiasHidden, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BiasVisible, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Weight, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NetworkType, OSTREAM, INDENT);
      return OSTREAM;
    }

  } // namespace model
} // namespace bcl
