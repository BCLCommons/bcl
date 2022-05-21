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
#include "model/bcl_model_train_restricted_boltzmann_machine_layer.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_const_reference.h"
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    //! the default input range for neural network transfer functions
    const math::Range< float> TrainRestrictedBoltzmannMachineLayer::s_DefaultInputRange( 0, 1);

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> TrainRestrictedBoltzmannMachineLayer::s_Instance
    (
      GetObjectInstances().AddInstance( new TrainRestrictedBoltzmannMachineLayer())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    TrainRestrictedBoltzmannMachineLayer::TrainRestrictedBoltzmannMachineLayer() :
      m_StochasticStepCount( 0),
      m_NumberFeaturesSeen( 0)
    {
    }

    //! @brief constructor from an untrained network
    //! @param LAYER the layer that this trainer will handle
    //! @param STOCHASTIC_STEP_COUNT helps control generalizability, see note above on m_StochasticStepCount
    TrainRestrictedBoltzmannMachineLayer::TrainRestrictedBoltzmannMachineLayer
    (
      RestrictedBoltzmannMachineLayer &LAYER,
      const size_t &STOCHASTIC_STEP_COUNT
    ) :
      m_Layer( LAYER),
      m_StochasticStepCount( STOCHASTIC_STEP_COUNT),
      m_WeightGradient( LAYER.GetNumberInputNeurons(), LAYER.GetNumberHiddenNeurons(), float( 0)),
      m_BiasVisibleGradient( LAYER.GetNumberInputNeurons(), float( 0)),
      m_BiasHiddenGradient( LAYER.GetNumberHiddenNeurons(), float( 0)),
      m_NoiseVisibleGradient( LAYER.GetNumberInputNeurons(), float( 0)),
      m_NoiseHiddenGradient( LAYER.GetNumberHiddenNeurons(), float( 0)),
      m_Hidden( LAYER.GetNumberHiddenNeurons(), float( 0)),
      m_VisibleReconstructed( LAYER.GetNumberInputNeurons(), float( 0)),
      m_HiddenReconstructed( LAYER.GetNumberHiddenNeurons(), float( 0)),
      m_VisibleSample( LAYER.GetNumberInputNeurons(), float( 0)),
      m_HiddenSample( LAYER.GetNumberHiddenNeurons(), float( 0)),
      m_NumberFeaturesSeen( 0)
    {
    }

    //! copy constructor
    TrainRestrictedBoltzmannMachineLayer *TrainRestrictedBoltzmannMachineLayer::Clone() const
    {
      return new TrainRestrictedBoltzmannMachineLayer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TrainRestrictedBoltzmannMachineLayer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief core function : takes a feature and uses it to update the associated gradients
    //! @param INPUT the input feature of interest
    //! @return reference to the hidden layer output
    linal::VectorConstReference< float> TrainRestrictedBoltzmannMachineLayer::Train( const linal::VectorConstInterface< float> &INPUT)
    {
      BCL_Assert
      (
        INPUT.GetSize() == m_Layer->GetNumberInputNeurons(),
        "INPUT.GetSize() does not match the size of the layer!"
      );

      const float *const itr_bias_vg_end( m_BiasVisibleGradient.End());
      const float *const itr_bias_hg_end( m_BiasHiddenGradient.End());
      const float *const itr_weight_gradient_end( m_WeightGradient.End());
      const float *const itr_hid_recon_end( m_HiddenReconstructed.End());

      // perform 1 round of Gibbs sampling
      m_Layer->GibbsSample
      (
        INPUT,
        m_Hidden,
        m_VisibleReconstructed,
        m_HiddenReconstructed,
        m_StochasticStepCount,
        m_HiddenSample,
        m_VisibleSample
      );

      // Update gradients
      // Update Visible bias gradient with FEATURE - visible_reconstructed
      {
        const float *itr_vis( INPUT.Begin());
        for
        (
          float *itr_bias_vg( m_BiasVisibleGradient.Begin()), *itr_vis_recon( m_VisibleReconstructed.Begin());
          itr_bias_vg != itr_bias_vg_end;
          ++itr_vis, ++itr_bias_vg, ++itr_vis_recon
        )
        {
          *itr_bias_vg += *itr_vis - *itr_vis_recon;
        }
      }
      {
        const float *itr_hidden( m_Hidden.Begin());
        // Update Hidden bias gradient hidden - hidden_reconstructed
        for
        (
          float *itr_bias_hg( m_BiasHiddenGradient.Begin()), *itr_hid_recon( m_HiddenReconstructed.Begin());
          itr_bias_hg != itr_bias_hg_end;
          ++itr_hidden, ++itr_bias_hg, ++itr_hid_recon
        )
        {
          *itr_bias_hg += *itr_hidden - *itr_hid_recon;
        }
      }
      // noise update
      if( m_Layer->m_NetworkType == RestrictedBoltzmannMachineLayer::e_StochasticSigmoid)
      {
        const float *itr_vis_recon( m_VisibleReconstructed.Begin()), *itr_vis_recon_end( m_VisibleReconstructed.End());
        const float *itr_vis( INPUT.Begin());
        for
        (
          float *itr_grad_noise( m_NoiseVisibleGradient.Begin());
          itr_vis_recon != itr_vis_recon_end;
          ++itr_vis, ++itr_vis_recon, ++itr_grad_noise
        )
        {
          *itr_grad_noise += math::Sqr( *itr_vis) - math::Sqr( *itr_vis_recon);
        }
        const float *itr_hid_recon( m_HiddenReconstructed.Begin()), *itr_hidden( m_Hidden.Begin());
        for
        (
          float *itr_grad_noise( m_NoiseHiddenGradient.Begin());
          itr_hid_recon != itr_hid_recon_end;
          ++itr_hidden, ++itr_hid_recon, ++itr_grad_noise
        )
        {
          *itr_grad_noise += math::Sqr( *itr_hidden) - math::Sqr( *itr_hid_recon);
        }
      }
      {
        // Update weight gradient with hidden' * FEATURE - hidden_reconstructed' * visible_reconstructed;
        float *itr_weight_gradient( m_WeightGradient.Begin());
        for
        (
          const float *itr_visible( INPUT.Begin()), *itr_visible_recon( m_VisibleReconstructed.Begin());
          itr_weight_gradient != itr_weight_gradient_end;
          ++itr_visible, ++itr_visible_recon
        )
        {
          const float *itr_hidden( m_Hidden.Begin());
          for
          (
            const float *itr_hid_recon( m_HiddenReconstructed.Begin());
            itr_hid_recon != itr_hid_recon_end;
            ++itr_hid_recon, ++itr_hidden, ++itr_weight_gradient
          )
          {
            *itr_weight_gradient += *itr_hidden * *itr_visible - *itr_hid_recon * *itr_visible_recon;
          }
        }
      }
      ++m_NumberFeaturesSeen;
      return linal::VectorConstReference< float>( m_Hidden);
    }

    //! @brief core function : takes a feature and a result and uses it to update the associated gradients
    //! @param INPUT the input feature of interest
    //! @param OUTPUT the output feature of interest
    //! @return reference to reconstructed output
    linal::VectorConstReference< float> TrainRestrictedBoltzmannMachineLayer::Train
    (
      const linal::VectorConstInterface< float> &INPUT,
      const linal::VectorConstInterface< float> &OUTPUT
    )
    {
      BCL_Assert
      (
        INPUT.GetSize() + OUTPUT.GetSize() == m_Layer->GetNumberInputNeurons(),
        "INPUT.GetSize() + OUTPUT.GetSize() does not match the size of the layer!"
      );

      const float *const itr_bias_hg_end( m_BiasHiddenGradient.End());
      const float *const itr_hid_recon_end( m_HiddenReconstructed.End());

      // perform 1 round of Gibbs sampling
      m_Layer->GibbsSample
      (
        INPUT,
        m_Hidden,
        m_VisibleReconstructed,
        m_HiddenReconstructed,
        m_StochasticStepCount,
        m_HiddenSample,
        m_VisibleSample
      );

      // Update gradients
      // Update Visible bias gradient with FEATURE - visible_reconstructed
      {
        float *itr_bias_vg( m_BiasVisibleGradient.Begin()), *itr_vis_recon( m_VisibleReconstructed.Begin());
        for
        (
          const float *itr_vis( INPUT.Begin()), *itr_vis_end( INPUT.End());
          itr_vis != itr_vis_end;
          ++itr_vis, ++itr_bias_vg, ++itr_vis_recon
        )
        {
          *itr_bias_vg += *itr_vis - *itr_vis_recon;
        }
        for
        (
          const float *itr_out( OUTPUT.Begin()), *itr_out_end( OUTPUT.End());
          itr_out != itr_out_end;
          ++itr_out, ++itr_bias_vg, ++itr_vis_recon
        )
        {
          *itr_bias_vg += *itr_out - *itr_vis_recon;
        }
      }
      {
        const float *itr_hidden( m_Hidden.Begin());
        // Update Hidden bias gradient hidden - hidden_reconstructed
        for
        (
          float *itr_bias_hg( m_BiasHiddenGradient.Begin()), *itr_hid_recon( m_HiddenReconstructed.Begin());
          itr_bias_hg != itr_bias_hg_end;
          ++itr_hidden, ++itr_bias_hg, ++itr_hid_recon
        )
        {
          *itr_bias_hg += *itr_hidden - *itr_hid_recon;
        }
      }
      // noise update
      if( m_Layer->m_NetworkType == RestrictedBoltzmannMachineLayer::e_StochasticSigmoid)
      {
        {
          const float *itr_vis_recon( m_VisibleReconstructed.Begin());
          float *itr_grad_noise( m_NoiseVisibleGradient.Begin());
          for
          (
            const float *itr_vis( INPUT.Begin()), *itr_vis_end( INPUT.End());
            itr_vis != itr_vis_end;
            ++itr_vis, ++itr_vis_recon, ++itr_grad_noise
          )
          {
            *itr_grad_noise += math::Sqr( *itr_vis) - math::Sqr( *itr_vis_recon);
          }
          for
          (
            const float *itr_out( OUTPUT.Begin()), *itr_out_end( OUTPUT.End());
            itr_out != itr_out_end;
            ++itr_out, ++itr_vis_recon, ++itr_grad_noise
          )
          {
            *itr_grad_noise += math::Sqr( *itr_out) - math::Sqr( *itr_vis_recon);
          }
        }
        const float *itr_hid_recon( m_HiddenReconstructed.Begin()), *itr_hidden( m_Hidden.Begin());
        for
        (
          float *itr_grad_noise( m_NoiseHiddenGradient.Begin());
          itr_hid_recon != itr_hid_recon_end;
          ++itr_hidden, ++itr_hid_recon, ++itr_grad_noise
        )
        {
          *itr_grad_noise += math::Sqr( *itr_hidden) - math::Sqr( *itr_hid_recon);
        }
      }
      {
        // Update weight gradient with hidden' * FEATURE - hidden_reconstructed' * visible_reconstructed;
        float *itr_weight_gradient( m_WeightGradient.Begin());
        const float *itr_visible_recon( m_VisibleReconstructed.Begin());
        for
        (
          const float *itr_visible( INPUT.Begin()), *itr_visible_end( INPUT.End());
          itr_visible != itr_visible_end;
          ++itr_visible, ++itr_visible_recon
        )
        {
          const float *itr_hidden( m_Hidden.Begin());
          for
          (
            const float *itr_hid_recon( m_HiddenReconstructed.Begin());
            itr_hid_recon != itr_hid_recon_end;
            ++itr_hid_recon, ++itr_hidden, ++itr_weight_gradient
          )
          {
            *itr_weight_gradient += *itr_hidden * *itr_visible - *itr_hid_recon * *itr_visible_recon;
          }
        }
        for
        (
          const float *itr_out( OUTPUT.Begin()), *itr_out_end( OUTPUT.End());
          itr_out != itr_out_end;
          ++itr_out, ++itr_visible_recon
        )
        {
          const float *itr_hidden( m_Hidden.Begin());
          for
          (
            const float *itr_hid_recon( m_HiddenReconstructed.Begin());
            itr_hid_recon != itr_hid_recon_end;
            ++itr_hid_recon, ++itr_hidden, ++itr_weight_gradient
          )
          {
            *itr_weight_gradient += *itr_hidden * *itr_out - *itr_hid_recon * *itr_visible_recon;
          }
        }
      }
      ++m_NumberFeaturesSeen;

      return linal::VectorConstReference< float>( OUTPUT.GetSize(), m_VisibleReconstructed.Begin() + INPUT.GetSize());
    }

    //! @brief Update layer to update the associated RestrictedBoltzmannMachineLayer's parameters
    //! @param WEIGHT_COST relative cost of weight
    //! @param WEIGHT_UPDATE update object for the weight
    //! @param VISIBLE_BIAS_UPDATE update object for the visible layer's bias
    //! @param HIDDEN_BIAS_UPDATE update object for the hidden layer's bias
    //! @param VISIBLE_NOISE_UPDATE update object for the visible layer's noise (for StochasticSigmoid RBMs)
    //! @param HIDDEN_NOISE_UPDATE update object for the hidden layer's noise (for StochasticSigmoid RBMs)
    void TrainRestrictedBoltzmannMachineLayer::UpdateLayer
    (
      const float &WEIGHT_COST,
      NeuralNetworkUpdateWeightsInterface &WEIGHT_UPDATE,
      NeuralNetworkUpdateWeightsInterface &VISIBLE_BIAS_UPDATE,
      NeuralNetworkUpdateWeightsInterface &HIDDEN_BIAS_UPDATE,
      NeuralNetworkUpdateWeightsInterface &VISIBLE_NOISE_UPDATE,
      NeuralNetworkUpdateWeightsInterface &HIDDEN_NOISE_UPDATE
    )
    {
      // divide all weight vectors by the number of feature seen
      m_BiasHiddenGradient /= float( m_NumberFeaturesSeen);
      m_BiasVisibleGradient /= float( m_NumberFeaturesSeen);
      m_NoiseHiddenGradient /= float( m_NumberFeaturesSeen);
      m_NoiseVisibleGradient /= float( m_NumberFeaturesSeen);
      m_WeightGradient /= float( m_NumberFeaturesSeen);

      const float *const itr_weight_gradient_end( m_WeightGradient.End());
      // weight decay and update of change weight
      // This code tries to just do an efficient version of the following:
      // weight_gradient -= WEIGHT_COST * m_Weight;
      // m_Weight += change_weight_gradient;
      for
      (
        float *itr_weight_grad( m_WeightGradient.Begin()), *itr_weight( m_Layer->m_Weight.Begin());
        itr_weight_grad != itr_weight_gradient_end;
        ++itr_weight_grad, ++itr_weight
      )
      {
        // weight decay on the gradient
        *itr_weight_grad -= *itr_weight * WEIGHT_COST;
      }
      // update weights
      WEIGHT_UPDATE( m_Layer->m_Weight.Begin(), m_WeightGradient.Begin());

      // update bias
      VISIBLE_BIAS_UPDATE( m_Layer->m_BiasVisible.Begin(), m_BiasVisibleGradient.Begin());
      HIDDEN_BIAS_UPDATE( m_Layer->m_BiasHidden.Begin(), m_BiasHiddenGradient.Begin());

      // update noise for stochastic sigmoidal network
      if( m_Layer->m_NetworkType == RestrictedBoltzmannMachineLayer::e_StochasticSigmoid)
      {
        m_NoiseVisibleGradient /= m_Layer->m_SlopeVisible;
        m_NoiseVisibleGradient /= m_Layer->m_SlopeVisible;
        m_NoiseHiddenGradient /= m_Layer->m_SlopeHidden;
        m_NoiseHiddenGradient /= m_Layer->m_SlopeHidden;

        VISIBLE_NOISE_UPDATE( m_Layer->m_SlopeVisible.Begin(), m_NoiseVisibleGradient.Begin());
        HIDDEN_NOISE_UPDATE( m_Layer->m_SlopeHidden.Begin(), m_NoiseHiddenGradient.Begin());
      }
    }

    //! @brief Reset this object, to prepare for the next batch of data
    void TrainRestrictedBoltzmannMachineLayer::Reset()
    {
      m_WeightGradient.SetZero();
      m_BiasVisibleGradient = float( 0.0);
      m_NoiseVisibleGradient = float( 0.0);
      m_BiasHiddenGradient = float( 0.0);
      m_BiasVisibleGradient = float( 0.0);
      m_NumberFeaturesSeen = 0;
    }

    //! @brief Accumulate changes from a different trainer of the same layer (used for reduction after thread-completion)
    //! @param OTHER the other trainer's whose changes should be added to this layers
    void TrainRestrictedBoltzmannMachineLayer::AccumulateChangesFrom( const TrainRestrictedBoltzmannMachineLayer &OTHER)
    {
      m_WeightGradient += OTHER.m_WeightGradient;
      m_BiasVisibleGradient += OTHER.m_BiasVisibleGradient;
      m_BiasHiddenGradient += OTHER.m_BiasHiddenGradient;
      m_NoiseVisibleGradient += OTHER.m_NoiseVisibleGradient;
      m_NoiseHiddenGradient += OTHER.m_NoiseHiddenGradient;
      m_NumberFeaturesSeen += OTHER.m_NumberFeaturesSeen;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read TrainRestrictedBoltzmannMachineLayer from std::istream
    std::istream &TrainRestrictedBoltzmannMachineLayer::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( *m_Layer, ISTREAM);
      io::Serialize::Read( m_StochasticStepCount, ISTREAM);
      io::Serialize::Read( m_WeightGradient, ISTREAM);
      io::Serialize::Read( m_BiasVisibleGradient, ISTREAM);
      io::Serialize::Read( m_BiasHiddenGradient, ISTREAM);
      io::Serialize::Read( m_NoiseVisibleGradient, ISTREAM);
      io::Serialize::Read( m_NoiseHiddenGradient, ISTREAM);
      io::Serialize::Read( m_NumberFeaturesSeen, ISTREAM);

      // reset the members that do not store state
      m_HiddenSample = m_HiddenReconstructed = m_Hidden = m_BiasVisibleGradient;
      m_VisibleReconstructed = m_VisibleSample = m_BiasVisibleGradient;
      return ISTREAM;
    }

    //! write TrainRestrictedBoltzmannMachineLayer into std::ostream
    std::ostream &TrainRestrictedBoltzmannMachineLayer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( *m_Layer, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_StochasticStepCount, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_WeightGradient, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BiasVisibleGradient, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_BiasHiddenGradient, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NoiseVisibleGradient, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NoiseHiddenGradient, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NumberFeaturesSeen, OSTREAM, INDENT);
      return OSTREAM;
    }

  } // namespace model
} // namespace bcl
