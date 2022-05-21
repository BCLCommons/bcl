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
#include "model/bcl_model_neural_network_selective_backpropagation_tolerance.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_matrix_const_reference.h"
#include "math/bcl_math_running_average_sd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> NeuralNetworkSelectiveBackpropagationTolerance::s_Instance
    (
      util::Enumerated< NeuralNetworkSelectiveBackpropagationInterface>::AddInstance
      (
        new NeuralNetworkSelectiveBackpropagationTolerance()
      )
    );

    //! @brief default constructor
    NeuralNetworkSelectiveBackpropagationTolerance::NeuralNetworkSelectiveBackpropagationTolerance() :
      m_ResultStartID( 0),
      m_ResultEndID( 1000),
      m_Tolerance( 0.0),
      m_Noise( 0.0),
      m_ReduceErrorByTolerance( false),
      m_ResultsSize( 1000),
      m_NumberThreads( 0)
    {
    }

    //! @brief copy constructor
    //! @return a new NeuralNetworkSelectiveBackpropagationTolerance copied from this instance
    NeuralNetworkSelectiveBackpropagationTolerance *NeuralNetworkSelectiveBackpropagationTolerance::Clone() const
    {
      return new NeuralNetworkSelectiveBackpropagationTolerance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &NeuralNetworkSelectiveBackpropagationTolerance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &NeuralNetworkSelectiveBackpropagationTolerance::GetAlias() const
    {
      static const std::string s_name( "Tolerant");
      return s_name;
    }

    //! @brief Initialize this object with the rescaled dataset and other important parameters
    //! @param RESCALED_DATA the already-rescaled training dataset
    //! @param OBJECTIVE the actual objective function
    //! @param THREADS # of threads that may be accessing this object at once
    void NeuralNetworkSelectiveBackpropagationTolerance::Initialize
    (
      const descriptor::Dataset &RESCALED_DATA,
      const ObjectiveFunctionInterface &OBJECTIVE,
      const size_t &NUMBER_THREADS
    )
    {
      // initialization
      if( m_ResultEndID > RESCALED_DATA.GetResultSize())
      {
        m_ResultEndID = RESCALED_DATA.GetResultSize();
      }
      m_ResultsSize = m_ResultEndID - m_ResultStartID;
      m_NumberThreads = NUMBER_THREADS;

      // set up sizes of all internally-held members
      m_RescaledNoise = linal::Vector< float>( m_ResultsSize, m_Noise);

      // regression task: Compute standard deviation over the selected columns
      m_RescaledTolerance = linal::Vector< float>( m_ResultsSize, m_Tolerance);
      const util::SiPtr< const RescaleFeatureDataSet> results_scaling( RESCALED_DATA.GetResultsPtr()->GetScaling());
      const float rescale_to_range_width( results_scaling->GetRange().GetWidth());
      for( size_t result( 0); result < m_ResultsSize; ++result)
      {
        const float rescaled_range_ratio( results_scaling->GetRescaleRanges()( result).GetWidth() / rescale_to_range_width);
        m_RescaledTolerance( result) /= rescaled_range_ratio;
        m_RescaledNoise( result) /= rescaled_range_ratio;
      }
    } // Initialize

    //! @brief select whether or not to backpropagate the current feature, and edit the error vector, if desired
    //! @param PREDICTION the prediction that was made by the neural network
    //! @param ERROR Reference to the error vector, (should already have been set to RESULT - PREDICTION)
    //! @param FEATURE_ID id of this feature in the dataset
    //! @param THREAD_ID id of this thread
    //! @return true if the feature should be backpropagated
    bool NeuralNetworkSelectiveBackpropagationTolerance::ShouldBackpropagate
    (
      const linal::VectorConstInterface< float> &PREDICTION,
      linal::VectorInterface< float> &ERROR,
      const size_t &FEATURE_ID,
      const size_t &THREAD_ID
    )
    {
      bool backprop_should_continue( false);
      // regression target
      for( size_t result_offset( 0); result_offset < m_ResultsSize; ++result_offset)
      {
        // get the actual index
        float &er( ERROR( result_offset));
        const float tolerance( m_RescaledTolerance( result_offset));
        if( math::Absolute( er) <= tolerance)
        {
          er = 0.0;
        }
        else
        {
          backprop_should_continue = true;
          if( m_ReduceErrorByTolerance)
          {
            backprop_should_continue = true;
            if( er > 0.0)
            {
              er -= tolerance;
            }
            else
            {
              er += tolerance;
            }
          }
          // add noise, if desirable
          if( m_RescaledNoise( result_offset))
          {
            er += random::GetGlobalRandom().RandomGaussian( 0.0, m_RescaledNoise( result_offset));
          }
        }
      }

      // return true if any part of ERROR remains to be backpropagated
      return backprop_should_continue;
    }

    //! @brief finalize the current round; occurs only after all threads were already joined
    void NeuralNetworkSelectiveBackpropagationTolerance::FinalizeRound()
    {
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer NeuralNetworkSelectiveBackpropagationTolerance::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Always back-propagating features that are mispredicted beyond the tolerance. "
        "Can also add noise to the results"
      );
      parameters.AddInitializer
      (
        "begin",
        "result columns to use this triager for",
        io::Serialization::GetAgent( &m_ResultStartID),
        "0"
      );
      parameters.AddInitializer
      (
        "end",
        "1+max result column to consider",
        io::Serialization::GetAgent( &m_ResultEndID),
        "1000"
      );
      parameters.AddInitializer
      (
        "noise",
        "Std of noise to add to errors"
        "This is related to the uncertainty in the experimental values; though even with perfect experimental values, it "
        "should be non-zero for training general models",
        io::Serialization::GetAgent( &m_Noise),
        "0.0"
      );
      parameters.AddInitializer
      (
        "tolerance",
        "If < than this amount of error above the cutoff, do not backpropagate under any circumstances",
        io::Serialization::GetAgent( &m_Tolerance),
        "0"
      );
      parameters.AddInitializer
      (
        "reduce error by tolerance",
        "Whether to reduce the error for backpropagated features by the tolerance when the feature is backpropagated",
        io::Serialization::GetAgent( &m_ReduceErrorByTolerance),
        "True"
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
