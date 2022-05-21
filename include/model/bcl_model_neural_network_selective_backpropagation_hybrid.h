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

#ifndef BCL_MODEL_NEURAL_NETWORK_SELECTIVE_BACKPROPAGATION_HYBRID_H_
#define BCL_MODEL_NEURAL_NETWORK_SELECTIVE_BACKPROPAGATION_HYBRID_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_neural_network_selective_backpropagation_interface.h"
#include "linal/bcl_linal_matrix_const_reference.h"
#include "math/bcl_math_linear_function.h"
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class NeuralNetworkSelectiveBackpropagationHybrid
    //! @brief Balances # of features backpropagated on each side of the cutoff; preferentially selects features that
    //!        are most sensitive and nearest the cutoff
    //!
    //! @see @link example_model_neural_network_selective_backpropagation_hybrid.cpp @endlink
    //! @author mendenjl
    //! @date Aug 05, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API NeuralNetworkSelectiveBackpropagationHybrid :
      public NeuralNetworkSelectiveBackpropagationInterface
    {
    private:

    //////////
    // data //
    //////////

      //! Parameters
      bool   m_PureClassification;   //!< Set for pure classification tasks, e.g. to true to pretend all actives are at the upper limit, and all inactives at the lower
      size_t m_ResultStartID;     //!< Min result to consider
      size_t m_ResultEndID;       //!< 1+Max result to consider
      float  m_EnrichmentCutoff;  //!< Cutoff for enrichment optimization
      float  m_StabilityRatio;    //!< Ratio of FP+FN backpropagated due to incorrectness to those backproppd due to instability
      float  m_SensitivityWeight; //!< Weighting of false positive/negative function vs. sensitivity
      float  m_FalsePositiveBonus; //!< Bonus for being a false positive
      float  m_FalseNegativeBonus; //!< Bonus for being a false negative
      bool   m_StartWithAll;       //!< True to backpropagate all features in the first round

      //! Constant data (after Initialize is called)
      size_t                m_ResultsSize;             //!< Results size, cached
      size_t                m_DatasetSize;             //!< Dataset size, cached
      size_t                m_NumberThreads;           //!< Number of threads
      bool                  m_Parity;                  //!< True if using enrichment cutoff and actives are above the cutoff
      linal::Vector< float> m_ScaledCutoffs;           //!< Rescaled cutoffs, one for each result
      linal::Vector< float> m_EffectiveCutoffs;        //!< Effective cutoff (= scaled cutoff, unless pure classification, then 0.5)
      linal::MatrixConstReference< float> m_Results;   //!< Results dataset
      linal::Vector< size_t> m_AboveCutoffCount;       //!< # of features (for each result) at or above the cutoff
      linal::Vector< size_t> m_BelowCutoffCount;       //!< # of features (for each result) below the cutoff

      //! Round - dependent variables (constant during a round, updated during FinalizeRound)
      //! Error multipliers for high/low values
      math::RunningAverage< linal::Matrix< float> > m_ChangeAverager; //!< Average change, weighted by round #
      linal::Matrix< float>  m_AverageChange;
      linal::Vector< float>  m_MaxDelta;
      size_t                 m_RoundNumber;
      linal::Vector< float>  m_ThresholdBPHigh;
      linal::Vector< float>  m_ThresholdBPLow;
      linal::Vector< float>  m_SensitivityScalingHigh; //!< Average change for each result / average distance to cutoff
      linal::Vector< float>  m_SensitivityScalingLow; //!< Average change for each result / average distance to cutoff

      //! Matrices changed during round, but which are not thread dependent
      linal::Matrix< float>   m_PredictedValue;            //!< Stores predicted values
      linal::Matrix< float>   m_ChangeLastRound;           //!< Delta from last round
      linal::Matrix< float>   m_BackPropPreference;        //!< Backpropagation probabilities for all features

      //! Thread-local variables
      //! # of incorrect results computed by each thread last turn
      storage::Vector< linal::Vector< size_t> > m_IncorrectTallyAboveCutoff;
      storage::Vector< linal::Vector< size_t> > m_IncorrectTallyBelowCutoff;
      storage::Vector< linal::Vector< size_t> > m_BackpropagatedTallyAboveCutoff;
      storage::Vector< linal::Vector< size_t> > m_BackpropagatedTallyBelowCutoff;
      storage::Vector< storage::Vector< float> > m_ThreadMaxDelta;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      NeuralNetworkSelectiveBackpropagationHybrid();

      //! @brief copy constructor
      //! @return a new NeuralNetworkSelectiveBackpropagationHybrid copied from this instance
      NeuralNetworkSelectiveBackpropagationHybrid *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Initialize this object with the rescaled dataset and other important parameters
      //! @param RESCALED_DATA the already-rescaled training dataset
      //! @param OBJECTIVE the actual objective function
      //! @param THREADS # of threads that may be accessing this object at once
      void Initialize
      (
        const descriptor::Dataset &RESCALED_DATA,
        const ObjectiveFunctionInterface &OBJECTIVE,
        const size_t &NUMBER_THREADS = size_t( 1)
      );

      //! @brief select whether or not to backpropagate the current feature, and edit the error vector, if desired
      //! @param PREDICTION the prediction that was made by the neural network
      //! @param ERROR Reference to the error vector, (should already have been set to RESULT - PREDICTION)
      //! @param FEATURE_ID id of this feature in the dataset
      //! @param THREAD_ID id of this thread
      //! @return true if the feature should be backpropagated
      bool ShouldBackpropagate
      (
        const linal::VectorConstInterface< float> &PREDICTION,
        linal::VectorInterface< float> &ERROR,
        const size_t &FEATURE_ID,
        const size_t &THREAD_ID
      );

      //! @brief finalize the current round; occurs only after all threads were already joined
      void FinalizeRound();

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class NeuralNetworkSelectiveBackpropagationHybrid

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_NEURAL_NETWORK_SELECTIVE_BACKPROPAGATION_HYBRID_H_
