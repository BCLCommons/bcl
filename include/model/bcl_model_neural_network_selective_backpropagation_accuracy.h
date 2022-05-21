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

#ifndef BCL_MODEL_NEURAL_NETWORK_SELECTIVE_BACKPROPAGATION_ACCURACY_H_
#define BCL_MODEL_NEURAL_NETWORK_SELECTIVE_BACKPROPAGATION_ACCURACY_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_neural_network_selective_backpropagation_interface.h"
#include "math/bcl_math_linear_function.h"
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class NeuralNetworkSelectiveBackpropagationAccuracy
    //! @brief Selective backpropagation for online balancing for classification problems; agnostic to side of cutoff
    //!
    //! @see @link example_model_neural_network_selective_backpropagation_accuracy.cpp @endlink
    //! @author mendenjl
    //! @date Jun 11, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API NeuralNetworkSelectiveBackpropagationAccuracy :
      public NeuralNetworkSelectiveBackpropagationInterface
    {
    private:

    //////////
    // data //
    //////////

      //! Parameters
      size_t m_ResultStartID;  //!< Min result to consider
      size_t m_ResultEndID;    //!< 1+Max result to consider
      float  m_IgnoreTruePredictionsBelow;   //!< Saturation boundary below cutoff for accurate predictions
      float  m_IgnoreTruePredictionsAbove;   //!< Saturation boundary above cutoff for accurate predictions
      float  m_IgnoreFalsePredictionsBelow;  //!< Saturation boundary below cutoff for inaccurate predictions
      float  m_IgnoreFalsePredictionsAbove;  //!< Saturation boundary above cutoff for inaccurate predictions
      float  m_SaturationUpperFalse;  //!< Saturation boundary above cutoff for inaccurate predictions
      float  m_SaturationLowerFalse;  //!< Saturation boundary below cutoff for inaccurate predictions
      bool   m_PureClassification; //!< True if training a pure classifier
      bool   m_BalanceError;       //!< True to try to balance error based on cutoff side
      bool   m_FlatErrorFunction;  //!< True to flatten the error function outside the saturation region

      size_t m_RoundsBeforeRefinementStage1; //!< # of rounds before stage 1 of refinement

      //! Constant data (after Initialize is called)
      size_t                m_ResultsSize;    //!< Results size, cached
      size_t                m_NumberThreads;  //!< Number of threads
      linal::Vector< float> m_ScaledCutoffs;  //!< Rescaled cutoffs, one for each result
      util::ShPtr< FeatureDataSet< float> > m_Results; //!< Results dataset

      size_t                m_RoundNumber;

      //! Inaccuracy last round; used to alter the probability of backpropagation if a prediction is made closer to the
      //! cutoff than the actual result was.  Updated after each iteration, using m_IncorrectTally
      linal::Vector< float>   m_LastRoundInaccuraciesAboveCutoff;
      linal::Vector< float>   m_LastRoundInaccuraciesBelowCutoff;
      linal::Vector< float>   m_LastRoundRelativeErrorAboveCutoff;
      linal::Vector< float>   m_LastRoundRelativeErrorBelowCutoff;

      //! Thread-local variables
      //! # of incorrect results computed by each thread last turn
      storage::Vector< linal::Vector< size_t> >    m_IncorrectAboveCutoffTally;
      storage::Vector< linal::Vector< size_t> >    m_IncorrectBelowCutoffTally;
      storage::Vector< linal::Vector< size_t> >    m_CorrectAboveCutoffTally;
      storage::Vector< linal::Vector< size_t> >    m_CorrectBelowCutoffTally;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      NeuralNetworkSelectiveBackpropagationAccuracy();

      //! @brief copy constructor
      //! @return a new NeuralNetworkSelectiveBackpropagationAccuracy copied from this instance
      NeuralNetworkSelectiveBackpropagationAccuracy *Clone() const;

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

    }; // class IterateResilientProagation

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_NEURAL_NETWORK_SELECTIVE_BACKPROPAGATION_ACCURACY_H_
