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

#ifndef BCL_MODEL_SCORE_DATASET_NEURAL_NETWORK_INPUT_SENSITIVITY_H
#define BCL_MODEL_SCORE_DATASET_NEURAL_NETWORK_INPUT_SENSITIVITY_H

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_retrieve_interface.h"
#include "bcl_model_score_dataset_interface.h"
#include "bcl_model_score_derivative_ensemble.h"
#include "linal/bcl_linal_vector.h"
#include "sched/bcl_sched_mutex.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ScoreDatasetNeuralNetworkInputSensitivity
    //! @brief calculates the sensitivity of an objective function to a model to given inputs
    //!
    //! @author mendenjl
    //! @see @link example_model_score_dataset_neural_network_input_sensitivity.cpp @endlink
    //! @date Jul 01, 2013
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ScoreDatasetNeuralNetworkInputSensitivity :
      public ScoreDatasetInterface
    {
    private:

    //////////
    // data //
    //////////

      //! Model storage; to retrieve models and descriptors
      util::Implementation< RetrieveInterface> m_ModelStorage;

      //! Optional key of the desired model
      std::string m_Key;

      //! pointer to vector of models used by this property
      util::SiPtr< const storage::Vector< util::OwnPtr< NeuralNetwork> > > m_Models;

      //! Final objective function
      mutable util::Implementation< ObjectiveFunctionInterface> m_Objective;

      //! scorer for the derivatives calculated by this function
      ScoreDerivativeEnsemble m_Scorer;

      //! Map from initializer string/type to the shptrvector of model interfaces; saves loading potentially
      //! gigantic models repetitively
      static storage::Map< std::string, storage::Vector< util::OwnPtr< NeuralNetwork> > > s_Models;

      //! Data used by thread
      mutable util::SiPtr< const descriptor::Dataset>    m_DatasetPtr;    //!< Dataset pointer
      mutable util::SiPtr< const linal::Vector< float> > m_DescriptorStd; //!< Standard deviation of the descriptor
      mutable storage::Vector< math::RunningAverage< linal::Vector< float> > > m_AveDescriptorScores; //!< Scores for each column
      mutable size_t                                     m_FeatureNumber;
      mutable sched::Mutex                               m_FeatureNumberMutex; //!< Mutex for changing feature-number

      mutable size_t                                     m_FeaturesToComputeSize; //!< Size of m_FeaturesToCompute, cached
      mutable storage::Vector< linal::Matrix< char> >    m_PredictionClassifications; //!< Actual predictions of the models (

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new ScoreDatasetNeuralNetworkInputSensitivity
      ScoreDatasetNeuralNetworkInputSensitivity *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief score a given dataset
      //! @param DATASET dataset of interest
      //! @return scores of the dataset
      linal::Vector< float> Score( const descriptor::Dataset &DATASET) const;

      //! @brief determine the side of the actual result
      //! @param ACTUAL actual / experimental value
      //! @param CUTOFF the cutoff value
      //! @return a vector containing size_ts 0/1/undefined depending on whether the model was bad, good, or unknown
      static storage::Vector< size_t> GetCutoffSides
      (
        const linal::VectorConstInterface< float> &ACTUAL,
        const float &CUTOFF
      );

      //! @brief determine the performance of the model ensemble on a particular feature
      //! @param PREDICTION_CLASSIFICATIONS vector of model prediction classifications
      //! @param FEATURE_NR feature number of interest
      //! @return a vector containing strings with P|p|N|n|\0 for whether each result is a TP, FP, TN, FN< or NA for each model
      static storage::Vector< std::string> PartitionModels
      (
        const storage::Vector< linal::Matrix< char> > &PREDICTION_CLASSIFICATIONS,
        const size_t &FEATURE_NR
      );

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write errors out to
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief Get the next feature for a thread to make predictions for
      //! @param LAST_FEATURE last feature that this thread computed
      //! @return the next feature for a thread to make predictions for
      size_t GetNextFeatureForPrediction( const size_t &LAST_FEATURE) const;

      //! @brief run a thread to compute the kernel for all features with ID = THREAD_ID % n_threads
      //! @param THREAD_ID id of the thread to run (0-indexed)
      void RunThread( const size_t &THREAD_ID) const;

    }; // class ScoreDatasetNeuralNetworkInputSensitivity

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_SCORE_DATASET_NEURAL_NETWORK_INPUT_SENSITIVITY_H

