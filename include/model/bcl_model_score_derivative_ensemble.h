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

#ifndef BCL_MODEL_SCORE_DERIVATIVE_ENSEMBLE_H
#define BCL_MODEL_SCORE_DERIVATIVE_ENSEMBLE_H

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_retrieve_interface.h"
#include "bcl_model_score_dataset_interface.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_running_average.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_running_min_max.h"
#include "sched/bcl_sched_mutex.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ScoreDerivativeEnsemble
    //! @brief calculates the saliency of all input features by taking the sum of absolute weights emanating from each input neuron
    //!
    //! @author mendenjl
    //! @see @link example_model_score_derivative_ensemble.cpp @endlink
    //! @date May 22, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ScoreDerivativeEnsemble :
      public util::SerializableInterface
    {
    private:

    //////////
    // data //
    //////////

      //! weight of sign consistency
      float m_ConsistencyWeight;

      //! Weight of the consistency of the best models
      float m_ConsistencyBestWeight;

      //! weight for square of weights
      float m_WeightSquareWeight;

      //! weight for average weight absolute
      float m_AverageWeightAbsWeight;

      //! weight for raw average score
      //! This score should never be used for feature selection, but is sometimes useful in analysis of the directional
      //! influence of various descriptors
      float m_RawAverageWeight;

      //! weight of utility score
      float m_UtilityWeight;

      //! weight of all scores added during each Score call, as opposed to during AddUtilityScore
      float m_NonUtilityWeight;

      //! flag to balance (weight columns by frequency of class)
      bool m_Balance;

      //! flag to use categorical-based balancing (appropriate when using CategoricalMax as the objective function)
      bool m_UseCategorical;

      //! Weighting for result columns; must be same size as # of results, or empty for default weights
      linal::Vector< float> m_ResultColumnWeighting;

      //! column weights for TP
      mutable linal::Vector< float> m_WeightsTruePositives;

      //! column weights for FP
      mutable linal::Vector< float> m_WeightsFalsePositives;

      //! column weights for TN
      mutable linal::Vector< float> m_WeightsTrueNegatives;

      //! column weights for FN
      mutable linal::Vector< float> m_WeightsFalseNegatives;

      //! Mutex for changing m_AveInfluence or any of the other variables in this class
      mutable sched::Mutex m_Mutex;

      //! for each result track utility, weighted by class weight
      //! Utility is defined as Abs Sum Correct, weighted by class weight
      mutable storage::Vector< storage::Vector< math::RunningAverage< float> > > m_ConsistencyP;
      mutable storage::Vector< storage::Vector< math::RunningAverage< float> > > m_ConsistencyN;
      mutable storage::Vector< storage::Vector< math::RunningAverage< float> > > m_ConsistencySomeTP;
      mutable storage::Vector< storage::Vector< math::RunningAverage< float> > > m_ConsistencySomeFN;
      mutable storage::Vector< storage::Vector< math::RunningAverage< float> > > m_ConsistencySomeTN;
      mutable storage::Vector< storage::Vector< math::RunningAverage< float> > > m_ConsistencySomeFP;
      mutable storage::Vector< storage::Vector< math::RunningAverage< float> > > m_AveWeightSomeTP;
      mutable storage::Vector< storage::Vector< math::RunningAverage< float> > > m_AveWeightSomeFN;
      mutable storage::Vector< storage::Vector< math::RunningAverage< float> > > m_AveWeightSomeTN;
      mutable storage::Vector< storage::Vector< math::RunningAverage< float> > > m_AveWeightSomeFP;

      mutable linal::Vector< float> m_PChanceConsistency;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ScoreDerivativeEnsemble();

      //! @brief Clone function
      //! @return pointer to new ScoreDerivativeEnsemble
      ScoreDerivativeEnsemble *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief test whether or not this score requires model scoring
      //! @return true if this scores uses model scoring
      bool GetUsesModelScores() const;

      //! @brief test whether or not this score balances based on a cutoff
      //! @return true if this score balances based on a cutoff
      bool GetDoesBalance() const
      {
        return m_Balance;
      }

      //! @brief initialize balancing (weighting of columns based on distribution)
      //! @param MODEL_CLASSIFICATIONS classifications of the model
      //! @param NR_DESCRIPTORS number of descriptor columns
      void InitializeBalancing
      (
        const storage::Vector< linal::Matrix< char> > &MODEL_CLASSIFICATIONS,
        const size_t &NR_DESCRIPTORS
      ) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief score a given vector of matrices; this version does not allow for balancing
      //! @param MODEL_DESCRIPTOR_DERIVATIVES(x)(y,z) corresponds to the derivative of feature y on output z for model x
      //! @param PREDICTION_CLASS PpNn\0 for each model, indicating whether the result was a TP,FP,TN,FN, or NA for each model
      //! @return an agglomerated score
      linal::Vector< float> Score
      (
        const storage::Vector< linal::Matrix< float> > &MODEL_DESCRIPTOR_DERIVATIVES,
        const storage::Vector< std::string> &PREDICTION_CLASS
      ) const;

      //! @brief score a given vector of matrices and add it to a running average with the appropriate weight
      //! @param MODEL_DESCRIPTOR_DERIVATIVES(x)(y,z) corresponds to the derivative of feature y on output z for model x
      //! @param PREDICTION_CLASS PpNn\0 for each model, indicating whether the result was a TP,FP,TN,FN, or NA for each model
      //! @param AVERAGES running average of vector with one value for each descriptor column
      //! This overloaded version allows for balancing
      void Score
      (
        const storage::Vector< linal::Matrix< float> > &MODEL_DESCRIPTOR_DERIVATIVES,
        const storage::Vector< std::string> &PREDICTION_CLASS,
        math::RunningAverage< linal::Vector< float> > &AVERAGES
      ) const;

      //! @brief add the utility score to the overall score
      void AddUtilityScore
      (
        math::RunningAverage< linal::Vector< float> > &AVERAGES
      ) const;

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

    }; // class ScoreDerivativeEnsemble

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_SCORE_DERIVATIVE_ENSEMBLE_H

