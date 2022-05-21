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

// include forward header of this class
#include "model/bcl_model_score_derivative_ensemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_running_average_sd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    //! @brief default constructor
    ScoreDerivativeEnsemble::ScoreDerivativeEnsemble() :
      m_ConsistencyWeight( 0.0),
      m_ConsistencyBestWeight( 0.0),
      m_WeightSquareWeight( 0.0),
      m_AverageWeightAbsWeight( 0.0),
      m_RawAverageWeight( 0.0),
      m_UtilityWeight( 0.0),
      m_Balance( true),
      m_UseCategorical( true)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ScoreDerivativeEnsemble
    ScoreDerivativeEnsemble *ScoreDerivativeEnsemble::Clone() const
    {
      return new ScoreDerivativeEnsemble( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ScoreDerivativeEnsemble::GetClassIdentifier() const
    {
      return GetStaticClassName< ScoreDerivativeEnsemble>();
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ScoreDerivativeEnsemble::GetAlias() const
    {
      static const std::string s_Name( "ScoreDerivativeEnsemble");
      return s_Name;
    }

    //! @brief test whether or not this score requires model scoring
    //! @return true if this scores uses model scoring
    bool ScoreDerivativeEnsemble::GetUsesModelScores() const
    {
      return m_ConsistencyBestWeight || m_UtilityWeight;
    }

    //! @brief initialize balancing (weighting of columns based on prior distribution)
    //! @param MATRIX results matrix
    //! @param NR_DESCRIPTORS number of descriptor columns
    void ScoreDerivativeEnsemble::InitializeBalancing
    (
      const storage::Vector< linal::Matrix< char> > &MODEL_CLASSIFICATIONS,
      const size_t &NR_DESCRIPTORS
    ) const
    {
      const size_t result_size( MODEL_CLASSIFICATIONS( 0).GetNumberCols());
      const size_t n_features( MODEL_CLASSIFICATIONS( 0).GetNumberRows());
      m_WeightsTruePositives = m_WeightsFalsePositives = m_WeightsTrueNegatives = m_WeightsFalseNegatives =
        linal::Vector< float>( result_size, float( 1.0) / float( result_size));
      m_AveWeightSomeTP = m_AveWeightSomeTN = m_AveWeightSomeFP = m_AveWeightSomeFN
        = storage::Vector< storage::Vector< math::RunningAverage< float> > >
          (
            result_size,
            storage::Vector< math::RunningAverage< float> >( NR_DESCRIPTORS)
          );
      m_ConsistencySomeTP = m_ConsistencySomeTN = m_ConsistencySomeFP = m_ConsistencySomeFN = m_ConsistencyP = m_ConsistencyN
        = storage::Vector< storage::Vector< math::RunningAverage< float> > >
          (
            result_size,
            storage::Vector< math::RunningAverage< float> >( NR_DESCRIPTORS)
          );
      const size_t n_models( MODEL_CLASSIFICATIONS.GetSize());
      storage::Vector< linal::Vector< size_t> > count_models_correct_p
      (
        result_size,
        linal::Vector< size_t>( size_t( n_models + 1), size_t( 0))
      );
      storage::Vector< linal::Vector< size_t> > count_models_correct_n
      (
        result_size,
        linal::Vector< size_t>( size_t( n_models + 1), size_t( 0))
      );

      // iterate over each model's predictions
      // translate the predictions from m_Predictions into the predictions for this model
      linal::Vector< size_t> positive_counts( result_size, size_t( 0));
      linal::Vector< size_t> true_positive_counts( result_size, size_t( 0));
      linal::Vector< size_t> false_positive_counts( result_size, size_t( 0));
      linal::Vector< size_t> true_negative_counts( result_size, size_t( 0));
      linal::Vector< size_t> false_negative_counts( result_size, size_t( 0));

      for( size_t datum( 0); datum < n_features; ++datum)
      {
        for( size_t result( 0); result < result_size; ++result)
        {
          size_t count_correct( 0);
          const bool is_positive
          (
            MODEL_CLASSIFICATIONS( 0)( datum, result) == 'P'
            || MODEL_CLASSIFICATIONS( 0)( datum, result) == 'n'
          );
          const char target_char( is_positive ? 'P' : 'N');
          for( size_t model( 0); model < n_models; ++model)
          {
            if( MODEL_CLASSIFICATIONS( model)( datum, result) == target_char)
            {
              ++count_correct;
            }
          }
          if( is_positive)
          {
            ++count_models_correct_p( result)( count_correct);
            positive_counts( result) += n_models;
            true_positive_counts( result) += count_correct;
            false_negative_counts( result) += n_models - count_correct;
          }
          else
          {
            ++count_models_correct_n( result)( count_correct);
            true_negative_counts( result) += count_correct;
            false_positive_counts( result) += n_models - count_correct;
          }
        }
      }

      BCL_MessageStd( "# models correct histogram for positives: " + util::Format()( count_models_correct_p));
      BCL_MessageStd( "# models correct histogram for negatives: " + util::Format()( count_models_correct_n));

      m_PChanceConsistency = linal::Vector< float>( n_models + 1);
      for( size_t model_tally( 0); model_tally <= n_models; ++model_tally)
      {
        math::RunningAverage< float> ave;
        double divisor( math::Pow( 2.0, model_tally - 1.0));
        for( size_t i( 0), end( model_tally / 2); i <= end; ++i)
        {
          double prob( 0);
          // determine the probability of randomly choosing i particular models or less
          for( size_t j( 0); j <= i; ++j)
          {
            prob += math::BinomialCoefficient( model_tally, j);
          }
          prob /= divisor;
          ave += std::max( 0.0, 1.0 - prob);
        }
        m_PChanceConsistency( model_tally) = ave.GetAverage();
      }

      // if not balancing or weighting non-equally, just continue
      if( !m_Balance && !m_UseCategorical)
      {
        return;
      }

      // for balancing, balance the influence of P vs N, and separately T vs F, so
      // P' = N'
      // T' = F'
      // where P' = x1 * TP + x2 * FN, N' = x3 * TN + x4 * FP, T' = x1 * TP + x3 * TN, F' = x2 * FN + x4 * FP
      // to properly constraint the equations, also require that
      // P = TP + FN = P'
      // x1 * TP = x3 * TN
      // This yields the solutions:
      // x1 = P / ( 2 * TP)
      // x2 = P / ( 2 * FN)
      // x3 = P / ( 2 * TN)
      // x4 = P / ( 2 * FP)
      if( m_Balance)
      {
        for( size_t j( 0); j < result_size; ++j)
        {
          // get half the positives
          const float half_positives( 0.5 * positive_counts( j));
          if( true_positive_counts( j))
          {
            m_WeightsTruePositives( j) = half_positives / float( true_positive_counts( j));
          }
          if( false_negative_counts( j))
          {
            m_WeightsFalseNegatives( j) = half_positives / float( false_negative_counts( j));
          }
          if( true_negative_counts( j))
          {
            m_WeightsTrueNegatives( j) = half_positives / float( true_negative_counts( j));
          }
          if( false_positive_counts( j))
          {
            m_WeightsFalsePositives( j) = half_positives / float( false_positive_counts( j));
          }
        }
      }
      if( m_UseCategorical && result_size > size_t( 1))
      {
        size_t sum_counts( positive_counts.Sum());
        size_t excess( sum_counts % n_features);

        // if there is a nearly integral # of counts per feature and there are least 1 point per feature
        if( sum_counts >= n_features && !excess)
        {
          for( size_t j( 0); j < result_size; ++j)
          {
            // get half the positives
            const float positive_frequency( float( positive_counts( j)) / float( n_features));
            m_WeightsTruePositives( j) *= positive_frequency;
            m_WeightsFalseNegatives( j) *= positive_frequency;
            m_WeightsTrueNegatives( j) *= positive_frequency;
            m_WeightsFalsePositives( j) *= positive_frequency;
          }
        }
        else
        {
          BCL_MessageCrt( "Requested categorical balancing but it appears invalid for this dataset!");
        }
      }
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief score a given vector of matrices
    //! @param MODEL_DESCRIPTOR_DERIVATIVES(x)(y,z) corresponds to the derivative of feature y on output z for model x
    //! @param PREDICTION_CLASS binary (0,1) score for each model; used for calculating f-score measure only
    //! @return an agglomerated score
    linal::Vector< float> ScoreDerivativeEnsemble::Score
    (
      const storage::Vector< linal::Matrix< float> > &MODEL_DESCRIPTOR_DERIVATIVES,
      const storage::Vector< std::string> &PREDICTION_CLASS
    ) const
    {
      math::RunningAverage< linal::Vector< float> > averages;
      this->Score( MODEL_DESCRIPTOR_DERIVATIVES, PREDICTION_CLASS, averages);
      return averages.GetAverage();
    }

    //! @brief score a given vector of matrices and add it to a running average with the appropriate weight
    //! @param MODEL_DESCRIPTOR_DERIVATIVES(x)(y,z) corresponds to the derivative of feature y on output z for model x
    //! @param PREDICTION_CLASS PpNn\0 for each model, indicating whether the result was a TP,FP,TN,FN, or NA for each model
    //! @param AVERAGES running average of vector with one value for each descriptor column
    //! This overloaded version allows for balancing
    void ScoreDerivativeEnsemble::Score
    (
      const storage::Vector< linal::Matrix< float> > &MODEL_DESCRIPTOR_DERIVATIVES,
      const storage::Vector< std::string> &PREDICTION_CLASS,
      math::RunningAverage< linal::Vector< float> > &AVERAGES
    ) const
    {
      const size_t n_models( MODEL_DESCRIPTOR_DERIVATIVES.GetSize());

      // handle trivial case where no models or no derivatives were given
      if( n_models == size_t( 0) || MODEL_DESCRIPTOR_DERIVATIVES.GetSize() == size_t( 0))
      {
        return;
      }

      const size_t feature_size( MODEL_DESCRIPTOR_DERIVATIVES( 0).GetNumberRows());
      const size_t result_size( MODEL_DESCRIPTOR_DERIVATIVES( 0).GetNumberCols());
      BCL_Assert
      (
        m_ResultColumnWeighting.IsEmpty() || m_ResultColumnWeighting.GetSize() == result_size,
        "Need same number of result column weights as outputs or none at all!"
      );

      // create multiplier for each row, based on WAS_ABOVE_CUTOFF
      linal::Vector< float> multiplier( result_size, float( 1.0));
      for( size_t i( 0); i < result_size; ++i)
      {
        size_t n_correct( 0);
        const char type( PREDICTION_CLASS( 0)[ i]);
        const bool is_positive( type == 'P' || type == 'n' || type == '\0');
        for( size_t model( 0); model < n_models; ++model)
        {
          if( !islower( PREDICTION_CLASS( model)[ i]))
          {
             ++n_correct;
          }
        }
        if( n_correct == 0)
        {
          // no models got this result correct; set weight to 0
          multiplier( i) = 0;
        }
        else if( is_positive)
        {
          multiplier( i) = m_WeightsTruePositives( i);
        }
        else
        {
          multiplier( i) = m_WeightsTrueNegatives( i);
        }
        if( !m_ResultColumnWeighting.IsEmpty())
        {
          multiplier( i) *= m_ResultColumnWeighting( i);
        }
      }

      float multiplier_sum( multiplier.Sum());
      if( !multiplier_sum)
      {
        // no weight for this class; likely all results were wrong
        return;
      }
      multiplier /= multiplier_sum;

      math::RunningAverage< linal::Matrix< float> > weights_stats;

      linal::Matrix< size_t> counts_above_zero( feature_size, result_size, size_t( 0));
      linal::Matrix< float> abs_derivative;
      math::RunningAverage< linal::Matrix< float> > abs_weights_stats;

      // for each model
      // 1. Obtain its matrix of weights
      // 2. Compute the average absolute value for each feature weight
      for( size_t model_number( 0); model_number < n_models; ++model_number)
      {
        // make a reference to the derivative vector
        const linal::Matrix< float> &derivative_matrix( MODEL_DESCRIPTOR_DERIVATIVES( model_number));
        weights_stats += derivative_matrix;
        abs_derivative = derivative_matrix;

        // update the count of values above zero
        for( size_t feature( 0); feature < feature_size; ++feature)
        {
          for( size_t result( 0); result < result_size; ++result)
          {
            if( derivative_matrix( feature, result) > 0.0)
            {
              ++counts_above_zero( feature, result);
            }
            else
            {
              abs_derivative( feature, result) = -abs_derivative( feature, result);
            }
          }
        }
        abs_weights_stats += abs_derivative;
      }

      linal::Matrix< float> abs_weights_ave( weights_stats.GetAverage());
      math::Absolute( abs_weights_ave);

      linal::Matrix< float> sign_purity( feature_size, result_size);

      const float n_modelsf( n_models);

      // Sum the matrices
      linal::Matrix< float> combined_measures( feature_size, result_size, float( 0.0));

      if( m_RawAverageWeight)
      {
        combined_measures += weights_stats.GetAverage();
        combined_measures *= m_RawAverageWeight;
      }

      if( m_ConsistencyWeight)
      {
        for( size_t i( 0); i < feature_size; ++i)
        {
          for( size_t j( 0); j < result_size; ++j)
          {
            // calculate the consistency of signs for this weight, scaled between 0 and 1
            // at 0, half the models have the opposite effect as the other half
            sign_purity( i, j) = std::max( counts_above_zero( i, j), n_models - counts_above_zero( i, j));
          }
        }
        // normalize sign purity, scaled between 0 and 1
        // at 0, half the models have the opposite effect as the other half
        sign_purity -= float( 0.75) * n_modelsf;
        sign_purity /= n_modelsf;
        sign_purity *= m_ConsistencyWeight * float( 4.0);
        combined_measures += sign_purity;
      }

      // rescale the weights matrix between between 0-1
      if( m_WeightSquareWeight)
      {
        math::RunningAverage< linal::Matrix< float> > weights_sq_stats;

        linal::Matrix< float> square_weight( feature_size, result_size);

        // for each model
        // 1. Obtain its matrix of weights
        // 2. Compute the average absolute value for each feature weight
        for( size_t model_number( 0); model_number < n_models; ++model_number)
        {
          // make a reference to the derivative vector
          const linal::Matrix< float> &derivative_matrix( MODEL_DESCRIPTOR_DERIVATIVES( model_number));

          // make a copy of the matrix, so that we can square it
          square_weight = derivative_matrix;

          // update the count of values above zero
          for( size_t feature( 0); feature < feature_size; ++feature)
          {
            for( size_t result( 0); result < result_size; ++result)
            {
              square_weight( feature, result) *= square_weight( feature, result);
            }
          }
          weights_sq_stats += square_weight;
        }

        square_weight = weights_sq_stats.GetAverage();
        for( size_t i( 0); i < feature_size; ++i)
        {
          for( size_t j( 0); j < result_size; ++j)
          {
            // calculate the consistency of signs for this weight, scaled between 0 and 1
            // at 0, half the models have the opposite effect as the other half
            square_weight( i, j) = math::Sqrt( square_weight( i, j));
          }
        }
        RescaleFeatureDataSet rescale_weight( square_weight, math::Range< float>( -1.0, 1.0), RescaleFeatureDataSet::e_AveStd);
        rescale_weight.RescaleMatrix( square_weight);
        square_weight *= m_WeightSquareWeight;
        combined_measures += square_weight;
      }

      // rescale the ave weights matrix between 0-1
      if( m_AverageWeightAbsWeight)
      {
        linal::Matrix< float> tmp_weight_ave( abs_weights_ave);
        RescaleFeatureDataSet rescale_weight( tmp_weight_ave, math::Range< float>( -1.0, 1.0), RescaleFeatureDataSet::e_AveStd);
        rescale_weight.RescaleMatrix( tmp_weight_ave); // comment out if you want raw rescaled changes
        tmp_weight_ave *= m_AverageWeightAbsWeight;
        combined_measures += tmp_weight_ave;
      }

      // model-performance based scores
      if( GetUsesModelScores())
      {
        counts_above_zero.SetZero();
        // take stats of the better and worse matrices
        storage::Vector< std::string>::const_iterator itr_good( PREDICTION_CLASS.Begin());
        linal::Matrix< size_t> counts_above_zero_worse( feature_size, result_size);
        linal::Vector< size_t> count_good_models( result_size, size_t( 0));
        linal::Vector< size_t> count_bad_models( result_size, n_models);
        for( size_t model_number( 0); model_number < n_models; ++model_number, ++itr_good)
        {
          std::string::const_iterator itr_result_good( itr_good->begin());
          const linal::Matrix< float> &derivative_matrix( MODEL_DESCRIPTOR_DERIVATIVES( model_number));
          linal::Matrix< float> derivatives_transposed( derivative_matrix.Transposed());
          for( size_t result( 0); result < result_size; ++result, ++itr_result_good)
          {
            if( !islower( *itr_result_good))
            {
              ++count_good_models( result);
              --count_bad_models( result);
              for( size_t feature( 0); feature < feature_size; ++feature)
              {
                if( derivative_matrix( feature, result) > 0.0)
                {
                  ++counts_above_zero( feature, result);
                }
              }
            }
            else if( m_UtilityWeight)
            {
              for( size_t feature( 0); feature < feature_size; ++feature)
              {
                if( derivative_matrix( feature, result) > 0.0)
                {
                  ++counts_above_zero_worse( feature, result);
                }
              }
            }
          }
        }

        sign_purity.SetZero();
        for( size_t i( 0); i < feature_size; ++i)
        {
          for( size_t j( 0); j < result_size; ++j)
          {
            const size_t n_better_models( count_good_models( j));
            if( n_better_models >= 2)
            {
              const float offset( 0.5 * n_better_models);
              // calculate the consistency of signs for this weight, scaled between 0 and 1
              // at 0, half the models have the opposite effect as the other half
              sign_purity( i, j) = ( std::max( counts_above_zero( i, j), n_better_models - counts_above_zero( i, j)) - offset) / offset;
            }
          }
        }
        linal::Matrix< float> worse_sign_purity( feature_size, result_size, float( 0.0));
        for( size_t i( 0); i < feature_size; ++i)
        {
          for( size_t j( 0); j < result_size; ++j)
          {
            const size_t n_worse_models( count_bad_models( j));
            if( n_worse_models >= 2)
            {
              const float offset( 0.5 * n_worse_models);
              // calculate the consistency of signs for this weight, scaled between 0 and 1
              // at 0, half the models have the opposite effect as the other half
              worse_sign_purity( i, j) = ( std::max( counts_above_zero_worse( i, j), n_worse_models - counts_above_zero_worse( i, j)) - offset) / offset;
            }
          }
        }

        // take stats of the better and worse matrices
        storage::Vector< math::RunningAverage< linal::Vector< float> > > better_weights_abs_stats( result_size);
        storage::Vector< math::RunningAverage< linal::Vector< float> > > worse_weights_abs_stats( result_size);
        itr_good = PREDICTION_CLASS.Begin();
        for( size_t model_number( 0); model_number < n_models; ++model_number, ++itr_good)
        {
          std::string::const_iterator itr_result_good( itr_good->begin());
          const linal::Matrix< float> &derivative_matrix( MODEL_DESCRIPTOR_DERIVATIVES( model_number));
          linal::Matrix< float> derivatives_transposed( derivative_matrix.Transposed());
          math::Absolute( derivatives_transposed);
          for( size_t result( 0); result < result_size; ++result, ++itr_result_good)
          {
            // check for P or N, which indicate TP & TN, respectively
            if( !islower( *itr_result_good))
            {
              better_weights_abs_stats( result) += derivatives_transposed.GetRow( result);
            }
            else
            {
              worse_weights_abs_stats( result) += derivatives_transposed.GetRow( result);
            }
          }
        }

        if( m_ConsistencyBestWeight)
        {
          for( size_t j( 0); j < result_size; ++j)
          {
            const size_t count_good( count_good_models( j));
            if( count_good < 2)
            {
              continue;
            }
            const bool is_positive( PREDICTION_CLASS( 0)[ j] == 'P' || PREDICTION_CLASS( 0)[ j] == 'n');
            const float multiplier_this_result( multiplier( j) * m_PChanceConsistency( count_good));
            storage::Vector< math::RunningAverage< float> > &consistency_array( is_positive ? m_ConsistencyP( j) : m_ConsistencyN( j));
            const linal::Vector< float> &better_model_weight_aves( better_weights_abs_stats( j).GetAverage());
            m_Mutex.Lock();
            for( size_t feature( 0); feature < feature_size; ++feature)
            {
              consistency_array( feature).AddWeightedObservation
              (
                sign_purity( feature, j),
                multiplier_this_result * better_model_weight_aves( feature)
              );
            }
            m_Mutex.Unlock();
          }
        }
        if( m_UtilityWeight)
        {
          for( size_t j( 0); j < result_size; ++j)
          {
            const size_t count_good( count_good_models( j));
            const size_t count_bad( count_bad_models( j));
            if( count_good && count_bad)
            {
              const linal::Vector< float> &better_model_weight_aves( better_weights_abs_stats( j).GetAverage());
              const linal::Vector< float> &worse_model_weight_aves( worse_weights_abs_stats( j).GetAverage());
              const float good_feature_multiplier( m_PChanceConsistency( count_good));
              const float bad_feature_multiplier( m_PChanceConsistency( count_bad));

              const bool is_positive( PREDICTION_CLASS( 0)[ j] == 'P' || PREDICTION_CLASS( 0)[ j] == 'n');

              storage::Vector< math::RunningAverage< float> > &weight_array_true( is_positive ? m_AveWeightSomeTP( j) : m_AveWeightSomeTN( j));
              storage::Vector< math::RunningAverage< float> > &weight_array_false( is_positive ? m_AveWeightSomeFN( j) : m_AveWeightSomeFP( j));
              storage::Vector< math::RunningAverage< float> > &consistency_array_true( is_positive ? m_ConsistencySomeTP( j) : m_ConsistencySomeTN( j));
              storage::Vector< math::RunningAverage< float> > &consistency_array_false( is_positive ? m_ConsistencySomeFN( j) : m_ConsistencySomeFP( j));

              m_Mutex.Lock();
              for( size_t feature( 0); feature < feature_size; ++feature)
              {
                consistency_array_true( feature).AddWeightedObservation
                (
                  sign_purity( feature, j),
                  good_feature_multiplier * better_model_weight_aves( feature)
                );
                weight_array_true( feature).AddWeightedObservation( better_model_weight_aves( feature), count_good);
                consistency_array_false( feature).AddWeightedObservation
                (
                  worse_sign_purity( feature, j),
                  bad_feature_multiplier * worse_model_weight_aves( feature)
                );
                weight_array_false( feature).AddWeightedObservation( worse_model_weight_aves( feature), count_bad);
              }
              m_Mutex.Unlock();
            }
          }
        }
      }
      if( m_NonUtilityWeight)
      {
        linal::Vector< float> norm( feature_size, float( 0));
        for( size_t i( 0); i < feature_size; ++i)
        {
          norm( i) = linal::VectorConstReference< float>( result_size, combined_measures[ i]) * multiplier;
        }

        AVERAGES.AddWeightedObservation( norm, multiplier_sum * m_NonUtilityWeight);
      }
    }

    //! @brief add the utility score to the overall score
    void ScoreDerivativeEnsemble::AddUtilityScore
    (
      math::RunningAverage< linal::Vector< float> > &AVERAGES
    ) const
    {
      if( !m_UtilityWeight && !m_ConsistencyBestWeight)
      {
        return;
      }

      const size_t n_features( m_ConsistencySomeTP( 0).GetSize()), n_results( m_ConsistencySomeTP.GetSize());
      linal::Vector< float> utility( n_features, float( 0.0));
      linal::Vector< float> consistency_best( n_features, float( 0.0));

      const double original_weight( AVERAGES.GetWeight());

      linal::Vector< float> result_weights( n_results);
      for( size_t j( 0); j < n_results; ++j)
      {
        result_weights( j) = m_WeightsTruePositives( j) + m_WeightsFalsePositives( j) + m_WeightsTrueNegatives( j) + m_WeightsTruePositives( j);
      }
      result_weights *= float( n_results) / result_weights.Sum();
      if( m_UtilityWeight)
      {
        //math::RunningAverage< float> average;
        math::RunningMinMax< float> minmax;
        math::RunningAverage< float> ave_ctp, ave_cfp, ave_ctn, ave_cfn;
        for( size_t i( 0); i < n_features; ++i)
        {
          for( size_t j( 0); j < n_results; ++j)
          {
            ave_ctp += m_ConsistencySomeTP( j)( i);
            ave_cfp += m_ConsistencySomeFP( j)( i);
            ave_ctn += m_ConsistencySomeTN( j)( i);
            ave_cfn += m_ConsistencySomeFN( j)( i);
          }
        }
        if( ave_ctp.GetAverage() < 1.0e-8)
        {
          ave_ctp.Reset();
          ave_ctp += 1.0e-8;
        }
        if( ave_cfp.GetAverage() < 1.0e-8)
        {
          ave_cfp.Reset();
          ave_cfp += 1.0e-8;
        }
        if( ave_ctn.GetAverage() < 1.0e-8)
        {
          ave_ctn.Reset();
          ave_ctn += 1.0e-8;
        }
        if( ave_cfn.GetAverage() < 1.0e-8)
        {
          ave_cfn.Reset();
          ave_cfn += 1.0e-8;
        }
        const float alpha( ave_ctp.GetAverage() / ave_cfn.GetAverage());
        const float beta( ave_ctp.GetAverage() / ave_cfp.GetAverage());
        const float chi( ave_ctp.GetAverage() / ave_ctn.GetAverage());
        for( size_t i( 0); i < n_features; ++i)
        {
          for( size_t j( 0); j < n_results; ++j)
          {
            float consistency_add( 0);
            const float consistency_tp( m_ConsistencySomeTP( j)( i).GetAverage());
            const float consistency_fp( beta * m_ConsistencySomeFP( j)( i).GetAverage());
            const float consistency_tn( chi * m_ConsistencySomeTN( j)( i).GetAverage());
            const float consistency_fn( alpha * m_ConsistencySomeFN( j)( i).GetAverage());
            if( consistency_tp > consistency_fn)
            {
              consistency_add += ( consistency_tp - consistency_fn) * m_AveWeightSomeTP( j)( i).GetAverage();
            }
            else
            {
              consistency_add += ( consistency_tp - consistency_fn) * m_AveWeightSomeFN( j)( i).GetAverage();
            }
            if( consistency_tn > consistency_fp)
            {
              consistency_add += ( consistency_tn - consistency_fp) * m_AveWeightSomeTN( j)( i).GetAverage();
            }
            else
            {
              consistency_add += ( consistency_tn - consistency_fp) * m_AveWeightSomeFP( j)( i).GetAverage();
            }
            consistency_add *= result_weights( j);
            utility( i) += consistency_add;
          }
          minmax += utility( i);
        }
        utility *= float( n_features) / float( n_results);
        if( m_NonUtilityWeight)
        {
          AVERAGES.AddWeightedObservation( utility, original_weight * m_UtilityWeight / m_NonUtilityWeight);
        }
        else
        {
          AVERAGES.AddWeightedObservation( utility, m_UtilityWeight);
        }
      }
      if( m_ConsistencyBestWeight)
      {
        for( size_t i( 0); i < n_features; ++i)
        {
          for( size_t j( 0); j < n_results; ++j)
          {
            consistency_best( i) += result_weights( j) * ( m_ConsistencyN( j)( i).GetAverage() + m_ConsistencyP( j)( i).GetAverage()) / 2.0;
          }
        }
        if( m_NonUtilityWeight)
        {
          AVERAGES.AddWeightedObservation( consistency_best, original_weight * m_ConsistencyBestWeight / m_NonUtilityWeight);
        }
        else
        {
          AVERAGES.AddWeightedObservation( consistency_best, m_ConsistencyBestWeight);
        }
      }
      for( size_t i( 0); i < n_features; ++i)
      {
        for( size_t j( 0); j < n_results; ++j)
        {
          BCL_MessageVrb
          (
            "Feature: " + util::Format()( i) + " Result: " + util::Format()( j)
            + " utility: " + util::Format()( utility( i))
            + " consistency-best: " + util::Format()( consistency_best( i))
          );
        }
      }
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write errors out to
    bool ScoreDerivativeEnsemble::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      const double sum
      (
        m_ConsistencyWeight + m_ConsistencyBestWeight
        + m_WeightSquareWeight + m_AverageWeightAbsWeight
        + m_RawAverageWeight + m_UtilityWeight
      );
      m_ConsistencyWeight /= sum;
      m_ConsistencyBestWeight /= sum;
      m_WeightSquareWeight /= sum;
      m_AverageWeightAbsWeight /= sum;
      m_RawAverageWeight /= sum;
      m_UtilityWeight /= sum;
      m_NonUtilityWeight = 1.0 - m_UtilityWeight - m_ConsistencyBestWeight;
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ScoreDerivativeEnsemble::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates a weighted sum of desired measures based on derivatives or pseudo-derivatives of a model"
      );

      parameters.AddInitializer
      (
        "consistency",
        "Change the weight of sign consistency on the final result. Higher values place more emphasis on columns that have the same effect on the "
        "output across the models",
        io::Serialization::GetAgentWithRange( &m_ConsistencyWeight, 0.0, 1.0),
        "0.0"
      );
      parameters.AddInitializer
      (
        "consistency best",
        "Change the weight of sign consistency for the best models on the final result",
        io::Serialization::GetAgentWithRange( &m_ConsistencyBestWeight, 0.0, 1.0),
        "0.0"
      );
      parameters.AddInitializer
      (
        "square",
        "Change the weight of the average (weight influence ^ 2) (and rescaled) on the score output",
        io::Serialization::GetAgentWithRange( &m_WeightSquareWeight, 0.0, 1.0),
        "0.0"
      );
      parameters.AddInitializer
      (
        "absolute",
        "Change the weight of the average (weight) (and rescaled) on the score output",
        io::Serialization::GetAgentWithRange( &m_AverageWeightAbsWeight, 0.0, 1.0),
        "0.0"
      );
      parameters.AddInitializer
      (
        "utility",
        "Change the weight of the utility function of input sensitivity.  The utility function considers whether models "
        "that gave incorrect output for a particular feature consistently used the feature more than the models that performed "
        "well on that particular feature",
        io::Serialization::GetAgentWithRange( &m_UtilityWeight, 0.0, 1.0),
        "0.0"
      );
      parameters.AddInitializer
      (
        "average",
        " This score should never be used for feature selection, but is sometimes useful in analysis of the directional "
        "influence of descriptors.  Likewise, it should always be used by itself alone",
        io::Serialization::GetAgentWithRange( &m_RawAverageWeight, 0.0, 1.0),
        "0.0"
      );
      parameters.AddInitializer
      (
        "balance",
        "True to weight columns by 1-column frequency",
        io::Serialization::GetAgent( &m_Balance),
        "True"
      );
      parameters.AddInitializer
      (
        "categorical",
        "True to have a constant weight per feature, divided equally amongst columns above vs below the cutoff",
        io::Serialization::GetAgent( &m_UseCategorical),
        "False"
      );
      parameters.AddOptionalInitializer
      (
        "result weights",
        "vector of weights; 1 value per result column in the dataset, or empty to weight all result columns equally",
        io::Serialization::GetAgent( &m_ResultColumnWeighting)
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl

