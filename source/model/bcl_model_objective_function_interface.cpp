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
#include "model/bcl_model_objective_function_interface.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_running_average_sd.h"
#include "model/bcl_model_feature_data_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    //! @brief Goal as string
    //! @param TYPE the goal
    //! @return the string for the goal
    const std::string &ObjectiveFunctionInterface::GetGoalName( const Goal &TYPE)
    {
      static const std::string s_names[] =
      {
        "Classification",     // Objective is based on ranking or a threshold
        "RankClassification", // Objective is based on ranking
        "Regression",         // Objective is based on minimizing error
        "Other",              // Objective is based on some other criteria
        GetStaticClassName< Goal>()
      };
      return s_names[ TYPE];
    }

    //! @brief get the threshold, for classification type objectives
    //! @return the threshold, for classification type objectives
    float ObjectiveFunctionInterface::GetThreshold() const
    {
      return util::GetUndefined< float>();
    }

    //! @brief get the parity, for rank classification type objectives
    //! @return the parity, for rank classification type objectives
    //! @note this has meaning only for rank-classification objectives; true means the objective is most interested
    //!       in prediction of values higher than the threshold, false means below threshold
    bool ObjectiveFunctionInterface::GetRankingParity() const
    {
      return true;
    }

    //! @brief get the desired hit rate (e.g. fraction of positive predictions)
    //! @return the desired hit rate, for rank classification type objectives
    //! @note this has meaning only for rank-classification objectives; 0.01 means, for example, that only the top 1%
    //!       of values will be considered
    float ObjectiveFunctionInterface::GetDesiredHitRate() const
    {
      return util::GetUndefined< float>();
    }

    //! @brief classify. Obtain a matrix of with N|n for all predicted negatives, P|p for all predicted positives, \0 for all non-predicted values
    //!        case indicates whether the prediction was true (upper) or false (lower)
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return matrix with PpNn\0 values indicating TP,FP,TN,FN,NA
    linal::Matrix< char> ObjectiveFunctionInterface::GetFeaturePredictionClassifications
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      const size_t dataset_size( EXPERIMENTAL.GetNumberFeatures());
      const size_t result_size( EXPERIMENTAL.GetFeatureSize());
      linal::Matrix< char> prediction_classification( dataset_size, result_size, '\0');
      const float threshold( GetThreshold());
      const float desired_hit_rate( GetDesiredHitRate());
      const bool ranking_parity( GetRankingParity());
      if( GetGoalType() == e_RankClassification)
      {
        storage::Vector< float> column( dataset_size, 0.0);
        for( size_t j( 0); j < result_size; ++j)
        {
          float prediction_cutoff( 0.0);
          if( util::IsDefined( desired_hit_rate))
          {
            for( size_t i( 0); i < dataset_size; ++i)
            {
              column( i) = PREDICTED( i)( j);
            }
            if( ranking_parity)
            {
              column.Sort( std::greater< float>());
            }
            else
            {
              column.Sort( std::less< float>());
            }

            // get the prediction cutoff
            prediction_cutoff =
              column
              (
                std::min
                (
                  size_t( column.GetSize() - 1),
                  size_t( column.GetSize() * desired_hit_rate)
                )
              );
          }
          else
          {
            column.Resize( 0);
            // auc or similar ranking.  Just choose the cutoff at which precision = 50%
            if( ranking_parity)
            {
              for( size_t i( 0); i < dataset_size; ++i)
              {
                if( EXPERIMENTAL( i)( j) >= threshold)
                {
                  column.PushBack( PREDICTED( i)( j));
                }
              }
            }
            else
            {
              for( size_t i( 0); i < dataset_size; ++i)
              {
                if( EXPERIMENTAL( i)( j) <= threshold)
                {
                  column.PushBack( PREDICTED( i)( j));
                }
              }
            }
            column.Sort( std::less< float>());

            // get the prediction cutoff
            prediction_cutoff = column( size_t( column.GetSize() * 0.5));
          }

          if( ranking_parity)
          {
            for( size_t i( 0); i < dataset_size; ++i)
            {
              const bool experimental_above_cutoff( EXPERIMENTAL( i)( j) >= threshold);
              if( ( PREDICTED( i)( j) >= prediction_cutoff) == experimental_above_cutoff)
              {
                // true
                prediction_classification( i, j) = experimental_above_cutoff ? 'P' : 'N';
              }
              else
              {
                // false
                prediction_classification( i, j) = experimental_above_cutoff ? 'n' : 'p';
              }
            }
          }
          else
          {
            for( size_t i( 0); i < dataset_size; ++i)
            {
              const bool experimental_below_cutoff( EXPERIMENTAL( i)( j) <= threshold);
              if( ( PREDICTED( i)( j) <= prediction_cutoff) == experimental_below_cutoff)
              {
                prediction_classification( i, j) = experimental_below_cutoff ? 'P' : 'N';
              }
              else
              {
                prediction_classification( i, j) = experimental_below_cutoff ? 'n' : 'p';
              }
            }
          }
        }
      }
      else if( GetGoalType() == e_Classification)
      {
        // for classification purposes, only use Pn
        for( size_t i( 0); i < dataset_size; ++i)
        {
          for( size_t j( 0); j < result_size; ++j)
          {
            prediction_classification( i, j)
              = char( ( EXPERIMENTAL( i)( j) >= threshold) == ( PREDICTED( i)( j) >= threshold) ? 'P' : 'n');
          }
        }
      }
      else
      {
        // for each column, compute the standard deviation of the deviations from experimental values
        // declare predictions incorrect that are more than 1 standard deviations from the mean deviation for this model
        for( size_t j( 0); j < result_size; ++j)
        {
          math::RunningAverageSD< float> result_ave_deviation;
          for( size_t i( 0); i < dataset_size; ++i)
          {
            if( util::IsDefined( EXPERIMENTAL( i)( j)))
            {
              result_ave_deviation += math::Absolute( PREDICTED( i)( j) - EXPERIMENTAL( i)( j));
            }
          }
          const float threshold_deviation( result_ave_deviation.GetAverage() + result_ave_deviation.GetStandardDeviation());
          for( size_t i( 0); i < dataset_size; ++i)
          {
            prediction_classification( i, j)
              = char( !util::IsDefined( EXPERIMENTAL( i)( j)) || math::Absolute( PREDICTED( i)( j) - EXPERIMENTAL( i)( j)) <= threshold_deviation ? 'P' : 'n');
          }
        }
      }
      return prediction_classification;
    }

  } // namespace model
} // namespace bcl
