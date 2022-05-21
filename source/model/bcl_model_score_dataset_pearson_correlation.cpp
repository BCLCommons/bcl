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
#include "model/bcl_model_score_dataset_pearson_correlation.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_running_average_sd.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ScoreDatasetPearsonCorrelation::s_Instance
    (
      util::Enumerated< ScoreDatasetInterface>::AddInstance( new ScoreDatasetPearsonCorrelation( false))
    );
    const util::SiPtr< const util::ObjectInterface> ScoreDatasetPearsonCorrelation::s_AutoInstance
    (
      util::Enumerated< ScoreDatasetInterface>::AddInstance( new ScoreDatasetPearsonCorrelation( true))
    );

    //! @brief Clone function
    //! @return pointer to new ScoreDatasetPearsonCorrelation
    ScoreDatasetPearsonCorrelation *ScoreDatasetPearsonCorrelation::Clone() const
    {
      return new ScoreDatasetPearsonCorrelation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ScoreDatasetPearsonCorrelation::GetClassIdentifier() const
    {
      return GetStaticClassName< ScoreDatasetPearsonCorrelation>();
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ScoreDatasetPearsonCorrelation::GetAlias() const
    {
      static const std::string s_Name( "PearsonCorrelation"), s_AutoName( "PearsonCorrelationNonRedundancy");
      return m_AutoCorrelation ? s_AutoName : s_Name;
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

    //! @brief score a given dataset
    //! @param DATASET dataset of interest
    //! @return scores of the dataset
    linal::Vector< float> ScoreDatasetPearsonCorrelation::Score( const descriptor::Dataset &DATASET) const
    {
      // determine # of features, results, and dataset size
      const size_t feature_size( DATASET.GetFeatureSize());
      const size_t result_size( m_AutoCorrelation ? feature_size : DATASET.GetResultSize());
      const size_t dataset_size( DATASET.GetSize());
      const linal::MatrixConstReference< float> features( DATASET.GetFeaturesReference());
      const linal::MatrixConstReference< float> results( m_AutoCorrelation ? features : DATASET.GetResultsReference());

      // handle empty dataset
      if( dataset_size == size_t( 0))
      {
        return linal::Vector< float>( feature_size, float( 0.0));
      }

      // compute the averages for each feature and result column
      math::RunningAverageSD< linal::Vector< float> > mean_feature, mean_result;
      // compute averages
      for( size_t counter( 0); counter < dataset_size; ++counter)
      {
        // average all features and results
        mean_feature += features.GetRow( counter);
      }
      // compute average results
      if( !m_AutoCorrelation)
      {
        for( size_t counter( 0); counter < dataset_size; ++counter)
        {
          mean_result += results.GetRow( counter);
        }
      }
      else
      {
        mean_result = mean_feature;
      }
      // for multiple result columns, return a score of the mean correlation for all feature_result pairings
      storage::Vector< storage::Vector< math::RunningAverage< float> > >
        mean_correlation( feature_size, storage::Vector< math::RunningAverage< float> >( result_size));

      linal::Vector< float>
        feature_temp( mean_feature.GetAverage()),
        result_temp( mean_result.GetAverage());

      // compute correlations
      for( size_t counter( 0); counter < dataset_size; ++counter)
      {
        // compute the difference from the mean
        feature_temp = features.GetRow( counter);
        result_temp = results.GetRow( counter);
        feature_temp -= mean_feature.GetAverage();
        result_temp -= mean_result.GetAverage();
        // BCL_Debug( counter);
        for( size_t feature_index( 0); feature_index < feature_size; ++feature_index)
        {
          for( size_t result_index( m_AutoCorrelation ? feature_index + 1 : 0); result_index < result_size; ++result_index)
          {
            mean_correlation( feature_index)( result_index) += feature_temp( feature_index) * result_temp( result_index);
          }
        }
      }

      // find the max absolute correlation for each feature
      linal::Vector< float> max_abs_correlations( feature_size, float( m_AutoCorrelation ? 0.0 : -1.0));
      for( size_t feature_index( 0); feature_index < feature_size; ++feature_index)
      {
        float max_abs_correlation( max_abs_correlations( feature_index));
        float std_feature( mean_feature.GetStandardDeviation()( feature_index));
        if( std_feature <= float( 0.0) || !util::IsDefined( std_feature))
        {
          if( m_AutoCorrelation)
          {
            max_abs_correlations( feature_index) = 1.0;
          }
          continue;
        }
        for( size_t result_index( m_AutoCorrelation ? feature_index + 1 : 0); result_index < result_size; ++result_index)
        {
          float std_result( mean_result.GetStandardDeviation()( result_index));
          if( std_result <= float( 0.0) || !util::IsDefined( std_result))
          {
            // skip results that did not vary
            continue;
          }
          float abs_correlation
          (
            math::Absolute( mean_correlation( feature_index)( result_index)) / ( std_feature * std_result)
          );

          // due to numeric truncation and roundoff, abs correlation may come out slightly above 1
          if( abs_correlation > 1.0)
          {
            abs_correlation = 1.0;
          }
          if( m_AutoCorrelation && abs_correlation > max_abs_correlations( result_index))
          {
            max_abs_correlations( result_index) = abs_correlation;
          }
          max_abs_correlation = std::max( abs_correlation, max_abs_correlation);
        }
        max_abs_correlations( feature_index) = max_abs_correlation;
      }
      if( m_AutoCorrelation)
      {
        max_abs_correlations *= -1.0;
        max_abs_correlations += float( 1.0);
      }

      // return f-scores for every feature column in dataset
      return max_abs_correlations;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ScoreDatasetPearsonCorrelation::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes absolute pearson correlation between features and results."
        "If multiple results are present, takes the maximum absolute correlation."
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
