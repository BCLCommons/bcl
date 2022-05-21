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
#include "model/bcl_model_score_dataset_non_redundant.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
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
    const util::SiPtr< const util::ObjectInterface> ScoreDatasetNonRedundant::s_Instance
    (
      util::Enumerated< ScoreDatasetInterface>::AddInstance( new ScoreDatasetNonRedundant())
    );

    //! @brief default constructor
    ScoreDatasetNonRedundant::ScoreDatasetNonRedundant() :
      m_ZScoreTolerance( 0.1),
      m_MaxOutliers( 3),
      m_AllowOneConstant( false),
      m_MinSpan( std::numeric_limits< float>::epsilon()),
      m_MinStd( std::numeric_limits< float>::epsilon()),
      m_MinAbsCrossCorrelation( 0.95)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ScoreDatasetNonRedundant
    ScoreDatasetNonRedundant *ScoreDatasetNonRedundant::Clone() const
    {
      return new ScoreDatasetNonRedundant( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ScoreDatasetNonRedundant::GetClassIdentifier() const
    {
      return GetStaticClassName< ScoreDatasetNonRedundant>();
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ScoreDatasetNonRedundant::GetAlias() const
    {
      static const std::string s_Name( "NonRedundant");
      return s_Name;
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
    linal::Vector< float> ScoreDatasetNonRedundant::Score( const descriptor::Dataset &DATASET) const
    {
      // determine # of features, results, and dataset size
      const size_t feature_size( DATASET.GetFeatureSize());
      const size_t dataset_size( DATASET.GetSize());
      const linal::MatrixConstReference< float> features( DATASET.GetFeaturesReference());

      linal::Vector< float> is_non_redundant( feature_size, float( 1.0));
      // handle empty dataset
      if( dataset_size == size_t( 0) || feature_size <= size_t( 1))
      {
        return is_non_redundant;
      }

      // create rescale object, z-score and for min-max
      // This is an easy way of getting the min/max/ave/std for each column
      RescaleFeatureDataSet rescaler_ave_std
      (
        features,
        math::Range< float>( 0.0, 1.0),
        RescaleFeatureDataSet::e_AveStd
      );
      RescaleFeatureDataSet rescaler_min_max
      (
        features,
        math::Range< float>( 0.0, 1.0),
        RescaleFeatureDataSet::e_MinMax
      );

      // Get the feature labels for display to the user
      const FeatureLabelSet feature_labels( DATASET.GetFeaturesPtr()->GetFeatureLabelSet()->SplitFeatureLabelSet( true));

      // first, handle constant descriptors. Only one constant descriptor is allowed, and then only if its constant
      // value is non-zero
      bool have_constant_non_zero_descriptor( false);
      size_t constant_feature_index( util::GetUndefined< size_t>());
      for( size_t feature_number( 0); feature_number < feature_size; ++feature_number)
      {
        const math::Range< float> &range( rescaler_min_max.GetRescaleRanges()( feature_number));
        const float std_span( rescaler_ave_std.GetRescaleRanges()( feature_number).GetWidth() * 0.25);

        if( range.GetWidth() < m_MinSpan || std_span < m_MinStd)
        {
          // constant feature
          is_non_redundant( feature_number) = 0.0;
          if
          (
            m_AllowOneConstant && !have_constant_non_zero_descriptor
            && math::Absolute( range.GetMiddle()) > 0.01
            && range.GetMin() == range.GetMax()
          )
          {
            // first non-zero, constant feature; it'll still have a 0 score to simplify the checks
            // in the main redundancy check loop. Afterwards, it'll be set back to 1.0
            have_constant_non_zero_descriptor = true;
            constant_feature_index = feature_number;
          }
          else
          {
            BCL_MessageStd
            (
              feature_labels.GetMemberLabels()( feature_number).ToString()
              + " is a constant " + util::Format()( rescaler_ave_std.GetRescaleRanges()( feature_number).GetMiddle())
            );
          }
        }
      }

      // for each feature
      const math::Range< float> rescale_to( -1.0, 1.0);
      for( size_t feature_number( 1); feature_number < feature_size; ++feature_number)
      {
        if( is_non_redundant( feature_number) == size_t( 0))
        {
          // skip constant descriptors
          continue;
        }

        // get the relevant range parameters for this descriptor
        const math::Range< float> ave_std_range_a( rescaler_ave_std.GetRescaleRanges()( feature_number));

        float range_a( ave_std_range_a.GetWidth());

        // construct a linear function to rescale this descriptor
        const math::LinearFunction a_to_zscore( ave_std_range_a, rescale_to);

        // for each other descriptor before this one that was non-redundant
        for( size_t feature_number_b( 0); feature_number_b < feature_number; ++feature_number_b)
        {
          if( is_non_redundant( feature_number_b) == size_t( 0))
          {
            // skip constant descriptors
            continue;
          }

          // get the relevant range parameters for this descriptor
          const math::Range< float> ave_std_range_b( rescaler_ave_std.GetRescaleRanges()( feature_number_b));

          // test Std(A)/Range(A) = Std(B)/Range(B)
          float range_b( ave_std_range_b.GetWidth());

          // construct a linear function to rescale this descriptor
          const math::LinearFunction b_to_zscore( ave_std_range_b, rescale_to);

          size_t n_outliers_pos( 0), n_outliers_neg( 0);
          math::RunningAverage< float> average_product;
          // run through dataset, find whether the columns disagree by more than zscore_tolerance at any point
          for( size_t row_index( 0); row_index < dataset_size; ++row_index)
          {
            const float *row( features[ row_index]);
            const float a_zscore( a_to_zscore( row[ feature_number]));
            const float b_zscore( b_to_zscore( row[ feature_number_b]));
            if( n_outliers_pos <= m_MaxOutliers && math::Absolute( a_zscore - b_zscore) > m_ZScoreTolerance)
            {
              ++n_outliers_pos;
              if( n_outliers_pos > m_MaxOutliers && n_outliers_neg > m_MaxOutliers)
              {
                break;
              }
            }
            if( n_outliers_neg <= m_MaxOutliers && math::Absolute( a_zscore + b_zscore) > m_ZScoreTolerance)
            {
              ++n_outliers_neg;
              if( n_outliers_pos > m_MaxOutliers && n_outliers_neg > m_MaxOutliers)
              {
                break;
              }
            }
            average_product += a_zscore * b_zscore;
          }
          const float rsq( std::min( math::Absolute( 4.0 * average_product.GetAverage()), 1.0));
          if( ( n_outliers_pos <= m_MaxOutliers || n_outliers_neg <= m_MaxOutliers) && rsq >= m_MinAbsCrossCorrelation)
          {
            size_t redundant( feature_number), kept( feature_number_b);
            if
            (
              feature_labels.GetMemberLabels()( kept).GetValue() != "Partial"
              && feature_labels.GetMemberLabels()( redundant).GetValue() == "Partial"
            )
            {
              // prefer to keep the partial over a scalar
              // this results in fewer descriptors being calculated
              std::swap( redundant, kept);
            }
            if( ave_std_range_a.GetMin() >= 0.0 && ave_std_range_b.GetMin() >= 0.0)
            {
              range_b = ave_std_range_b.GetMiddle();
              range_a = ave_std_range_a.GetMiddle();
            }
            BCL_MessageStd
            (
              feature_labels.GetMemberLabels()( redundant).ToString()
              + " is redundant with " + feature_labels.GetMemberLabels()( kept).ToString()
              + std::string( n_outliers_pos < n_outliers_neg ? " Multiplier +" : " Multiplier -")
              + util::Format()( redundant != feature_number ? range_b / range_a : range_a / range_b)
              + " Cross correlation: " + util::Format()( rsq)
            );
            is_non_redundant( redundant) = 0.0;
            break;
          }
        }
      }

      if( have_constant_non_zero_descriptor)
      {
        is_non_redundant( constant_feature_index) = 1.0;
      }

      // return f-scores for every feature column in dataset
      return is_non_redundant;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ScoreDatasetNonRedundant::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "1 for each descriptor that is not a linear function of a previously seen descriptor, 0 otherwise. "
        "Also returns 0 for descriptors that are always 0"
      );
      parameters.AddInitializer
      (
        "max outliers",
        "Maximum number of data points that can violate the z-score tolerance for a particular descriptor before declaring "
        "a pair of descriptors non-redundant",
        io::Serialization::GetAgent( &m_MaxOutliers),
        "3"
      );
      parameters.AddInitializer
      (
        "tol",
        "Maximum difference in z-score between two redundant descriptors",
        io::Serialization::GetAgent( &m_ZScoreTolerance),
        "0.1"
      );
      parameters.AddInitializer
      (
        "allow a constant",
        "allow up to one constant, non-zero, descriptor in a dataset. "
        "This is useful for keeping an offset for linear regression",
        io::Serialization::GetAgent( &m_AllowOneConstant),
        "False"
      );
      parameters.AddInitializer
      (
        "min span",
        "Minimum range of values that a descriptor needs to span to be considered non-constant",
        io::Serialization::GetAgent( &m_MinSpan),
        util::Format()( std::numeric_limits< float>::epsilon())
      );
      parameters.AddInitializer
      (
        "min std",
        "Minimum standard deviation of values that a descriptor needs to span to be considered non-constant",
        io::Serialization::GetAgent( &m_MinStd),
        util::Format()( std::numeric_limits< float>::epsilon())
      );
      parameters.AddInitializer
      (
        "min rsq",
        "Minimum r-squared for declaration of feature as redundant",
        io::Serialization::GetAgent( &m_MinAbsCrossCorrelation),
        "0.95"
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
