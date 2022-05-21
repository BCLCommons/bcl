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
#include "model/bcl_model_score_dataset_f_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
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
    const util::SiPtr< const util::ObjectInterface> ScoreDatasetFScore::s_Instance
    (
      util::Enumerated< ScoreDatasetInterface>::AddInstance( new ScoreDatasetFScore())
    );

    //! @brief Clone function
    //! @return pointer to new ScoreDatasetFScore
    ScoreDatasetFScore *ScoreDatasetFScore::Clone() const
    {
      return new ScoreDatasetFScore( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ScoreDatasetFScore::GetClassIdentifier() const
    {
      return GetStaticClassName< ScoreDatasetFScore>();
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ScoreDatasetFScore::GetAlias() const
    {
      static const std::string s_Name( "FScore");
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

    //! @brief calculate f-score values for every column in combined active and inactive dataset
    //!                                     (X_avg_pos - X_avg)^2 + (X_avg_neg - X_avg)^2
    //! f(i)= -------------------------------------------------------------------------------------------------------
    //!       1/(n_pos - 1) SUM k=1 to n_pos (X_k_pos - X_avg)^2 + 1/(n_neg - 1) SUM k=1 to n_neg (X_k_neg - X_avg)^2
    //!
    //! @param FEATURES matrix of features
    //! @param RESULTS matrix of results
    //! @return scores of the dataset
    //! @brief score a given dataset
    //! @param DATASET dataset of interest
    //! @return scores of the dataset
    linal::Vector< float> ScoreDatasetFScore::Score( const descriptor::Dataset &DATASET) const
    {
      // determine # of features, results, and dataset size
      const size_t feature_size( DATASET.GetFeatureSize());
      const size_t result_size( DATASET.GetResultSize());
      const size_t dataset_size( DATASET.GetSize());
      const linal::MatrixConstReference< float> features( DATASET.GetFeaturesReference());
      const linal::MatrixConstReference< float> results( DATASET.GetResultsReference());

      // handle empty dataset
      if( dataset_size == size_t( 0))
      {
        return linal::Vector< float>( feature_size, float( 0.0));
      }

      // for multiple result columns, return a score of the mean f-score for all result columns
      math::RunningAverage< linal::Vector< float> > mean_fscores( linal::Vector< float>( result_size, 0.0));

      for( size_t result_index( 0); result_index < result_size; ++result_index)
      {
        // setup running average / variancess for the actives, inactives, and variance
        // initialize them all with a feature vector so that the running average will be for vectors of the right size
        math::RunningAverageSD< linal::Vector< float> > actives_stats, inactives_stats;

        // compute averages
        for( size_t counter( 0); counter < dataset_size; ++counter)
        {
          // check if the result is above the cutoff, e.g. nominally active
          if( results[ counter][ result_index] > m_Cutoff)
          {
            // add the feature to the averages for the actives
            actives_stats += linal::VectorConstReference< float>( feature_size, features[ counter]);
          }
          else
          {
            // add the feature to the averages for the inactives
            inactives_stats += linal::VectorConstReference< float>( feature_size, features[ counter]);
          }
        }

        if( actives_stats.GetWeight() == double( 0.0) || inactives_stats.GetWeight() == double( 0.0))
        {
          BCL_MessageCrt
          (
            "Ignoring result column " + util::Format()( result_index) + " because it contained no values "
            + ( actives_stats.GetWeight() == double( 0.0) ? " above " : " below ") + util::Format()( m_Cutoff)
          );
          continue;
        }

        math::RunningAverage< linal::Vector< float> > averages_total( actives_stats.GetAverage());

        averages_total.AddWeightedObservation( actives_stats.GetAverage(), actives_stats.GetWeight());
        averages_total.AddWeightedObservation( inactives_stats.GetAverage(), inactives_stats.GetWeight());

        // calculate parts of f-score, beginning with the numerator

        // f-score initialized with variance between the average actives and average totals
        linal::Vector< float> fscore( linal::SqrVector( actives_stats.GetAverage() - averages_total.GetAverage()));

        // add the variance between the averages inactives and averages total
        fscore += linal::SqrVector( inactives_stats.GetAverage() - averages_total.GetAverage());

        // n / (n - 1) * sample-variance = population variance
        //! @see http://en.wikipedia.org/wiki/Bessel%27s_correction
        // the fscore formula uses the population variance, but we have computed the sample variance, so we need to
        // multiply the values by these ratios, for the actives and inactives, respectively
        const float sample_population_variance_ratio_actives
        (
          actives_stats.GetWeight() / ( actives_stats.GetWeight() - 1.0)
        );
        const float sample_population_variance_ratio_inactives
        (
          inactives_stats.GetWeight() / ( inactives_stats.GetWeight() - 1.0)
        );

        // get an iterator over the fscore vector
        float *itr_fscore( fscore.Begin());

        // iterate over all standard deviations of all columns
        for
        (
          const float *itr_var_actives( actives_stats.GetVariance().Begin()),
                      *itr_var_inactives( inactives_stats.GetVariance().Begin()),
                      *itr_var_actives_end( actives_stats.GetVariance().End());
          itr_var_actives != itr_var_actives_end;
          ++itr_var_actives, ++itr_var_inactives, ++itr_fscore
        )
        {
          // denominator
          float denominator( *itr_var_actives * sample_population_variance_ratio_actives);
          denominator += *itr_var_inactives * sample_population_variance_ratio_inactives;
          *itr_fscore /= denominator;
          if( !util::IsDefined( *itr_fscore))
          {
            *itr_fscore = 0.0;
          }
        }
        mean_fscores += fscore;
      }

      // return f-scores for every feature column in dataset
      return mean_fscores.GetAverage();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ScoreDatasetFScore::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "uses fscore, see Yi-Wei Chen and Chih-Jen Lin. Combining SVMs with various feature selection strategies."
        "In Feature extraction, foundations and applications. Springer, 2006."
      );

      parameters.AddInitializer
      (
        "cutoff",
        "result value that separates actives from inactives",
        io::Serialization::GetAgent( &m_Cutoff),
        "0"
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
