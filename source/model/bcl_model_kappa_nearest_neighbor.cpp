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
#include "model/bcl_model_kappa_nearest_neighbor.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_running_average.h"
#include "storage/bcl_storage_triplet.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    //! the default input range
    const math::Range< float> KappaNearestNeighbor::s_DefaultInputRange( 0, 1);

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> KappaNearestNeighbor::s_Instance
    (
      GetObjectInstances().AddInstance( new KappaNearestNeighbor())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    KappaNearestNeighbor::KappaNearestNeighbor() :
      m_TrainingData(),
      m_Kappa( 1)
    {
    }

    //! @brief constructor from training data, query data, and rescale functions
    //! @param TRAINING_DATA data to train the NeuralNetwork on
    //! @param KAPPA kappa value for number of nearest neighbors to consider
    KappaNearestNeighbor::KappaNearestNeighbor
    (
      const util::ShPtr< descriptor::Dataset> &TRAINING_DATA,
      const size_t KAPPA
    ) :
      m_TrainingData( TRAINING_DATA),
      m_Kappa( KAPPA)
    {
    }

    //! @brief copy constructor
    KappaNearestNeighbor *KappaNearestNeighbor::Clone() const
    {
      return new KappaNearestNeighbor( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &KappaNearestNeighbor::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &KappaNearestNeighbor::GetAlias() const
    {
      static const std::string s_Name( "KappaNearestNeighbor");
      return s_Name;
    }

    //! @brief Set the scaling of a feature set according to the model
    //! @param FEATURES feature set of interest
    //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
    //!       when operator() is called
    void KappaNearestNeighbor::Rescale( FeatureDataSet< float> &FEATURE) const
    {
      if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_TrainingData->GetFeaturesPtr()->GetScaling())
      {
        FEATURE.DeScale();
        FEATURE.Rescale( *m_TrainingData->GetFeaturesPtr()->GetScaling());
      }
    }

    //! @brief performs the kappa nearest neighbor algorithm
    //! @param QUERY_VECTOR math::Vector< float> of features which to the neighbors will be found
    //! @param KAPPA number of neigbors
    //! @param TRAINING_FEATURES that from the kappa nearest neighbors are selected
    //! @return StorageSet of kappa nearest neighbors and the corresponding distance from QUERY_VECTOR
    storage::Vector< storage::Pair< float, size_t> > KappaNearestNeighbor::FindWithoutRescaling
    (
      const FeatureReference< float> &QUERY_VECTOR,
      const size_t KAPPA,
      const FeatureDataSet< float> &TRAINING_FEATURES
    )
    {
      BCL_Assert( TRAINING_FEATURES.GetNumberFeatures() > 0, "kNN requires training data!");

      // storing size of data set
      const size_t training_size( TRAINING_FEATURES.GetNumberFeatures());

      BCL_Message( util::Message::e_Debug, "training data size " + util::Format()( training_size));

      // make a set with keys being the kappa lowest distances and place in the training data from which they came
      // If we only used the distance as the key in the map, then if two training data had the same distance from they
      // query vector, the second one would overwrite the value of the first.  With the indice as the second value,
      // we will keep both values in the set
      storage::Set< storage::Pair< float, size_t> >
        kappa_nearest_neighbor_distance;

      // keep track of the kappa-th smallest distance seen so far
      float cutoff_distance( std::numeric_limits< float>::infinity());

      // store the end of the query vector
      const float *itr_query_end( QUERY_VECTOR.End());

      // loop over reference data
      for( size_t itr( 0); itr < training_size; ++itr)
      {
        float distance( 0.0);
        // store norm of difference vector
        for
        (
          const float *itr_data( TRAINING_FEATURES( itr).Begin()), *itr_query( QUERY_VECTOR.Begin());
          distance < cutoff_distance && itr_query != itr_query_end;
          ++itr_data, ++itr_query
        )
        {
          distance += math::Sqr( *itr_data - *itr_query);
        }

        // if the new distance is >= to the cutoff distance, then the output value won't be considered, so continue
        if( distance >= cutoff_distance)
        {
          continue;
        }

        distance = math::Sqrt( distance);

        // if the map has fewer distances than KAPPA, go ahead and add the distance to the map
        kappa_nearest_neighbor_distance.Insert( storage::Pair< float, size_t>( distance, itr));

        // if the map is now larger than kappa, remove the last element
        if( kappa_nearest_neighbor_distance.GetSize() > KAPPA)
        {
          storage::Set< storage::Pair< float, size_t> >::iterator
            itr_end( kappa_nearest_neighbor_distance.End()), new_itr_end;

          // go to the last element
          new_itr_end = --itr_end;

          // new_itr_end needs to go one forward
          --new_itr_end;

          // remove the last element
          kappa_nearest_neighbor_distance.RemoveElement( itr_end);

          // store the new cutoff distance, which is the distance of the last element in the map
          cutoff_distance = math::Sqr( new_itr_end->First());
        }
      }
      //copy into vector
      return storage::Vector< storage::Pair< float, size_t> >( kappa_nearest_neighbor_distance.Begin(), kappa_nearest_neighbor_distance.End());
    }

    //! @brief performs the kappa nearest neighbor algorithm
    //! @param QUERY_VECTORS math::Vector< float> of features
    //!        and known output of query point
    //! @param RESULT_STORAGE pointer to vector or matrix where result should be stored; must be preallocated
    storage::Vector< storage::Triplet< float, size_t, size_t> > KappaNearestNeighbor::FindWithoutRescaling
    (
      const FeatureDataSet< float> &QUERY_VECTORS,
      const size_t KAPPA,
      const FeatureDataSet< float> &TRAINING_FEATURES
    )
    {
      BCL_Assert( TRAINING_FEATURES.GetNumberFeatures() > 0, "kNN requires training data!");

      // storing size of data set
      const size_t training_size( TRAINING_FEATURES.GetNumberFeatures());

      BCL_Message( util::Message::e_Debug, "training data size " + util::Format()( training_size));

      // make a set with keys being the kappa lowest distances and place in the training data from which they came
      // If we only used the distance as the key in the map, then if two training data had the same distance from they
      // query vector, the second one would overwrite the value of the first.  With the indice as the second value,
      // we will keep both values in the set
      storage::Set< storage::Triplet< float, size_t, size_t> >
        kappa_nearest_neighbor_distance;

      // keep track of the kappa-th smallest distance seen so far
      float cutoff_distance( std::numeric_limits< float>::infinity());
      size_t current_size( 0);
      for( size_t query_nr( 0), query_sz( QUERY_VECTORS.GetNumberFeatures()); query_nr < query_sz; ++query_nr)
      {
        // store the end of the query vector
        const float *itr_query_end( QUERY_VECTORS( query_nr).End());

        // loop over reference data
        for( size_t itr( 0); itr < training_size; ++itr)
        {
          float distance( 0.0);
          // store norm of difference vector
          for
          (
            const float *itr_data( TRAINING_FEATURES[ itr]), *itr_query( QUERY_VECTORS[ query_nr]);
            distance < cutoff_distance && itr_query != itr_query_end;
            ++itr_data, ++itr_query
          )
          {
            distance += math::Sqr( *itr_data - *itr_query);
          }

          // if the new distance is >= to the cutoff distance, then the output value won't be considered, so continue
          if( distance >= cutoff_distance)
          {
            continue;
          }

          distance = math::Sqrt( distance);

          // if the map has fewer distances than KAPPA, go ahead and add the distance to the map
          kappa_nearest_neighbor_distance.Insert( storage::Triplet< float, size_t, size_t>( distance, query_nr, itr));

          // if the map is now larger than kappa, remove the last element
          if( current_size + 1 > KAPPA)
          {
            storage::Set< storage::Triplet< float, size_t, size_t> >::iterator
              itr_end( kappa_nearest_neighbor_distance.End()), new_itr_end;

            // go to the last element
            new_itr_end = --itr_end;

            // new_itr_end needs to go one forward
            --new_itr_end;

            // remove the last element
            kappa_nearest_neighbor_distance.RemoveElement( itr_end);

            // store the new cutoff distance, which is the distance of the last element in the map
            cutoff_distance = math::Sqr( new_itr_end->First());
          }
          else
          {
            ++current_size;
          }
        }
      }
      //copy into vector
      return storage::Vector< storage::Triplet< float, size_t, size_t> >( kappa_nearest_neighbor_distance.Begin(), kappa_nearest_neighbor_distance.End());
    }

  //////////////
  // operator //
  //////////////

    //! @brief predict result with model using a NOT rescaled feature vector
    //! @param FEATURE not rescaled feature vector
    //! @return predicted result vector using a model
    FeatureDataSet< float> KappaNearestNeighbor::PredictWithoutRescaling( const FeatureDataSetInterface< float> &FEATURE) const
    {
      const size_t num_features( FEATURE.GetNumberFeatures());
      linal::Matrix< float> predicted_results( num_features, m_TrainingData->GetResultSize());
      for( size_t feature_number( 0); feature_number < num_features; ++feature_number)
      {
        EuclideanDistance( FEATURE( feature_number), predicted_results[ feature_number]);
      }

      return FeatureDataSet< float>( predicted_results);
    }

    //! @brief predict result with model using a rescaled feature vector
    //! @param FEATURE normalized or rescaled feature vector
    //! @return predicted result vector using a model
    FeatureDataSet< float> KappaNearestNeighbor::operator()( const FeatureDataSetInterface< float> &FEATURE) const
    {
      // handle the case where rescaling is necessary
      if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_TrainingData->GetFeaturesPtr()->GetScaling())
      {
        FeatureDataSet< float> feature( FEATURE);
        feature.Rescale( *m_TrainingData->GetFeaturesPtr()->GetScaling());
        return PredictWithoutRescaling( feature).DeScale();
      }

      // data is already rescaled
      return PredictWithoutRescaling( FEATURE).DeScale();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read KappaNearestNeighbor from std::istream
    std::istream &KappaNearestNeighbor::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_TrainingData,  ISTREAM);
      io::Serialize::Read( m_Kappa,         ISTREAM);

      // end
      return ISTREAM;
    }

    //! write KappaNearestNeighbor into std::ostream
    std::ostream &KappaNearestNeighbor::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_TrainingData,  OSTREAM) << '\n';
      io::Serialize::Write( m_Kappa,         OSTREAM);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief performs the kappa nearest neighbor algorithm
    //! @param QUERY_VECTOR linal::Vector< float> of features
    //!        and known output of query point
    //! @param RESULT_STORAGE pointer to vector or matrix where result should be stored; must be preallocated
    void KappaNearestNeighbor::EuclideanDistance
    (
      const FeatureReference< float> &QUERY_VECTOR,
      float *RESULT_STORAGE
    ) const
    {
      BCL_Assert( m_TrainingData->GetSize() > 0, "kNN requires training data!");

      // storing size of data set
      const size_t training_size( m_TrainingData->GetSize());
      const FeatureDataSetInterface< float> &training_features( *m_TrainingData->GetFeaturesPtr());

      BCL_MessageDbg( "training data size " + util::Format()( training_size));

      // make a set with keys being the kappa lowest distances and place in the training data from which they came
      // If we only used the distance as the key in the map, then if two training data had the same distance from they
      // query vector, the second one would overwrite the value of the first.  With the indice as the second value,
      // we will keep both values in the set
      storage::Set< storage::Pair< float, size_t> >
        kappa_nearest_neighbor_distance;

      // keep track of the kappa-th smallest distance seen so far
      float cutoff_distance( std::numeric_limits< float>::infinity());

      // store the end of the query vector
      const float *itr_query_end( QUERY_VECTOR.End());

      // loop over reference data
      for( size_t itr( 0); itr < training_size; ++itr)
      {
        float distance( 0.0);
        // store norm of difference vector
        for
        (
          const float *itr_data( training_features( itr).Begin()), *itr_query( QUERY_VECTOR.Begin());
          distance < cutoff_distance && itr_query != itr_query_end;
          ++itr_data, ++itr_query
        )
        {
          distance += math::Sqr( *itr_data - *itr_query);
        }

        // if the new distance is >= to the cutoff distance, then the output value won't be considered, so continue
        if( distance >= cutoff_distance)
        {
          continue;
        }

        distance = math::Sqrt( distance);

        // if the map has fewer distances than m_Kappa, go ahead and add the distance to the map
        kappa_nearest_neighbor_distance.Insert( storage::Pair< float, size_t>( distance, itr));

        // if the map is now larger than kappa, remove the last element
        if( kappa_nearest_neighbor_distance.GetSize() > m_Kappa)
        {
          storage::Set< storage::Pair< float, size_t> >::iterator
            itr_end( kappa_nearest_neighbor_distance.End()), new_itr_end;

          // go to the last element
          new_itr_end = --itr_end;

          // new_itr_end needs to go one forward
          --new_itr_end;

          // remove the last element
          kappa_nearest_neighbor_distance.RemoveElement( itr_end);

          // store the new cutoff distance, which is the distance of the last element in the map
          cutoff_distance = math::Sqr( new_itr_end->First());
        }
      }

      // setup activities sum
      const FeatureDataSetInterface< float> &training_results( *m_TrainingData->GetResultsPtr());
      math::RunningAverage< linal::Vector< float> > weighted_ave_activity( training_results( 0));

      // check for distances of zero; if zero distances are present, average them and return
      // to avoid nearly - zero distances that result in overflow, consider any distance < 1e-16 to be equivalent
      storage::Set< storage::Pair< float, size_t> >::const_iterator
        itr_act( kappa_nearest_neighbor_distance.Begin()),
        itr_act_end( kappa_nearest_neighbor_distance.End());

      if( itr_act->First() < float( 1.0e-16))
      {
        for( ; itr_act != itr_act_end && itr_act->First() < float( 1.0e-16); ++itr_act)
        {
          weighted_ave_activity += training_results( itr_act->Second());
        }
      }
      else
      {
        // no distances of essentially zero, so add up the kappa nearest neighbors, weighted by distance
        for( ; itr_act != itr_act_end; ++itr_act)
        {
          // weight each output by distance ^ -1
          const float weight( 1.0 / itr_act->First());

          // add the weighted activity to the weighted activities sum
          weighted_ave_activity.AddWeightedObservation( training_results( itr_act->Second()), weight);
        }
      }

      std::copy( weighted_ave_activity.GetAverage().Begin(), weighted_ave_activity.GetAverage().End(), RESULT_STORAGE);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer KappaNearestNeighbor::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "see http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm"
      );

      parameters.AddInitializer
      (
        "kappa",
        "number of nearest neighbors to use for computing average",
        io::Serialization::GetAgentWithMin( &m_Kappa, size_t( 1)),
        "3"
      );

      return parameters;
    }

  } // namespace model
} // namespace bcl
