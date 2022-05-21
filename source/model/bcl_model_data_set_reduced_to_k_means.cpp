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
#include "model/bcl_model_data_set_reduced_to_k_means.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_running_average.h"
#include "model/bcl_model_training_schedule.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> DataSetReducedToKMeans::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new DataSetReducedToKMeans())
    );

    //! @brief default constructor
    //! @param NUMBER_CLUSTERS the desired number of clusters
    //! @param MAX_RECLUSTER_ATTEMPTS the maximum number of times to try improving the clusters
    DataSetReducedToKMeans::DataSetReducedToKMeans
    (
      const size_t &NUMBER_CLUSTERS,
      const size_t &MAX_RECLUSTER_ATTEMPTS,
      const bool &AUTOSCALE
    ) :
      m_NumberClusters( NUMBER_CLUSTERS),
      m_MaxReclusterAttempts( MAX_RECLUSTER_ATTEMPTS),
      m_Autoscale( AUTOSCALE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new DataSetReducedToKMeans
    DataSetReducedToKMeans *DataSetReducedToKMeans::Clone() const
    {
      return new DataSetReducedToKMeans( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &DataSetReducedToKMeans::GetClassIdentifier() const
    {
      return GetStaticClassName< DataSetReducedToKMeans>();
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &DataSetReducedToKMeans::GetAlias() const
    {
      static const std::string s_alias( "KMeans");
      return s_alias;
    }

    //! @brief return the maximum number of clusters that will be returned
    //! @return the maximum number of clusters that will be returned
    //! the actual number of clusters will be no larger than the size of the data set
    size_t DataSetReducedToKMeans::GetNumberClusters() const
    {
      return m_NumberClusters;
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    //! @return the code / label for the feature (1st part) of the data set
    void DataSetReducedToKMeans::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectFeatures( CODE);
      m_Retriever->SelectFeatures( CODE);
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @return the code / label for the result (2nd part) of the data set
    void DataSetReducedToKMeans::SelectResults( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectResults( CODE);
      m_Retriever->SelectResults( CODE);
    }

    //! @brief Set which id columns to retrieve
    //! @param CODE the id column names to retrieve
    void DataSetReducedToKMeans::SelectIds( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectIds( CODE);
      m_Retriever->SelectIds( CODE);
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet DataSetReducedToKMeans::GetFeatureLabelsWithSizes() const
    {
      return m_Retriever->GetFeatureLabelsWithSizes();
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet DataSetReducedToKMeans::GetResultCodeWithSizes() const
    {
      return m_Retriever->GetResultCodeWithSizes();
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet DataSetReducedToKMeans::GetIdCodeWithSizes() const
    {
      return m_Retriever->GetIdCodeWithSizes();
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > DataSetReducedToKMeans::GetNumberPartitionsAndIds() const
    {
      return m_Retriever->GetNumberPartitionsAndIds();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief compute the square of the distance between two vectors, so long as it is no greater than LIMIT
    //! @param VEC_A, VEC_B the two vectors
    //! @param LIMIT the maximum distance of interest
    //! @return the distance, so long as it is less than LIMIT, otherwise, a value larger than LIMIT
    float DataSetReducedToKMeans::LimitedSquareDistance
    (
      const linal::VectorConstInterface< float> &VEC_A,
      const linal::VectorConstInterface< float> &VEC_B,
      const float &LIMIT
    )
    {
      BCL_Assert( VEC_A.GetSize() == VEC_B.GetSize(), "Vectors must be of the same size to compute the distance");

      // get the beginning of each vector
      const float *vec_a_begin( VEC_A.Begin()), *vec_b_begin( VEC_B.Begin());

      // compute the distance, but stop if limit is reached
      float distance_squared( 0.0);
      for( size_t index( 0), size( VEC_A.GetSize()); distance_squared < LIMIT && index < size; ++index)
      {
        distance_squared += math::Sqr( vec_a_begin[ index] - vec_b_begin[ index]);
      }

      return distance_squared;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief generate dataset, reduced to the desired # of cluster centers
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset>
      DataSetReducedToKMeans::GenerateDataSet()
    {
      return operator ()( m_Retriever->GenerateDataSet());
    }

    //! @brief reduce a data set
    //! @param DATA the data set
    //! @return the reduced data set
    util::ShPtr< descriptor::Dataset>
      DataSetReducedToKMeans::operator ()
      (
        const util::ShPtr< descriptor::Dataset> &DATA,
        const TrainingSchedule &BALANCING
      ) const
    {
      // make references to the feature and result matrices
      const linal::MatrixConstReference< float> features( DATA->GetFeaturesReference());
      const linal::MatrixConstReference< float> results( DATA->GetResultsReference());

      // cache some commonly needed numbers
      const size_t data_set_size( features.GetNumberRows());
      const size_t feature_size( features.GetNumberCols());
      const size_t result_size( results.GetNumberCols());

      // there are fewer features than m_NumberClusters, then just return a copy of the data set
      if( data_set_size <= m_NumberClusters)
      {
        return DATA;
      }

      TrainingSchedule balancer( BALANCING);
      // initialize the balancer w/ dummy balancing if no balancing was specified
      if( !balancer.GetSize())
      {
        balancer.Setup( *DATA->GetResultsPtr(), 0.0);
      }

      // make objects to rescale the input/output, and rescaled data, if desired
      util::ShPtr< linal::Matrix< float> > sp_feature_matrix;

      // make objects to rescale the input/output, and rescaled data, if desired
      RescaleFeatureDataSet rescale_feature_result;

      // make an object to hold si-ptrs to all the feature-results
      storage::Vector< storage::VectorND< 2, linal::VectorConstReference< float> > > feature_results( data_set_size);

      // perform the scaling functions
      if( m_Autoscale)
      {
        sp_feature_matrix = util::CloneToShPtr( DATA->GetFeaturesPtr()->GetRawMatrix());

        linal::Matrix< float> &features_rescaled( *sp_feature_matrix);

        rescale_feature_result = RescaleFeatureDataSet( *sp_feature_matrix, math::Range< float>( 0, 1));
        rescale_feature_result.RescaleMatrix( features_rescaled);

        for( size_t index( 0); index < data_set_size; ++index)
        {
          feature_results( index) =
            storage::VectorND< 2, linal::VectorConstReference< float> >
            (
              linal::VectorConstReference< float>( feature_size, features_rescaled[ index]),
              linal::VectorConstReference< float>( result_size, results[ index])
            );
        }
      }
      else
      {
        // just add the feature-results to the si-ptr vector
        for( size_t index( 0); index < data_set_size; ++index)
        {
          feature_results( index) =
            storage::VectorND< 2, linal::VectorConstReference< float> >
            (
              linal::VectorConstReference< float>( feature_size, features[ index]),
              linal::VectorConstReference< float>( result_size, results[ index])
            );
        }
      }

      // initialize the cluster index vector, which stores the cluster id of each feature
      storage::Vector< size_t> cluster_indices( data_set_size, util::GetUndefined< size_t>());

      // initialize the clusters
      storage::Vector< math::RunningAverage< linal::Vector< float> > > clusters
      (
        m_NumberClusters,
        math::RunningAverage< linal::Vector< float> >( linal::Vector< float>( feature_size, 0.0))
      );

      // initialize the clusters with the first m_NumberClusters features
      for( size_t i( 0); i < m_NumberClusters; ++i)
      {
        clusters( i) += feature_results( i).First();

        // set the cluster indices up to the associated cluster, since we know that these are the closest values
        // this isn't necessary, but it speeds up the search for the closest node later on
        cluster_indices( i) = i;
      }

      // refine clusters for m_MaxReclusterAttempts or until they do not change, whichever comes first
      size_t clustering_attempt_number( 0);
      while
      (
        clustering_attempt_number < m_MaxReclusterAttempts
        && RefineClusters( feature_results, cluster_indices, clusters, balancer)
      )
      {
        ++clustering_attempt_number;
      }
      BCL_MessageVrb( "# reclusterings: " + util::Format()( clustering_attempt_number));

      // make an object to hold the clusters that are found
      // There is no sensical way to ID the clusters, so leave the ID columns empty
      util::ShPtr< descriptor::Dataset> sp_clusters;
      if
      (
        DATA->GetFeaturesPtr()->GetFeatureLabelSet().IsDefined()
        && DATA->GetResultsPtr()->GetFeatureLabelSet().IsDefined()
        && DATA->GetIdsPtr()->GetFeatureLabelSet().IsDefined()
      )
      {
        sp_clusters =
          util::ShPtr< descriptor::Dataset>
          (
            new descriptor::Dataset
            (
              m_NumberClusters,
              *DATA->GetFeaturesPtr()->GetFeatureLabelSet(),
              *DATA->GetResultsPtr()->GetFeatureLabelSet(),
              *DATA->GetIdsPtr()->GetFeatureLabelSet()
            )
          );
      }
      else
      {
        // initialize using just the basic sizes
        sp_clusters =
          util::ShPtr< descriptor::Dataset>
          (
            new descriptor::Dataset
            (
              m_NumberClusters,
              feature_size,
              result_size,
              DATA->GetIdSize()
            )
          );
      }
      linal::MatrixReference< float> clustered_features( sp_clusters->GetFeaturesReference());
      linal::MatrixReference< float> clustered_results( sp_clusters->GetResultsReference());

        // just get the results
      for( size_t itr( 0); itr < m_NumberClusters; ++itr)
      {
        clustered_features.ReplaceRow( itr, clusters( itr).GetAverage());
      }

      // descale if necessary
      if( m_Autoscale)
      {
        rescale_feature_result.DeScaleMatrix( clustered_features);
      }

      // reset the clusters so they can be used to average the results
      for( size_t itr( 0); itr < m_NumberClusters; ++itr)
      {
        clusters( itr) = math::RunningAverage< linal::Vector< float> >( linal::Vector< float>( result_size, 0.0));
      }

      // now get the average results for each cluster
      for( size_t itr( 0); itr < data_set_size; ++itr)
      {
        clusters( cluster_indices( itr)) += feature_results( itr).Second();
      }

      // add the results to the clustered_features_results
      for( size_t i( 0); i < m_NumberClusters; ++i)
      {
        clustered_results.ReplaceRow( i, clusters( i).GetAverage());
      }

      // return a new data set initialized with the cluster centers
      return sp_clusters;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer DataSetReducedToKMeans::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription
      (
        "Create a dataset using the k-means average of another dataset"
      );
      member_data.AddInitializer
      (
        "dataset",
        "dataset retriever to call to get the entire data set",
        io::Serialization::GetAgent( &m_Retriever)
      );
      member_data.AddInitializer
      (
        "clusters",
        "number of clusters for k-means clustering",
        io::Serialization::GetAgentWithMin( &m_NumberClusters, size_t( 1))
      );

      member_data.AddInitializer
      (
        "max steps",
        "maximum number of times to refine the clusters before accepting the solution",
        io::Serialization::GetAgentWithMin( &m_MaxReclusterAttempts, size_t( 1)),
        "10"
      );

      member_data.AddInitializer
      (
        "autoscale",
        "whether to normalize / denormalize the data set before / after clustering",
        io::Serialization::GetAgent( &m_Autoscale),
        "0"
      );

      return member_data;
    }

    storage::Vector< size_t> DataSetReducedToKMeans::IdentifyClusterIndices( util::ShPtr< descriptor::Dataset> &REDUCED_MATRIX) const
    {
      // make references to the feature and result matrices
      const linal::MatrixConstReference< float> features( REDUCED_MATRIX->GetFeaturesReference());
      const linal::MatrixConstReference< float> results( REDUCED_MATRIX->GetResultsReference());

      // cache some commonly needed numbers
      const size_t data_set_size( features.GetNumberRows());
      const size_t feature_size( features.GetNumberCols());
      const size_t result_size( results.GetNumberCols());

      // initialize the cluster index vector, which stores the cluster id of each feature
      storage::Vector< size_t> cluster_indices( data_set_size, util::GetUndefined< size_t>());

      // initialize the clusters
      storage::Vector< math::RunningAverage< linal::Vector< float> > > clusters
      (
        m_NumberClusters,
        math::RunningAverage< linal::Vector< float> >( linal::Vector< float>( feature_size, 0.0))
      );

      // make objects to rescale the input/output, and rescaled data, if desired
      util::ShPtr< linal::Matrix< float> > sp_feature_matrix;

      // make objects to rescale the input/output, and rescaled data, if desired
      RescaleFeatureDataSet rescale_feature_result;

      // make an object to hold si-ptrs to all the feature-results
      storage::Vector< storage::VectorND< 2, linal::VectorConstReference< float> > > feature_results( data_set_size);

      // perform the scaling functions
      if( m_Autoscale)
      {
        sp_feature_matrix = util::CloneToShPtr( REDUCED_MATRIX->GetFeaturesPtr()->GetRawMatrix());

        linal::Matrix< float> &features_rescaled( *sp_feature_matrix);

        rescale_feature_result = RescaleFeatureDataSet( *sp_feature_matrix, math::Range< float>( 0, 1));
        rescale_feature_result.RescaleMatrix( features_rescaled);

        for( size_t index( 0); index < data_set_size; ++index)
        {
          feature_results( index) =
            storage::VectorND< 2, linal::VectorConstReference< float> >
            (
              linal::VectorConstReference< float>( feature_size, features_rescaled[ index]),
              linal::VectorConstReference< float>( result_size, results[ index])
            );
        }
      }
      else
      {
        // just add the feature-results to the si-ptr vector
        for( size_t index( 0); index < data_set_size; ++index)
        {
          feature_results( index) =
            storage::VectorND< 2, linal::VectorConstReference< float> >
            (
              linal::VectorConstReference< float>( feature_size, features[ index]),
              linal::VectorConstReference< float>( result_size, results[ index])
            );
        }
      }

      // initialize the clusters with the first m_NumberClusters features
      for( size_t i( 0); i < m_NumberClusters; ++i)
      {
        clusters( i) += feature_results( i).First();

        // set the cluster indices up to the associated cluster, since we know that these are the closest values
        // this isn't necessary, but it speeds up the search for the closest node later on
        cluster_indices( i) = i;
      }

      // refine clusters for m_MaxReclusterAttempts or until they do not change, whichever comes first
      size_t clustering_attempt_number( 0);
      while
        (
            clustering_attempt_number < m_MaxReclusterAttempts
            && RefineClusters( feature_results, cluster_indices, clusters, TrainingSchedule())
        )
      {
        ++clustering_attempt_number;
      }
      return cluster_indices;
    }

    //! @brief refines the current set of clusters
    //! @param FEATURES the features to cluster
    //! @param FEATURE_CLUSTER_INDICES the index of the closest mean in the clusters last round, for each feature
    //! @param CLUSTERS the current vector of clusters, which will be updated
    //! @return true if the clusters changed
    bool DataSetReducedToKMeans::RefineClusters
    (
      const storage::Vector< storage::VectorND< 2, linal::VectorConstReference< float> > > &FEATURES,
      storage::Vector< size_t> &FEATURE_CLUSTER_INDICES,
      storage::Vector< math::RunningAverage< linal::Vector< float> > > &CLUSTERS,
      const TrainingSchedule &TRAINING_SCHEDULE
    ) const
    {
      size_t nr_features( TRAINING_SCHEDULE.GetSize());
      bool cluster_changed( false);

      // find the closest centroid for each feature
      for( size_t itr( 0); itr < nr_features; ++itr)
      {
        // get a reference to the current feature
        const linal::VectorConstInterface< float> &feature( FEATURES( TRAINING_SCHEDULE( itr)).First());

        // maintain the shortest distance found between this feature and any centroid, as well as the id of the centroid
        float best_distance( std::numeric_limits< float>::infinity());
        size_t best_cluster_id( std::numeric_limits< size_t>::max());

        // usually only a few features change clusters, so first get the distance to the cluster that this
        // feature was last associated with.  this way, our best_distance is likely to be the shortest, in which
        // case the limited square distance calculation in the for loop below can stop early because it will be
        // the other centroids will mostly be much farther away than the closest centroid last time
        if( FEATURE_CLUSTER_INDICES( TRAINING_SCHEDULE( itr)) < m_NumberClusters)
        {
          // distance found
          best_distance =
            LimitedSquareDistance
            (
              feature,
              CLUSTERS( FEATURE_CLUSTER_INDICES( TRAINING_SCHEDULE( itr))).GetAverage(),
              best_distance
            );
          best_cluster_id = FEATURE_CLUSTER_INDICES( TRAINING_SCHEDULE( itr));
        }

        // look in the other clusters, to see if any of them are closer than the best-distance away
        for( size_t itr_cluster( 0); itr_cluster < m_NumberClusters; ++itr_cluster)
        {
          // ignore the cluster that was closest last round, since we already have its distance
          if( itr_cluster == FEATURE_CLUSTER_INDICES( TRAINING_SCHEDULE( itr)))
          {
            continue;
          }

          // compute the distance to the cluster with id itr_cluster, limited by the best_distance seen so far
          const float distance( LimitedSquareDistance( feature, CLUSTERS( itr_cluster).GetAverage(), best_distance));

          // if the new distance is better than the previous best, update the best distance
          if( distance < best_distance)
          {
            best_distance = distance;
            best_cluster_id = itr_cluster;
          }
        }

        // check whether the best cluster id changed from the last round
        if( best_cluster_id != FEATURE_CLUSTER_INDICES( TRAINING_SCHEDULE( itr)))
        {
          // yep, so update the cluster_changed bool and the closest cluster indice for this feature
          cluster_changed = true;
          FEATURE_CLUSTER_INDICES( TRAINING_SCHEDULE( itr)) = best_cluster_id;
        }
      }

      // reset all the clusters, so that we can find the new centroids
      for( size_t i( 0); i < m_NumberClusters; ++i)
      {
        CLUSTERS( i).Reset();
      }

      // add the features to the nearest cluster's running average to compute the new clusters
      for( size_t i( 0); i < nr_features; ++i)
      {
        CLUSTERS( FEATURE_CLUSTER_INDICES( TRAINING_SCHEDULE( i))) += FEATURES( TRAINING_SCHEDULE( i)).First();
      }

      return cluster_changed;
    }

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t DataSetReducedToKMeans::GetNominalSize() const
    {
      return std::min( m_Retriever->GetNominalSize(), m_NumberClusters);
    }

  } // namespace model
} // namespace bcl
