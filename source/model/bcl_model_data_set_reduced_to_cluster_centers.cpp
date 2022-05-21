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
#include "model/bcl_model_data_set_reduced_to_cluster_centers.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> DataSetReducedToClusterCenters::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new DataSetReducedToClusterCenters())
    );

    //! @brief get the name associated with a particular selection
    //! @param SELECTION the selection type of interest
    //! @return the name associated with that selection type
    const std::string &DataSetReducedToClusterCenters::GetSelectionName( const Selection &SELECTION)
    {
      static const std::string s_names[] =
      {
        "First",
        "Center",
        "Max",
        "Min",
        GetStaticClassName< Selection>()
      };
      return s_names[ SELECTION];
    }

    //! @brief constructor from the rmsd
    //! @param RMSD the rmsd to use in the clustering algorithm
    //! @param RANGE the min and max # of clusters that would be acceptable
    //! @param MAX_STEPS_TO_REACH_CLUSTER_SIZE the maximum # of attempts to change the rmsd to get the # of clusters into the range
    //! @param AUTOSCALE true if the data need to be normalized before clustering
    DataSetReducedToClusterCenters::DataSetReducedToClusterCenters
    (
      const double &RMSD,
      const math::Range< size_t> &RANGE,
      const size_t &MAX_STEPS_TO_REACH_CLUSTER_SIZE,
      const bool &AUTOSCALE,
      const Selection &SELECTION
    ) :
      m_RMSD( RMSD),
      m_ClusterSizeRange( RANGE),
      m_MaxSteps( MAX_STEPS_TO_REACH_CLUSTER_SIZE),
      m_Autoscale( AUTOSCALE),
      m_Selection( SELECTION)
    {
      BCL_Assert( m_RMSD >= 0.0 && m_RMSD <= 1.0, "RMSD must be between 0.0 and 1.0");
    }

    //! @brief Clone function
    //! @return pointer to new DataSetReducedToClusterCenters
    DataSetReducedToClusterCenters *DataSetReducedToClusterCenters::Clone() const
    {
      return new DataSetReducedToClusterCenters( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &DataSetReducedToClusterCenters::GetClassIdentifier() const
    {
      return GetStaticClassName< DataSetReducedToClusterCenters>();
    }

    //! @brief return the rmsd parameter
    //! @return the rmsd parameter
    double DataSetReducedToClusterCenters::GetRMSDParameter() const
    {
      return m_RMSD;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &DataSetReducedToClusterCenters::GetAlias() const
    {
      static const std::string s_Name( "Downsample");
      return s_Name;
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    //! @return the code / label for the feature (1st part) of the data set
    void DataSetReducedToClusterCenters::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectFeatures( CODE);
      m_Retriever->SelectFeatures( CODE);
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @return the code / label for the result (2nd part) of the data set
    void DataSetReducedToClusterCenters::SelectResults( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectResults( CODE);
      m_Retriever->SelectResults( CODE);
    }

    //! @brief Set which id columns to retrieve
    //! @param CODE the id column names to retrieve
    void DataSetReducedToClusterCenters::SelectIds( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectIds( CODE);
      m_Retriever->SelectIds( CODE);
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > DataSetReducedToClusterCenters::GetNumberPartitionsAndIds() const
    {
      return m_Retriever->GetNumberPartitionsAndIds();
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet DataSetReducedToClusterCenters::GetFeatureLabelsWithSizes() const
    {
      return m_Retriever->GetFeatureLabelsWithSizes();
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet DataSetReducedToClusterCenters::GetResultCodeWithSizes() const
    {
      return m_Retriever->GetResultCodeWithSizes();
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet DataSetReducedToClusterCenters::GetIdCodeWithSizes() const
    {
      return m_Retriever->GetIdCodeWithSizes();
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief generate dataset, reduced to the desired # of cluster centers
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset>
      DataSetReducedToClusterCenters::GenerateDataSet()
    {
      // get the base data set and use the clustering algorithm on it
      util::ShPtr< descriptor::Dataset> changeable( m_Retriever->GenerateDataSet());
      this->operator()( changeable);
      return changeable;
    }

    //! @brief reduce a data set
    //! @param DATA the data set
    //! @return the reduced data set
    void DataSetReducedToClusterCenters::operator ()( util::ShPtr< descriptor::Dataset> DATA) const
    {
      // make references to the matrices
      const linal::MatrixConstReference< float> features( DATA->GetFeaturesReference());

      // cache some commonly needed numbers
      const size_t data_set_size( DATA->GetSize());

      // make objects to rescale the input/output, and rescaled data, if desired
      util::ShPtr< linal::Matrix< float> > sp_feature_matrix;

      // simple pointer to the feature matrix that will be used; this will be changed to the autoscaled matrix
      // if that matrix is selected
      util::SiPtr< const linal::MatrixConstInterface< float> > si_ptr_feature_matrix( features);

      if( m_Autoscale)
      {
        // make objects to rescale the input/output, and rescaled data, if desired
        RescaleFeatureDataSet rescale_features( DATA->GetFeaturesReference(), math::Range< float>( 0.0, 1.0));

        // clone the features
        sp_feature_matrix = util::ShPtr< linal::Matrix< float> >( new linal::Matrix< float>( features));
        si_ptr_feature_matrix = sp_feature_matrix;
        rescale_features.RescaleMatrix( *sp_feature_matrix);
      }

      // determine the actual min and max # of clusters, which can be no larger than the number of feature results
      const size_t min_number_clusters( std::min( data_set_size, m_ClusterSizeRange.GetMin()));
      const size_t max_number_clusters( std::min( data_set_size, m_ClusterSizeRange.GetMax()));

      // make an object to hold the representative feature indices
      storage::Vector< size_t> representatives;

      double min_rmsd( 0.0), max_rmsd( 2.0 * m_RMSD);

      // keep track of how many steps have been performed so far
      size_t number_steps( 0);

      // track cluster ids
      storage::Vector< size_t> cluster_ids;

      // determine the max rmsd by doubling max_rmsd until the number of clusters is below the maximum number
      while( number_steps < m_MaxSteps)
      {
        representatives = GetRepresentativeFeatures( *si_ptr_feature_matrix, max_rmsd, cluster_ids);
        if( representatives.GetSize() > max_number_clusters)
        {
          min_rmsd = max_rmsd;
          max_rmsd *= 2.0;
        }
        else
        {
          // the max rmsd is large enough, so stop
          break;
        }
        ++number_steps;
      }

      // now perform a binary search for the right rmsd
      for( ; number_steps < m_MaxSteps; ++number_steps)
      {
        const double rmsd_guess( ( min_rmsd + max_rmsd) / 2.0);
        representatives = GetRepresentativeFeatures( *si_ptr_feature_matrix, rmsd_guess, cluster_ids);

        if( representatives.GetSize() > max_number_clusters)
        {
          min_rmsd = rmsd_guess * ( 1.0 + std::numeric_limits< double>::epsilon());
        }
        else if( representatives.GetSize() < min_number_clusters)
        {
          max_rmsd = rmsd_guess / ( 1.0 + std::numeric_limits< double>::epsilon());
        }
        else
        {
          break;
        }
      }

      if( m_Selection == e_Center)
      {
        representatives = GetClusterCenters( cluster_ids, *si_ptr_feature_matrix);
      }
      else if( m_Selection == e_Max)
      {
        storage::Vector< linal::VectorConstReference< float> > cluster_maxes( representatives.GetSize());
        const linal::MatrixReference< float> ref_results( DATA->GetResultsReference());
        for( size_t i( 0), n_rows( ref_results.GetNumberRows()); i < n_rows; ++i)
        {
          // account for id offset of 1
          const size_t effective_id( cluster_ids( i) - 1);
          if
          (
            cluster_maxes( effective_id).GetSize() == size_t( 0)
            || cluster_maxes( effective_id) < ref_results.GetRow( i)
          )
          {
            cluster_maxes( effective_id) = ref_results.GetRow( i);
            representatives( effective_id) = i;
          }
        }
      }
      else if( m_Selection == e_Min)
      {
        storage::Vector< linal::VectorConstReference< float> > cluster_mins( representatives.GetSize());
        const linal::MatrixReference< float> ref_results( DATA->GetResultsReference());
        for( size_t i( 0), n_rows( ref_results.GetNumberRows()); i < n_rows; ++i)
        {
          const size_t effective_id( cluster_ids( i) - 1);
          if
          (
            cluster_mins( effective_id).GetSize() == size_t( 0)
            || ref_results.GetRow( i) < cluster_mins( effective_id)
          )
          {
            cluster_mins( effective_id) = ref_results.GetRow( i);
            representatives( effective_id) = i;
          }
        }
      }

      // select the desired features
      BCL_MessageStd( "Representatives: " + util::Format()( representatives));
      DATA->KeepRows( representatives);
    }

    //! @brief reduce a data set
    //! @param DATA the data set
    //! @return the reduced data set
    util::ShPtr< descriptor::Dataset> DataSetReducedToClusterCenters::operator()( const descriptor::Dataset &DATA) const
    {
      util::ShPtr< descriptor::Dataset> changeable( DATA.HardCopy());
      this->operator()( changeable);
      return changeable;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get a set of features that cover the same space as FEATURES at resolution RMSD
    //! @param FEATURES the set of features to consider
    //! @param RMSD the rmsd; vectors closer than this rmsd will be considered representative of each other
    //! @param CLUSTER_IDS vector that will store cluster ids (1-offset), used for refinement if m_ClusterIDs was set
    //! @return indices of features that cover the same space as FEATURES at resolution RMSD
    storage::Vector< size_t> DataSetReducedToClusterCenters::GetRepresentativeFeatures
    (
      const linal::MatrixConstInterface< float> &FEATURES,
      const double &RMSD,
      storage::Vector< size_t> &CLUSTER_IDS
    ) const
    {
      const size_t data_set_size( FEATURES.GetNumberRows());
      const size_t feature_size( FEATURES.GetNumberCols());

      // given the rmsd cutoff, we can get a minor speed improvement by computing the maximum square deviation
      // between the data set feature vectors, because then there is no need to take a square root or divide by the number
      // of features in the end
      const double mean_square_deviation( RMSD * RMSD);
      const double max_square_deviation( mean_square_deviation * feature_size);

      // keep track of how many points have not yet been clustered
      size_t unclustered_points( data_set_size);

      BCL_MessageDbg( "rmsd: " + util::Format()( max_square_deviation));

      // make a vector that holds 0 if a datum is already used in a cluster
      CLUSTER_IDS = storage::Vector< size_t>( data_set_size, size_t( 0));

      // make an object to hold the clusters that are found
      storage::Vector< size_t> cluster_centers;

      // allocate enough memory so that the vector never needs reallocation
      cluster_centers.AllocateMemory( data_set_size);
      for( size_t index( 0); index < data_set_size; ++index)
      {
        // ignore datums that are already in the cluster
        if( CLUSTER_IDS( index))
        {
          continue;
        }

        // add this data point to the set of clusters
        cluster_centers.PushBack( index);

        // if there are already too many clusters, stop
        if( cluster_centers.GetSize() > m_ClusterSizeRange.GetMax())
        {
          break;
        }
        CLUSTER_IDS( index) = cluster_centers.GetSize();
        --unclustered_points;

        // get a reference on the feature vector of the current data set vector, which defines the center of the cluster
        const float *cluster_vector_begin( FEATURES[ index]);

        // get the end of the current data set vector
        const float *cluster_vector_end( cluster_vector_begin + feature_size);

        // look at all other data that have not yet been placed in clusters, set their datum in cluster to 1 if
        // they are within the rmsd of this datum
        for( size_t index_b( index + 1); index_b < data_set_size; ++index_b)
        {
          // ignore datums that are already in the cluster
          if( CLUSTER_IDS( index_b))
          {
            continue;
          }

          // initialize a sum of the square deviation of the feature at index_b from the cluster vector
          double square_deviation( 0.0);

          // compute the square deviation, but stop if we get over the cutoff
          // since there will typically be many clusters, we usually will exceed the cutoff very quickly, so
          // short-circuiting is important here
          for
          (
            const float *itr_cluster_vec( cluster_vector_begin), *itr_current_vec( FEATURES[ index_b]);
            square_deviation < max_square_deviation && itr_cluster_vec != cluster_vector_end;
            ++itr_cluster_vec, ++itr_current_vec
          )
          {
            square_deviation += math::Sqr( *itr_cluster_vec - *itr_current_vec);
          }

          BCL_MessageDbg( " square_deviation between " + util::Format()( index) + " " + util::Format()( index_b) + " = " + util::Format()( square_deviation));

          // mark the datum at index_b as being already in a cluster so that we don't use it to form a new cluster
          if( square_deviation < max_square_deviation)
          {
            CLUSTER_IDS( index_b) = cluster_centers.GetSize();
            --unclustered_points;
          }
        }
        BCL_MessageVrb
        (
          "Point: " + util::Format()( index)
          + " # cluster centers: " + util::Format()( cluster_centers.GetSize())
          + " # unclustered points: " + util::Format()( unclustered_points)
        );

        // if there are no longer enough clusters to satisfy the min, then stop
        if( cluster_centers.GetSize() + unclustered_points < m_ClusterSizeRange.GetMin())
        {
          break;
        }
      }

      BCL_MessageStd( " RMSD: " + util::Format()( RMSD) + " #clusters: " + util::Format()( cluster_centers.GetSize()));

      return cluster_centers;
    }

    //! @brief refine the id set to be the center-most element of each cluster
    //! @param CLUSTER_IDS vector that will store cluster ids (1-offset), used for refinement if m_ClusterIDs was set
    //! @param FEATURES set of features to consider
    //! @return the updated feature ids
    storage::Vector< size_t> DataSetReducedToClusterCenters::GetClusterCenters
    (
      storage::Vector< size_t> &CLUSTER_IDS,
      const linal::MatrixConstInterface< float> &FEATURES
    ) const
    {
      storage::Vector< size_t> cluster_centers;
      const size_t data_set_size( FEATURES.GetNumberRows());
      const size_t feature_size( FEATURES.GetNumberCols());

      for( size_t current_cluster_id( 1);; ++current_cluster_id)
      {
        // get all the rows in this cluster
        storage::Vector< size_t> current_cluster_ids;
        storage::Vector< linal::VectorConstReference< float> > cluster_members;
        for( size_t row( 0); row < data_set_size; ++row)
        {
          if( CLUSTER_IDS( row) == current_cluster_id)
          {
            current_cluster_ids.PushBack( row);
            cluster_members.PushBack( linal::VectorConstReference< float>( feature_size, FEATURES[ row]));
          }
        }

        const size_t number_cluster_members( current_cluster_ids.GetSize());
        if( number_cluster_members == 0)
        {
          break;
        }
        // track sums of msd's for each row
        linal::Vector< float> msd_sum( number_cluster_members, float( 0.0));
        size_t min_rmsd_id( 0);
        for( size_t row_a( 0); row_a < number_cluster_members; ++row_a)
        {
          for( size_t row_b( row_a + 1); row_b < number_cluster_members; ++row_b)
          {
            // initialize a sum of the square deviation of the feature at index_b from the cluster vector
            const double rmsd( linal::Distance( cluster_members( row_a), cluster_members( row_b)));

            // add the squared deviation to the msd sum for each cluster point
            msd_sum( row_a) += rmsd;
            msd_sum( row_b) += rmsd;
          }

          // track the lowest-rmsd instance
          if( msd_sum( row_a) < msd_sum( min_rmsd_id))
          {
            min_rmsd_id = row_a;
          }
        }

        // get the actual feature id
        cluster_centers.PushBack( current_cluster_ids( min_rmsd_id));
      }

      return cluster_centers;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer DataSetReducedToClusterCenters::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription
      (
        "Downsample the dataset while maintaining near maximal coverage of feature space"
      );
      member_data.AddInitializer
      (
        "dataset",
        "dataset retriever to call to get the entire data set",
        io::Serialization::GetAgent( &m_Retriever)
      );

      member_data.AddInitializer
      (
        "rmsd",
        "features less than this root mean squared distance apart are initially eligible to be part of the same cluster",
        io::Serialization::GetAgentWithRange( &m_RMSD, 0.0, 1.0)
      );

      member_data.AddInitializer
      (
        "acceptable cluster size range",
        "range of # of desired points in the resulting data set.\n"
        "RMSD will be adjusted to try to get a number of representative features in this range",
        io::Serialization::GetAgent( &m_ClusterSizeRange)
      );

      member_data.AddInitializer
      (
        "max steps",
        "maximum number of times to change rmsd to try to downsample the data set to within the acceptable size of the cluster range",
        io::Serialization::GetAgent( &m_MaxSteps),
        "1"
      );

      member_data.AddInitializer
      (
        "autoscale",
        "whether to normalize / denormalize the data set before / after clustering",
        io::Serialization::GetAgent( &m_Autoscale),
        "0"
      );

      member_data.AddInitializer
      (
        "select",
        "Use First for the fastest selection of cluster-center-member (the first found in the dataset). "
        "Use Center to the center-most feature of each cluster (making the algorithm order independent and generally more precise)"
        "note that this will result in much slower performance if the number of clusters is small"
        "Use Max (Min) to choose the feature with the Max (or Min) result.  If multiple results, the min/max of each column is selected",
        io::Serialization::GetAgent( &m_Selection),
        "First"
      );

      return member_data;
    }

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t DataSetReducedToClusterCenters::GetNominalSize() const
    {
      return std::min( m_Retriever->GetNominalSize(), m_ClusterSizeRange.GetMax());
    }

  } // namespace model
} // namespace bcl
