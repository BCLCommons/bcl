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
#include "model/bcl_model_kohonen_network_average.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "sched/bcl_sched_mutex.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_tertiary_function_job_with_data.h"
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> KohonenNetworkAverage::s_Instance
    (
      GetObjectInstances().AddInstance( new KohonenNetworkAverage())
    );

    //! @brief predict result with model using a NOT rescaled feature vector
    //! @param FEATURES rescaled features
    //! @return predicted result vector using a model
    FeatureDataSet< float> KohonenNetworkAverage::PredictWithoutRescaling
    (
      const FeatureDataSetInterface< float> &FEATURES
    ) const
    {
      storage::Vector< size_t> winning_nodes;
      storage::Vector< sched::Mutex> muteces;

      const size_t number_features( FEATURES.GetNumberFeatures());

      GetWinningNodeIndices
      (
        FEATURES,
        math::Range< size_t>( 0, FEATURES.GetNumberFeatures()),
        storage::Vector< size_t>(),
        winning_nodes,
        muteces
      );

      // create a results matrix
      linal::Matrix< float> results( FEATURES.GetNumberFeatures(), m_Codebook( 0).GetResultVector().GetSize());

      for( size_t feature_id( 0); feature_id < number_features; ++feature_id)
      {
        const linal::Vector< float> &result_vector( m_Codebook( winning_nodes( feature_id)).GetResultVector());
        std::copy( result_vector.Begin(), result_vector.End(), results[ feature_id]);
      }

      return FeatureDataSet< float>( results);
    }

    //! @brief predict result with model using a rescaled feature vector
    //! @param FEATURE normalized or rescaled feature vector
    //! @return predicted result vector using a model
    FeatureDataSet< float> KohonenNetworkAverage::operator()( const FeatureDataSetInterface< float> &FEATURE) const
    {
      // handle the case where rescaling is necessary
      if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_RescaleInput)
      {
        FeatureDataSet< float> feature( FEATURE);
        feature.Rescale( *m_RescaleInput);
        return PredictWithoutRescaling( feature).DeScale();
      }

      // data is already rescaled
      return PredictWithoutRescaling( FEATURE).DeScale();
    }

    //! @brief set the nodes up according to the map dimensions
    //! @param NORMALIZED_DATA data that can be used to set the initial positions of the nodes
    void KohonenNetworkAverage::InitializeNodes
    (
      const descriptor::Dataset &NORMALIZED_DATA
    )
    {
      // create a list with all the positions
      storage::List< linal::Vector< float> > positions;
      const size_t dimension( m_MapDimensions.GetSize());
      BCL_Assert( dimension, "Kohonen networks cannot be 0D!")

      // create an initial position vector
      linal::Vector< double> position( dimension, 0.0);

      while( position( 0) < m_MapDimensions( 0))
      {
        // add the current position to the positions list, converting it to float along the way
        positions.PushBack( linal::Vector< float>( position.Begin(), position.End()));

        // walk from the end of position back to the beginning to find the next point
        for
        (
          double *itr_position( position[ dimension - 1]),
                 *itr_position_end( position.Begin()),
                 *itr_map_dim( m_MapDimensions[ dimension - 1]);
          ; // condition in loop
          --itr_position, --itr_map_dim
        )
        {
          *itr_position += double( 1.0);
          if( *itr_position >= *itr_map_dim && itr_position != itr_position_end) // wrap around
          {
            *itr_position -= *itr_map_dim;
          }
          else
          {
            break;
          }
        }
      }

      // create the code book
      m_Codebook = storage::Vector< KohonenNode>( positions.GetSize());

      // iterate over the nodes and setup their positions and feature/results properly
      storage::List< linal::Vector< float> >::iterator itr_position( positions.Begin());

      // if no data was given, set up the network to have all 0 features and results
      if( NORMALIZED_DATA.IsEmpty())
      {
        // if there is a valid rescale function, use it to determine the correct dimensions
        const size_t feature_size
        (
          m_RescaleInput.IsDefined()
          ? m_RescaleInput->GetSize()
          : 0
        );
        const size_t result_size( 0);

        // create 0-initialized vectors for the average feature and result for each node
        linal::Vector< float> initial_feature( feature_size, float( 0.0)), initial_result( result_size, float( 0.0));

        for
        (
          storage::Vector< KohonenNode>::iterator itr( m_Codebook.Begin()), itr_end( m_Codebook.End());
          itr != itr_end;
          ++itr, ++itr_position
        )
        {
          *itr = KohonenNode( *itr_position, initial_feature, initial_result);
        }
      }
      else // data was given, use it to set up the initial positions
      {
        size_t data_number( 0);
        const size_t number_data( NORMALIZED_DATA.GetSize());
        // use the positions from the normalized data vector
        for
        (
          storage::Vector< KohonenNode>::iterator itr( m_Codebook.Begin()), itr_end( m_Codebook.End());
          itr != itr_end;
          ++itr, ++itr_position, ++data_number
        )
        {
          const size_t wrapped_data_number( data_number % number_data);
          *itr =
            KohonenNode
            (
              *itr_position,
              NORMALIZED_DATA.GetFeaturesPtr()->operator()( wrapped_data_number),
              NORMALIZED_DATA.GetResultsPtr()->operator()( wrapped_data_number)
            );
        }
      }
    }

    //! @brief reset sets the weight of all nodes to 0
    void KohonenNetworkAverage::Reset()
    {
      for
      (
        storage::Vector< KohonenNode>::iterator itr( m_Codebook.Begin()), itr_end( m_Codebook.End());
        itr != itr_end;
        ++itr
      )
      {
        // reset the weight of each node (positions, features, and result averages remain the same)
        itr->Reset();
      }
    }

    //! @brief gets the winning node
    //! @param VECTOR the sample of data from the training data
    //! @param HINT a hint for the index of the winning node; doubles speed on average if set to previous winner
    //! @return the winning node, paired with the error
    size_t KohonenNetworkAverage::GetIndexOfWinningNode
    (
      const linal::VectorConstInterface< float> &VECTOR,
      const size_t &HINT
    ) const
    {
      BCL_Assert( HINT < m_Codebook.GetSize(), "Hint node doesn't exist!");
      size_t winner( HINT);

      // Go through all code vectors
      double best_square_distance( linal::SquareDistance( m_Codebook( HINT).GetFeatureVector(), VECTOR));

      for( size_t counter( 0), number_codes( m_Codebook.GetSize()); counter < number_codes; ++counter)
      {
        if( HINT == counter)
        {
          // HINT was already used
          continue;
        }
        // get a reference to the current feature
        const linal::Vector< float> &feature( m_Codebook( counter).GetFeatureVector());
        BCL_Assert( feature.GetSize() == VECTOR.GetSize(), "Different vector sizes!");

        // compute the square distance between VECTOR and the current node's feature
        // because computing the distance is costly, only go far enough to see whether the current distance is longer
        // than the shortest distance
        double square_distance( 0.0);
        for
        (
          const float *itr_feature( feature.Begin()), *itr_feature_end( feature.End()), *itr_vector( VECTOR.Begin());
          square_distance < best_square_distance && itr_feature != itr_feature_end;
          ++itr_feature, ++itr_vector
        )
        {
          square_distance += math::Sqr( *itr_feature - *itr_vector);
        }

        // If distance is smaller than previous distances keep the new node
        if( square_distance < best_square_distance)
        {
          winner = counter;
          best_square_distance = square_distance;
        }
      }

      // return best node position
      return winner;
    }

    //! @brief adds each node from the other network to this network. The other network is probably a partial network.
    //! @param OTHER another Kohonen network to add from
    //! @return the resulting network
    KohonenNetworkAverage &KohonenNetworkAverage::operator +=( const KohonenNetworkAverage &OTHER)
    {
      storage::Vector< KohonenNode>::const_iterator from( OTHER.m_Codebook.Begin());
      for
      (
        storage::Vector< KohonenNode>::iterator to( m_Codebook.Begin()), to_end( m_Codebook.End());
        to != to_end;
        ++to, ++from
      )
      {
        *to += *from;
      }

      return *this;
    }

    //! @brief adds each node from the other network to this network. The other network is probably a partial network.
    //! @param OTHER another Kohonen network to add from
    //! @return the resulting network
    KohonenNetworkAverage &KohonenNetworkAverage::AddNetworkWithWeight
    (
      const KohonenNetworkAverage &OTHER,
      const float &WEIGHT
    )
    {
      storage::Vector< KohonenNode>::const_iterator from( OTHER.m_Codebook.Begin());
      for
      (
        storage::Vector< KohonenNode>::iterator to( m_Codebook.Begin()), to_end( m_Codebook.End());
        to != to_end;
        ++to, ++from
      )
      {
        to->AddNodeWithWeight( *from, WEIGHT);
      }

      return *this;
    }

    //! @brief give centroid values to nodes with no associated nearest features
    KohonenNetworkAverage &KohonenNetworkAverage::FixEmptyNodes()
    {
      for
      (
        storage::Vector< KohonenNode>::iterator to( m_Codebook.Begin()), to_end( m_Codebook.End());
        to != to_end;
        ++to
      )
      {
        if( to->GetWeight() > 0.0)
        {
          continue;
        }
        KohonenNode nearby;
        for
        (
          storage::Vector< KohonenNode>::const_iterator itr( m_Codebook.Begin()), itr_end( m_Codebook.End());
          itr != itr_end;
          ++itr
        )
        {
          if( itr->GetWeight())
          {
            nearby.MapData
            (
              itr->GetFeatureVector(),
              itr->GetResultVector(),
              1.0 / linal::SquareDistance( to->GetPosition(), itr->GetPosition())
            );
          }
        }
        to->MapData
        (
          nearby.GetFeatureVector(),
          nearby.GetResultVector(),
          0.0
        );
      }
      return *this;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &KohonenNetworkAverage::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_RescaleInput, ISTREAM);
      io::Serialize::Read( m_Codebook, ISTREAM);
      io::Serialize::Read( m_MapDimensions, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &KohonenNetworkAverage::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_RescaleInput, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Codebook, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MapDimensions, OSTREAM, INDENT);
      return OSTREAM;
    }

    //! @brief thread worker function for finding the winning indices for a data set with the network
    //! @param DATA the data to operate on
    //! @param RANGE the designated range of the training data to work on
    //! @param WINNING_NODES a vector where the winning node indices will be stored
    void KohonenNetworkAverage::GetWinningNodeIndicesHelper
    (
      const FeatureDataSetInterface< float> &DATA,
      const storage::Pair< math::Range< size_t>, util::SiPtr< const storage::Vector< size_t> > > &RANGE_AND_ORDER,
      storage::Pair
      <
        util::SiPtr< storage::Vector< size_t> >,
        util::SiPtr< storage::Vector< sched::Mutex> >
      > &WINNING_NODES_AND_MUTECES
    ) const
    {
      // get the actual range
      const math::Range< size_t> &range( RANGE_AND_ORDER.First());
      const storage::Vector< size_t> &order( *RANGE_AND_ORDER.Second());
      storage::Vector< size_t> &winning_indices( *WINNING_NODES_AND_MUTECES.First());
      storage::Vector< sched::Mutex> &muteces( *WINNING_NODES_AND_MUTECES.Second());

      if( !muteces.IsEmpty())
      {
        // threaded case
        // associate the data with the winning nodes in this range.
        for
        (
          size_t data_index( range.GetMin()), data_index_last( range.GetMax());
          data_index < data_index_last;
          ++data_index
        )
        {
          // get the real node
          const size_t node( order( data_index));
          // get the previous winning node
          muteces( node).Lock();
          const size_t prev_winning_node( winning_indices( node));
          muteces( node).Unlock();
          // find the newest winning node
          const size_t winning_node( GetIndexOfWinningNode( DATA( node), prev_winning_node));
          muteces( node).Lock();
          winning_indices( node) = winning_node;
          muteces( node).Unlock();
        }
      }
      else
      {
        // serial case
        for
        (
          size_t data_index( range.GetMin()), data_index_last( range.GetMax());
          data_index < data_index_last;
          ++data_index
        )
        {
          // get the real node
          const size_t node( order( data_index));
          // update the winning index
          winning_indices( node) = GetIndexOfWinningNode( DATA( node), winning_indices( node));
        }
      }

    }

    //! @brief determines the winning node indices of a normalized data set
    //! @param NORMALIZED_DATA all the training data, already normalized
    //! @param RANGE the range of training data to determine the winning indices for
    //! @param ORDER the order in which to visit the nodes
    //! @param PREVIOUS_WINNERS previous rounds node winners, will be updated
    //! @note this function is threaded for performance
    void KohonenNetworkAverage::GetWinningNodeIndices
    (
      const FeatureDataSetInterface< float> &NORMALIZED_DATA,
      const math::Range< size_t> &RANGE,
      const storage::Vector< size_t> &ORDER,
      storage::Vector< size_t> &PREVIOUS_WINNERS,
      storage::Vector< sched::Mutex> &MUTECES
    ) const
    {
      // construct the standardized range so that we don't have to worry about borders
      const math::Range< size_t> std_range( RANGE.StandardizeRange());

      storage::Vector< size_t> order;
      if( ORDER.IsEmpty())
      {
        order = storage::CreateIndexVector( NORMALIZED_DATA.GetNumberFeatures());
      }
      const storage::Vector< size_t> &order_ref( ORDER.IsEmpty() ? order : ORDER);

      const size_t max_data( std::min( std_range.GetMax(), order_ref.GetSize()));

      // construct a new range limited to the range of data available in normalized data
      math::Range< size_t> range
      (
        math::RangeBorders::e_LeftClosed, std_range.GetMin(), max_data, math::RangeBorders::e_RightOpen
      );

      // if there is no data, return an empty vector
      if( range.IsEmpty())
      {
        return;
      }

      // size of dataset
      const size_t data_set_size( range.GetWidth());

      // group for scheduler
      const size_t group_id( 0);

      const size_t number_of_available_cpus( std::min( sched::GetNumberCPUs(), data_set_size));

      // data interval for jobs
      const size_t interval( data_set_size / number_of_available_cpus);

      // calculate how many threads will operate on interval + 1 features
      const size_t number_threads_with_interval_plus_one( data_set_size % number_of_available_cpus);

      // scheduler for managing concurrent jobs
      util::ShPtrVector< sched::JobInterface> schedule;

      // number of jobs
      size_t number_jobs( number_of_available_cpus);

      schedule.AllocateMemory( number_jobs);

      storage::Vector< storage::Pair< math::Range< size_t>, util::SiPtr< const storage::Vector< size_t> > > > ranges;
      ranges.AllocateMemory( number_jobs);
      if( PREVIOUS_WINNERS.GetSize() != NORMALIZED_DATA.GetNumberFeatures())
      {
        PREVIOUS_WINNERS.Resize( NORMALIZED_DATA.GetNumberFeatures(), size_t( 0));
      }

      // create ranges
      for
      (
        size_t job_number( 0), end_position( range.GetMin());
        job_number < number_jobs;
        ++job_number
      )
      {
        // save the previous end position as the current start position
        const size_t start_position( end_position);

        // add the interval
        end_position += interval;

        // add one if necessary
        if( job_number < number_threads_with_interval_plus_one)
        {
          ++end_position;
        }

        // push back this range
        ranges.PushBack
        (
          storage::Pair< math::Range< size_t>, util::SiPtr< const storage::Vector< size_t> > >
          (
            math::Range< size_t>( start_position, end_position),
            util::SiPtr< const storage::Vector< size_t> >( order_ref)
          )
        );
      }

      storage::Pair
      <
        util::SiPtr< storage::Vector< size_t> >,
        util::SiPtr< storage::Vector< sched::Mutex> >
      > winners_and_muteces( PREVIOUS_WINNERS, MUTECES);

      // create jobs
      for
      (
        storage::Vector< storage::Pair< math::Range< size_t>, util::SiPtr< const storage::Vector< size_t> > > >::const_iterator
          itr_ranges( ranges.Begin()), itr_ranges_end( ranges.End());
        itr_ranges != itr_ranges_end;
        ++itr_ranges
      )
      {
        // create and pushback jobs into scheduler
        schedule.PushBack
        (
          util::ShPtr< sched::JobInterface>
          (
            new sched::TertiaryFunctionJobWithData
            <
              const FeatureDataSetInterface< float>,
              const storage::Pair< math::Range< size_t>, util::SiPtr< const storage::Vector< size_t> > >,
              storage::Pair
              <
                util::SiPtr< storage::Vector< size_t> >,
                util::SiPtr< storage::Vector< sched::Mutex> >
              >,
              void,
              KohonenNetworkAverage
            >
            (
              group_id,
              *this,
              &KohonenNetworkAverage::GetWinningNodeIndicesHelper,
              NORMALIZED_DATA,
              *itr_ranges,
              winners_and_muteces,
              sched::JobInterface::e_READY,
              NULL
            )
          )
        );

        // submit this particular job to start a thread
        sched::GetScheduler().RunJob( schedule.LastElement());
      }

      // join all jobs ( need to be finished)
      for( size_t job_id( 0); job_id < number_jobs; ++job_id)
      {
        sched::GetScheduler().Join( schedule( job_id));
      }
    } // Calibrate

    //! @brief Set the scaling of a feature set according to the model
    //! @param FEATURES feature set of interest
    //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
    //!       when operator() is called
    void KohonenNetworkAverage::Rescale( FeatureDataSet< float> &FEATURE) const
    {
      if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_RescaleInput)
      {
        FEATURE.DeScale();
        FEATURE.Rescale( *m_RescaleInput);
      }
    }

  } // namespace model
} // namespace bcl
