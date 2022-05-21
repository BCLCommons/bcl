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
#include "model/bcl_model_approximator_kohonen_network.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "model/bcl_model_data_set_reduced_to_k_means.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_tertiary_function_job_with_data.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ApproximatorKohonenNetwork::s_Instance
    (
      util::Enumerated< ApproximatorBase>::AddInstance( new ApproximatorKohonenNetwork())
    );

    //! @brief Kernel as string
    //! @param KERNEL the kernel
    //! @return the string for the kernel
    const std::string &ApproximatorKohonenNetwork::GetKernelName( const ApproximatorKohonenNetwork::NeighborKernel &KERNEL)
    {
      static const std::string s_Names[] =
      {
        "Gaussian",
        "Bubble",
        GetStaticClassName< ApproximatorKohonenNetwork::NeighborKernel>()
      };
      return s_Names[ size_t( KERNEL)];
    }

    //! @brief Initializer as string
    //! @param INITIALIZER the initialization
    //! @return the string for the initializer
    const std::string &ApproximatorKohonenNetwork::GetInitializationName( const Initialization &INITIALIZER)
    {
      static const std::string s_Names[] =
      {
        "FirstVectors",           // Initial vectors taken from the start of the training data
        "RandomlyChosenVectors",  // Initial vectors randomly chosen from the training data
        "RandomlyChosenElements", // Initial vectors composed of randomly chosen values for each element of the vector
        "KMeans",
        GetStaticClassName< ApproximatorKohonenNetwork::Initialization>()
      };
      return s_Names[ size_t( INITIALIZER)];
    }

    //! @brief default constructor
    ApproximatorKohonenNetwork::ApproximatorKohonenNetwork() :
      m_Network(),
      m_Length( 1),
      m_Radius( 0.0),
      m_CurrentRadius( 0.0),
      m_UpdateEveryNthFeature( 0),
      m_NeighborKernel( e_Gaussian),
      m_Initializer( e_FirstVectors),
      m_NoTrack( false),
      m_RescaleType( RescaleFeatureDataSet::e_MinMax),
      m_Cutoff( 0.5)
    {
    }

    //! @brief constructor from all necessary parameters
    //! @param MAP_DIMENSIONS dimensions of the map
    //! @param INITIAL_LENGTH how many iterations to train for (if applicable).
    //! @param INITAL_RADIUS the initial neighborhood radius
    //! @param OBJECTIVE_FUNCTION the objective function from the approximator framework
    //! @param UPDATE_EVERY_NTH_FEATURE update the nodes after seeing this many features
    //! @param NEIGHBOR_KERNEL the neighbor kernel type
    ApproximatorKohonenNetwork::ApproximatorKohonenNetwork
    (
      const linal::Vector< double> &MAP_DIMENSIONS,
      const size_t &INITIAL_LENGTH,
      const float &INITAL_RADIUS,
      const util::ShPtr< ObjectiveFunctionWrapper> &OBJECTIVE_FUNCTION,
      const size_t &UPDATE_EVERY_NTH_FEATURE,
      const NeighborKernel &NEIGHBOR_KERNEL
    ) :
      ApproximatorBase( OBJECTIVE_FUNCTION),
      m_Network( MAP_DIMENSIONS),
      m_Length( INITIAL_LENGTH),
      m_Radius( INITAL_RADIUS),
      m_CurrentRadius( INITAL_RADIUS),
      m_UpdateEveryNthFeature( UPDATE_EVERY_NTH_FEATURE),
      m_NeighborKernel( NEIGHBOR_KERNEL),
      m_Initializer( e_FirstVectors),
      m_NoTrack( false),
      m_RescaleType( RescaleFeatureDataSet::e_MinMax),
      m_Cutoff( 0.5)
    {
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ApproximatorKohonenNetwork::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ApproximatorKohonenNetwork::GetAlias() const
    {
      static const std::string s_alias( "Kohonen");
      return s_alias;
    }

    //! @brief set training data set for a specific iterate in approximator framework
    //! @param DATA training data set
    void ApproximatorKohonenNetwork::SetTrainingData
    (
      util::ShPtr< descriptor::Dataset> &DATA
    )
    {
      m_CurrentRadius = m_Radius;
      m_CurrentIteration = 0;
      BCL_Assert( !DATA->IsEmpty(), "no training data loaded");

      const size_t dataset_size( DATA->GetResultsPtr()->GetNumberFeatures());
      m_TrainingData = DATA;
      m_TrainingData->GetFeatures().Rescale( math::Range< float>( 0.0, 1.0), m_RescaleType);

      BCL_MessageStd( "Setting up training data with " + util::Format()( DATA->GetSize()) + " points");

      // Setup the initial network's positions
      m_Network = KohonenNetworkAverage( m_Network.GetMapDimensions(), GetRescaleFeatureDataSet());

      m_Schedule.Setup( DATA->GetResults(), m_Cutoff);

      // choose the initial feature vectors
      m_Network.InitializeNodes( *m_TrainingData);
      if( m_Initializer == e_FirstVectors)
      {
        // map is already correctly initialized
      }
      else if( m_Initializer == e_RandomlyChosenElements || m_Initializer == e_RandomlyChosenVectors)
      {
        // create a new feature dataset with randomly selected vectors or elements
        // determine the # of nodes in the code book
        const size_t map_size( m_Network.GetCodeBook().GetSize());

        // Because the results are not used for node assignment, their initial values do not matter, so long
        // as they are the right length, thus, the initial results can always be taken from the first N vectors, where
        // N is the map size
        linal::Matrix< float> first_results( m_TrainingData->GetResultsPtr()->GetMatrix( 0, map_size));

        // determine the size of each feature
        const size_t feature_size( m_TrainingData->GetFeatureSize());

        // create a matrix to store the selected features
        linal::Matrix< float> selected_features( map_size, feature_size);

        if( m_Initializer == e_RandomlyChosenVectors)
        {
          // copy random vectors from the training data into each row of the selected features
          for( size_t node_id( 0); node_id < map_size; ++node_id)
          {
            // choose the training data point index to use
            const size_t chosen_feature_index( random::GetGlobalRandom().Random( dataset_size - 1));

            // get a feature reference to that feature
            const FeatureReference< float> chosen_feature
            (
              m_TrainingData->GetFeaturesPtr()->operator()( chosen_feature_index)
            );

            // copy that feature into the selected features matrix
            std::copy( chosen_feature.Begin(), chosen_feature.End(), selected_features[ node_id]);
          }
        }
        else // if( m_Initializer == e_RandomlyChosenElements)
        {
          // copy random elements from the training data into each row of the selected features
          for( size_t node_id( 0); node_id < map_size; ++node_id)
          {
            for( size_t element_id( 0); element_id < feature_size; ++element_id)
            {
              // choose the training data point index to use
              const size_t chosen_feature_index( random::GetGlobalRandom().Random( dataset_size - 1));

              // get a feature reference to that feature
              const FeatureReference< float> chosen_feature
              (
                m_TrainingData->GetFeaturesPtr()->operator()( chosen_feature_index)
              );

              // copy that feature into the element of the selected feature matrix
              selected_features( node_id, element_id) = chosen_feature( element_id);
            }
          }
        }

        // create a new feature data set with those matrices
        const descriptor::Dataset seed_positions( selected_features, first_results);

        // use that data set to set up the initial features
        m_Network.InitializeNodes( seed_positions);
      }
      else if( m_Initializer == e_KMeans)
      {
        // create a new feature dataset with randomly selected vectors or elements
        // determine the # of nodes in the code book
        const size_t map_size( m_Network.GetCodeBook().GetSize());
        DataSetReducedToKMeans kmeans( map_size, 4 * map_size, false);
        m_Network.InitializeNodes( *kmeans( DATA, m_Schedule));
      }

      size_t features_per_update( m_UpdateEveryNthFeature);
      // if the m_UpdateEveryNthFeature was 0, then just use the data set size
      if( m_UpdateEveryNthFeature == size_t( 0))
      {
        features_per_update = m_Schedule.GetSize();
      }

      // set up the training ranges for each thread
      m_DataSetRanges.Reset();

      // calculate the # of weigh updates per run through the data set
      const size_t number_epochs( m_Schedule.GetSize() / features_per_update);

      m_DataSetRanges.AllocateMemory( number_epochs);

      // if m_TrainingData->GetSize() is not a multiple of features_per_update, distribute the extra features
      // over the initial number_epochs_with_extra_feature epochs
      const size_t number_epochs_with_extra_feature( m_Schedule.GetSize() % features_per_update);

      for( size_t epoch_number( 0), end_index( 0); epoch_number < number_epochs; ++epoch_number)
      {
        const size_t start_index( end_index);

        // add the # of features that will be examined in this epoch
        end_index += features_per_update + size_t( epoch_number < number_epochs_with_extra_feature);

        m_DataSetRanges.PushBack
        (
          math::Range< size_t>
          (
            math::RangeBorders::e_LeftClosed,
            start_index,
            end_index,
            math::RangeBorders::e_RightOpen
          )
        );
      }

      SetupThreadRanges();

      m_LastClosestNodes.Resize( dataset_size);
      m_LastClosestNodes.SetAllElements( 0);
    }

    //! @brief construct a model from the current iterate
    //! @return shptr to the new model interface
    util::ShPtr< Interface> ApproximatorKohonenNetwork::GetCurrentModel() const
    {
      return util::ShPtr< Interface>( m_Network.Clone());
    }

    //! @brief returns the current approximation
    //! @return current argument result pair
    const util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> > ApproximatorKohonenNetwork::GetCurrentApproximation() const
    {
      util::ShPtr< Interface> new_network( m_Network.Clone());
      util::ShPtr< storage::Pair< util::ShPtr< Interface>, float> > return_value
      (
        new storage::Pair< util::ShPtr< Interface>, float>( new_network, m_ObjectiveFunction->operator ()( new_network))
      );
      return return_value;
    }

    //! @brief conducts the next approximation step and stores the approximation
    void ApproximatorKohonenNetwork::Next()
    {
      m_Schedule.Next();
      // assert when there is no interval
      BCL_Assert( m_TrainingData.IsDefined() && !m_TrainingData->IsEmpty(), "no training data given!");

      // stay at the last point if length is reached
      if( m_CurrentRadius < 0.4)
      {
        m_CurrentIteration = m_Length - 1;
      }
      else
      {
        m_CurrentRadius = std::max( m_CurrentRadius * ( 1.0 - m_CurrentIteration / ( 4.0 * m_Length)), 0.4);
      }

      util::ShPtrVector< sched::JobInterface> all_jobs( m_NumberThreadsNodes);

      const size_t group_id( 0);

      // store the id of each thread for use by the jobs
      storage::Vector< size_t> thread_ids( m_NumberThreadsNodes);
      for( size_t thread_number( 0); thread_number < m_NumberThreadsNodes; ++thread_number)
      {
        thread_ids( thread_number) = thread_number;
      }

      // create a const reference to the features of the dataset
      const FeatureDataSetInterface< float> &features( *m_TrainingData->GetFeaturesPtr());

      // if the training schedule includes balancing, it is necessary to have a mutex to protect from
      // concurrent access to m_PreviousWinners
      if( m_Schedule.IsBalanced())
      {
        m_Muteces.Resize( features.GetNumberFeatures());
      }

      // run through all the data once, updating nodes after every m_UpdateEveryNthFeature features
      for
      (
        size_t epoch_number( 0), number_epochs( m_DataSetRanges.GetSize());
        epoch_number < number_epochs;
        ++epoch_number
      )
      {
        // get the winning indices for this segment of the data set (threaded)
        m_Network.GetWinningNodeIndices
        (
          features,
          m_DataSetRanges( epoch_number),
          m_Schedule.GetOrder(),
          m_LastClosestNodes,
          m_Muteces
        );

        // run through the data set once, updating nodes after every epoch (m_UpdateEveryNthFeature features, data set size by default)
        for( size_t thread_number( 0); thread_number < m_NumberThreadsNodes; ++thread_number)
        {
          all_jobs( thread_number)
            = util::ShPtr< sched::JobInterface>
              (
                new sched::TertiaryFunctionJobWithData
                <
                  const size_t,
                  const size_t,
                  const size_t,
                  void,
                  ApproximatorKohonenNetwork
                >
                (
                  group_id,
                  *this,
                  &ApproximatorKohonenNetwork::TrainThread,
                  thread_ids( thread_number),
                  m_DataSetRanges( epoch_number).GetMin(),
                  m_DataSetRanges( epoch_number).GetMax(),
                  sched::JobInterface::e_READY,
                  NULL
                )
              );
          sched::GetScheduler().RunJob( all_jobs( thread_number));
        }

        // join the first job
        sched::GetScheduler().Join( all_jobs( 0));

        // join remaining jobs
        for( size_t job_id( 1); job_id < m_NumberThreadsNodes; ++job_id)
        {
          sched::GetScheduler().Join( all_jobs( job_id));

          // add their results to the first new network
          m_NewNetworks( 0) += m_NewNetworks( job_id);
        }

        m_NewNetworks.FirstElement().FixEmptyNodes();
        m_Network = m_NewNetworks( 0);
        m_Network.FixEmptyNodes();
      }

      m_CurrentIteration++;

      if( !m_NoTrack)
      {
        this->GetTracker().Track( GetCurrentApproximation());
      }
    }

    //! @brief thread helper function for iterate, works on a subset of the training data, computes partial
    //!        adjustments to the network, then adds them to the final result.
    //! @param RANGE_ID the partition of the training data to work on
    //! @param WINNING_INDICES indices of the winning nodes for the data being operated on
    //! @param START_DATA_NUMBER the index in the training data for the current epoch
    void ApproximatorKohonenNetwork::TrainThread
    (
      const size_t &RANGE_ID,
      const size_t &START_DATA_NUMBER,
      const size_t &END_DATA_NUMBER
    )
    {
      // get the network where results will be placed
      KohonenNetworkAverage &network( m_NewNetworks( RANGE_ID));

      // reset the weights in that network
      network.Reset();

      // get the range of nodes that will be iterated over
      const math::Range< size_t> &range( m_NodeRanges( RANGE_ID));

      // node count will be a histogram of nodes inside the range
      storage::Vector< size_t> node_counts( range.GetWidth(), size_t( 0));

      size_t number_winners( 0), last_winner( 0);

      // walk over all the winning indices
      for
      (
        storage::Vector< size_t>::const_iterator
          itr_winners( m_Schedule.GetOrder().Begin() + START_DATA_NUMBER),
          itr_winners_end( m_Schedule.GetOrder().Begin() + END_DATA_NUMBER);
        itr_winners != itr_winners_end;
        ++itr_winners
      )
      {
        const size_t winning_node( m_LastClosestNodes( *itr_winners));
        // check if this thread is responsible for this range
        if( range.IsWithin( winning_node))
        {
          // get a reference on the appropriate member of node count
          size_t &node_count( node_counts( winning_node - range.GetMin()));

          // if this was not previously a winning node, mark it as such
          if( node_count == 0)
          {
            node_count = 1;
            ++number_winners;
            last_winner = winning_node;
          }

          // map the data, no need for learning rate with batch iterate
          network.GetCodeBook()( winning_node).MapData
          (
            m_TrainingData->GetFeaturesPtr()->operator()( *itr_winners),
            m_TrainingData->GetResultsPtr()->operator()( *itr_winners)
          );
        }
      }

      if( number_winners == 0)
      {
        // no winners in this range, nothing to do
      }
      // optimization for when there is only one node that was chosen as a winner
      if( number_winners == 1)
      {
        // overlap is impossible, so adapt network directly

        // use a partial vector adapt function and the neighbor adapt function on the network
        AdaptNeighbors( network, network.GetCodeBook()( last_winner), m_CurrentRadius);
      }
      else
      {
        // need a swap network to avoid changing the results of this network
        KohonenNetworkAverage &network_swap( m_NewNetworksSwap( RANGE_ID));
        network_swap.Reset();

        for( size_t counter( 0), size( node_counts.GetSize()); counter < size; ++counter)
        {
          if( node_counts( counter) > 0)
          {
            // get the index of the node to update
            const size_t node_index( counter + range.GetMin());
            // use a partial vector adapt function and the neighbor adapt function on the network
            AdaptNeighbors( network_swap, network.GetCodeBook()( node_index), m_CurrentRadius);
          }
        }

        // swap the networks' codebooks
        network.GetCodeBook().InternalData().swap( network_swap.GetCodeBook().InternalData());
      }
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ApproximatorKohonenNetwork::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_Network, ISTREAM);
      io::Serialize::Read( m_Length, ISTREAM);
      io::Serialize::Read( m_CurrentIteration, ISTREAM);
      io::Serialize::Read( m_Radius, ISTREAM);
      io::Serialize::Read( m_TrainingData, ISTREAM);
      io::Serialize::Read( m_NeighborKernel, ISTREAM);
      io::Serialize::Read( m_Initializer, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ApproximatorKohonenNetwork::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_Network, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Length, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CurrentIteration, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Radius, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TrainingData, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NeighborKernel, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Initializer, OSTREAM, INDENT);
      return OSTREAM;
    }

    //! @brief setup the data ranges and node ranges for training
    //! This should only be called while setting the training data or after reading
    void ApproximatorKohonenNetwork::SetupThreadRanges()
    {
      // set up the training ranges for each thread
      m_NodeRanges.Reset();

      const size_t number_nodes( m_Network.GetCodeBook().GetSize());

      // determine the # of threads for splitting up nodes
      m_NumberThreadsNodes = std::min( number_nodes, sched::GetNumberCPUs());

      BCL_MessageStd( "Set up node ranges with # threads: " + util::Format()( m_NumberThreadsNodes));

      const size_t nodes_per_thread( number_nodes / m_NumberThreadsNodes);
      const size_t number_threads_with_additional_node( number_nodes % m_NumberThreadsNodes);
      m_NodeRanges.AllocateMemory( m_NumberThreadsNodes);

      for( size_t thread_number( 0), end_node( 0); thread_number < m_NumberThreadsNodes; ++thread_number)
      {
        // save the old end node #
        const size_t old_end_node( end_node);

        // determine the # of nodes that will be handled by this thread
        const size_t nodes_this_thread( nodes_per_thread + size_t( thread_number < number_threads_with_additional_node));

        end_node += nodes_this_thread;

        m_NodeRanges.PushBack
        (
          math::Range< size_t>
          (
            math::RangeBorders::e_LeftClosed, old_end_node, end_node, math::RangeBorders::e_RightOpen
          )
        );
      }

      m_NewNetworks = m_NewNetworksSwap = storage::Vector< KohonenNetworkAverage>( m_NumberThreadsNodes, m_Network);
    }

    //! @brief spreads the adaptation from a particular node to other nodes within RADIUS
    //! @param NETWORK the network to adapt, usually a copy to avoid conflicting changes
    //! @param NODE_TO_SPREAD is the node whose feature and result vectors should be propagated
    //! @param RADIUS the current neighborhood radius
    void ApproximatorKohonenNetwork::AdaptNeighbors
    (
      KohonenNetworkAverage &NETWORK,
      const KohonenNode &NODE_TO_SPREAD,
      const float &RADIUS
    ) const
    {
      const linal::Vector< float> &position( NODE_TO_SPREAD.GetPosition());

      const double node_weight( NODE_TO_SPREAD.GetWeight());

      const double max_radius( m_NeighborKernel == e_Gaussian ? 2.0 * RADIUS : RADIUS);
      for
      (
        storage::Vector< KohonenNode>::iterator
          itr( NETWORK.GetCodeBook().Begin()), itr_end( NETWORK.GetCodeBook().End());
        itr != itr_end;
        ++itr
      )
      {
        // the distance between two map nodes
        const float distance( linal::Distance( position, itr->GetPosition()));

        // only if learning will occur do we adapt the neighbor
        if( distance <= max_radius)
        {
          // the learning rate to be applied
          float learning_rate( node_weight);
          if( m_NeighborKernel == e_Gaussian)
          {
            learning_rate *= std::exp( float( -0.5 * math::Sqr( distance / RADIUS)));
          }
          // if m_function was bubble, then learning rate should == node_weight, which it already is

          if( position == itr->GetPosition())
          {
            itr->MapData( NODE_TO_SPREAD.GetFeatureVector(), NODE_TO_SPREAD.GetResultVector(), 0.8 * learning_rate);
          }
          else
          {
            itr->MapData( NODE_TO_SPREAD.GetFeatureVector(), NODE_TO_SPREAD.GetResultVector(), 0.2 * learning_rate);
          }
        }
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ApproximatorKohonenNetwork::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "A kohonen-network based predictor. See http://en.wikipedia.org/wiki/Self-organizing_map"
      );
      parameters.Merge( m_Schedule.GetSerializer());
      parameters.AddInitializer
      (
        "objective function",
        "function that evaluates the model after each batch step",
        io::Serialization::GetAgent( &m_ObjectiveFunction->GetImplementation()),
        "RMSD"
      );
      parameters.AddInitializer
      (
        "map dimensions",
        "size of each dimension, grid spacing for each node will always be 1.0",
        io::Serialization::GetAgentContainerWithCheck
        (
          &m_Network.GetMapDimensions(),
          io::Serialization::GetAgentWithRange( double( 0.0), double( 10000.0))
        )
      );
      parameters.AddInitializer
      (
        "steps per update",
        "# of features seen between each update of nodes (set to 0 to use the size of the training data set)",
        io::Serialization::GetAgent( &m_UpdateEveryNthFeature),
        "0"
      );
      parameters.AddInitializer
      (
        "length",
        "# of iterations it takes for radius to decrease to 0",
        io::Serialization::GetAgent( &m_Length),
        "10"
      );
      parameters.AddInitializer
      (
        "radius",
        "initial radius; larger radii are often better, at the expense of training time, up to about half the distance from one end of the map to the other",
        io::Serialization::GetAgent( &m_Radius),
        "10"
      );
      parameters.AddInitializer
      (
        "neighbor kernel",
        "Determines neighborhood type.  With guassian the influence decreases with distance (within the radii), with bubble it remains constant",
        io::Serialization::GetAgent( &m_NeighborKernel),
        "Bubble"
      );
      parameters.AddInitializer
      (
        "initializer",
        "Determines how the map is initialized.  FirstVectors chooses the first N vectors to populate the map, other methods allow selection of random training vectors or random training vector elements",
        io::Serialization::GetAgent( &m_Initializer),
        "FirstVectors"
      );
      parameters.AddInitializer
      (
        "scaling",
        "Type of input scaling. Normally AveStd works best, but MinMax and None may also be used in some circumstances",
        io::Serialization::GetAgent( &m_RescaleType),
        "MinMax"
      );
      parameters.AddInitializer
      (
        "cutoff",
        "Cutoff between actives and inactives. Needed if balancing",
        io::Serialization::GetAgent( &m_Cutoff),
        "0.5"
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
