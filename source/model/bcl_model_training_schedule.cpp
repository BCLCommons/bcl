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
#include "model/bcl_model_training_schedule.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_reflecting.h"
#include "model/bcl_model_retrieve_data_set_from_delimited_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from members
    TrainingSchedule::TrainingSchedule
    (
      const bool &SHUFFLE,
      const bool &BALANCE,
      const size_t &MAX_REPEATS,
      const float &TARGET_OVERSAMPLING
    ) :
      m_Balance( BALANCE),
      m_Shuffle( SHUFFLE),
      m_BalanceMaxRepeatedFeatures( MAX_REPEATS),
      m_BalanceMaxOversampling( TARGET_OVERSAMPLING)
    {
    }

    //! @brief Clone function
    //! @return pointer to new TrainingSchedule
    TrainingSchedule *TrainingSchedule::Clone() const
    {
      return new TrainingSchedule( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &TrainingSchedule::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &TrainingSchedule::GetAlias() const
    {
      static const std::string s_name( "TrainingScheduler");
      return s_name;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief Setup the balancer
    //! @param RESULTS results feature data set
    //! @param CUTOFF cutoff from the objective function that separates classes
    void TrainingSchedule::Setup
    (
      const FeatureDataSetInterface< float> &RESULTS,
      const float &CUTOFF
    )
    {
      const linal::MatrixConstReference< float> results( RESULTS.GetMatrix());

      const size_t results_size( RESULTS.GetFeatureSize());
      const size_t dataset_size( RESULTS.GetNumberFeatures());

      if( util::IsDefined( CUTOFF))
      {
        // Rescaled cutoffs, one for each result
        linal::Vector< float> scaled_cutoffs( results_size, CUTOFF);

        // compute the rescaled cutoffs
        if( RESULTS.IsRescaled())
        {
          const util::SiPtr< const RescaleFeatureDataSet> results_scaling( RESULTS.GetScaling());
          for( size_t result( 0); result < results_size; ++result)
          {
            scaled_cutoffs( result) = results_scaling->RescaleValue( result, CUTOFF);
          }
        }

        const size_t bitsize( sizeof( size_t) * 8);
        const size_t n_sizets( ( results_size - 1) / bitsize + 1);
        storage::Vector< size_t> result_class( n_sizets, size_t( 0));
        m_ResultClass = storage::Vector< size_t>( RESULTS.GetNumberFeatures(), size_t( 0));
        // for each result, compute the result class
        size_t n_classes_so_far( 0);
        storage::Map< storage::Vector< size_t>, size_t> class_to_index;
        storage::Vector< storage::Vector< size_t> > classes;
        for( size_t feature( 0); feature < dataset_size; ++feature)
        {
          result_class.SetAllElements( 0);
          linal::VectorConstReference< float> result_row( results.GetRow( feature));
          for( size_t result_index( 0); result_index < results_size; ++result_index)
          {
            if( result_row( result_index) >= scaled_cutoffs( result_index))
            {
              result_class( result_index / bitsize) |= ( 1 << ( result_index % bitsize));
            }
          }
          storage::Map< storage::Vector< size_t>, size_t>::const_iterator itr_map( class_to_index.Find( result_class));
          if( itr_map == class_to_index.End())
          {
            class_to_index[ result_class] = m_ResultClass( feature) = n_classes_so_far++;
            classes.PushBack( result_class);
          }
          else
          {
            m_ResultClass( feature) = itr_map->second;
          }
        }
        // accumulate all classes
        m_PeerFeatures = storage::Vector< storage::Vector< size_t> >( n_classes_so_far);
        for( size_t feature( 0); feature < dataset_size; ++feature)
        {
          m_PeerFeatures( m_ResultClass( feature)).PushBack( feature);
        }
        std::stringstream info_stream;
        info_stream << "Found " << n_classes_so_far << " classes of training points. Counts per class: ";
        for( size_t class_id( 0); class_id < n_classes_so_far; ++class_id)
        {
          info_stream << m_PeerFeatures( class_id).GetSize() << ' ';
        }
        info_stream << ".  Class Binary IDs: ";
        for( size_t class_id( 0); class_id < n_classes_so_far; ++class_id)
        {
          for( size_t result_index( 0); result_index < results_size; ++result_index)
          {
            info_stream << int( bool( classes( class_id)( result_index / bitsize) & ( 1 << ( result_index % bitsize))));
          }
          info_stream << ' ';
        }
        BCL_MessageStd( info_stream.str());
      }

      m_Order.Resize( dataset_size);
      for( size_t i( 0); i < dataset_size; ++i)
      {
        m_Order( i) = i;
      }

      if( m_Balance)
      {
        // get scaling and cutoff information
        if( !util::IsDefined( CUTOFF))
        {
          BCL_MessageCrt( "Cannot balance a dataset when a non-classification type objective is used!");
          return;
        }
        // create reflecting iterators for each class and determine the size of the most populated class
        const size_t n_classes( m_PeerFeatures.GetSize());
        BCL_Assert( n_classes > 1, "Only a single class of features was found; balancing is impossible!");
        size_t max_class_size( m_PeerFeatures.FirstElement().GetSize());
        storage::Vector< iterate::Reflecting< const size_t> > class_iterators;
        for( size_t class_id( 0); class_id < n_classes; ++class_id)
        {
          class_iterators.PushBack
          (
            iterate::Reflecting< const size_t>( m_PeerFeatures( class_id).Begin(), m_PeerFeatures( class_id).End())
          );
          max_class_size = std::max( max_class_size, m_PeerFeatures( class_id).GetSize());
        }
        size_t second_most_common_class_size( 0);
        for( size_t class_id( 0); class_id < n_classes; ++class_id)
        {
          if( m_PeerFeatures( class_id).GetSize() != max_class_size)
          {
            second_most_common_class_size = std::max( second_most_common_class_size, m_PeerFeatures( class_id).GetSize());
          }
        }

        size_t max_target_size( std::min( size_t( max_class_size * m_BalanceMaxOversampling + 1), max_class_size));
        storage::Vector< size_t> max_this_class( n_classes, max_class_size);
        float oversampling_factor
        (
          std::min( float( m_BalanceMaxRepeatedFeatures), float( max_target_size) / float( second_most_common_class_size))
        );
        if( oversampling_factor > 1.0)
        {
          for( size_t class_id( 0); class_id < n_classes; ++class_id)
          {
            if( m_PeerFeatures( class_id).GetSize() < max_class_size)
            {
              max_this_class( class_id) =
                std::min( m_PeerFeatures( class_id).GetSize() * oversampling_factor, float( max_target_size));
            }
          }
          m_Order.Reset();
          m_Order.AllocateMemory( n_classes * max_class_size);
          for( size_t class_feature( 0); class_feature < max_class_size; ++class_feature)
          {
            for( size_t class_id( 0); class_id < n_classes; ++class_id)
            {
              if( class_feature <= max_this_class( class_id))
              {
                m_Order.PushBack( *class_iterators( class_id));
                ++class_iterators( class_id);
              }
            }
          }
        }
      }
      //m_Order.Shuffle();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer TrainingSchedule::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Chooses the order in which training examples are visited, optionally to balance classes"
      );
      parameters.AddInitializer
      (
        "shuffle",
        "primarily for non-batch update; if true, shuffle the order or data points between each run through the data",
        io::Serialization::GetAgent( &m_Shuffle),
        "False"
      );
      parameters.AddInitializer
      (
        "balance",
        "Whether to automatically balance each class (as defined by the objective function's cutoff, if applicable)",
        io::Serialization::GetAgent( &m_Balance),
        "False"
      );
      parameters.AddInitializer
      (
        "balance max repeats",
        "Applies only if balance=True; absolute maximum number of times that a feature can be repeated in order to reach"
        "the targeted ratio of positives to negatives",
        io::Serialization::GetAgent( &m_BalanceMaxRepeatedFeatures),
        "1000000"
      );
      parameters.AddInitializer
      (
        "balance target ratio",
        "Applies only if balance=True; target ratio between most-common and underrepresented class in the dataset "
        "achieved by data replication; to simulate normal balancing; this should be 1, but smaller values may yield "
        "more general models",
        io::Serialization::GetAgentWithRange( &m_BalanceMaxOversampling, float( 0.0), float( 1.0)),
        "1"
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
