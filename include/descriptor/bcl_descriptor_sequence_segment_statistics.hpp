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

// include header of this class
#include "bcl_descriptor_sequence_segment_statistics.h"
// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
#include "linal/bcl_linal_vector_operations.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief default constructor from whether to be springy at the end and whether to do a cumulative average
    template< typename t_DataType>
    SequenceSegmentStatistics< t_DataType>::SequenceSegmentStatistics() :
      m_ConditionalDescriptor(),
      m_Descriptor(),
      m_Statistics(),
      m_StatisticsSizes(),
      m_StatisticsSizePerCondition( 0),
      m_Conditions()
    {
    }

    //! @brief Clone function
    //! @return pointer to new SequenceSegmentStatistics
    template< typename t_DataType>
    SequenceSegmentStatistics< t_DataType> *SequenceSegmentStatistics< t_DataType>::Clone() const
    {
      return new SequenceSegmentStatistics( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &SequenceSegmentStatistics< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &SequenceSegmentStatistics< t_DataType>::GetAlias() const
    {
      static const std::string s_window_name( "SequenceSegmentStatistics");
      return s_window_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t SequenceSegmentStatistics< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return m_StatisticsSizePerCondition * ( 1 + m_Conditions.GetSize());
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

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer SequenceSegmentStatistics< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes statistics of segments across the sequence -- regions of the sequence that have the same value of a descriptor, which is "
        "termed the condition. An additional descriptor can be averaged over the elements of the segments and/or the "
        "sequence. Segment statistics for particular conditions of interest can also be computed"
      );
      parameters.AddInitializer
      (
        "descriptor",
        "descriptor to take statistics of",
        io::Serialization::GetAgent( &m_Descriptor)
      );
      parameters.AddInitializer
      (
        "condition",
        "descriptor that defines the segments; adjacent elements in the sequence with the same value of this descriptor"
        " are grouped into a segment",
        io::Serialization::GetAgent( &m_ConditionalDescriptor)
      );
      parameters.AddOptionalInitializer
      (
        "statistics",
        "List of statistics to calculate. If no statistics are given, all of them will be added",
        io::Serialization::GetAgent( &m_Statistics)
      );
      parameters.AddOptionalInitializer
      (
        "conditions",
        "values returned by the condition descriptor to compute separate sequence segment statistics for. "
        "Overall segment statistics will be returned regardless of the setting for this",
        io::Serialization::GetAgent( &m_Conditions)
      );
      return parameters;
    } // GetSerializer

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType>
    void SequenceSegmentStatistics< t_DataType>::Calculate
    (
      linal::VectorReference< float> &STORAGE
    )
    {
      // create segments over the sequence
      SegmentFinder segment_finder;
      for( Iterator< t_DataType> itr( this->GetCurrentObject()->GetIterator()); itr.NotAtEnd(); ++itr)
      {
        segment_finder.ConsiderNextElementOfSequence
        (
          m_ConditionalDescriptor->operator()( itr),
          m_Descriptor->operator()( itr)
        );
      }
      segment_finder.Finalize();

      size_t position_in_storage( 0);
      {
        storage::Vector< size_t>::const_iterator itr_stat_size( m_StatisticsSizes.Begin());
        // walk over the segment statistics; get values for each of them
        for
        (
          storage::Vector< SegmentFinder::StatisticEnum>::const_iterator
            itr( m_Statistics.Begin()), itr_end( m_Statistics.End());
          itr != itr_end;
          ++itr, ++itr_stat_size
        )
        {
          linal::VectorReference< float> store( STORAGE.CreateSubVectorReference( *itr_stat_size, position_in_storage));
          segment_finder.CopyStatisticConsideringAllSegments( store, *itr);
          position_in_storage += *itr_stat_size;
        }
      }

      // iterate over conditions, collect statistics for them too
      for
      (
        storage::Vector< linal::Vector< float> >::const_iterator
          itr_condition( m_Conditions.Begin()), itr_condition_end( m_Conditions.End());
        itr_condition != itr_condition_end;
        ++itr_condition
      )
      {
        storage::Vector< size_t>::const_iterator itr_stat_size( m_StatisticsSizes.Begin());
        // walk over the segment statistics; get values for each of them
        for
        (
          storage::Vector< SegmentFinder::StatisticEnum>::const_iterator
            itr( m_Statistics.Begin()), itr_end( m_Statistics.End());
          itr != itr_end;
          ++itr, ++itr_stat_size
        )
        {
          linal::VectorReference< float> store( STORAGE.CreateSubVectorReference( *itr_stat_size, position_in_storage));
          segment_finder.CopyStatisticConsideringSegmentsWithCondition( store, *itr, *itr_condition);
          position_in_storage += *itr_stat_size;
        }
      }
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > SequenceSegmentStatistics< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_ConditionalDescriptor, &m_Descriptor + 1);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool SequenceSegmentStatistics< t_DataType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( m_Descriptor->GetType().GetDimension() != size_t( 1))
      {
        ERR_STREAM << "Expected an element-wise descriptor for descriptor";
        return false;
      }
      if( m_ConditionalDescriptor->GetType().GetDimension() != size_t( 1))
      {
        ERR_STREAM << "Expected an element-wise descriptor for conditioned descriptor";
        return false;
      }
      m_StatisticsSizePerCondition = 0;
      // handle special case where no descriptors were given; in this case, add all statistics
      if( m_Statistics.IsEmpty())
      {
        for( SegmentFinder::StatisticEnum statistic; statistic < SegmentFinder::s_NumberStatistics; ++statistic)
        {
          m_Statistics.PushBack( statistic);
        }
      }

      const size_t condition_size( m_ConditionalDescriptor->GetSizeOfFeatures());
      const size_t descriptor_size( m_Descriptor->GetSizeOfFeatures());

      // walk over the selected statistics; partition them into window and raw statistics
      for
      (
        storage::Vector< SegmentFinder::StatisticEnum>::const_iterator
          itr( m_Statistics.Begin()), itr_end( m_Statistics.End());
        itr != itr_end;
        ++itr
      )
      {
        // determine the number of values in the descriptor
        const size_t stat_size
        (
          SegmentFinder::GetIsStatisticOfCondition( *itr)
          ? condition_size
          : SegmentFinder::GetIsStatisticOfDescriptor( *itr)
            ? descriptor_size
            : size_t( 1)
        );
        m_StatisticsSizes.PushBack( stat_size);
        m_StatisticsSizePerCondition += stat_size;
      }

      return true;
    }

  } // namespace descriptor
} // namespace bcl
