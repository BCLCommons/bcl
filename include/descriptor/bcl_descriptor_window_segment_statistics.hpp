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
#include "bcl_descriptor_window_segment_statistics.h"
// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
#include "iterate/bcl_iterate_reflecting.h"
#include "linal/bcl_linal_vector_operations.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief default constructor from whether to be springy at the end and whether to do a cumulative average
    template< typename t_DataType>
    WindowSegmentStatistics< t_DataType>::WindowSegmentStatistics() :
      m_ConditionalDescriptor(),
      m_ValueDescriptor(),
      m_SegmentsRadius( 0),
      m_IgnoreUndefinedConditionSegments( false),
      m_Statistics(),
      m_SegmentFinder(),
      m_SegmentWindowValueCount( 0)
    {
    }

    //! @brief Clone function
    //! @return pointer to new WindowSegmentStatistics
    template< typename t_DataType>
    WindowSegmentStatistics< t_DataType> *WindowSegmentStatistics< t_DataType>::Clone() const
    {
      return new WindowSegmentStatistics( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &WindowSegmentStatistics< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &WindowSegmentStatistics< t_DataType>::GetAlias() const
    {
      static const std::string s_window_name( "WindowSegmentStatistics");
      return s_window_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t WindowSegmentStatistics< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return ( 2 * m_SegmentsRadius + 1) * m_SegmentWindowValueCount;
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
    io::Serializer WindowSegmentStatistics< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "computes statistics over segments -- regions of the sequence that have the same value of a descriptor, which is "
        "termed the condition. An additional descriptor can be averaged over the elements of the segments and/or the "
        "sequence"
      );
      parameters.AddInitializer
      (
        "radius",
        "maximum number of segments to consider in each direction beyond the current segment",
        io::Serialization::GetAgent( &m_SegmentsRadius)
      );
      parameters.AddInitializer
      (
        "descriptor",
        "descriptor to take statistics of",
        io::Serialization::GetAgent( &m_ValueDescriptor)
      );
      parameters.AddInitializer
      (
        "condition",
        "descriptor that defines the segments; adjacent elements in the sequence with the same value of this descriptor"
        " are grouped into a segment",
        io::Serialization::GetAgent( &m_ConditionalDescriptor)
      );
      parameters.AddInitializer
      (
        "ignore undefined elements",
        "If set, elements for which the conditional descriptor is not defined will be ignored (not counted towards any "
        "statistics or segment sizes). The default behavior is that undefined values are treated like other numbers, so "
        "regions where the condition is undefined will normally be distinct segments from the surrounding regions. This "
        "setting causes the elements to be explicitly ignored",
        io::Serialization::GetAgent( &m_IgnoreUndefinedConditionSegments),
        "False"
      );
      parameters.AddOptionalInitializer
      (
        "statistics",
        "List of statistics to calculate. If no statistics are given, all of them will be added",
        io::Serialization::GetAgent( &m_Statistics)
      );
      return parameters;
    } // GetSerializer

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    template< typename t_DataType>
    void WindowSegmentStatistics< t_DataType>::Calculate
    (
      const iterate::Generic< const t_DataType> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      if( m_SegmentFinder.IsEmpty())
      {
        RunSegmentFinder();
      }
      if( m_SegmentFinder.IsEmpty())
      {
        // no segments
        STORAGE = float( 0);
        return;
      }

      size_t element_position( ELEMENT.GetPosition());
      if( m_IgnoreUndefinedConditionSegments)
      {
        element_position = m_SegmentFinderIndex( element_position);
        if( !util::IsDefined( element_position))
        {
          // segment with undefined condition value, return undefined for all values
          STORAGE = util::GetUndefined< float>();
          return;
        }
      }

      SegmentFinder::const_iterator itr_seg( m_SegmentFinder.SegmentForPosition( element_position));
      size_t position_in_storage( 0);
      // walk to the left
      iterate::Reflecting< const SegmentInfo> itr_segment_left
      (
        m_SegmentFinder.SegmentBegin(),
        m_SegmentFinder.SegmentEnd(),
        itr_seg
      );

      // iterate over segments to the left
      for
      (
        size_t number_segments_left( 0);
        number_segments_left <= m_SegmentsRadius;
        ++number_segments_left, --itr_segment_left
      )
      {
        storage::Vector< size_t>::const_iterator itr_stat_size( m_StatisticsSizes.Begin());
        // walk over the segment window statistics; get values for each of them
        for
        (
          storage::Vector< SegmentInfo::StatisticEnum>::const_iterator
            itr( m_Statistics.Begin()), itr_end( m_Statistics.End());
          itr != itr_end;
          ++itr, ++itr_stat_size
        )
        {
          linal::VectorReference< float> store( STORAGE.CreateSubVectorReference( *itr_stat_size, position_in_storage));
          itr_segment_left->CopyStatistic( store, *itr, element_position);
          position_in_storage += *itr_stat_size;
        }
      }
      // walk right
      iterate::Reflecting< const SegmentInfo> itr_segment_right
      (
        m_SegmentFinder.SegmentBegin(),
        m_SegmentFinder.SegmentEnd(),
        itr_seg
      );
      ++itr_segment_right;
      for
      (
        size_t number_segments_right( 1);
        number_segments_right <= m_SegmentsRadius;
        ++number_segments_right, ++itr_segment_right
      )
      {
        storage::Vector< size_t>::const_iterator itr_stat_size( m_StatisticsSizes.Begin());
        // walk over the segment window statistics; get values for each of them
        for
        (
          storage::Vector< SegmentInfo::StatisticEnum>::const_iterator
            itr( m_Statistics.Begin()), itr_end( m_Statistics.End());
          itr != itr_end;
          ++itr, ++itr_stat_size
        )
        {
          linal::VectorReference< float> store( STORAGE.CreateSubVectorReference( *itr_stat_size, position_in_storage));
          itr_segment_right->CopyStatistic( store, *itr, element_position);
          position_in_storage += *itr_stat_size;
        }
      }
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > WindowSegmentStatistics< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_ConditionalDescriptor, &m_ValueDescriptor + 1);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    template< typename t_DataType>
    void WindowSegmentStatistics< t_DataType>::SetObjectHook()
    {
      if( !m_SegmentFinder.IsEmpty())
      {
        m_SegmentFinder.Reset();
        m_SegmentFinderIndex.Reset();
      }
    }

    //! @brief run the segment finder over the sequence
    template< typename t_DataType>
    void WindowSegmentStatistics< t_DataType>::RunSegmentFinder()
    {
      if( m_IgnoreUndefinedConditionSegments)
      {
        // keep a vector map (m_SegmentFinderIndex) from the index in the sequence to the index in the segment finder
        m_SegmentFinderIndex.Reset();
        m_SegmentFinderIndex.AllocateMemory( this->GetCurrentObject()->GetSize());
        size_t defined_conditions( 0);
        for( Iterator< t_DataType> itr( this->GetCurrentObject()->GetIterator()); itr.NotAtEnd(); ++itr)
        {
          linal::VectorConstReference< float> condition( m_ConditionalDescriptor->operator()( itr));
          if( condition.IsDefined())
          {
            m_SegmentFinder.ConsiderNextElementOfSequence( condition, m_ValueDescriptor->operator()( itr));
            m_SegmentFinderIndex.PushBack( defined_conditions++);
          }
          else
          {
            m_SegmentFinderIndex.PushBack( util::GetUndefined< size_t>());
          }
        }
      }
      else
      {
        for( Iterator< t_DataType> itr( this->GetCurrentObject()->GetIterator()); itr.NotAtEnd(); ++itr)
        {
          m_SegmentFinder.ConsiderNextElementOfSequence
          (
            m_ConditionalDescriptor->operator()( itr),
            m_ValueDescriptor->operator()( itr)
          );
        }
      }
      m_SegmentFinder.Finalize();
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    template< typename t_DataType>
    bool WindowSegmentStatistics< t_DataType>::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERR_STREAM
    )
    {
      if( m_ValueDescriptor->GetType().GetDimension() != size_t( 1))
      {
        ERR_STREAM << "Expected an element-wise descriptor for descriptor";
        return false;
      }
      if( m_ConditionalDescriptor->GetType().GetDimension() != size_t( 1))
      {
        ERR_STREAM << "Expected an element-wise descriptor for conditioned descriptor";
        return false;
      }
      m_SegmentWindowValueCount = 0;
      // handle special case where no descriptors were given; in this case, add all statistics
      if( m_Statistics.IsEmpty())
      {
        for( SegmentInfo::StatisticEnum statistic; statistic < SegmentInfo::s_NumberStatistics; ++statistic)
        {
          m_Statistics.PushBack( statistic);
        }
      }

      const size_t condition_size( m_ConditionalDescriptor->GetSizeOfFeatures());
      const size_t descriptor_size( m_ValueDescriptor->GetSizeOfFeatures());
      m_StatisticsSizes.AllocateMemory( m_Statistics.GetSize());

      // walk over the selected statistics; partition them into window and raw statistics
      for
      (
        storage::Vector< SegmentInfo::StatisticEnum>::const_iterator
          itr( m_Statistics.Begin()), itr_end( m_Statistics.End());
        itr != itr_end;
        ++itr
      )
      {
        // determine the number of values in the descriptor
        const size_t stat_size
        (
          *itr == SegmentInfo::e_Condition
          ? condition_size
          : *itr == SegmentInfo::e_DescriptorAve || *itr == SegmentInfo::e_DescriptorSD
            || *itr == SegmentInfo::e_DescriptorMin || *itr == SegmentInfo::e_DescriptorMax
            || *itr == SegmentInfo::e_DescriptorFirst || *itr == SegmentInfo::e_DescriptorLast
            ? descriptor_size
            : size_t( 1)
        );
        m_StatisticsSizes.PushBack( stat_size);
        m_SegmentWindowValueCount += stat_size;
      }

      return true;
    }

  } // namespace descriptor
} // namespace bcl
