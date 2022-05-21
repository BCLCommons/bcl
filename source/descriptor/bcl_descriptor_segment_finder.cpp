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
#include "descriptor/bcl_descriptor_segment_finder.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief Statistic as string
    //! @param STATISTIC the statistic
    //! @return the Statistic as string
    const std::string &SegmentFinder::GetStatisticString( const Statistic &STATISTIC)
    {
      static const std::string s_names[ s_NumberStatistics + 1] =
      {
        "NumSegments",
        "FractionSegments",
        "FractionElements",
        "LengthSegmentsAve",
        "LengthSegmentsSD",
        "LengthSegmentsElementwiseAve",
        "LengthSegmentsElementwiseSD",
        "DescriptorSegmentsAve",
        "DescriptorSegmentsSD",
        "DescriptorSegmentsElementwiseAve",
        "DescriptorSegmentsElementwiseSD",
        "ConditionSegmentsAve",
        "ConditionSegmentsSD",
        "ConditionSegmentsElementwiseAve",
        "ConditionSegmentsElementwiseSD",
        GetStaticClassName< SegmentInfo>()
      };

      return s_names[ STATISTIC];
    }

    //! @brief test whether the statistic returns a value for every value in the condition
    //! @param STATISTIC the statistic of interest
    //! @return true if the statistic returns a value for every value in the condition
    bool SegmentFinder::GetIsStatisticOfCondition( const Statistic &STATISTIC)
    {
      return STATISTIC == e_ConditionSegmentsAve
             || STATISTIC == e_ConditionSegmentsSD
             || STATISTIC == e_ConditionSegmentsElementwiseAve
             || STATISTIC == e_ConditionSegmentsElementwiseSD;
    }

    //! @brief test whether the statistic returns a value for every value in the statistic descriptor
    //! @param STATISTIC the statistic of interest
    //! @return true if the statistic returns a value for every value in the statistic descriptor
    bool SegmentFinder::GetIsStatisticOfDescriptor( const Statistic &STATISTIC)
    {
      return STATISTIC == e_DescriptorSegmentsAve
             || STATISTIC == e_DescriptorSegmentsSD
             || STATISTIC == e_DescriptorSegmentsElementwiseAve
             || STATISTIC == e_DescriptorSegmentsElementwiseSD;
    }

    //! @brief function that compares two linal::Vectors, ignoring effectively considering undefined values to be FLT_MAX
    //! @param A, B the vectors to compare
    //! @return A < B, respective of undefined values
    //! This is necessary because maps require that their keys have weak ordering; which is not possible if the values
    //! are undefined
    bool SegmentFinder::VectorLessThanRespectUndefined::operator()
    (
      const linal::Vector< float> &A,
      const linal::Vector< float> &B
    ) const
    {
      if( A.GetSize() != B.GetSize())
      {
        return A.GetSize() < B.GetSize();
      }
      // both vectors have undefined values, need to compare them element by element
      for
      (
        const float *itr_a( A.Begin()), *itr_b( B.Begin()), *itr_a_end( A.End());
        itr_a != itr_a_end;
        ++itr_a, ++itr_b
      )
      {
        const bool a_def( util::IsDefined( *itr_a)), b_def( util::IsDefined( *itr_b));
        if( a_def != b_def)
        {
          return a_def;
        }
        else if( a_def && *itr_a != *itr_b)
        {
          return *itr_a < *itr_b;
        }
        // both undefined, just continue
      }
      // same pattern of undefineds
      return false;
    }

    //! @brief simpler test for ==
    //! @param A, B the vectors to compare
    //! @return A == B, respective of undefined values
    bool SegmentFinder::VectorLessThanRespectUndefined::Equal( const linal::Vector< float> &A, const linal::Vector< float> &B)
    {
      if( A.GetSize() != B.GetSize())
      {
        return false;
      }
      // both vectors have undefined values, need to compare them element by element
      for
      (
        const float *itr_a( A.Begin()), *itr_b( B.Begin()), *itr_a_end( A.End());
        itr_a != itr_a_end;
        ++itr_a, ++itr_b
      )
      {
        const bool a_def( util::IsDefined( *itr_a)), b_def( util::IsDefined( *itr_b));
        if( a_def != b_def || ( a_def && *itr_a != *itr_b))
        {
          return false;
        }
      }
      // same pattern of undefineds
      return true;
    }

    //! @brief default constructor
    SegmentFinder::SegmentFinder() :
      m_UnconditionalSequenceStatistics(),
      m_CurrentConditionToSequenceStatisticsIterator( m_ConditionToSequenceStatistics.End())
    {
      m_UnconditionalSequenceStatistics.Third() = 0;
    }

    //! @brief copy constructor
    SegmentFinder::SegmentFinder( const SegmentFinder &PARENT) :
      m_UnconditionalSequenceStatistics(),
      m_CurrentConditionToSequenceStatisticsIterator( m_ConditionToSequenceStatistics.End())
    {
      m_UnconditionalSequenceStatistics.Third() = 0;
      BCL_Assert
      (
        PARENT.GetNumberSegments() == size_t( 0),
        "Cannot copy a non-empty SegmentFinder because it contains objects that hold pointers to internally-held data"
      );
    }

    //! @brief Clone function
    //! @return pointer to new SegmentFinder
    SegmentFinder *SegmentFinder::Clone() const
    {
      return new SegmentFinder( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SegmentFinder::GetClassIdentifier() const
    {
      return GetStaticClassName< SegmentFinder>();
    }

    //! @brief get an iterator to the first segment
    //! @return an iterator to the first segment
    SegmentFinder::const_iterator SegmentFinder::SegmentBegin() const
    {
      return m_SequenceSegments.Begin();
    }

    //! @brief get an iterator beyond the last segment
    //! @return an iterator beyond the last segment
    SegmentFinder::const_iterator SegmentFinder::SegmentEnd() const
    {
      return m_SequenceSegments.End();
    }

    //! @brief get an iterator for the element at a particular position in the original sequence
    //! @param POSITION position in the original sequence to get an iterator for
    //! @return an iterator for the element at a particular position in the original sequence
    SegmentFinder::const_iterator SegmentFinder::SegmentForPosition( const size_t &POSITION) const
    {
      return m_IteratorsForSequence[ POSITION];
    }

    //! @brief test whether there all segments' statistics have been accumulated
    //! @return true if the sequence has been finalized
    bool SegmentFinder::IsFinalized() const
    {
      return m_CurrentConditionToSequenceStatisticsIterator == m_ConditionToSequenceStatistics.End();
    }

    //! @brief get the number of segments found so far
    //! @return the number of segments found so far
    size_t SegmentFinder::GetNumberSegments() const
    {
      return m_SequenceSegments.GetSize();
    }

    //! @brief test whether the segment finder is empty
    //! @return true if the segment finder is empty
    bool SegmentFinder::IsEmpty() const
    {
      return m_SequenceSegments.IsEmpty();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add information for the next element in the sequence to this class
    //! @param CONDITION current value of the condition
    //! @param DESCRIPTOR current value of the descriptor
    void SegmentFinder::ConsiderNextElementOfSequence
    (
      const linal::VectorConstInterface< float> &CONDITION,
      const linal::VectorConstInterface< float> &DESCRIPTOR
    )
    {
      // handle the case that the condition changed
      if
      (
        m_SequenceSegments.IsEmpty()
        || !VectorLessThanRespectUndefined::Equal( m_CurrentConditionToSequenceStatisticsIterator->first, CONDITION)
      )
      {
        AddEmptySegment( CONDITION);
      }
      SegmentInfo &current_segment( m_SequenceSegments.LastElement());
      m_IteratorsForSequence.push_back( m_SequenceSegments.Last());
      current_segment.PushBack( DESCRIPTOR);
    }

    //! @brief finalize; flushes most recent segment's statistics into the various maps
    void SegmentFinder::Finalize()
    {
      if( IsFinalized())
      {
        // nothing to do; already finalized
        return;
      }

      // get the current element
      const SegmentInfo &current_segment( m_SequenceSegments.LastElement());

      // add unweighted statistics
      const float length( current_segment.GetLength());
      const linal::Vector< float> &average( current_segment.GetDescriptorAverageSD().GetAverage());
      m_CurrentConditionToSequenceStatisticsIterator->second.First().First() += average;
      m_CurrentConditionToSequenceStatisticsIterator->second.Second().First() += length;
      m_CurrentConditionToSequenceStatisticsIterator->second.First().Second().AddWeightedObservation( average, length);
      m_CurrentConditionToSequenceStatisticsIterator->second.Second().Second().AddWeightedObservation( length, length);

      // update the overall sequence statistics information, unweighted
      const linal::Vector< float> &condition( current_segment.GetCondition());
      m_UnconditionalSequenceStatistics.First().First() += average;
      m_UnconditionalSequenceStatistics.First().Second() += condition;
      m_UnconditionalSequenceStatistics.Second() += length;

      // update the weighted sequence statistics
      m_UnconditionalSequenceElementwiseStatistics.First().First().AddWeightedObservation( average, length);
      m_UnconditionalSequenceElementwiseStatistics.First().Second().AddWeightedObservation( condition, length);
      m_UnconditionalSequenceElementwiseStatistics.Second().AddWeightedObservation( length, length);

      // set the iterator to the end; this indicates that the segment finder is currently fully updated
      m_CurrentConditionToSequenceStatisticsIterator = m_ConditionToSequenceStatistics.End();
    }

    //! @brief reset -- resets this class so that it forgets about all the segments
    void SegmentFinder::Reset()
    {
      m_SequenceSegments.Reset();
      m_ConditionToSequenceStatistics.Reset();
      m_UnconditionalSequenceStatistics.First().First().Reset();
      m_UnconditionalSequenceStatistics.First().Second().Reset();
      m_UnconditionalSequenceStatistics.Second().Reset();
      m_UnconditionalSequenceStatistics.Third() = 0;
      m_IteratorsForSequence.clear();
      m_UnconditionalSequenceElementwiseStatistics.First().First().Reset();
      m_UnconditionalSequenceElementwiseStatistics.First().Second().Reset();
      m_UnconditionalSequenceElementwiseStatistics.Second().Reset();
      m_CurrentConditionToSequenceStatisticsIterator = m_ConditionToSequenceStatistics.End();
    }

    //! @brief copy a statistic into a vector
    //! @param STORAGE storage for the statistic
    //! @param STATISTIC the statistic desired
    void SegmentFinder::CopyStatisticConsideringAllSegments
    (
      linal::VectorInterface< float> &REFERENCE,
      const Statistic &STATISTIC
    ) const
    {
      switch( STATISTIC)
      {
        case e_NumSegments:
          REFERENCE( 0) = m_UnconditionalSequenceStatistics.Third();
          break;
        case e_FractionElements:
        case e_FractionSegments:
          REFERENCE( 0) = 1.0;
          break;
        case e_LengthSegmentsAve:
          REFERENCE( 0) = m_UnconditionalSequenceStatistics.Second().GetAverage();
          break;
        case e_LengthSegmentsSD:
          REFERENCE( 0) = m_UnconditionalSequenceStatistics.Second().GetStandardDeviation();
          break;
        case e_LengthSegmentsElementwiseAve:
          REFERENCE( 0) = m_UnconditionalSequenceElementwiseStatistics.Second().GetAverage();
          break;
        case e_LengthSegmentsElementwiseSD:
          REFERENCE( 0) = m_UnconditionalSequenceElementwiseStatistics.Second().GetStandardDeviation();
          break;
        case e_DescriptorSegmentsAve:
          std::copy
          (
            m_UnconditionalSequenceStatistics.First().First().GetAverage().Begin(),
            m_UnconditionalSequenceStatistics.First().First().GetAverage().End(),
            REFERENCE.Begin()
          );
          break;
        case e_DescriptorSegmentsSD:
          std::copy
          (
            m_UnconditionalSequenceStatistics.First().First().GetStandardDeviation().Begin(),
            m_UnconditionalSequenceStatistics.First().First().GetStandardDeviation().End(),
            REFERENCE.Begin()
          );
          break;
        case e_DescriptorSegmentsElementwiseAve:
          std::copy
          (
            m_UnconditionalSequenceElementwiseStatistics.First().First().GetAverage().Begin(),
            m_UnconditionalSequenceElementwiseStatistics.First().First().GetAverage().End(),
            REFERENCE.Begin()
          );
          break;
        case e_DescriptorSegmentsElementwiseSD:
          std::copy
          (
            m_UnconditionalSequenceElementwiseStatistics.First().First().GetStandardDeviation().Begin(),
            m_UnconditionalSequenceElementwiseStatistics.First().First().GetStandardDeviation().End(),
            REFERENCE.Begin()
          );
          break;
        case e_ConditionSegmentsAve:
          std::copy
          (
            m_UnconditionalSequenceStatistics.First().Second().GetAverage().Begin(),
            m_UnconditionalSequenceStatistics.First().Second().GetAverage().End(),
            REFERENCE.Begin()
          );
          break;
        case e_ConditionSegmentsSD:
          std::copy
          (
            m_UnconditionalSequenceStatistics.First().Second().GetStandardDeviation().Begin(),
            m_UnconditionalSequenceStatistics.First().Second().GetStandardDeviation().End(),
            REFERENCE.Begin()
          );
          break;
        case e_ConditionSegmentsElementwiseAve:
          std::copy
          (
            m_UnconditionalSequenceElementwiseStatistics.First().Second().GetAverage().Begin(),
            m_UnconditionalSequenceElementwiseStatistics.First().Second().GetAverage().End(),
            REFERENCE.Begin()
          );
          break;
        case e_ConditionSegmentsElementwiseSD:
          std::copy
          (
            m_UnconditionalSequenceElementwiseStatistics.First().Second().GetStandardDeviation().Begin(),
            m_UnconditionalSequenceElementwiseStatistics.First().Second().GetStandardDeviation().End(),
            REFERENCE.Begin()
          );
          break;
        case s_NumberStatistics:
        default:
          BCL_Exit( "Invalid statistic: " + GetStatisticString( STATISTIC), -1);
          break;
      }
    }

    //! @brief copy a statistic into a vector
    //! @param STORAGE storage for the statistic
    //! @param STATISTIC the statistic desired
    //! @param CONDITION the condition of interest
    void SegmentFinder::CopyStatisticConsideringSegmentsWithCondition
    (
      linal::VectorInterface< float> &REFERENCE,
      const Statistic &STATISTIC,
      const linal::VectorConstInterface< float> &CONDITION
    ) const
    {
      // handle the trivial case where the statistic is only dependent on the condition
      // since this function only considers the condition value passed in, the average and standard deviation are trivial
      if( GetIsStatisticOfCondition( STATISTIC))
      {
        if( STATISTIC == e_ConditionSegmentsAve || STATISTIC == e_ConditionSegmentsElementwiseAve)
        {
          std::copy( CONDITION.Begin(), CONDITION.End(), REFERENCE.Begin());
        }
        else // if( STATISTIC == e_ConditionSegmentsSD || STATISTIC == e_ConditionSegmentsElementwiseSD)
        {
          REFERENCE = float( 0.0);
        }
        return;
      }
      // locate the statistic in the map
      storage::Map
      <
        linal::Vector< float>,
        storage::Triplet
        <
          storage::VectorND< 2, math::RunningAverageSD< linal::Vector< float> > >,
          storage::VectorND< 2, math::RunningAverageSD< float> >,
          size_t
        >,
        VectorLessThanRespectUndefined
      >::const_iterator itr( m_ConditionToSequenceStatistics.Find( CONDITION));
      if( itr == m_ConditionToSequenceStatistics.End())
      {
        // the condition has never been encountered
        // set all elements to 0.0, unless the statistic is descriptor-specific, in which case, use the sequence averages
        if( GetIsStatisticOfDescriptor( STATISTIC))
        {
          CopyStatisticConsideringAllSegments( REFERENCE, STATISTIC);
        }
        else
        {
          REFERENCE = float( 0.0);
        }
        return;
      }

      // from here on, STATISTIC is not a condition descriptor, and statistics are known
      const storage::Triplet
      <
        storage::VectorND< 2, math::RunningAverageSD< linal::Vector< float> > >,
        storage::VectorND< 2, math::RunningAverageSD< float> >,
        size_t
      > &stats( itr->second);
      switch( STATISTIC)
      {
        case e_NumSegments:
          REFERENCE( 0) = stats.Third();
          break;
        case e_FractionSegments:
          REFERENCE( 0) = float( stats.Third()) / float( m_UnconditionalSequenceStatistics.Third());
          break;
        case e_FractionElements:
          // use the weight from the residue-weighted length for the local condition divided by the weight from the
          // residue weighted length for the overall sequence.
          REFERENCE( 0) = float( stats.Second().Second().GetWeight()) // weight for the residue-weighted length
                          /
                          float( m_UnconditionalSequenceElementwiseStatistics.Second().GetWeight());
          break;
        case e_LengthSegmentsAve:
          REFERENCE( 0) = stats.Second().First().GetAverage();
          break;
        case e_LengthSegmentsSD:
          REFERENCE( 0) = stats.Second().First().GetStandardDeviation();
          break;
        case e_LengthSegmentsElementwiseAve:
          REFERENCE( 0) = stats.Second().Second().GetAverage();
          break;
        case e_LengthSegmentsElementwiseSD:
          REFERENCE( 0) = stats.Second().Second().GetStandardDeviation();
          break;
        case e_DescriptorSegmentsAve:
          std::copy
          (
            stats.First().First().GetAverage().Begin(),
            stats.First().First().GetAverage().End(),
            REFERENCE.Begin()
          );
          break;
        case e_DescriptorSegmentsSD:
          std::copy
          (
            stats.First().First().GetStandardDeviation().Begin(),
            stats.First().First().GetStandardDeviation().End(),
            REFERENCE.Begin()
          );
          break;
        case e_DescriptorSegmentsElementwiseAve:
          std::copy
          (
            stats.First().Second().GetAverage().Begin(),
            stats.First().Second().GetAverage().End(),
            REFERENCE.Begin()
          );
          break;
        case e_DescriptorSegmentsElementwiseSD:
          std::copy
          (
            stats.First().Second().GetStandardDeviation().Begin(),
            stats.First().Second().GetStandardDeviation().End(),
            REFERENCE.Begin()
          );
          break;
        default:
          BCL_Exit( "Invalid statistic: " + GetStatisticString( STATISTIC), -1);
          break;
      }
    }

    //! @brief get a statistic as a vector
    //! @param STATISTIC the statistic desired
    linal::Vector< float> SegmentFinder::GetStatisticConsideringAllSegments( const Statistic &STATISTIC) const
    {
      linal::Vector< float> value( GetStatisticSize( STATISTIC), float( 0));
      CopyStatisticConsideringAllSegments( value, STATISTIC);
      return value;
    }

    //! @brief get a statistic as a vector
    //! @param STATISTIC the statistic desired
    //! @param CONDITION the condition of interest
    linal::Vector< float> SegmentFinder::GetStatisticConsideringSegmentsWithCondition
    (
      const Statistic &STATISTIC,
      const linal::VectorConstInterface< float> &CONDITION
    ) const
    {
      linal::Vector< float> value( GetStatisticSize( STATISTIC), float( 0));
      CopyStatisticConsideringSegmentsWithCondition( value, STATISTIC, CONDITION);
      return value;
    }

    //! @brief get the size of the statistic
    //! @param STATISTIC the actual statistic of interest
    //! @return the size of the statistic
    size_t SegmentFinder::GetStatisticSize( const Statistic &STATISTIC) const
    {
      size_t expected_size( 1);
      if( GetIsStatisticOfCondition( STATISTIC))
      {
        expected_size = m_UnconditionalSequenceStatistics.First().Second().GetAverage().GetSize();
      }
      else if( GetIsStatisticOfDescriptor( STATISTIC))
      {
        expected_size = m_UnconditionalSequenceStatistics.First().First().GetAverage().GetSize();
      }
      return expected_size;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator; allowed only if RHS is empty (SegmentInfo has si-ptrs to objects held by *this)
    //! @param RHS the segment finder of interest
    SegmentFinder &SegmentFinder::operator=( const SegmentFinder &RHS)
    {
      if( this == &RHS)
      {
        return *this;
      }
      BCL_Assert
      (
        RHS.GetNumberSegments() == size_t( 0),
        "Cannot copy a non-empty SegmentFinder because it contains objects that hold pointers to internally-held data"
      );
      Reset();
      return *this;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SegmentFinder::Read( std::istream &ISTREAM)
    {
      BCL_Exit( "SegmentFinder cannot be read; it contains iterators", -1);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SegmentFinder::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_SequenceSegments, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ConditionToSequenceStatistics, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UnconditionalSequenceStatistics, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UnconditionalSequenceElementwiseStatistics, OSTREAM, INDENT);
      return OSTREAM;
    }

    //! @brief append an empty segment
    //! @param CONDITION current value of the condition
    void SegmentFinder::AddEmptySegment( const linal::VectorConstInterface< float> &CONDITION)
    {
      // finalize if necessary to end the current segment
      Finalize();

      // find the condition in the map; or insert it if it is not yet there
      m_CurrentConditionToSequenceStatisticsIterator = m_ConditionToSequenceStatistics.Find( CONDITION);
      if( m_CurrentConditionToSequenceStatisticsIterator == m_ConditionToSequenceStatistics.End())
      {
        // insert the new condition into the map construction
        m_CurrentConditionToSequenceStatisticsIterator =
          m_ConditionToSequenceStatistics.Insert
          (
            std::make_pair
            (
              linal::Vector< float>( CONDITION),
              storage::Triplet
              <
                storage::VectorND< 2, math::RunningAverageSD< linal::Vector< float> > >,
                storage::VectorND< 2, math::RunningAverageSD< float> >,
                size_t
              >
              (
                storage::VectorND< 2, math::RunningAverageSD< linal::Vector< float> > >(),
                storage::VectorND< 2, math::RunningAverageSD< float> >(),
                0
              )
            )
          ).first;
      }

      // Pushback a new segment with the correct information
      m_SequenceSegments.PushBack
      (
        SegmentInfo
        (
          !m_SequenceSegments.IsEmpty() ? m_SequenceSegments.LastElement().GetEndPosition() + 1 : 0,
          util::SiPtr< const size_t>( m_UnconditionalSequenceStatistics.Third()),
          util::SiPtr< const size_t>( m_CurrentConditionToSequenceStatisticsIterator->second.Third()),
          util::SiPtr< const linal::Vector< float> >( m_CurrentConditionToSequenceStatisticsIterator->first)
        )
      );

      // update counts of segments, conditional then unconditional
      ++m_CurrentConditionToSequenceStatisticsIterator->second.Third();
      ++m_UnconditionalSequenceStatistics.Third();
    }

  } // namespace descriptor
} // namespace bcl
