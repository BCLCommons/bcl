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
#include "descriptor/bcl_descriptor_segment_info.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief Statistic as string
    //! @param STATISTIC the statistic
    //! @return the Statistic as string
    const std::string &SegmentInfo::GetStatisticString( const Statistic &STATISTIC)
    {
      static const std::string s_names[ s_NumberStatistics + 1] =
      {
        "Condition",
        "Length",
        "DescriptorAve",
        "DescriptorSD",
        "DescriptorMin",
        "DescriptorMax",
        "DescriptorFirst",
        "DescriptorLast",
        "SegmentsBefore",
        "SegmentsAfter",
        "SegmentsTillNearestEnd",
        "ElementsFromSegmentStart",
        "ElementsFromSegmentEnd",
        "ElementsTillNearestSegmentEnd",
        "ConditionalSegmentsBefore",
        "ConditionalSegmentsAfter",
        "ConditionalSegmentsTillNearestEnd",
        GetStaticClassName< SegmentInfo>()
      };

      return s_names[ STATISTIC];
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor; not used but needed so that this object can appear in a list
    SegmentInfo::SegmentInfo() :
      m_StartPosition( 0),
      m_EndPosition( size_t( -1)),
      m_SegmentId( 0),
      m_ConditionalSegmentId( 0),
      m_DescriptorStatistics()
    {
    }

    //! @brief constructor from members for an empty segment
    //! @param START_POS index in the original sequence of the segment start
    //! @param SEGMENT_COUNT_REF pointer to a persistent count of segments in the sequence
    //! @param SEGMENT_CONDITION_COUNT_REF pointer to a persistent count of segments in the sequence with the same condition
    //! @param CONDITION_REF reference to the condition; not owned by this object
    SegmentInfo::SegmentInfo
    (
      const size_t &START_POS,
      const util::SiPtr< const size_t> &SEGMENT_COUNT_REF,
      const util::SiPtr< const size_t> &SEGMENT_CONDITION_COUNT_REF,
      const util::SiPtr< const linal::Vector< float> > &CONDITION_REF
    ) :
      m_StartPosition( START_POS),
      m_EndPosition( START_POS - 1),
      m_SegmentId( *SEGMENT_COUNT_REF),
      m_ConditionalSegmentId( *SEGMENT_CONDITION_COUNT_REF),
      m_DescriptorStatistics(),
      m_DescriptorLimits(),
      m_ConditionPtr( CONDITION_REF),
      m_ConditionNumSegmentsPtr( SEGMENT_CONDITION_COUNT_REF),
      m_OverallNumSegmentsPtr( SEGMENT_COUNT_REF)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SegmentInfo
    SegmentInfo *SegmentInfo::Clone() const
    {
      return new SegmentInfo( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SegmentInfo::GetClassIdentifier() const
    {
      return GetStaticClassName< SegmentInfo>();
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get the first id
    //! @return the first id
    size_t SegmentInfo::GetStartPosition() const
    {
      return m_StartPosition;
    }

    //! @brief get the last id
    //! @return the last id
    size_t SegmentInfo::GetEndPosition() const
    {
      return m_EndPosition;
    }

    //! @brief get the last id
    //! @return the last id
    size_t SegmentInfo::GetLength() const
    {
      return m_EndPosition - m_StartPosition + 1;
    }

    //! @brief get the id/index of this segment
    //! @return the id/index of this segment
    size_t SegmentInfo::GetSegmentId() const
    {
      return m_SegmentId;
    }

    //! @brief get the id/index of this segment relative to other segments with the same conditional value
    //! @return the id/index of this segment relative to other segments with the same conditional value
    size_t SegmentInfo::GetConditionalSegmentId() const
    {
      return m_ConditionalSegmentId;
    }

    //! @brief get the average of the averaged descriptor for this segment
    //! @return the average of the averaged descriptor for this segment
    const math::RunningAverageSD< linal::Vector< float> > &SegmentInfo::GetDescriptorAverageSD() const
    {
      return m_DescriptorStatistics;
    }

    //! @brief get the limits of the descriptor for this segment
    //! @return the limits of the descriptor for this segment
    const math::RunningMinMax< linal::Vector< float> > &SegmentInfo::GetDescriptorLimits() const
    {
      return m_DescriptorLimits;
    }

    //! @brief get the condition of the conditional descriptor for this segment
    //! @return he condition of the conditional descriptor for this segment
    const linal::VectorConstInterface< float> &SegmentInfo::GetCondition() const
    {
      return *m_ConditionPtr;
    }

    //! @brief copy a statistic into a vector
    //! @param STATISTIC the statistic desired
    //! @param STORAGE storage for the statistic
    //! @param POSITION current position in the original sequence
    void SegmentInfo::CopyStatistic
    (
      linal::VectorInterface< float> &REFERENCE,
      const Statistic &STATISTIC,
      const size_t &POSITION
    ) const
    {
      switch( STATISTIC)
      {
        case e_Condition:
          std::copy( GetCondition().Begin(), GetCondition().End(), REFERENCE.Begin());
          break;
        case e_Length:
          REFERENCE( 0) = GetLength();
          break;
        case e_DescriptorAve:
          std::copy
          (
            m_DescriptorStatistics.GetAverage().Begin(),
            m_DescriptorStatistics.GetAverage().End(),
            REFERENCE.Begin()
          );
          break;
        case e_DescriptorSD:
          std::copy
          (
            m_DescriptorStatistics.GetStandardDeviation().Begin(),
            m_DescriptorStatistics.GetStandardDeviation().End(),
            REFERENCE.Begin()
          );
          break;
        case e_DescriptorMin:
          std::copy
          (
            m_DescriptorLimits.GetMin().Begin(),
            m_DescriptorLimits.GetMin().End(),
            REFERENCE.Begin()
          );
          break;
        case e_DescriptorMax:
          std::copy
          (
            m_DescriptorLimits.GetMax().Begin(),
            m_DescriptorLimits.GetMax().End(),
            REFERENCE.Begin()
          );
          break;
        case e_DescriptorFirst:
          std::copy
          (
            m_DescriptorFirst.Begin(),
            m_DescriptorFirst.End(),
            REFERENCE.Begin()
          );
          break;
        case e_DescriptorLast:
          std::copy
          (
            m_DescriptorLast.Begin(),
            m_DescriptorLast.End(),
            REFERENCE.Begin()
          );
          break;
        case e_SegmentsBefore:
          REFERENCE( 0) = m_SegmentId;
          break;
        case e_SegmentsAfter:
          REFERENCE( 0) = *m_OverallNumSegmentsPtr - m_SegmentId - 1;
          break;
        case e_SegmentsTillNearestEnd:
          REFERENCE( 0) = float( *m_OverallNumSegmentsPtr) - float( m_SegmentId) - 1.0;
          REFERENCE( 0) = std::min( REFERENCE( 0), float( m_SegmentId));
          break;
        case e_ElementsFromSegmentStart:
          REFERENCE( 0) = float( POSITION) - float( m_StartPosition);
          break;
        case e_ElementsFromSegmentEnd:
          REFERENCE( 0) = float( m_EndPosition) - float( POSITION);
          break;
        case e_ElementsTillNearestSegmentEnd:
          REFERENCE( 0) = math::Absolute( float( m_EndPosition) - float( POSITION));
          REFERENCE( 0) = std::min( REFERENCE( 0), float( math::Absolute( float( POSITION) - float( m_StartPosition))));
          break;
        case e_ConditionalSegmentsBefore:
          REFERENCE( 0) = m_ConditionalSegmentId;
          break;
        case e_ConditionalSegmentsAfter:
          REFERENCE( 0) = *m_ConditionNumSegmentsPtr - m_ConditionalSegmentId - 1;
          break;
        case e_ConditionalSegmentsTillNearestEnd:
          REFERENCE( 0) =
            std::min( m_ConditionalSegmentId, size_t( *m_ConditionNumSegmentsPtr - m_ConditionalSegmentId - 1));
          break;
        case s_NumberStatistics:
        default:
          BCL_Exit( "Invalid statistic: " + GetStatisticString( STATISTIC), -1);
          break;
      }
    }

    //! @brief get a statistic as a vector
    //! @param STATISTIC the statistic desired
    //! @param POSITION current position in the original sequence
    linal::Vector< float> SegmentInfo::GetStatistic( const Statistic &STATISTIC, const size_t &POSITION) const
    {
      linal::Vector< float> value( GetStatisticSize( STATISTIC), float( 0));
      CopyStatistic( value, STATISTIC, POSITION);
      return value;
    }

    //! @brief get the size of the statistic
    //! @param STATISTIC the actual statistic of interest
    //! @return the size of the statistic
    size_t SegmentInfo::GetStatisticSize( const Statistic &STATISTIC) const
    {
      size_t expected_size( 1);
      if( STATISTIC == e_Condition)
      {
        expected_size = GetCondition().GetSize();
      }
      else if
      (
        STATISTIC == e_DescriptorAve || STATISTIC == e_DescriptorSD
        || STATISTIC == e_DescriptorMin || STATISTIC == e_DescriptorMax
        || STATISTIC == e_DescriptorFirst || STATISTIC == e_DescriptorLast
      )
      {
        expected_size = m_DescriptorStatistics.GetAverage().GetSize();
      }
      return expected_size;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add another element to the segment info
    //! @param DESCRIPTOR the actual descriptor value
    void SegmentInfo::PushBack( const linal::VectorConstInterface< float> &DESCRIPTOR)
    {
      ++m_EndPosition;
      m_DescriptorStatistics += DESCRIPTOR;
      m_DescriptorLimits += DESCRIPTOR;
      if( m_EndPosition == m_StartPosition)
      {
        m_DescriptorFirst = DESCRIPTOR;
      }
      m_DescriptorLast = DESCRIPTOR;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SegmentInfo::Read( std::istream &ISTREAM)
    {
      BCL_Exit( "SegmentInfo contains iterators and, as such, cannot be copied", -1);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SegmentInfo::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_StartPosition, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_EndPosition, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SegmentId, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ConditionalSegmentId, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DescriptorStatistics, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DescriptorLimits, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ConditionPtr, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ConditionNumSegmentsPtr, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_OverallNumSegmentsPtr, OSTREAM, INDENT);
      return OSTREAM;
    }

  } // namespace descriptor
} // namespace bcl
