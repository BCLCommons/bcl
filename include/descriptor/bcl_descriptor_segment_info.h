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

#ifndef BCL_DESCRIPTOR_SEGMENT_INFO_H_
#define BCL_DESCRIPTOR_SEGMENT_INFO_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_running_min_max.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SegmentInfo
    //! @brief contains information about each segment in the sequence
    //!
    //! @see @link example_descriptor_segment_info.cpp @endlink
    //! @author mendenjl
    //! @date Mar 05, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SegmentInfo :
      public util::ObjectInterface
    {
    public:

      enum Statistic
      {
        e_Condition,                      //!< Value of the condition for the segment
        e_Length,                         //!< Length of the segment
        e_DescriptorAve,                  //!< Average value of the descriptor for this segment
        e_DescriptorSD,                   //!< Standard deviation of the descriptor for this segment
        e_DescriptorMin,                  //!< Minimum value of the descriptor for this segment
        e_DescriptorMax,                  //!< Maximum value of the descriptor for this segment
        e_DescriptorFirst,                //!< Value of the descriptor for the first element in the segment
        e_DescriptorLast,                 //!< Value of the descriptor for the last element in the segment
        e_SegmentsBefore,                 //!< number of segments before this one in the sequence
        e_SegmentsAfter,                  //!< number of segments after this one in the sequence
        e_SegmentsTillNearestEnd,         //!< number of segments before the nearest of the two ends
        e_ElementsFromSegmentStart,       //!< number of elements before this one in the segment
        e_ElementsFromSegmentEnd,         //!< number of elements after this one in the segment
        e_ElementsTillNearestSegmentEnd,  //!< number of elements before the nearest of the two ends
        e_ConditionalSegmentsBefore,      //!< number of segments before this one with the same condition
        e_ConditionalSegmentsAfter,       //!< number of segments after this one with the same condition
        e_ConditionalSegmentsTillNearestEnd, //!< number of segments till the nearest end with the same condition
        s_NumberStatistics
      };

      //! @brief Statistic as string
      //! @param STATISTIC the statistic
      //! @return the Statistic as string
      static const std::string &GetStatisticString( const Statistic &STATISTIC);

      //! simplifies the usage of the UnconditionalStatistic enum of this class
      typedef util::WrapperEnum< Statistic, &GetStatisticString, s_NumberStatistics> StatisticEnum;

    private:
    //////////
    // data //
    //////////

      size_t m_StartPosition;          //!< Position that starts this segment
      size_t m_EndPosition;            //!< Final position for this segment
      size_t m_SegmentId;              //!< Index of the segment overall
      size_t m_ConditionalSegmentId;   //!< Index of the condition in the segment
      math::RunningAverageSD< linal::Vector< float> > m_DescriptorStatistics; //!< descriptor stats over this segment
      math::RunningMinMax< linal::Vector< float> >    m_DescriptorLimits;     //!< descriptor segmental min and max
      linal::Vector< float>                           m_DescriptorFirst; //!< descriptor value for the first element
      linal::Vector< float>                           m_DescriptorLast; //!< descriptor value for the last element

      //! Pointer to condition
      util::SiPtr< const linal::Vector< float> > m_ConditionPtr;

      //! Pointer to size_t representing # of segments with the same condition
      util::SiPtr< const size_t> m_ConditionNumSegmentsPtr;

      //! Pointer to size_t representing # of segments in the sequence
      util::SiPtr< const size_t> m_OverallNumSegmentsPtr;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor; not used but needed so that this object can appear in a list
      SegmentInfo();

      //! @brief constructor from members for an empty segment
      //! @param START_POS index in the original sequence of the segment start
      //! @param SEGMENT_COUNT_REF pointer to a persistent count of segments in the sequence
      //! @param SEGMENT_CONDITION_COUNT_REF pointer to a persistent count of segments in the sequence with the same condition
      //! @param CONDITION_REF reference to the condition; not owned by this object
      SegmentInfo
      (
        const size_t &START_POS,
        const util::SiPtr< const size_t> &SEGMENT_COUNT_REF,
        const util::SiPtr< const size_t> &SEGMENT_CONDITION_COUNT_REF,
        const util::SiPtr< const linal::Vector< float> > &CONDITION_REF
      );

      //! @brief Clone function
      //! @return pointer to new SegmentInfo
      SegmentInfo *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get the first id
      //! @return the first id
      size_t GetStartPosition() const;

      //! @brief get the last id
      //! @return the last id
      size_t GetEndPosition() const;

      //! @brief get the last id
      //! @return the last id
      size_t GetLength() const;

      //! @brief get the id/index of this segment
      //! @return the id/index of this segment
      size_t GetSegmentId() const;

      //! @brief get the id/index of this segment relative to other segments with the same conditional value
      //! @return the id/index of this segment relative to other segments with the same conditional value
      size_t GetConditionalSegmentId() const;

      //! @brief get the average of the averaged descriptor for this segment
      //! @return the average of the averaged descriptor for this segment
      const math::RunningAverageSD< linal::Vector< float> > &GetDescriptorAverageSD() const;

      //! @brief get the limits of the descriptor for this segment
      //! @return the limits of the descriptor for this segment
      const math::RunningMinMax< linal::Vector< float> > &GetDescriptorLimits() const;

      //! @brief get the condition of the conditional descriptor for this segment
      //! @return he condition of the conditional descriptor for this segment
      const linal::VectorConstInterface< float> &GetCondition() const;

      //! @brief copy a statistic into a vector
      //! @param STORAGE storage for the statistic
      //! @param STATISTIC the statistic desired
      //! @param POSITION current position in the original sequence
      void CopyStatistic
      (
        linal::VectorInterface< float> &REFERENCE,
        const Statistic &STATISTIC,
        const size_t &POSITION = size_t( 0)
      ) const;

      //! @brief get a statistic as a vector
      //! @param STATISTIC the statistic desired
      //! @param POSITION current position in the original sequence
      linal::Vector< float> GetStatistic
      (
        const Statistic &STATISTIC,
        const size_t &POSITION = size_t( 0)
      ) const;

      //! @brief get the size of the statistic
      //! @param STATISTIC the actual statistic of interest
      //! @return the size of the statistic
      size_t GetStatisticSize( const Statistic &STATISTIC) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief add another element to the segment info
      //! @param DESCRIPTOR the actual descriptor value
      void PushBack( const linal::VectorConstInterface< float> &DESCRIPTOR);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class SegmentInfo

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_SEGMENT_INFO_H_
