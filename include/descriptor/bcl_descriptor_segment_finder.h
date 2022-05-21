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

#ifndef BCL_DESCRIPTOR_SEGMENT_FINDER_H_
#define BCL_DESCRIPTOR_SEGMENT_FINDER_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_segment_info.h"
#include "linal/bcl_linal_vector.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_running_average_sd.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_triplet.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SegmentFinder
    //! @brief Computes information about the current and nearby segments
    //!        A segment is a region of consecutive elements of a sequence which have the same values returned by a user
    //!        defined conditional descriptor.  The user can also select a descriptor to take mean/sd for over this
    //!        segment if they desire. Information about adjacent segments, as well as overall sequence statistics may
    //!        also be obtained. This class is similar to WindowConditional average but differs in that it looks at
    //!        segments, which are not relative to any particular amino acid, rather than length-restricted windows.
    //!
    //! @tparam t_DataType type of the argument for which the description is going to be generated
    //!
    //! @see @link example_descriptor_segment_finder.cpp @endlink
    //! @author mendenjl
    //! @date Mar 05, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SegmentFinder :
      public util::ObjectInterface
    {
    public:

      //! These statistics are for the whole sequence
      enum Statistic
      {
        e_NumSegments,                      //!< number of segments in the sequence
        e_FractionSegments,                 //!< fraction of segments with a condition value (1 if no condition is given)
        e_FractionElements,                 //!< fraction of elements with a condition value (1 if no condition is given)
        e_LengthSegmentsAve,                //!< average length of segments in the sequence
        e_LengthSegmentsSD,                 //!< standard deviation of length of segments in the sequence
        e_LengthSegmentsElementwiseAve,     //!< average length of segments in the sequence, weighted by # residues in segment
        e_LengthSegmentsElementwiseSD,      //!< standard deviation of length of segments in the sequence, weighted by # residues in segment
        e_DescriptorSegmentsAve,            //!< average of descriptor over all segments
        e_DescriptorSegmentsSD,             //!< standard deviation of descriptor over all segments
        e_DescriptorSegmentsElementwiseAve, //!< average of descriptor weighted by # residues in segment
        e_DescriptorSegmentsElementwiseSD,  //!< standard deviation of descriptor weighted by # residues in segment
        e_ConditionSegmentsAve,             //!< average of condition value over all segments
        e_ConditionSegmentsSD,              //!< standard deviation of condition value over all segments
        e_ConditionSegmentsElementwiseAve,  //!< average of condition value over all segments, weighted by # residues in segment
        e_ConditionSegmentsElementwiseSD,   //!< standard deviation of condition value over all segments, weighted by # residues in segment
        s_NumberStatistics
      };

      //! @brief Statistic as string
      //! @param STATISTIC the statistic
      //! @return the Statistic as string
      static const std::string &GetStatisticString( const Statistic &STATISTIC);

      //! simplifies the usage of the UnconditionalStatistic enum of this class
      typedef util::WrapperEnum< Statistic, &GetStatisticString, s_NumberStatistics> StatisticEnum;

      //! @brief test whether the statistic returns a value for every value in the condition
      //! @param STATISTIC the statistic of interest
      //! @return true if the statistic returns a value for every value in the condition
      static bool GetIsStatisticOfCondition( const Statistic &STATISTIC);

      //! @brief test whether the statistic returns a value for every value in the statistic descriptor
      //! @param STATISTIC the statistic of interest
      //! @return true if the statistic returns a value for every value in the statistic descriptor
      static bool GetIsStatisticOfDescriptor( const Statistic &STATISTIC);

    private:

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class VectorLessThanRespectUndefined
      //! @brief Class used internally to compare vectors properly, even if they have undefined values
      //!        undefined values only compare == to each other in this class
      //!        C++ ISO standard dictates that NaNs normally compare == to nothing; but descriptors sometimes use NaNs
      //!        to indicate that the descriptor could not be calculated
      //! @remarks example unnecessary
      //! @author mendenjl
      //! @date Jun 19, 2014
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      struct VectorLessThanRespectUndefined
      {
        //! @brief function that compares two linal::Vectors, ignoring effectively considering undefined values to be FLT_MAX
        //! @param A, B the vectors to compare
        //! @return A < B, respective of undefined values
        //! This is necessary because maps require that their keys have weak ordering; which is not possible if the values
        //! are undefined
        bool operator()( const linal::Vector< float> &A, const linal::Vector< float> &B) const;

        //! @brief simpler test for ==
        //! @param A, B the vectors to compare
        //! @return A == B, respective of undefined values
        static bool Equal( const linal::Vector< float> &A, const linal::Vector< float> &B);
      };

    //////////
    // data //
    //////////

      //! list of segment info for the current sequence
      storage::List< SegmentInfo> m_SequenceSegments;

      //! Statistics for every condition
      //! Key is the conditional value
      //! Value.First().First() -> Stats for descriptor, all segments weighted equally
      //! Value.First().Second() -> Stats for descriptor, segments weighted by length
      //! Value.Second().First() -> Stats for length, all segments weighted equally
      //! Value.Second().Second() -> Stats for length, segments weighted by length
      //! Value.Third() -> Number of segments with the given condition
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
      > m_ConditionToSequenceStatistics;

      //! unconditional sequence statistics, all segments considered equally
      //! First().First() is stats for the descriptor
      //! First().Second() is stats for the condition
      //! Second() is stats for the length
      //! Third() is number of segments in the sequence
      storage::Triplet
      <
        storage::VectorND< 2, math::RunningAverageSD< linal::Vector< float> > >,
        math::RunningAverageSD< float>,
        size_t
      > m_UnconditionalSequenceStatistics;

      //! unconditional sequence statistics, weighted by segment size
      //! First().First() is stats for the descriptor
      //! First().Second() is stats for the condition
      //! Second() is stats for the length
      storage::Pair
      <
        storage::VectorND< 2, math::RunningAverageSD< linal::Vector< float> > >,
        math::RunningAverageSD< float>
      > m_UnconditionalSequenceElementwiseStatistics;

      //! vector; index is iterator position; value is iterator to the segment info for that position
      std::vector< storage::List< SegmentInfo>::const_iterator> m_IteratorsForSequence;

      //! iterator on the statistics for the current condition
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
      >::iterator m_CurrentConditionToSequenceStatisticsIterator;

    public:

      typedef storage::List< SegmentInfo>::const_iterator const_iterator;

      //! instances of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SegmentFinder();

      //! @brief copy constructor
      SegmentFinder( const SegmentFinder &PARENT);

      //! @brief Clone function
      //! @return pointer to new SegmentFinder
      SegmentFinder *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get an iterator to the first segment
      //! @return an iterator to the first segment
      const_iterator SegmentBegin() const;

      //! @brief get an iterator beyond the last segment
      //! @return an iterator beyond the last segment
      const_iterator SegmentEnd() const;

      //! @brief get an iterator for the element at a particular position in the original sequence
      //! @param POSITION position in the original sequence to get an iterator for
      //! @return an iterator for the element at a particular position in the original sequence
      const_iterator SegmentForPosition( const size_t &POSITION) const;

      //! @brief test whether there all segments' statistics have been accumulated
      //! @return true if the sequence has been finalized
      bool IsFinalized() const;

      //! @brief get the number of segments found so far
      //! @return the number of segments found so far
      size_t GetNumberSegments() const;

      //! @brief test whether the segment finder is empty
      //! @return true if the segment finder is empty
      bool IsEmpty() const;

      //! @brief copy a statistic into a vector
      //! @param STORAGE storage for the statistic
      //! @param STATISTIC the statistic desired
      void CopyStatisticConsideringAllSegments
      (
        linal::VectorInterface< float> &REFERENCE,
        const Statistic &STATISTIC
      ) const;

      //! @brief copy a statistic into a vector
      //! @param STORAGE storage for the statistic
      //! @param STATISTIC the statistic desired
      //! @param CONDITION the condition of interest
      void CopyStatisticConsideringSegmentsWithCondition
      (
        linal::VectorInterface< float> &REFERENCE,
        const Statistic &STATISTIC,
        const linal::VectorConstInterface< float> &CONDITION
      ) const;

      //! @brief get a statistic as a vector
      //! @param STATISTIC the statistic desired
      linal::Vector< float> GetStatisticConsideringAllSegments( const Statistic &STATISTIC) const;

      //! @brief get a statistic as a vector
      //! @param STATISTIC the statistic desired
      //! @param CONDITION the condition of interest
      linal::Vector< float> GetStatisticConsideringSegmentsWithCondition
      (
        const Statistic &STATISTIC,
        const linal::VectorConstInterface< float> &CONDITION
      ) const;

      //! @brief get the size of the statistic
      //! @param STATISTIC the actual statistic of interest
      //! @return the size of the statistic
      size_t GetStatisticSize( const Statistic &STATISTIC) const;

    ////////////////
    // operations //
    ////////////////

      //! @brief add information for the next element in the sequence to this class
      //! @param CONDITION current value of the condition
      //! @param DESCRIPTOR current value of the descriptor
      void ConsiderNextElementOfSequence
      (
        const linal::VectorConstInterface< float> &CONDITION,
        const linal::VectorConstInterface< float> &DESCRIPTOR
      );

      //! @brief finalize; flushes most recent segment's statistics into the various maps
      void Finalize();

      //! @brief reset -- resets this class so that it forgets about all the segments
      void Reset();

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment operator; allowed only if RHS is empty (SegmentInfo has si-ptrs to objects held by *this)
      //! @param RHS the segment finder of interest
      SegmentFinder &operator=( const SegmentFinder &RHS);

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

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief append an empty segment
      //! @param CONDITION current value of the condition
      void AddEmptySegment( const linal::VectorConstInterface< float> &CONDITION);

    }; // class SegmentFinder

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_SEGMENT_FINDER_H_
