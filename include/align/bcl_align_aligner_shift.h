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

#ifndef BCL_ALIGN_ALIGNER_SHIFT_H_
#define BCL_ALIGN_ALIGNER_SHIFT_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_alignment_interface.h"
#include "bcl_align_alignment_leaf.h"
#include "bcl_align_pairwise_aligner_interface.h"
#include "bcl_align_sequence.h"
#include "function/bcl_function_unary_interface.h"
#include "math/bcl_math_comparisons.h"
#include "util/bcl_util_binary_function_interface.h"
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignerShift
    //! @brief aligns two sequences by systematically creating shifted alignments
    //! @details The alignments are created by systematically trying all relative shifted aligments, e.g.
    //! seq1: ABCD
    //! seq2: abcd
    //! align1: ---ABCD
    //!         abcd---
    //! align2: --ABCD
    //!         abcd--
    //! .....
    //! align7: ABCD---
    //!         ---abcd
    //! and the alignment with the highest score is returned
    //!
    //! @tparam t_Member type of members in assignments within the alignments
    //!
    //! @see @link example_align_aligner_shift.cpp @endlink
    //! @author woetzen
    //! @date Mar 14, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class AlignerShift :
      public PairwiseAlignerInterface< t_Member>
    {

    private:

    //////////
    // data //
    //////////

      //! minimum number of assignment without gap
      size_t m_MinNumberAssignmentsWithoutGap;

      //! maximum number of gaps
      size_t m_MaxNumberGapsLeft;
      size_t m_MaxNumberGapsRight;

      //! score a complete alignment
      util::ShPtr< function::UnaryInterface< const AlignmentInterface< t_Member>, double> > m_AlignmentScore;

      //! comparison function for score
      util::SiPtr< const util::BinaryFunctionInterface< double, double, bool> > m_ScoreComparison;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AlignerShift() :
        m_MinNumberAssignmentsWithoutGap( 1),
        m_MaxNumberGapsLeft( util::GetUndefined< int>()),
        m_MaxNumberGapsRight( util::GetUndefined< int>()),
        m_AlignmentScore(),
        m_ScoreComparison( &( **math::Comparisons< double>::GetEnums().e_Less))
      {
      }

      //! @brief Clone function
      //! @return pointer to new AlignerShift< t_Member>
      AlignerShift< t_Member> *Clone() const
      {
        return new AlignerShift< t_Member>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief set the minimum number of assignments if non gap members
      //! @param MIN_NUMBER_ASSIGNMENTS_WITHOUT_GAP number of assignments within alignment without gap - or the minimum number of assigned members
      void SetMinNumberAssignmentsWithoutGap( const size_t MIN_NUMBER_ASSIGNMENTS_WITHOUT_GAP)
      {
        m_MinNumberAssignmentsWithoutGap = MIN_NUMBER_ASSIGNMENTS_WITHOUT_GAP;
      }

      //! @brief set the maximum number of gaps
      //! @param MAX_NUMBER_GAPS_LEFT max number of gaps within alignment on left side
      //! @param MAX_NUMBER_GAPS_RIGHT max number of gaps within alignment on right side
      void SetMaxNumberGaps( const size_t MAX_NUMBER_GAPS_LEFT, const size_t MAX_NUMBER_GAPS_RIGHT)
      {
        m_MaxNumberGapsLeft = MAX_NUMBER_GAPS_LEFT;
        m_MaxNumberGapsRight = MAX_NUMBER_GAPS_RIGHT;
      }

      //! @brief SetScoringFunction sets the score::Assignment scoring function
      //! @param SCORE is the score::Assignment to be used for scoring an assignment
      void SetScoringFunction( const score::AssignmentWithGap< t_Member> &SCORE)
      {

      }

      //! @brief set scoring function for complete alignment
      //! @param SP_SCORE the scoring function for an entire elignment
      void SetScoringFunction( const util::ShPtr< function::UnaryInterface< const AlignmentInterface< t_Member>, double> > &SP_SCORE)
      {
        m_AlignmentScore = SP_SCORE;
      }

      //! @brief set the comparison function for a better score
      //! @param SP_COMPARISON the comparison function to determine if a score is better than another
      void SetScoreComparisonFunction( const util::SiPtr< const util::BinaryFunctionInterface< double, double, bool> > &SP_COMPARISON)
      {
        m_ScoreComparison = SP_COMPARISON;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Align is the function which does the actually aligns data together
      //! @param ALIGNMENT_TOP the first alignment to be aligned (no const, it must go in AlignmentNode's m_ChildAlignments)
      //! @param ALIGNMENT_BOTTOM the second alignment to be aligned
      //! @return returns a pair of the Alignment and a double which is the score
      storage::Pair< AlignmentNode< t_Member>, double> AlignPair
      (
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_TOP,
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_BOTTOM
      ) const
      {
        storage::Pair< AlignmentNode< t_Member>, double> best_alignment;
        best_alignment.Second() = util::GetUndefined< double>();

        const std::set< int> shifts
        (
          PossibleShifts
          (
            *ALIGNMENT_TOP, *ALIGNMENT_BOTTOM,
            m_MinNumberAssignmentsWithoutGap, m_MaxNumberGapsLeft, m_MaxNumberGapsRight
          )
        );

        // iterate through all possible shifts
        for( std::set< int>::const_iterator itr( shifts.begin()), itr_end( shifts.end()); itr != itr_end; ++itr)
        {
          storage::Pair< AlignmentNode< t_Member>, double> current_alignment
          (
            AlignmentFromShift( ALIGNMENT_TOP, ALIGNMENT_BOTTOM, *itr),
            double( 0.0)
          );

          // score the current alignment and update best
          current_alignment.Second() = m_AlignmentScore->operator()( current_alignment.First());
          if( !util::IsDefined( best_alignment.Second()))
          {
            best_alignment = current_alignment;
          }
          else if
          (
               util::IsDefined( current_alignment.Second())
            && m_ScoreComparison->operator()( current_alignment.Second(), best_alignment.Second())
          )
          {
            best_alignment = current_alignment;
          }
        }

        // end
        return best_alignment;
      }

      //! @brief generate all possible alignments for two Pairs of sequences
      //! @param ALIGNMENT_LEFT_TOP left top alignment
      //! @param ALIGNMENT_LEFT_BOTTOM left bottom alignment aligned to left top to make up the left part
      //! @param ALIGNMENT_RIGHT_TOP right top alignment
      //! @param ALIGNMENT_RIGHT_BOTTOM right bottom alignment aligned to right top to make up the right part
      //! @param MAX_NUMBER_GAPS_RIGHT number of gaps right to left alignment (first part) at fusion point
      //! @param MAX_NUMBER_GAPS_LEFT number of gaps left to right alignment (second part) at fusion point
      //! @param INSERT_SEPARATOR insert a separator, i.e. an assignment of gaps
      //! @return all possible fuses alignments between left and right
      storage::Vector< AlignmentNode< t_Member> > AlignmentsPairPair
      (
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_LEFT_TOP,
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_LEFT_BOTTOM,
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_RIGHT_TOP,
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_RIGHT_BOTTOM,
        const size_t MAX_NUMBER_GAPS_RIGHT, const size_t MAX_NUMBER_GAPS_LEFT,
        const bool INSERT_SEPARATOR
      ) const
      {
        // left
        const std::set< int> shifts_left
        (
          PossibleShifts
          (
            *ALIGNMENT_LEFT_TOP, *ALIGNMENT_LEFT_BOTTOM,
            m_MinNumberAssignmentsWithoutGap, m_MaxNumberGapsLeft,
            MAX_NUMBER_GAPS_RIGHT
          )
        );
        // right
        const std::set< int> shifts_right
        (
          PossibleShifts
          (
            *ALIGNMENT_RIGHT_TOP, *ALIGNMENT_RIGHT_BOTTOM,
            m_MinNumberAssignmentsWithoutGap, MAX_NUMBER_GAPS_LEFT,
            m_MaxNumberGapsRight
          )
        );

        storage::Vector< AlignmentNode< t_Member> > alignments;
        alignments.AllocateMemory( shifts_left.size() * shifts_right.size());

        // iterate through all possible shifts
        for( std::set< int>::const_iterator itr1( shifts_left.begin()), itr1_end( shifts_left.end()); itr1 != itr1_end; ++itr1)
        {
          // create alignment for left sequences
          AlignmentNode< t_Member> alignment_left
          (
            AlignmentFromShift( ALIGNMENT_LEFT_TOP, ALIGNMENT_LEFT_BOTTOM, *itr1)
          );
          for( std::set< int>::const_iterator itr2( shifts_right.begin()), itr2_end( shifts_right.end()); itr2 != itr2_end; ++itr2)
          {
            // create alignment for right seuqences
            AlignmentNode< t_Member> alignment_right
            (
              AlignmentFromShift( ALIGNMENT_RIGHT_TOP, ALIGNMENT_RIGHT_BOTTOM, *itr2)
            );

            // sequence for the top and bottom
            util::ShPtr< Sequence< t_Member> > sequence_top( new Sequence< t_Member>( *ALIGNMENT_LEFT_TOP->GetSequences().FirstElement(), ALIGNMENT_LEFT_TOP->GetSequences().FirstElement()->GetSequenceId()));
            util::ShPtr< Sequence< t_Member> > sequence_bottom( new Sequence< t_Member>( *ALIGNMENT_LEFT_BOTTOM->GetSequences().FirstElement(), ALIGNMENT_LEFT_BOTTOM->GetSequences().FirstElement()->GetSequenceId()));
            sequence_top->Append( *ALIGNMENT_RIGHT_TOP->GetSequences().FirstElement());
            sequence_bottom->Append( *ALIGNMENT_RIGHT_BOTTOM->GetSequences().FirstElement());

            // alignment for top and bottom
            util::ShPtr< AlignmentInterface< t_Member> > align_leaf_a( new AlignmentLeaf< t_Member>( sequence_top));
            util::ShPtr< AlignmentInterface< t_Member> > align_leaf_b( new AlignmentLeaf< t_Member>( sequence_bottom));

            // generated alignment from sequences
            AlignmentNode< t_Member> current_alignment( align_leaf_a, align_leaf_b);
            // fuse alignments
            current_alignment.Append( alignment_left.GetAssignments());

            // separator so that left and right side can be identified easier
            if( INSERT_SEPARATOR)
            {
              current_alignment.Append
              (
                util::ShPtr< Assignment< t_Member> >
                (
                  new Assignment< t_Member>( util::SiPtr< const t_Member>(), util::SiPtr< const t_Member>())
                )
              );
            }
            current_alignment.Append( alignment_right.GetAssignments());

            // insert alignment to all possible alignments
            alignments.PushBack( current_alignment);
          }
        }

        return alignments;
      }

      //! @brief align pair of pairs
      //! @param ALIGNMENT_LEFT_TOP left top alignment
      //! @param ALIGNMENT_LEFT_BOTTOM left bottom alignment aligned to left top to make up the left part
      //! @param ALIGNMENT_RIGHT_TOP right top alignment
      //! @param ALIGNMENT_RIGHT_BOTTOM right bottom alignment aligned to right top to make up the right part
      //! @param MAX_NUMBER_GAPS_RIGHT number of gaps right to left alignment (first part) at fusion point
      //! @param MAX_NUMBER_GAPS_LEFT number of gaps left to right alignment (second part) at fusion point
      //! @param INSERT_SEPARATOR insert a separator, i.e. an assignment of gaps
      //! @return best fuse alignment between left and right and its score
      storage::Pair< AlignmentNode< t_Member>, double> AlignPairPair
      (
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_LEFT_TOP,
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_LEFT_BOTTOM,
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_RIGHT_TOP,
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_RIGHT_BOTTOM,
        const size_t MAX_NUMBER_GAPS_RIGHT, const size_t MAX_NUMBER_GAPS_LEFT,
        const bool INSERT_SEPARATOR
      ) const
      {
        // generate all possible alignment
        const storage::Vector< AlignmentNode< t_Member> > all_alignments
        (
          AlignmentsPairPair
          (
            ALIGNMENT_LEFT_TOP,
            ALIGNMENT_LEFT_BOTTOM,
            ALIGNMENT_RIGHT_TOP,
            ALIGNMENT_RIGHT_BOTTOM,
            MAX_NUMBER_GAPS_RIGHT, MAX_NUMBER_GAPS_LEFT,
            INSERT_SEPARATOR
          )
        );

        storage::Pair< AlignmentNode< t_Member>, double> best_alignment;
        best_alignment.Second() = util::GetUndefined< double>();

        // iterate through all possible alignments and determine the best scoring
        for
        (
          typename storage::Vector< AlignmentNode< t_Member> >::const_iterator
            itr( all_alignments.Begin()), itr_end( all_alignments.End());
          itr != itr_end;
          ++itr
        )
        {
          // score the current alignment and update best
          const double current_score( m_AlignmentScore->operator()( *itr));

          // first alignment with defined score is the best alignment
          if( !util::IsDefined( best_alignment.Second()))
          {
            best_alignment.First() = *itr;
            best_alignment.Second() = current_score;
          }
          // current alignment scores better than the best alignment -> update
          else if( util::IsDefined( current_score) && m_ScoreComparison->operator()( current_score, best_alignment.Second()))
          {
            best_alignment.First() = *itr;
            best_alignment.Second() = current_score;
          }
        }

        return best_alignment;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_MinNumberAssignmentsWithoutGap, ISTREAM);
        io::Serialize::Read( m_MaxNumberGapsLeft             , ISTREAM);
        io::Serialize::Read( m_MaxNumberGapsRight            , ISTREAM);
        io::Serialize::Read( m_AlignmentScore                , ISTREAM);

        // return the stream
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_MinNumberAssignmentsWithoutGap, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_MaxNumberGapsLeft             , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_MaxNumberGapsRight            , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_AlignmentScore                , OSTREAM, INDENT);

        // return the stream
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief calculate the possible shift for given alignments
      //! @param ALIGNMENT_TOP alignment on top
      //! @param ALIGNMENT_BOTTOM alignment on bottom
      //! @param MIN_NUMBER_ASSIGNMENTS_WITHOUT_GAP minimal number of assignments that any shift should have
      //! @param MAX_NUMBER_GAPS_LEFT maximal number of gaps on the left, that should be left
      //! @param MAX_NUMBER_GAPS_RIGHT maximal number of gaps that should be left on the right
      //! @return set of shifts as integers
      static std::set< int> PossibleShifts
      (
        const AlignmentInterface< t_Member> &ALIGNMENT_TOP,
        const AlignmentInterface< t_Member> &ALIGNMENT_BOTTOM,
        const size_t MIN_NUMBER_ASSIGNMENTS_WITHOUT_GAP,
        const size_t MAX_NUMBER_GAPS_LEFT,
        const size_t MAX_NUMBER_GAPS_RIGHT
      )
      {
        std::set< int> possible_shifts;

        const size_t size_top( ALIGNMENT_TOP.GetSize());
        const size_t size_bottom( ALIGNMENT_BOTTOM.GetSize());
        const size_t smallest_size( std::min( size_top, size_bottom));

        // if either size is too small, return
        if( smallest_size < MIN_NUMBER_ASSIGNMENTS_WITHOUT_GAP)
        {
          return possible_shifts;
        }

        // gaps left and right
        const int size_difference( int( size_top) - int( size_bottom));
        const size_t max_gaps_top( size_bottom - MIN_NUMBER_ASSIGNMENTS_WITHOUT_GAP);
        const size_t max_gaps_bottom( size_top - MIN_NUMBER_ASSIGNMENTS_WITHOUT_GAP);

        const int max_gaps_top_right( int( std::min( max_gaps_top, MAX_NUMBER_GAPS_RIGHT)));
        const int max_gaps_bottom_right( int( std::min( max_gaps_bottom, MAX_NUMBER_GAPS_RIGHT)));

        const int max_gaps_top_left( int( std::min( max_gaps_top, MAX_NUMBER_GAPS_LEFT)));
        const int max_gaps_bottom_left( int( std::min( max_gaps_bottom, MAX_NUMBER_GAPS_LEFT)));

        // left
        std::set< int> possible_shifts_left;
        for( int shift( -max_gaps_bottom_left); shift <= max_gaps_top_left; ++shift)
        {
          possible_shifts_left.insert( shift);
        }
        // right
        std::set< int> possible_shifts_right;
        for( int shift( -max_gaps_top_right); shift <= max_gaps_bottom_right; ++shift)
        {
          possible_shifts_right.insert( shift - size_difference);
        }

        // use only shifts, that fulfill conditions on both sides
        std::set_intersection
        (
          possible_shifts_left.begin(), possible_shifts_left.end(),
          possible_shifts_right.begin(), possible_shifts_right.end(),
          std::inserter( possible_shifts, possible_shifts.begin())
        );

        // end
        return possible_shifts;
      }

      //! @brief create alignment node for given pair of alignments and shift
      //! a shift of 3 results in an alignment:
      //! ---ABC
      //! abcd--
      //! while a shift of -2 results an an alignment:
      //! ABC---
      //! --abcd
      //! @param ALIGNMENT_TOP the first alignment to be aligned (no const, it must go in AlignmentNode's m_ChildAlignments)
      //! @param ALIGNMENT_BOTTOM the second alignment to be aligned
      //! @param SHIFT the shift between the top and the bottom alignment desired; positive means gaps on top, negative means gaps on bottom
      //! @return the Alignment with gaps according to given shift
      static AlignmentNode< t_Member> AlignmentFromShift
      (
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_TOP,
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_BOTTOM,
        const int SHIFT
      )
      {
        AlignmentNode< t_Member> current_alignment( ALIGNMENT_TOP, ALIGNMENT_BOTTOM);

        typename AlignmentInterface< t_Member>::const_iterator
          itr_a( ALIGNMENT_TOP->GetAssignments().Begin()), itr_a_end( ALIGNMENT_TOP->GetAssignments().End());
        typename AlignmentInterface< t_Member>::const_iterator
          itr_b( ALIGNMENT_BOTTOM->GetAssignments().Begin()), itr_b_end( ALIGNMENT_BOTTOM->GetAssignments().End());

        // left gap alignments; gaps on top
        for( int i( 0); i < SHIFT && itr_b != itr_b_end; ++i, ++itr_b)
        {
          current_alignment.Append
          (
            util::ShPtr< Assignment< t_Member> >( new Assignment< t_Member>( util::SiPtr< const t_Member>(), **itr_b))
          );
        }
        // left gap alignments; gaps on bottom
        for( int i( 0); i < -SHIFT && itr_a != itr_a_end; ++i, ++itr_a)
        {
          current_alignment.Append
          (
            util::ShPtr< Assignment< t_Member> >( new Assignment< t_Member>( **itr_a, util::SiPtr< const t_Member>()))
          );
        }

        // member alignments
        for( ; itr_a != itr_a_end && itr_b != itr_b_end; ++itr_a, ++itr_b)
        {
          current_alignment.Append
          (
            util::ShPtr< Assignment< t_Member> >( new Assignment< t_Member>( **itr_a, **itr_b))
          );
        }

        // right gap alignments
        for( ; itr_b != itr_b_end; ++itr_b)
        {
          current_alignment.Append
          (
            util::ShPtr< Assignment< t_Member> >( new Assignment< t_Member>( util::SiPtr< const t_Member>(), **itr_b)
            )
          );
        }
        for( ; itr_a != itr_a_end; ++itr_a)
        {
          current_alignment.Append
          (
            util::ShPtr< Assignment< t_Member> >( new Assignment< t_Member>( **itr_a, util::SiPtr< const t_Member>()))
          );
        }

        // end
        return current_alignment;
      }

    }; // template class AlignerShift

    // instantiate s_Instance
    template< typename t_Member>
    const util::SiPtr< const util::ObjectInterface> AlignerShift< t_Member>::s_Instance
    (
      GetObjectInstances().AddInstance( new AlignerShift< t_Member>())
    );

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_ALIGNER_SHIFT_H_ 
