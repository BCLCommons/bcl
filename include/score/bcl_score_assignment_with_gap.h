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

#ifndef BCL_SCORE_ASSIGNMENT_WITH_GAP_H_
#define BCL_SCORE_ASSIGNMENT_WITH_GAP_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_score_assignment_gap_simple.h"
#include "align/bcl_align_assignment.h"
#include "function/bcl_function_binary_interface.h"
#include "function/bcl_function_unary_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AssignmentWithGap
    //! @brief This is a template class for scoring assignments
    //!
    //! @see @link example_score_assignment_with_gap.cpp @endlink
    //! @author meilerj
    //! @date 21.11.2005
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class AssignmentWithGap :
      public function::UnaryInterface< const align::Assignment< t_Member>, double>
    {
    private:

    //////////
    // data //
    //////////

      //! Function to score pair of t_GroupMember
      util::ShPtr< function::BinaryInterface< const t_Member, const t_Member, double> > m_ScorePair;

      //! Function to score enclosed gaps
      util::ShPtr< AssignmentGapSimple> m_ScoreGapEnclosedOpen;
      //! Function to score enclosed gaps
      util::ShPtr< AssignmentGapSimple> m_ScoreGapEnclosedExtend;
      //! Function to score boundary gaps
      util::ShPtr< AssignmentGapSimple> m_ScoreGapBoundaryOpen;
      //! Function to score boundary gaps
      util::ShPtr< AssignmentGapSimple> m_ScoreGapBoundaryExtend;

      //! bool to indicate a trivial scoring function to skip the scoring
      bool m_TrivialScoringFunction;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructer
      AssignmentWithGap() :
        m_ScorePair(),
        m_ScoreGapEnclosedOpen( new AssignmentGapSimple( 0.0)),
        m_ScoreGapEnclosedExtend( new AssignmentGapSimple( 0.0)),
        m_ScoreGapBoundaryOpen( new AssignmentGapSimple( 0.0)),
        m_ScoreGapBoundaryExtend( new AssignmentGapSimple( 0.0)),
        m_TrivialScoringFunction( true)
      {
      }

      //! @brief construct from scorepair function and optional parameters for handling gap score
      //! @param SCORE_PAIR util::ShPtr to function to score an assigned pair of amino acids
      //! @param ENCLOSED_SINGLE_GAP   penalty for single   gaps within sequence assignment
      //! @param ENCLOSED_MULTIPLE_GAP penalty for multiple gaps within sequence assignment
      //! @param BOUNDARY_SINGLE_GAP   penalty for single   gaps at boundary of sequence assignment
      //! @param BOUNDARY_MULTIPLE_GAP penalty for multiple gaps at boundary of sequence assignment
      AssignmentWithGap
      (
        const util::ShPtr< function::BinaryInterface< const t_Member, const t_Member, double> > &SCORE_PAIR,
        const double ENCLOSED_SINGLE_GAP = 0.0, const double ENCLOSED_MULTIPLE_GAP = 0.0,
        const double BOUNDARY_SINGLE_GAP = 0.0, const double BOUNDARY_MULTIPLE_GAP = 0.0
      ) :
        m_ScorePair( SCORE_PAIR),
        m_ScoreGapEnclosedOpen( new AssignmentGapSimple( ENCLOSED_SINGLE_GAP)),
        m_ScoreGapEnclosedExtend( new AssignmentGapSimple( ENCLOSED_MULTIPLE_GAP)),
        m_ScoreGapBoundaryOpen( new AssignmentGapSimple( BOUNDARY_SINGLE_GAP)),
        m_ScoreGapBoundaryExtend( new AssignmentGapSimple( BOUNDARY_MULTIPLE_GAP)),
        m_TrivialScoringFunction( false)
      {
        if
        (
          !m_ScorePair.IsDefined() && ENCLOSED_SINGLE_GAP == 0.0 && ENCLOSED_MULTIPLE_GAP == 0.0
          && BOUNDARY_SINGLE_GAP == 0.0 && BOUNDARY_MULTIPLE_GAP == 0.0
        )
        {
          m_TrivialScoringFunction = true;
        }
      }

      //! virtual copy constructor
      AssignmentWithGap< t_Member> *Clone() const
      {
        return new AssignmentWithGap< t_Member>( *this);
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

    ///////////////
    // operators //
    ///////////////

      //! score enclosed gap
      double ScoreGapEnclosedOpen( const size_t SIZE = 1) const
      {
        return m_ScoreGapEnclosedOpen->operator()( SIZE);
      }

      //! score enclosed gap
      double ScoreGapEnclosedExtend( const size_t SIZE = 1) const
      {
        return m_ScoreGapEnclosedExtend->operator()( SIZE);
      }

      //! score enclosed gap
      double ScoreGapBoundaryOpen( const size_t SIZE = 1) const
      {
        return m_ScoreGapBoundaryOpen->operator()( SIZE);
      }

      //! score enclosed gap
      double ScoreGapBoundaryExtend( const size_t SIZE = 1) const
      {
        return m_ScoreGapBoundaryExtend->operator()( SIZE);
      }

      //! score pair
      double operator()
      (
        const util::SiPtr< const t_Member> &MEMBER_A,
        const util::SiPtr< const t_Member> &MEMBER_B
      ) const
      {
        // check for gap
        if( !MEMBER_A.IsDefined() || !MEMBER_B.IsDefined())
        {
          if( !MEMBER_A.IsDefined() && !MEMBER_B.IsDefined())
          {
            return 0.0;
          }

          return m_ScoreGapEnclosedExtend->GetPenalty();
        }

        // if a score was provided
        if( m_ScorePair.IsDefined())
        {
          // score
          return m_ScorePair->operator()( *MEMBER_A, *MEMBER_B);
        }
        else
        {
          return double( 0.0);
        }
      }

      //! score assignment
      double operator()( const align::Assignment< t_Member> &ASSIGNMENT) const;

      //! score pair of assignments
      double operator()
      (
        const align::Assignment< t_Member> &ASSIGNMENT_A,
        const align::Assignment< t_Member> &ASSIGNMENT_B
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        // read members
        io::Serialize::Read( m_ScorePair             , ISTREAM);
        io::Serialize::Read( m_ScoreGapEnclosedOpen  , ISTREAM);
        io::Serialize::Read( m_ScoreGapEnclosedExtend, ISTREAM);
        io::Serialize::Read( m_ScoreGapBoundaryOpen  , ISTREAM);
        io::Serialize::Read( m_ScoreGapBoundaryExtend, ISTREAM);
        io::Serialize::Read( m_TrivialScoringFunction, ISTREAM);

        // end
        return ISTREAM;
      }

      //! write to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members
        io::Serialize::Write( m_ScorePair             , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_ScoreGapEnclosedOpen  , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_ScoreGapEnclosedExtend, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_ScoreGapBoundaryOpen  , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_ScoreGapBoundaryExtend, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_TrivialScoringFunction, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class AssignmentGap

    //! score assignment
    template< typename t_Member>
    double AssignmentWithGap< t_Member>::operator()( const align::Assignment< t_Member> &ASSIGNMENT) const
    {
      // check if the scoring function is trivial and size big enough to score
      if( m_TrivialScoringFunction || ASSIGNMENT.GetMembers().GetSize() <= 1)
      {
        return 0.0;
      }

      // initialize score
      double score( 0.0);

      const util::SiPtrList< const t_Member> &members( ASSIGNMENT.GetMembers());

      // loop over all pairs of t_GroupMember within "big_group"
      for
      (
        typename util::SiPtrList< const t_Member>::const_iterator
          member_itr_a( members.Begin()), member_itr_end( members.End());
        member_itr_a != member_itr_end;
        ++member_itr_a
      )
      {
        // initialize member_itr_b before the loop, it must be increased first
        typename util::SiPtrList< const t_Member>::const_iterator member_itr_b( member_itr_a);
        ++member_itr_b;
        for( ; member_itr_b != member_itr_end; ++member_itr_b)
        {
          score += operator()( *member_itr_a, *member_itr_b);
        }
      }

      // normalize
      score /= double( members.GetSize() * ( members.GetSize() - 1) / 2);

      // return score
      return score;
    }

    //! score pair of assignments
    template< typename t_Member>
    double AssignmentWithGap< t_Member>::operator()
    (
      const align::Assignment< t_Member> &ASSIGNMENT_A,
      const align::Assignment< t_Member> &ASSIGNMENT_B
    ) const
    {
      const util::SiPtrList< const t_Member> &members_a( ASSIGNMENT_A.GetMembers());
      const util::SiPtrList< const t_Member> &members_b( ASSIGNMENT_B.GetMembers());

      // check size
      if( members_a.IsEmpty() || members_b.IsEmpty())
      {
        return 0;
      }

      // initialize score
      double score( 0);

      // loop over all pairs
      for
      (
        typename util::SiPtrList< const t_Member>::const_iterator
          member_a_itr( members_a.Begin()), member_a_itr_end( members_a.End()),
          member_b_itr_end( members_b.End()); // initialize end itr for ASSIGNMENT_B, it does not change
        member_a_itr != member_a_itr_end;
        ++member_a_itr
      )
      {
        for
        (
          typename util::SiPtrList< const t_Member>::const_iterator member_b_itr( members_b.Begin());
          member_b_itr != member_b_itr_end;
          ++member_b_itr
        )
        {
          score += operator()( *member_a_itr, *member_b_itr);
        }
      }

      // normalize
      score /= double( members_a.GetSize() * members_b.GetSize());

      // return score
      return score;
    }

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_ASSIGNMENT_WITH_GAP_H_
