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

#ifndef BCL_ALIGN_ALIGNER_WORDBASED_H_
#define BCL_ALIGN_ALIGNER_WORDBASED_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_alignment_hit.h"
#include "bcl_align_alignment_leaf.h"
#include "bcl_align_alignment_node.h"
#include "bcl_align_multiple_aligner_interface.h"
#include "bcl_align_pairwise_aligner_interface.h"
#include "bcl_align_sequence_interface.h"
#include "bcl_align_word_generator_subsequences.h"
#include "score/bcl_score_assignment_with_gap.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignerWordbased
    //! @brief calculates an alignment using a word-based alignment algorithm
    //! @details This class calculates an alignment with a alignment algorithm derived from ideas described in blast
    //!
    //! @see @link example_align_aligner_wordbased.cpp @endlink
    //! @author heinzes1
    //! @date Mar 24, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class AlignerWordbased :
      public PairwiseAlignerInterface< t_Member>
    {

    private:

    //////////
    // data //
    //////////

      score::AssignmentWithGap< t_Member> m_ScoreAssignment; //!< ScoreAssignment object

      size_t m_WordLength; //!< word length
      size_t m_MaximalExtensionWithoutImprovement; //!< max number of AAs tested for an extension w/o score improvement

      static const size_t s_DefaultWordLength = 3; //!< default word length
      static const size_t s_DefaultMaximalExtensionWithoutImprovement = 20; //!< default max extension w/o score improvement

    public:

      static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< single instance of that class

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, optionally taking assignment score, word length and max extension length
      //! @param SCORE_ASSIGNMENT the assignment scoring function to score
      //! @param WORD_LENGTH the length of a word
      //! @param MAXIMAL_EXTENSION_WITHOUT_IMPROVEMENT the maximal extension length without an improvement
      AlignerWordbased
      (
        const score::AssignmentWithGap< t_Member> &SCORE_ASSIGNMENT = score::AssignmentWithGap< t_Member>(),
        const size_t WORD_LENGTH = s_DefaultWordLength,
        const size_t MAXIMAL_EXTENSION_WITHOUT_IMPROVEMENT = s_DefaultMaximalExtensionWithoutImprovement
      ) :
        m_ScoreAssignment( SCORE_ASSIGNMENT),
        m_WordLength( WORD_LENGTH),
        m_MaximalExtensionWithoutImprovement( MAXIMAL_EXTENSION_WITHOUT_IMPROVEMENT)
      {
      }

      //! @brief Clone function
      //! @return pointer to new AlignerWordbased
      AlignerWordbased< t_Member> *Clone() const
      {
        return new AlignerWordbased< t_Member>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief SetScoringFunction sets the score::Assignment scoring function
      //! @param SCORE is the score::Assignment to be used for scoring an assignment
      void SetScoringFunction( const score::AssignmentWithGap< t_Member> &SCORE)
      {
        m_ScoreAssignment = SCORE;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief AlignPair is the function which does the actually aligns data together
      //! @param QUERY_ALIGNMENT first alignment to be aligned (no const, must go in AlignmentNode's m_ChildAlignments)
      //! @param ALIGNMENT second alignment to be aligned
      //! @return returns a pair of the Alignment and a double which is the score
      storage::Pair< AlignmentNode< t_Member>, double> AlignPair
      (
        util::ShPtr< AlignmentInterface< t_Member> > &QUERY_ALIGNMENT,
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT
      ) const;

      //! @brief AlignPairwise is the function which does the actually aligns data together
      //! @param QUERY_ALIGNMENT the alignment to be aligned
      //! @param ALIGNMENT_LIST the list of alignments ALIGNMENT is aligned with
      //! @return returns a list of pairs of the Alignment and doubles which are the score
      storage::List< storage::Pair< AlignmentNode< t_Member>, double> > AlignPairwise
      (
        util::ShPtr< AlignmentInterface< t_Member> > &QUERY_ALIGNMENT,
        util::ShPtrList< AlignmentInterface< t_Member> > &ALIGNMENT_LIST
      ) const;

      //! @brief GenerateHits is the function which generates the hits between the alignments
      //! @param QUERY_ALIGNMENT the alignment to be used to generate words
      //! @param ALIGNMENT_LIST the list of alignments ALIGNMENT is that are searched for those words
      //! @return returns a vector of hits
      storage::List< AlignmentHit< t_Member> > GenerateHits
      (
        util::ShPtr< AlignmentInterface< t_Member> > &QUERY_ALIGNMENT,
        util::ShPtrList< AlignmentInterface< t_Member> > &ALIGNMENT_LIST
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read member
        io::Serialize::Read( m_WordLength     , ISTREAM);
        io::Serialize::Read( m_ScoreAssignment, ISTREAM);

        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_WordLength     , OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_ScoreAssignment, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief scans the given alignment for hits (i.e. high scoring segments aligning a word and an alignment part)
      //! @param ALIGNMENT the alignment which is scanned
      //! @param QUERY_WORDS the vector of words which the alignment is scanned for
      //! @return vector of hits
      storage::List< AlignmentHit< t_Member> > ScanAlignmentForHits
      (
        const util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT,
        const storage::Vector< AlignmentWord< t_Member> > &QUERY_WORDS
      ) const;

      //! @brief returns the assignment build from the two assignments pointed at by the itrs
      //! @param first itr (different itr types possible, fwd itrs, rev itrs etc.)
      //! @param second itr
      //! @return ShPtr to an assignment
      template< typename t_Iterator>
      util::ShPtr< Assignment< t_Member> > GetAssignment( t_Iterator &ITR_A, t_Iterator &ITR_B) const;

      //! @brief extends the given hits in both directions
      //! @param HITS the vector of hits that are extended (no copies are made)
      void ExtendHits( storage::List< AlignmentHit< t_Member> > &HITS) const;

      //! @brief extents the given hit alignment forward using the information from the hit
      //! @param HIT hit from which the hit alignment was generated; contains information about where to start extending
      void ExtendHitForward( AlignmentHit< t_Member> &HIT) const;

      //! @brief extents the given hit alignment backward using the information from the hit
      //! @param HIT hit from which the hit alignment was generated; contains information about where to start extending
      void ExtendHitBackward( AlignmentHit< t_Member> &HIT) const;

    }; // template class AlignmentEngineWordbased

  ///////////////////////////////////
  // public operations definitions //
  ///////////////////////////////////

    //! @brief AlignPair is the function which does the actually aligns data together
    //! @param QUERY_ALIGNMENT first alignment to be aligned (no const, must go in AlignmentNode's m_ChildAlignments)
    //! @param ALIGNMENT second alignment to be aligned
    //! @return returns a pair of the Alignment and a double which is the score
    template< typename t_Member>
    storage::Pair< AlignmentNode< t_Member>, double> AlignerWordbased< t_Member>::AlignPair
    (
      util::ShPtr< AlignmentInterface< t_Member> > &QUERY_ALIGNMENT,
      util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT
    ) const
    {
      // add second alignment to an alignment list and use AlignPairwise()
      util::ShPtrList< AlignmentInterface< t_Member> > alignment_list;
      alignment_list.Append( ALIGNMENT);
      return AlignPairwise( QUERY_ALIGNMENT, alignment_list).FirstElement();
    }

    //! @brief AlignPairwise is the function which does the actually aligns data together
    //! @param QUERY_ALIGNMENT the alignment to be aligned
    //! @param ALIGNMENT_LIST the list of alignments ALIGNMENT is aligned with
    //! @return returns a list of pairs of the Alignment and doubles which are the score
    template< typename t_Member>
    storage::List< storage::Pair< AlignmentNode< t_Member>, double> > AlignerWordbased< t_Member>::AlignPairwise
    (
      util::ShPtr< AlignmentInterface< t_Member> > &QUERY_ALIGNMENT,
      util::ShPtrList< AlignmentInterface< t_Member> > &ALIGNMENT_LIST
    ) const
    {
      storage::List< AlignmentHit< t_Member> > hits( GenerateHits( QUERY_ALIGNMENT, ALIGNMENT_LIST));
      ExtendHits( hits);

      storage::List< storage::Pair< AlignmentNode< t_Member>, double> > alignment_and_score_pairs;
      for
      (
        typename storage::List< AlignmentHit< t_Member> >::const_iterator itr( hits.Begin()), itr_end( hits.End());
        itr != itr_end;
        ++itr
      )
      {
        AlignmentNode< t_Member> alignment( *itr);
        double score( alignment.Score( m_ScoreAssignment));
        alignment_and_score_pairs.Append( storage::Pair< AlignmentNode< t_Member>, double>( alignment, score));
      }

      return alignment_and_score_pairs;
    }

    //! @brief GenerateHits is the function which generates the hits between the alignments
    //! @param QUERY_ALIGNMENT the alignment to be used to generate words
    //! @param ALIGNMENT_LIST the list of alignments ALIGNMENT is that are searched for those words
    //! @return returns a list of hits
    template< typename t_Member>
    storage::List< AlignmentHit< t_Member> > AlignerWordbased< t_Member>::GenerateHits
    (
      util::ShPtr< AlignmentInterface< t_Member> > &QUERY_ALIGNMENT,
      util::ShPtrList< AlignmentInterface< t_Member> > &ALIGNMENT_LIST
    ) const
    {
      if( ALIGNMENT_LIST.IsEmpty()) // no alignment possible
      {
        return storage::List< AlignmentHit< t_Member> >();
      }

      // generate a list/vector of high scoring words for the first alignment
      const WordGeneratorSubsequences< t_Member> generator;
      const storage::Vector< AlignmentWord< t_Member> > query_words( generator.Generate( QUERY_ALIGNMENT, m_WordLength));

      // scan alignment list for hits
      storage::List< AlignmentHit< t_Member> > hits;
      for
      (
        typename util::ShPtrList< AlignmentInterface< t_Member> >::iterator
          itr( ALIGNMENT_LIST.Begin()),
          itr_end( ALIGNMENT_LIST.End());
        itr != itr_end;
        ++itr
      )
      {
        hits.Append( ScanAlignmentForHits( *itr, query_words));
      }

      return hits;
    }

    //! @brief scans the given alignment for hits (i.e. high scoring segments aligning a word and an alignment part)
    //! @param ALIGNMENT the alignment which is scanned
    //! @param QUERY_WORDS the vector of words which the alignment is scanned for
    //! @return vector of hits
    template< typename t_Member>
    storage::List< AlignmentHit< t_Member> > AlignerWordbased< t_Member>::ScanAlignmentForHits
    (
      const util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT,
      const storage::Vector< AlignmentWord< t_Member> > &QUERY_WORDS
    ) const
    {
      storage::List< AlignmentHit< t_Member> > hits;

      // create two itr on alignment in distance of word length
      typename AlignmentInterface< t_Member>::const_iterator
        itr_alignment_hit_begin( ALIGNMENT->GetAssignments().Begin()),
        itr_alignment_hit_end( itr_alignment_hit_begin);
      storage::AdvanceIterator( itr_alignment_hit_end, ALIGNMENT->GetAssignments().End(), m_WordLength);

      // check first assignment separately, then second..end in the loop; iterate over words in wordlist
      for
      (
        typename storage::Vector< AlignmentWord< t_Member> >::const_iterator itr_query_word( QUERY_WORDS.Begin());
        itr_query_word != QUERY_WORDS.End();
        ++itr_query_word
      )
      {
        // clone query word in ShPtr, create second word; create hit from both words; add hit to the vector of hits
        util::ShPtr< AlignmentWord< t_Member> >
          alignment_query_word( util::ShPtr< AlignmentWord< t_Member> >( itr_query_word->Clone())),
          alignment_word( new AlignmentWord< t_Member>( ALIGNMENT, itr_alignment_hit_begin, itr_alignment_hit_end));
        AlignmentHit< t_Member> hit( alignment_query_word, alignment_word);
        hits.PushBack( hit);
      }

      // iterate over assignment in alignment, end when itr_alignment_hit_end itr reaches end of alignment
      do
      {
        ++itr_alignment_hit_begin, ++itr_alignment_hit_end;

        // iterate over words in wordlist
        for
        (
          typename storage::Vector< AlignmentWord< t_Member> >::const_iterator itr_query_word( QUERY_WORDS.Begin());
          itr_query_word != QUERY_WORDS.End();
          ++itr_query_word
        )
        {
          // clone query word in ShPtr, create second word; create hit from both words; add hit to the vector of hits
          util::ShPtr< AlignmentWord< t_Member> >
            alignment_query_word( util::ShPtr< AlignmentWord< t_Member> >( itr_query_word->Clone())),
            alignment_word( new AlignmentWord< t_Member>( ALIGNMENT, itr_alignment_hit_begin, itr_alignment_hit_end));
          AlignmentHit< t_Member> hit( alignment_query_word, alignment_word);
          hits.PushBack( hit);
        }
      }
      while( itr_alignment_hit_end != ALIGNMENT->GetAssignments().End());

      return hits;
    }

    //! @brief extends the given hits in both directions
    //! @param HITS the list of hits
    template< typename t_Member>
    void AlignerWordbased< t_Member>::ExtendHits( storage::List< AlignmentHit< t_Member> > &HITS) const
    {
      // loop over all hits and extend them, first forward, then backward
      for
      (
        typename storage::List< AlignmentHit< t_Member> >::iterator itr( HITS.Begin()), itr_end( HITS.End());
        itr != itr_end;
        ++itr
      )
      {
        ExtendHitForward( *itr); // forward and backward extension
        ExtendHitBackward( *itr);
      }
    }

    //! @brief returns the assignment build from the two assignments pointed at by the itrs
    //! @param first itr
    //! @param second itr
    //! @return ShPtr to an assignment
    template< typename t_Member> template< typename t_Iterator>
    util::ShPtr< Assignment< t_Member> > AlignerWordbased< t_Member>::GetAssignment
    (
      t_Iterator &ITR_A,
      t_Iterator &ITR_B
    ) const
    {
      util::ShPtr< Assignment< t_Member> > assignment( new Assignment< t_Member>()); // create new assignment
      assignment->Append( ( **ITR_A).GetMembers()); // add members for each child alignment
      assignment->Append( ( **ITR_B).GetMembers());
      return assignment;
    }

    //! @brief extents the given hit alignment forward using the information from the hit
    //! @param HIT hit from which the hit alignment was generated; contains information about where to start extending
    template< typename t_Member>
    void AlignerWordbased< t_Member>::ExtendHitForward( AlignmentHit< t_Member> &HIT) const
    {
      double maximal_score( 0.0), current_score( 0.0); // to reset remaining_max_extension
      size_t remaining_max_extension( m_MaximalExtensionWithoutImprovement); // stop extending if no improvement found
      util::ShPtrList< Assignment< t_Member> > assignment_list; // to collect all assignment tested for this extension

      const util::ShPtr< AlignmentInterface< t_Member> >
        word_a( HIT.GetChildAlignments().FirstElement()),
        word_b( HIT.GetChildAlignments().LastElement());
      typename AlignmentInterface< t_Member>::const_iterator
        word_a_end( HIT.GetEndItrA()),
        alignment_a_end( word_a->GetChildAlignments().FirstElement()->GetAssignments().End()),
        word_b_end( HIT.GetEndItrB()),
        alignment_b_end( word_b->GetChildAlignments().FirstElement()->GetAssignments().End());

      // try to extend the hit, go at most until the end of either alignment or max extension w/o score improvement
      while( word_a_end != alignment_a_end && word_b_end != alignment_b_end && remaining_max_extension > 0)
      {
        // create new assignment and add it to assignment list
        util::ShPtr< Assignment< t_Member> > assignment( GetAssignment( word_a_end, word_b_end));
        assignment_list.Append( assignment);

        // score new assignment
        current_score += m_ScoreAssignment( *assignment);
        if( current_score > maximal_score) // if score for new list is higher/better
        {
          maximal_score = current_score; // then update max score
          HIT.Append( assignment_list); // append list to alignment and then remove everything from the list
          assignment_list.Reset();
          remaining_max_extension = m_MaximalExtensionWithoutImprovement; // reset remaining number of AAs to extend
        }
        else
        {
          --remaining_max_extension; // otherwise decrease the remaining number of AAs to extend
        }

        ++word_a_end, ++word_b_end; // move both itrs forward
      }
    } // ExtendHitForward

    //! @brief extents the given hit alignment backward using the information from the hit
    //! @param HIT hit from which the hit alignment was generated; contains information about where to start extending
    template< typename t_Member>
    void AlignerWordbased< t_Member>::ExtendHitBackward( AlignmentHit< t_Member> &HIT) const
    {
      double maximal_score( 0.0), current_score( 0.0); // to reset remaining_max_extension
      size_t remaining_max_extension( m_MaximalExtensionWithoutImprovement); // stop extending if no improvement found
      util::ShPtrList< Assignment< t_Member> > assignment_list; // to collect all assignment tested for this extension

      const util::ShPtr< AlignmentInterface< t_Member> >
        word_a( HIT.GetChildAlignments().FirstElement()),
        word_b( HIT.GetChildAlignments().LastElement());
      typename AlignmentInterface< t_Member>::const_reverse_iterator
        word_a_begin( HIT.GetBeginItrA()),
        alignment_a_end( word_a->GetChildAlignments().FirstElement()->GetAssignments().ReverseEnd()),
        word_b_begin( HIT.GetBeginItrB()),
        alignment_b_end( word_b->GetChildAlignments().FirstElement()->GetAssignments().ReverseEnd());

      // try to extend the hit, go at most until the end of either alignment or max extension w/o score improvement
      while( word_a_begin != alignment_a_end && word_b_begin != alignment_b_end && remaining_max_extension > 0)
      {
        // create new assignment and add it to assignment list
        util::ShPtr< Assignment< t_Member> > assignment( GetAssignment( word_a_begin, word_b_begin));
        assignment_list.Append( assignment);

        // score new assignment
        current_score += m_ScoreAssignment( *assignment);
        if( current_score > maximal_score) // if score for new list is higher/better
        {
          maximal_score = current_score; // then update max score
          HIT.Prepend( assignment_list); // append list to alignment and then remove everything from the list
          assignment_list.Reset();
          remaining_max_extension = m_MaximalExtensionWithoutImprovement; // reset remaining number of AAs to extend
        }
        else
        {
          --remaining_max_extension; // otherwise decrease the remaining number of AAs to extend
        }

        ++word_a_begin, ++word_b_begin; // move both (reverse) itrs forward
      }
    } // ExtendHitBackward

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_ALIGNER_WORDBASED_H_
