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

#ifndef BCL_ALIGN_ALIGNER_PROGRESSIVE_H_
#define BCL_ALIGN_ALIGNER_PROGRESSIVE_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_alignment_leaf.h"
#include "bcl_align_multiple_aligner_interface.h"
#include "bcl_align_pairwise_aligner_classes.h"
#include "bcl_align_pairwise_aligner_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignerProgressive
    //! @brief This class creates an alignment from multiple sequences by progressive pairwise alignments 
    //!
    //! @tparam t_Member type of object that the Assignment stores
    //!
    //! @see @link example_align_aligner_progressive.cpp @endlink
    //! @author heinzes1
    //! @date Jan 25, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class AlignerProgressive :
      public MultipleAlignerInterface< t_Member>
    {
      typedef storage::Pair< util::ShPtr< AlignmentInterface< t_Member> >, double> alignment_pair_type;

    //////////
    // data //
    //////////

      util::ShPtr< PairwiseAlignerInterface< t_Member> > m_PairwiseAligner; //!< aligner to do the pairwise alignments
      score::AssignmentWithGap< t_Member>                m_ScoreAssignment; //!< ScoreAssignment object
      bool                                               m_FullLinkage;     //!< whether or not each sequence is aligned with every other

    public:

      static const util::SiPtr< const util::ObjectInterface> s_Instance; //!< single instance of that class

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor taking a PairwiseSequenceAligner and a bool whether or not do full linkage alignment
      //! @param PAIRWISE_ALIGNER the aligner to calculate pairwise alignments
      //! @param SCORE the scoring function to use
      //! @param FULL_LINKAGE whether or not calculate full linkage alignments
      AlignerProgressive
      ( 
        const util::ShPtr< PairwiseAlignerInterface< t_Member> > &PAIRWISE_ALIGNER =
          GetPairwiseAlignerClasses< t_Member>().e_AlignerDynamicProgramming->HardCopy(),
        const score::AssignmentWithGap< t_Member> &SCORE = score::AssignmentWithGap< t_Member>(),
        bool FULL_LINKAGE = false
      ) :
        m_PairwiseAligner( PAIRWISE_ALIGNER),
        m_ScoreAssignment( SCORE),
        m_FullLinkage( FULL_LINKAGE)
      {
        m_PairwiseAligner->SetScoringFunction( SCORE);
      }

      //! @brief Clone function
      //! @return pointer to new EngineMerge< t_Member>
      AlignerProgressive< t_Member> *Clone() const
      {
        return new AlignerProgressive< t_Member>( *this);
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

      //! @brief SetScoringFunction sets the score::Assignment scoring function
      //! @param SCORE is the score::Assignment to be used for scoring an assignment
      void SetScoringFunction( const score::AssignmentWithGap< t_Member> &SCORE)
      {
        m_ScoreAssignment = SCORE;
        m_PairwiseAligner->SetScoringFunction( SCORE);
      }

      //! @brief SetPairwiseAligner sets the aligner used to calculate the pairwise alignments
      //! @param ALIGNER a ShPtr to a PairwiseSequenceAligner
      void SetPairwiseAligner( const util::ShPtr< PairwiseAlignerInterface< t_Member> > &PAIRWISE_ALIGNER)
      {
        m_PairwiseAligner = PAIRWISE_ALIGNER;
      }

      //! @brief SetFullLinkage sets bool m_FullLinkage
      //! @param FULL_LINKAGE is bool to set whether or not full linkage alignments are calculated
      void SetFullLinkage( const bool FULL_LINKAGE)
      {
        m_FullLinkage = FULL_LINKAGE;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Align is the function which does the actually aligns data together
      //! @param ALIGNMENTS is the list of Alignments to be aligned
      //! @return returns a pair of the Alignment and a double which is the score
      storage::Pair< AlignmentNode< t_Member>, double> AlignMultiple
      (
        util::ShPtrList< AlignmentInterface< t_Member> > &ALIGNMENTS
      ) const
      {
        if( ALIGNMENTS.GetSize() < 2) // no alignment possible
        {
          return storage::Pair< AlignmentNode< t_Member>, double>( AlignmentNode< t_Member>(), util::GetUndefinedDouble());
        }

        // create align_pairs as Vector of Pairs of Alignments and scores and allocate the necessary memory
        storage::Vector< alignment_pair_type> align_pairs;
        align_pairs.AllocateMemory( ALIGNMENTS.GetSize());

        // iterate through ALIGNMENTS to create an align_pair from each alignment
        for
        (
          typename util::ShPtrList< AlignmentInterface< t_Member> >::const_iterator
            align_itr( ALIGNMENTS.Begin()), align_itr_end( ALIGNMENTS.End());
          align_itr != align_itr_end;
          ++align_itr
        )
        {
          align_pairs.PushBack( alignment_pair_type( *align_itr, util::GetUndefinedDouble()));
        }

        // build start alignments
        BCL_MessageStd( "Start multiple sequence alignment");

        size_t round( 0); // create size_t round and initialize to zero

        // align sequences while align_pairs has at least two Alignments
        while( align_pairs.GetSize() > 1)
        {
          // increase round
          round++;

          // output messages
          BCL_MessageStd( "Current   alignment : " + util::Format().W( 3)( round));
          BCL_MessageStd( "Remaining alignments: " + util::Format().W( 3)( align_pairs.GetSize() - 1));

          // create align_pair and best_align_pair as Pairs of Alignment and double
          alignment_pair_type align_pair, best_align_pair;
          // set the score of best_align_pair to an undefined double
          best_align_pair.Second() = util::GetUndefinedDouble();

          // create iterators best_index and best_index_b, and initialize to the end of align_pairs
          typename std::vector< alignment_pair_type>::iterator
            best_index( align_pairs.End()),
            best_index_b( align_pairs.End());

          // iterate over align_pairs
          typename std::vector< alignment_pair_type>::iterator
            align_pair_itr( align_pairs.Begin()),
            align_pair_itr_end( align_pairs.End());

          for( ; align_pair_itr != align_pair_itr_end; ++align_pair_itr)
          {
            // iterate over align_pairs to compare all Alignments to the Alignment denoted by align_pair_itr
            typename std::vector< alignment_pair_type>::iterator align_pair_itr_b( align_pair_itr + 1);
            for( ; align_pair_itr_b != align_pair_itr_end; ++align_pair_itr_b)
            {
              // if m_FullLinkage is true or any the depth of alignments denoted by align_pair_itr or
              // align_pair_itr_b is the same as the current round
              if( m_FullLinkage || align_pair_itr->First()->GetDepth() == round || align_pair_itr_b->First()->GetDepth() == round)
              {
                // set align_pair to alignment and score from aligning align_pair_itr and align_pair_itr_b
                // first convert ShPtr<AlignmentNode> to ShPtr<AlignmentInterface> for use in AlignPair
                storage::Pair< AlignmentNode< t_Member>, double>
                  result( m_PairwiseAligner->AlignPair( align_pair_itr->First(), align_pair_itr_b->First()));
                align_pair =
                  alignment_pair_type( util::ShPtr< AlignmentInterface< t_Member> >( result.First().Clone()), result.Second());

                // if score of best_alignment_pair is undefined, or the score of align_pair is better than the
                // score of best_align_pair
                if( !util::IsDefined( best_align_pair.Second()) || align_pair.Second() > best_align_pair.Second())
                {
                  // set best_index, best_index_b, and best_align_pair
                  best_index = align_pair_itr;
                  best_index_b = align_pair_itr_b;
                  best_align_pair = align_pair;
                }
              }
            }
          }

          // after calculating the alignment of each alignment to every other alignment, replace alignment and score
          // denoted by best_index with best_align_pair and delete alignment denoted by best_index_b
          BCL_Assert( best_index != align_pairs.End() && best_index_b != align_pairs.End(), "No alignment found");
          *best_index = best_align_pair;
          align_pairs.RemoveElement( best_index_b);
        }

        // return last remaining alignment and score, but convert type i.e. to no ShPtr to alignment
        alignment_pair_type result( align_pairs.FirstElement());
        return storage::Pair< AlignmentNode< t_Member>, double>( *result.First(), result.Second());
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
        return ISTREAM; // return the stream
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM; // return the stream
      }

    }; // template class AlignerProgressive

    // instantiate s_Instance
    template< typename t_Member>
    const util::SiPtr< const util::ObjectInterface> AlignerProgressive< t_Member>::s_Instance
    (
      GetObjectInstances().AddInstance( new AlignerProgressive< t_Member>())
    );

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_ALIGNER_PROGRESSIVE_H_
