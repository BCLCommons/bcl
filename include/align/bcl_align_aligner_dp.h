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

#ifndef BCL_ALIGN_ALIGNER_DP_H_
#define BCL_ALIGN_ALIGNER_DP_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_alignment_leaf.h"
#include "bcl_align_pairwise_aligner_interface.h"
#include "linal/bcl_linal_matrix.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignerDP
    //! @brief This is a class which uses dynamic programming in order to calculate an Alignment.
    //! @details  This algorithm is taken from the following publication:
    //! Optimal sequence alignment using affine gap costs. Altschul SF, Erickson BW.
    //! Bull Math Biol. 1986; 48(5-6):603-16.
    //!
    //! @tparam t_Member the type of object that the Assignment stores
    //!
    //! @see @link example_align_aligner_dp.cpp @endlink
    //! @author heinzes1
    //! @date Jan 20, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class AlignerDP :
      public PairwiseAlignerInterface< t_Member>
    {

    private:

    //////////
    // data //
    //////////

      score::AssignmentWithGap< t_Member> m_ScoreAssignment; //!< scoring function object

      double m_ScoreGapEnclosedOpen; //!< score enclosed open gap
      double m_ScoreGapEnclosedExtend; //!< score enclosed extend gap
      double m_ScoreGapBoundaryOpen; //!< score boundary open gap
      double m_ScoreGapBoundaryExtend; //!< score enclosed gap

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new AlignerDP< t_Member>
      AlignerDP< t_Member> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief SetScoringFunction sets the score::Assignment scoring function
      //! @param SCORE is the score::Assignment to be used for scoring an assignment
      void SetScoringFunction( const score::AssignmentWithGap< t_Member> &SCORE);

    ////////////////
    // operations //
    ////////////////

      //! @brief Align is the function which does the actually aligns data together
      //! @param ALIGNMENT_A the first alignment to be aligned (no const, it must go in AlignmentNode's m_ChildAlignments)
      //! @param ALIGNMENT_B the second alignment to be aligned
      //! @return pair of the Alignment and a double which is the score
      storage::Pair< AlignmentNode< t_Member>, double> AlignPair
      (
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_A,
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_B
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief computes the score from individual assignments
      //! @param ALIGNMENT is the alignment whose score is calculated
      //! @return a double which is the score
      double CalculateScore( const AlignmentInterface< t_Member> &ALIGNMENT) const;

    }; // template class AlignerDP

    //! @brief Clone function
    //! @return pointer to new AlignerDP< t_Member>
    template< typename t_Member> AlignerDP< t_Member> *AlignerDP< t_Member>::Clone() const
    {
      return new AlignerDP< t_Member>( *this);
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_Member> const std::string &AlignerDP< t_Member>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief SetScoringFunction sets the score::Assignment scoring function
    //! @param SCORE is the score::Assignment to be used for scoring an assignment
    template< typename t_Member> void AlignerDP< t_Member>::SetScoringFunction
    (
      const score::AssignmentWithGap< t_Member> &SCORE
    )
    {
      m_ScoreAssignment = SCORE;
      m_ScoreGapEnclosedOpen = m_ScoreAssignment.ScoreGapEnclosedOpen();
      m_ScoreGapEnclosedExtend = m_ScoreAssignment.ScoreGapEnclosedExtend();
      m_ScoreGapBoundaryOpen = m_ScoreAssignment.ScoreGapBoundaryOpen();
      m_ScoreGapBoundaryExtend = m_ScoreAssignment.ScoreGapBoundaryExtend();
    }

    //! @brief Align is the function which does the actually aligns data together
    //! @param ALIGNMENT_A the first alignment to be aligned (no const, it must go in AlignmentNode's m_ChildAlignments)
    //! @param ALIGNMENT_B the second alignment to be aligned
    //! @return pair of the Alignment and a double which is the score
    template< typename t_Member> storage::Pair< AlignmentNode< t_Member>, double> AlignerDP< t_Member>::AlignPair
    (
      util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_A,
      util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_B
    ) const
    {
      // if an alignment is empty, return empty Alignment and undefined score
      if( !ALIGNMENT_A.IsDefined() || !ALIGNMENT_B.IsDefined() || ALIGNMENT_A->IsEmpty() || ALIGNMENT_B->IsEmpty())
      {
        return storage::Pair< AlignmentNode< t_Member>, double>( AlignmentNode< t_Member>(), util::GetUndefinedDouble());
      }

      // define bit operations for better readability
      #define SET_BIT(bits, which_bit) ((bits) |= (which_bit))
      #define CLEAR_BIT(bits, which_bit) ((bits) &= ~(which_bit))
      #define IS_BIT_SET(bits, which_bit) ((bits) & (which_bit))

      // create matrices
      const size_t rows( ALIGNMENT_A->GetSize() + 1), cols( ALIGNMENT_B->GetSize() + 1);
      linal::Matrix< double> p_matrix( rows, cols, 0.0); // vertical_score_matrix
      linal::Matrix< double> q_matrix( rows, cols, 0.0); // horizontal_score_matrix
      linal::Matrix< double> r_matrix( rows, cols, 0.0); // diagonal_score_matrix
      linal::Matrix< size_t> path_matrix( rows + 1, cols + 1, ( size_t) 0); // path matrix, 7 bits of path information for each position

      // path bits
      enum
      {
        align_vertical_a   = 1 << 0, // 1
        align_horizontal_b = 1 << 1, // 2
        align_diagonal_c   = 1 << 2, // 4
        align_vertical_d   = 1 << 3, // 8
        align_vertical_e   = 1 << 4, // 16
        align_horizontal_f = 1 << 5, // 32
        align_horizontal_g = 1 << 6  // 64
      };

      // initialize matrices
      for( size_t row( 0); row < rows; ++row) // fill first column (col=0)
      {
        q_matrix( row, 0) = -std::numeric_limits< double>::max(); // use -max, larger numbers indicate better score
        r_matrix( row, 0) = row == 0 ? 0 : m_ScoreGapBoundaryOpen + ( row - 1) * m_ScoreGapBoundaryExtend;
      }
      for( size_t col( 0); col < cols; ++col) // fill first row (row=0)
      {
        p_matrix( 0, col) = -std::numeric_limits< double>::max(); // use -max, larger numbers indicate better score
        r_matrix( 0, col) = col == 0 ? 0 : m_ScoreGapBoundaryOpen + ( col - 1) * m_ScoreGapBoundaryExtend;
      }
      SET_BIT( path_matrix( rows, cols), align_diagonal_c);

      // cost assignment: fill matrix interior
      // initialize non-changing values before the loops
      const util::ShPtrList< Assignment< t_Member> > &assignments_a( ALIGNMENT_A->GetAssignments());
      const util::ShPtrList< Assignment< t_Member> > &assignments_b( ALIGNMENT_B->GetAssignments());

      typename AlignmentLeaf< t_Member>::const_iterator itr_row( assignments_a.Begin());
      for( size_t row( 0); row < rows; ++row)
      {
        typename AlignmentLeaf< t_Member>::const_iterator itr_col( assignments_b.Begin());
        for( size_t col( 0); col < cols; ++col)
        {
          // calculate vertical edge (gap) scores
          if( row > 0) // p_matrix( 0, col) is already initialized
          {
            const double vertical_extend_gap_score
            (
              row == rows - 1 || col == 0 || col == cols - 1 ?
                p_matrix( row - 1, col) + m_ScoreGapBoundaryExtend : p_matrix( row - 1, col) + m_ScoreGapEnclosedExtend
            );
            const double vertical_open_gap_score
            (
              row == rows - 1 || col == 0 || col == cols - 1 ?
                r_matrix( row - 1, col) + m_ScoreGapBoundaryOpen : r_matrix( row - 1, col) + m_ScoreGapEnclosedOpen
            );
            p_matrix( row, col) = std::max( vertical_extend_gap_score, vertical_open_gap_score);

            // set vertical edge path bits
            if( p_matrix( row, col) == vertical_extend_gap_score)
            {
              SET_BIT( path_matrix( row - 1, col), align_vertical_d);
            }
            else //if( p_matrix( row, col) == vertical_open_gap_score)
            {
              SET_BIT( path_matrix( row - 1, col), align_vertical_e);
            }
          }

          // calculate horizontal edge (gap) scores
          if( col > 0) // q_matrix( row, 0) is already initialized
          {
            const double horizontal_extend_gap_score
            (
              row == 0 || row == rows - 1 || col == cols - 1 ?
                q_matrix( row, col - 1) + m_ScoreGapBoundaryExtend : q_matrix( row, col - 1) + m_ScoreGapEnclosedExtend
            );
            const double horizontal_open_gap_score
            (
              row == 0 || row == rows - 1 || col == cols - 1 ?
                r_matrix( row, col - 1) + m_ScoreGapBoundaryOpen : r_matrix( row, col - 1) + m_ScoreGapEnclosedOpen
            );
            q_matrix( row, col) = std::max( horizontal_extend_gap_score, horizontal_open_gap_score);

            // set horizontal edge path bits
            if( q_matrix( row, col) == horizontal_extend_gap_score)
            {
              SET_BIT( path_matrix( row, col - 1), align_horizontal_f);
            }
            else //if( q_matrix( row, col) == horizontal_open_gap_score)
            {
              SET_BIT( path_matrix( row, col - 1), align_horizontal_g);
            }
          }

          // determine the minimum score for current position
          if( row > 0 || col > 0) // r_matrix is already initialize at (0,0)
          {
            const double diagonal_score
            (
              row > 0 && col > 0 ?
                r_matrix( row - 1, col - 1) + m_ScoreAssignment( **itr_row, **itr_col) :
                -std::numeric_limits< double>::max() // if one of row and col is 0, assign -max
            );

            r_matrix( row, col) = std::max( p_matrix( row, col), std::max( q_matrix( row, col), diagonal_score));

            // set path bits
            if( r_matrix( row, col) == p_matrix( row, col))
            {
              SET_BIT( path_matrix( row, col), align_vertical_a);
            }
            else if( r_matrix( row, col) == q_matrix( row, col))
            {
              SET_BIT( path_matrix( row, col), align_horizontal_b);
            }
            else // if( r_matrix( row - 1, col - 1) == diagonal_score)
            {
              SET_BIT( path_matrix( row, col), align_diagonal_c);
            }
          }

          if( col > 0)
          {
            ++itr_col; // only move itr after col > 0
          }
        }

        if( row > 0)
        {
          ++itr_row; // only move itr after row > 0
        }
      }

      // path/edge assignment: the meaning of the bits in the path_matrix changes with this loop, see reference
      // accessing row+1 or col+1 is not going out of bounds as path_matrix has size (rows+1, cols+1)
      for( int row( rows - 1); row >= 0; --row) // use int here to allow for checking >= 0
      {
        for( int col( cols - 1); col >= 0; --col)
        {
          // if no optimal path is passing through the current node, remove edges V(row, col), H(row, col), D(row, col)
          if
          (
            (
              !IS_BIT_SET( path_matrix( row + 1, col), align_vertical_a) ||
              !IS_BIT_SET( path_matrix( row, col), align_vertical_e)
            ) &&
            (
              !IS_BIT_SET( path_matrix( row, col + 1), align_horizontal_b) ||
              !IS_BIT_SET( path_matrix( row, col), align_horizontal_g)
            ) &&
            ( !IS_BIT_SET( path_matrix( row + 1, col + 1), align_diagonal_c))
          )
          {
            // set bits to 0
            CLEAR_BIT( path_matrix( row, col), align_vertical_a | align_horizontal_b | align_diagonal_c);
          }

          // if no optimal path exists though the current node (none of these bits are set), go to next cell in matrix
          if
          (
            !IS_BIT_SET( path_matrix( row + 1, col), align_vertical_a) &&
            !IS_BIT_SET( path_matrix( row, col + 1), align_horizontal_b) &&
            !IS_BIT_SET( path_matrix( row + 1, col + 1), align_diagonal_c)
          )
          {
            continue;
          }

          // if V(row+1, col) is an optimal path and requires V(row, col),
          // determine if an optimal path through V(row+1, col) must use V(row, col)
          if
          (
            IS_BIT_SET( path_matrix( row + 1, col), align_vertical_a) &&
            IS_BIT_SET( path_matrix( row, col), align_vertical_d)
          )
          {
            // set bit d to (1 - bit e) i.e. to (! bit e)
            if( IS_BIT_SET( path_matrix( row, col), align_vertical_e))
            {
              CLEAR_BIT( path_matrix( row + 1, col), align_vertical_d);
            }
            else
            {
              SET_BIT( path_matrix( row + 1, col), align_vertical_d);
            }

            // set bit e to (1 - bit a)
            if( IS_BIT_SET( path_matrix( row, col), align_vertical_a))
            {
              CLEAR_BIT( path_matrix( row, col), align_vertical_e);
            }
            else
            {
              SET_BIT( path_matrix( row, col), align_vertical_e);
            }

            // set bit a
            SET_BIT( path_matrix( row, col), align_vertical_a);
          }
          else
          {
            // set bit d and bit e to 0
            CLEAR_BIT( path_matrix( row + 1, col), align_vertical_d);
            CLEAR_BIT( path_matrix( row, col), align_vertical_e);
          }

          // if H(row, col+1) is an optimal path and requires H(row, col),
          // determine if an optimal path through H(row, col+1) must use H(row, col)
          if
          (
            IS_BIT_SET( path_matrix( row, col + 1), align_horizontal_b) &&
            IS_BIT_SET( path_matrix( row, col), align_horizontal_f)
          )
          {
            // set bit f to (1 - bit g) i.e. to (! bit g)
            if( IS_BIT_SET( path_matrix( row, col), align_horizontal_g))
            {
              CLEAR_BIT( path_matrix( row, col + 1), align_horizontal_f);
            }
            else
            {
              SET_BIT( path_matrix( row, col + 1), align_horizontal_f);
            }

            // set bit g to (1 - bit b)
            if( IS_BIT_SET( path_matrix( row, col), align_horizontal_b))
            {
              CLEAR_BIT( path_matrix( row, col), align_horizontal_g);
            }
            else
            {
              SET_BIT( path_matrix( row, col), align_horizontal_g);
            }

            // set bit b
            SET_BIT( path_matrix( row, col), align_horizontal_b);
          }
          else
          {
            // set bit f and bit g to 0
            CLEAR_BIT( path_matrix( row, col + 1), align_horizontal_f);
            CLEAR_BIT( path_matrix( row, col), align_horizontal_g);
          }
        }
      }

      const double score( r_matrix( rows - 1, cols - 1)); // get the score
      storage::Pair< AlignmentNode< t_Member>, double>
        result( AlignmentNode< t_Member>( ALIGNMENT_A, ALIGNMENT_B), score);

      // initialize for traceback: itrs to point to current assignment in ALIGNMENT_A/B, sizes to insert gaps correctly
      typename AlignmentLeaf< t_Member>::const_reverse_iterator r_itr_row( ALIGNMENT_A->GetAssignments().ReverseBegin());
      typename AlignmentLeaf< t_Member>::const_reverse_iterator r_itr_col( ALIGNMENT_B->GetAssignments().ReverseBegin());
      const size_t row_alignment_size( ALIGNMENT_A->GetDepth());
      const size_t col_alignment_size( ALIGNMENT_B->GetDepth());

      // trace alignment back through path matrix
      size_t row( rows - 1), col( cols - 1);
      while( row != 0 || col != 0)
      {
        util::ShPtr< Assignment< t_Member> > new_assignment( new Assignment< t_Member>()); // create new assignment

        if( IS_BIT_SET( path_matrix( row, col), align_diagonal_c))
        {
          new_assignment->Append( ( **r_itr_row).GetMembers());
          new_assignment->Append( ( **r_itr_col).GetMembers());
          --row, --col, ++r_itr_row, ++r_itr_col; // reverse itr need ++ to move forward i.e. backward on the sequence
        }
        else if( IS_BIT_SET( path_matrix( row, col), align_vertical_a))
        {
          new_assignment->Append( ( **r_itr_row).GetMembers());
          new_assignment->Append( util::SiPtrList< const t_Member>( col_alignment_size)); // gaps for col
          --row, ++r_itr_row;
        }
        else if( IS_BIT_SET( path_matrix( row, col), align_horizontal_b))
        {
          new_assignment->Append( util::SiPtrList< const t_Member>( row_alignment_size)); // gaps for row
          new_assignment->Append( ( **r_itr_col).GetMembers());
          --col, ++r_itr_col;
        }
        else if
        (
          IS_BIT_SET( path_matrix( row + 1, col), align_vertical_d)
          && IS_BIT_SET( path_matrix( row + 1, col), align_vertical_a)
        )
        {
          new_assignment->Append( ( **r_itr_row).GetMembers());
          new_assignment->Append( util::SiPtrList< const t_Member>( col_alignment_size)); // gaps for col
          --row, ++r_itr_row;
        }
        else if
        (
          row > 0 && IS_BIT_SET( path_matrix( row - 1, col), align_vertical_e)
          && IS_BIT_SET( path_matrix( row - 1, col), align_vertical_a)
        )
        {
          new_assignment->Append( ( **r_itr_row).GetMembers());
          new_assignment->Append( util::SiPtrList< const t_Member>( col_alignment_size)); // gaps for col
          --row, ++r_itr_row;
        }
        else if
        (
          IS_BIT_SET( path_matrix( row, col + 1), align_horizontal_f)
          && IS_BIT_SET( path_matrix( row, col + 1), align_horizontal_b)
        )
        {
          new_assignment->Append( util::SiPtrList< const t_Member>( row_alignment_size)); // gaps for row
          new_assignment->Append( ( **r_itr_col).GetMembers());
          --col, ++r_itr_col;
        }
        else if( col > 0) // if( col > 0 && IS_BIT_SET( path_matrix( row, col - 1), align_horizontal_g) && IS_BIT_SET( path_matrix( row, col - 1), align_horizontal_b))
        {
          new_assignment->Append( util::SiPtrList< const t_Member>( row_alignment_size)); // gaps for row
          new_assignment->Append( ( **r_itr_col).GetMembers());
          --col, ++r_itr_col;
        }
        else
        {
          BCL_Exit( "Cannot determine path!", 1);
        }

        result.First().Prepend( new_assignment); // new_assignment cannot be empty here (b/c exit), so just prepend
      }

      // undef to only allow local usage
      #undef SET_BIT
      #undef CLEAR_BIT
      #undef IS_BIT_SET

      // rescore if multiple sequence alignment
      if( !result.First().IsEmpty() && ( *result.First().GetAssignments().Begin())->GetMembers().GetSize() > 2)
      {
        result.Second() = CalculateScore( result.First());
      }

      return result;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_Member> std::istream &AlignerDP< t_Member>::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_ScoreAssignment, ISTREAM); // read member
      return ISTREAM; // return the stream
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    template< typename t_Member> std::ostream &AlignerDP< t_Member>::Write
    (
      std::ostream &OSTREAM,
      const size_t INDENT
    ) const
    {
      io::Serialize::Write( m_ScoreAssignment, OSTREAM, INDENT); // write member
      return OSTREAM; // return the stream
    }

    //! @brief computes the score from individual assignments
    //! @param ALIGNMENT is the alignment whose score is calculated
    //! @return a double which is the score
    template< typename t_Member> double AlignerDP< t_Member>::CalculateScore
    (
      const AlignmentInterface< t_Member> &ALIGNMENT
    ) const
    {
      double score( 0);

      //iterate over all
      typename AlignmentLeaf< t_Member>::const_iterator assignment_itr( ALIGNMENT.GetAssignments().Begin());
      typename AlignmentLeaf< t_Member>::const_iterator assignment_itr_end( ALIGNMENT.GetAssignments().End());
      for( ; assignment_itr != assignment_itr_end; ++assignment_itr)
      {
        score += m_ScoreAssignment( **assignment_itr);
      }

      return score;
    }

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_ALIGNER_DP_H_
