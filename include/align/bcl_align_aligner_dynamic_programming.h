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

#ifndef BCL_ALIGN_ALIGNER_DYNAMIC_PROGRAMMING_H_
#define BCL_ALIGN_ALIGNER_DYNAMIC_PROGRAMMING_H_

// include the namespace header
#include "bcl_align.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_align_alignment_leaf.h"
#include "bcl_align_alignment_node.h"
#include "bcl_align_handler_pir.h"
#include "bcl_align_pairwise_aligner_interface.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "score/bcl_score_assignment_with_gap.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace align
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlignerDynamicProgramming
    //! @brief This is a class which uses dynamic programming in order to calculate an Alignment.
    //!
    //! @tparam t_Member the type of object that the Assignment stores
    //!
    //! @see @link example_align_aligner_dynamic_programming.cpp @endlink
    //! @author alexanns, meilerj
    //! @date Nov 2, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_Member>
    class AlignerDynamicProgramming :
      public PairwiseAlignerInterface< t_Member>
    {
    private:

    //////////
    // data //
    //////////

      score::AssignmentWithGap< t_Member> m_ScoreAssignment; //!< ScoreAssignment object

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, optionally taking a score::Assignment
      //! @param SCORE_ASSIGNMENT is the score::Assignment object which will initialize "m_ScoreAssignment"
      AlignerDynamicProgramming
      (
        const score::AssignmentWithGap< t_Member> &SCORE_ASSIGNMENT = score::AssignmentWithGap< t_Member>()
      ) :
        m_ScoreAssignment( SCORE_ASSIGNMENT)
      {
      }

      //! @brief Clone is the virtual copy constructor
      AlignerDynamicProgramming< t_Member> *Clone() const
      {
        return new AlignerDynamicProgramming< t_Member>( *this);
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
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief AlignPair is the overwritten function of PairwiseSequenceAligner doing the alignment of two Alignments
      //! @param ALIGNMENT_VERTICAL is the first alignment to be aligned
      //! @param ALIGNMENT_HORIZONTAL is the second alignment to be aligned
      //! @return returns a pair of the Alignment and a double which is the score
      storage::Pair< AlignmentNode< t_Member>, double> AlignPair
      (
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_VERTICAL,
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_HORIZONTAL
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief write container to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_ScoreAssignment, OSTREAM, INDENT);

        return OSTREAM;
      }

      //! @brief read container from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read in member
        io::Serialize::Read( m_ScoreAssignment, ISTREAM);

        return ISTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief FillScoreAndGapMatrices calculates the scoring matrix and any gaps associated with cells of the matrix
      //! @param ALIGNMENT_VERTICAL the first alignment involved in the Alignment
      //! @param ALIGNMENT_HORIZONTAL the second alignment involved in the Alignment
      //! @param SCORE_MATRIX the matrix which will hold all the scores calculated
      //! @param HORIZONTAL_GAP_MATRIX matrix to hold the horizontal gaps which are associated with each cell
      //! @param VERTICAL_GAP_MATRIX matrix to hold the vertical gaps which are associated with each cell
      void FillScoreAndGapMatrices
      (
        const AlignmentInterface< t_Member> &ALIGNMENT_VERTICAL,
        const AlignmentInterface< t_Member> &ALIGNMENT_HORIZONTAL,
        linal::Matrix< double> &SCORE_MATRIX,
        linal::MatrixInterface< size_t> &HORIZONTAL_GAP_MATRIX,
        linal::MatrixInterface< size_t> &VERTICAL_GAP_MATRIX
      ) const;

      //! @brief FillInitialBoundaryScoresAndGaps calculates the scores and gaps for the upper and left matrix sides
      //! @param HORIZONTAL_DIRECTION_MAX_LENGTH
      //! @param VERTICAL_DIRECTION_MAX_LENGTH
      //! @param ITR_ROW iterator to first element of "ALIGNMENT_VERTICAL"
      //! @param ITR_COLUMN iterator to first element of "ALIGNMENT_HORIZONTAL"
      //! @param SCORE_MATRIX the matrix which will hold all the scores calculated
      //! @param HORIZONTAL_GAP_MATRIX matrix to hold the horizontal gaps which are associated with each cell
      //! @param VERTICAL_GAP_MATRIX matrix to hold the vertical gaps which are associated with each cell
      //! @param ROW_MATCH_MAX_SCORE Vector which holds best score of any matched cell in each row
      //! @param ROW_MATCH_MAX_SCORE_POSITION Vector holds the position in each row that the match_max_score occurs
      //! @param COLUMN_MATCH_MAX_SCORE Vector which holds best score of any matched cell in each column
      //! @param COLUMN_MATCH_MAX_SCORE_POSITION  Vector holds position in each column that the match_max_score occurs
      void FillInitialBoundaryScoresAndGaps
      (
        const typename AlignmentLeaf< t_Member>::const_iterator &ITR_ROW,
        const typename AlignmentLeaf< t_Member>::const_iterator &ITR_COLUMN,
        const size_t &HORIZONTAL_DIRECTION_MAX_LENGTH,
        const size_t &VERTICAL_DIRECTION_MAX_LENGTH,
        linal::Matrix< double> &SCORE_MATRIX,
        linal::MatrixInterface< size_t> &HORIZONTAL_GAP_MATRIX,
        linal::MatrixInterface< size_t> &VERTICAL_GAP_MATRIX,
        storage::Vector< double> &ROW_MATCH_MAX_SCORE,
        storage::Vector< size_t> &ROW_MATCH_MAX_SCORE_POSITION,
        storage::Vector< double> &COLUMN_MATCH_MAX_SCORE,
        storage::Vector< size_t> &COLUMN_MATCH_MAX_SCORE_POSITION
      ) const;

      //! @brief DeterminePathAndScore calculates best score for a cell and path it must come from to get that score
      //! @param ROW current row in matrix
      //! @param COLUMN current column in matrix
      //! @param ITR_ROW current element of VERTICAL_ALIGNMENT
      //! @param ITR_COLUMN current element of HORIZONTAL_ALIGNMENT
      //! @param SCORE_MATRIX the matrix which will hold all the scores calculated
      //! @param HORIZONTAL_GAP_MATRIX matrix to hold the horizontal gaps which are associated with each cell
      //! @param VERTICAL_GAP_MATRIX matrix to hold the vertical gaps which are associated with each cell
      //! @param ROW_MATCH_MAX_SCORE Vector which holds best score of any matched cell in each row
      //! @param ROW_MATCH_MAX_SCORE_POSITION Vector holds the position in each row that the match_max_score occurs
      //! @param COLUMN_MATCH_MAX_SCORE Vector which holds best score of any matched cell in each column
      //! @param COLUMN_MATCH_MAX_SCORE_POSITION  Vector holds position in each column that the match_max_score occurs
      void DeterminePathAndScore
      (
        const size_t &ROW,
        const size_t &COLUMN,
        const typename AlignmentLeaf< t_Member>::const_iterator &ITR_ROW,
        const typename AlignmentLeaf< t_Member>::const_iterator &ITR_COLUMN,
        linal::Matrix< double> &SCORE_MATRIX,
        linal::MatrixInterface< size_t> &HORIZONTAL_GAP_MATRIX,
        linal::MatrixInterface< size_t> &VERTICAL_GAP_MATRIX,
        storage::Vector< double> &ROW_MATCH_MAX_SCORE,
        storage::Vector< size_t> &ROW_MATCH_MAX_SCORE_POSITION,
        storage::Vector< double> &COLUMN_MATCH_MAX_SCORE,
        storage::Vector< size_t> &COLUMN_MATCH_MAX_SCORE_POSITION
      ) const;

      //! @brief SetPathAndScore decides where a cell is coming from based on the possible pathway scores
      //! @param DIAGONAL_PATH_SCORE score cell would acquire coming from diagonal path
      //! @param HORIZONTAL_PATH_SCORE score cell would acquire coming from horizontal neighbor
      //! @param VERTICAL_PATH_SCORE score cell would acquire coming from vertical neighbor
      //! @param PATH_ROW_MAX_SCORE score cell would acquire coming from cell with maximum score in row
      //! @param PATH_COLUMN_MAX_SCORE score cell would acquire coming from cell with maximum score in column
      //! @param ROW current row in matrix
      //! @param COLUMN current column in matrix
      //! @param SCORE_MATRIX the matrix which will hold all the scores calculated
      //! @param HORIZONTAL_GAP_MATRIX matrix to hold the horizontal gaps which are associated with each cell
      //! @param VERTICAL_GAP_MATRIX matrix to hold the vertical gaps which are associated with each cell
      //! @param ROW_MATCH_MAX_SCORE Vector which holds best score of any matched cell in each row
      //! @param ROW_MATCH_MAX_SCORE_POSITION Vector holds the position in each row that the match_max_score occurs
      //! @param COLUMN_MATCH_MAX_SCORE Vector which holds best score of any matched cell in each column
      //! @param COLUMN_MATCH_MAX_SCORE_POSITION  Vector holds position in each column that the match_max_score occurs
      void SetPathAndScore
      (
        const double &DIAGONAL_PATH_SCORE,
        const double &HORIZONTAL_PATH_SCORE,
        const double &VERTICAL_PATH_SCORE,
        const double &PATH_ROW_MAX_SCORE,
        const double &PATH_COLUMN_MAX_SCORE,
        const size_t &ROW,
        const size_t &COLUMN,
        linal::Matrix< double> &SCORE_MATRIX,
        linal::MatrixInterface< size_t> &HORIZONTAL_GAP_MATRIX,
        linal::MatrixInterface< size_t> &VERTICAL_GAP_MATRIX,
        storage::Vector< double> &ROW_MATCH_MAX_SCORE,
        storage::Vector< size_t> &ROW_MATCH_MAX_SCORE_POSITION,
        storage::Vector< double> &COLUMN_MATCH_MAX_SCORE,
        storage::Vector< size_t> &COLUMN_MATCH_MAX_SCORE_POSITION
      ) const;

      //! @brief TraceBack builds the new alignment
      //! @param ALIGNMENT_VERTICAL the first alignment involved in the Alignment
      //! @param ALIGNMENT_HORIZONTAL the second alignment involved in the Alignment
      //! @param SCORE_MATRIX the matrix which holds all the scores calculated
      //! @param HORIZONTAL_GAP_MATRIX matrix holds the horizontal gaps which are associated with each cell
      //! @param VERTICAL_GAP_MATRIX matrix holds the vertical gaps which are associated with each cell
      //! @return returns the new alignment and its score in a storage::Pair
      storage::Pair< AlignmentNode< t_Member>, double> TraceBack
      (
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_VERTICAL,
        util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_HORIZONTAL,
        const linal::Matrix< double> &SCORE_MATRIX,
        const linal::MatrixConstInterface< size_t> &HORIZONTAL_GAP_MATRIX,
        const linal::MatrixConstInterface< size_t> &VERTICAL_GAP_MATRIX
      ) const;

      //! @brief BuildInteriorOfAlignment builds new alignment until either the upper or left side of matrix is reached
      //! @param ROW index of last row in the matrix
      //! @param COLUMN index of last column in matrix
      //! @param ALIGNMENT_VERTICAL the first alignment involved in the Alignment
      //! @param ALIGNMENT_HORIZONTAL the second alignment involved in the Alignment
      //! @param REV_ROW_ITR reverse iterator on reverse begin of ALIGNMENT_VERTICAL
      //! @param REV_COLUMN_ITR reverse iterator on reverse begin of ALIGNMENT_HORIZONTAL
      //! @param ALIGNMENT new alignment which is being built
      //! @param SCORE_MATRIX the matrix which holds all the scores calculated
      //! @param HORIZONTAL_GAP_MATRIX matrix holds the horizontal gaps which are associated with each cell
      //! @param VERTICAL_GAP_MATRIX matrix holds the vertical gaps which are associated with each cell
      void BuildInteriorOfAlignment
      (
        size_t &ROW,
        size_t &COLUMN,
        const AlignmentInterface< t_Member> &ALIGNMENT_VERTICAL,
        const AlignmentInterface< t_Member> &ALIGNMENT_HORIZONTAL,
        typename AlignmentLeaf< t_Member>::const_reverse_iterator &REV_ROW_ITR,
        typename AlignmentLeaf< t_Member>::const_reverse_iterator &REV_COLUMN_ITR,
        AlignmentNode< t_Member> &ALIGNMENT,
        const linal::Matrix< double> &SCORE_MATRIX,
        const linal::MatrixConstInterface< size_t> &HORIZONTAL_GAP_MATRIX,
        const linal::MatrixConstInterface< size_t> &VERTICAL_GAP_MATRIX
      ) const;

      //! @brief AddGapsToBeginningOfAlignment adds remaining assignments after finished with BuildInteriorOfAlignment
      //! @param ROW current row in matrix
      //! @param COLUMN current column in matrix
      //! @param REV_ROW_ITR iterator for current position of REV_ROW_ITR after BuildInteriorOfAlignment
      //! @param REV_ROW_ITR_END iterator at reverse end of ALIGNMENT_VERTICAL
      //! @param REV_COLUMN_ITR iterator for current position of REV_COLUMN_ITR after BuildInteriorOfAlignment
      //! @param REV_COLUMN_ITR_END iterator at reverse end of ALIGNMENT_HORIZONTAL
      //! @param ALIGNMENT new alignment which is being built
      //! @param HORIZONTAL_GAP_MATRIX matrix holds the horizontal gaps which are associated with each cell
      //! @param VERTICAL_GAP_MATRIX matrix holds the vertical gaps which are associated with each cell
      void AddGapsToBeginningOfAlignment
      (
        size_t &ROW,
        size_t &COLUMN,
        typename AlignmentLeaf< t_Member>::const_reverse_iterator &REV_ROW_ITR,
        typename AlignmentLeaf< t_Member>::const_reverse_iterator &REV_ROW_ITR_END,
        typename AlignmentLeaf< t_Member>::const_reverse_iterator &REV_COLUMN_ITR,
        typename AlignmentLeaf< t_Member>::const_reverse_iterator &REV_COLUMN_ITR_END,
        AlignmentNode< t_Member> &ALIGNMENT,
        const linal::MatrixConstInterface< size_t> &HORIZONTAL_GAP_MATRIX,
        const linal::MatrixConstInterface< size_t> &VERTICAL_GAP_MATRIX
      ) const;

      //! @brief AddAssignmentFromDiagonalPath adds an assignment of assigned elements to the new alignment
      //! @param ITERATOR_ALIGNMENT_VERTICAL iterator to current element of ALIGNMENT_VERTICAL being assigned
      //! @param ITERATOR_ALIGNMENT_HORIZONTAL iterator to current element of ALIGNMENT_HORIZONTAL being assigned
      //! @param ALIGNMENT is the new alignment which is being built
      void AddAssignmentFromDiagonalPath
      (
        const typename AlignmentLeaf< t_Member>::const_reverse_iterator &ITERATOR_ALIGNMENT_VERTICAL,
        const typename AlignmentLeaf< t_Member>::const_reverse_iterator &ITERATOR_ALIGNMENT_HORIZONTAL,
        AlignmentNode< t_Member> &ALIGNMENT
      ) const;

      //! @brief AddAssignmentFromHorizontalPath adds a gapped assignment
      //! @param ITERATOR_ALIGNMENT_VERTICAL iterator to current element of ALIGNMENT_VERTICAL; not being assigned
      //! @param ITERATOR_ALIGNMENT_HORIZONTAL itr to current element of ALIGNMENT_HORIZONTAL being assigned with gaps
      //! @param ALIGNMENT is the new alignment which is being built
      void AddAssignmentFromHorizontalPath
      (
        const typename AlignmentLeaf< t_Member>::const_reverse_iterator &ITERATOR_ALIGNMENT_VERTICAL,
        const typename AlignmentLeaf< t_Member>::const_reverse_iterator &ITERATOR_ALIGNMENT_HORIZONTAL,
        AlignmentNode< t_Member> &ALIGNMENT
      ) const;

      //! @brief AddAssignmentFromVerticalPath adds a gapped assignment
      //! @param ITERATOR_ALIGNMENT_VERTICAL itr to current element of ALIGNMENT_VERTICAL being assigned with gaps
      //! @param ITERATOR_ALIGNMENT_HORIZONTAL iterator to current element of ALIGNMENT_HORIZONTAL; not being assigned
      //! @param ALIGNMENT is the new alignment which is being built
      void AddAssignmentFromVerticalPath
      (
        const typename AlignmentLeaf< t_Member>::const_reverse_iterator &ITERATOR_ALIGNMENT_VERTICAL,
        const typename AlignmentLeaf< t_Member>::const_reverse_iterator &ITERATOR_ALIGNMENT_HORIZONTAL,
        AlignmentNode< t_Member> &ALIGNMENT
      ) const;

      //! @brief UpdateScore computes the Alignment score from individual Assignments
      //! @param ALIGNMENT is the Alignment whose score is calculated
      //! @return returns a double which is the score of ALIGNMENT
      double UpdateScore( AlignmentInterface< t_Member> &ALIGNMENT) const;

    }; // template class AlignerDynamicProgramming

  ///////////////////////////////////
  // public operations definitions //
  ///////////////////////////////////

    //! @brief AlignPair performs the aligning of two alignments
    //! @param ALIGNMENT_VERTICAL the first Alignment
    //! @param ALIGNMENT_HORIZONTAL the second Alignment which will be added to ALIGNMENT_VERTICAL
    //! @return returns a pair of the alignment of the two sequences and its score
    template< typename t_Member>
    storage::Pair< AlignmentNode< t_Member>, double> AlignerDynamicProgramming< t_Member>::AlignPair
    (
      util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_VERTICAL,
      util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_HORIZONTAL
    ) const
    {
      // if ALIGNMENT_VERTICAL or ALIGNMENT_HORIZONTAL are empty, return empty Alignment and undefined score
      if( !ALIGNMENT_VERTICAL.IsDefined() || !ALIGNMENT_HORIZONTAL.IsDefined() || ALIGNMENT_VERTICAL->IsEmpty() || ALIGNMENT_HORIZONTAL->IsEmpty())
      {
        return storage::Pair< AlignmentNode< t_Member>, double>( AlignmentNode< t_Member>(), util::GetUndefinedDouble());
      }

      // create vertical_length and horizontal_length
      size_t vertical_length( ALIGNMENT_VERTICAL->GetSize());
      size_t horizontal_length( ALIGNMENT_HORIZONTAL->GetSize());

      // create score_matrix and initialize to dimensions of vertical_length and horizontal_length
      linal::Matrix< double> score_matrix( vertical_length, horizontal_length);

      // horizontal_gap_matrix/vertical_gap_matrix keep track of how many horizontal/vertical gaps led to a given
      // (row, column) of score_matrix
      linal::Matrix< size_t> horizontal_gap_matrix( vertical_length, horizontal_length);
      linal::Matrix< size_t> vertical_gap_matrix( vertical_length, horizontal_length);

      // fill in score_matrix, horizontal_gap_matrix, and vertical_gap_matrix
      FillScoreAndGapMatrices( *ALIGNMENT_VERTICAL, *ALIGNMENT_HORIZONTAL, score_matrix, horizontal_gap_matrix, vertical_gap_matrix);

      // do TraceBack to get the alignment of ALIGNMENT_VERTICAL and ALIGNMENT_HORIZONTAL and score
      storage::Pair< AlignmentNode< t_Member>, double> alignment_and_score
      (
        TraceBack( ALIGNMENT_VERTICAL, ALIGNMENT_HORIZONTAL, score_matrix, horizontal_gap_matrix, vertical_gap_matrix)
      );

      // rescore if multiple sequence alignment
      if( !alignment_and_score.First().IsEmpty() && ( *alignment_and_score.First().GetAssignments().Begin())->GetMembers().GetSize() > 2)
      {
        alignment_and_score.Second() = UpdateScore( alignment_and_score.First());
      }

      return alignment_and_score;
    } // AlignPair

  ////////////////////////////////////
  // private operations definitions //
  ////////////////////////////////////

    //! @brief FillScoreAndGapMatrices calculates the scoring matrix and any gaps associated with cells of the matrix
    //! @param ALIGNMENT_VERTICAL the first alignment involved in the Alignment
    //! @param ALIGNMENT_HORIZONTAL the second alignment involved in the Alignment
    //! @param SCORE_MATRIX the matrix which will hold all the scores calculated
    //! @param HORIZONTAL_GAP_MATRIX matrix to hold the horizontal gaps which are associated with each cell
    //! @param VERTICAL_GAP_MATRIX matrix to hold the vertical gaps which are associated with each cell
    template< typename t_Member>
    void AlignerDynamicProgramming< t_Member>::FillScoreAndGapMatrices
    (
      const AlignmentInterface< t_Member> &ALIGNMENT_VERTICAL,
      const AlignmentInterface< t_Member> &ALIGNMENT_HORIZONTAL,
      linal::Matrix< double> &SCORE_MATRIX,
      linal::MatrixInterface< size_t> &HORIZONTAL_GAP_MATRIX,
      linal::MatrixInterface< size_t> &VERTICAL_GAP_MATRIX
    ) const
    {
      // create size_t length_vertical and length_horizontal
      const size_t length_vertical( ALIGNMENT_VERTICAL.GetSize());
      const size_t length_horizontal( ALIGNMENT_HORIZONTAL.GetSize());

      // row_match_max_score is, for a row, the max score of any cell which comes from a diagonal pathway
      storage::Vector< double> row_match_max_score( length_vertical, util::GetUndefinedDouble());
      // row_match_max_score_position is the position where row_match_max_score occurs in the row
      storage::Vector< size_t> row_match_max_score_position( length_vertical);

      // column_match_max_score is, for a column, the max score of any cell which comes from a diagonal pathway
      storage::Vector< double> column_match_max_score( length_horizontal, util::GetUndefinedDouble());
      // column_match_max_score_position is the position where column_match_max_score occurs in the column
      storage::Vector< size_t> column_match_max_score_position( length_horizontal);

      // fill in the upper and left edges (boundaries) of the score matrix
      FillInitialBoundaryScoresAndGaps
      (
        ALIGNMENT_VERTICAL.GetAssignments().Begin(),
        ALIGNMENT_HORIZONTAL.GetAssignments().Begin(),
        length_horizontal,
        length_vertical,
        SCORE_MATRIX,
        HORIZONTAL_GAP_MATRIX,
        VERTICAL_GAP_MATRIX,
        row_match_max_score,
        row_match_max_score_position,
        column_match_max_score,
        column_match_max_score_position
      );

      // prepare for filling interior of "SCORE_MATRIX"
      size_t row( 1); //< starting position in vertical direction
      size_t column( 1); //< starting position in horizontal direction

      // iterate through vertical direction of "SCORE_MATRIX"
      typename AlignmentLeaf< t_Member>::const_iterator itr_row( ++ALIGNMENT_VERTICAL.GetAssignments().Begin());
      typename AlignmentLeaf< t_Member>::const_iterator itr_row_end( ALIGNMENT_VERTICAL.GetAssignments().End());
      for( ; row != length_vertical; ++itr_row, ++row)
      {
        // make sure itr_row does not reach the end of ALIGNMENT_VERTICAL before row reaches length_vertical
        BCL_Assert( itr_row != itr_row_end, "Row iterator reached end, but row index did not");

        // iterate through horizontal direction of SCORE_MATRIX
        typename AlignmentLeaf< t_Member>::const_iterator itr_column( ++ALIGNMENT_HORIZONTAL.GetAssignments().Begin());
        typename AlignmentLeaf< t_Member>::const_iterator itr_column_end( ALIGNMENT_HORIZONTAL.GetAssignments().End());
        for( ; column != length_horizontal; ++itr_column, ++column)
        {
          // make sure itr_column does not reach the end of ALIGNMENT_HORIZONTAL before column reaches length_horizontal
          BCL_Assert( itr_column != itr_column_end, "Column iterator reached end, but column index did not");

          // determine path the cell (row, column) comes from and its score
          // also determine the gap matrices which tell how many gaps (if any) have led to cell (row, column)
          DeterminePathAndScore
          (
            row,
            column,
            itr_row,
            itr_column,
            SCORE_MATRIX,
            HORIZONTAL_GAP_MATRIX,
            VERTICAL_GAP_MATRIX,
            row_match_max_score,
            row_match_max_score_position,
            column_match_max_score,
            column_match_max_score_position
          );
        }

        // reset column in preparation of going to next row
        column = 1;
      }
    } // FillScoreAndGapMatrices

    //! @brief FillInitialBoundaryScoresAndGaps calculates the scores and gaps for the upper and left matrix sides
    //! @param HORIZONTAL_DIRECTION_MAX_LENGTH
    //! @param VERTICAL_DIRECTION_MAX_LENGTH
    //! @param ITR_ROW iterator to first element of "ALIGNMENT_VERTICAL"
    //! @param ITR_COLUMN iterator to first element of "ALIGNMENT_HORIZONTAL"
    //! @param SCORE_MATRIX the matrix which will hold all the scores calculated
    //! @param HORIZONTAL_GAP_MATRIX matrix to hold the horizontal gaps which are associated with each cell
    //! @param VERTICAL_GAP_MATRIX matrix to hold the vertical gaps which are associated with each cell
    //! @param ROW_MATCH_MAX_SCORE Vector which holds best score of any matched cell in each row
    //! @param ROW_MATCH_MAX_SCORE_POSITION Vector holds the position in each row that the match_max_score occurs
    //! @param COLUMN_MATCH_MAX_SCORE Vector which holds best score of any matched cell in each column
    //! @param COLUMN_MATCH_MAX_SCORE_POSITION  Vector holds position in each column that the match_max_score occurs
    template< typename t_Member>
    void AlignerDynamicProgramming< t_Member>::FillInitialBoundaryScoresAndGaps
    (
      const typename AlignmentLeaf< t_Member>::const_iterator &ITR_ROW,
      const typename AlignmentLeaf< t_Member>::const_iterator &ITR_COLUMN,
      const size_t &HORIZONTAL_DIRECTION_MAX_LENGTH,
      const size_t &VERTICAL_DIRECTION_MAX_LENGTH,
      linal::Matrix< double> &SCORE_MATRIX,
      linal::MatrixInterface< size_t> &HORIZONTAL_GAP_MATRIX,
      linal::MatrixInterface< size_t> &VERTICAL_GAP_MATRIX,
      storage::Vector< double> &ROW_MATCH_MAX_SCORE,
      storage::Vector< size_t> &ROW_MATCH_MAX_SCORE_POSITION,
      storage::Vector< double> &COLUMN_MATCH_MAX_SCORE,
      storage::Vector< size_t> &COLUMN_MATCH_MAX_SCORE_POSITION
    ) const
    {
      // initialize bool column_origin_found to false; this is for origins in the left boundary which is column 0
      bool column_origin_found( false);
      // initialize bool row_origin_found to false; this is for origins in the upper boundary which is row 0
      bool row_origin_found( false);

      // top left corner cell (0, 0)
      // initialize assign_score to assignment score of ITR_ROW and ITR_COLUMN (row 0, column 0)
      double assign_score( m_ScoreAssignment( **ITR_ROW, **ITR_COLUMN));
      // initialize gap_score to score of opening a gap on the boundary
      double gap_score( m_ScoreAssignment.ScoreGapBoundaryOpen( 1));

      if( assign_score >= gap_score) // assign ITR_ROW and ITR_COLUMN (row 0, column 0) together
      {
        SCORE_MATRIX( 0, 0) = assign_score;
        HORIZONTAL_GAP_MATRIX( 0, 0) = VERTICAL_GAP_MATRIX( 0, 0) = 0;
        ROW_MATCH_MAX_SCORE( 0) = COLUMN_MATCH_MAX_SCORE( 0) = assign_score;
        ROW_MATCH_MAX_SCORE_POSITION( 0) = COLUMN_MATCH_MAX_SCORE_POSITION( 0) = 0;
        column_origin_found = row_origin_found = true;
        BCL_MessageDbg( "Cell ( 0, 0) origin found");
      }
      else if( assign_score < gap_score) // make a gap
      {
        SCORE_MATRIX( 0, 0) = gap_score;
        HORIZONTAL_GAP_MATRIX( 0, 0) = VERTICAL_GAP_MATRIX( 0, 0) = 1;
      }
      else // NaN?
      {
        BCL_Exit
        (
          "Unable to determine (0, 0) while trying to align\n"
            + util::Format()( *ITR_ROW) + "\nand\n" + util::Format()( *ITR_COLUMN)
            + "\nassign_score=" + util::Format()( assign_score) + ", gap_score=" + util::Format()( gap_score),
          -1
        );
      }

      // create const_iterator itr_row and itr_column and initialize to ITR_ROW and ITR_COLUMN respectively
      typename AlignmentLeaf< t_Member>::const_iterator itr_row( ITR_ROW), itr_column( ITR_COLUMN);

      // cell (1, 0)
      // move itr_row to next position
      ++itr_row;
      // set assign_score to score of assigning row 1 with column 0 plus opening a boundary gap
      assign_score = m_ScoreAssignment( **itr_row, **itr_column) + m_ScoreAssignment.ScoreGapBoundaryOpen( 1);
      // if cell (0, 0) is assigned: gap at cell (1, 0) is enclosed, gap_score is score of (0, 0) + opening an enclosed gap
      if( !HORIZONTAL_GAP_MATRIX( 0, 0) && !VERTICAL_GAP_MATRIX( 0, 0))
      {
        gap_score = SCORE_MATRIX( 0, 0) + m_ScoreAssignment.ScoreGapEnclosedOpen( 1);
      }
      else // cell (0, 0) is not assigned: gap at cell (1, 0) is boundary extension gap
      {
        gap_score = SCORE_MATRIX( 0, 0) + m_ScoreAssignment.ScoreGapBoundaryExtend( 1);
      }

      if( assign_score >= gap_score) // assign (1, 0) together
      {
        SCORE_MATRIX( 1, 0) = assign_score;
        HORIZONTAL_GAP_MATRIX( 1, 0) = VERTICAL_GAP_MATRIX( 1, 0) = 0;
        ROW_MATCH_MAX_SCORE( 1) = assign_score;
        ROW_MATCH_MAX_SCORE_POSITION( 1) = 0;
        column_origin_found = true;
        BCL_MessageDbg( "Cell ( 1, 0) origin found");
      }
      else if( assign_score < gap_score) // assign score is not larger than gap score, make gap
      {
        SCORE_MATRIX( 1, 0) = gap_score;
        HORIZONTAL_GAP_MATRIX( 1, 0) = 0;
        VERTICAL_GAP_MATRIX( 1, 0) = 1;
      }
      else
      {
        BCL_Exit
        (
          "Unable to determine (1, 0) while trying to align\n"
            + util::Format()( *itr_row) + "\nand\n" + util::Format()( *itr_column)
            + "\nassign_score=" + util::Format()( assign_score) + ", gap_score=" + util::Format()( gap_score),
          -1
        );
      }

      // cell (0, 1)
      // step column forward, set row to beginning
      ++itr_column;
      itr_row = ITR_ROW;
      assign_score = m_ScoreAssignment( **itr_row, **itr_column) + m_ScoreAssignment.ScoreGapBoundaryOpen( 1);
      // if cell (0, 0) is assigned: gap at cell (0, 1) is enclosed, gap_score is score of (0, 0) + opening an enclosed gap
      if( !HORIZONTAL_GAP_MATRIX( 0, 0) && !VERTICAL_GAP_MATRIX( 0, 0))
      {
        gap_score = SCORE_MATRIX( 0, 0) + m_ScoreAssignment.ScoreGapEnclosedOpen( 1);
      }
      else // cell (0, 0) is not assigned: gap at cell (0, 1) is boundary extension gap
      {
        gap_score = SCORE_MATRIX( 0, 0) + m_ScoreAssignment.ScoreGapBoundaryExtend( 1);
      }

      if( assign_score >= gap_score)
      {
        SCORE_MATRIX( 0, 1) = assign_score;
        HORIZONTAL_GAP_MATRIX( 0, 1) = VERTICAL_GAP_MATRIX( 0, 1) = 0;
        COLUMN_MATCH_MAX_SCORE( 1) = assign_score;
        COLUMN_MATCH_MAX_SCORE_POSITION( 1) = 0;
        row_origin_found = true;
        BCL_MessageDbg( "( 0, 1) origin found");
      }
      else if( assign_score < gap_score) // assign score is not larger than gap score
      {
        SCORE_MATRIX( 0, 1) = gap_score;
        HORIZONTAL_GAP_MATRIX( 0, 1) = 1;
        VERTICAL_GAP_MATRIX( 0, 1) = 0;
      }
      else
      {
        BCL_Exit
        (
          "Unable to determine (0, 1) while trying to align\n"
            + util::Format()( *itr_row) + "\nand\n" + util::Format()( *itr_column)
            + "\nassign_score=" + util::Format()( assign_score) + ", gap_score=" + util::Format()( gap_score),
          -1
        );
      }

      // remaining boundary gaps in vertical direction
      // set itr_row to third row, itr_column to first column
      ++++itr_row;
      itr_column = ITR_COLUMN;
      for( size_t row( 2); row != VERTICAL_DIRECTION_MAX_LENGTH; ++row, ++itr_row)
      {
        double path_column_max_score( util::GetUndefined< double>()); // score of coming from max score in column 0
        double assign_score
        (
          m_ScoreAssignment( **itr_row, **itr_column)
          + m_ScoreAssignment.ScoreGapBoundaryOpen( 1)
          + m_ScoreAssignment.ScoreGapBoundaryExtend( row - 1)
        );
        double gap_score( 0.0);
        // Better off arriving from cell above (which would be a open or extend enclosed gap penalty, depending
        // on if it is an origin (e.g. has gap matrices of 0 and 0))? Or begin own origin which means must subtract
        // appropriate boundary gap extension and open penalties?
        if( !HORIZONTAL_GAP_MATRIX( row - 1, 0) && !VERTICAL_GAP_MATRIX( row - 1, 0)) // cell above is an origin
        {
          gap_score = SCORE_MATRIX( row - 1, 0) + m_ScoreAssignment.ScoreGapEnclosedOpen( 1);
        }
        else if( !HORIZONTAL_GAP_MATRIX( row - 1, 0) && VERTICAL_GAP_MATRIX( row - 1, 0)) // cell above is no origin
        {
          // gap score is the score of the cell above plus an enclosed gap extension
          gap_score = SCORE_MATRIX( row - 1, 0) + m_ScoreAssignment.ScoreGapEnclosedExtend( 1);
        }
        else
        {
          BCL_Exit( "Cannot determine how to score vertical boundary gap for row " + util::Format()( row), -1);
        }

        // if there is a max score defined for column zero then must calculate the score of coming from that cell
        if( util::IsDefined( COLUMN_MATCH_MAX_SCORE( 0)))
        {
          // calculate score of coming from cell with the max score within the column 0
          path_column_max_score = COLUMN_MATCH_MAX_SCORE( 0) + m_ScoreAssignment.ScoreGapEnclosedOpen( 1)
            + m_ScoreAssignment.ScoreGapEnclosedExtend( row - COLUMN_MATCH_MAX_SCORE_POSITION( 0) - 1);
        }
        // true if path_column_max_score is defined, then need to take it into account when determining path
        if( util::IsDefined( path_column_max_score))
        {
          // true if assign_score is largest
          if( assign_score >= gap_score && assign_score >= path_column_max_score)
          {
            SCORE_MATRIX( row, 0) = assign_score;
            HORIZONTAL_GAP_MATRIX( row, 0) = 0;
            VERTICAL_GAP_MATRIX( row, 0) = 0;
            ROW_MATCH_MAX_SCORE( row) = assign_score;
            ROW_MATCH_MAX_SCORE_POSITION( row) = 0;

            if // true if assign_score is a new maximum score in column 0
            (
              ( SCORE_MATRIX( row, 0) + m_ScoreAssignment.ScoreGapEnclosedOpen( 1)) >=
              ( path_column_max_score + m_ScoreAssignment.ScoreGapEnclosedExtend( 1))
            )
            {
              COLUMN_MATCH_MAX_SCORE( 0) = assign_score;
              COLUMN_MATCH_MAX_SCORE_POSITION( 0) = row;
            }
            column_origin_found = true;
            BCL_MessageDbg( "Cell (" + util::Format()( row) + ", 0) origin found");
          }
          // true if gap_score is largest
          else if( gap_score >= assign_score && gap_score >= path_column_max_score)
          {
            SCORE_MATRIX( row, 0) = gap_score;
            HORIZONTAL_GAP_MATRIX( row, 0) = 0;
            VERTICAL_GAP_MATRIX( row, 0) = VERTICAL_GAP_MATRIX( row - 1, 0) + 1;
          }
          // true if path_column_max_score is largest
          else if( path_column_max_score >= assign_score && path_column_max_score >= gap_score)
          {
            SCORE_MATRIX( row, 0) = path_column_max_score;
            HORIZONTAL_GAP_MATRIX( row, 0) = 0;
            VERTICAL_GAP_MATRIX( row, 0) = row - COLUMN_MATCH_MAX_SCORE_POSITION( 0);
          }
          else
          {
            BCL_Exit
            (
              "Unable to determine (" + util::Format()( row) + ", 0) while trying to align\n"
                + util::Format()( *itr_row) + "\nand\n" + util::Format()( *itr_column)
                + "\nassign_score=" + util::Format()( assign_score) + ", gap_score=" + util::Format()( gap_score)
                + "path_column_max_score=" + util::Format()( path_column_max_score),
              -1
            );
          }
        }
        else // path_column_max_score is not defined, so can not take it into account when determining path
        {
          BCL_MessageDbg( "( " + util::Format()( row) + ", 0): path_column_max_score NOT used");

          // true if assign_score is largest
          if( assign_score >= gap_score)
          {
            SCORE_MATRIX( row, 0) = assign_score;
            HORIZONTAL_GAP_MATRIX( row, 0) = VERTICAL_GAP_MATRIX( row, 0) = 0;
            ROW_MATCH_MAX_SCORE( row) = assign_score;
            ROW_MATCH_MAX_SCORE_POSITION( row) = 0;
            // since max score for column zero was not previously defined, can automatically set it here
            COLUMN_MATCH_MAX_SCORE( 0) = assign_score;
            COLUMN_MATCH_MAX_SCORE_POSITION( 0) = row;
            column_origin_found = true;
            BCL_MessageDbg( "Cell (" + util::Format()( row) + ", 0) origin found");
          }
          // true if gap_score is largest
          else if( gap_score >= assign_score)
          {
            SCORE_MATRIX( row, 0) = gap_score;
            HORIZONTAL_GAP_MATRIX( row, 0) = 0;
            VERTICAL_GAP_MATRIX( row, 0) = VERTICAL_GAP_MATRIX( row - 1, 0) + 1;
          }
          else
          {
            BCL_Exit
            (
              "Unable to determine (" + util::Format()( row) + ", 0) while trying to align\n"
                + util::Format()( *itr_row) + "\nand\n" + util::Format()( *itr_column)
                + "\nassign_score=" + util::Format()( assign_score) + ", gap_score=" + util::Format()( gap_score)
                + "path_column_max_score=" + util::Format()( path_column_max_score),
              -1
            );
          }
        }
      }

      // boundary gaps in horizontal direction
      ++++itr_column;
      itr_row = ITR_ROW;
      for( size_t column( 2); column != HORIZONTAL_DIRECTION_MAX_LENGTH; ++column, ++itr_column)
      {
        double path_row_max_score( util::GetUndefined< double>()); // score of coming from max score in row 0
        double assign_score
        (
          m_ScoreAssignment( ( **itr_row), ( **itr_column))
          + m_ScoreAssignment.ScoreGapBoundaryOpen( 1)
          + m_ScoreAssignment.ScoreGapBoundaryExtend( column - 1)
        );
        double gap_score( 0.0);
        // Better off arriving from cell next door (which would be a open or extend enclosed gap penalty, depending
        // on if it is an origin (e.g. as gap matrices of 0 and 0))? Or begin own origin which means must subtract
        // appropriate boundary gap extension and open penalties?
        if( !HORIZONTAL_GAP_MATRIX( 0, column - 1) && !VERTICAL_GAP_MATRIX( 0, column - 1)) // neighbor is an origin
        {
          gap_score = SCORE_MATRIX( 0, column - 1) + m_ScoreAssignment.ScoreGapEnclosedOpen( 1);
        }
        // neighbor is not an origin
        else if( HORIZONTAL_GAP_MATRIX( 0, column - 1) && !VERTICAL_GAP_MATRIX( 0, column - 1))
        {
          gap_score = SCORE_MATRIX( 0, column - 1) + m_ScoreAssignment.ScoreGapEnclosedExtend( 1);
        }
        else
        {
          BCL_Exit
          (
            "Cannot determine how to score horizontal boundary gap for column" + util::Format()( column)
            + "previous horizontal gap is " + util::Format()( HORIZONTAL_GAP_MATRIX( 0, column - 1))
            + " and previous vertical gap is " + util::Format()( VERTICAL_GAP_MATRIX( 0, column - 1)),
            -1
          );
        }

        // if there is a max score defined for row zero then must calculate the score of coming from that cell
        if( util::IsDefined( ROW_MATCH_MAX_SCORE( 0)))
        {
          // calculate score of coming from cell with the max score within the row 0
          path_row_max_score = ROW_MATCH_MAX_SCORE( 0) + m_ScoreAssignment.ScoreGapEnclosedOpen( 1)
            + m_ScoreAssignment.ScoreGapEnclosedExtend( column - ROW_MATCH_MAX_SCORE_POSITION( 0) - 1);
        }
        // true if path_column_max_score is defined, then need to take it into account when determining path
        if( util::IsDefined( path_row_max_score))
        {
          // true if assign_score is largest
          if( assign_score >= gap_score && assign_score >= path_row_max_score)
          {
            SCORE_MATRIX( 0, column) = assign_score;
            HORIZONTAL_GAP_MATRIX( 0, column) = VERTICAL_GAP_MATRIX( 0, column) = 0;
            COLUMN_MATCH_MAX_SCORE( column) = assign_score;
            COLUMN_MATCH_MAX_SCORE_POSITION( column) = 0;

            if // true if assign_score is a new maximum score in row 0
            (
              ( SCORE_MATRIX( 0, column) + m_ScoreAssignment.ScoreGapEnclosedOpen( 1)) >=
              ( path_row_max_score + m_ScoreAssignment.ScoreGapEnclosedExtend( 1))
            )
            {
              ROW_MATCH_MAX_SCORE( 0) = assign_score;
              ROW_MATCH_MAX_SCORE_POSITION( 0) = column;
            }
            row_origin_found = true;
            BCL_MessageDbg( "Cell (0, " + util::Format()( column) + ") origin found");
          }
          // true if gap_score is largest
          else if( gap_score >= assign_score && gap_score >= path_row_max_score)
          {
            SCORE_MATRIX( 0, column) = gap_score;
            HORIZONTAL_GAP_MATRIX( 0, column) = HORIZONTAL_GAP_MATRIX( 0, column - 1) + 1;
            VERTICAL_GAP_MATRIX( 0, column) = 0;
          }
          // true if path_row_max_score is largest
          else if( path_row_max_score >= assign_score && path_row_max_score >= gap_score)
          {
            SCORE_MATRIX( 0, column) = path_row_max_score;
            HORIZONTAL_GAP_MATRIX( 0, column) = column - ROW_MATCH_MAX_SCORE_POSITION( 0);
            VERTICAL_GAP_MATRIX( 0, column) = 0;
          }
          else
          {
            BCL_Exit
            (
              "Unable to determine (0, " + util::Format()( column) + ") while trying to align\n"
                + util::Format()( *itr_row) + "\nand\n" + util::Format()( *itr_column)
                + "\nassign_score=" + util::Format()( assign_score) + ", gap_score=" + util::Format()( gap_score)
                + "path_row_max_score=" + util::Format()( path_row_max_score),
              -1
            );
          }
        }
        else // path_row_max_score is not defined so don't need to take it into account
        {
          BCL_MessageDbg( "( 0, " + util::Format()( column) + "): path_row_max_score NOT used");
          if( assign_score >= gap_score)
          {
            SCORE_MATRIX( 0, column) = assign_score;
            HORIZONTAL_GAP_MATRIX( 0, column) = VERTICAL_GAP_MATRIX( 0, column) = 0;
            COLUMN_MATCH_MAX_SCORE( column) = assign_score;
            COLUMN_MATCH_MAX_SCORE_POSITION( column) = 0;
            // since max score for row zero was not previously defined, can automatically set it here
            ROW_MATCH_MAX_SCORE( 0) = assign_score;
            ROW_MATCH_MAX_SCORE_POSITION( 0) = column;
            row_origin_found = true;
            BCL_MessageDbg( "Cell ( 0, " + util::Format()( column) + ") origin found");
          }
          else if( assign_score < gap_score)
          {
            SCORE_MATRIX( 0, column) = gap_score;
            HORIZONTAL_GAP_MATRIX( 0, column) = HORIZONTAL_GAP_MATRIX( 0, column - 1) + 1;
            VERTICAL_GAP_MATRIX( 0, column) = 0;
          }
          else
          {
            BCL_Exit
            (
              "Unable to determine (0, " + util::Format()( column) + ") while trying to align\n"
                + util::Format()( *itr_row) + "\nand\n" + util::Format()( *itr_column)
                + "\nassign_score=" + util::Format()( assign_score) + ", gap_score=" + util::Format()( gap_score)
                + "path_row_max_score=" + util::Format()( path_row_max_score),
              -1
            );
          }
        }
      }

      BCL_Assert( column_origin_found, "No origin was found in column 0 (the boundary column)");
      BCL_Assert( row_origin_found, "No origin was found in row 0 (the boundary row)");
    } // FillInitialBoundaryScoresAndGaps

    //! @brief DeterminePathAndScore calculates best score for a cell and path it must come from to get that score
    //! @param ROW current row in matrix
    //! @param COLUMN current column in matrix
    //! @param ITR_ROW current element of VERTICAL_ALIGNMENT
    //! @param ITR_COLUMN current element of HORIZONTAL_ALIGNMENT
    //! @param SCORE_MATRIX the matrix which will hold all the scores calculated
    //! @param HORIZONTAL_GAP_MATRIX matrix to hold the horizontal gaps which are associated with each cell
    //! @param VERTICAL_GAP_MATRIX matrix to hold the vertical gaps which are associated with each cell
    //! @param ROW_MATCH_MAX_SCORE Vector which holds best score of any matched cell in each row
    //! @param ROW_MATCH_MAX_SCORE_POSITION Vector holds the position in each row that the match_max_score occurs
    //! @param COLUMN_MATCH_MAX_SCORE Vector which holds best score of any matched cell in each column
    //! @param COLUMN_MATCH_MAX_SCORE_POSITION  Vector holds position in each column that the match_max_score occurs
    template< typename t_Member>
    void AlignerDynamicProgramming< t_Member>::DeterminePathAndScore
    (
      const size_t &ROW,
      const size_t &COLUMN,
      const typename AlignmentLeaf< t_Member>::const_iterator &ITR_ROW,
      const typename AlignmentLeaf< t_Member>::const_iterator &ITR_COLUMN,
      linal::Matrix< double> &SCORE_MATRIX,
      linal::MatrixInterface< size_t> &HORIZONTAL_GAP_MATRIX,
      linal::MatrixInterface< size_t> &VERTICAL_GAP_MATRIX,
      storage::Vector< double> &ROW_MATCH_MAX_SCORE,
      storage::Vector< size_t> &ROW_MATCH_MAX_SCORE_POSITION,
      storage::Vector< double> &COLUMN_MATCH_MAX_SCORE,
      storage::Vector< size_t> &COLUMN_MATCH_MAX_SCORE_POSITION
    ) const
    {
      // diagonal path score
      double path_diagonal_score( SCORE_MATRIX( ROW - 1, COLUMN - 1) + m_ScoreAssignment( **ITR_ROW, **ITR_COLUMN));
      double path_horizontal_score; // horizontal path score
      double path_vertical_score; // vertical path score
      double path_row_max_score( util::GetUndefinedDouble()); // score of coming from max score in row
      double path_column_max_score( util::GetUndefinedDouble()); // score of coming from max score in column

      // true if score is being calculated for either the right or bottom edge of matrix
      // if true then need to use boundary gaps for score calculations
      if( ROW == SCORE_MATRIX.GetNumberRows() - 1 || COLUMN == SCORE_MATRIX.GetNumberCols() - 1)
      {
        // true if there is a maximum score defined for ROW; doensn't make sense to calculate a score from a
        // something that is undefined
        if( util::IsDefined( ROW_MATCH_MAX_SCORE( ROW)))
        {
          path_row_max_score = ROW_MATCH_MAX_SCORE( ROW) + m_ScoreAssignment.ScoreGapBoundaryOpen( 1)
            + m_ScoreAssignment.ScoreGapBoundaryExtend( COLUMN - ROW_MATCH_MAX_SCORE_POSITION( ROW) - 1);
        }
        // true if there is a maximum score defined for COLUMN; doensn't make sense to calculate a score from a
        // something that is undefined
        if( util::IsDefined( COLUMN_MATCH_MAX_SCORE( COLUMN)))
        {
          path_column_max_score = COLUMN_MATCH_MAX_SCORE( COLUMN) + m_ScoreAssignment.ScoreGapBoundaryOpen( 1)
            + m_ScoreAssignment.ScoreGapBoundaryExtend( ROW - COLUMN_MATCH_MAX_SCORE_POSITION( COLUMN) - 1);
        }
        if( HORIZONTAL_GAP_MATRIX( ROW, COLUMN - 1)) // true if horizontal neighbor cell is a horizontal gap
        {
          path_horizontal_score = SCORE_MATRIX( ROW, COLUMN - 1) + m_ScoreAssignment.ScoreGapBoundaryExtend( 1);
        }
        else // horizontal neighbor cell is not a horizontal gap
        {
          path_horizontal_score = SCORE_MATRIX( ROW, COLUMN - 1) + m_ScoreAssignment.ScoreGapBoundaryOpen( 1);
        }
        if( VERTICAL_GAP_MATRIX( ROW - 1, COLUMN)) // true if vertical neighbor cell is already a vertical gap
        {
          path_vertical_score = SCORE_MATRIX( ROW - 1, COLUMN) + m_ScoreAssignment.ScoreGapBoundaryExtend( 1);
        }
        else // vertical neighbor cell is not a gap
        {
          path_vertical_score = SCORE_MATRIX( ROW - 1, COLUMN) + m_ScoreAssignment.ScoreGapBoundaryOpen( 1);
        }
      }
      else // don't use boundary gaps
      {
        // true if there is a maximum score defined for ROW; doensn't make sense to calculate a score from a
        // something that is undefined
        if( util::IsDefined( ROW_MATCH_MAX_SCORE( ROW)))
        {
          path_row_max_score = ROW_MATCH_MAX_SCORE( ROW) + m_ScoreAssignment.ScoreGapEnclosedOpen( 1)
            + m_ScoreAssignment.ScoreGapEnclosedExtend( COLUMN - ROW_MATCH_MAX_SCORE_POSITION( ROW) - 1);
        }
        // true if there is a maximum score defined for COLUMN; doensn't make sense to calculate a score from a
        // something that is undefined
        if( util::IsDefined( COLUMN_MATCH_MAX_SCORE( COLUMN)))
        {
          path_column_max_score = COLUMN_MATCH_MAX_SCORE( COLUMN) + m_ScoreAssignment.ScoreGapEnclosedOpen( 1)
            + m_ScoreAssignment.ScoreGapEnclosedExtend( ROW - COLUMN_MATCH_MAX_SCORE_POSITION( COLUMN) - 1);
        }
        if( HORIZONTAL_GAP_MATRIX( ROW, COLUMN - 1)) // true if horizontal neighbor cell is a horizontal gap
        {
          path_horizontal_score = SCORE_MATRIX( ROW, COLUMN - 1) + m_ScoreAssignment.ScoreGapEnclosedExtend( 1);
        }
        else // horizontal neighbor cell is not a horizontal gap
        {
          path_horizontal_score = SCORE_MATRIX( ROW, COLUMN - 1) + m_ScoreAssignment.ScoreGapEnclosedOpen( 1);
        }
        if( VERTICAL_GAP_MATRIX( ROW - 1, COLUMN)) // true if vertical neighbor cell is already a vertical gap
        {
          path_vertical_score = SCORE_MATRIX( ROW - 1, COLUMN) + m_ScoreAssignment.ScoreGapEnclosedExtend( 1);
        }
        else // vertical neighbor cell is not a gap
        {
          path_vertical_score = SCORE_MATRIX( ROW - 1, COLUMN) + m_ScoreAssignment.ScoreGapEnclosedOpen( 1);
        }
      }

      // make sure path_row_max_score was calculated or not calculated in correct fashion
      BCL_Assert
      (
        util::IsDefined( ROW_MATCH_MAX_SCORE( ROW)) == util::IsDefined( path_row_max_score),
        "ROW_MATCH_MAX_SCORE(ROW)=" + util::Format()( util::IsDefined( ROW_MATCH_MAX_SCORE( ROW)))
          + " but path_row_max_score=" + util::Format()( util::IsDefined( path_row_max_score))
          + " for row " + util::Format()( ROW)
      );
      // make sure path_column_max_score was calculated or not calculated in correct fashion
      BCL_Assert
      (
        util::IsDefined( COLUMN_MATCH_MAX_SCORE( COLUMN)) == util::IsDefined( path_column_max_score),
        "COLUMN_MATCH_MAX_SCORE(COLUMN)=" + util::Format()( util::IsDefined( COLUMN_MATCH_MAX_SCORE( COLUMN)))
          + " but path_column_max_score=" + util::Format()( util::IsDefined( path_column_max_score))
          + " for column" + util::Format()( COLUMN)
      );

      // set the path which gives the maximum score and sets the score and gap matrices
      SetPathAndScore
      (
        path_diagonal_score,
        path_horizontal_score,
        path_vertical_score,
        path_row_max_score,
        path_column_max_score,
        ROW,
        COLUMN,
        SCORE_MATRIX,
        HORIZONTAL_GAP_MATRIX,
        VERTICAL_GAP_MATRIX,
        ROW_MATCH_MAX_SCORE,
        ROW_MATCH_MAX_SCORE_POSITION,
        COLUMN_MATCH_MAX_SCORE,
        COLUMN_MATCH_MAX_SCORE_POSITION
      );
    } // DeterminePathAndScore

    //! @brief SetPathAndScore decides where a cell is coming from based on the possible pathway scores
    //! @param DIAGONAL_PATH_SCORE score cell would acquire coming from diagonal path
    //! @param HORIZONTAL_PATH_SCORE score cell would acquire coming from horizontal neighbor
    //! @param VERTICAL_PATH_SCORE score cell would acquire coming from vertical neighbor
    //! @param PATH_ROW_MAX_SCORE score cell would acquire coming from cell with maximum score in row
    //! @param PATH_COLUMN_MAX_SCORE score cell would acquire coming from cell with maximum score in column
    //! @param ROW current row in matrix
    //! @param COLUMN current column in matrix
    //! @param SCORE_MATRIX the matrix which will hold all the scores calculated
    //! @param HORIZONTAL_GAP_MATRIX matrix to hold the horizontal gaps which are associated with each cell
    //! @param VERTICAL_GAP_MATRIX matrix to hold the vertical gaps which are associated with each cell
    //! @param ROW_MATCH_MAX_SCORE Vector which holds best score of any matched cell in each row
    //! @param ROW_MATCH_MAX_SCORE_POSITION Vector holds the position in each row that the match_max_score occurs
    //! @param COLUMN_MATCH_MAX_SCORE Vector which holds best score of any matched cell in each column
    //! @param COLUMN_MATCH_MAX_SCORE_POSITION  Vector holds position in each column that the match_max_score occurs
    template< typename t_Member>
    void AlignerDynamicProgramming< t_Member>::SetPathAndScore
    (
      const double &DIAGONAL_PATH_SCORE,
      const double &HORIZONTAL_PATH_SCORE,
      const double &VERTICAL_PATH_SCORE,
      const double &PATH_ROW_MAX_SCORE,
      const double &PATH_COLUMN_MAX_SCORE,
      const size_t &ROW,
      const size_t &COLUMN,
      linal::Matrix< double> &SCORE_MATRIX,
      linal::MatrixInterface< size_t> &HORIZONTAL_GAP_MATRIX,
      linal::MatrixInterface< size_t> &VERTICAL_GAP_MATRIX,
      storage::Vector< double> &ROW_MATCH_MAX_SCORE,
      storage::Vector< size_t> &ROW_MATCH_MAX_SCORE_POSITION,
      storage::Vector< double> &COLUMN_MATCH_MAX_SCORE,
      storage::Vector< size_t> &COLUMN_MATCH_MAX_SCORE_POSITION
    ) const
    {
      // create doubles for holding what will be used for the penalty of opening and extending a gap
      double score_gap_open, score_gap_extend;

      // true if need to use boundary gap scoring
      if( ROW == SCORE_MATRIX.GetNumberRows() - 1 || COLUMN == SCORE_MATRIX.GetNumberCols() - 1)
      {
        score_gap_open = m_ScoreAssignment.ScoreGapBoundaryOpen( 1);
        score_gap_extend = m_ScoreAssignment.ScoreGapBoundaryExtend( 1);
      }
      else // use enclosed gap scoring
      {
        score_gap_open = m_ScoreAssignment.ScoreGapEnclosedOpen( 1);
        score_gap_extend = m_ScoreAssignment.ScoreGapEnclosedExtend( 1);
      }

      // true if both ROW_MATCH_MAX_SCORE( ROW) and COLUMN_MATCH_MAX_SCORE( COLUMN) are defined
      if( util::IsDefined( ROW_MATCH_MAX_SCORE( ROW)) && util::IsDefined( COLUMN_MATCH_MAX_SCORE( COLUMN)))
      {
        if // is DIAGONAL_PATH_SCORE largest?
        (
          ( DIAGONAL_PATH_SCORE >= HORIZONTAL_PATH_SCORE)
          && ( DIAGONAL_PATH_SCORE >= VERTICAL_PATH_SCORE)
          && ( DIAGONAL_PATH_SCORE >= PATH_ROW_MAX_SCORE)
          && ( DIAGONAL_PATH_SCORE >= PATH_COLUMN_MAX_SCORE)
        )
        {
          SCORE_MATRIX( ROW, COLUMN) = DIAGONAL_PATH_SCORE;
          HORIZONTAL_GAP_MATRIX( ROW, COLUMN) = VERTICAL_GAP_MATRIX( ROW, COLUMN) = 0;

          // true if DIAGONAL_PATH_SCORE is a new maximum score in ROW
          if( ( SCORE_MATRIX( ROW, COLUMN) + score_gap_open) >= ( PATH_ROW_MAX_SCORE + score_gap_extend))
          {
            ROW_MATCH_MAX_SCORE( ROW) = DIAGONAL_PATH_SCORE;
            ROW_MATCH_MAX_SCORE_POSITION( ROW) = COLUMN;
          }
          // true if DIAGONAL_PATH_SCORE is a new maximum score in COLUMN
          if( ( SCORE_MATRIX( ROW, COLUMN) + score_gap_open) >= ( PATH_COLUMN_MAX_SCORE + score_gap_extend))
          {
            COLUMN_MATCH_MAX_SCORE( COLUMN) = DIAGONAL_PATH_SCORE;
            COLUMN_MATCH_MAX_SCORE_POSITION( COLUMN) = ROW;
          }
        }
        else if // is HORIZONTAL_PATH_SCORE largest?
        (
          ( HORIZONTAL_PATH_SCORE >= DIAGONAL_PATH_SCORE)
          && ( HORIZONTAL_PATH_SCORE >= VERTICAL_PATH_SCORE)
          && ( HORIZONTAL_PATH_SCORE >= PATH_ROW_MAX_SCORE)
          && ( HORIZONTAL_PATH_SCORE >= PATH_COLUMN_MAX_SCORE)
        )
        {
          SCORE_MATRIX( ROW, COLUMN) = HORIZONTAL_PATH_SCORE;
          HORIZONTAL_GAP_MATRIX( ROW, COLUMN) = 1;
          VERTICAL_GAP_MATRIX( ROW, COLUMN) = 0;
        }
        else if // is VERTICAL_PATH_SCORE largest?
        (
          ( VERTICAL_PATH_SCORE >= DIAGONAL_PATH_SCORE)
          && ( VERTICAL_PATH_SCORE >= HORIZONTAL_PATH_SCORE)
          && ( VERTICAL_PATH_SCORE >= PATH_ROW_MAX_SCORE)
          && ( VERTICAL_PATH_SCORE >= PATH_COLUMN_MAX_SCORE)
        )
        {
          SCORE_MATRIX( ROW, COLUMN) = VERTICAL_PATH_SCORE;
          HORIZONTAL_GAP_MATRIX( ROW, COLUMN) = 0;
          VERTICAL_GAP_MATRIX( ROW, COLUMN) = 1;
        }
        else if // is PATH_ROW_MAX_SCORE largest?
        (
          ( PATH_ROW_MAX_SCORE >= DIAGONAL_PATH_SCORE)
          && ( PATH_ROW_MAX_SCORE >= HORIZONTAL_PATH_SCORE)
          && ( PATH_ROW_MAX_SCORE >= VERTICAL_PATH_SCORE)
          && ( PATH_ROW_MAX_SCORE >= PATH_COLUMN_MAX_SCORE)
        )
        {
          SCORE_MATRIX( ROW, COLUMN) = PATH_ROW_MAX_SCORE;
          HORIZONTAL_GAP_MATRIX( ROW, COLUMN) = COLUMN - ROW_MATCH_MAX_SCORE_POSITION( ROW);
          VERTICAL_GAP_MATRIX( ROW, COLUMN) = 0;
        }
        else if // is PATH_COLUMN_MAX_SCORE largest?
        (
          ( PATH_COLUMN_MAX_SCORE >= DIAGONAL_PATH_SCORE)
          && ( PATH_COLUMN_MAX_SCORE >= HORIZONTAL_PATH_SCORE)
          && ( PATH_COLUMN_MAX_SCORE >= VERTICAL_PATH_SCORE)
          && ( PATH_COLUMN_MAX_SCORE >= PATH_ROW_MAX_SCORE)
        )
        {
          SCORE_MATRIX( ROW, COLUMN) = PATH_COLUMN_MAX_SCORE;
          HORIZONTAL_GAP_MATRIX( ROW, COLUMN) = 0;
          VERTICAL_GAP_MATRIX( ROW, COLUMN) = ROW - COLUMN_MATCH_MAX_SCORE_POSITION( COLUMN);
        }
        else
        {
          BCL_Exit( "A.) Unable to determine path while filling matrix.", -1);
        }
      }
      // true if ROW_MATCH_MAX_SCORE( ROW) is defined but not COLUMN_MATCH_MAX_SCORE( COLUMN)
      else if( util::IsDefined( ROW_MATCH_MAX_SCORE( ROW)) && !util::IsDefined( COLUMN_MATCH_MAX_SCORE( COLUMN)))
      {
        if // is DIAGONAL_PATH_SCORE largest?
        (
          ( DIAGONAL_PATH_SCORE >= HORIZONTAL_PATH_SCORE)
          && ( DIAGONAL_PATH_SCORE >= VERTICAL_PATH_SCORE)
          && ( DIAGONAL_PATH_SCORE >= PATH_ROW_MAX_SCORE)
        )
        {
          SCORE_MATRIX( ROW, COLUMN) = DIAGONAL_PATH_SCORE;
          HORIZONTAL_GAP_MATRIX( ROW, COLUMN) = VERTICAL_GAP_MATRIX( ROW, COLUMN) = 0;
          COLUMN_MATCH_MAX_SCORE( COLUMN) = DIAGONAL_PATH_SCORE;
          COLUMN_MATCH_MAX_SCORE_POSITION( COLUMN) = ROW;
          // true if DIAGONAL_PATH_SCORE is a new maximum score in ROW
          if( ( SCORE_MATRIX( ROW, COLUMN) + score_gap_open) >= ( PATH_ROW_MAX_SCORE + score_gap_extend))
          {
            ROW_MATCH_MAX_SCORE( ROW) = DIAGONAL_PATH_SCORE;
            ROW_MATCH_MAX_SCORE_POSITION( ROW) = COLUMN;
          }
        }
        else if // is HORIZONTAL_PATH_SCORE largest?
        (
          ( HORIZONTAL_PATH_SCORE >= DIAGONAL_PATH_SCORE)
          && ( HORIZONTAL_PATH_SCORE >= VERTICAL_PATH_SCORE)
          && ( HORIZONTAL_PATH_SCORE >= PATH_ROW_MAX_SCORE)
        )
        {
          SCORE_MATRIX( ROW, COLUMN) = HORIZONTAL_PATH_SCORE;
          HORIZONTAL_GAP_MATRIX( ROW, COLUMN) = 1;
          VERTICAL_GAP_MATRIX( ROW, COLUMN) = 0;
        }
        else if // is VERTICAL_PATH_SCORE largest?
        (
          ( VERTICAL_PATH_SCORE >= DIAGONAL_PATH_SCORE)
          && ( VERTICAL_PATH_SCORE >= HORIZONTAL_PATH_SCORE)
          && ( VERTICAL_PATH_SCORE >= PATH_ROW_MAX_SCORE)
        )
        {
          SCORE_MATRIX( ROW, COLUMN) = VERTICAL_PATH_SCORE;
          HORIZONTAL_GAP_MATRIX( ROW, COLUMN) = 0;
          VERTICAL_GAP_MATRIX( ROW, COLUMN) = 1;
        }
        else if // is PATH_ROW_MAX_SCORE largest?
        (
          ( PATH_ROW_MAX_SCORE >= DIAGONAL_PATH_SCORE)
          && ( PATH_ROW_MAX_SCORE >= HORIZONTAL_PATH_SCORE)
          && ( PATH_ROW_MAX_SCORE >= VERTICAL_PATH_SCORE)
        )
        {
          SCORE_MATRIX( ROW, COLUMN) = PATH_ROW_MAX_SCORE;
          HORIZONTAL_GAP_MATRIX( ROW, COLUMN) = COLUMN - ROW_MATCH_MAX_SCORE_POSITION( ROW);
          VERTICAL_GAP_MATRIX( ROW, COLUMN) = 0;
        }
        else
        {
          BCL_Exit( "B.) Unable to determine path while filling matrix.", -1);
        }
      }
      // true if ROW_MATCH_MAX_SCORE( ROW) is undefined and COLUMN_MATCH_MAX_SCORE( COLUMN) is defined
      else if( !util::IsDefined( ROW_MATCH_MAX_SCORE( ROW)) && util::IsDefined( COLUMN_MATCH_MAX_SCORE( COLUMN)))
      {
        if // is DIAGONAL_PATH_SCORE largest?
        (
          ( DIAGONAL_PATH_SCORE >= HORIZONTAL_PATH_SCORE)
          && ( DIAGONAL_PATH_SCORE >= VERTICAL_PATH_SCORE)
          && ( DIAGONAL_PATH_SCORE >= PATH_COLUMN_MAX_SCORE)
        )
        {
          SCORE_MATRIX( ROW, COLUMN) = DIAGONAL_PATH_SCORE;
          HORIZONTAL_GAP_MATRIX( ROW, COLUMN) = VERTICAL_GAP_MATRIX( ROW, COLUMN) = 0;
          ROW_MATCH_MAX_SCORE( ROW) = DIAGONAL_PATH_SCORE;
          ROW_MATCH_MAX_SCORE_POSITION( ROW) = COLUMN;

          // true if DIAGONAL_PATH_SCORE is a new maximum score in COLUMN
          if( ( SCORE_MATRIX( ROW, COLUMN) + score_gap_open) >= ( PATH_COLUMN_MAX_SCORE + score_gap_extend))
          {
            COLUMN_MATCH_MAX_SCORE( COLUMN) = DIAGONAL_PATH_SCORE;
            COLUMN_MATCH_MAX_SCORE_POSITION( COLUMN) = ROW;
          }
        }
        else if // is HORIZONTAL_PATH_SCORE largest?
        (
          ( HORIZONTAL_PATH_SCORE >= DIAGONAL_PATH_SCORE)
          && ( HORIZONTAL_PATH_SCORE >= VERTICAL_PATH_SCORE)
          && ( HORIZONTAL_PATH_SCORE >= PATH_COLUMN_MAX_SCORE)
        )
        {
          SCORE_MATRIX( ROW, COLUMN) = HORIZONTAL_PATH_SCORE;
          HORIZONTAL_GAP_MATRIX( ROW, COLUMN) = 1;
          VERTICAL_GAP_MATRIX( ROW, COLUMN) = 0;
        }
        else if // is VERTICAL_PATH_SCORE largest?
        (
          ( VERTICAL_PATH_SCORE >= DIAGONAL_PATH_SCORE)
          && ( VERTICAL_PATH_SCORE >= HORIZONTAL_PATH_SCORE)
          && ( VERTICAL_PATH_SCORE >= PATH_COLUMN_MAX_SCORE)
        )
        {
          SCORE_MATRIX( ROW, COLUMN) = VERTICAL_PATH_SCORE;
          HORIZONTAL_GAP_MATRIX( ROW, COLUMN) = 0;
          VERTICAL_GAP_MATRIX( ROW, COLUMN) = 1;
        }
        else if // is PATH_COLUMN_MAX_SCORE largest?
        (
          ( PATH_COLUMN_MAX_SCORE >= DIAGONAL_PATH_SCORE)
          && ( PATH_COLUMN_MAX_SCORE >= HORIZONTAL_PATH_SCORE)
          && ( PATH_COLUMN_MAX_SCORE >= VERTICAL_PATH_SCORE)
        )
        {
          SCORE_MATRIX( ROW, COLUMN) = PATH_COLUMN_MAX_SCORE;
          HORIZONTAL_GAP_MATRIX( ROW, COLUMN) = 0;
          VERTICAL_GAP_MATRIX( ROW, COLUMN) = ROW - COLUMN_MATCH_MAX_SCORE_POSITION( COLUMN);
        }
        else
        {
          BCL_Exit( "C.) Unable to determine path while filling matrix.", -1);
        }
      }
      // true if both ROW_MATCH_MAX_SCORE( ROW) and COLUMN_MATCH_MAX_SCORE( COLUMN) are undefined
      else if( !util::IsDefined( ROW_MATCH_MAX_SCORE( ROW)) && !util::IsDefined( COLUMN_MATCH_MAX_SCORE( COLUMN)))
      {
        // is DIAGONAL_PATH_SCORE largest?
        if( ( DIAGONAL_PATH_SCORE >= HORIZONTAL_PATH_SCORE) && ( DIAGONAL_PATH_SCORE >= VERTICAL_PATH_SCORE))
        {
          SCORE_MATRIX( ROW, COLUMN) = DIAGONAL_PATH_SCORE;
          HORIZONTAL_GAP_MATRIX( ROW, COLUMN) = VERTICAL_GAP_MATRIX( ROW, COLUMN) = 0;
          ROW_MATCH_MAX_SCORE( ROW) = DIAGONAL_PATH_SCORE;
          ROW_MATCH_MAX_SCORE_POSITION( ROW) = COLUMN;
          COLUMN_MATCH_MAX_SCORE( COLUMN) = DIAGONAL_PATH_SCORE;
          COLUMN_MATCH_MAX_SCORE_POSITION( COLUMN) = ROW;
        }
        // is HORIZONTAL_PATH_SCORE largest?
        else if( ( HORIZONTAL_PATH_SCORE >= DIAGONAL_PATH_SCORE) && ( HORIZONTAL_PATH_SCORE >= VERTICAL_PATH_SCORE))
        {
          SCORE_MATRIX( ROW, COLUMN) = HORIZONTAL_PATH_SCORE;
          HORIZONTAL_GAP_MATRIX( ROW, COLUMN) = 1;
          VERTICAL_GAP_MATRIX( ROW, COLUMN) = 0;
        }
        // is VERTICAL_PATH_SCORE largest?
        else if( ( VERTICAL_PATH_SCORE >= DIAGONAL_PATH_SCORE) && ( VERTICAL_PATH_SCORE >= HORIZONTAL_PATH_SCORE))
        {
          SCORE_MATRIX( ROW, COLUMN) = VERTICAL_PATH_SCORE;
          HORIZONTAL_GAP_MATRIX( ROW, COLUMN) = 0;
          VERTICAL_GAP_MATRIX( ROW, COLUMN) = 1;
        }
        else
        {
          BCL_Exit( "D.) Unable to determine path while filling matrix.", -1);
        }
      }
      else
      {
        BCL_Exit( "Cannot determine defined state of row_match_max_score and column_match_max_score", -1);
      }
    } // SetPathAndScore

    //! @brief TraceBack builds the new alignment
    //! @param ALIGNMENT_VERTICAL the first alignment involved in the Alignment
    //! @param ALIGNMENT_HORIZONTAL the second alignment involved in the Alignment
    //! @param SCORE_MATRIX the matrix which holds all the scores calculated
    //! @param HORIZONTAL_GAP_MATRIX matrix holds the horizontal gaps which are associated with each cell
    //! @param VERTICAL_GAP_MATRIX matrix holds the vertical gaps which are associated with each cell
    //! @return returns the new alignment and its score in a storag::Pair
    template< typename t_Member>
    storage::Pair< AlignmentNode< t_Member>, double> AlignerDynamicProgramming< t_Member>::TraceBack
    (
      util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_VERTICAL,
      util::ShPtr< AlignmentInterface< t_Member> > &ALIGNMENT_HORIZONTAL,
      const linal::Matrix< double> &SCORE_MATRIX,
      const linal::MatrixConstInterface< size_t> &HORIZONTAL_GAP_MATRIX,
      const linal::MatrixConstInterface< size_t> &VERTICAL_GAP_MATRIX
    ) const
    {
      // new_alignment will be alignment of ALIGNMENT_VERTICAL and ALIGNMENT_HORIZONTAL
      AlignmentNode< t_Member> new_alignment( ALIGNMENT_VERTICAL, ALIGNMENT_HORIZONTAL);

      // initialize to score of alignment of ALIGNMENT_VERTICAL and ALIGNMENT_HORIZONTAL
      double score( SCORE_MATRIX( ALIGNMENT_VERTICAL->GetSize() - 1, ALIGNMENT_HORIZONTAL->GetSize() - 1));

      // create reverse_iterators to reverse begin (last element) of ALIGNMENT_VERTICAL and ALIGNMENT_HORIZONTAL
      typename AlignmentLeaf< t_Member>::const_reverse_iterator
        rev_itr_row( ALIGNMENT_VERTICAL->GetAssignments().ReverseBegin()), // reverse iterator to traverse vertically
        rev_itr_row_reverse_end( ALIGNMENT_VERTICAL->GetAssignments().ReverseEnd()), // reverse iterator at reverse end of vertical
        rev_itr_column( ALIGNMENT_HORIZONTAL->GetAssignments().ReverseBegin()), // reverse iterator to traverse horizontally
        rev_itr_column_reverse_end( ALIGNMENT_HORIZONTAL->GetAssignments().ReverseEnd()); // reverse iterator at reverse end of horizontal

      size_t row( ALIGNMENT_VERTICAL->GetSize() - 1); // position in vertical direction
      size_t column( ALIGNMENT_HORIZONTAL->GetSize() - 1); // position in horizontal direction

      // build new_alignment until beginning of ALIGNMENT_VERTICAL or ALIGNMENT_HORIZONTAL is reached
      BuildInteriorOfAlignment
      (
        row,
        column,
        *ALIGNMENT_VERTICAL,
        *ALIGNMENT_HORIZONTAL,
        rev_itr_row,
        rev_itr_column,
        new_alignment,
        SCORE_MATRIX,
        HORIZONTAL_GAP_MATRIX,
        VERTICAL_GAP_MATRIX
      );

      // row and/or column has reached zero so fill in remaining gaps/assignment
      AddGapsToBeginningOfAlignment
      (
        row,
        column,
        rev_itr_row,
        rev_itr_row_reverse_end,
        rev_itr_column,
        rev_itr_column_reverse_end,
        new_alignment,
        HORIZONTAL_GAP_MATRIX,
        VERTICAL_GAP_MATRIX
      );

      // return Pair of new_alignment and score
      return storage::Pair< AlignmentNode< t_Member>, double>( new_alignment, score);
    } // TraceBack

    //! @brief BuildInteriorOfAlignment builds new alignment until either the upper or left side of matrix is reached
    //! @param ROW index of last row in the matrix
    //! @param COLUMN index of last column in matrix
    //! @param ALIGNMENT_VERTICAL the first alignment involved in the Alignment
    //! @param ALIGNMENT_HORIZONTAL the second alignment involved in the Alignment
    //! @param REV_ROW_ITR reverse iterator on reverse begin of ALIGNMENT_VERTICAL
    //! @param REV_COLUMN_ITR reverse iterator on reverse begin of ALIGNMENT_HORIZONTAL
    //! @param ALIGNMENT new alignment which is being built
    //! @param SCORE_MATRIX the matrix which holds all the scores calculated
    //! @param HORIZONTAL_GAP_MATRIX matrix holds the horizontal gaps which are associated with each cell
    //! @param VERTICAL_GAP_MATRIX matrix holds the vertical gaps which are associated with each cell
    template< typename t_Member>
    void AlignerDynamicProgramming< t_Member>::BuildInteriorOfAlignment
    (
      size_t &ROW,
      size_t &COLUMN,
      const AlignmentInterface< t_Member> &ALIGNMENT_VERTICAL,
      const AlignmentInterface< t_Member> &ALIGNMENT_HORIZONTAL,
      typename AlignmentLeaf< t_Member>::const_reverse_iterator &REV_ROW_ITR,
      typename AlignmentLeaf< t_Member>::const_reverse_iterator &REV_COLUMN_ITR,
      AlignmentNode< t_Member> &ALIGNMENT,
      const linal::Matrix< double> &SCORE_MATRIX,
      const linal::MatrixConstInterface< size_t> &HORIZONTAL_GAP_MATRIX,
      const linal::MatrixConstInterface< size_t> &VERTICAL_GAP_MATRIX
    ) const
    {
      typename AlignmentLeaf< t_Member>::const_reverse_iterator
        rev_itr_end_alignment_vertical( --ALIGNMENT_VERTICAL.GetAssignments().ReverseEnd()), // first element in alignment
        rev_itr_end_alignment_horizontal( --ALIGNMENT_HORIZONTAL.GetAssignments().ReverseEnd()); // first element

      while( ROW > 0 && COLUMN > 0) // when false, only gaps need to be added
      {
        // make sure iterators do not reach beginning before indices do
        BCL_Assert
        (
          REV_ROW_ITR != rev_itr_end_alignment_vertical && REV_COLUMN_ITR != rev_itr_end_alignment_horizontal,
          "Iterator has reached a beginning boundary of matrix but row/column has not"
        );
        // true if no gaps led to this cell (came from diagonal path); need to add match to alignment
        if( !HORIZONTAL_GAP_MATRIX( ROW, COLUMN) && !VERTICAL_GAP_MATRIX( ROW, COLUMN))
        {
          AddAssignmentFromDiagonalPath( REV_ROW_ITR, REV_COLUMN_ITR, ALIGNMENT);
          --ROW;
          --COLUMN;
          ++REV_ROW_ITR;
          ++REV_COLUMN_ITR;
        }
        // true if horizontal gaps led to this cell; need to add horizontal gaps to alignment
        else if( HORIZONTAL_GAP_MATRIX( ROW, COLUMN) && !VERTICAL_GAP_MATRIX( ROW, COLUMN))
        {
          // add number of horizontal gaps that led to this cell as indicated by HORIZONTAL_GAP_MATRIX( ROW, COLUMN)
          size_t number_gaps( 0), number_gaps_to_fill( HORIZONTAL_GAP_MATRIX( ROW, COLUMN));
          for( ; number_gaps < number_gaps_to_fill; ++number_gaps)
          {
            AddAssignmentFromHorizontalPath( REV_ROW_ITR, REV_COLUMN_ITR, ALIGNMENT);
            --COLUMN;
            ++REV_COLUMN_ITR;
          }
        }
        // true if vertical gaps led to this cell; need to add vertical gaps to alignment
        else if( !HORIZONTAL_GAP_MATRIX( ROW, COLUMN) && VERTICAL_GAP_MATRIX( ROW, COLUMN))
        {
          // add number of vertical gaps that led to this cell as indicated by VERTICAL_GAP_MATRIX( ROW, COLUMN)
          size_t number_gaps( 0), number_gaps_to_fill( VERTICAL_GAP_MATRIX( ROW, COLUMN));
          for( ; number_gaps < number_gaps_to_fill; ++number_gaps)
          {
            AddAssignmentFromVerticalPath( REV_ROW_ITR, REV_COLUMN_ITR, ALIGNMENT);
            --ROW;
            ++REV_ROW_ITR;
          }
        }
        else
        {
          BCL_Exit( "Traceback does not know where to go", -1);
        }
      }
    } // BuildInteriorOfAlignment

    //! @brief AddGapsToBeginningOfAlignment adds remaining assignments after finished with BuildInteriorOfAlignment
    //! @param ROW current row in matrix
    //! @param COLUMN current column in matrix
    //! @param REV_ROW_ITR iterator for current position of REV_ROW_ITR after BuildInteriorOfAlignment
    //! @param REV_ROW_ITR_END iterator at reverse end of ALIGNMENT_VERTICAL
    //! @param REV_COLUMN_ITR iterator for current position of REV_COLUMN_ITR after BuildInteriorOfAlignment
    //! @param REV_COLUMN_ITR_END iterator at reverse end of ALIGNMENT_HORIZONTAL
    //! @param ALIGNMENT new alignment which is being built
    //! @param HORIZONTAL_GAP_MATRIX matrix holds the horizontal gaps which are associated with each cell
    //! @param VERTICAL_GAP_MATRIX matrix holds the vertical gaps which are associated with each cell
    template< typename t_Member>
    void AlignerDynamicProgramming< t_Member>::AddGapsToBeginningOfAlignment
    (
      size_t &ROW,
      size_t &COLUMN,
      typename AlignmentLeaf< t_Member>::const_reverse_iterator &REV_ROW_ITR,
      typename AlignmentLeaf< t_Member>::const_reverse_iterator &REV_ROW_ITR_END,
      typename AlignmentLeaf< t_Member>::const_reverse_iterator &REV_COLUMN_ITR,
      typename AlignmentLeaf< t_Member>::const_reverse_iterator &REV_COLUMN_ITR_END,
      AlignmentNode< t_Member> &ALIGNMENT,
      const linal::MatrixConstInterface< size_t> &HORIZONTAL_GAP_MATRIX,
      const linal::MatrixConstInterface< size_t> &VERTICAL_GAP_MATRIX
    ) const
    {
      BCL_MessageDbg( "Add gaps: start row=" + util::Format()( ROW) + ", column=" + util::Format()( COLUMN));

      // true if ROW and COLUMN are zero; means we are at the upper left cell of matrix
      if( !ROW && !COLUMN)
      {
        // add match if first cell is assigned together
        if( !HORIZONTAL_GAP_MATRIX( ROW, COLUMN) && !VERTICAL_GAP_MATRIX( ROW, COLUMN))
        {
          AddAssignmentFromDiagonalPath( REV_ROW_ITR, REV_COLUMN_ITR, ALIGNMENT);
        }
        // true if first cell is not assigned together
        else if( HORIZONTAL_GAP_MATRIX( ROW, COLUMN) && VERTICAL_GAP_MATRIX( ROW, COLUMN))
        {
          // assign each first element of row and column to a gap
          AddAssignmentFromVerticalPath( REV_ROW_ITR, REV_COLUMN_ITR, ALIGNMENT);
          AddAssignmentFromHorizontalPath( REV_ROW_ITR, REV_COLUMN_ITR, ALIGNMENT);
        }
        else
        {
          BCL_Exit( "Cannot add origin cell to alignment", -1);
        }
      }
      // true if ROW has value but COLUMN does not; elements of ROW will be assigned with gaps
      else if( ROW && !COLUMN)
      {
        // create bool origin_set for keeping track if first match at start of alignment has occurred
        bool origin_set( false);

        // while all elements of row have not been assigned yet
        while( REV_ROW_ITR != REV_ROW_ITR_END)
        {
          // if REV_ROW_ITR and REV_COLUMN_ITR are assigned and the origin has not been set yet
          if( !HORIZONTAL_GAP_MATRIX( ROW, COLUMN) && !VERTICAL_GAP_MATRIX( ROW, COLUMN) && !origin_set)
          {
            AddAssignmentFromDiagonalPath( REV_ROW_ITR, REV_COLUMN_ITR, ALIGNMENT);
            origin_set = true;
          }
          // if REV_ROW_ITR and REV_COLUMN_ITR are not assigned or if the orgin has been set
          // when origin has been set the rest of REV_ROW_ITRs will be assigned with gaps
          else if( ( !HORIZONTAL_GAP_MATRIX( ROW, COLUMN) && VERTICAL_GAP_MATRIX( ROW, COLUMN)) || origin_set)
          {
            AddAssignmentFromVerticalPath( REV_ROW_ITR, REV_COLUMN_ITR, ALIGNMENT);
            // add number of vertical gaps that led to this cell as indicated by VERTICAL_GAP_MATRIX( ROW, COLUMN)
            size_t number_gaps( 1), number_gaps_to_fill( VERTICAL_GAP_MATRIX( ROW, COLUMN));
            for( ; number_gaps < number_gaps_to_fill; ++number_gaps)
            {
              --ROW;
              ++REV_ROW_ITR;
              AddAssignmentFromVerticalPath( REV_ROW_ITR, REV_COLUMN_ITR, ALIGNMENT);
            }
          }
          // true if the upper left corner of the matrix has been reached and no origins occurred on the way there
          // and the first of each sequence should not be assigned together
          else if
          (
            !ROW && !COLUMN && HORIZONTAL_GAP_MATRIX( ROW, COLUMN) && VERTICAL_GAP_MATRIX( ROW, COLUMN) && !origin_set
          )
          {
            // assign each first element of row and column to a gap
            AddAssignmentFromVerticalPath( REV_ROW_ITR, REV_COLUMN_ITR, ALIGNMENT);
            AddAssignmentFromHorizontalPath( REV_ROW_ITR, REV_COLUMN_ITR, ALIGNMENT);
          }
          else
          {
            BCL_Exit
            (
              "Cannot determine if row=" + util::Format()( ROW) + " is an origin;\n"
                + "column=" + util::Format()( COLUMN) + " (should be zero).\n"
                + "horizontal_gap_matrix( " + util::Format()( ROW) + ", " + util::Format()( COLUMN) + ")="
                + util::Format()( HORIZONTAL_GAP_MATRIX( ROW, COLUMN))
                + "; vertical_gap_matrix( " + util::Format()( ROW) + ", " + util::Format()( COLUMN) + ")="
                + util::Format()( VERTICAL_GAP_MATRIX( ROW, COLUMN)),
              -1
            );
          }

          ++REV_ROW_ITR;
          --ROW;
        }
      }
      // true if COLUMN has value but ROW does not; elements of COLUMN will be assigned with gaps
      else if( !ROW && COLUMN)
      {
        // create bool origin_set for keeping track if first match at start of alignment has occured
        bool origin_set( false);

        // while all elements of column have not been assigned yet
        while( REV_COLUMN_ITR != REV_COLUMN_ITR_END)
        {
          // if REV_ROW_ITR and REV_COLUMN_ITR are assigned and the origin has not been set yet
          if( !HORIZONTAL_GAP_MATRIX( ROW, COLUMN) && !VERTICAL_GAP_MATRIX( ROW, COLUMN) && !origin_set)
          {
            AddAssignmentFromDiagonalPath( REV_ROW_ITR, REV_COLUMN_ITR, ALIGNMENT);
            origin_set = true;
          }
          // if REV_ROW_ITR and REV_COLUMN_ITR are not assigned or if the origin has been set
          // when origin has been set the rest of REV_ROW_ITRs will be assigned with gaps
          else if( ( HORIZONTAL_GAP_MATRIX( ROW, COLUMN) && !VERTICAL_GAP_MATRIX( ROW, COLUMN)) || origin_set)
          {
            AddAssignmentFromHorizontalPath( REV_ROW_ITR, REV_COLUMN_ITR, ALIGNMENT);
            // add number of horizontal gaps that led to this cell as indicated by HORIZONTAL_GAP_MATRIX( ROW, COLUMN)
            size_t number_gaps( 1), number_gaps_to_fill( HORIZONTAL_GAP_MATRIX( ROW, COLUMN));
            for( ; number_gaps < number_gaps_to_fill; ++number_gaps)
            {
              ++REV_COLUMN_ITR;
              --COLUMN;
              AddAssignmentFromHorizontalPath( REV_ROW_ITR, REV_COLUMN_ITR, ALIGNMENT);
            }
          }
          // true if the upper left corner of the matrix has been reached and no origins occurred on the way there
          // and the first of each sequence should not be assigned together
          else if
          (
            !ROW && !COLUMN && HORIZONTAL_GAP_MATRIX( ROW, COLUMN) && VERTICAL_GAP_MATRIX( ROW, COLUMN) && !origin_set
          )
          {
            // assign each first element of row and column to a gap
            AddAssignmentFromVerticalPath( REV_ROW_ITR, REV_COLUMN_ITR, ALIGNMENT);
            AddAssignmentFromHorizontalPath( REV_ROW_ITR, REV_COLUMN_ITR, ALIGNMENT);
          }
          else
          {
            BCL_Exit( "Cannot determine if column " + util::Format()( COLUMN) + " is an origin", -1);
          }

          ++REV_COLUMN_ITR;
          --COLUMN;
        }
      }
      else
      {
        BCL_Exit( "Cannot add gap to beginning of alignment: row=" + util::Format()( ROW) + ", column=" + util::Format()( COLUMN), -1);
      }
    } // AddGapsToBeginningOfAlignment

    //! @brief AddAssignmentFromDiagonalPath adds an assignment of assigned elements to the new alignment
    //! @param ITERATOR_ALIGNMENT_VERTICAL iterator to current element of ALIGNMENT_VERTICAL being assigned
    //! @param ITERATOR_ALIGNMENT_HORIZONTAL iterator to current element of ALIGNMENT_HORIZONTAL being assigned
    //! @param ALIGNMENT is the new alignment which is being built
    template< typename t_Member>
    void AlignerDynamicProgramming< t_Member>::AddAssignmentFromDiagonalPath
    (
      const typename AlignmentLeaf< t_Member>::const_reverse_iterator &ITERATOR_ALIGNMENT_VERTICAL,
      const typename AlignmentLeaf< t_Member>::const_reverse_iterator &ITERATOR_ALIGNMENT_HORIZONTAL,
      AlignmentNode< t_Member> &ALIGNMENT
    ) const
    {
      // create SharedPointer new_assignment to a new Assignment and append vertical_members
      util::ShPtr< Assignment< t_Member> > new_assignment( new Assignment< t_Member>());
      // add the members from ITERATOR_ALIGNMENT_VERTICAL
      new_assignment->Append( ( **ITERATOR_ALIGNMENT_VERTICAL).GetMembers());
      // add the members from ITERATOR_ALIGNMENT_HORIZONTAL
      new_assignment->Append( ( **ITERATOR_ALIGNMENT_HORIZONTAL).GetMembers());

      // insert new_assignment into ALIGNMENT
      ALIGNMENT.Prepend( new_assignment);
    } // AddAssignmentFromDiagonalPath

    //! @brief AddAssignmentFromHorizontalPath adds a gapped assignment
    //! @param ITERATOR_ALIGNMENT_VERTICAL iterator to current element of ALIGNMENT_VERTICAL; not being assigned
    //! @param ITERATOR_ALIGNMENT_HORIZONTAL itr to current element of ALIGNMENT_HORIZONTAL being assigned with gaps
    //! @param ALIGNMENT is the new alignment which is being built
    template< typename t_Member>
    void AlignerDynamicProgramming< t_Member>::AddAssignmentFromHorizontalPath
    (
      const typename AlignmentLeaf< t_Member>::const_reverse_iterator &ITERATOR_ALIGNMENT_VERTICAL,
      const typename AlignmentLeaf< t_Member>::const_reverse_iterator &ITERATOR_ALIGNMENT_HORIZONTAL,
      AlignmentNode< t_Member> &ALIGNMENT
    ) const
    {
      // create SharedPointer new_assignment to new Assignment and append horizontal_members
      util::ShPtr< Assignment< t_Member> > new_assignment( new Assignment< t_Member>());
      // add gaps equal to the number of members from ITERATOR_ALIGNMENT_VERTICAL
      size_t vertical_size( ( **ITERATOR_ALIGNMENT_VERTICAL).GetMembers().GetSize());
      new_assignment->Append( util::SiPtrList< const t_Member>( vertical_size));
      // add the members from ITERATOR_ALIGNMENT_HORIZONTAL
      new_assignment->Append( ( **ITERATOR_ALIGNMENT_HORIZONTAL).GetMembers());

      // insert new_assignment into ALIGNMENT
      ALIGNMENT.Prepend( new_assignment);
    } // AddAssignmentFromHorizontalPath

    //! @brief AddAssignmentFromVerticalPath adds a gapped assignment
    //! @param ITERATOR_ALIGNMENT_VERTICAL itr to current element of ALIGNMENT_VERTICAL being assigned with gaps
    //! @param ITERATOR_ALIGNMENT_HORIZONTAL iterator to current element of ALIGNMENT_HORIZONTAL; not being assigned
    //! @param ALIGNMENT is the new alignment which is being built
    template< typename t_Member>
    void AlignerDynamicProgramming< t_Member>::AddAssignmentFromVerticalPath
    (
      const typename AlignmentLeaf< t_Member>::const_reverse_iterator &ITERATOR_ALIGNMENT_VERTICAL,
      const typename AlignmentLeaf< t_Member>::const_reverse_iterator &ITERATOR_ALIGNMENT_HORIZONTAL,
      AlignmentNode< t_Member> &ALIGNMENT
    ) const
    {
      // create SharedPointer new_assignment to new Assignment constructed from vertical_members
      util::ShPtr< Assignment< t_Member> > new_assignment( new Assignment< t_Member>());
      // add the members from ITERATOR_ALIGNMENT_VERTICAL
      new_assignment->Append( ( **ITERATOR_ALIGNMENT_VERTICAL).GetMembers());
      // add gaps equal to the number of members from ITERATOR_ALIGNMENT_HORIZONTAL
      size_t horizontal_size( ( **ITERATOR_ALIGNMENT_HORIZONTAL).GetMembers().GetSize());
      new_assignment->Append( util::SiPtrList< const t_Member>( horizontal_size));

      // insert new_assignment into ALIGNMENT
      ALIGNMENT.Prepend( new_assignment);
    } // AddAssignmentFromVerticalPath

    //! @brief UpdateScore computes the Alignment score from individual Assignments
    //! @param ALIGNMENT is the Alignment whose score is calculated
    //! @return returns a double which is the score of ALIGNMENT
    template< typename t_Member>
    double AlignerDynamicProgramming< t_Member>::UpdateScore( AlignmentInterface< t_Member> &ALIGNMENT) const
    {
      double score( 0);

      //iterate over all
      typename AlignmentLeaf< t_Member>::const_iterator assign_itr( ALIGNMENT.GetAssignments().Begin());
      typename AlignmentLeaf< t_Member>::const_iterator assign_itr_end( ALIGNMENT.GetAssignments().End());
      for( ; assign_itr != assign_itr_end; ++assign_itr)
      {
        score += m_ScoreAssignment( **assign_itr);
      }

      return score;
    } // UpdateScore

  } // namespace align
} // namespace bcl

#endif // BCL_ALIGN_ALIGNER_DYNAMIC_PROGRAMMING_H_
