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
#include "contact/bcl_contact_calculate_correlations_sm.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_alignment_interface.h"
#include "biol/bcl_biol_aa_base.h"
#include "math/bcl_math_running_average_sd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

  //////////
  // data //
  //////////

    //! Default is e_BLOSUM_62
    const score::AAAssignmentBLOSUM CalculateCorrelationsSM::s_SimilarityMatrix( score::AAAssignmentBLOSUM::e_BLOSUM_62);

    // TODO: Use assign score so you have a general score with type set as parameter

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CalculateCorrelationsSM::s_Instance
    (
      GetObjectInstances().AddInstance( new CalculateCorrelationsSM())
    );

    const double CalculateCorrelationsSM::s_GapCutoff( 0.1);

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CalculateCorrelationsSM::CalculateCorrelationsSM() :
      m_WeightMap() // holds all the weights for all pairs of sequences k and l
    {
    }

    //! @brief Clone function
    //! @return pointer to new CalculateCorrelationsSM
    CalculateCorrelationsSM *CalculateCorrelationsSM::Clone() const
    {
      return new CalculateCorrelationsSM( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CalculateCorrelationsSM::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief calculates a correlation matrix if given an AlignmentInterface of type bio::AABase
    //! @details
    //! @param ALIGNMENT_INTERFACE is the MSA representation from which the CM is calculated
    //! @return Returns a correlation matrix of dimensions N by N where N is the length of the MSA, containing doubles

    // Currently returns nan for perfectly conserved positions
    CorrelationMatrix CalculateCorrelationsSM::operator()( const align::AlignmentInterface< biol::AABase> &ALIGNMENT_INTERFACE) const
    {
      // Acquire N length of MSA, number of sequences K from AlignmentInterface
      const size_t msa_length( ALIGNMENT_INTERFACE.GetSize());
      const size_t num_sequences( ALIGNMENT_INTERFACE.GetDepth());

      // Create vector to hold matrices with similarity values for each column of the MSA
      storage::Vector< linal::Matrix< double> > similarity_matrix_vector;
      similarity_matrix_vector.AllocateMemory( msa_length);

      // Create vector to hold all statistics as saved via RunningAverageSD< double> indexed by position pos_index
      storage::Vector< math::RunningAverageSD< double> > stat_mean_sd_vector;

      // Get the list of assignments
      const util::ShPtrList< align::Assignment< biol::AABase> > &assignments( ALIGNMENT_INTERFACE.GetAssignments());

      // Set sequence position index that keeps track of position in alignment
      size_t pos_index( 0);

      // Itr Alignment Interface across all nodes and get assignments
      for
      (
        util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator
          assignment_itr( assignments.Begin()),
          assignment_itr_end( assignments.End());
        assignment_itr != assignment_itr_end;
        ++assignment_itr, ++pos_index
      )
      {
        // TODO: MOVE TO CONSTRUCTOR and SKIP ALLOCATION OF MEMORY to line 98
        // Create matrix of doubles to store similarity values to be calculated at this "column" of assignments and pushback into vector of matrices to keep track
        // TODO: Initialize all members of the matrix to zero, do I just add to the third field to accomplish that?
        similarity_matrix_vector.PushBack( linal::Matrix< double>( num_sequences, num_sequences));

        // Get SiPtrLst of members
        const util::SiPtrList< const biol::AABase> &members( ( *assignment_itr)->GetMembers());

        // Create DataSetStatisc
        math::RunningAverageSD< double> stat_mean_col_current;

        // Set sequence position k_index, and l_index
        size_t k_index( 0);
        size_t l_index( 0);

        for
        (
          util::SiPtrList< const biol::AABase>::const_iterator
            k_assign_itr( members.Begin()),
            assign_itr_end( members.End());
          k_assign_itr != assign_itr_end;
          ++k_assign_itr, ++k_index
        )
        {
          // Set l index and l iterator such that sequences aren't compared to self and all pairs are accounted for once
          l_index = k_index + 1;
          // Set l_assign_itr such that it is != to k_assign_itr
          util::SiPtrList< const biol::AABase>::const_iterator l_assign_itr( k_assign_itr);
          ++l_assign_itr;
          for( ; l_assign_itr != assign_itr_end; ++l_assign_itr, ++l_index)
          {
            // Calculate similarity value for Sikl using s_SimilarityMatrix
            // this code is ridiculous
            // Check for null pointers since GAPs are stored in alignments as null pointers and not AABases of type GAP
            const util::SiPtr< const biol::AABase> &aa_k_siptr( *k_assign_itr);
            const biol::AAType type_k( aa_k_siptr.IsDefined() ? aa_k_siptr->GetType() : biol::GetAATypes().GAP);
            // Check for null pointers since GAPs are stored in alignments as null pointers and not AABases of type GAP
            const util::SiPtr< const biol::AABase> &aa_l_siptr( *l_assign_itr);
            const biol::AAType type_l( aa_l_siptr.IsDefined() ? aa_l_siptr->GetType() : biol::GetAATypes().GAP);
            const double similarity_kl
            (
              s_SimilarityMatrix.Probability
              (
                s_SimilarityMatrix.e_BLOSUM_62,
                type_k,
                type_l
              )
            );

            // for current stat_mean_col_current add the value so that later it can computer mean and standard deviation
            stat_mean_col_current += similarity_kl;

            // add similarity value for Sikl to position k,l in the current similarity matrix
            similarity_matrix_vector( pos_index)( k_index, l_index) = similarity_kl;

            // add to weight if non-identical elements present in k and l in position k,l in the weights map keyed by
            // k,l
            const storage::Pair< size_t, size_t> temp_pair( k_index, l_index);
            // check to see if this pair is already in the m_WeightMap if no, insert pair along with current count for non-identical elements
            double non_identical_count( ( type_k != type_l));
            // TODO: Change map to symmetric matrix
            // TODO: check logic
            if( m_WeightMap.Has( temp_pair) && non_identical_count) // Only executed if non_identical_count is != 0
            {
              m_WeightMap[ temp_pair] += non_identical_count;
            }
            else
            {
              m_WeightMap.InsertElement( storage::Pair< storage::Pair< size_t, size_t>, double>( temp_pair, non_identical_count));
            }
          }
        }
        // pushback most recently finished statistics object into statistics vector
        stat_mean_sd_vector.PushBack( stat_mean_col_current);
      } // end alignment iterator for loop

      // TODO: Sum to 1 function in math matrix, change this and maybe make a linal::symmetric matrix

      // Normalize weights in m_WeightsMap such that they sum to 1
      // First iterate through Map and sum up all values in weights_sum
      double weights_sum( 0.0);
      for
      (
        storage::Map< storage::Pair< size_t, size_t>, double>::const_iterator
          map_itr( m_WeightMap.Begin()),
          map_itr_end( m_WeightMap.End());
        map_itr != map_itr_end;
        ++map_itr
      )
      {
        weights_sum += map_itr->second;
      }

      // Iterate through Map again, now saving the values as value/total sum of values to normalize
      for
      (
        storage::Map< storage::Pair< size_t, size_t>, double>::iterator
          map_itr( m_WeightMap.Begin()),
          map_itr_end( m_WeightMap.End());
        map_itr != map_itr_end;
        ++map_itr
      )
      {
        map_itr->second /= weights_sum;
//        m_WeightMap[ map_itr->first] = ( m_WeightMap[ map_itr->first] / weights_sum);
      }

      // Create CorrelationMatrix to fill up with values
      CorrelationMatrix correlation_matrix( msa_length);

      // Iterate through all pairs of k l summing values of formula for all pairings i j
      for( size_t pos_i( 0); pos_i < msa_length; ++pos_i)
      {
        for( size_t pos_j( pos_i + 1); pos_j < msa_length; ++pos_j)
        {
          // Keep track of sum of all pairs of kl
          double sum_for_ij( 0.0);

          // Iterate through all possible pairs of k and l applying the formula below
          // For all pairs of i and j where i != j
          // Correlation at Rij =
          // 1/N^2 *
          // Sum(for all pairs of sequences k and l where k != l)
          // (Weight(k,l)  *  ( similarity at i between k and l - average similarity at i) *
          // ( similarity at j between k and l - average similarity at j)) divided by
          // (standard deviation of similarities at column i)*(standard deviation of similarities at j)
          for( size_t seq_k( 0); seq_k < num_sequences; ++seq_k)
          {
            for( size_t seq_l( seq_k + 1); seq_l < num_sequences; ++seq_l)
            {
              // TODO: statistics function for this score, if not there implement and call here
              sum_for_ij += ( GetWeight( seq_k, seq_l)
                  * ( similarity_matrix_vector( pos_i)( seq_k, seq_l) - stat_mean_sd_vector( pos_i).GetAverage())
                  * ( similarity_matrix_vector( pos_j)( seq_k, seq_l) - stat_mean_sd_vector( pos_j).GetAverage()))
                  / ( stat_mean_sd_vector( pos_i).GetStandardDeviation() * stat_mean_sd_vector( pos_j).GetStandardDeviation());
            }
          }
          // Divide the sum of all pairs kl by the total sequence length squared and save that to r_ij
          correlation_matrix( pos_i, pos_j) = ( sum_for_ij)/( math::Sqr( msa_length));
          // Correlation matrix is diagonally symmetrical
          correlation_matrix( pos_j, pos_i) = correlation_matrix( pos_i, pos_j);
        }
      }
      // Return filled in correlation_matrix
      return correlation_matrix;
    }

    //! @brief Checks to see if param CalculateCorrelationsSM object is equal to this one
    //! @param CALCULATE_CORRELATIONS reference to a CalculateCorrelationsSM object
    //! @return bool true for all members are equal and false otherwise
    bool CalculateCorrelationsSM::operator==( const CalculateCorrelationsSM &CALCULATE_CORRELATIONS) const
    {
      if
      (
        s_GapCutoff == CALCULATE_CORRELATIONS.s_GapCutoff &&
          s_SimilarityMatrix.GetTableType() == CALCULATE_CORRELATIONS.s_SimilarityMatrix.GetTableType() &&
          m_WeightMap.GetSize() == CALCULATE_CORRELATIONS.m_WeightMap.GetSize()
      )
      {
        for
        (
          storage::Map< storage::Pair< size_t, size_t>, double>::const_iterator
            map_itr_a( m_WeightMap.Begin()),
            map_itr_b( CALCULATE_CORRELATIONS.m_WeightMap.Begin()),
            map_itr_end( m_WeightMap.End());
          map_itr_a != map_itr_end;
          ++map_itr_a, ++map_itr_b
        )
        {
          if
          (
            map_itr_a->first.First() != map_itr_b->first.First() ||
              map_itr_a->first.Second() != map_itr_b->first.Second() ||
              map_itr_a->second != map_itr_b->second
          )
          {
            return false;
          }
        }
        return true;
      }
      // else
      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CalculateCorrelationsSM::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_WeightMap, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CalculateCorrelationsSM::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_WeightMap, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Takes two sequences from the MSA given to CalculateWeights and returns the weight value for that pair
    //! @param K_SEQ_INDEX from the MSA for sequence k
    //! @param L_SEQ_INDEX from the MSA for sequence l
    //! @return Returns a double weight value
    double CalculateCorrelationsSM::GetWeight( const size_t K_SEQ_INDEX, const size_t L_SEQ_INDEX) const
    {
      // check to make sure that K_SEQ_INDEX is < L_SEQ_INDEX and if not swap them and use them as the keys for map lookup
      if( K_SEQ_INDEX < L_SEQ_INDEX)
      {
        return m_WeightMap.GetValue( storage::Pair< size_t, size_t>( K_SEQ_INDEX, L_SEQ_INDEX));
      }

      return m_WeightMap.GetValue( storage::Pair< size_t, size_t>( L_SEQ_INDEX, K_SEQ_INDEX));
    }

  } // namespace contact
} // namespace bcl
