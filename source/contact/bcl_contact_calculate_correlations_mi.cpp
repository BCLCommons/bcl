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
#include "contact/bcl_contact_calculate_correlations_mi.h"

// includes from bcl - sorted alphabetically
#include "align/bcl_align_alignment_interface.h"
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

    //! @brief ConstAATypePair used to make code more readable
    typedef storage::Pair< biol::AAType, biol::AAType> AATypePair;

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CalculateCorrelationsMI::s_Instance
    (
      GetObjectInstances().AddInstance( new CalculateCorrelationsMI())
    );

    const size_t CalculateCorrelationsMI::s_LogBase = 20;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CalculateCorrelationsMI::CalculateCorrelationsMI()
    {
    }

    //! @brief Clone function
    //! @return pointer to new CalculateCorrelationsMI
    CalculateCorrelationsMI *CalculateCorrelationsMI::Clone() const
    {
      return new CalculateCorrelationsMI( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CalculateCorrelationsMI::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    CorrelationMatrix CalculateCorrelationsMI::operator ()( const align::AlignmentInterface< biol::AABase> &ALIGNMENT_INTERFACE) const
    {
      // Acquire N length of MSA and number of sequences from ALIGNMENT_INTERFACE
      const size_t msa_length = ALIGNMENT_INTERFACE.GetSize();
      const size_t num_sequences = ALIGNMENT_INTERFACE.GetDepth();

      // Create a vector to store all the probabilities for base type x at position i
      storage::Vector< storage::Map< biol::AAType, double> > probability_vector( msa_length);

      // Create a matrix to hold each map of pair probabilities
      storage::SymmetricMatrix
      <
        storage::Map< AATypePair, double>
      >
        pair_probability_map_matrix( msa_length);

      // Get the list of assignments
      const util::ShPtrList< align::Assignment< biol::AABase> > assignments( ALIGNMENT_INTERFACE.GetAssignments());

      // Create vector to hold iterators pointing to each assignment list in the Alignment Interface
      // Using std::vector because storage::Vector tries to write out iterators and can't/not designed to
      std::vector< util::SiPtrList< const biol::AABase>::const_iterator> assignments_vector;
      assignments_vector.reserve( msa_length);

      // Iterate through assignment iterator vector setting each to beginning of corresponding members list
      util::ShPtrList< align::Assignment< biol::AABase> >::const_iterator
        assignments_itr( assignments.Begin()),
        assignments_itr_end( assignments.End());
      // Save end iterator for first position to keep track of how far to go
      util::SiPtrList< const biol::AABase>::const_iterator
        first_members_itr_end( ( *assignments_itr)->GetMembers().End());
      for
      (
        ;
        assignments_itr != assignments_itr_end;
        ++assignments_itr
      )
      {
        // Add the iterator for the current member list to the assignments_vector
        assignments_vector.push_back( ( *assignments_itr)->GetMembers().Begin());
      }

      // Iterate through each row represented by a vector of iterators pointing each to members in the same sequence
      // TODO: Use IS EMPTY?
      while( assignments_vector.front() != first_members_itr_end)
      {
        // Set sequence position index i which keeps track of positions in the alignment
        size_t i_index( 0);

        // Iterate through each position i
        for
        (
          std::vector< util::SiPtrList< const biol::AABase>::const_iterator>::iterator
            assignments_vector_itr( assignments_vector.begin()),
            assignments_vector_itr_end( assignments_vector.end());
          assignments_vector_itr != assignments_vector_itr_end;
          ++assignments_vector_itr, ++i_index
        )
        {
          // Save frequency data for type at this i in this row
          // Get current SiPtr and check to see if is NULL and set type appropriately
          util::SiPtr< const biol::AABase> aa_siptr( **assignments_vector_itr);
          biol::AAType temp_type_i( aa_siptr.IsDefined() ? aa_siptr->GetType() : biol::GetAATypes().GAP);

          // Check for key in map and insert if not present otherwise, increment count for type
          if( probability_vector( i_index).Has( temp_type_i))
          {
            ++( probability_vector( i_index)[ temp_type_i]); // #################### Not sure if this works!!!!!!
          }
          else
          {
            probability_vector( i_index).InsertElement
            (
              storage::Pair< biol::AAType, double>( temp_type_i, double( 1.0))
            );
          }

          // TODO: Change temp type stuff to a symmetric matrix using enums

          // Set sequence position index j which keeps track of positions in the alignment
          size_t j_index( 0);
          // Iterate through each position j where j is always > i
          std::vector< util::SiPtrList< const biol::AABase>::const_iterator>::iterator
            j_assignments_vector_itr( assignments_vector_itr);
          ++j_assignments_vector_itr;
          for( ; j_assignments_vector_itr != assignments_vector_itr_end; ++j_assignments_vector_itr, ++j_index)
          {
            // Get current SiPtr and check to see if is NULL and set type appropriately
            util::SiPtr< const biol::AABase> aa_siptr( **j_assignments_vector_itr);
            const biol::AAType temp_type_j( aa_siptr.IsDefined() ? aa_siptr->GetType() : biol::GetAATypes().GAP);

            // Create type ordered type pair to search map from types i and j
            AATypePair temp_pair;
            if( temp_type_i <= temp_type_j)
            {
              temp_pair = AATypePair( temp_type_i, temp_type_j);
            }
            else
            {
              temp_pair = AATypePair( temp_type_j, temp_type_i);
            }

            // Save a new map into pair_probability_map_matrix to keep track of x,y pair frequencies at position i,j
            if( pair_probability_map_matrix( i_index, j_index).IsEmpty())
            {
              pair_probability_map_matrix( i_index, j_index) = storage::Map< AATypePair, double >();
              pair_probability_map_matrix( i_index, j_index).InsertElement
              (
                storage::Pair< AATypePair, double>( temp_pair, double( 1.0))
              );
            }
            else
            {

              // Check to see if map at this i,j already has temp_pairing, increment x,y (x must always be <= y) or
              // add with count of 1 if not present yet
              if( pair_probability_map_matrix( i_index, j_index).Has( temp_pair))
              {
                ++( pair_probability_map_matrix( i_index, j_index)[ temp_pair]);
              }
              else
              {
                pair_probability_map_matrix( i_index, j_index).InsertElement
                  (
                    storage::Pair< AATypePair, double>( temp_pair, double( 1.0))
                  );
              }
            }
          } // End of iterator for position j
          // Increment position i iterator within the vector once done with all comparisons
          ( *assignments_vector_itr) = ++( *assignments_vector_itr);

        } // End of iterator for position i
      } // End of while loop which iterates through all rows in the ALIGNMENT_INTERFACE

      ////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////
      // Iterate through all counts and change all at current position into probabilities by dividing them by
      // num_sequences for all probabilities in the probability vector
      for( size_t pos_index( 0); pos_index < msa_length; ++pos_index)
      {
        for
        (
          storage::Map< biol::AAType, double>::iterator
            count_itr( probability_vector( pos_index).Begin()),
            count_itr_end( probability_vector( pos_index).End());
          count_itr != count_itr_end;
          ++count_itr
        )
        {
          probability_vector( pos_index)[ ( *count_itr).first] /= num_sequences;
        }
      }

      // Iterate through all positions in the matrix of pair probabilities and change all the probabilities to
      // decimal percentages by dividing by total num_sequences
      for( size_t i_index( 0); i_index < msa_length; ++i_index)
      {
        size_t j_index( i_index);
        ++j_index;
        for( ; j_index < msa_length; ++j_index)
        {
          for
          (
            storage::Map< AATypePair, double>::iterator
              map_itr( pair_probability_map_matrix( i_index, j_index).Begin()),
              map_itr_end( pair_probability_map_matrix( i_index, j_index).End());
            map_itr != map_itr_end;
            ++map_itr
          )
          {
            pair_probability_map_matrix( i_index, j_index)[ map_itr->first] /= num_sequences;
          }
        }
      }

      // Create correlation matrix to be filled in and returned
      CorrelationMatrix correlation_matrix;
      // Iterate through SymmetricMatrix with Maps of all probability pairs xi,yj- skips duplicate half across diagonal
      for( size_t i_index( 0); i_index < pair_probability_map_matrix.GetSize(); ++i_index)
      {
        // Set double sum_ij that keeps track of the sum at position i,j
        double sum_ij( 0.0);
        for( size_t j_index( 0); j_index <= i_index; ++j_index)
        {
          // Iterate through all frequency pairs at this i,j positionPredictorCorrelations
          for
          (
            storage::Map< AATypePair, double>::const_iterator
              pair_freq_itr( pair_probability_map_matrix( i_index, j_index).Begin()),
              pair_freq_itr_end( pair_probability_map_matrix( i_index, j_index).End());
            pair_freq_itr != pair_freq_itr_end;
            ++pair_freq_itr
          )
          {
            // Get types x and y at positions i and j respectively
            const biol::AAType &x_type( pair_freq_itr->first.First());
            const biol::AAType &y_type( pair_freq_itr->first.Second());

            //###### Should be possible to skip this for the sake of a slight speedup assuming that the above code
            //###### ensuring pairs are stored such that i <= j is not tampered with...
            // Create temp_pair for lookup such that in pair x,y x <= y
//            AATypePair temp_type_pair;
//            if( x_type <= y_type)
//            {
//              temp_type_pair = AATypePair( x_type, y_type);
//            }
//            else
//            {
//              temp_type_pair = AATypePair( y_type, x_type);
//            }
            // Get frequency for (xi,yj)
            const double prob_xi_yj( pair_freq_itr->second);

            // Get frequency for xi and for yj independently
            const double prob_xi( probability_vector( i_index).GetValue( x_type));
            const double prob_yj( probability_vector( j_index).GetValue( y_type));

            // Calculate for log_b(a) by calculating log(a)/log(b)
            const double log_result( log( ( prob_xi_yj) / ( prob_xi * prob_yj)) / log( double( s_LogBase)));
            sum_ij += prob_xi_yj * log_result;
          }
          // Save -sum_ij to correlation matrix
          correlation_matrix( i_index, j_index) = 0 - sum_ij;
        }
      }

      //!**********************************************************

//      // Iterate through probabilityies for all pairings xi,yj
//      for
//      (
//          storage::Vector< storage::Map< biol::AAType, double> >::const_iterator
//            i_prob_itr( probability_vector.Begin()),
//            prob_itr_end( probability_vector.End());
//          i_prob_itr != prob_itr_end;
//          ++i_prob_itr, ++i_index
//      )
//      {
//        // Set variable keeping track of the sum at position i,j
//        double sum_ij( 0.0);
//        // Set j_index as one position advanced from i_index
//        size_t j_index( i_index);
//        ++j_index;
//        // Set j position iterator one advanced from the i position iterator
//        storage::Vector< storage::Map< biol::AAType, double> >::const_iterator j_prob_iter = i_prob_itr;
//        ++j_prob_iter;
//        for
//        ( ; j_prob_iter != prob_itr_end; ++j_prob_iter, ++j_index)
//        {
//          // Set the first iterator for all the types in the probability map
//          for
//          (
//            storage::Map< biol::AAType, double>::const_iterator
//              k_type_itr( probability_vector(i_index).Begin()),
//              type_itr_end( probability_vector( i_index).End());
//            k_type_itr != type_itr_end;
//            ++k_type_itr
//          )
//          {
//            // Set the second iterator for all the types present in the probability map and start it one position advanced
//            storage::Map< biol::AAType, double>::const_iterator l_type_itr( k_type_itr);
//            ++l_type_itr;
//            for
//            ( ; l_type_itr != type_itr_end; ++l_type_itr)
//            {
//              double prob_xi( probability_vector( i_index).GetValue( k_type_itr->first));
//              double prob_yj( probability_vector( j_index).GetValue( l_type_itr->first));
//              // Create temp_pair for lookup such that in pair x,y x <= y
//              AATypePair temp_type_pair;
//              if ( k_type_itr->first <= l_type_itr->first)
//              {
//                temp_type_pair = AATypePair( k_type_itr->first, l_type_itr->first);
//              }
//              else
//              {
//                temp_type_pair = AATypePair( l_type_itr->first, k_type_itr->first);
//              }
//              double prob_xi_yj( pair_probability_map_matrix( i_index, j_index)[ temp_type_pair]);
//              // Calculate for log_b(a) by calculating log_2(a)/log_2(b)
//              double log_result( log2(( prob_xi_yj) / ( prob_xi * prob_yj)) / log2( double( 20.0)));
//              sum_ij +=  prob_xi_yj * log_result;
//            }
//          }
//        }
//
//        // Save the negative of the sum as the correlation for pair i,j in the correlation matrix
//        correlation_matrix( i_index, j_index) = sum_ij;
//      }

      return correlation_matrix;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CalculateCorrelationsMI::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CalculateCorrelationsMI::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace contact
} // namespace bcl
