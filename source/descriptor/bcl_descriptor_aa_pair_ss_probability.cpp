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
#include "descriptor/bcl_descriptor_aa_pair_ss_probability.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_running_average_sd.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AAPairSSProbability::s_Instance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AAPairSSProbability()
      )
    );

    // Map from initializer string/type to the shptrvector of model interfaces; saves loading potentially
    // gigantic models repetitively
    storage::Map< util::ObjectDataLabel, storage::Vector< linal::Matrix< float> > > AAPairSSProbability::s_ProbabilityMapStorage =
      storage::Map< util::ObjectDataLabel, storage::Vector< linal::Matrix< float> > >();

    //  Mutex for access to s_ProbabilityMapStorage
    sched::Mutex AAPairSSProbability::s_ProbabilityMapMutex = sched::Mutex();

    //! @brief SSInfoType as string
    //! @param SS_INFO_TYPE the message level
    //! @return the SSInfoType as string
    const std::string &AAPairSSProbability::GetSSInfoTypeString( const SSInfoType &SS_INFO_TYPE)
    {
      static const std::string s_names [] =
      {
        "Helix",         // P(central residue is in a helix)
        "Strand",        // P(central residue is in a strand)
        "Coil",          // P(central residue is in a coil)
        "ToHelix",       // P(SS transition to helix over the stencil)
        "ToStrand",      // P(SS transition to strand over the stencil)
        "ToCoil",        // P(SS transition to coil over the stencil)
        "FromHelix",     // P(SS transition from helix over the stencil)
        "FromStrand",    // P(SS transition from strand over the stencil)
        "FromCoil",      // P(SS transition from coil over the stencil)
        GetStaticClassName< SSInfoType>()
      };
      return s_names[ SS_INFO_TYPE];
    }

    //! @brief SSInfoType directionality (relative to the central residue that the descriptor considers)
    //! @param SS_INFO_TYPE the type of interest
    //! @return -1 if the left is considered, 1 if the right is considered, 0 if both are considered
    int AAPairSSProbability::GetDirectionality( const SSInfoType &SS_INFO_TYPE)
    {
      static const int s_directions [] =
      {
        0,  // P(central residue is in a helix)
        0,  // P(central residue is in a strand)
        0,  // P(central residue is in a coil)
        -1, // P(SS transition to helix over the stencil)
        -1, // P(SS transition to strand over the stencil)
        -1, // P(SS transition to coil over the stencil)
        1,  // P(SS transition from helix over the stencil)
        1,  // P(SS transition from strand over the stencil)
        1,  // P(SS transition from coil over the stencil)
        0
      };
      return s_directions[ SS_INFO_TYPE];
    }

    //! @brief Get the column offset in the input file for the given type
    //! @param SS_INFO_TYPE the type of interest
    //! @return Column offset (4 for helical types, 7 for strand types, 10 for coil types)
    size_t AAPairSSProbability::GetColumnOffset( const SSInfoType &SS_INFO_TYPE)
    {
      static const int s_offset [] =
      {
        0,  // P(central residue is in a helix)
        1,  // P(central residue is in a strand)
        2, // P(central residue is in a coil)
        0,  // P(SS transition to helix over the stencil)
        1,  // P(SS transition to strand over the stencil)
        2, // P(SS transition to coil over the stencil)
        0,  // P(SS transition from helix over the stencil)
        1,  // P(SS transition from strand over the stencil)
        2, // P(SS transition from coil over the stencil)
        0
      };
      return s_offset[ SS_INFO_TYPE];
    }

    namespace
    {
      //! @brief create the aa type map
      storage::Vector< size_t> CreateAATypeMap()
      {
        storage::Vector< size_t> mapping( size_t( 256), util::GetUndefined< size_t>());
        for( size_t type_id( 0); type_id < biol::AATypes::s_NumberStandardAATypes; ++type_id)
        {
          mapping( int( biol::AAType( type_id)->GetOneLetterCode())) = type_id;
        }
        mapping( int( 'X')) = biol::AATypes::s_NumberStandardAATypes;
        return mapping;
      }
    }

    //! @brief get a vector that takes as input any letter, and returns the corresponding AA type index
    //! @return the vector
    const storage::Vector< size_t> &AAPairSSProbability::GetAATypeId()
    {
      static const storage::Vector< size_t> s_mapping( CreateAATypeMap());
      return s_mapping;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new BaseElement
    AAPairSSProbability *AAPairSSProbability::Clone() const
    {
      return new AAPairSSProbability( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AAPairSSProbability::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AAPairSSProbability::GetAlias() const
    {
      static const std::string s_name( "AAPairSSProbability");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AAPairSSProbability::GetNormalSizeOfFeatures() const
    {
      return 1;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief set up a types vector and find the significant types
    //! @param TYPES types vector of interest
    //! @param VALUES vector of all significant values; should be zeroed before calling this function
    //! @param SIGNIFICANT_TYPES an initially empty vector; will contain all significant indices of TYPES
    void FindSignificantTypes
    (
      const linal::VectorConstInterface< float> &TYPES,
      linal::Vector< float> &VALUES,
      storage::Vector< size_t> &SIGNIFICANT_TYPES
    )
    {
      for( size_t aa_type_val( 0); aa_type_val < biol::AATypes::s_NumberStandardAATypes; ++aa_type_val)
      {
        if( !math::EqualWithinAbsoluteTolerance( 0.0, TYPES( aa_type_val), 0.0001))
        {
          SIGNIFICANT_TYPES.PushBack( aa_type_val);
          VALUES( aa_type_val) = TYPES( aa_type_val);
        }
      }
      // if there were no significant types, add undefined type
      if( SIGNIFICANT_TYPES.IsEmpty())
      {
        SIGNIFICANT_TYPES.PushBack( biol::AATypes::s_NumberStandardAATypes);
        VALUES( biol::AATypes::s_NumberStandardAATypes) = 1.0;
      }
      else
      {
        VALUES /= linal::AbsSum( VALUES);
      }
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT: the element of interest
    //! @param STORAGE storage for the descriptor
    void AAPairSSProbability::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      const int directionality( GetDirectionality());
      storage::Vector< size_t> types( m_Probabilities->GetSize(), biol::AATypes::s_NumberStandardAATypes);

      // reset internal arrays
      for( size_t i( 0), n_arrays( m_SignificantTypes.GetSize()); i < n_arrays; ++i)
      {
        m_SignificantTypes( i).Reset();
        m_TypeCalculatorOutput( i) = float( 0.0);
      }
      FindSignificantTypes
      (
        m_TypeCalculator->operator ()( Iterator< biol::AABase>( ELEMENT)),
        m_TypeCalculatorOutput.LastElement(),
        m_SignificantTypes.LastElement()
      );

      // determine the types
      storage::Vector< storage::Vector< size_t> >::iterator itr_sig_types( m_SignificantTypes.Begin());
      storage::Vector< linal::Vector< float> >::iterator itr_type_values( m_TypeCalculatorOutput.Begin());
      if( directionality <= 0)
      {
        // copy the iterator
        Iterator< biol::AABase> itr_element( ELEMENT);
        for
        (
          size_t stencil_id( 0), stencil_size( m_Stencil.GetSize());
          stencil_id < stencil_size;
          ++stencil_id, ++itr_type_values, ++itr_sig_types
        )
        {
          const int target_pdb_id( ELEMENT->GetPdbID() - m_Stencil( stencil_id));
          while( itr_element.GetPosition() && itr_element( 0)->GetPdbID() > target_pdb_id)
          {
            --itr_element;
          }
          FindSignificantTypes( m_TypeCalculator->operator ()( itr_element), *itr_type_values, *itr_sig_types);
        }
      }
      if( directionality >= 0)
      {
        // copy the iterator
        Iterator< biol::AABase> itr_element( ELEMENT);
        for
        (
          size_t stencil_id( 0), stencil_size( m_Stencil.GetSize());
          stencil_id < stencil_size;
          ++stencil_id, ++itr_type_values, ++itr_sig_types
        )
        {
          const int target_pdb_id( ELEMENT->GetPdbID() + m_Stencil( stencil_id));
          while( itr_element.NotAtEnd() && itr_element( 0)->GetPdbID() < target_pdb_id)
          {
            ++itr_element;
          }
          if( !itr_element.NotAtEnd())
          {
            --itr_element;
          }
          FindSignificantTypes( m_TypeCalculator->operator ()( itr_element), *itr_type_values, *itr_sig_types);
        }
      }

      // now the types array has been loaded up, all that remains is to calculate the probabilities
      // this needs to be handled differently based on directionality
      // The directional cases are actually simplest since we just take a raw average, so handle those first
      const linal::Vector< float> &primary_type_values( m_TypeCalculatorOutput.LastElement());
      const storage::Vector< size_t> &primary_significant_type( m_SignificantTypes.LastElement());
      itr_sig_types = m_SignificantTypes.Begin();
      itr_type_values = m_TypeCalculatorOutput.Begin();
      storage::Vector< linal::Matrix< float> >::const_iterator
        itr_matrices( m_Probabilities->Begin()), itr_matrices_end( m_Probabilities->End());
      math::RunningAverageSD< float> average_value_pp, average_value_nn, average_value_pn, average_value_np;
      for
      (
        ;
        itr_matrices != itr_matrices_end;
        ++itr_sig_types, ++itr_type_values, ++itr_matrices
      )
      {
        for
        (
          storage::Vector< size_t>::const_iterator
            itr_psig_type( primary_significant_type.Begin()), itr_psig_type_end( primary_significant_type.End());
          itr_psig_type != itr_psig_type_end;
          ++itr_psig_type
        )
        {
          const size_t ptype_id( *itr_psig_type);
          const float pweight( primary_type_values( ptype_id));
          math::RunningAverage< float> ave_this_type_p, ave_this_type_n;
          for
          (
            storage::Vector< size_t>::const_iterator
              itr_ssig_type( itr_sig_types->Begin()), itr_ssig_type_end( itr_sig_types->End());
            itr_ssig_type != itr_ssig_type_end;
            ++itr_ssig_type
          )
          {
            const size_t stype_id( *itr_ssig_type);
            const float sweight( itr_type_values->operator()( stype_id));
            ( sweight > 0.0 ? ave_this_type_p : ave_this_type_n).AddWeightedObservation
            (
              itr_matrices->operator()( ptype_id, stype_id),
              math::Absolute( sweight)
            );
          }
          if( pweight > 0.0)
          {
            if( ave_this_type_p.GetWeight())
            {
              average_value_pp.AddWeightedObservation( ave_this_type_p.GetAverage(), pweight);
            }
            if( ave_this_type_n.GetWeight())
            {
              average_value_pn.AddWeightedObservation( ave_this_type_n.GetAverage(), pweight);
            }
          }
          else
          {
            if( ave_this_type_p.GetWeight())
            {
              average_value_np.AddWeightedObservation( ave_this_type_p.GetAverage(), -pweight);
            }
            if( ave_this_type_n.GetWeight())
            {
              average_value_nn.AddWeightedObservation( ave_this_type_n.GetAverage(), -pweight);
            }
          }
        }
      }
      float total_average( 0);
      if( average_value_pp.GetWeight())
      {
        total_average += average_value_pp.GetAverage();
      }
      if( average_value_nn.GetWeight())
      {
        total_average -= average_value_nn.GetAverage();
      }
      // handle the case that all blast values were at or below zero.  This should typically happen only with
      // unnatural AA types
      if( !average_value_pp.GetWeight() || !average_value_nn.GetWeight())
      {
        STORAGE( 0) = average_value_pp.GetAverage() + average_value_nn.GetAverage();
      }
      else
      {
        // estimate overall energy of the configuration
        STORAGE( 0)
          = 4 * average_value_pp.GetAverage() - average_value_pn.GetAverage()
            - average_value_np.GetAverage() - average_value_nn.GetAverage();
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAPairSSProbability::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes aa-type pairwise probabilities"
      );
      parameters.AddInitializer
      (
        "filename",
        "file that contains the statistics to read",
        io::Serialization::GetAgentInputFilename( &m_Filename)
      );
      parameters.AddInitializer
      (
        "probability type",
        "type of probability to consider",
        io::Serialization::GetAgent( &m_Type)
      );
      parameters.AddInitializer
      (
        "aa type",
        "AA type descriptor to consider (must return 20 values per AA)",
        io::Serialization::GetAgent( &m_TypeCalculator),
        "AAType"
      );
      parameters.AddInitializer
      (
        "stencil",
        "AA relative distances to consider (not counting the central residue)",
        io::Serialization::GetAgentContainerWithCheck
        (
          &m_Stencil,
          io::Serialization::GetAgentWithRange( size_t( 1), size_t( 6))
        ),
        "(1,5)"
      );
      parameters.AddOptionalInitializer
      (
        "force alignment",
        "Force the alignment from the type implied by probability type.  For example, ToHelix normally implies Left "
        "because the amino acid(s) to the left would normally control a To type transition; Helix defaults to Center, "
        "FromHelix defaults to Right",
        io::Serialization::GetAgent( &m_ForceAlignment)
      );
      return parameters;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool AAPairSSProbability::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      if( !m_TypeCalculator.IsDefined())
      {
        ERR_STREAM << "A valid AA type descriptor is required";
        return false;
      }
      else if( m_TypeCalculator->GetSizeOfFeatures() != size_t( 20))
      {
        ERR_STREAM << "A descriptor that returns 20 values is required for the AA type descriptor";
        return false;
      }

      // sort the stencil, in case the user left it out of order
      std::sort( m_Stencil.Begin(), m_Stencil.End());

      // get the expected # of matrices (twice the size of the stencil for bidirectional stencil types, same size for others)
      const int directionality( GetDirectionality());
      const size_t expected_number_matrices( m_Stencil.GetSize() * ( directionality ? 1 : 2));

      // one extra row/column in matrix for undefined AA type
      const size_t number_standard_types( biol::AATypes::s_NumberStandardAATypes);
      const size_t number_aa_types( number_standard_types + 1);

      // create temporary vectors
      m_TypeCalculatorOutput
        = storage::Vector< linal::Vector< float> >
          (
            expected_number_matrices + 1,
            linal::Vector< float>( number_aa_types)
          );
      m_SignificantTypes = storage::Vector< storage::Vector< size_t> >
          (
            expected_number_matrices + 1,
            storage::Vector< size_t>( number_aa_types)
          );

      util::ObjectDataLabel this_label( this->GetLabel());
      s_ProbabilityMapMutex.Lock();
      storage::Map< util::ObjectDataLabel, storage::Vector< linal::Matrix< float> > >::const_iterator
        itr_prob( s_ProbabilityMapStorage.Find( this_label));
      if( itr_prob != s_ProbabilityMapStorage.End())
      {
        // file was already loaded and normalized, just use it directly
        m_Probabilities = util::ToSiPtr( itr_prob->second);
        s_ProbabilityMapMutex.Unlock();
        return true;
      }

      // create the object in the map with the correct size
      storage::Vector< linal::Matrix< float> > &probabilities( s_ProbabilityMapStorage[ this_label]);

      probabilities
        = storage::Vector< linal::Matrix< float> >
          (
            expected_number_matrices,
            linal::Matrix< float>( number_aa_types, number_aa_types, float( 0.0))
          );
      // attach the SiPtr to the newly-created matrices
      m_Probabilities = util::ToSiPtr( probabilities);

      // Load the file
      io::IFStream statistics_file;
      io::File::MustOpenIFStream( statistics_file, m_Filename);
      storage::Vector< storage::Vector< std::string> > tokenized_lines
      (
        util::SplittedStringLineListFromIStream( statistics_file)
      );
      io::File::CloseClearFStream( statistics_file);

      // get the number of columns to determine the proper bins
      const size_t n_columns( tokenized_lines.FirstElement().GetSize());

      // determine whether there is a distance column
      const bool has_distance_column
      (
        n_columns > size_t( 3)
        && util::LengthOfIntegerType( tokenized_lines.FirstElement()( 2))
           == tokenized_lines.FirstElement()( 2).size()
      );

      // determine whether there is a distance column
      const bool has_counts_column
      (
        n_columns > size_t( 3 + has_distance_column)
        && util::LengthOfIntegerType( tokenized_lines.FirstElement()( 2 + has_distance_column))
           == tokenized_lines.FirstElement()( 2 + has_distance_column).size()
      );

      const bool is_from_type( GetDirectionality( m_Type) < 0);
      const bool is_to_type( GetDirectionality( m_Type) > 0);

      // determine whether there are transition columns
      const bool has_transition_columns( int( n_columns - has_distance_column - has_counts_column) >= 11);
      if( !has_transition_columns && GetDirectionality( m_Type))
      {
        ERR_STREAM << "Given file lacked transition state columns, but they are required for this probability type";
        return false;
      }

      // place statistics from each line into the proper matrices
      bool had_negative_values( false);

      // the general file format is
      // AAType1 AAType2 Distance Count HHProb XHProb HXProb SSProb XSProb SXProb CCProb XCProb CXProb (possibly other columns)
      const size_t type_type_column( has_distance_column + has_counts_column + 2 + ( 1 + 2 * has_transition_columns) * GetColumnOffset( m_Type));
      const size_t to_type_column( has_transition_columns ? type_type_column + 1 : util::GetUndefined< size_t>());
      const size_t from_type_column( has_transition_columns ? type_type_column + 2 : util::GetUndefined< size_t>());

      storage::Vector< size_t> stencil_mapping( size_t( 7), util::GetUndefined< size_t>());
      for( size_t distance( 0); distance < m_Stencil.GetSize(); ++distance)
      {
        stencil_mapping( m_Stencil( distance) - 1) = distance;
      }

      // get the AA type mapping for fast lookups of aa's from one letter code
      const storage::Vector< size_t> &one_letter_code_to_aa_type( GetAATypeId());

      // vector to hold all counts info
      storage::Vector< linal::Matrix< size_t> > counts
      (
        expected_number_matrices,
        linal::Matrix< size_t>( number_aa_types, number_aa_types, size_t( 0))
      );

      for
      (
        storage::Vector< storage::Vector< std::string> >::const_iterator
          itr_line( tokenized_lines.Begin()), itr_line_end( tokenized_lines.End());
        itr_line != itr_line_end;
        ++itr_line
      )
      {
        // check that the line is of sufficient size
        if( itr_line->GetSize() < size_t( type_type_column + 1))
        {
          continue;
        }

        const storage::Vector< std::string> &tokens( *itr_line);

        // get the distance
        const size_t distance
        (
          has_distance_column
          ? util::ConvertStringToNumericalValue< size_t>( tokens( 2)) - 1
          : 0
        );

        // test whether the distance will be used
        const size_t array_pos_left( stencil_mapping( distance));
        if( distance > size_t( 6) || array_pos_left >= expected_number_matrices)
        {
          continue;
        }
        const size_t array_pos_right( array_pos_left + m_Stencil.GetSize());

        // get the aa type ids
        const size_t first_aa_type_id( one_letter_code_to_aa_type( int( tokens( 0)[ 0])));
        const size_t second_aa_type_id( one_letter_code_to_aa_type( int( tokens( 1)[ 0])));

        // get the count
        const size_t count
        (
          has_counts_column
          ? util::ConvertStringToNumericalValue< size_t>( tokens( 3))
          : size_t( 1000)
        );

        // get the probability that the residue maintains the ss type of interest (strand/helix/coil)
        const float prob_maintain_type( util::ConvertStringToNumericalValue< float>( tokens( type_type_column)));

        // get the probability that the residue changes the ss type to the type of interest
        const float prob_to_type
        (
          has_transition_columns
          ? util::ConvertStringToNumericalValue< float>( tokens( to_type_column))
          : 0.0
        );

        // get the probability that the residue changes the ss type from the type of interest
        const float prob_from_type
        (
          has_transition_columns
          ? util::ConvertStringToNumericalValue< float>( tokens( from_type_column))
          : 0.0
        );
        had_negative_values = had_negative_values || prob_maintain_type < 0.0 || prob_to_type < 0.0 || prob_from_type < 0.0;

        if( util::IsDefined( second_aa_type_id) && util::IsDefined( first_aa_type_id))
        {
          // add the values to the matrices
          // for the left matrix, invert the values.  In this way, the row will always refer to the amino acid type of
          // interest, and the column refers to the target amino acid type
          if( !directionality)
          {
            probabilities( array_pos_left)( second_aa_type_id, first_aa_type_id) += prob_maintain_type + prob_to_type;
            probabilities( array_pos_right)( first_aa_type_id, second_aa_type_id) += prob_maintain_type + prob_from_type;
            counts( array_pos_left)( second_aa_type_id, first_aa_type_id) += count;
            counts( array_pos_right)( first_aa_type_id, second_aa_type_id) += count;
          }
          else if( directionality < 0)
          {
            probabilities( array_pos_left)( second_aa_type_id, first_aa_type_id) += is_to_type ? prob_to_type : is_from_type ? prob_from_type : prob_maintain_type;
            counts( array_pos_left)( second_aa_type_id, first_aa_type_id) += count;
          }
          else if( directionality > 0)
          {
            probabilities( array_pos_left)( first_aa_type_id, second_aa_type_id) += is_to_type ? prob_to_type : is_from_type ? prob_from_type : prob_maintain_type;
            counts( array_pos_left)( first_aa_type_id, second_aa_type_id) += count;
          }
        }
        else if( !util::IsDefined( second_aa_type_id) && !util::IsDefined( first_aa_type_id))
        {
          // both types are wildcards, should not happen
          BCL_Exit( "Undefined aa types: " + tokens( 0) + " " + tokens( 1), -1);
        }
        else if( util::IsDefined( first_aa_type_id))
        {
          for( size_t wildcard( 0); wildcard < number_standard_types; ++wildcard)
          {
            // add the values to the matrices
            // for the left matrix, invert the values.  In this way, the row will always refer to the amino acid type of
            // interest, and the column refers to the target amino acid type
            if( !directionality)
            {
              probabilities( array_pos_left)( wildcard, first_aa_type_id) += prob_maintain_type + prob_to_type;
              probabilities( array_pos_right)( first_aa_type_id, wildcard) += prob_maintain_type + prob_from_type;
            }
            else if( directionality < 0)
            {
              probabilities( array_pos_left)( wildcard, first_aa_type_id) += is_to_type ? prob_to_type : is_from_type ? prob_from_type : prob_maintain_type;
            }
            else if( directionality > 0)
            {
              probabilities( array_pos_left)( first_aa_type_id, wildcard) += is_to_type ? prob_to_type : is_from_type ? prob_from_type : prob_maintain_type;
            }
          }
        }
        else
        {
          for( size_t wildcard( 0); wildcard < number_standard_types; ++wildcard)
          {
            // add the values to the matrices
            // for the left matrix, invert the values.  In this way, the row will always refer to the amino acid type of
            // interest, and the column refers to the target amino acid type
            if( !directionality)
            {
              probabilities( array_pos_left)( second_aa_type_id, wildcard) += prob_maintain_type + prob_to_type;
              probabilities( array_pos_right)( wildcard, second_aa_type_id) += prob_maintain_type + prob_from_type;
            }
            else if( directionality < 0)
            {
              probabilities( array_pos_left)( second_aa_type_id, wildcard) += is_to_type ? prob_to_type : is_from_type ? prob_from_type : prob_maintain_type;
            }
            else if( directionality > 0)
            {
              probabilities( array_pos_left)( wildcard, second_aa_type_id) += is_to_type ? prob_to_type : is_from_type ? prob_from_type : prob_maintain_type;
            }
          }
        }

      }

      if( !had_negative_values)
      {
        // add pseudo-counts
        const size_t pseudocount( 20);
        for( size_t aatype_id( 0); aatype_id < number_standard_types; ++aatype_id)
        {
          // depending on where the statistics were collected, they may be heavily biased.
          // Therefore, as a heuristic to reduce noise, add a weighted pseudo-count of 10, using the AAType-UnknownType probabilities
          for( size_t stencil_pos( 0); stencil_pos < expected_number_matrices; ++stencil_pos)
          {
            linal::Matrix< float> &stencil_probs( probabilities( stencil_pos));
            linal::Matrix< size_t> &stencil_counts( counts( stencil_pos));
            for( size_t aatype_id_b( 0); aatype_id_b < number_standard_types; ++aatype_id_b)
            {
              const size_t count( stencil_counts( aatype_id, aatype_id_b));
              float &stencil_prob_ab( stencil_probs( aatype_id, aatype_id_b));
              const float stencil_prob_ax( stencil_probs( number_standard_types, aatype_id_b));
              const float stencil_prob_xb( stencil_probs( aatype_id, number_standard_types));
              const float ave_stencil_prob_axb( 0.5 * ( stencil_prob_ax + stencil_prob_xb));
              stencil_prob_ab
                = log( ( stencil_prob_ab * count + pseudocount * ave_stencil_prob_axb + 0.5) / float( count + pseudocount));
            }
          }
        }
        // TODO Consider normalization for cases e_Helix, e_Strand, and e_Coil, specifically like in the Dunbrack
        // paper, ala S in P(X,Y) * P(Y,Z) / (S * P(Y))
        // For now, just take roots
      }
      // BCL_Debug( *m_Probabilities);

      s_ProbabilityMapMutex.Unlock();
      return true;
    }

    //! @brief get the effective directionality (considering implicit direction of m_Type, and m_ForceAlignment)
    //! @return -1 for left, 0 for center, 1 for right
    int AAPairSSProbability::GetDirectionality() const
    {
      return m_ForceAlignment == s_NumberOfWindowAlignments
             ? GetDirectionality( m_Type)
             : m_ForceAlignment == e_Left
               ? -1
               : m_ForceAlignment == e_Right ? 1 : 0;
    }

  } // namespace descriptor
} // namespace bcl
