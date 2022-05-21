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
#include "descriptor/bcl_descriptor_aa_triplet_lookup.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_reflecting.h"
#include "math/bcl_math_running_average.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AATripletLookup::s_Instance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AATripletLookup()
      )
    );

    // Map from initializer string/type to the shptrvector of model interfaces; saves loading potentially
    // gigantic models repetitively
    storage::Map< std::string, AATripletLookup::TripletSSPropensityType> AATripletLookup::s_ProbabilitiesStorage =
      storage::Map< std::string, AATripletLookup::TripletSSPropensityType>();

    //  Mutex for access to s_ProbabilityMapStorage
    sched::Mutex AATripletLookup::s_TripletProbabilitiesMutex = sched::Mutex();

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
        return mapping;
      }

      //! @brief create string that can be used to go from aa type index to one letter code
      std::string CreateAAOneLetterString()
      {
        std::string one_letter_codes;
        for
        (
          biol::AATypes::const_iterator itr( biol::GetAATypes().Begin()), itr_end( biol::GetAATypes().End());
          itr != itr_end;
          ++itr
        )
        {
          one_letter_codes += ( *itr)->GetOneLetterCode();
        }
        return one_letter_codes;
      }

      //! @brief Convert a hash to an actual sequence string
      //! @param HASH the hash to convert to a string
      //! @param SEQ_STORAGE storage space for the sequence string
      void ConvertHashToString( size_t HASH, std::string &SEQ_STORAGE)
      {
        static const std::string s_aa_one_letter_codes( CreateAAOneLetterString());
        for( int i( 0); i < 3; ++i)
        {
          SEQ_STORAGE[ i] = s_aa_one_letter_codes[ size_t( HASH % biol::AATypes::s_NumberStandardAATypes)];
          HASH /= size_t( biol::AATypes::s_NumberStandardAATypes);
        }
      }
    }
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new BaseElement
    AATripletLookup *AATripletLookup::Clone() const
    {
      return new AATripletLookup( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AATripletLookup::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AATripletLookup::GetAlias() const
    {
      static const std::string s_name( "AATripletLookup");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AATripletLookup::GetNormalSizeOfFeatures() const
    {
      return 3;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT: the element of interest
    //! @param STORAGE storage for the descriptor
    void AATripletLookup::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      storage::Vector< size_t> aa_index_window( size_t( 3));
      storage::Vector< size_t> aa_index_window_pssm( size_t( 3), util::GetUndefined< size_t>());
      storage::Vector< float> aa_index_window_pssm_weight( size_t( 3), 0.0);
      // get a window of size 3 with the nearby AA types
      {
        storage::Vector< biol::AAType> aa_type_window( size_t( 3));
        iterate::Reflecting< const biol::AABase> reflecting_itr( ELEMENT);
        --reflecting_itr;
        for( size_t i( 0), max_i( 3); i < max_i; ++i, ++reflecting_itr)
        {
          aa_type_window( i) = reflecting_itr->GetType();
          if( reflecting_itr->GetType().GetIndex() > size_t( 20))
          {
            continue;
          }
          if( !m_EffTypeCalculator.IsDefined())
          {
            // user does not want to use blast profile for claculation
            continue;
          }
          Iterator< biol::AABase> descriptor_iterator( reflecting_itr);
          linal::VectorConstReference< float> blast_probability( m_EffTypeCalculator->operator()( descriptor_iterator));
          if( blast_probability( aa_type_window( i)) < m_BlastProbabilityThreshold)
          {
            continue;
          }
          size_t aa_blast_index( 0);
          size_t max_blast_index( util::GetUndefined< size_t>());
          double max_blast_value( m_BlastProbabilityThreshold);
          const size_t actual_index( aa_type_window( i).GetIndex());
          for
          (
            linal::VectorConstReference< float>::const_iterator
              itr_prob( blast_probability.Begin()), itr_prob_end( blast_probability.End());
            itr_prob != itr_prob_end;
            ++itr_prob, ++aa_blast_index
          )
          {
            if( aa_blast_index == actual_index)
            {
              continue;
            }
            if( *itr_prob > max_blast_value)
            {
              max_blast_value = *itr_prob;
              max_blast_index = aa_blast_index;
            }
          }
          if( util::IsDefined( max_blast_index))
          {
            aa_index_window_pssm( i) = max_blast_index;
            aa_index_window_pssm_weight( i) = max_blast_value;
          }
        }
        // set undefined types to leucine, the most common amino acid, and get the index of all natural aas
        for( size_t i( 0); i < 3; ++i)
        {
          if( !aa_type_window( i).IsDefined())
          {
            aa_type_window( i) = biol::GetAATypes().LEU;
          }
          else if( !aa_type_window( i)->IsNaturalAminoAcid())
          {
            // try to get the parent type
            aa_type_window( i) = aa_type_window( i)->GetParentType();
            if( !aa_type_window( i)->IsNaturalAminoAcid())
            {
              // if the parent type is also unnatural, set to leucine
              aa_type_window( i) = biol::GetAATypes().LEU;
            }
          }
          aa_index_window( i) = aa_type_window( i).GetIndex();
        }
      }

      // determine the ids for each member of the triplet
      size_t triplet_id
      (
        aa_index_window( 0) + 20 * aa_index_window( 1) + 400 * aa_index_window( 2)
      );
      math::RunningAverage< linal::Vector< float> > average_result;
      average_result += m_TripletProbabilitiesPtr->operator()( triplet_id);

      for( size_t blast_position( 0), max_blast_pos( 3); blast_position < max_blast_pos; ++blast_position)
      {
        if( !util::IsDefined( aa_index_window_pssm( blast_position)))
        {
          continue;
        }
        std::swap( aa_index_window_pssm( blast_position), aa_index_window( blast_position));
        triplet_id = aa_index_window( 0) + 20 * aa_index_window( 1) + 400 * aa_index_window( 2);
        average_result.AddWeightedObservation
        (
          m_TripletProbabilitiesPtr->operator()( triplet_id),
          aa_index_window_pssm_weight( blast_position)
        );
        std::swap( aa_index_window_pssm( blast_position), aa_index_window( blast_position));
      }

      STORAGE.CopyValues( average_result.GetAverage());
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to forward calls from SetObject on to all internal implementations
    //! If InjectDimensions is set to true, then internal objects will also be called with SetDimension, whenever
    //! that function is called
    iterate::Generic< Base< biol::AABase, float> > AATripletLookup::GetInternalDescriptors()
    {
      if( m_EffTypeCalculator.IsDefined())
      {
        return iterate::Generic< Base< biol::AABase, float> >( &m_EffTypeCalculator, &m_EffTypeCalculator + 1);
      }
      return iterate::Generic< Base< biol::AABase, float> >();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AATripletLookup::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes lookups on an optionally smoothed aa triplet table"
      );
      parameters.AddInitializer
      (
        "filename",
        "file that contains the aa triplet HSC probabilities.  An example row in this file should look like: "
        "CPA 1.34 1.535 0.12 \nThe meaning of each number is shown on the annotated table below\n"
        "CPARK -- AA triplet \n"
        "0 -- # of times the C (in CPARK) is in a helix\n"
        "1 -- # of times the C (in CPARK) is in a strand\n"
        "5 -- # of times the C (in CPARK) is in a coil\n"
        "0 -- # of times the P (in CPARK) is in a helix\n"
        "and so forth for the remaining positions in the triplet.\n"
        "Note: This file will be overwritten if it is incomplete with a complete table, which will be computed using "
        "the transitions file",
        io::Serialization::GetAgentInputFilename( &m_SequencesFilename)
      );
      parameters.AddOptionalInitializer
      (
        "blast threshold",
        "number of steps which to smooth out SS types from triplet neighbors, using the transition matrices",
        io::Serialization::GetAgent( &m_BlastProbabilityThreshold)
      );
      parameters.AddOptionalInitializer
      (
        "blast type",
        "Blast or AA type descriptor to consider (must return 20 values per AA). Defaults to AABlastProbability",
        io::Serialization::GetAgent( &m_TypeCalculator)
      );
      return parameters;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    namespace
    {
      //! @brief a helper class that allows iteration between AA triplet hashes that are separated by a single AA
      class AATripletHashIterator
      {
      private:

        size_t m_Increment;   //!< Amount to increment during each ++
        size_t m_AATypeIndex; //!< Current AA type index
        size_t m_Hash;        //!< Current hash value
        static const size_t s_NumberAATypes = 20;

        //! @brief get the increment for a particular position in the hash
        static size_t GetIncrement( const size_t &POSITION)
        {
          static const size_t s_Increments[ 5] = { 1, 20, 400, 8000};
          return s_Increments[ POSITION];
        }

      public:

        //! @brief default constructor
        AATripletHashIterator() :
          m_Increment( 1),
          m_AATypeIndex( 0),
          m_Hash( 0)
        {
        }

        //! @brief HASH the hashed triplet value to consider
        //! @param POSITION the position (0-3) within the triplet to increment
        AATripletHashIterator( const size_t &HASH, const size_t &POSITION) :
          m_Increment( GetIncrement( POSITION)),
          m_AATypeIndex( ( HASH / m_Increment) % s_NumberAATypes),
          m_Hash( HASH)
        {
        }

        //! @brief get the current aa type index
        size_t GetAATypeIndex() const
        {
          return m_AATypeIndex;
        }

        //! @brief set the aa type index
        //! @param TYPE_INDEX the new type index
        void SetAATypeIndex( const size_t &TYPE_INDEX)
        {
          if( TYPE_INDEX > m_AATypeIndex)
          {
            m_Hash += ( TYPE_INDEX - m_AATypeIndex) * m_Increment;
          }
          else
          {
            m_Hash -= ( m_AATypeIndex - TYPE_INDEX) * m_Increment;
          }
          m_AATypeIndex = TYPE_INDEX;
        }

        //! @brief conversion to size_t (the current hash value)
        operator const size_t &() const
        {
          return m_Hash;
        }

        //! @brief increment operator
        AATripletHashIterator &operator++()
        {
          ++m_AATypeIndex;
          m_Hash += m_Increment;
          return *this;
        }

        //! @brief decrement operator
        AATripletHashIterator &operator--()
        {
          --m_AATypeIndex;
          m_Hash -= m_Increment;
          return *this;
        }

        //! @brief restart the iterator at the beginning of the sequence
        AATripletHashIterator &Restart()
        {
          m_Hash -= m_AATypeIndex * m_Increment;
          m_AATypeIndex = 0;
          return *this;
        }

        //! @brief get an iterator at the end of the sequence
        AATripletHashIterator End() const
        {
          AATripletHashIterator itr;
          itr.m_Hash = m_Hash + ( s_NumberAATypes - m_AATypeIndex) * m_Increment;
          itr.m_Increment = m_Increment;
          itr.m_AATypeIndex = s_NumberAATypes;
          return itr;
        }

        //! @brief equality operator
        bool operator ==( const AATripletHashIterator &OTHER) const
        {
          return m_Hash == OTHER.m_Hash;
        }

        //! @brief less than operator
        bool operator <( const AATripletHashIterator &OTHER) const
        {
          return m_Hash < OTHER.m_Hash;
        }

        //! @brief inequality operator
        bool operator !=( const AATripletHashIterator &OTHER) const
        {
          return m_Hash != OTHER.m_Hash;
        }
      };
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool AATripletLookup::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      if( m_TypeCalculator.IsDefined() && m_TypeCalculator->GetSizeOfFeatures() != size_t( 20))
      {
        ERR_STREAM << "A descriptor that returns 20 values is required for the AA type descriptor";
        return false;
      }
      m_EffTypeCalculator = m_TypeCalculator;
      if( !m_TypeCalculator.IsDefined() && m_BlastProbabilityThreshold < float( 1.0))
      {
        m_EffTypeCalculator = "AABlastProbability";
      }

      // lock the mutex before looking the map
      s_TripletProbabilitiesMutex.Lock();

      storage::Map< std::string, AATripletLookup::TripletSSPropensityType>::iterator
        itr_map( s_ProbabilitiesStorage.Find( m_SequencesFilename));

      if( itr_map != s_ProbabilitiesStorage.End())
      {
        s_TripletProbabilitiesMutex.Unlock();
        m_TripletProbabilitiesPtr = util::SiPtr< const TripletSSPropensityType>( itr_map->second);
        return true;
      }

      itr_map = s_ProbabilitiesStorage.Insert( std::make_pair( m_SequencesFilename, TripletSSPropensityType())).first;
      itr_map->second.Resize( 8000);
      LoadInitialPropensities( itr_map->second);
      m_TripletProbabilitiesPtr = util::SiPtr< TripletSSPropensityType>( itr_map->second);
      s_TripletProbabilitiesMutex.Unlock();

      return true;
    }

    //! @brief load the initial propensities vector
    //! @param TRIPLETS_PROB storage for all triplets loaded
    void AATripletLookup::LoadInitialPropensities( TripletSSPropensityType &TRIPLETS_PROB)
    {
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_SequencesFilename);

      static const storage::Vector< size_t> s_aa_index_from_char( CreateAATypeMap());

      // An example row in this file should look like:
      // CPARK 0 1 5 0 1 5 3 0 5 3 2 1 4 5 6
      // The meaning of each number is shown on the annotated table below
      // CPARK -- AA triplet
      // 0 -- # of times the C (in CPARK) is in a helix
      // 1 -- # of times the C (in CPARK) is in a strand
      // 5 -- # of times the C (in CPARK) is in a coil
      // 0 -- # of times the P (in CPARK) is in a helix
      // and so forth for the remaining positions in the triplet.
      // Floating point numbers can also be used for the values
      std::string sequence( size_t( 3), ' ');
      static const size_t s_radii[ 5] = { 1, 20, 400, 8000, 160000};
      while( input.good())
      {
        std::string line;
        std::getline( input, line);

        // test whether the first character is a space.  Normally it should not be
        if( !line.empty() && isspace( line[ 0]))
        {
          line.erase( 0, line.find_first_not_of( " \t\n"));
        }

        // comment line or empty line, both allowed
        if( line.size() < 10 || line[ 0] == '#')
        {
          continue;
        }

        // get the triplet index
        size_t triplet_index( 0);
        for( int i( 0); i < 3; ++i)
        {
          triplet_index += s_aa_index_from_char( int( line[ i])) * s_radii[ i];
        }

        // copy the first 4 characters into the sequence
        std::copy( line.begin(), line.begin() + 3, sequence.begin());

        // strip off the sequence from the line
        line.erase( 0, 4);

        // convert the line into the numbers
        storage::Vector< float> ss_counts( util::SplitStringToNumerical< float>( line, " \t"));

        // ensure the correct number of entries exist on the line
        BCL_Assert
        (
          ss_counts.GetSize() == size_t( 3),
          "Line should have contained 3 entries but instead contained: " + util::Format()( ss_counts.GetSize())
          + " line was: " + line
        );

        Propensity &propensity( TRIPLETS_PROB( triplet_index));
        BCL_Assert( propensity.Max() == propensity.Min(), "Repeated row for sequence: " + sequence);
        propensity( 0) = ss_counts( 0);
        propensity( 1) = ss_counts( 1);
        propensity( 2) = ss_counts( 2);
      }

      io::File::CloseClearFStream( input);
    }

  } // namespace descriptor
} // namespace bcl
