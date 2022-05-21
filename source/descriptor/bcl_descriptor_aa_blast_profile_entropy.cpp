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
#include "descriptor/bcl_descriptor_aa_blast_profile_entropy.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_running_average.h"
#include "score/bcl_score.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AABlastProfileEntropy::s_Instance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AABlastProfileEntropy( false)
      )
    );

    // Map from initializer string/type to the shptrvector of model interfaces; saves loading potentially
    // gigantic models repetitively
    storage::Map< std::string, AABlastProfileEntropy::t_AAEntropyStorage> AABlastProfileEntropy::s_EntropyMapStorage =
      storage::Map< std::string, AABlastProfileEntropy::t_AAEntropyStorage>();

    //  Mutex for access to s_EntropyMapStorage
    sched::Mutex AABlastProfileEntropy::s_EntropyMapMutex = sched::Mutex();

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @param INITIALIZE true if it is desired that this class be fully initialized; should not be set to true unless
    //!        using this class within the code
    AABlastProfileEntropy::AABlastProfileEntropy( const bool &INITIALIZE) :
      m_Filename( score::Score::AddHistogramPath( "blast_profile_entropies.txt")),
      m_MutantWeight( 1.0),
      m_NativeDominance( false),
      m_NumberEntropies( 6)
    {
      if( INITIALIZE)
      {
        ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
      }
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    AABlastProfileEntropy *AABlastProfileEntropy::Clone() const
    {
      return new AABlastProfileEntropy( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AABlastProfileEntropy::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AABlastProfileEntropy::GetAlias() const
    {
      static const std::string s_name( "AABlastProfileEntropy");
      return s_name;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief core function, computes the ss-specific blast entropy for the blast profile of the given AABase
    //! @param ELEMENT iterator to the AA of interest
    //! @return vector containing the values, one for each AA
    AABlastProfileEntropy::t_Entropy AABlastProfileEntropy::ComputeBlastSSEntropy
    (
      const iterate::Generic< const biol::AABase> &ELEMENT
    )
    {
      // compute the blast descriptor
      Iterator< biol::AABase> itr_blast( ELEMENT);
      linal::VectorConstReference< float> blast_profile( m_TypeCalculator->operator()( itr_blast));

      math::RunningAverage< t_Entropy> average_entropy;

      const int threshold_profile( std::min( blast_profile.Max(), float( 0.0)));
      t_AAEntropyStorage::const_iterator itr_entropy_storage( m_EntropiesPtr->Begin());

      const biol::AAType aa_type( ELEMENT->GetType());
      const bool is_likely_mutant
      (
        !aa_type->IsNaturalAminoAcid()
        || blast_profile( aa_type.GetIndex()) < threshold_profile
      );
      linal::VectorConstReference< float>::const_iterator
        itr_profile_native( blast_profile.Begin() + aa_type.GetIndex());

      // iterate through the AA types to find the dominant types
      for
      (
        linal::VectorConstReference< float>::const_iterator
          itr_profile( blast_profile.Begin()), itr_profile_end( blast_profile.End());
        itr_profile != itr_profile_end;
        ++itr_profile, ++itr_entropy_storage
      )
      {
        const int profile_value( *itr_profile);
        if( m_NativeDominance)
        {
          if( itr_profile != itr_profile_native)
          {
            continue;
          }
        }
        else if( profile_value < threshold_profile)
        {
          continue;
        }
        storage::Vector< t_BlastProfileBinEntropies>::const_iterator
          itr_entropy_target_storage( itr_entropy_storage->Begin());
        const float weight( profile_value > 5 ? 2.0 : 1.0);
        for
        (
          linal::VectorConstReference< float>::const_iterator itr_target( blast_profile.Begin());
          itr_target != itr_profile_end;
          ++itr_target, ++itr_entropy_target_storage
        )
        {
          const int target_value( *itr_profile);
          const int target_value_bin( GetBin( target_value));
          average_entropy.AddWeightedObservation
          (
            itr_entropy_target_storage->operator()( target_value_bin),
            weight
          );
        }
      }

      t_Entropy average( average_entropy.GetAverage());
      if( is_likely_mutant)
      {
        average *= m_MutantWeight;
      }
      if( average.GetSize() != m_NumberEntropies)
      {
        average = t_Entropy( m_NumberEntropies, util::GetUndefined< float>());
      }
      return average;
    }

    //! @brief calculate the descriptors
    //! @param ELEMENT: the element of interest
    //! @param STORAGE storage for the descriptor
    void AABlastProfileEntropy::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      STORAGE.CopyValues( ComputeBlastSSEntropy( ELEMENT));
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AABlastProfileEntropy::GetSerializer() const
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
        io::Serialization::GetAgentInputFilename( &m_Filename),
        score::Score::AddHistogramPath( "blast_profile_entropies.txt")
      );
      parameters.AddInitializer
      (
        "mutant multiplier",
        "multiplier to use if the blast profile suggests that the AA is a mutant, not normally present at the position",
        io::Serialization::GetAgentWithRange( &m_MutantWeight, 0.0, 1.0),
        "1.0"
      );
      parameters.AddInitializer
      (
        "native dominance",
        "true if the only dominant type to consider is the native",
        io::Serialization::GetAgent( &m_NativeDominance),
        "False"
      );
      parameters.AddInitializer
      (
        "blast type",
        "Blast or AA type descriptor to consider (must return 20 values per AA). Defaults to AABlastProfile",
        io::Serialization::GetAgent( &m_TypeCalculator),
        "AABlastProfile"
      );
      return parameters;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AABlastProfileEntropy::GetNormalSizeOfFeatures() const
    {
      return m_NumberEntropies;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to forward calls from SetObject on to all internal implementations
    //! If InjectDimensions is set to true, then internal objects will also be called with SetDimension, whenever
    //! that function is called
    iterate::Generic< Base< biol::AABase, float> > AABlastProfileEntropy::GetInternalDescriptors()
    {
      if( m_TypeCalculator.IsDefined())
      {
        return iterate::Generic< Base< biol::AABase, float> >( &m_TypeCalculator, &m_TypeCalculator + 1);
      }
      return iterate::Generic< Base< biol::AABase, float> >();
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool AABlastProfileEntropy::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      if( m_TypeCalculator.IsDefined() && m_TypeCalculator->GetSizeOfFeatures() != size_t( 20))
      {
        ERR_STREAM << "A descriptor that returns 20 values is required for the AA type descriptor";
        return false;
      }
      if( !m_TypeCalculator.IsDefined())
      {
        m_TypeCalculator = "AABlastProfile";
      }

      s_EntropyMapMutex.Lock();
      storage::Map< std::string, t_AAEntropyStorage>::const_iterator
        itr_entropy( s_EntropyMapStorage.Find( m_Filename));
      if( itr_entropy != s_EntropyMapStorage.End())
      {
        // file was already loaded and normalized, just use it directly
        m_EntropiesPtr = util::ToSiPtr( itr_entropy->second);
        m_NumberEntropies = m_EntropiesPtr->operator ()( 0)( 0)( 0).GetSize();
        s_EntropyMapMutex.Unlock();
        return true;
      }

      // create the object in the map with the correct size
      t_AAEntropyStorage &entropies( s_EntropyMapStorage[ m_Filename]);

      entropies.Resize
      (
        size_t( 20),
        storage::Vector< t_BlastProfileBinEntropies>( size_t( 20))
      );
      // attach the SiPtr to the newly-created matrices
      m_EntropiesPtr = util::ToSiPtr( entropies);

      // Load the file
      io::IFStream statistics_file;
      io::File::MustOpenIFStream( statistics_file, m_Filename);
      storage::Vector< storage::Vector< std::string> > tokenized_lines
      (
        util::SplittedStringLineListFromIStream( statistics_file)
      );
      io::File::CloseClearFStream( statistics_file);

      // get the number of columns to determine the proper bins
      BCL_Assert
      (
        tokenized_lines.GetSize() >= size_t( 5 * 20 * 20),
        "File contained insufficient lines to contain all statistics"
      );

      size_t number_data_lines( 0);
      m_NumberEntropies = 0;
      for
      (
        storage::Vector< storage::Vector< std::string> >::const_iterator
          itr_line( tokenized_lines.Begin()), itr_line_end( tokenized_lines.End());
        itr_line != itr_line_end;
        ++itr_line
      )
      {
        const storage::Vector< std::string> tokens( *itr_line);
        if( tokens.IsEmpty() || tokens( 0)[ 0] == '#')
        {
          // empty or comment line, continue
          continue;
        }
        ++number_data_lines;
        BCL_Assert
        (
          tokens.GetSize() > size_t( 3),
          "Incorrect number of entries on line; should be > 3, but was: " + util::Format()( tokens.GetSize())
        );
        BCL_Assert
        (
          tokens( 0).size() == size_t( 1) && tokens( 1).size() == size_t( 1) && tokens( 2).size() == size_t( 1),
          "First two tokens on each line should have just had an AA letter and 1-digit blast bin on them, not "
          + tokens( 0) + " " + tokens( 1) + " " + tokens( 2)
        );
        if( !m_NumberEntropies)
        {
          m_NumberEntropies = tokens.GetSize() - 3;
        }
        else
        {
          BCL_Assert( m_NumberEntropies + 3 == tokens.GetSize(), "Inconsistent # of entropies across file");
        }
        const size_t aa_dominant( biol::GetAATypes().AATypeFromOneLetterCode( tokens( 0)[ 0]).GetIndex());
        const size_t aa_target( biol::GetAATypes().AATypeFromOneLetterCode( tokens( 1)[ 0]).GetIndex());
        const size_t blast_bin( tokens( 2)[ 0] - '0');

        // get the entropy object for the designated dominant, target aa types and blast bin
        t_Entropy &entropy( entropies( aa_dominant)( aa_target)( blast_bin));
        entropy = linal::Vector< float>( m_NumberEntropies, float( 0.0));

        for( size_t entropy_n( 0); entropy_n < m_NumberEntropies; ++entropy_n)
        {
          entropy( 0) = util::ConvertStringToNumericalValue< float>( tokens( entropy_n + 3));
        }
      }

      BCL_Assert
      (
        number_data_lines == size_t( 5 * 20 * 20),
        "File contained wrong number of statistics lines, should be 2000"
      );

      s_EntropyMapMutex.Unlock();
      return true;
    }

    //! @brief get the bin for a particular blast profile value
    //! @param BLAST_PROFILE the blast profile value
    //! @return the bin for the blast profile value
    int AABlastProfileEntropy::GetBin( const int &BP_VALUE) const
    {
      // equivalent to (but faster than) max(0,min(4,(BP_VALUE + 9)/4))
      static const int s_bins[ 31] =
      {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, // -15 - -6
        1, 1, 1, 1, // -5 - -2
        2, 2, 2, 2, // -1 - 2
        3, 3, 3, 3, // 3 - 6
        4, 4, 4, 4, 4, 4, 4, 4, 4 // 7 - 15
      };
      return s_bins[ BP_VALUE + 15];
    }

  } // namespace descriptor
} // namespace bcl
