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
#include "descriptor/bcl_descriptor_aa_distance_from_aa_type.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> AADistanceFromAAType::s_BlastInstance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AADistanceFromAAType( true)
      )
    );
    const util::SiPtr< const util::ObjectInterface> AADistanceFromAAType::s_Instance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AADistanceFromAAType( false)
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AADistanceFromAAType::AADistanceFromAAType( const bool &USE_BLAST_PROBABILITY) :
      m_HaveLocated( false),
      m_UseBlastProbability( USE_BLAST_PROBABILITY),
      m_MaxDistance( 10),
      m_MaxToFind( 1)
    {
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    AADistanceFromAAType *AADistanceFromAAType::Clone() const
    {
      return new AADistanceFromAAType( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AADistanceFromAAType::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AADistanceFromAAType::GetAlias() const
    {
      static const std::string s_name( "AA_DistanceFromAAType"), s_blast_name( "AA_DistanceFromAABlastProfileType");
      return m_UseBlastProbability ? s_blast_name : s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AADistanceFromAAType::GetNormalSizeOfFeatures() const
    {
      return m_MaxToFind;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT_A the element of interest
    //! @param STORAGE storage for the descriptor
    void AADistanceFromAAType::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // if we have not yet located the AAs for this protein, do so now
      if( !m_HaveLocated)
      {
        FindAAs();
      }

      STORAGE = m_MaxDistance + 1;

      // if there are no AAs of the specified type, return
      if( m_SeqIDsOfType.IsEmpty())
      {
        return;
      }

      const int seq_id( ELEMENT->GetData()->GetSeqID());
      storage::Set< int>::const_iterator set_itr( m_SeqIDsOfType.LowerBound( seq_id));
      storage::Set< int>::const_iterator set_begin( m_SeqIDsOfType.Begin());
      storage::Set< int>::const_iterator set_end( m_SeqIDsOfType.End());
      if( set_itr == m_SeqIDsOfType.End())
      {
        --set_itr;
      }
      storage::Set< int>::const_iterator original_set_itr( set_itr);

      // center
      storage::Vector< size_t> distances;

      size_t n_found_right( 0), n_found_left( 0);
      if( m_Direction != e_Right)
      {
        // reverse until the first amino acid that is to the left or equal to the current AA
        if( *set_itr > seq_id)
        {
          if( set_itr != set_begin)
          {
            --set_itr;
          }
        }

        // find the first m_MaxToFind instances to the left of this instance
        while( 1)
        {
          const size_t seq_distance( seq_id - *set_itr);
          if( seq_distance > m_MaxDistance)
          {
            break;
          }
          distances.PushBack( seq_distance);
          if( set_itr == set_begin || ++n_found_left == m_MaxToFind)
          {
            break;
          }
          --set_itr;
        }
      }
      set_itr = original_set_itr;
      if( m_Direction != e_Left)
      {
        // forward until the first amino acid that is to the right of the current AA
        if( *set_itr <= seq_id)
        {
          ++set_itr;
        }

        // find the first m_MaxToFind instances to the left of this instance
        while( set_itr != set_end && n_found_right < m_MaxToFind)
        {
          const size_t seq_distance( *set_itr - seq_id);
          if( seq_distance > m_MaxDistance)
          {
            break;
          }
          distances.PushBack( seq_distance);
          ++n_found_right;
          ++set_itr;
        }
      }
      if( n_found_left && n_found_right)
      {
        distances.Sort( std::less< size_t>());
      }

      // put the distances into the storage
      for( size_t dist( 0), max_dist( std::min( m_MaxToFind, distances.GetSize())); dist < max_dist; ++dist)
      {
        STORAGE( dist) = distances( dist);
      }

      return;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AADistanceFromAAType::SetObjectHook()
    {
      m_SeqIDsOfType.Reset();
      m_HaveLocated = false;
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference AADistanceFromAAType::GetNormalCachePreference() const
    {
      return e_PreferCache;
    }

    //! @brief function to find the aas of the given type in the specified sequence
    void AADistanceFromAAType::FindAAs()
    {
      // Get protein model
      util::SiPtr< const assemble::ProteinModelWithCache> sp_protein_model( this->GetCurrentObject());

      // get an iterator on it
      iterate::Generic< const biol::AABase> aa_itr( sp_protein_model->GetIterator());

      // find all the AAs of the specified type
      if( m_UseBlastProbability)
      {
        for( ; aa_itr.NotAtEnd(); ++aa_itr)
        {
          if( aa_itr->GetBlastProfile().GetProbabilities()( m_Type.GetIndex()) > 0.05)
          {
            m_SeqIDsOfType.Insert( aa_itr->GetSeqID());
          }
        }
      }
      else
      {
        for( ; aa_itr.NotAtEnd(); ++aa_itr)
        {
          if( aa_itr->GetType() == m_Type)
          {
            m_SeqIDsOfType.Insert( aa_itr->GetSeqID());
          }
        }
      }

      // set the bool
      m_HaveLocated = true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AADistanceFromAAType::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        m_UseBlastProbability
        ? "The sequence distance to an amino acid whose blast probability is at least 0.05"
        : "The sequence distance to an amino acid with the specified type"
      );

      parameters.AddInitializer
      (
        "type",
        "type of amino acid to find",
        io::Serialization::GetAgent( &m_Type)
      );

      parameters.AddInitializer
      (
        "direction",
        "directions in which to search for the aa",
        io::Serialization::GetAgent( &m_Direction),
        "Center"
      );

      parameters.AddInitializer
      (
        "max",
        "maximum distance to search; if the type is not found, returns this # + 1",
        io::Serialization::GetAgentWithMin( &m_MaxDistance, size_t( 1))
      );

      parameters.AddInitializer
      (
        "size",
        "# of AAs of the specified type to locate; e.g. 2 to find the nearest two matches in the given direction",
        io::Serialization::GetAgentWithMin( &m_MaxToFind, size_t( 1)),
        "1"
      );

      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
