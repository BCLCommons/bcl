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
#include "descriptor/bcl_descriptor_aa_sse_info.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    //! @brief SSEInfoType as string
    //! @param SSE_INFO_TYPE the message level
    //! @return the SSEInfoType as string
    const std::string &AASSEInfo::GetSSEInfoTypeString( const SSEInfoType &SSE_INFO_TYPE)
    {
      static const std::string s_types[] =
      {
          "AA_SSESize", "AA_SSEID", "AA_PositionInSSE", "AA_DistanceFromSSECenter",
        GetStaticClassName< SSEInfoType>()
      };

      return s_types[ SSE_INFO_TYPE];
    }

    //! @brief SSEInfoType as string
    //! @param SSE_INFO_TYPE the message level
    //! @return the SSEInfoType as string
    const std::string &GetSSEInfoDescription( const AASSEInfo::SSEInfoType &SSE_INFO_TYPE)
    {
      static const std::string s_types[] =
      {
        "Size of the AA's SSE",
        "Index of the AA's SSE",
        "Normalized [0-1] AA's position in the SSE",
        "Normalized [0,1] AA's distance from the center of the SSE",
        GetStaticClassName< AASSEInfo::SSEInfoType>()
      };

      return s_types[ SSE_INFO_TYPE];
    }

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> AASSEInfo::s_AASSESizeInstance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AASSEInfo( e_AASSESize)
      )
    );
    const util::SiPtr< const util::ObjectInterface> AASSEInfo::s_AASSEIDInstance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AASSEInfo( e_AASSEID)
      )
    );
    const util::SiPtr< const util::ObjectInterface> AASSEInfo::s_AAPositionInSSEInstance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AASSEInfo( e_AAPositionInSSE)
      )
    );
    const util::SiPtr< const util::ObjectInterface> AASSEInfo::s_AADistanceFromSSECenterInstance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AASSEInfo( e_AADistanceFromSSECenter)
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AASSEInfo::AASSEInfo( const SSEInfoTypeEnum &SSE_INFO_TYPE_ENUM) :
      m_SSEInfoTypeEnum( SSE_INFO_TYPE_ENUM),
      m_Method( sspred::GetMethods().e_OCTOPUS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    AASSEInfo *AASSEInfo::Clone() const
    {
      return new AASSEInfo( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AASSEInfo::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AASSEInfo::GetAlias() const
    {
      return m_SSEInfoTypeEnum.GetString();
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AASSEInfo::GetNormalSizeOfFeatures() const
    {
      return 1;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT_A the element of interest
    //! @param STORAGE storage for the descriptor
    void AASSEInfo::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      // if the sequence maps are empty, this indicates that this object is now operating over a new sequence,
      // so it is necessary to reload the files
      if( m_SSESizeVec.IsEmpty())
      {
        LoadFiles();
      }

      const int seq_id( ELEMENT->GetData()->GetSeqID());
      storage::Map< int, size_t>::const_iterator map_itr( m_SeqIDToIndexMap.Find( seq_id));

      if( map_itr == m_SeqIDToIndexMap.End())
      {
        STORAGE( 0) = util::GetUndefined< float>();
        return;
      }

      const size_t vec_index( map_itr->second);

      switch( m_SSEInfoTypeEnum)
      {
        case e_AASSESize: // SSE Size
          STORAGE( 0) = m_SSESizeVec( vec_index);
          break;

        case e_AASSEID: // SSE ID
          STORAGE( 0) = m_SSEIDVec( vec_index);
          break;

        case e_AAPositionInSSE: // SSE Position
          STORAGE( 0) = m_PositionIndexVec( vec_index) / float( m_SSESizeVec( vec_index));
          break;

        case e_AADistanceFromSSECenter: // SSE Position from center
          STORAGE( 0) = 2.0 * std::abs( ( m_PositionIndexVec( vec_index) / float( m_SSESizeVec( vec_index)) - 0.5));
          break;

        default:
          STORAGE( 0) = util::GetUndefined< float>();
          break;
      }

      return;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AASSEInfo::SetObjectHook()
    {
      m_SeqIDToIndexMap.Reset();
      m_SSESizeVec.Reset();
      m_SSEIDVec.Reset();
      m_PositionIndexVec.Reset();
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference AASSEInfo::GetNormalCachePreference() const
    {
      return e_PreferCache;
    }

    //! @brief function to load files; should only be called the first time Calculate is called with a new sequence
    //! since the results are often in the cache
    void AASSEInfo::LoadFiles()
    {
      // Get protein model
      util::SiPtr< const assemble::ProteinModelWithCache> sp_protein_model( this->GetCurrentObject());

      // Check to make sure there is only a single chain in the given PDB file
      BCL_Assert
      (
        sp_protein_model->GetChains().GetSize() == 1,
        "More than one chain in the given protein model: "
         + util::Format()( sp_protein_model->GetProteinModelData()->GetData( assemble::ProteinModelData::e_PDBFile))
      );

      // Create a string to take all SSE single character predictions
      std::string ss_pred_str( "");

      // Store the highest seq ID seen
      size_t def_seq_ids_index( 0);

      util::ShPtr< biol::AASequence> sp_seq( sp_protein_model->GetChains().FirstElement()->GetSequence());

      BCL_MessageDbg( "Preparing to read in one state SS predictions...");

      // Iterate through the entire protein sequence checking seq IDs and SS single character predictions
      const bool method_requires_coordinates
      (
        m_Method == sspred::GetMethods().e_PDB
        || m_Method == sspred::GetMethods().e_PALSSE
        || m_Method == sspred::GetMethods().e_Stride
        || m_Method == sspred::GetMethods().e_DSSP
        || m_Method == sspred::GetMethods().e_StrideDSSP
        || m_Method == sspred::GetMethods().e_MAHSSMI
      );
      for( biol::AASequence::const_iterator itr( sp_seq->Begin()), itr_end( sp_seq->End()); itr < itr_end; ++itr)
      {
        const int temp_seq_id( ( *itr)->GetData()->GetSeqID());
        if( util::IsDefined( temp_seq_id))
        {
          if( !( *itr)->GetCenter().IsDefined() && method_requires_coordinates)
          {
            continue;
          }
          ss_pred_str += ( *itr)->GetSSPrediction( m_Method)->GetOneStateSSPrediction()->GetOneLetterCode();

          // Add entry to map
          m_SeqIDToIndexMap.Insert( storage::Pair< int, size_t>( temp_seq_id, def_seq_ids_index));

          ++def_seq_ids_index;
        }
      }
      const size_t max_vec_size( def_seq_ids_index + 1);

      BCL_MessageDbg( "Done reading in one state SS predictions");

      // Set vector sizes based on max seq ID
      m_SSESizeVec.AllocateMemory( max_vec_size);
      m_SSEIDVec.AllocateMemory( max_vec_size);
      m_PositionIndexVec.AllocateMemory( max_vec_size);

      // Keep track of current index position in string
      size_t index( 0);
      size_t sse_id( 0);

      // Iterate through string setting maximum SSE sizes
      for
      (
        std::string::const_iterator itr( ss_pred_str.begin()), itr_end( ss_pred_str.end());
        itr < itr_end;
        ++sse_id
      )
      {
        size_t next_index( std::min( ss_pred_str.size(), ss_pred_str.find_first_not_of( *itr, index)));
        size_t temp_size( next_index - index);

        if( m_SSEInfoTypeEnum == e_AAPositionInSSE || m_SSEInfoTypeEnum == e_AADistanceFromSSECenter)
        {
          // the desired size is 1 less, due to the normalization
          temp_size = std::max( temp_size, size_t( 2)) - 1;
        }

        // Iterate internally until the current SSE is finished
        for( size_t pos( 0); index < next_index; ++itr, ++index, ++pos)
        {
          m_SSESizeVec.PushBack( temp_size);
          m_SSEIDVec.PushBack( sse_id);
          m_PositionIndexVec.PushBack( pos);
        }
      }

      // DEBUG print out the final vectors so I can check them manually!
      BCL_MessageDbg( "SSE Sizes vector: " + util::Format()( m_SSESizeVec));
      BCL_MessageDbg( "SSE ID vector: " + util::Format()( m_SSEIDVec));
      BCL_MessageDbg( "AA in SSE position index vector: " + util::Format()( m_PositionIndexVec));
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AASSEInfo::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        GetSSEInfoDescription( m_SSEInfoTypeEnum) + " based on predictions or actual secondary structure. "
      );
      parameters.AddInitializer
      (
        "method",
        "secondary structure prediction or analysis method to use to obtain result",
        io::Serialization::GetAgent( &m_Method),
        "OCTOPUS"
      );

      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
