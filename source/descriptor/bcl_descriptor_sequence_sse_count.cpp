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
#include "descriptor/bcl_descriptor_sequence_sse_count.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> SequenceSSECount::s_SSECountsInstance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new SequenceSSECount( false, false)
      )
    );
    const util::SiPtr< const util::ObjectInterface> SequenceSSECount::s_MembraneSSECountsInstance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new SequenceSSECount( true, false)
      )
    );
    const util::SiPtr< const util::ObjectInterface> SequenceSSECount::s_SolutionSSECountsInstance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new SequenceSSECount( false, true)
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SequenceSSECount::SequenceSSECount( const bool &ONLY_MEMBRANE_COUNTS, const bool &ONLY_SOLUTION_COUNTS) :
      m_MinHelixSize( 0),
      m_MinStrandSize( 0),
      m_MinCoilSize( 0),
      m_OnlyMembraneCounts( ONLY_MEMBRANE_COUNTS),
      m_OnlySolutionCounts( ONLY_SOLUTION_COUNTS)
    {
      if( m_OnlyMembraneCounts && m_OnlySolutionCounts)
      {
        // if both only parameters are set to true, set them both to false since that just means ignore environment
        m_OnlyMembraneCounts = m_OnlySolutionCounts = false;
      }
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    SequenceSSECount *SequenceSSECount::Clone() const
    {
      return new SequenceSSECount( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SequenceSSECount::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &SequenceSSECount::GetAlias() const
    {
      static const std::string s_name( "HelixStrandCoilCounts"),
                               s_membrane_name( "HelixStrandCoilCountsInMembrane"),
                               s_solution_name( "HelixStrandCoilCountsInSolution");
      return m_OnlyMembraneCounts ? s_membrane_name : m_OnlySolutionCounts ? s_solution_name : s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t SequenceSSECount::GetNormalSizeOfFeatures() const
    {
      return 3;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    void SequenceSSECount::Calculate
    (
      linal::VectorReference< float> &STORAGE
    )
    {
      if( m_OnlyMembraneCounts || m_OnlySolutionCounts)
      {
        const biol::EnvironmentType required_env
        (
          m_OnlyMembraneCounts
          ? biol::GetEnvironmentTypes().e_MembraneCore
          : biol::GetEnvironmentTypes().e_Solution
        );

        // iterate over all residues
        for
        (
          iterate::Generic< const biol::AABase> itr( this->GetCurrentObject()->GetIterator());
          itr.NotAtEnd();
          // iteration in loop
        )
        {
          // access the environment type
          util::SiPtr< const sspred::MethodInterface> si_env( itr->GetSSPrediction( m_EnvironmentMethod));

          // check whether the environment type is defined and the desired value
          if( !si_env.IsDefined() || si_env->GetOneStateTMPrediction() != required_env)
          {
            ++itr;
            continue;
          }

          // access the ss-type
          util::SiPtr< const sspred::MethodInterface> si_method( itr->GetSSPrediction( m_SSMethod));

          // handle undefined predictions
          if( !si_method.IsDefined())
          {
            ++itr;
            continue;
          }

          // record the starting pdb id
          const int pdb_id_start( itr->GetPdbID());

          // record the starting ss type
          const biol::SSType ss_type( si_method->GetOneStateSSPrediction());

          // initialize the ending pdb id
          int pdb_id_end( pdb_id_start);

          // continue until the end of the sse or environment (possibly also the end of the protein)
          for( ++itr; itr.NotAtEnd(); ++itr)
          {
            // access the ss-type
            util::SiPtr< const sspred::MethodInterface> si_next_method( itr->GetSSPrediction( m_SSMethod));

            // handle undefined predictions and change in sse
            if( !si_next_method.IsDefined() || si_next_method->GetOneStateSSPrediction() != ss_type)
            {
              break;
            }

            // access the env-type
            util::SiPtr< const sspred::MethodInterface> si_next_env( itr->GetSSPrediction( m_EnvironmentMethod));

            // handle undefined predictions and change in sse
            if( !si_next_env.IsDefined() || si_next_env->GetOneStateTMPrediction() != required_env)
            {
              break;
            }

            // update the pdb_id of the end; since we are still on the same sse
            pdb_id_end = itr->GetPdbID();
          }

          // check sizes
          const size_t sse_size( pdb_id_end - pdb_id_start + 1);
          if( sse_size >= GetMinSSESize( ss_type))
          {
            // min size met; add 1 to the storage vector
            STORAGE( ss_type->GetIndex() <= size_t( 2) ? ss_type->GetIndex() : size_t( 2)) += 1.0;
          }
        }
      }
      else
      {
        // iterate over all residues
        for
        (
          iterate::Generic< const biol::AABase> itr( this->GetCurrentObject()->GetIterator());
          itr.NotAtEnd();
          // iteration in loop
        )
        {
          // access the ss-type
          util::SiPtr< const sspred::MethodInterface> si_method( itr->GetSSPrediction( m_SSMethod));

          // handle undefined predictions
          if( !si_method.IsDefined())
          {
            ++itr;
            continue;
          }

          // record the starting pdb id
          const int pdb_id_start( itr->GetPdbID());

          // record the starting ss type
          const biol::SSType ss_type( si_method->GetOneStateSSPrediction());

          // initialize the ending pdb id
          int pdb_id_end( pdb_id_start);

          // continue until the end of the sse (possibly also the end of the protein)
          for( ++itr; itr.NotAtEnd(); ++itr)
          {
            // access the ss-type
            util::SiPtr< const sspred::MethodInterface> si_next_method( itr->GetSSPrediction( m_SSMethod));

            // handle undefined predictions
            if( !si_next_method.IsDefined() || si_next_method->GetOneStateSSPrediction() != ss_type)
            {
              break;
            }

            // update the pdb_id of the end; since we are still on the same sse
            pdb_id_end = itr->GetPdbID();
          }

          // check sizes
          const size_t sse_size( pdb_id_end - pdb_id_start + 1);
          if( sse_size >= GetMinSSESize( ss_type))
          {
            // min size met; add 1 to the storage vector
            STORAGE( ss_type->GetIndex() <= size_t( 2) ? ss_type->GetIndex() : size_t( 2)) += 1.0;
          }
        }
      }
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference SequenceSSECount::GetNormalCachePreference() const
    {
      return e_PreferCache;
    }

    //! @brief get the minimum sse size for a particular ss-type
    //! @param TYPE the type of sse to get the minimum sse size for
    size_t SequenceSSECount::GetMinSSESize( const biol::SSType &TYPE) const
    {
      return TYPE == biol::GetSSTypes().COIL
             ? m_MinCoilSize
             : TYPE == biol::GetSSTypes().STRAND
               ? m_MinStrandSize
               : m_MinHelixSize;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SequenceSSECount::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Counts of each each secondary structure element (helices, strands, and coils) "
        + std::string( m_OnlyMembraneCounts ? " in the membrane" : m_OnlySolutionCounts ? " in solution" : "")
      );
      parameters.AddInitializer
      (
        "method",
        "secondary structure prediction or analysis method for which counts should be made",
        io::Serialization::GetAgent( &m_SSMethod)
      );
      if( m_OnlyMembraneCounts || m_OnlySolutionCounts)
      {
        parameters.AddInitializer
        (
          "membrane method",
          "membrane environment prediction or analysis method that should define the membrane",
          io::Serialization::GetAgent( &m_EnvironmentMethod)
        );
      }
      parameters.AddInitializer
      (
        "min helix size",
        "minimum number of residues in a helix to require before declaring it a helix",
        io::Serialization::GetAgent( &m_MinHelixSize),
        "0"
      );
      parameters.AddInitializer
      (
        "min strand size",
        "minimum number of residues in a helix to require before declaring it a helix",
        io::Serialization::GetAgent( &m_MinStrandSize),
        "0"
      );
      parameters.AddInitializer
      (
        "min coil size",
        "minimum number of residues in a helix to require before declaring it a helix",
        io::Serialization::GetAgent( &m_MinCoilSize),
        "0"
      );

      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
