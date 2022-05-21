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
#include "descriptor/bcl_descriptor_mutation_nearest_in_sequence.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MutationNearestInSequence::s_Instance
    (
      util::Enumerated< Base< biol::Mutation, float> >::AddInstance
      (
        new MutationNearestInSequence()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutationNearestInSequence::MutationNearestInSequence() :
      m_IncludeDistances( false),
      m_IncludeDirection( false),
      m_MaxToFind( 1),
      m_SkipUndefined( false)
    {
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    MutationNearestInSequence *MutationNearestInSequence::Clone() const
    {
      return new MutationNearestInSequence( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutationNearestInSequence::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &MutationNearestInSequence::GetAlias() const
    {
      static const std::string s_name( "NearestMutants");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t MutationNearestInSequence::GetNormalSizeOfFeatures() const
    {
      return m_MaxToFind * ( m_Descriptor->GetSizeOfFeatures() + ( m_IncludeDirection ? 1 : 0) + ( m_IncludeDistances ? 1 : 0));
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< biol::Mutation, float> > MutationNearestInSequence::GetInternalDescriptors()
    {
      return iterate::Generic< Base< biol::Mutation, float> >( &m_Descriptor, &m_Descriptor + 1);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    void MutationNearestInSequence::RecalculateImpl
    (
      const Iterator< biol::Mutation> &ITR,
      linal::VectorReference< float> &STORAGE
    )
    {
      // if we have not yet located the AAs for this protein, do so now
      if( m_Mutations.IsEmpty())
      {
        FindMutants( ITR);
      }

      const int seq_id( ITR( 0)->GetResidueNumber());
      const size_t n_features( m_Descriptor->GetSizeOfFeatures());
      auto set_itr_next( m_Mutations.Find( *ITR( 0)));
      auto set_begin( m_Mutations.Begin());
      auto set_end( m_Mutations.End());

      auto original_set_itr( set_itr_next);
      ++set_itr_next;
      // if the mutation is the first one, there is no need to look backwards
      size_t pos( 0), n( 0);
      auto set_itr_prev( original_set_itr);
      bool can_go_down( original_set_itr != set_begin), can_go_up( set_itr_next != set_end);
      if( can_go_down)
      {
        --set_itr_prev;
      }
      for( ; n < m_MaxToFind; ++n)
      {
        bool go_up( false);
        if( can_go_up)
        {
          if( can_go_down)
          {
            go_up = set_itr_next->first.GetResidueNumber() - original_set_itr->first.GetResidueNumber()
                     < original_set_itr->first.GetResidueNumber() - set_itr_prev->first.GetResidueNumber();
          }
          else
          {
            go_up = true;
          }
        }
        else if( !can_go_down)
        {
          break;
        }
        if( go_up)
        {
          if( m_IncludeDistances)
          {
            STORAGE( pos++) = set_itr_next->first.GetResidueNumber() - original_set_itr->first.GetResidueNumber();
          }
          if( m_IncludeDirection)
          {
            STORAGE( pos++) = set_itr_next->first.GetResidueNumber() == original_set_itr->first.GetResidueNumber() ? 0.0 : 1.0;
          }
          const linal::VectorConstReference< float> result( ( *m_Descriptor)( set_itr_next->second));
          if( m_SkipUndefined && !result.IsDefined())
          {
            if( m_IncludeDistances)
            {
              --pos;
            }
            if( m_IncludeDirection)
            {
              --pos;
            }
            --n;
          }
          else
          {
            BCL_MessageStd
            (
              "using mutant " + set_itr_next->second( 0)->ToString() +
              "(next) to describe " + ITR( 0)->ToString()
            );
            STORAGE.CreateSubVectorReference( n_features, pos).CopyValues( result);
            pos += n_features;
          }
          ++set_itr_next;
          can_go_up = set_itr_next != set_end;
        }
        else // if( go_down)
        {
          if( m_IncludeDistances)
          {
            STORAGE( pos++) = original_set_itr->first.GetResidueNumber() - set_itr_prev->first.GetResidueNumber();
          }
          if( m_IncludeDirection)
          {
            STORAGE( pos++) = set_itr_prev->first.GetResidueNumber() == original_set_itr->first.GetResidueNumber() ? 0.0 : -1.0;
          }

          const linal::VectorConstReference< float> result( ( *m_Descriptor)( set_itr_prev->second));
          if( m_SkipUndefined && !result.IsDefined())
          {
            if( m_IncludeDistances)
            {
              --pos;
            }
            if( m_IncludeDirection)
            {
              --pos;
            }
            --n;
          }
          else
          {
            BCL_MessageStd
            (
              "using mutant " + set_itr_prev->second( 0)->ToString() +
              "(prev) to describe " + ITR( 0)->ToString()
            );
            STORAGE.CreateSubVectorReference( n_features, pos).CopyValues( result);
            pos += n_features;
          }
          if( set_itr_prev == set_begin)
          {
            can_go_down = false;
          }
          else
          {
            --set_itr_prev;
          }
        }
      }

      return;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void MutationNearestInSequence::SetObjectHook()
    {
      m_Mutations.Reset();
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference MutationNearestInSequence::GetNormalCachePreference() const
    {
      return e_PreferCache;
    }

    //! @brief function to find the aas of the given type in the specified sequence
    void MutationNearestInSequence::FindMutants( const Iterator< biol::Mutation> &ITR)
    {
      for( Iterator< biol::Mutation> itr( ITR.Begin()->Begin()); itr.NotAtEnd(); ++itr)
      {
        m_Mutations.Insert( std::make_pair( *( itr( 0)), itr));
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutationNearestInSequence::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "retrieve properties of the nearest mutations to the current one");

      parameters.AddInitializer
      (
        "distance",
        "whether to write the direction before the descriptor and direction, if requested)",
        io::Serialization::GetAgent( &m_IncludeDistances),
        "False"
      );

      parameters.AddInitializer
      (
        "direction",
        "whether to write the direction (before the descriptors value and after the distance, if requested)",
        io::Serialization::GetAgent( &m_IncludeDirection),
        "False"
      );

      parameters.AddInitializer
      (
        "count",
        "# of nearby mutations to return",
        io::Serialization::GetAgentWithMin( &m_MaxToFind, size_t( 1)),
        "1"
      );

      parameters.AddInitializer
      (
        "descriptor",
        "Descriptor to retrieve for the nearest mutations",
        io::Serialization::GetAgent( &m_Descriptor)
      );

      parameters.AddInitializer
      (
        "skip undefined",
        "set this to true to skip descriptors that are undefined",
        io::Serialization::GetAgent( &m_SkipUndefined),
        "False"
      );

      return parameters;
    }

    //! @brief return the type of symmetry this descriptor has
    //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
    Type::Symmetry MutationNearestInSequence::GetSymmetry() const
    {
      return m_Descriptor->GetType().GetSymmetry();
    }

    //! @brief return whether this descriptor is valid if repeated elements are given
    //! @return true if this descriptor is valid if repeated elements are given
    //! This will be the case if the descriptor may have a legitimate value for A-A
    bool MutationNearestInSequence::ConsiderRepeatedElements() const
    {
      return m_Descriptor->ConsiderRepeatedElements();
    }

    //! @brief get the normal dimension for this descriptor
    //! @return the normal dimension for this descriptor
    size_t MutationNearestInSequence::GetNormalDimension() const
    {
      return m_Descriptor->GetType().GetDimension();
    }

  } // namespace descriptor
} // namespace bcl
