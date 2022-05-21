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
#include "assemble/bcl_assemble_voxel_grid_mutation.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "descriptor/bcl_descriptor_mutation_nearest_spatially.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MutationNearestSpatially::s_Instance
    (
      util::Enumerated< Base< biol::Mutation, float> >::AddInstance
      (
        new MutationNearestSpatially()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutationNearestSpatially::MutationNearestSpatially() :
      m_IncludeDistances( false),
      m_IncludeSeqDistances( false),
      m_IncludeDirection( false),
      m_IncludeSeqDirection( false),
      m_SkipUndefined( false),
      m_MaxToFind( 1)
    {
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    MutationNearestSpatially *MutationNearestSpatially::Clone() const
    {
      return new MutationNearestSpatially( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutationNearestSpatially::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &MutationNearestSpatially::GetAlias() const
    {
      static const std::string s_name( "NearestMutantsSpatially");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t MutationNearestSpatially::GetNormalSizeOfFeatures() const
    {
      return m_MaxToFind *
             (
               m_Descriptor->GetSizeOfFeatures()
               + ( m_IncludeDirection ? 3 : 0)
               + ( m_IncludeDistances ? 1 : 0)
               + ( m_IncludeSeqDirection ? 1 : 0)
               + ( m_IncludeSeqDistances ? 1 : 0)
             );
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< biol::Mutation, float> > MutationNearestSpatially::GetInternalDescriptors()
    {
      return iterate::Generic< Base< biol::Mutation, float> >( &m_Descriptor, &m_Descriptor + 1);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    void MutationNearestSpatially::RecalculateImpl
    (
      const Iterator< biol::Mutation> &ITR,
      linal::VectorReference< float> &STORAGE
    )
    {
      // if we have not yet located the AAs for this protein, do so now
      storage::Vector< storage::Triplet< double, Iterator< biol::Mutation>, util::SiPtr< const linal::Vector3D> > > vec_dist_iterator;
      vec_dist_iterator.AllocateMemory( m_MaxToFind);
      assemble::VoxelGridMutation vgm( 4.0, false, true, false);
      auto positions( vgm.ExtractPositions( *ITR( 0)));
      for( auto itr_pos( positions.Begin()), itr_pos_end( positions.End()); itr_pos != itr_pos_end; ++itr_pos)
      {
        for( Iterator< biol::Mutation> itr_restart( ITR.Begin()->Begin()); itr_restart.NotAtEnd(); ++itr_restart)
        {
          if( itr_restart.GetPosition() == ITR.GetPosition())
          {
            continue;
          }
          const double sqdist
          (
            linal::SquareDistance
            (
              **itr_pos,
              itr_restart( 0)->GetAAs()( 0)->GetCA().GetCenter()
            )
          );
          if( !util::IsDefined( sqdist))
          {
            continue;
          }
          if( vec_dist_iterator.GetSize() < m_MaxToFind || sqdist < vec_dist_iterator( m_MaxToFind - 1).First())
          {
            if( m_SkipUndefined && !( *m_Descriptor)( itr_restart).IsDefined())
            {
              continue;
            }
            if( vec_dist_iterator.GetSize() < m_MaxToFind)
            {
              vec_dist_iterator.PushBack
              (
                storage::Triplet< double, Iterator< biol::Mutation>, util::SiPtr< const linal::Vector3D> >
                (
                  sqdist,
                  itr_restart,
                  *itr_pos
                )
              );
              if( vec_dist_iterator.GetSize() == m_MaxToFind)
              {
                vec_dist_iterator.Sort( std::less< storage::Triplet< double, Iterator< biol::Mutation>, util::SiPtr< const linal::Vector3D> > >());
              }
            }
            else if( sqdist < vec_dist_iterator( m_MaxToFind - 1).First())
            {
              size_t lpos( m_MaxToFind - 1);
              while( lpos && sqdist < vec_dist_iterator( lpos - 1).First())
              {
                vec_dist_iterator( lpos) = vec_dist_iterator( lpos - 1);
                --lpos;
              }
              vec_dist_iterator( lpos) =
                storage::Triplet< double, Iterator< biol::Mutation>, util::SiPtr< const linal::Vector3D> >( sqdist, itr_restart, *itr_pos);
            }
          }
        }
      }

      const size_t n_features( m_Descriptor->GetSizeOfFeatures());
      size_t pos( 0);
      for
      (
        auto itr_closest( vec_dist_iterator.Begin()), itr_closest_end( vec_dist_iterator.End());
        itr_closest != itr_closest_end;
        ++itr_closest
      )
      {
        if( m_IncludeDistances)
        {
          STORAGE( pos++) = math::Sqrt( std::max( itr_closest->First(), 0.0));
        }
        if( m_IncludeSeqDistances)
        {
          STORAGE( pos++) = math::Absolute( ITR( 0)->GetResidueNumber() - itr_closest->Second()( 0)->GetResidueNumber());
        }
        if( m_IncludeDirection)
        {
          if( itr_closest->Second()( 0)->GetResidueNumber() == ITR( 0)->GetResidueNumber())
          {
            STORAGE( pos++) = 0.0;
            STORAGE( pos++) = 0.0;
            STORAGE( pos++) = 0.0;
          }
          else
          {
            linal::Vector3D direction
            (
              linal::UnitVector
              (
                itr_closest->Second()( 0)->GetAAs()( 0)->GetCenter(),
                *itr_closest->Third()
              )
            );
            STORAGE( pos++) = direction( 0);
            STORAGE( pos++) = direction( 1);
            STORAGE( pos++) = direction( 2);
          }
        }
        if( m_IncludeSeqDirection)
        {
          STORAGE( pos++)
            = ITR( 0)->GetResidueNumber() == itr_closest->Second()( 0)->GetResidueNumber()
              ? 0.0
              : ITR( 0)->GetResidueNumber() < itr_closest->Second()( 0)->GetResidueNumber()
                ? 1.0
                : -1.0;
        }
        const linal::VectorConstReference< float> result( ( *m_Descriptor)( itr_closest->Second()));
        STORAGE.CreateSubVectorReference( n_features, pos).CopyValues( result);
        pos += n_features;
      }

      return;
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference MutationNearestSpatially::GetNormalCachePreference() const
    {
      return e_PreferCache;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutationNearestSpatially::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "retrieve properties of the nearest mutations to the current one");

      parameters.AddInitializer
      (
        "distance",
        "whether to write the distance before the descriptor and direction, if requested)",
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
        "seqdistance",
        "whether to write the sequence distance before the descriptor and direction, if requested)",
        io::Serialization::GetAgent( &m_IncludeSeqDistances),
        "False"
      );

      parameters.AddInitializer
      (
        "seqdirection",
        "whether to write the sequence direction (before the descriptors value and after the distance, if requested)",
        io::Serialization::GetAgent( &m_IncludeSeqDirection),
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
        "skip undefined",
        "Set to true to skip nan values if the descriptor returns them. The default behavior would be for this "
        "descriptor to forward the NaN, which would result in the feature being skipped",
        io::Serialization::GetAgent( &m_SkipUndefined),
        "False"
      );

      parameters.AddInitializer
      (
        "descriptor",
        "Descriptor to retrieve for the nearest mutations",
        io::Serialization::GetAgent( &m_Descriptor)
      );

      return parameters;
    }

    //! @brief return the type of symmetry this descriptor has
    //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
    Type::Symmetry MutationNearestSpatially::GetSymmetry() const
    {
      return m_Descriptor->GetType().GetSymmetry();
    }

    //! @brief return whether this descriptor is valid if repeated elements are given
    //! @return true if this descriptor is valid if repeated elements are given
    //! This will be the case if the descriptor may have a legitimate value for A-A
    bool MutationNearestSpatially::ConsiderRepeatedElements() const
    {
      return m_Descriptor->ConsiderRepeatedElements();
    }

    //! @brief get the normal dimension for this descriptor
    //! @return the normal dimension for this descriptor
    size_t MutationNearestSpatially::GetNormalDimension() const
    {
      return m_Descriptor->GetType().GetDimension();
    }

  } // namespace descriptor
} // namespace bcl
