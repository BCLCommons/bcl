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
#include "descriptor/bcl_descriptor_aa_hbond_neighbor.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_aa_dssp_info.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> AAHbondNeighbor::s_Instance
    (
      util::Enumerated< Base< biol::AABase, float> >::AddInstance
      (
        new AAHbondNeighbor()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor, can take a descriptor to take for the hbonded neighbor
    //! @param IMPL the descriptor of interest
    AAHbondNeighbor::AAHbondNeighbor( const util::Implementation< Base< biol::AABase, float> > &IMPL) :
      m_Descriptor( IMPL)
    {
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    AAHbondNeighbor *AAHbondNeighbor::Clone() const
    {
      return new AAHbondNeighbor( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AAHbondNeighbor::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AAHbondNeighbor::GetAlias() const
    {
      static const std::string s_name( "AAHbondNeighbor");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AAHbondNeighbor::GetNormalSizeOfFeatures() const
    {
      return m_Descriptor->GetNormalSizeOfFeatures();
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    void AAHbondNeighbor::RecalculateImpl( const Iterator< biol::AABase> &ITR, linal::VectorReference< float> &STORAGE)
    {
      // if tv is empty, this indicates that this object is now operating over a new sequence,
      // so it is necessary to regenerate the mapping
      if( m_PositionMapping.IsEmpty())
      {
        AADSSPInfo neighbor_finder( AADSSPInfo::e_MaxHBondNeighborOffset);
        neighbor_finder.SetObject( *this->GetCurrentObject());
        Iterator< biol::AABase> itr( ITR);
        storage::Map< int, size_t> pdb_ids_to_position;
        for( itr.GotoPosition( 0); itr.NotAtEnd(); ++itr)
        {
          pdb_ids_to_position[ itr( 0)->GetPdbID()] = itr.GetPosition();
        }
        m_PositionMapping.AllocateMemory( itr.GetSize());
        for( itr.GotoPosition( 0); itr.NotAtEnd(); ++itr)
        {
          int neighbor_pdb_id( itr( 0)->GetPdbID() + neighbor_finder( itr)( 0));
          storage::Map< int, size_t>::const_iterator itr_neighbor( pdb_ids_to_position.Find( neighbor_pdb_id));
          if( itr_neighbor == pdb_ids_to_position.End())
          {
            m_PositionMapping.PushBack( util::GetUndefined< size_t>());
          }
          else
          {
            m_PositionMapping.PushBack( itr_neighbor->second);
          }
        }
      }

      Iterator< biol::AABase> itr( ITR);
      const size_t current_position( ITR.GetPosition());
      if( m_PositionMapping( current_position) >= itr.GetSize())
      {
        STORAGE = util::GetUndefined< float>();
        return;
      }
      itr.GotoPosition( m_PositionMapping( current_position));

      STORAGE.CopyValues( m_Descriptor->operator ()( itr));
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void AAHbondNeighbor::SetObjectHook()
    {
      m_PositionMapping.Reset();
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< biol::AABase, float> > AAHbondNeighbor::GetInternalDescriptors()
    {
      return iterate::Generic< Base< biol::AABase, float> >( &m_Descriptor, &m_Descriptor + 1);
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AAHbondNeighbor::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Retrieve a descriptor from the highest energy hbonded neighbor of an aa according to the dssp file"
      );
      parameters.AddInitializer
      (
        "",
        "Descriptor to retrieve from the neighbor",
        io::Serialization::GetAgent( &m_Descriptor)
      );

      return parameters;
    }

    //! @brief return the type of symmetry this descriptor has
    //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
    Type::Symmetry AAHbondNeighbor::GetSymmetry() const
    {
      return m_Descriptor->GetType().GetSymmetry();
    }

    //! @brief return whether this descriptor is valid if repeated elements are given
    //! @return true if this descriptor is valid if repeated elements are given
    //! This will be the case if the descriptor may have a legitimate value for A-A
    bool AAHbondNeighbor::ConsiderRepeatedElements() const
    {
      return m_Descriptor->ConsiderRepeatedElements();
    }

    //! @brief get the normal dimension for this descriptor
    //! @return the normal dimension for this descriptor
    size_t AAHbondNeighbor::GetNormalDimension() const
    {
      return m_Descriptor->GetType().GetDimension();
    }

  } // namespace descriptor
} // namespace bcl
