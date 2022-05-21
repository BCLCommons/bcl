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
#include "descriptor/bcl_descriptor_mutation_aa_property.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_protein_mutation_set.h"
#include "chemistry/bcl_chemistry_fragment_configuration_shared.h"
#include "io/bcl_io_serialization.h"
#include "iterate/bcl_iterate_generic.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> MutationAAProperty::s_Instance
    (
      util::Enumerated< Base< biol::Mutation, float> >::AddInstance( new MutationAAProperty( false))
    );
    const util::SiPtr< const util::ObjectInterface> MutationAAProperty::s_NativeInstance
    (
      util::Enumerated< Base< biol::Mutation, float> >::AddInstance( new MutationAAProperty( true))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from exposure measure
    //! @param CONSIDER_SIDE_CHAIN whether to consider side chain atoms in the calculation
    //! @param CONSIDER_BACK_BONE whether to consider back bone atoms in the calculation
    //! @param PROPERTY property to calculate
    MutationAAProperty::MutationAAProperty
    (
      const bool &NATIVE,
      const util::Implementation< Base< biol::AABase, float> > &PROPERTY
    ) :
      m_Property( PROPERTY),
      m_UseNative( NATIVE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MutationAAProperty
    MutationAAProperty *MutationAAProperty::Clone() const
    {
      return new MutationAAProperty( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutationAAProperty::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &MutationAAProperty::GetAlias() const
    {
      static const std::string s_mut_name( "MutantAADescriptor"),
                               s_nat_name( "NativeAADescriptor");
      return m_UseNative ? s_nat_name : s_mut_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t MutationAAProperty::GetNormalSizeOfFeatures() const
    {
      return m_Property->GetNormalSizeOfFeatures();
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MutationAAProperty::Calculate
    (
      const iterate::Generic< const biol::Mutation> &ELEMENT,
      linal::VectorReference< float> &STORAGE
    )
    {
      util::SiPtr< const assemble::ProteinModelWithCache> pmwc
      (
        m_UseNative ? m_MutationsMap->GetValue( biol::Mutation()) : m_MutationsMap->GetValue( *ELEMENT)
      );
      if( !m_UseNative)
      {
        m_Property->SetObject( *pmwc);
      }

      // this is a very slow, but reliable way of finding the right AA for the iterator. If larger (10-100k) mutation datasets
      // become available, we should definitely move to a predictive search
      Iterator< biol::AABase> itr_aa( pmwc->GetIterator());
      for( ; itr_aa.NotAtEnd(); ++itr_aa)
      {
        if( itr_aa( 0)->GetSeqID() == ELEMENT->GetResidueNumber())
        {
          break;
        }
      }
      if( !itr_aa.NotAtEnd())
      {
        STORAGE = util::GetUndefined< float>();
      }
      else
      {
        STORAGE.CopyValues( m_Property->operator ()( itr_aa));
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutationAAProperty::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Uses AA descriptors on individual residues of the "
        + std::string( m_UseNative ? "native/wild-type" : "mutant") + " protein"
      );
      serializer.AddInitializer
      (
        "",
        "aa property desired",
        io::Serialization::GetAgent( &m_Property)
      );

      return serializer;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void MutationAAProperty::SetObjectHook()
    {
      util::SiPtr< const biol::ProteinMutationSet> mut_set( this->GetCurrentObject());
      m_ModelWMutations = mut_set;
      static sched::Mutex s_mutex;
      s_mutex.Lock();
      static
        storage::Map
        <
          util::SiPtr< const biol::ProteinMutationSet>,
          storage::Map< biol::Mutation, util::ShPtr< assemble::ProteinModelWithMutations> >
        > s_map;
      storage::Map< biol::Mutation, util::ShPtr< assemble::ProteinModelWithMutations> > &map
      (
        s_map[ m_ModelWMutations]
      );
      if( map.IsEmpty())
      {
        const bool req_coordinates( mut_set->GetNativeType().GetRequiresCoordinates());
        for( auto itr( mut_set->GetIterator()); itr.NotAtEnd(); ++itr)
        {
          map.Insert
          (
            std::make_pair
            (
              *itr,
              util::ShPtr< assemble::ProteinModelWithMutations>
              (
                new assemble::ProteinModelWithMutations
                (
                  assemble::ProteinModel::HardCopy( mut_set->GetNativeType()),
                  req_coordinates
                )
              )
            )
          ).first->second->Mutate( *itr);
        }
        // add native type
        map.Insert
        (
          std::make_pair
          (
            biol::Mutation(),
            util::ShPtr< assemble::ProteinModelWithMutations>
            (
              new assemble::ProteinModelWithMutations
              (
                assemble::ProteinModel::HardCopy( mut_set->GetNativeType()),
                req_coordinates
              )
            )
          )
        );
      }
      m_MutationsMap = util::ToSiPtr( map);
      s_mutex.Unlock();
      if( m_UseNative)
      {
        m_Property->SetObject( *m_MutationsMap->GetValue( biol::Mutation()));
      }
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference MutationAAProperty::GetNormalCachePreference() const
    {
      // never cache this descriptor for mutants, since their values will change every switch of the mutant
      return e_PreferCache;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MutationAAProperty::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      m_Property->SetDimension( 1);
      return true;
    }

  } // namespace descriptor
} // namespace bcl
