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
#include "descriptor/bcl_descriptor_mutation_aa_density.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "biol/bcl_biol_protein_mutation_set.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MutationAADensity::s_Instance
    (
      util::Enumerated< Base< biol::Mutation, float> >::AddInstance
      (
        new MutationAADensity()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutationAADensity::MutationAADensity() :
      m_IncludeWeight( false),
      m_TransitionWidth( 0.5),
      m_CutoffDistance( 5.0)
    {
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    MutationAADensity *MutationAADensity::Clone() const
    {
      return new MutationAADensity( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutationAADensity::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &MutationAADensity::GetAlias() const
    {
      static const std::string s_name( "MutantAAPropertyDensity");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t MutationAADensity::GetNormalSizeOfFeatures() const
    {
      return m_Descriptor->GetSizeOfFeatures() + ( m_IncludeWeight ? 1 : 0);
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< biol::Mutation, float> > MutationAADensity::GetInternalDescriptors()
    {
      return iterate::Generic< Base< biol::Mutation, float> >();
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    void MutationAADensity::RecalculateImpl
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
      auto pairs( m_VoxelGrid.GetNeighbors( ITR( 0)->GetAAs()( 0)->GetCenter(), m_CutoffDistance));
      math::RunningAverage< linal::Vector< float> > ave_descriptor;
      const size_t n_features( m_Descriptor->GetSizeOfFeatures());
      size_t pos( 0);
      for( auto itr( pairs.Begin()), itr_end( pairs.End()); itr != itr_end; ++itr)
      {
        BCL_Debug( itr->First()->GetIdentification());
        const double weight( m_Transition( itr->Second()));
        ave_descriptor.AddWeightedObservation
        (
          ( *m_Descriptor)( m_Mutations.GetValue( std::make_pair( itr->First()->GetSeqID(), itr->First()->GetChainID()))),
          weight
        );
      }
      if( m_IncludeWeight)
      {
        STORAGE( pos++) = ave_descriptor.GetWeight();
      }
      std::copy( ave_descriptor.GetAverage().Begin(), ave_descriptor.GetAverage().End(), STORAGE.Begin() + 1);

      return;
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference MutationAADensity::GetNormalCachePreference() const
    {
      return e_PreferCache;
    }

    //! @brief function to find the aas of the given type in the specified sequence
    void MutationAADensity::FindMutants( const Iterator< biol::Mutation> &ITR)
    {
      util::SiPtr< const biol::ProteinMutationSet> si_pmc( this->GetCurrentObject());

      // Get the filename
      auto itr_aas( si_pmc->GetNativeType().GetIterator());
      util::SiPtrVector< const biol::AABase> aas;
      aas.AllocateMemory( itr_aas.GetSize());

      // Remove the last extension
      Iterator< biol::AABase> desc_itr( itr_aas);
      for( ; itr_aas.NotAtEnd(); ++itr_aas, ++desc_itr)
      {
        aas.PushBack( *itr_aas);
        m_Mutations.Insert
        (
          std::make_pair
          (
            storage::Pair< int, char>( itr_aas->GetSeqID(), itr_aas->GetChainID()),
            desc_itr
          )
        );
      }
      m_VoxelGrid.SetObjects( aas);
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MutationAADensity::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      m_VoxelGrid = assemble::VoxelGridAA( m_CutoffDistance);
      m_Transition = math::TrigonometricTransition( m_CutoffDistance - m_TransitionWidth, 1.0, m_CutoffDistance, 0.0);
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutationAADensity::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "retrieve properties of the nearest mutations to the current one");

      parameters.AddInitializer
      (
        "weight",
        "whether to include the weight (number of mutations that went into the computation)",
        io::Serialization::GetAgent( &m_IncludeWeight),
        "False"
      );

      parameters.AddInitializer
      (
        "transition",
        "width, in angstroms, of the transition between mutations receiving full weight and receiving 0 weight",
        io::Serialization::GetAgentWithMin( &m_TransitionWidth, float( 0.0)),
        "0.5"
      );

      parameters.AddInitializer
      (
        "cutoff",
        "maximum distance to consider",
        io::Serialization::GetAgentWithMin( &m_CutoffDistance, float( 0.0))
      );

      parameters.AddInitializer
      (
        "descriptor",
        "Descriptor to retrieve for the nearest mutations",
        io::Serialization::GetAgent( &m_Descriptor)
      );

      return parameters;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void MutationAADensity::SetObjectHook()
    {
      m_Mutations.Reset();
      util::SiPtr< const biol::ProteinMutationSet> si_pmc( this->GetCurrentObject());

      // Get the filename
      m_Descriptor->SetObject( si_pmc->GetNativeType());
    }

    //! @brief return the type of symmetry this descriptor has
    //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
    Type::Symmetry MutationAADensity::GetSymmetry() const
    {
      return m_Descriptor->GetType().GetSymmetry();
    }

    //! @brief return whether this descriptor is valid if repeated elements are given
    //! @return true if this descriptor is valid if repeated elements are given
    //! This will be the case if the descriptor may have a legitimate value for A-A
    bool MutationAADensity::ConsiderRepeatedElements() const
    {
      return m_Descriptor->ConsiderRepeatedElements();
    }

    //! @brief get the normal dimension for this descriptor
    //! @return the normal dimension for this descriptor
    size_t MutationAADensity::GetNormalDimension() const
    {
      return m_Descriptor->GetType().GetDimension();
    }

  } // namespace descriptor
} // namespace bcl
