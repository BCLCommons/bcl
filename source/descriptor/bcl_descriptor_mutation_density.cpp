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
#include "descriptor/bcl_descriptor_mutation_density.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model_with_cache.h"
#include "biol/bcl_biol_protein_mutation_set.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_running_average_sd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MutationDensity::s_Instance
    (
      util::Enumerated< Base< biol::Mutation, float> >::AddInstance
      (
        new MutationDensity()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutationDensity::MutationDensity() :
      m_IncludeWeight( false),
      m_IncludeStd( false),
      m_TransitionWidth( 0.5),
      m_CutoffDistance( 5.0),
      m_DefaultValue(),
      m_PseudocountWeight( 1.0)
    {
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    MutationDensity *MutationDensity::Clone() const
    {
      return new MutationDensity( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutationDensity::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &MutationDensity::GetAlias() const
    {
      static const std::string s_name( "MutantPropertyDensity");
      return s_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t MutationDensity::GetNormalSizeOfFeatures() const
    {
      return m_Descriptor->GetSizeOfFeatures() * ( m_IncludeStd ? 2 : 1) + ( m_IncludeWeight ? 1 : 0);
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    iterate::Generic< Base< biol::Mutation, float> > MutationDensity::GetInternalDescriptors()
    {
      return iterate::Generic< Base< biol::Mutation, float> >( &m_Descriptor, &m_Descriptor + 1);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    void MutationDensity::RecalculateImpl
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
      if( m_VoxelGrid.ExtractPositions( *ITR( 0)).IsEmpty())
      {
        STORAGE = util::GetUndefined< float>();
        return;
      }
      auto pairs( m_VoxelGrid.GetNeighbors( ITR( 0)->GetAAs()( 0)->GetCenter(), m_CutoffDistance));
      math::RunningAverageSD< linal::Vector< float> > ave_descriptor;
      size_t pos( 0);
      for( auto itr( pairs.Begin()), itr_end( pairs.End()); itr != itr_end; ++itr)
      {
        if( *itr->First() == *ITR( 0))
        {
          continue;
        }
        const double weight( m_Transition( itr->Second()));
        auto result( ( *m_Descriptor)( m_Mutations.GetValue( itr->First()->GetAAs()( 0))));
        if( result.IsDefined())
        {
          ave_descriptor.AddWeightedObservation( result, weight);
        }
      }
      if( m_IncludeWeight)
      {
        STORAGE( pos++) = ave_descriptor.GetWeight();
      }
      if( !m_DefaultValue.IsEmpty() && m_PseudocountWeight)
      {
        ave_descriptor.AddWeightedObservation( m_DefaultValue, m_PseudocountWeight);
      }
      std::copy( ave_descriptor.GetAverage().Begin(), ave_descriptor.GetAverage().End(), STORAGE.Begin() + pos);
      pos += ave_descriptor.GetAverage().GetSize();
      if( m_IncludeStd)
      {
        std::copy( ave_descriptor.GetStandardDeviation().Begin(), ave_descriptor.GetStandardDeviation().End(), STORAGE.Begin() + pos);
      }
      return;
    }

    //! @brief hook that derived classes can override to add behavior after every time SetObject is called
    void MutationDensity::SetObjectHook()
    {
      m_Mutations.Reset();
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference MutationDensity::GetNormalCachePreference() const
    {
      return e_PreferCache;
    }

    //! @brief function to find the aas of the given type in the specified sequence
    void MutationDensity::FindMutants( const Iterator< biol::Mutation> &ITR)
    {
      util::SiPtr< const biol::ProteinMutationSet> si_pmc( this->GetCurrentObject());

      // Get the filename
      util::SiPtrVector< const biol::Mutation> muts;
      muts.AllocateMemory( si_pmc->GetSize());

      // Remove the last extension
      for( Iterator< biol::Mutation> itr( ITR.Begin()->Begin()); itr.NotAtEnd(); ++itr)
      {
        muts.PushBack( util::ToSiPtr( *itr( 0)));
      }
      m_VoxelGrid.SetObjects( muts);
      for( Iterator< biol::Mutation> itr( ITR.Begin()->Begin()); itr.NotAtEnd(); ++itr)
      {
        for( auto itr_aa( itr( 0)->GetAAs().Begin()), itr_aa_end( itr( 0)->GetAAs().End()); itr_aa != itr_aa_end; ++itr_aa)
        {
          BCL_Assert
          (
            m_Mutations.Insert( std::make_pair( *itr_aa, itr)).second,
            "Need to hard copy AAs...hmm " + itr( 0)->ToString()
            + " current map: " + util::Format()( m_Mutations.GetKeys())
          );
        }
      }
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MutationDensity::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      m_VoxelGrid = assemble::VoxelGridMutation( m_CutoffDistance);
      m_Transition = math::TrigonometricTransition( m_CutoffDistance - m_TransitionWidth, m_CutoffDistance, 1.0, 0.0);
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MutationDensity::GetSerializer() const
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
        "std",
        "whether to include the standard deviation after the main descriptors",
        io::Serialization::GetAgent( &m_IncludeStd),
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

      parameters.AddInitializer
      (
        "default",
        "Default value, applied with pseudocount weight to get the final density average",
        io::Serialization::GetAgent( &m_DefaultValue),
        "(0.0)"
      );
      parameters.AddInitializer
      (
        "pseudocount",
        "Amount of pseudocounts to add into the density calculation",
        io::Serialization::GetAgent( &m_PseudocountWeight),
        "1.0"
      );
      return parameters;
    }

    //! @brief return the type of symmetry this descriptor has
    //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
    Type::Symmetry MutationDensity::GetSymmetry() const
    {
      return m_Descriptor->GetType().GetSymmetry();
    }

    //! @brief return whether this descriptor is valid if repeated elements are given
    //! @return true if this descriptor is valid if repeated elements are given
    //! This will be the case if the descriptor may have a legitimate value for A-A
    bool MutationDensity::ConsiderRepeatedElements() const
    {
      return m_Descriptor->ConsiderRepeatedElements();
    }

    //! @brief get the normal dimension for this descriptor
    //! @return the normal dimension for this descriptor
    size_t MutationDensity::GetNormalDimension() const
    {
      return m_Descriptor->GetType().GetDimension();
    }

  } // namespace descriptor
} // namespace bcl
