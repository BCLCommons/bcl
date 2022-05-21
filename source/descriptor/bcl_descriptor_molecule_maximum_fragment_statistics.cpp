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
#include "descriptor/bcl_descriptor_molecule_maximum_fragment_statistics.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

// Uncomment the next line to perform profiling in this class
//#define BCL_PROFILE_MoleculeMaximumFragmentStatistics

namespace bcl
{
  namespace descriptor
  {
    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MoleculeMaximumFragmentStatistics::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new MoleculeMaximumFragmentStatistics()
      )
    );

    //! @brief default constructor
    MoleculeMaximumFragmentStatistics::MoleculeMaximumFragmentStatistics()
    {
    }

    //! @brief constructor from implementation and desired statistics
    MoleculeMaximumFragmentStatistics::MoleculeMaximumFragmentStatistics
    (
      const CheminfoProperty &DESCRIPTOR
    ) :
      m_Descriptor( DESCRIPTOR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MoleculeMaximumFragmentStatistics
    MoleculeMaximumFragmentStatistics *MoleculeMaximumFragmentStatistics::Clone() const
    {
      return new MoleculeMaximumFragmentStatistics( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoleculeMaximumFragmentStatistics::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &MoleculeMaximumFragmentStatistics::GetAlias() const
    {
      static std::string s_alias( "MolecularMaxFragmentStatistics");
      return s_alias;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t MoleculeMaximumFragmentStatistics::GetNormalSizeOfFeatures() const
    {
      return m_Descriptor->GetSizeOfFeatures();
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

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool MoleculeMaximumFragmentStatistics::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      m_Descriptor->SetDimension( 0);
      if( m_Maximum.IsDefined())
      {
        m_Maximum->SetDimension( 0);
        if( m_Maximum->GetNormalSizeOfFeatures() != size_t( 1))
        {
          ERR_STREAM << "weight must be a molecular property that returns a single, molecular value";
          return false;
        }
      }
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeMaximumFragmentStatistics::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Statistics of a descriptor across the series of fragments split off the original molecule"
      );
      parameters.AddInitializer
      (
        "",
        "The descriptor to compute the statistics of across the series of fragments",
        io::Serialization::GetAgent( &m_Descriptor)
      );
      parameters.AddOptionalInitializer
      (
        "maximize",
        "Descriptor that can be defined to weight the values of the primary descriptor. "
        "Primarily used for average and standard deviation calculations, but also influences Min and Max, "
        "because if the weight descriptor is defined and is <= 0, the value is ignored for both min/max and ave/std",
        io::Serialization::GetAgent( &m_Maximum)
      );
      parameters.AddInitializer
      (
        "splitter",
        "The method of splitting the molecule into a series of fragments",
        io::Serialization::GetAgent( &m_Splitter)
      );

      return parameters;
    }

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    void MoleculeMaximumFragmentStatistics::Calculate( linal::VectorReference< float> &STORAGE)
    {
      #ifdef BCL_PROFILE_MoleculeMaximumFragmentStatistics
      static util::Stopwatch s_split( "Splitting", util::Time( 1, 0), util::Message::e_Standard, true, false);
      s_split.Start();
      #endif
      util::SiPtr< const chemistry::ConformationInterface> mol_ptr
      (
        Base< chemistry::AtomConformationalInterface, float>::GetCurrentObject()
      );
      chemistry::FragmentEnsemble ensemble( m_Splitter->operator()( *mol_ptr));
      #ifdef BCL_PROFILE_MoleculeMaximumFragmentStatistics
      s_split.Stop();
      static util::Stopwatch s_descriptor( "Descriptor Calculation", util::Time( 1, 0), util::Message::e_Standard, true, false);
      s_descriptor.Start();
      #endif

      m_Descriptor->SetObject( *Base< chemistry::AtomConformationalInterface, float>::GetCurrentObject());

      util::SiPtr< const chemistry::FragmentComplete> max_conformation;
      double max_prediction( -std::numeric_limits< double>::max());
      #ifdef BCL_PROFILE_MoleculeMaximumFragmentStatistics
      BCL_Debug( ensemble.GetMolecules().GetSize());
      #endif
      for
      (
        chemistry::FragmentEnsemble::const_iterator
          itr_mol( ensemble.Begin()), itr_mol_end( ensemble.End());
        itr_mol != itr_mol_end;
        ++itr_mol
      )
      {
        m_Maximum->SetObject( *Base< chemistry::AtomConformationalInterface, float>::GetCurrentObject());
        float prediction( 1.0);

        m_Maximum->SetObject( *itr_mol);
        Iterator< chemistry::AtomConformationalInterface> itr( m_Descriptor->GetType(), *itr_mol);
        prediction = m_Maximum->operator()( itr).First();

        if( prediction > max_prediction)
        {
          max_prediction = prediction;
          max_conformation = util::SiPtr< const chemistry::FragmentComplete>( *itr_mol);
        }
      }

      if( max_conformation.IsDefined())
      {
        m_Descriptor->SetObject( *max_conformation);
        Iterator< chemistry::AtomConformationalInterface> itr( m_Descriptor->GetType(), *max_conformation);
        BCL_Assert( itr.GetSize() == size_t( 1), "Should have been a molecular descriptor");
        linal::VectorConstReference< float> result( m_Descriptor->operator()( itr));
        STORAGE.CopyValues( result);
      }
      else
      {
        STORAGE = util::GetUndefined< float>();
      }
      #ifdef BCL_PROFILE_MoleculeMaximumFragmentStatistics
      s_descriptor.Stop();
      #endif
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference MoleculeMaximumFragmentStatistics::GetNormalCachePreference() const
    {
      return e_PreferCache;
    }

  } // namespace descriptor
} // namespace bcl
