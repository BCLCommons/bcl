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
#include "descriptor/bcl_descriptor_molecule_fragment_statistics.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_running_statistics.h"

// external includes - sorted alphabetically

// Uncomment the next line to perform profiling in this class
//#define BCL_PROFILE_MoleculeFragmentStatistics
#ifdef BCL_PROFILE_MoleculeFragmentStatistics
#include "util/bcl_util_stopwatch.h"
#endif

namespace bcl
{
  namespace descriptor
  {
    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MoleculeFragmentStatistics::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new MoleculeFragmentStatistics()
      )
    );

    //! @brief Statistic as string
    //! @param STAT the statistic
    //! @return the string for the stat
    const std::string &MoleculeFragmentStatistics::GetStatisticName( const Statistic &STAT)
    {
      static const std::string s_names[] =
      {
        "Min",
        "Max",
        "Mean",
        "StDev",
        "Sum",
        GetStaticClassName< Statistic>()
      };
      return s_names[ STAT];
    }

    //! @brief default constructor
    MoleculeFragmentStatistics::MoleculeFragmentStatistics()
    {
    }

    //! @brief constructor from implementation and desired statistics
    MoleculeFragmentStatistics::MoleculeFragmentStatistics
    (
      const CheminfoProperty &DESCRIPTOR,
      const storage::Vector< StatisticEnum> &STATS
    ) :
      m_Statistics( STATS),
      m_Descriptor( DESCRIPTOR)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MoleculeFragmentStatistics
    MoleculeFragmentStatistics *MoleculeFragmentStatistics::Clone() const
    {
      return new MoleculeFragmentStatistics( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoleculeFragmentStatistics::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &MoleculeFragmentStatistics::GetAlias() const
    {
      static std::string s_alias( "MolecularFragmentStatistics");
      return s_alias;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t MoleculeFragmentStatistics::GetNormalSizeOfFeatures() const
    {
      return m_Statistics.GetSize() * m_Descriptor->GetSizeOfFeatures();
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
    bool MoleculeFragmentStatistics::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      m_Descriptor->SetDimension( 0);
      if( m_Weighting.IsDefined())
      {
        m_Weighting->SetDimension( 0);
        if( m_Weighting->GetNormalSizeOfFeatures() != size_t( 1))
        {
          ERR_STREAM << "weight must be a molecular property that returns a single, molecular value";
          return false;
        }
      }
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeFragmentStatistics::GetSerializer() const
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
        "weight",
        "Descriptor that can be defined to weight the values of the primary descriptor. "
        "Primarily used for average and standard deviation calculations, but also influences Min and Max, "
        "because if the weight descriptor is defined and is <= 0, the value is ignored for both min/max and ave/std",
        io::Serialization::GetAgent( &m_Weighting)
      );
      parameters.AddInitializer
      (
        "splitter",
        "The method of splitting the molecule into a series of fragments",
        io::Serialization::GetAgent( &m_Splitter)
      );
      parameters.AddInitializer
      (
        "statistics",
        "Statistics desired",
        io::Serialization::GetAgent( &m_Statistics)
      );
      return parameters;
    }

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    void MoleculeFragmentStatistics::Calculate( linal::VectorReference< float> &STORAGE)
    {
      #ifdef BCL_PROFILE_MoleculeFragmentStatistics
      static util::Stopwatch s_split( "Splitting", util::Time( 1, 0), util::Message::e_Standard, true, false);
      s_split.Start();
      #endif
      util::SiPtr< const chemistry::ConformationInterface> mol_ptr
      (
        Base< chemistry::AtomConformationalInterface, float>::GetCurrentObject()
      );
      chemistry::FragmentEnsemble ensemble( m_Splitter->operator()( *mol_ptr));
      #ifdef BCL_PROFILE_MoleculeFragmentStatistics
      s_split.Stop();
      static util::Stopwatch s_descriptor( "Descriptor Calculation", util::Time( 1, 0), util::Message::e_Standard, true, false);
      s_descriptor.Start();
      #endif

      m_Descriptor->SetObject( *Base< chemistry::AtomConformationalInterface, float>::GetCurrentObject());

      if( m_Weighting.IsDefined())
      {
        m_Weighting->SetObject( *Base< chemistry::AtomConformationalInterface, float>::GetCurrentObject());
      }

      math::RunningAverageSD< linal::Vector< float> > ave_sd;
      math::RunningMinMax< linal::Vector< float> > min_max;

      const bool defined_weight( m_Weighting.IsDefined());
      bool have_first_observation( false);
      #ifdef BCL_PROFILE_MoleculeFragmentStatistics
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
        m_Descriptor->SetObject( *itr_mol);
        Iterator< chemistry::AtomConformationalInterface> itr( m_Descriptor->GetType(), *itr_mol);
        BCL_Assert( itr.GetSize() == size_t( 1), "Should have been a molecular descriptor");
        linal::VectorConstReference< float> result( m_Descriptor->operator()( itr));
        if( !result.IsDefined())
        {
          continue;
        }
        float weight( 1.0);
        if( defined_weight)
        {
          m_Weighting->SetObject( *itr_mol);
          weight = m_Weighting->operator()( itr).First();
          #ifdef BCL_PROFILE_MoleculeFragmentStatistics
          BCL_Debug( weight);
          #endif
          if( !util::IsDefined( weight))
          {
            continue;
          }
          else if( weight <= 0.0)
          {
            if( !have_first_observation)
            {
              ave_sd.AddWeightedObservation( result, 0.0);
              min_max = math::RunningMinMax< linal::Vector< float> >( result);
              min_max.Reset();
              have_first_observation = true;
            }
            continue;
          }
          have_first_observation = true;
        }
        ave_sd.AddWeightedObservation( result, weight);
        min_max += result;
      }
      if( ave_sd.GetAverage().GetSize())
      {
        size_t position( 0), sz( ave_sd.GetAverage().GetSize());
        for
        (
          storage::Vector< StatisticEnum>::const_iterator itr( m_Statistics.Begin()), itr_end( m_Statistics.End());
          itr != itr_end;
          ++itr, position += sz
        )
        {
          linal::VectorReference< float> ref( STORAGE.CreateSubVectorReference( sz, position));
          switch( itr->GetEnum())
          {
            case e_Min: ref.CopyValues( min_max.GetMin()); break;
            case e_Max: ref.CopyValues( min_max.GetMax()); break;
            case e_Mean: ref.CopyValues( ave_sd.GetAverage()); break;
            case e_StandardDeviation: ref.CopyValues( ave_sd.GetStandardDeviation()); break;
            case e_Sum: ref.CopyValues( ave_sd.GetAverage()); ref *= float( ave_sd.GetWeight()); break;
            default: BCL_Exit( "Undefined statistic", -1); break;
          }
        }
      }
      else
      {
        STORAGE = util::GetUndefined< float>();
        size_t position( 0), sz( STORAGE.GetSize() / m_Statistics.GetSize());
        for
        (
          storage::Vector< StatisticEnum>::const_iterator itr( m_Statistics.Begin()), itr_end( m_Statistics.End());
          itr != itr_end;
          ++itr, position += sz
        )
        {
          if( itr->GetEnum() == e_Sum)
          {
            linal::VectorReference< float> ref( STORAGE.CreateSubVectorReference( sz, position));
            ref = float( 0.0);
          }
        }
      }
      #ifdef BCL_PROFILE_MoleculeFragmentStatistics
      s_descriptor.Stop();
      #endif
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference MoleculeFragmentStatistics::GetNormalCachePreference() const
    {
      return e_PreferCache;
    }

  } // namespace descriptor
} // namespace bcl
