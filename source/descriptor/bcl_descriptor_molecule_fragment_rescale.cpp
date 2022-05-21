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
#include "descriptor/bcl_descriptor_molecule_fragment_rescale.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_running_statistics.h"
#include "sdf/bcl_sdf_mdl_handler.h"
// external includes - sorted alphabetically

// Uncomment the next line to perform profiling in this class
//#define BCL_PROFILE_MoleculeFragmentRescale
#ifdef BCL_PROFILE_MoleculeFragmentRescale
#include "util/bcl_util_stopwatch.h"
#endif

namespace bcl
{
  namespace descriptor
  {
    // add each of the possible instances to the enumerated instances
    const util::SiPtr< const util::ObjectInterface> MoleculeFragmentRescale::s_Instance
    (
      util::Enumerated< Base< chemistry::AtomConformationalInterface, float> >::AddInstance
      (
        new MoleculeFragmentRescale()
      )
    );

    //! @brief Statistic as string
    //! @param STAT the statistic
    //! @return the string for the stat
    const std::string &MoleculeFragmentRescale::GetStatisticName( const Statistic &STAT)
    {
      static const std::string s_names[] =
      {
        "MinMax",
        "ZScore",
        "SubtractMean",
        "DivideMean",
        GetStaticClassName< Statistic>()
      };
      return s_names[ STAT];
    }

    //! @brief default constructor
    MoleculeFragmentRescale::MoleculeFragmentRescale() : m_CacheConformations( true)
    {
    }

    //! @brief constructor from implementation and desired statistics
    MoleculeFragmentRescale::MoleculeFragmentRescale
    (
      const CheminfoProperty &DESCRIPTOR,
      const storage::Vector< StatisticEnum> &STATS
    ) :
      m_Statistics( STATS),
      m_Descriptor( DESCRIPTOR),
      m_CacheConformations( true)
    {
    }

    //! @brief Clone function
    //! @return pointer to new MoleculeFragmentRescale
    MoleculeFragmentRescale *MoleculeFragmentRescale::Clone() const
    {
      return new MoleculeFragmentRescale( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MoleculeFragmentRescale::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &MoleculeFragmentRescale::GetAlias() const
    {
      static std::string s_alias( "MolecularFragmentRescale");
      return s_alias;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t MoleculeFragmentRescale::GetNormalSizeOfFeatures() const
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
    bool MoleculeFragmentRescale::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      m_Descriptor->SetDimension( 0);
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeFragmentRescale::GetSerializer() const
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
      parameters.AddInitializer
      (
        "splitter",
        "The method of splitting the molecule into a series of fragments",
        io::Serialization::GetAgent( &m_Splitter)
      );
      parameters.AddInitializer
      (
        "rescale",
        "Rescalings to perform",
        io::Serialization::GetAgent( &m_Statistics)
      );
      parameters.AddInitializer
      (
        "cache",
        "Whether to cache conformation ensembles for each conformation",
        io::Serialization::GetAgent( &m_CacheConformations),
        "1"
      );
      return parameters;
    }

    //! @brief calculate the descriptors
    //! @param STORAGE storage for the descriptor
    void MoleculeFragmentRescale::Calculate( linal::VectorReference< float> &STORAGE)
    {
      #ifdef BCL_PROFILE_MoleculeFragmentRescale
      static util::Stopwatch s_split( "Splitting", util::Time( 1, 0), util::Message::e_Standard, true, false);
      s_split.Start();
      #endif
      util::SiPtr< const chemistry::ConformationInterface> mol_ptr
      (
        Base< chemistry::AtomConformationalInterface, float>::GetCurrentObject()
      );
      static storage::Map< std::string, chemistry::FragmentEnsemble> s_ensembles; // cache ensembles internally
      static sched::Mutex s_mutex;
      chemistry::FragmentEnsemble ensemble_store;
      std::string full_hash;
      if( m_CacheConformations)
      {
        std::string hash
        (
          sdf::MdlHandler::CreateConfigurationalHashString( mol_ptr->GetAtomInfo(), mol_ptr->GetBondInfo(), chemistry::ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness)
        );
        full_hash = hash + m_Splitter->GetString();
        s_mutex.Lock();
        if( s_ensembles[ full_hash].IsEmpty())
        {
          s_mutex.Unlock();
          ensemble_store = m_Splitter->operator()( *mol_ptr);
          s_mutex.Lock();
          if( s_ensembles[ full_hash].IsEmpty())
          {
            s_ensembles[ full_hash] = ensemble_store;
          }
        }
      }
      else
      {
        ensemble_store = m_Splitter->operator()( *mol_ptr);
      }
      chemistry::FragmentEnsemble &ensemble
      (
        m_CacheConformations
        ? s_ensembles[ full_hash]
        : ensemble_store
      );
      if( m_CacheConformations)
      {
        s_mutex.Unlock();
      }

      #ifdef BCL_PROFILE_MoleculeFragmentRescale
      s_split.Stop();
      static util::Stopwatch s_descriptor( "Descriptor Calculation", util::Time( 1, 0), util::Message::e_Standard, true, false);
      s_descriptor.Start();
      #endif

      m_Descriptor->SetObject( *Base< chemistry::AtomConformationalInterface, float>::GetCurrentObject());
      Iterator< chemistry::AtomConformationalInterface> itr_descriptor( m_Descriptor->GetType(), *mol_ptr);
      linal::Vector< float> this_result( m_Descriptor->operator()( itr_descriptor));

      math::RunningAverageSD< linal::Vector< float> > ave_sd;
      math::RunningMinMax< linal::Vector< float> > min_max;

      ave_sd += this_result;
      min_max += this_result;

      bool have_first_observation( false);
      #ifdef BCL_PROFILE_MoleculeFragmentRescale
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
        ave_sd += result;
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
          if( *itr == e_MinMax)
          {
            for( size_t i( 0); i < sz; ++i)
            {
              ref( i) = min_max.GetRange()( i) > 1e-8 ? 2.0 * ( this_result( i) - min_max.GetMin()( i)) / min_max.GetRange()( i) - 1.0 : 0.0;
            }
          }
          else if( *itr == e_ZScore)
          {
            for( size_t i( 0); i < sz; ++i)
            {
              ref( i) = ave_sd.GetStandardDeviation()( i) > 1e-8 ? ( this_result( i) - ave_sd.GetAverage()( i)) / ave_sd.GetStandardDeviation()( i) : 0.0;
            }
          }
          else if( *itr == e_SubtractMean)
          {
            for( size_t i( 0); i < sz; ++i)
            {
              ref( i) = ( this_result( i) - ave_sd.GetAverage()( i));
            }
          }
          else // if( *itr == e_DivideMean)
          {
            for( size_t i( 0); i < sz; ++i)
            {
              ref( i) = math::Absolute( ave_sd.GetAverage()( i)) > 1e-8 ? this_result( i) / ave_sd.GetAverage()( i) : 1.0;
            }
          }
        }
      }
      else
      {
        STORAGE = util::GetUndefined< float>();
      }
      #ifdef BCL_PROFILE_MoleculeFragmentRescale
      s_descriptor.Stop();
      #endif
    }

    //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
    //! @return the cache preference, assuming this feature has its normal dimension setting
    CachePreference MoleculeFragmentRescale::GetNormalCachePreference() const
    {
      return e_PreferCache;
    }

  } // namespace descriptor
} // namespace bcl
