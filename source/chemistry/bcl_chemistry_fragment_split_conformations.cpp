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
#include "chemistry/bcl_chemistry_fragment_split_conformations.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_sample_conformations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentSplitConformations::s_Instance
    (
      util::Enumerated< FragmentSplitInterface>::AddInstance( new FragmentSplitConformations())
    );

    // cache map
    storage::Map< std::string, storage::Pair< storage::List< FragmentSplitConformations::t_SplitCacheTriplet>, sched::Mutex> >
      FragmentSplitConformations::s_ConformationCache =
        storage::Map< std::string, storage::Pair< storage::List< FragmentSplitConformations::t_SplitCacheTriplet>, sched::Mutex> >();

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    FragmentSplitConformations::FragmentSplitConformations() :
      m_UseCache( true)
    {
    }

    //! virtual copy constructor
    FragmentSplitConformations *FragmentSplitConformations::Clone() const
    {
      return new FragmentSplitConformations( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &FragmentSplitConformations::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentSplitConformations::GetAlias() const
    {
      static const std::string s_name_isolate( "SampleConformations");
      return s_name_isolate;
    }

    //! get the minimum size of a component of interest
    const size_t FragmentSplitConformations::GetMinSize() const
    {
      return 0;
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief operator returns an ensemble of conformations
    //! @param MOLECULE molecule of interest
    //! @return enseble of conformations
    FragmentEnsemble FragmentSplitConformations::operator()( const ConformationInterface &MOLECULE) const
    {
      if( m_UseCache)
      {
        // because this class is often called multiple times for the same molecule, the results from each molecule are
        // temporarily cached locally
        storage::Vector< AtomType> atom_types( MOLECULE.GetAtomTypesVector());
        storage::Vector< sdf::BondInfo> bond_info( MOLECULE.GetBondInfo());

        m_Mutex->Lock();
        size_t queue_size( 0);

        for
        (
          storage::List< t_SplitCacheTriplet>::iterator
            itr( m_LocalCache->Begin()), itr_end( m_LocalCache->End());
          itr != itr_end;
          ++itr, ++queue_size
        )
        {
          if( itr->First() == atom_types && itr->Second() == bond_info)
          {
            FragmentEnsemble &conf_ensemble( itr->Third());
            if( itr != m_LocalCache->Begin())
            {
              m_LocalCache->InternalData().splice( m_LocalCache->Begin(), m_LocalCache->InternalData(), itr);
            }
            m_Mutex->Unlock();
            return conf_ensemble;
          }
        }
        m_Mutex->Unlock();

        // sample the conformations
        FragmentEnsemble ensemble( m_SampleConformations( MOLECULE).First());

        // store the conformations in a queue
        m_Mutex->Lock();
        if( queue_size >= size_t( 32))
        {
          m_LocalCache->PopBack();
        }
        m_LocalCache->PushFront( t_SplitCacheTriplet( atom_types, bond_info, ensemble));
        m_Mutex->Unlock();
        return ensemble;
      }

      return m_SampleConformations( MOLECULE).First();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer FragmentSplitConformations::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.Merge( m_SampleConformations.GetCompleteSerializer());
      serializer.AddInitializer
      (
        "cache",
        "Whether to cache conformations. If 2+ descriptors or classes are creating conformers for the same molecule, "
        "then caching makes a lot of sense, but it can substantially increase memory overhead, especially if they are "
        "not needed",
        io::Serialization::GetAgent( &m_UseCache),
        "True"
      );
      return serializer;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentSplitConformations::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      // locate this variable in the cache
      static sched::Mutex s_mutex;
      if( !m_SampleConformations.ReadInitializerSuccessHook( LABEL, ERROR_STREAM))
      {
        return false;
      }
      if( m_UseCache)
      {
        s_mutex.Lock();
        storage::Pair< storage::List< t_SplitCacheTriplet>, sched::Mutex> &reference
        (
          s_ConformationCache[ m_SampleConformations.GetLabel().ToString()]
        );
        m_LocalCache = util::ToSiPtrNonConst( reference.First());
        m_Mutex = util::ToSiPtrNonConst( reference.Second());

        s_mutex.Unlock();
      }
      return true;
    }

  } // namespace chemistry
} // namespace bcl
