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
#include "chemistry/bcl_chemistry_configuration_set.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief construct from list of ConstitutionInterface
    //! @param MOLECULES an iterator to configuration of molecules
    ConfigurationSet::ConfigurationSet( iterate::Generic< const ConfigurationInterface> MOLECULES)
    {
      for( ; MOLECULES.NotAtEnd(); ++MOLECULES)
      {
        Insert( *MOLECULES);
      }
    }

    //! @brief Clone function
    //! @return pointer to new ConfigurationSet
    ConfigurationSet *ConfigurationSet::Clone() const
    {
      return new ConfigurationSet( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ConfigurationSet::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ConfigurationSet::Read( std::istream &ISTREAM)
    {
      //read
//      *this = ConfigurationSet( ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ConfigurationSet::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
//      SmallMoleculeFactory::WriteToMDLFile( *this, OSTREAM);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief find a particular configuration
    //! @param FRAGMENT a configuration to search for
    //! @return an iterator to the configuration in this set
    ConfigurationSet::const_iterator ConfigurationSet::Find
    (
      const ConfigurationInterface &FRAGMENT,
      util::SiPtr< storage::Vector< size_t> > ISOMORPHISM
    ) const
    {
      storage::Vector< size_t> isomorphism;
      if( !ISOMORPHISM.IsDefined())
      {
        // create a dummy isomorphism object that can usually be used to avoid the need for 2 isomorphism searches
        ISOMORPHISM = util::ToSiPtrNonConst( isomorphism);
      }

      // find the constitution
      ConstitutionSet::const_iterator itr( m_Constitutions.Find( FragmentConstitutionShared( FRAGMENT), ISOMORPHISM));
      if( itr == m_Constitutions.End())
      {
        return End();
      }

      // find the constitution in the map
      const ConfigurationSetSameConstitution &configurations_for_constitution( m_ConfigurationsMap.Find( *itr)->second);

      // return the result of calling find on the configuration set same constitution
      ConfigurationSet::const_iterator configurations_itr
      (
        configurations_for_constitution.Find( FRAGMENT, *ISOMORPHISM)
      );
      if( configurations_itr == configurations_for_constitution.End())
      {
        return End();
      }
      return m_IteratorMap.find( *configurations_itr)->second;
    }

    //! @brief append a molecule to SetConstitution
    //! @param FRAGMENT fragment constitution shared that needs to be added to set constitution
    //! @param ISOMORPHISM si ptr to vector that can be used to cache the isomorphism from another layer
    //! @return a pair containing iterator to configuration layer of FRAGMENT and whether configuration has been seen
    //!         earlier (false) or seen for the first time (true)
    std::pair< ConfigurationSet::const_iterator, bool> ConfigurationSet::Insert
    (
      const ConfigurationInterface &FRAGMENT,
      util::SiPtr< storage::Vector< size_t> > ISOMORPHISM
    )
    {
      storage::Vector< size_t> isomorphism;
      if( !ISOMORPHISM.IsDefined())
      {
        // create a dummy isomorphism object that can usually be used to avoid the need for 2 isomorphism searches
        ISOMORPHISM = util::ToSiPtrNonConst( isomorphism);
      }
      // insert into the constitution
      std::pair< ConstitutionSet::const_iterator, bool>
        constitution_insert( m_Constitutions.Insert( FragmentConstitutionShared( FRAGMENT), *ISOMORPHISM));

      // create a pair to hold the result of the insert on the appropriate configuration set same constitution
      std::pair< ConfigurationSet::const_iterator, bool> configuration_insert;

      // if a new constitution was found, create a new configuration set object for it and add it to the map
      if( constitution_insert.second)
      {
        configuration_insert =
          m_ConfigurationsMap.Insert
          (
            std::make_pair
            (
              util::ToSiPtr( **constitution_insert.first),
              ConfigurationSetSameConstitution( *constitution_insert.first)
            )
          ).first->second.Insert( FRAGMENT, *constitution_insert.first, *ISOMORPHISM);
      }
      else
      {
        // same constitution as something seen previously
        configuration_insert =
          m_ConfigurationsMap[ *constitution_insert.first].Insert( FRAGMENT, *constitution_insert.first, *ISOMORPHISM);
      }
      if( configuration_insert.second)
      {
        m_Configurations.PushBack( *configuration_insert.first);
        m_IteratorMap[ *configuration_insert.first] = m_Configurations.Last();
      }
      return configuration_insert;
    }

  } // namespace chemistry
} // namespace bcl
