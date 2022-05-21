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

#ifndef BCL_CHEMISTRY_CONFIGURATION_SET_H_
#define BCL_CHEMISTRY_CONFIGURATION_SET_H_

// include the namespace header
#include "bcl_chemistry.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_configuration_set_same_constitution.h"
#include "bcl_chemistry_constitution_set.h"
#include "bcl_chemistry_fragment_constitution_shared.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_hydrogens_handler.h"
#include "math/bcl_math_range.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ConfigurationSet
    //! @brief Container class for FragmentConstitutionShared objects
    //!
    //! @see @link example_chemistry_configuration_set.cpp @endlink
    //! @author kothiwsk, mendenjl
    //! @date Mar 09, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ConfigurationSet :
      public util::ObjectInterface
    {
    public:

    //////////////
    // typedefs //
    //////////////

      typedef util::ShPtrList< FragmentConfigurationShared>::const_iterator const_iterator;

    private:

    //////////
    // data //
    //////////

      ConstitutionSet m_Constitutions;

      //! map of constitution to various configurations
      storage::Map< util::SiPtr< const FragmentConstitutionShared>, ConfigurationSetSameConstitution> m_ConfigurationsMap;

      //! map from configurations of same constitution iterators to corresponding iterator on m_Configurations
      std::map< util::SiPtr< const FragmentConfigurationShared>, const_iterator> m_IteratorMap;

      //! All configurations that are known
      util::ShPtrList< FragmentConfigurationShared> m_Configurations;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ConfigurationSet()
      {}

      //! @brief construct from list of FragmentConstitutionShared
      //! @param MOLECULES an iterator to configuration of molecules
      ConfigurationSet( iterate::Generic< const ConfigurationInterface> MOLECULES);

      //! @brief Clone function
      //! @return pointer to new ConfigurationSet
      ConfigurationSet *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the constitution set
      //! @return the constitution set
      const ConstitutionSet &GetConstitutions() const
      {
        return m_Constitutions;
      }
      //! @brief return a map of constitution and related configurations
      //! @return a map of constitution and related configurations
      const storage::Map< util::SiPtr< const FragmentConstitutionShared>, ConfigurationSetSameConstitution> &GetConfigurationsMap() const
      {
        return m_ConfigurationsMap;
      }

      //! @brief return a list containing all configurations
      //! @return a list containing all configurations
      const util::ShPtrList< FragmentConfigurationShared> &GetConfigurations() const
      {
        return m_Configurations;
      }

      //! @brief returns the total number of configurations
      //! @return the total number of configurations
      size_t GetSize() const
      {
        return m_Configurations.GetSize();
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns const iterator begin to configurations
      //! @return const iterator begin to configurations
      const_iterator Begin() const
      {
        return m_Configurations.Begin();
      }

      //! @brief returns const iterator end to configurations
      //! @return const iterator end to configurations
      const_iterator End() const
      {
        return m_Configurations.End();
      }

      //! @brief find a particular configuration
      //! @param FRAGMENT a configuration to search for
      //! @param ISOMORPHISM si ptr to vector that can be used to cache the isomorphism from another layer
      //! @return a simple pointer to the configuration in the set
      const_iterator Find
      (
        const ConfigurationInterface &FRAGMENT,
        util::SiPtr< storage::Vector< size_t> > ISOMORPHISM = util::SiPtr< storage::Vector< size_t> >()
      ) const;

      //! @brief append a molecule to SetConstitution
      //! @param FRAGMENT fragment constitution shared that needs to be added to set constitution
      //! @param ISOMORPHISM si ptr to vector that can be used to cache the isomorphism from another layer
      //! @return a pair containing iterator to configuration layer of FRAGMENT and whether configuration has been seen
      //!         earlier (false) or seen for the first time (true)
      std::pair< const_iterator, bool> Insert
      (
        const ConfigurationInterface &FRAGMENT,
        util::SiPtr< storage::Vector< size_t> > ISOMORPHISM = util::SiPtr< storage::Vector< size_t> >()
      );

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    }; // class ConfigurationSet

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_CONFIGURATION_SET_H_
