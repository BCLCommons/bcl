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
#include "chemistry/bcl_chemistry_conformation_set_same_configuration.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_set_same_constitution.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief default constructor
    ConformationSetSameConstitution::ConformationSetSameConstitution()
    {
    }

    //! @brief constructor
    //! @param METHOD conformation comparison method to determine if conformations are same
    //! @param TOLERANCE tolerance value for comparing conformation values
    ConformationSetSameConstitution::ConformationSetSameConstitution
    (
      util::Implementation< ConformationComparisonInterface> METHOD,
      double TOLERANCE
    )
    :
      m_ConfomationComparer( METHOD),
      m_Tolerance( TOLERANCE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new ConformationSetSameConstitution
    ConformationSetSameConstitution *ConformationSetSameConstitution::Clone() const
    {
      return new ConformationSetSameConstitution( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ConformationSetSameConstitution::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    std::istream &ConformationSetSameConstitution::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    std::ostream &ConformationSetSameConstitution::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief append a molecule to SetConstitution
    //! @param FRAGMENT fragment conformation  that needs to be added to set constitution
    //! @param CONSTITUTION constitution to link the conformations with
    //! @param ISOMORPHISM isomorphism which if it was determined at constitution layer
    //! @return a pair of iterator to conformation and a bool, true indicating conformation insert was successful
    std::pair< ConformationSetSameConstitution::const_iterator, bool> ConformationSetSameConstitution::Insert
    (
      const ConformationInterface &FRAGMENT,
      const util::ShPtr< FragmentConstitutionShared> &CONSTITUTION,
      storage::Vector< size_t> &ISOMORPHISM
    )
    {
      // insert configuration of the given fragment into constitution set
      std::pair< ConfigurationSetSameConstitution::const_iterator, bool> configuration_insert
      (
        m_ConfigurationSet.Insert( FragmentConfigurationShared( FRAGMENT), CONSTITUTION, ISOMORPHISM)
      );

      // if configuration insertion wasn't successful that means configuration already exists, so insert the conformation
      // in the appropriate ConformationSetSameConfiguration set. Else create a new ConformationSetSameConfiguration and
      // insert in there
      std::pair< ConformationSetSameConfiguration::const_iterator, bool> conformation_insert;

      if( configuration_insert.second)
      {
        ConformationSetSameConfiguration &conformation_set_for_configuration
        (
          m_ConformationsMap[ *configuration_insert.first] =
              ConformationSetSameConfiguration( m_ConfomationComparer, m_Tolerance)
        );
        conformation_insert =
          conformation_set_for_configuration.Insert( FRAGMENT, *configuration_insert.first, ISOMORPHISM);
      }
      else
      {
        conformation_insert =
          m_ConformationsMap[ *configuration_insert.first].Insert( FRAGMENT, *configuration_insert.first, ISOMORPHISM);
      }

      // if conformation_insert was successfull in  ConformationSetSameConfiguration, insert conformation into
      // conformation container
      if( conformation_insert.second)
      {
        m_Conformations.PushBack( *conformation_insert.first);
      }

      return conformation_insert;
    }
  } // namespace chemistry
} // namespace bcl
