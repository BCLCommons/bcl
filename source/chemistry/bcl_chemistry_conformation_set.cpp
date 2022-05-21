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
#include "chemistry/bcl_chemistry_conformation_comparison_by_dihedral_bins.h"
#include "chemistry/bcl_chemistry_conformation_set.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief default constructor
    ConformationSet::ConformationSet()
    :
      m_ConfomationComparer( util::Implementation< ConformationComparisonInterface>( ConformationComparisonByDihedralBins( 30.0))),
      m_Tolerance( 1.0)
    {
    }

    //! @brief constructor
    //! @param METHOD conformation comparison method to determine if conformations are same
    //! @param TOLERANCE tolerance value for comparing conformation values
    ConformationSet::ConformationSet
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
    //! @return pointer to new ConformationSet
    ConformationSet *ConformationSet::Clone() const
    {
      return new ConformationSet( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ConformationSet::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the conformation comparer and tolerance
    void ConformationSet::Setup( util::Implementation< ConformationComparisonInterface> METHOD, double TOLERANCE)
    {
      m_Tolerance = TOLERANCE;
      m_ConfomationComparer = METHOD;
    }

    std::istream &ConformationSet::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    std::ostream &ConformationSet::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief append a molecule to configuration set
    //! @param FRAGMENT fragment conformation shared that needs to be added to configuration set
    //! @return a pair of iterator to conformation and a bool, true indicating conformation insert was successful
    std::pair< ConformationSet::const_iterator, bool> ConformationSet::Insert( const ConformationInterface &FRAGMENT)
    {
      // cache isomorphism between layers
      storage::Vector< size_t> isomorphism;

      // insert configuration in configuration set
      std::pair< ConfigurationSet::const_iterator, bool> configuration_insert
      (
        m_Configurations.Insert( FragmentConfigurationShared( FRAGMENT), isomorphism)
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
          conformation_set_for_configuration.Insert( FRAGMENT, *configuration_insert.first, isomorphism);
      }
      else
      {
        conformation_insert =
          m_ConformationsMap[ *configuration_insert.first].Insert( FRAGMENT, *configuration_insert.first, isomorphism);
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
