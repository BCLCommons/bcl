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
#include "chemistry/bcl_chemistry_fragment_complete.h"

// includes from bcl - sorted alphabetically
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief default constructor
    ConformationSetSameConfiguration::ConformationSetSameConfiguration()
    {
    }

    //! @brief constructor
    //! @param METHOD conformation comparison method to determine if conformations are same
    //! @param TOLERANCE tolerance value for comparing conformation values
    ConformationSetSameConfiguration::ConformationSetSameConfiguration
    (
      util::Implementation< ConformationComparisonInterface> METHOD,
      double TOLERANCE
    ) :
      m_ConfomationComparer( METHOD),
      m_Tolerance( TOLERANCE),
      m_Configuration(),
      m_Conformations()
    {
    }

    //! @brief Clone function
    //! @return pointer to new ConformationSetSameConfiguration
    ConformationSetSameConfiguration *ConformationSetSameConfiguration::Clone() const
    {
      return new ConformationSetSameConfiguration( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ConformationSetSameConfiguration::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ConformationSetSameConfiguration::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ConformationSetSameConfiguration::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief append a molecule to configuration set
    //! @param FRAGMENT fragment conformation shared that needs to be added to configuration set
    //! @param CONFIGURATION configuration of interest
    //! @param ISOMORPHISM isomorphism which if it was determined at configuration layer
    //! @return a pair of iterator to conformation and a bool, true indicating conformation insert was successful
    std::pair< ConformationSetSameConfiguration::const_iterator, bool>
      ConformationSetSameConfiguration::Insert
      (
        const ConformationInterface &FRAGMENT,
        util::ShPtr< FragmentConfigurationShared> CONFIGURATION,
        util::SiPtr< storage::Vector< size_t> > ISOMORPHISM
      )
    {
      // if the configuration of interest has been passed in for the first time then initialize member variables
      if( !m_Configuration.IsDefined())
      {
        // storage the atom and bond info to avoid recalculation
        const storage::Vector< sdf::AtomInfo> atom_info( FRAGMENT.GetAtomInfo());
        const storage::Vector< sdf::BondInfo> bond_info( FRAGMENT.GetBondInfo());
        if( CONFIGURATION.IsDefined())
        {
          m_Configuration = CONFIGURATION;
        }
        else
        {
          // create an atom vector from atom and bond info
          AtomVector< AtomComplete> atoms( atom_info, bond_info);
          if( ISOMORPHISM.IsDefined() && !ISOMORPHISM->IsEmpty())
          {
            atoms.Reorder( *ISOMORPHISM);
          }
          m_Configuration = util::ShPtr< FragmentConfigurationShared>
          (
            new FragmentConfigurationShared( FragmentComplete( atoms, FRAGMENT.GetName()))
          );
        }

        AtomVector< AtomConformationalShared> atoms_conformational( atom_info, bond_info);
        if( ISOMORPHISM.IsDefined() && !ISOMORPHISM->IsEmpty())
        {
          atoms_conformational.Reorder( *ISOMORPHISM);
        }
        util::ShPtr< FragmentConformationShared> shptr_fragment
        (
          new FragmentConformationShared( m_Configuration, atoms_conformational)
        );
        shptr_fragment->StoreProperties( FRAGMENT);

        // Push back fragment into conformation list
        m_Conformations.PushBack( shptr_fragment);
        return std::make_pair( m_Conformations.Last(), true);
      }

      if( !m_Conformations.IsEmpty() && m_Tolerance >= double( 1000))
      {
        return std::make_pair( m_Conformations.Last(), false);
      }

      // check if correct constitution is passed in. This is done in stages.
      // first assert if wrong molecule is passed in the set
      BCL_Assert
      (
        m_Configuration->GetNumberAtoms() == FRAGMENT.GetNumberAtoms()
        && m_Configuration->GetNumberBonds() == FRAGMENT.GetNumberBonds(),
        "Wrong size of molecule passed in"
      );
      // storage the atom and bond info to avoid recalculation
      const storage::Vector< sdf::AtomInfo> atom_info( FRAGMENT.GetAtomInfo());
      const storage::Vector< sdf::BondInfo> bond_info( FRAGMENT.GetBondInfo());

      AtomVector< AtomConformationalShared> atoms_ordered( atom_info, bond_info);

      // if there was a isomorphism required to match on the configurational level, use it to remap the atom vector
      if( ISOMORPHISM.IsDefined() && !ISOMORPHISM->IsEmpty())
      {
        atoms_ordered.Reorder( *ISOMORPHISM);
      }
      util::ShPtr< FragmentConformationShared> shptr_conformation
      (
        new FragmentConformationShared( m_Configuration, atoms_ordered)
      );

      // pushback conformation into conformation container. If comparison method is dihedral bins then gather statistics
      // else compare conformation with conformations that already exist in the set.
      for
      (
        util::ShPtrList< FragmentConformationShared>::const_iterator
          itr( m_Conformations.Begin()), itr_end( m_Conformations.End());
        itr != itr_end;
        ++itr
      )
      {
        if( ( *m_ConfomationComparer)( **itr, *shptr_conformation) <= m_Tolerance)
        {
          return std::make_pair( itr, false);
        }
      }
      shptr_conformation->StoreProperties( FRAGMENT);
      m_Conformations.PushBack( shptr_conformation);
      return std::make_pair( m_Conformations.Last(), true);
    }

  } // namespace chemistry
} // namespace bcl
