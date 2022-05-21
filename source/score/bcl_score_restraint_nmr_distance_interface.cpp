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
#include "score/bcl_score_restraint_nmr_distance_interface.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  /////////////////
  // data access //
  /////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief gets the sum of the number of bonds between each atom and the CB
    //! @param ASSIGNMENT AtomDistanceAssignment containing the two atoms
    //! @return the sum of the number of bonds between each atom and the CB
    size_t RestraintNMRDistanceInterface::GetTotalBondsFromCB( const restraint::AtomDistanceAssignment &ASSIGNMENT)
    {
      return GetBondsFromCB( ASSIGNMENT.GetAtomA().GetType()) + GetBondsFromCB( ASSIGNMENT.GetAtomB().GetType());
    }

    //! @brief gets the number of bonds from the cb for a given atom type, H and HA return 0 since their positions can
    //!        be determined
    //! @param ATOM_TYPE atom type
    //! @return the number of bonds from the cb for a given atom type
    size_t RestraintNMRDistanceInterface::GetBondsFromCB( const biol::AtomType &ATOM_TYPE)
    {
      // create static undefined distance
      static const size_t s_undefined_distance( util::GetUndefined< size_t>());

      // if the atom type is CB, H, or HA
      if
      (
        ATOM_TYPE == biol::GetAtomTypes().CB ||
        ATOM_TYPE == biol::GetAtomTypes().H ||
        ATOM_TYPE == biol::GetAtomTypes().HA ||
        ATOM_TYPE == biol::GetAtomTypes().HA2 ||
        ATOM_TYPE == biol::GetAtomTypes().HA3
      )
      {
        return 0;
      }

      // if this is a spin label
      if( ATOM_TYPE == biol::GetAtomTypes().O1)
      {
        return GetSpinLabelLength();
      }

      // if the atom is not a proton
      if( ATOM_TYPE->GetElementType() != chemistry::GetElementTypes().e_Hydrogen)
      {
        // return an undefined distance
        return s_undefined_distance;
      }

      // read in the second character to determine the side chain position (i.e. beta, gamma, delta, etc.)
      const char atom_pos( ATOM_TYPE.GetName()[ 1]);

      // initialize number of bonds between proton and CB
      size_t nr_bonds( s_undefined_distance);

      // use the atom_pos char to determine how many bonds away the proton is
      switch( atom_pos)
      {
        case 'A':
          nr_bonds = 2;
          break;
        case 'B':
          nr_bonds = 1;
          break;
        case 'G':
          nr_bonds = 2;
          break;
        case 'D':
          nr_bonds = 3;
          break;
        case 'E':
          nr_bonds = 4;
          break;
        case 'Z':
          nr_bonds = 5;
          break;
        case 'H':
          nr_bonds = 6;
          break;
        // proton position not determined so return distance will remain undefined
      }

      // end
      return nr_bonds;
    }

    //! @brief Gets the atom type (CB, H, or HA) from the given assignment
    //! @param RESTRAINT assignment containing the atoms
    //! @return the atom type
    biol::AtomType RestraintNMRDistanceInterface::GetAtomType( const restraint::AtomDistanceAssignment &RESTRAINT)
    {
      // store types
      const biol::AtomType &type_a( RESTRAINT.GetAtomA().GetType());
      const biol::AtomType &type_b( RESTRAINT.GetAtomB().GetType());

      if
      (
        type_a->GetElementType() == chemistry::GetElementTypes().e_Hydrogen
        || type_b->GetElementType() == chemistry::GetElementTypes().e_Hydrogen
      )
      {
        // if either atom type is H
        if( type_a == biol::GetAtomTypes().H || type_b == biol::GetAtomTypes().H)
        {
          return biol::GetAtomTypes().H;
        }

        // if either type is HA
        if
        (
          type_a == biol::GetAtomTypes().HA || type_a == biol::GetAtomTypes().HA2 || type_a == biol::GetAtomTypes().HA3 ||
          type_b == biol::GetAtomTypes().HA || type_b == biol::GetAtomTypes().HA2 || type_b == biol::GetAtomTypes().HA3
        )
        {
          return biol::GetAtomTypes().HA;
        }
      }

      // type is not H or HA so return CB
      return biol::GetAtomTypes().CB;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RestraintNMRDistanceInterface::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        " an NMR distance restraint score. "
        "Currently these are not enumerated, so there are no user-changeable parameters"
      );
      return serializer;
    }

    //! @brief Get the number of bonds that the spin label is from CB from the command line flag
    //! @return the number of bonds that the spin label is from CB from the command line flag
    size_t &RestraintNMRDistanceInterface::GetSpinLabelLength()
    {
      static size_t s_spin_label_length( 6);
      return s_spin_label_length;
    }

  } // namespace score

} // namespace bcl
