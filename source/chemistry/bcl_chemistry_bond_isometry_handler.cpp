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
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_vector.h"
#include "chemistry/bcl_chemistry_substituent_conformational.h"
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief Add E/Z isometry information to the bonds, given a conformation
    //! @param FRAGMENT Conformation upon which to add chirality information
    //! @param FORCE_RECALCULATION whether to relculate even if the bond isometry is already given
    void BondIsometryHandler::AddIsometryInformation
    (
      AtomVector< AtomComplete> &FRAGMENT,
      const bool &FORCE_RECALCULATION
    )
    {
      // if bond isometry must be recalculated, remove all bond isometry information first
      if( FORCE_RECALCULATION)
      {
        for
        (
          AtomVector< AtomComplete>::iterator itr_atom( FRAGMENT.Begin()), itr_atom_end( FRAGMENT.End());
          itr_atom != itr_atom_end;
          ++itr_atom
        )
        {
          for
          (
            storage::Vector< BondConformational>::const_iterator
              itr_bond( itr_atom->GetBonds().Begin()), itr_bond_end( itr_atom->GetBonds().End());
            itr_bond != itr_bond_end;
            ++itr_bond
          )
          {
            // force bonds with non-trivial isometry to have unknown isometry to force the reculculation
            if
            (
              itr_bond->GetBondType()->GetIsometry() != e_NonIsometric
              || itr_bond->GetBondType() == GetConfigurationalBondTypes().e_ConjugatedDoubleBond
              || itr_bond->GetBondType() == GetConfigurationalBondTypes().e_ConjugatedDoubleBondInRing
            )
            {
              // set the bond type; monodirection since all bonds will be iterated over
              itr_atom->SetBondTypeMonoDirectional
              (
                itr_bond->GetTargetAtom(),
                itr_bond->GetBondType()->WithIsometry( e_UnknownIsometry)
              );
            }
          }
        }
      }

      // for each bond in the molecule
      for
      (
        AtomVector< AtomComplete>::iterator itr_atom( FRAGMENT.Begin()), itr_atom_end( FRAGMENT.End());
        itr_atom != itr_atom_end;
        ++itr_atom
      )
      {
        // make a reference on atom_a
        AtomComplete &atom_a( *itr_atom);

        // keep track of whether we need to look for stereoisomers for this atom
        bool atom_could_have_isometric_bonds( true);

        // skip atoms with types that do not have 3 bonds
        if( atom_a.GetAtomType()->GetNumberBonds() != size_t( 3))
        {
          atom_could_have_isometric_bonds = false;
        }

        // skip atoms that do not have at least two explicit bonds
        if( atom_a.GetBonds().GetSize() < size_t( 2))
        {
          atom_could_have_isometric_bonds = false;
        }

        // make a reference on the bonds of atom a
        const storage::Vector< BondConformational> &bonds_a( atom_a.GetBonds());

        // walk over the bonds of this atom
        for
        (
          storage::Vector< BondConformational>::const_iterator itr_bond( bonds_a.Begin()), itr_bond_end( bonds_a.End());
          itr_bond != itr_bond_end;
          ++itr_bond
        )
        {
          // get a reference to the bonded atom
          AtomComplete &atom_b( FRAGMENT( itr_bond->GetTargetAtom()));

          // skip bonds that could not have isometry or for which isometry is obvious
          if( itr_bond->GetBondType()->GetIsometry() != e_UnknownIsometry)
          {
            continue;
          }

          // change stereoisometry of bonds that clearly have no isometry
          if( !atom_could_have_isometric_bonds)
          {
            atom_a.SetBondTypeTo( atom_b, itr_bond->GetBondType()->WithoutIsometry());
            continue;
          }

          // skip if the other atom type does not have exactly 3 bonds
          if( atom_b.GetAtomType()->GetNumberBonds() != size_t( 3))
          {
            atom_a.SetBondTypeTo( atom_b, itr_bond->GetBondType()->WithoutIsometry());
            continue;
          }

          // skip if the other atom type does not have at least two explicitly defined bonds
          if( atom_b.GetBonds().GetSize() < size_t( 2))
          {
            atom_a.SetBondTypeTo( atom_b, itr_bond->GetBondType()->WithoutIsometry());
            continue;
          }

          // make a reference on the bonds of atom b
          const storage::Vector< BondConformational> &bonds_b( atom_b.GetBonds());

          // find the highest priority substituent of the two atoms, excepting each other
          util::SiPtr< const AtomConformationalInterface> highest_priority_substituent_a;

          if( bonds_a.GetSize() == size_t( 2))
          {
            highest_priority_substituent_a = &bonds_a( itr_bond == bonds_a.Begin() ? 1 : 0).GetTargetAtom();

            // the explicit substituent of atom_a must be different from H for this bond to be isometric
            // otherwise, skip this bond
            if( highest_priority_substituent_a->GetElementType() == GetElementTypes().e_Hydrogen)
            {
              atom_a.SetBondTypeTo( atom_b, itr_bond->GetBondType()->WithoutIsometry());
              continue;
            }
          }
          else
          {
            // more general case, need to test which substituent has higher priority
            util::SiPtr< const AtomConformationalInterface> substituent_a_ptr, substituent_b_ptr;

            // get the index of this bond in atom a's bonds
            const size_t bond_index( itr_bond - bonds_a.Begin());

            substituent_a_ptr = bonds_a( bond_index == size_t( 0) ? 1 : 0).GetTargetAtom();
            substituent_b_ptr = bonds_a( bond_index == size_t( 2) ? 1 : 2).GetTargetAtom();

            SubstituentConformational substituent_a( *substituent_a_ptr);
            SubstituentConformational substituent_b( *substituent_b_ptr);

            if( substituent_a < substituent_b)
            {
              highest_priority_substituent_a = *substituent_a_ptr;
            }
            else if( substituent_b < substituent_a)
            {
              highest_priority_substituent_a = *substituent_b_ptr;
            }
            else
            {
              // substituents are equal, no isometry
              atom_a.SetBondTypeTo( atom_b, itr_bond->GetBondType()->WithoutIsometry());
              continue;
            }
          }

          // find the highest priority substituent of the two atoms, excepting each other
          util::SiPtr< const AtomConformationalInterface> highest_priority_substituent_b;

          // get the index of atom_a in atom_b's bonds
          storage::Vector< BondConformational>::const_iterator itr_bonds_b( atom_b.FindBondTo( atom_a));
          if( bonds_b.GetSize() == size_t( 2))
          {
            highest_priority_substituent_b = &bonds_b( itr_bonds_b == bonds_b.Begin() ? 1 : 0).GetTargetAtom();

            // the explicit substituent of atom_b must be different from H for this bond to be isometric
            // otherwise, skip this bond
            if( highest_priority_substituent_b->GetElementType() == GetElementTypes().e_Hydrogen)
            {
              atom_a.SetBondTypeTo( atom_b, itr_bond->GetBondType()->WithoutIsometry());
              continue;
            }
          }
          else
          {
            // more general case, need to test which substituent has higher priority
            util::SiPtr< const AtomConformationalInterface> substituent_a_ptr, substituent_b_ptr;

            // get the index of this bond in atom a's bonds
            const size_t bond_index( itr_bonds_b - bonds_b.Begin());

            substituent_a_ptr = bonds_b( bond_index == size_t( 0) ? 1 : 0).GetTargetAtom();
            substituent_b_ptr = bonds_b( bond_index == size_t( 2) ? 1 : 2).GetTargetAtom();

            SubstituentConformational substituent_a( *substituent_a_ptr);
            SubstituentConformational substituent_b( *substituent_b_ptr);

            if( substituent_a < substituent_b)
            {
              highest_priority_substituent_b = *substituent_a_ptr;
            }
            else if( substituent_b < substituent_a)
            {
              highest_priority_substituent_b = *substituent_b_ptr;
            }
            else
            {
              atom_a.SetBondTypeTo( atom_b, itr_bond->GetBondType()->WithoutIsometry());
              // substituents are equal, no isometry
              continue;
            }
          }

          // Determine the dihedral angles between the highest connected priority substituents
          const double dihedral_angle
          (
            math::Absolute
            (
              linal::Dihedral
              (
                highest_priority_substituent_a->GetPosition(),
                atom_a.GetPosition(),
                atom_b.GetPosition(),
                highest_priority_substituent_b->GetPosition()
              )
            )
          );

          // if the absolute dihedral angle is >= 90 degrees (pi/2), then the bond is E, otherwise it's Z
          if( !util::IsDefined( dihedral_angle))
          {
            // nearly 0 dihedral angle; isometry remains undefined, likely due to bad coordinates
          }
          else if( dihedral_angle >= math::g_Pi / 2.0) // if the absolute dihedral angle is >= 90 degrees (pi/2), then the bond is E, otherwise it's Z
          {
            atom_a.SetBondTypeTo( atom_b, itr_bond->GetBondType()->WithIsometry( e_EIsometry));
          }
          else
          {
            atom_a.SetBondTypeTo( atom_b, itr_bond->GetBondType()->WithIsometry( e_ZIsometry));
          }
        }
      }
    }

  } // namespace chemistry
} // namespace bcl
