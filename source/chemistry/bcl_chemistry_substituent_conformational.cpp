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
#include "chemistry/bcl_chemistry_substituent_conformational.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_conformational.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief Compare two atoms based on atomic mass (to distinguish isotopes when they are added) and, for atoms of equal mass, compares # H (counting valences)
    //! @param ATOM_A, ATOM_B the two SiPtr< const Atom> to compare
    //! @return true if ATOM_A has higher atomic number than ATOM_B
    bool AtomicMassIsGreaterOrHasMoreH
    (
      const util::SiPtr< const AtomConformationalInterface> &ATOM_A,
      const util::SiPtr< const AtomConformationalInterface> &ATOM_B
    )
    {
      const double mass_a( ATOM_A->GetElementType()->GetProperty( ElementTypeData::e_Mass));
      const double mass_b( ATOM_B->GetElementType()->GetProperty( ElementTypeData::e_Mass));
      if( mass_a > mass_b)
      {
        return true;
      }
      else if( mass_a < mass_b)
      {
        return false;
      }

      if( ATOM_A->GetAtomType()->IsGasteigerAtomType() && ATOM_B->GetAtomType()->IsGasteigerAtomType())
      {
        return ATOM_A->GetNumberValenceBonds() + ATOM_A->GetNumberCovalentlyBoundHydrogens()
               > ATOM_B->GetNumberValenceBonds() + ATOM_B->GetNumberCovalentlyBoundHydrogens();
      }
      return ATOM_A->GetNumberCovalentlyBoundHydrogens() > ATOM_B->GetNumberCovalentlyBoundHydrogens();
    }

    //! @brief test if one vector of SubstituentConformationals is less than another
    //! @param LHS, RHS the left and right hand vector of atoms with priority
    //! @return 1 if LHS > RHS, 0 if LHS = RHS, -1 if LHS < RHS
    int CompareSubstituentConformationalVectors
    (
      const util::SiPtrVector< const AtomConformationalInterface> &LHS,
      const util::SiPtrVector< const AtomConformationalInterface> &RHS
    )
    {
      // Compare the vectors of all bonded atoms, sorted in descending order of atomic mass
      for
      (
        util::SiPtrVector< const AtomConformationalInterface>::const_iterator
          itr_atoms_a( LHS.Begin()), itr_atoms_b( RHS.Begin()),
          itr_atoms_a_end( LHS.End()), itr_atoms_b_end( RHS.End());
        itr_atoms_a != itr_atoms_a_end && itr_atoms_b != itr_atoms_b_end;
        ++itr_atoms_a, ++itr_atoms_b
      )
      {
        if( AtomicMassIsGreaterOrHasMoreH( *itr_atoms_a, *itr_atoms_b))
        {
          return 1;
        }

        if( AtomicMassIsGreaterOrHasMoreH( *itr_atoms_b, *itr_atoms_a))
        {
          return -1;
        }
      }

      //If the compared atomic numbers are equal, check to see if there are still some left in only one vector
      if( LHS.GetSize() > RHS.GetSize())
      {
        return 1;
      }

      if( LHS.GetSize() < RHS.GetSize())
      {
        return -1;
      }

      // there is nothing to differentiate the vectors, so they have equal priority
      return 0;
    }

  //////////
  // data //
  //////////

    const util::SiPtr< const util::ObjectInterface> SubstituentConformational::s_Instance
    (
      GetObjectInstances().AddInstance( new SubstituentConformational())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SubstituentConformational::SubstituentConformational()
    {
    }

    //! @brief construct a SubstituentConformational from an ATOM
    SubstituentConformational::SubstituentConformational
    (
      const AtomConformationalInterface &ATOM
    )
    {
      m_Substituents.PushBack( util::SiPtrVector< const AtomConformationalInterface>( ATOM));
    }

    //! @brief virtual copy constructor
    SubstituentConformational *SubstituentConformational::Clone() const
    {
      return new SubstituentConformational( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &SubstituentConformational::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the atom at the root of this search
    //! @return the atom at the root of this search
    const util::SiPtr< const AtomConformationalInterface> &SubstituentConformational::GetRootAtom() const
    {
      return m_Substituents.FirstElement().FirstElement();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief test if one SubstituentConformational is less than another
    //! @param RHS right hand SubstituentConformationals
    //! @return true if this has higher priority than RHS
    bool SubstituentConformational::operator <( const SubstituentConformational &RHS) const
    {
      storage::List< util::SiPtrVector< const AtomConformationalInterface> >::const_iterator
        itr_lhs( m_Substituents.Begin()), itr_lhs_end( m_Substituents.End()),
        itr_rhs( RHS.m_Substituents.Begin()), itr_rhs_end( RHS.m_Substituents.End());

      // compare the known SubstituentConformationals of each SubstituentConformational
      for
      (
        ;
        itr_lhs != itr_lhs_end && itr_rhs != itr_rhs_end;
        ++itr_lhs, ++itr_rhs
      )
      {
        int comparison_result( CompareSubstituentConformationalVectors( *itr_lhs, *itr_rhs));
        if( comparison_result != 0)
        {
          return comparison_result > 0;
        }
      }

      // the LHS and RHS agree up to the points where they have been expanded

      // expand the RHS until it is as big as the LHS or it cannot be expanded any farther
      while( m_Substituents.GetSize() > RHS.m_Substituents.GetSize() && !RHS.IsComplete())
      {
        RHS.Expand();
        const int comparison_result( CompareSubstituentConformationalVectors( *itr_lhs++, RHS.m_Substituents.LastElement()));
        if( comparison_result != 0)
        {
          return comparison_result > 0;
        }
      }

      // do the same for the LHS
      while( RHS.m_Substituents.GetSize() > m_Substituents.GetSize() && !IsComplete())
      {
        Expand();
        const int comparison_result( CompareSubstituentConformationalVectors( m_Substituents.LastElement(), *itr_rhs++));
        if( comparison_result != 0)
        {
          return comparison_result > 0;
        }
      }

      // expand both sides simultaneously, now that they have the same number of SubstituentConformational vectors
      // and are still equal up to this point
      while( !RHS.IsComplete() && !IsComplete())
      {
        Expand();
        RHS.Expand();
        const int comparison_result( CompareSubstituentConformationalVectors( m_Substituents.LastElement(), RHS.m_Substituents.LastElement()));
        if( comparison_result != 0)
        {
          return comparison_result > 0;
        }
      }

      // one of the two sides is complete, so we need only compare the # of atoms seen to determine
      // if LHS belongs before RHS
      return m_ExpandedAtoms.GetSize() > RHS.m_ExpandedAtoms.GetSize();
    }

    //! @brief tell whether or not the SubstituentConformational is complete, e.g. whether it has any leaves that can be expanded
    //! @return true if there are
    bool SubstituentConformational::IsComplete() const
    {
      return m_Substituents.LastElement().IsEmpty();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SubstituentConformational::Read( std::istream &ISTREAM)
    {
      BCL_Exit( "Read is not available for " + GetClassIdentifier(), 1);

      // end
      return ISTREAM;
    }

    //! @brief read from std::ostream
    //! @param OSTREAM input stream
    //! @return ostream which was read from
    std::ostream &SubstituentConformational::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Substituents, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PriorGhosts, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ExpandedAtoms, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief add the SubstituentConformationals at the next level of depth
    //! @return the SubstituentConformationals at the next level of description
    const util::SiPtrVector< const AtomConformationalInterface> &SubstituentConformational::Expand() const
    {
      if( IsComplete())
      {
        return m_Substituents.LastElement();
      }

      // find the 'roots' of the breadth first search, so far as it has been conducted
      const util::SiPtrVector< const AtomConformationalInterface> &roots( m_Substituents.LastElement());

      // add a new vector for the subsituents' atoms we find at the next depth; initialize it with ghosted atoms from last time
      m_Substituents.PushBack( m_PriorGhosts);

      // get a reference on the new SubstituentConformationals
      util::SiPtrVector< const AtomConformationalInterface> &new_substituents( m_Substituents.LastElement());

      // reset the ghost vector
      m_PriorGhosts.Reset();

      for
      (
        util::SiPtrVector< const AtomConformationalInterface>::const_iterator
          itr_roots( roots.Begin()), itr_roots_end( roots.End());
        itr_roots != itr_roots_end;
        ++itr_roots
      )
      {
        // get the pointer to the atom
        const util::SiPtr< const AtomConformationalInterface> &si_atom( *itr_roots);

        // if the atom already had its connectivity expanded, do not try to expand it again
        if( m_ExpandedAtoms.Contains( *itr_roots))
        {
          continue;
        }

        // get a pointer to this atoms bonds
        const storage::Vector< BondConformational> &bonds( ( *itr_roots)->GetBonds());

        //  Go through all the bonded atoms to current root and add to growing list of atoms at current branch level.
        for
        (
          storage::Vector< BondConformational>::const_iterator itr_bonds( bonds.Begin()), itr_bonds_end( bonds.End());
          itr_bonds != itr_bonds_end;
          ++itr_bonds
        )
        {
          util::SiPtr< const AtomConformationalInterface> bonded_atom( itr_bonds->GetTargetAtom());
          if( bonded_atom->GetElementType() == GetElementTypes().e_Hydrogen)
          {
            continue;
          }
          new_substituents.PushBack( bonded_atom);
          const ConfigurationalBondType bond_type( itr_bonds->GetBondType());

          if
          (
            bond_type->GetBondData( ConfigurationalBondTypeData::e_IsAromatic)
            || bond_type->GetNumberOfElectrons() == 4
          )
          { // Based on rules of stereochemistry, a double-bonded atom is counted twice,
            // aromatic bonds are here considered the same as double bonds, thus aromatic rings have higher priority
            // than non-aromatic rings
            new_substituents.PushBack( bonded_atom);
            m_PriorGhosts.PushBack( si_atom);
          }
          else if( bond_type->GetNumberOfElectrons() == 6)
          {
            // and an atom bound by a triple bond is counted three times
            new_substituents.PushBack( bonded_atom);
            m_PriorGhosts.PushBack( si_atom);
            new_substituents.PushBack( bonded_atom);
            m_PriorGhosts.PushBack( si_atom);
          }
        }

        // add the atom to the list of atoms which have already been expanded
        m_ExpandedAtoms.Insert( si_atom);
      }

      // Sort the two vectors of all bonded atoms so that their atomic numbers can be compared pairwise.
      new_substituents.Sort( &AtomicMassIsGreaterOrHasMoreH);

      return new_substituents;
    }

  } // namespace chemistry
} // namespace bcl

