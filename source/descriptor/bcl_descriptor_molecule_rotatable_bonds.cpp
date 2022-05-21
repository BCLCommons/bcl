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
#include "descriptor/bcl_descriptor_molecule_rotatable_bonds.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_conformational.h"
#include "chemistry/bcl_chemistry_priority_dihedral_angles.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeRotatableBonds
    MoleculeRotatableBonds *MoleculeRotatableBonds::Clone() const
    {
      return new MoleculeRotatableBonds( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeRotatableBonds::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeRotatableBonds::GetAlias() const
    {
      static const std::string s_name( "NRotBond"), s_symname( "NRotBondSym");
      return m_SymmetryAware ? s_symname : s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeRotatableBonds::Calculate( linal::VectorReference< float> &STORAGE)
    {
      size_t n_rotatable( 0);
      if( !m_SymmetryAware)
      {
        for
        (
          iterate::Generic< const chemistry::AtomConformationalInterface> itr_atom( this->GetCurrentObject()->GetIterator());
          itr_atom.NotAtEnd();
          ++itr_atom
        )
        {
          // if atom is terminal it has only one other atom attached to it
          if( itr_atom->GetBonds().GetSize() <= 1)
          {
            continue;
          }

          // skip atoms that are trivially rotatable; e.g. all but one bond is to a hydrogen
          if( itr_atom->GetNumberCovalentlyBoundHydrogens() >= itr_atom->GetBonds().GetSize() - 1)
          {
            continue;
          }
          const bool is_mandatory_linear
          (
            itr_atom->GetAtomType()->GetNumberBonds() == 2
            && itr_atom->GetAtomType()->GetNumberElectronsInBonds() == 4
          );

          for
          (
            storage::Vector< chemistry::BondConformational>::const_iterator
              itr_bonds( itr_atom->GetBonds().Begin()), itr_bonds_end( itr_atom->GetBonds().End());
            itr_bonds != itr_bonds_end;
            ++itr_bonds
          )
          {
            const chemistry::AtomConformationalInterface &target_atom( itr_bonds->GetTargetAtom());
            if( &target_atom > &*itr_atom)
            {
              continue;
            }
            if
            (
              itr_bonds->GetBondType()->GetBondData( chemistry::ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness)
              != 1
            )
            {
              continue;
            }

            // count hydrogens around this atom
            const size_t h_count_neighbor( target_atom.GetNumberCovalentlyBoundHydrogens());

            // skip atoms that are trivially rotatable; e.g. all but one bond is to a hydrogen
            if( h_count_neighbor >= target_atom.GetBonds().GetSize() - 1)
            {
              continue;
            }

            if
            (
              is_mandatory_linear
              && target_atom.GetAtomType()->GetNumberBonds() == 2
              && target_atom.GetAtomType()->GetNumberElectronsInBonds() == 4
            )
            {
              // if both atoms are linear, then (usually) this is just a linear terminal segment.
              // If not (e.g. the linear bonds both lead to non-linear segments), then it is still true that rotation about
              // this bond is equivalent to rotation about the more distal rotatable bond, in which case this bond doesn't
              // contribute any additional degrees of freedom
              continue;
            }

            // the bond passes the list of rules to be rotatable
            ++n_rotatable;
          }
        }
      }
      else
      {
        util::SiPtr< const chemistry::ConformationInterface> mol_ptr( this->GetCurrentObject());
        storage::Vector< size_t> priorities( chemistry::PriorityDihedralAngles::GetPriority( *mol_ptr));
        for
        (
          iterate::Generic< const chemistry::AtomConformationalInterface> itr_atom( this->GetCurrentObject()->GetIterator());
          itr_atom.NotAtEnd();
          ++itr_atom
        )
        {
          // if atom is terminal it has only one other atom attached to it
          if( itr_atom->GetBonds().GetSize() <= 1)
          {
            continue;
          }

          // skip atoms that are trivially rotatable; e.g. all but one bond is to a hydrogen
          if( itr_atom->GetNumberCovalentlyBoundHydrogens() >= itr_atom->GetBonds().GetSize() - 1)
          {
            continue;
          }
          const bool is_mandatory_linear
          (
            itr_atom->GetAtomType()->GetNumberBonds() == 2
            && itr_atom->GetAtomType()->GetNumberElectronsInBonds() == 4
          );

          bool may_be_symmetric( false);
          size_t symmetric_bond_number( 0);
          if
          (
            itr_atom->GetAtomType()->GetNumberBonds() == size_t( 4)
            ||
            (
              itr_atom->GetAtomType()->GetHybridOrbitalType() == chemistry::GetHybridOrbitalTypes().e_SP2
              && itr_atom->GetAtomType()->GetNumberBonds() == size_t( 3)
            )
          )
          {
            may_be_symmetric = true;
            storage::Vector< size_t> neighbor_priorities;
            for
            (
              storage::Vector< chemistry::BondConformational>::const_iterator
                itr_bonds( itr_atom->GetBonds().Begin()), itr_bonds_end( itr_atom->GetBonds().End());
              itr_bonds != itr_bonds_end;
              ++itr_bonds
            )
            {
              const chemistry::AtomConformationalInterface &target_atom( itr_bonds->GetTargetAtom());
              neighbor_priorities.PushBack( priorities( mol_ptr->GetAtomIndex( target_atom)));
            }
            while( neighbor_priorities.GetSize() < itr_atom->GetAtomType()->GetNumberBonds())
            {
              neighbor_priorities.PushBack( priorities.GetSize() + size_t( 1));
            }
            if( neighbor_priorities( 0) == neighbor_priorities( 1))
            {
              if( neighbor_priorities( 0) == neighbor_priorities( 2))
              {
                symmetric_bond_number = util::GetUndefined< size_t>();
              }
              else
              {
                symmetric_bond_number = 2;
              }
            }
            else if( neighbor_priorities( 0) == neighbor_priorities( 2))
            {
              symmetric_bond_number = 1;
            }
            else if( neighbor_priorities( 1) == neighbor_priorities( 2))
            {
              symmetric_bond_number = 0;
            }
            else
            {
              may_be_symmetric = false;
            }
            if( may_be_symmetric && neighbor_priorities.GetSize() == size_t( 4))
            {
              if( !util::IsDefined( symmetric_bond_number))
              {
                if( neighbor_priorities( 0) != neighbor_priorities( 3))
                {
                  symmetric_bond_number = 3;
                }
              }
              else if( neighbor_priorities( !symmetric_bond_number) != neighbor_priorities( 3))
              {
                may_be_symmetric = false;
              }
            }
          }

          size_t bond_nr( 0);
          for
          (
            storage::Vector< chemistry::BondConformational>::const_iterator
              itr_bonds( itr_atom->GetBonds().Begin()), itr_bonds_end( itr_atom->GetBonds().End());
            itr_bonds != itr_bonds_end;
            ++itr_bonds, ++bond_nr
          )
          {
            const chemistry::AtomConformationalInterface &target_atom( itr_bonds->GetTargetAtom());
            if( &target_atom > &*itr_atom)
            {
              continue;
            }
            if
            (
              itr_bonds->GetBondType()->GetBondData( chemistry::ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness)
              != 1
            )
            {
              continue;
            }
            if( may_be_symmetric && ( !util::IsDefined( symmetric_bond_number) || bond_nr == symmetric_bond_number))
            {
              continue;
            }

            // count hydrogens around this atom
            const size_t h_count_neighbor( target_atom.GetNumberCovalentlyBoundHydrogens());

            // skip atoms that are trivially rotatable; e.g. all but one bond is to a hydrogen
            if( h_count_neighbor >= target_atom.GetBonds().GetSize() - 1)
            {
              continue;
            }

            if
            (
              is_mandatory_linear
              && target_atom.GetAtomType()->GetNumberBonds() == 2
              && target_atom.GetAtomType()->GetNumberElectronsInBonds() == 4
            )
            {
              // if both atoms are linear, then (usually) this is just a linear terminal segment.
              // If not (e.g. the linear bonds both lead to non-linear segments), then it is still true that rotation about
              // this bond is equivalent to rotation about the more distal rotatable bond, in which case this bond doesn't
              // contribute any additional degrees of freedom
              continue;
            }

            if
            (
              target_atom.GetAtomType()->GetNumberBonds() == size_t( 4)
              ||
              (
                target_atom.GetAtomType()->GetHybridOrbitalType() == chemistry::GetHybridOrbitalTypes().e_SP2
                && target_atom.GetAtomType()->GetNumberBonds() == size_t( 3)
              )
            )
            {
              storage::Vector< size_t> neighbor_priorities;
              for
              (
                storage::Vector< chemistry::BondConformational>::const_iterator
                  itr_bonds( target_atom.GetBonds().Begin()), itr_bonds_end( target_atom.GetBonds().End());
                itr_bonds != itr_bonds_end;
                ++itr_bonds
              )
              {
                const chemistry::AtomConformationalInterface &target_atomb( itr_bonds->GetTargetAtom());
                if( &target_atomb != &*itr_atom)
                {
                  neighbor_priorities.PushBack( priorities( mol_ptr->GetAtomIndex( target_atomb)));
                }
              }
              while( neighbor_priorities.GetSize() < target_atom.GetAtomType()->GetNumberBonds() - size_t( 1))
              {
                neighbor_priorities.PushBack( priorities.GetSize() + size_t( 1));
              }
              if
              (
                neighbor_priorities( 0) == neighbor_priorities( 1)
                &&
                ( neighbor_priorities.GetSize() == size_t( 2) || neighbor_priorities( 0) == neighbor_priorities( 2))
              )
              {
                continue;
              }
            }

            // the bond passes the list of rules to be rotatable
            ++n_rotatable;
          }
        }
      }

      STORAGE = n_rotatable;
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeRotatableBonds::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "calculates the number of non-trivially rotatable bonds");
      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
