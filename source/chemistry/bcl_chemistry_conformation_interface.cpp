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
#include "chemistry/bcl_chemistry_conformation_interface.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_conformation_comparison_by_dihedral_bins.h"
#include "chemistry/bcl_chemistry_fragment_split_unbridged_rings.h"
#include "linal/bcl_linal_matrix3x3.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_angle.h"
#include "math/bcl_math_linear_least_squares.h"
#include "math/bcl_math_running_average_sd.h"
#include "sdf/bcl_sdf_mdl_handler.h"
// external includes - sorted alphabetically

// Uncomment the following line to see a printout of ring planarity by type; useful for benchmarking aromatization schemes
//#define BCL_ShowAromaticRingPlanarityDetail

namespace bcl
{
  namespace chemistry
  {

  /////////////////
  // data access //
  /////////////////

    //! @brief return a number of hydrogens in the conformation
    //! @return the number of hydrogens in the conformation
    size_t ConformationInterface::GetNumberHydrogens() const
    {
      size_t number_h( 0);
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        number_h += itr_atoms->GetElementType() == GetElementTypes().e_Hydrogen;
      }
      return number_h;
    }

    //! @brief return the number of bonds with a particular property (e.g. is in ring, is aromatic, is isometric)
    //! @param DATA the data to retrieve for each bond
    //! @param VALUE the value to count, e.g. if ConfigurationalBondType->GetData( DATA) == VALUE, ++return
    //! @return the number of bonds with the property of interest
    size_t ConformationInterface::CountNonValenceBondsWithProperty
    (
      const ConfigurationalBondTypeData::Data &DATA,
      const size_t &VALUE
    ) const
    {
      size_t number_bonds( 0);
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        number_bonds += itr_atoms->CountNonValenceBondsWithProperty( DATA, VALUE);
      }
      number_bonds /= 2;
      return number_bonds;
    }

    //! @brief get the number of valences for the entire molecule or fragment
    //! @return the number of valences on the entire molecule or fragment
    size_t ConformationInterface::GetNumberValences() const
    {
      size_t number_valences( 0);

      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        number_valences += itr_atoms->GetNumberValenceBonds();
      }
      return number_valences;
    }

    //! @return a vector with the atom types in the molecule
    storage::Vector< AtomType> ConformationInterface::GetAtomTypesVector() const
    {
      storage::Vector< AtomType> atom_types;
      atom_types.AllocateMemory( GetNumberAtoms());

      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        atom_types.PushBack( itr_atoms->GetAtomType());
      }
      return atom_types;
    }

    //! @brief get a string with each atom type in the molecule
    //! @return a string with each atom type in the molecule
    std::string ConformationInterface::GetAtomTypesString() const
    {
      std::string atom_types;

      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        atom_types += itr_atoms->GetAtomType().GetName() + " ";
      }
      return atom_types;
    }

    //! @return a string with each element type in the molecule
    std::string ConformationInterface::GetElementTypesString() const
    {
      std::string atom_types;

      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        atom_types += itr_atoms->GetElementType()->GetChemicalSymbol() + " ";
      }
      return atom_types;
    }

    //! @brief get a string with each bond type in the molecule
    //! @return a string with each bond type in the molecule
    std::string ConformationInterface::GetBondTypesString() const
    {
      std::ostringstream bond_types;
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr( GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        const storage::Vector< BondConformational> &connected_atoms( itr->GetBonds());

        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_connected( connected_atoms.Begin()),
            itr_connected_end( connected_atoms.End());
          itr_connected != itr_connected_end;
          ++itr_connected
        )
        {
          if( &*itr < &itr_connected->GetTargetAtom())
          {
            bond_types << itr_connected->GetBondType().GetName() << ' ';
          }
        }
      }

      return bond_types.str();
    }

    //! @brief return a string containing the double bond isometry of the atoms of this molecule as a string
    std::string ConformationInterface::GetDoubleBondIsometryString() const
    {
      std::ostringstream bond_types;
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr( GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        const storage::Vector< BondConformational> &connected_atoms( itr->GetBonds());

        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_connected( connected_atoms.Begin()),
            itr_connected_end( connected_atoms.End());
          itr_connected != itr_connected_end;
          ++itr_connected
        )
        {
          if( &*itr < &itr_connected->GetTargetAtom() && itr_connected->GetBondType()->GetNumberOfElectrons() == size_t( 4))
          {
            bond_types << GetIsometryName( itr_connected->GetBondType()->GetIsometry()) << ' ';
          }
        }
      }

      return bond_types.str();
    }

    //! @brief get a string containing stereochemistry of atoms
    //! @brief return a string containing the stereochemistry of the atoms of this molecule as a string
    std::string ConformationInterface::GetChiralityString() const
    {
      std::ostringstream stereochem;
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr( GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        stereochem << GetChiralityName( itr->GetChirality()) << ' ';
      }

      return stereochem.str();
    }

    //! @brief returns a string with the sum formula
    //! @return a string with the sum formula
    std::string ConformationInterface::GetSumFormula() const
    {
      // vector of positive integers that count the number of atoms that appear in the molecule
      storage::Vector< size_t> formula( GetElementTypes().GetEnumCount(), size_t( 0));

      // number of unknown elements
      size_t number_unknown( 0);

      // iterate through all atoms of the molecule and count their appearance
      for( iterate::Generic< const AtomConformationalInterface> itr( GetAtomsIterator()); itr.NotAtEnd(); ++itr)
      {
        if( itr->GetElementType().IsDefined())
        {
          ++formula( itr->GetElementType().GetIndex());
        }
        else
        {
          ++number_unknown;
        }
      }

      // stringstream for the final sum formula
      std::stringstream sum_formula( "");

      // sum formula begins always with a carbon if there is one
      if( formula( GetElementTypes().e_Carbon) != 0)
      {
        sum_formula << GetElementTypes().e_Carbon->GetChemicalSymbol();

        // only print the number of atoms if it is greater than 1
        if( formula( GetElementTypes().e_Carbon) > 1)
        {
          sum_formula << formula( GetElementTypes().e_Carbon);
        }
      }

      // assemble sum formula with delimiter "_"
      for
      (
        ElementTypes::const_iterator itr( GetElementTypes().Begin()), itr_end( GetElementTypes().End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->IsDefined() && formula( *itr) != 0 && ( *itr) != GetElementTypes().e_Carbon)
        {
          sum_formula << ( *itr)->GetChemicalSymbol();

          // only print the number of atoms if it is greater than 1
          if( formula( *itr) > 1)
          {
            sum_formula << formula( *itr);
          }
        }
      }

      if( number_unknown != 0)
      {
        sum_formula << "X";

        // only print the number of atoms if it is greater than 1
        if( number_unknown > 1)
        {
          sum_formula << number_unknown;
        }
      }

      // return sum formula
      return sum_formula.str();
    }

    //! @brief assignment operator
    ConformationInterface &ConformationInterface::operator=( const ConformationInterface &C)
    {
      descriptor::SequenceInterface< AtomConformationalInterface>::operator=( C);
      HasPropertiesInterface< coord::MovableInterface>::operator=( C);
      m_ChangeSignal.Emit( *this);
      return *this;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief checks for 'bad' geometry, i.e., z-coord./all coords. == const.
    //! @return flag whether geometry is 'bad'
    storage::Vector< size_t> ConformationInterface::GetAtomsWithBadGeometry() const
    {
      storage::Set< size_t> atoms_bad_geo;
      iterate::Generic< const AtomConformationalInterface> atoms( this->GetAtomsIterator());
      // loop over all atoms in molecule, checking the following conditions for a molecule with valid geometry:
      // - all coordinates are defined
      // - atoms do not overlap (e.g. no closer than their average covalent radii)
      // - atoms with 4+ neighbors should have at least 1 neighbor w/ a non-0 z coordinate
      // - SP3 atoms with 3 neighbors should have at least 1 neighbor w/ a non-0 z coordinate

      // note: the warnings below are verbose since this function is called while loading in molecules to detect
      // the presence of 3d geometries and it is not necessarily an error when reading in configurations

      // first, check that all coordinates are defined, and obtain atomic radii for all atoms
      storage::Vector< double> atom_covalent_radii, atom_vdw_radii;
      atom_covalent_radii.AllocateMemory( atoms.GetSize());
      atom_vdw_radii.AllocateMemory( atoms.GetSize());
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( atoms.Begin());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        // check that coordinates are defined
        if( !itr_atoms->GetPosition().IsDefined())
        {
          // undefined coordinates; geometry is bad
          BCL_MessageVrb
          (
            "Bad Geometry! Atom with undefined position (type: " + itr_atoms->GetAtomType().GetName() + ")"
          );
          atoms_bad_geo.Insert( itr_atoms.GetPosition());
        }

        // look for overlapping atoms (defined as atoms closer than the average of the covalent radii)
        const double atom_type_covalent_radius( BondLengths::GetAverageCovalentRadius( *itr_atoms));
        const double covalent_radius
        (
          util::IsDefined( atom_type_covalent_radius)
          ? atom_type_covalent_radius
          : itr_atoms->GetElementType()->GetProperty( ElementTypeData::e_CovalentRadius)
        );
        atom_covalent_radii.PushBack( covalent_radius);

        double csd_vdw_radius( itr_atoms->GetAtomType()->GetAtomTypeProperty( AtomTypeData::e_VdWaalsRadiusCSD));

        // only keep vdw radii for atoms outside rings, since the vdw radii are incorrect otherwise
        // also, ignore the vdw radii for atoms with valences and hydrogens because hydrogen bonds shorten the distance
        // between atoms, rendering the notion of a proper vdw radius useless
        if( util::IsDefined( csd_vdw_radius))
        {
          if
          (
            itr_atoms->GetNumberValenceBonds()
            || itr_atoms->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, 1)
            || itr_atoms->GetNumberCovalentlyBoundHydrogens()
          )
          {
            csd_vdw_radius = util::GetUndefined< double>();
          }
        }
        atom_vdw_radii.PushBack( csd_vdw_radius);
      }

      storage::Vector< double>::const_iterator
        itr_atoms_cov_radius( atom_covalent_radii.Begin()), itr_atoms_vdw_radius( atom_vdw_radii.Begin());
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( atoms.Begin());
        itr_atoms.NotAtEnd();
        ++itr_atoms, ++itr_atoms_cov_radius, ++itr_atoms_vdw_radius
      )
      {
        // skip metals; their ionic bonds have much less regular distances and so the checks below are meaningless
        // for them
        if
        (
          itr_atoms->GetElementType() != GetElementTypes().e_Hydrogen
          && ( itr_atoms->GetElementType()->GetMainGroup() == 1 || itr_atoms->GetElementType()->GetMainGroup() == 2)
        )
        {
          continue;
        }

        // look for overlapping atoms (non-bonded atoms closer than their covalent radius sum; bonded atoms closer than
        // their average covalent radius)
        storage::Vector< double>::const_iterator
          itr_atoms_cov_radius_b( atom_covalent_radii.Begin()), itr_atoms_vdw_radius_b( atom_vdw_radii.Begin());
        for
        (
          iterate::Generic< const AtomConformationalInterface> itr_atoms_b( atoms.Begin());
          itr_atoms_b != itr_atoms;
          ++itr_atoms_b, ++itr_atoms_cov_radius_b, ++itr_atoms_vdw_radius_b
        )
        {
          // skip metals; their ionic bonds have much less regular distances and so the checks below are meaningless
          // for them
          if
          (
            itr_atoms_b->GetElementType() != GetElementTypes().e_Hydrogen
            && ( itr_atoms_b->GetElementType()->GetMainGroup() == 1 || itr_atoms_b->GetElementType()->GetMainGroup() == 2)
          )
          {
            continue;
          }

          // compute the bond length (nominally = the sum of the covalent radii)
          double bond_length( *itr_atoms_cov_radius + *itr_atoms_cov_radius_b);

          // compute the minimum distance for these atoms to be apart, assuming they meet all the requirements for the
          // vdw to be valid, as was checked earlier, except the requirement that they not have opposite charge
          const double vdw_min_distance
          (
            itr_atoms->GetCharge() != itr_atoms_b->GetCharge() && itr_atoms->GetCharge() && itr_atoms_b->GetCharge()
            ? util::GetUndefined< double>()
            : *itr_atoms_cov_radius_b + *itr_atoms_cov_radius
          );

          // skip atom types for which the vdw distance and the bond length is undefined
          if( !util::IsDefined( vdw_min_distance) && !util::IsDefined( bond_length))
          {
            continue;
          }

          // compute the actual distance between the two atoms
          const double distance( linal::Distance( itr_atoms->GetPosition(), itr_atoms_b->GetPosition()));

          double max_bond_distance_ratio( 1.25);
          // period 3 elements can have longer bonds due to d-orbital bonding; so allow a larger tolerance for them
          if( itr_atoms->GetElementType()->GetPeriod() > 2)
          {
            max_bond_distance_ratio += 0.125;
          }

          if( itr_atoms_b->GetElementType()->GetPeriod() > 2)
          {
            max_bond_distance_ratio += 0.125;
          }

          storage::Vector< BondConformational>::const_iterator itr_bond( itr_atoms->FindBondTo( *itr_atoms_b));
          const bool are_bonded( itr_bond != itr_atoms->GetBonds().End());
          if( are_bonded)
          {
            const double expected_bond_length
            (
              BondLengths::GetBondLength
              (
                itr_atoms->GetAtomType(),
                itr_bond->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromatic),
                itr_atoms_b->GetAtomType()
              )
            );
            if( util::IsDefined( expected_bond_length))
            {
              bond_length = expected_bond_length;
            }
          }

          // check whether the atoms are bonded
          if( are_bonded)
          {
            // bonded atom, check for bond lengths less than half the normal values
            if( util::IsDefined( bond_length))
            {
              if( distance < 0.75 * bond_length)
              {
                BCL_MessageVrb
                (
                  "Bad Geometry! Too short a bond between atoms of type: " + itr_atoms->GetAtomType().GetName()
                  + " and " + itr_atoms_b->GetAtomType().GetName() + ", distance: "
                  + util::Format()( distance) + " expected bond length: " + util::Format()( bond_length)
                );
                atoms_bad_geo.Insert( itr_atoms.GetPosition());
                atoms_bad_geo.Insert( itr_atoms_b.GetPosition());
              }
              else if( distance > max_bond_distance_ratio * bond_length)
              {
                BCL_MessageVrb
                (
                  "Bad Geometry! Too long a bond between atoms of type: " + itr_atoms->GetAtomType().GetName()
                  + " and " + itr_atoms_b->GetAtomType().GetName() + ", distance: "
                  + util::Format()( distance) + " expected bond length: " + util::Format()( bond_length)
                );
                atoms_bad_geo.Insert( itr_atoms.GetPosition());
                atoms_bad_geo.Insert( itr_atoms_b.GetPosition());
              }
            }
          }
          // if the distance is <= the expected bond length
          else if( util::IsDefined( bond_length) && distance <= bond_length)
          {
            BCL_MessageVrb
            (
              "Bad Geometry! Non-bonded atom (type: " + itr_atoms->GetAtomType().GetName()
              + ") within bonding distance of another atom (type: " + itr_atoms_b->GetAtomType().GetName() + "), distance: "
              + util::Format()( distance) + " expected bond length: " + util::Format()( bond_length)
            );
            atoms_bad_geo.Insert( itr_atoms.GetPosition());
            atoms_bad_geo.Insert( itr_atoms_b.GetPosition());
          }
          else if( util::IsDefined( vdw_min_distance) && distance < vdw_min_distance)
          {
            // check for
            BCL_MessageVrb
            (
              "Bad Geometry! Non-bonded atom (type: " + itr_atoms->GetAtomType().GetName()
              + ") inside van-der waals sphere of another atom (type: " + itr_atoms_b->GetAtomType().GetName() + "), distance: "
              + util::Format()( distance) + " vdw distance: " + util::Format()( vdw_min_distance)
            );
            atoms_bad_geo.Insert( itr_atoms.GetPosition());
            atoms_bad_geo.Insert( itr_atoms_b.GetPosition());
          }
        }

        // for the 3D check, skip atoms with a defined z coordinate
        if( !math::EqualWithinTolerance( itr_atoms->GetPosition().Z(), 0))
        {
          continue;
        }

        const storage::Vector< BondConformational> &connected_atoms( itr_atoms->GetBonds());

        // atoms with 0-2 bonds could have all neighboring atoms with a 0 z-coordinate
        if( connected_atoms.GetSize() < 3)
        {
          continue;
        }

        // atoms with 3 bonds could have all neighboring atoms with a 0 z-coordinate if the geometry is not SP3
        if
        (
          connected_atoms.GetSize() == 3
          && itr_atoms->GetAtomType()->GetHybridOrbitalType() != GetHybridOrbitalTypes().e_SP3
        )
        {
          continue;
        }

        bool had_3d_neighbor( false);
        // check that at least one neighbor has a non-zero z-coordinate
        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_connected( connected_atoms.Begin()),
            itr_connected_end( connected_atoms.End());
          itr_connected != itr_connected_end;
          ++itr_connected
        )
        {
          if( !math::EqualWithinTolerance( itr_connected->GetTargetAtom().GetPosition().Z(), 0))
          {
            had_3d_neighbor = true;
            break;
          }
        }
        if( !had_3d_neighbor)
        {
          const std::string &name( itr_atoms->GetAtomType().GetName());
          if( connected_atoms.GetSize() == 3)
          {
            BCL_MessageVrb
            (
              "Bad Geometry! SP3 atom with (type: " + name + "); all neighbors have 0 z-coordinate"
            );
          }
          else
          {
            BCL_MessageVrb
            (
              "Bad Geometry! Atom with 4 bonds (type: " + name + "); all neighbors have 0 z-coordinate"
            );
          }
          atoms_bad_geo.Insert( itr_atoms.GetPosition());
        }
      }
      return storage::Vector< size_t>( atoms_bad_geo.Begin(), atoms_bad_geo.End());
    }

    //! @brief checks for 'bad' geometry, i.e., z-coord./all coords. == const.
    //! @return flag whether geometry is 'bad'
    bool ConformationInterface::HasBadGeometry() const
    {
      iterate::Generic< const AtomConformationalInterface> atoms( this->GetAtomsIterator());
      // loop over all atoms in molecule, checking the following conditions for a molecule with valid geometry:
      // - all coordinates are defined
      // - atoms do not overlap (e.g. no closer than their average covalent radii)
      // - atoms with 4+ neighbors should have at least 1 neighbor w/ a non-0 z coordinate
      // - SP3 atoms with 3 neighbors should have at least 1 neighbor w/ a non-0 z coordinate

      // note: the warnings below are verbose since this function is called while loading in molecules to detect
      // the presence of 3d geometries and it is not necessarily an error when reading in configurations

      // first, check that all coordinates are defined, and obtain atomic radii for all atoms
      storage::Vector< double> atom_covalent_radii, atom_vdw_radii;
      atom_covalent_radii.AllocateMemory( atoms.GetSize());
      atom_vdw_radii.AllocateMemory( atoms.GetSize());
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( atoms.Begin());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        // check that coordinates are defined
        if( !itr_atoms->GetPosition().IsDefined())
        {
          // undefined coordinates; geometry is bad
          BCL_MessageVrb
          (
            "Bad Geometry! Atom with undefined position (type: " + itr_atoms->GetAtomType().GetName() + ")"
          );
          return true;
        }

        // look for overlapping atoms (defined as atoms closer than the average of the covalent radii)
        const double atom_type_covalent_radius( BondLengths::GetAverageCovalentRadius( *itr_atoms));
        const double covalent_radius
        (
          util::IsDefined( atom_type_covalent_radius)
          ? atom_type_covalent_radius
          : itr_atoms->GetElementType()->GetProperty( ElementTypeData::e_CovalentRadius)
        );
        atom_covalent_radii.PushBack( covalent_radius);

        double csd_vdw_radius( itr_atoms->GetAtomType()->GetAtomTypeProperty( AtomTypeData::e_VdWaalsRadiusCSD));

        // only keep vdw radii for atoms outside rings, since the vdw radii are incorrect otherwise
        // also, ignore the vdw radii for atoms with valences and hydrogens because hydrogen bonds shorten the distance
        // between atoms, rendering the notion of a proper vdw radius useless
        if( util::IsDefined( csd_vdw_radius))
        {
          if
          (
            itr_atoms->GetNumberValenceBonds()
            || itr_atoms->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, 1)
            || itr_atoms->GetNumberCovalentlyBoundHydrogens()
          )
          {
            csd_vdw_radius = util::GetUndefined< double>();
          }
        }
        atom_vdw_radii.PushBack( csd_vdw_radius);
      }

      storage::Vector< double>::const_iterator
        itr_atoms_cov_radius( atom_covalent_radii.Begin()), itr_atoms_vdw_radius( atom_vdw_radii.Begin());
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( atoms.Begin());
        itr_atoms.NotAtEnd();
        ++itr_atoms, ++itr_atoms_cov_radius, ++itr_atoms_vdw_radius
      )
      {
        // skip metals; their ionic bonds have much less regular distances and so the checks below are meaningless
        // for them
        if
        (
          itr_atoms->GetElementType() != GetElementTypes().e_Hydrogen
          && ( itr_atoms->GetElementType()->GetMainGroup() == 1 || itr_atoms->GetElementType()->GetMainGroup() == 2)
        )
        {
          continue;
        }

        // look for overlapping atoms (non-bonded atoms closer than their covalent radius sum; bonded atoms closer than
        // their average covalent radius)
        storage::Vector< double>::const_iterator
          itr_atoms_cov_radius_b( atom_covalent_radii.Begin()), itr_atoms_vdw_radius_b( atom_vdw_radii.Begin());
        for
        (
          iterate::Generic< const AtomConformationalInterface> itr_atoms_b( atoms.Begin());
          itr_atoms_b != itr_atoms;
          ++itr_atoms_b, ++itr_atoms_cov_radius_b, ++itr_atoms_vdw_radius_b
        )
        {
          // skip metals; their ionic bonds have much less regular distances and so the checks below are meaningless
          // for them
          if
          (
            itr_atoms_b->GetElementType() != GetElementTypes().e_Hydrogen
            && ( itr_atoms_b->GetElementType()->GetMainGroup() == 1 || itr_atoms_b->GetElementType()->GetMainGroup() == 2)
          )
          {
            continue;
          }

          // compute the bond length (nominally = the sum of the covalent radii)
          double bond_length( *itr_atoms_cov_radius + *itr_atoms_cov_radius_b);

          // compute the minimum distance for these atoms to be apart, assuming they meet all the requirements for the
          // vdw to be valid, as was checked earlier, except the requirement that they not have opposite charge
          const double vdw_min_distance
          (
            itr_atoms->GetCharge() != itr_atoms_b->GetCharge() && itr_atoms->GetCharge() && itr_atoms_b->GetCharge()
            ? util::GetUndefined< double>()
            : *itr_atoms_cov_radius_b + *itr_atoms_cov_radius
          );

          // skip atom types for which the vdw distance and the bond length is undefined
          if( !util::IsDefined( vdw_min_distance) && !util::IsDefined( bond_length))
          {
            continue;
          }

          // compute the actual distance between the two atoms
          const double distance( linal::Distance( itr_atoms->GetPosition(), itr_atoms_b->GetPosition()));

          double max_bond_distance_ratio( 1.25);
          // period 3 elements can have longer bonds due to d-orbital bonding; so allow a larger tolerance for them
          if( itr_atoms->GetElementType()->GetPeriod() > 2)
          {
            max_bond_distance_ratio += 0.125;
          }

          if( itr_atoms_b->GetElementType()->GetPeriod() > 2)
          {
            max_bond_distance_ratio += 0.125;
          }

          storage::Vector< BondConformational>::const_iterator itr_bond( itr_atoms->FindBondTo( *itr_atoms_b));
          const bool are_bonded( itr_bond != itr_atoms->GetBonds().End());
          if( are_bonded)
          {
            const double expected_bond_length
            (
              BondLengths::GetBondLength
              (
                itr_atoms->GetAtomType(),
                itr_bond->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromatic),
                itr_atoms_b->GetAtomType()
              )
            );
            if( util::IsDefined( expected_bond_length))
            {
              bond_length = expected_bond_length;
            }
          }

          // check whether the atoms are bonded
          if( are_bonded)
          {
            // bonded atom, check for bond lengths less than half the normal values
            if( util::IsDefined( bond_length))
            {
              if( distance < 0.75 * bond_length)
              {
                BCL_MessageVrb
                (
                  "Bad Geometry! Too short a bond between atoms of type: " + itr_atoms->GetAtomType().GetName()
                  + " and " + itr_atoms_b->GetAtomType().GetName() + ", distance: "
                  + util::Format()( distance) + " expected bond length: " + util::Format()( bond_length)
                );
                return true;
              }
              else if( distance > max_bond_distance_ratio * bond_length)
              {
                BCL_MessageVrb
                (
                  "Bad Geometry! Too long a bond between atoms of type: " + itr_atoms->GetAtomType().GetName()
                  + " and " + itr_atoms_b->GetAtomType().GetName() + ", distance: "
                  + util::Format()( distance) + " expected bond length: " + util::Format()( bond_length)
                );
                return true;
              }
            }
          }
          // if the distance is <= the expected bond length
          else if( util::IsDefined( bond_length) && distance <= bond_length)
          {
            BCL_MessageVrb
            (
              "Bad Geometry! Non-bonded atom (type: " + itr_atoms->GetAtomType().GetName()
              + ") within bonding distance of another atom (type: " + itr_atoms_b->GetAtomType().GetName() + "), distance: "
              + util::Format()( distance) + " expected bond length: " + util::Format()( bond_length)
            );
            return true;
          }
          else if( util::IsDefined( vdw_min_distance) && distance < vdw_min_distance)
          {
            // check for
            BCL_MessageVrb
            (
              "Bad Geometry! Non-bonded atom (type: " + itr_atoms->GetAtomType().GetName()
              + ") inside van-der waals sphere of another atom (type: " + itr_atoms_b->GetAtomType().GetName() + "), distance: "
              + util::Format()( distance) + " vdw distance: " + util::Format()( vdw_min_distance)
            );
            return true;
          }
        }

        // for the 3D check, skip atoms with a defined z coordinate
        if( !math::EqualWithinTolerance( itr_atoms->GetPosition().Z(), 0))
        {
          continue;
        }

        const storage::Vector< BondConformational> &connected_atoms( itr_atoms->GetBonds());

        // atoms with 0-2 bonds could have all neighboring atoms with a 0 z-coordinate
        if( connected_atoms.GetSize() < 3)
        {
          continue;
        }

        // atoms with 3 bonds could have all neighboring atoms with a 0 z-coordinate if the geometry is not SP3
        if
        (
          connected_atoms.GetSize() == 3
          && itr_atoms->GetAtomType()->GetHybridOrbitalType() != GetHybridOrbitalTypes().e_SP3
        )
        {
          continue;
        }

        bool had_3d_neighbor( false);
        // check that at least one neighbor has a non-zero z-coordinate
        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_connected( connected_atoms.Begin()),
            itr_connected_end( connected_atoms.End());
          itr_connected != itr_connected_end;
          ++itr_connected
        )
        {
          if( !math::EqualWithinTolerance( itr_connected->GetTargetAtom().GetPosition().Z(), 0))
          {
            had_3d_neighbor = true;
            break;
          }
        }
        if( !had_3d_neighbor)
        {
          const std::string &name( itr_atoms->GetAtomType().GetName());
          if( connected_atoms.GetSize() == 3)
          {
            BCL_MessageVrb
            (
              "Bad Geometry! SP3 atom with (type: " + name + "); all neighbors have 0 z-coordinate"
            );
          }
          else
          {
            BCL_MessageVrb
            (
              "Bad Geometry! Atom with 4 bonds (type: " + name + "); all neighbors have 0 z-coordinate"
            );
          }
          return true;
        }
      }

      //const double score( AtomClashScore( true)( *this));
      //if( score > 0.02)
      //{
      //  BCL_MessageVrb( "Molecule would have passed except for clash score: " + util::Format()( score));
      //  AtomClashScore( true, true)( *this);
      //  return true;
      //}

      return false;
    }

    //! @brief checks for non-gasteiger atom types
    //! @return flag whether any non-gasteiger types are present
    bool ConformationInterface::HasNonGasteigerAtomTypes() const
    {
      return ConformationInterface::HasNonGasteigerAtomTypes( GetAtomsIterator());
    }

    //! @brief checks for chiral centers
    //! @return flag whether any chiral centers are present
    bool ConformationInterface::HasChiralCenters() const
    {
      for( auto itr( GetAtomsIterator()); itr.NotAtEnd(); ++itr)
      {
        if( itr->GetChirality() != e_NonChiral)
        {
          return true;
        }
      }
      return false;
    }

    //! @brief checks for double bond isometry
    //! @return double bond isometry
    bool ConformationInterface::HasIsometry() const
    {
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr( GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        const storage::Vector< BondConformational> &connected_atoms( itr->GetBonds());

        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_connected( connected_atoms.Begin()),
            itr_connected_end( connected_atoms.End());
          itr_connected != itr_connected_end;
          ++itr_connected
        )
        {
          if( itr_connected->GetBondType()->GetIsometry() != e_NonIsometric)
          {
            return true;
          }
        }
      }
      return false;
    }

    //! @brief checks any arbitrary generic iterator over AtomConformationInterfcaes for non-gasteiger atom types
    //! @param ATOMS a generic iterator over the atoms
    //! @return flag whether any non-gasteiger types are present
    bool ConformationInterface::HasNonGasteigerAtomTypes( iterate::Generic< const AtomConformationalInterface> ATOMS)
    {
      for( ; ATOMS.NotAtEnd(); ++ATOMS)
      {
        if( !ATOMS->GetAtomType()->IsGasteigerAtomType())
        {
          return true;
        }
      }
      return false;
    }

    //! @brief get the center of the molecule, as defined by the positions
    linal::Vector3D ConformationInterface::GetCenter() const
    {
      math::RunningAverage< linal::Vector3D> averager;

      // iterate through the list of AtomsWithPosition and average their positions
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        if( itr_atoms->GetPosition().IsDefined())
        {
          averager += itr_atoms->GetPosition();
        }
      }

      return averager.GetAverage();
    }

    //! @brief return a SiPtrVector to linal::Vector3D containing atomic coordinates
    util::SiPtrVector< const linal::Vector3D> ConformationInterface::GetAtomCoordinates() const
    {
      util::SiPtrVector< const linal::Vector3D> positions( GetNumberAtoms());
      util::SiPtrVector< const linal::Vector3D>::iterator position_itr( positions.Begin());

      // iterate through the list of AtomsWithPosition and record their positions
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms, ++position_itr
      )
      {
        *position_itr = &itr_atoms->GetPosition();
      }

      return positions;
    }

    //! @brief return a SiPtrVector to linal::Vector3D containing atomic coordinates
    util::SiPtrVector< const linal::Vector3D> ConformationInterface::GetHeavyAtomCoordinates() const
    {
      util::SiPtrVector< const linal::Vector3D> positions;
      positions.AllocateMemory( GetNumberAtoms());

      // iterate through the list of AtomsWithPosition and record their positions
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        if( itr_atoms->GetElementType() != GetElementTypes().e_Hydrogen)
        {
          positions.PushBack( itr_atoms->GetPosition());
        }
      }

      return positions;
    }

    //! @brief calculate the bond angles
    //! @return a map from atom conformation interface in ascending order of atom position
    storage::Vector< double>
      ConformationInterface::GetBondAngles() const
    {
      // begin a vector of bond angles
      storage::Vector< double> bond_angles;

      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        const storage::Vector< BondConformational> &connected_atoms( itr_atoms->GetBonds());

        // get the number of neighbors this atom has
        const size_t number_neighbors( connected_atoms.GetSize());

        // 0-1 neighbors -> no bond angles
        if( number_neighbors < size_t( 2))
        {
          continue;
        }

        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_connected( connected_atoms.Begin()),
            itr_connected_end( connected_atoms.End());
          itr_connected != itr_connected_end;
          ++itr_connected
        )
        {
          // walk through the list of connected neighbor indices with two iterators,
          // thus if the neighbors are 1, 2, 3 for atom 0, the angles will be calculated in the following order:
          // 1 - 0 - 2
          // 1 - 0 - 3
          // 2 - 0 - 3
          for
          (
            storage::Vector< BondConformational>::const_iterator itr_connected_smaller( connected_atoms.Begin());
            itr_connected_smaller != itr_connected;
            ++itr_connected_smaller
          )
          {
            // add the bond angle between these atoms
            const double angle
            (
              linal::ProjAngle
              (
                itr_atoms->GetPosition(),
                itr_connected_smaller->GetTargetAtom().GetPosition(),
                itr_connected->GetTargetAtom().GetPosition()
              )
            );

            // add the bond angle to the vector
            bond_angles.PushBack( angle);
          }
        }
      }

      return bond_angles;
    }

    //! @brief calculate all the dihedral angles
    //! @return a vector containing dihedral angles, in ascending order of atom positions
    storage::Map< storage::VectorND< 5, size_t>, storage::Vector< double> > ConformationInterface::GetBondAnglesByType() const
    {
      // begin a vector of bond angles
      storage::Map< storage::VectorND< 5, size_t>, storage::Vector< double> > bond_angles;

      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( GetAtomsIterator());
        itr_atoms.NotAtEnd();
        ++itr_atoms
      )
      {
        const storage::Vector< BondConformational> &connected_atoms( itr_atoms->GetBonds());

        // get the number of neighbors this atom has
        const size_t number_neighbors( connected_atoms.GetSize());

        // 0-1 neighbors -> no bond angles
        if( number_neighbors < size_t( 2))
        {
          continue;
        }

        // allocate enough memory to hold the resulting bond angles
        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_connected( connected_atoms.Begin()),
            itr_connected_end( connected_atoms.End());
          itr_connected != itr_connected_end;
          ++itr_connected
        )
        {
          // compute the angle A-B-C, where
          // B is the position of the central atom
          // C is an atom that B is connected to
          // A is an atom also connected to B, but with higher index than C
          for
          (
            storage::Vector< BondConformational>::const_iterator itr_connected_smaller( connected_atoms.Begin());
            itr_connected_smaller != itr_connected;
            ++itr_connected_smaller
          )
          {
            // add the bond angle between these atoms
            const double angle
            (
              linal::ProjAngle
              (
                itr_atoms->GetPosition(),
                itr_connected_smaller->GetTargetAtom().GetPosition(),
                itr_connected->GetTargetAtom().GetPosition()
              )
            );

            // compose the vectornd5 containing AT1, BT12, AT2, BT23, AT3(AT = atom type, BT = bond type)
            storage::VectorND< 5, size_t> bond_type;
            bond_type( 0) = itr_connected->GetTargetAtom().GetAtomType().GetIndex();
            bond_type( 1) = itr_connected->GetBondType()->WithoutIsometry()->WithoutAromaticOrder();
            bond_type( 2) = itr_atoms->GetAtomType().GetIndex();
            bond_type( 3) = itr_connected_smaller->GetBondType()->WithoutIsometry()->WithoutAromaticOrder();
            bond_type( 4) = itr_connected_smaller->GetTargetAtom().GetAtomType().GetIndex();

            // order the bond type
            bool need_to_reverse( false);
            for( size_t i_forward( 0), i_reverse( 4); i_forward < i_reverse; ++i_forward, --i_reverse)
            {
              if( bond_type( i_forward) != bond_type( i_reverse))
              {
                need_to_reverse = bond_type( i_forward) > bond_type( i_reverse);
                break;
              }
            }
            if( need_to_reverse)
            {
              for( size_t i_forward( 0), i_reverse( 4); i_forward < i_reverse; ++i_forward, --i_reverse)
              {
                std::swap( bond_type( i_forward), bond_type( i_reverse));
              }
            }

            bond_angles[ bond_type].PushBack( angle);
          }
        }
      }
      return bond_angles;
    }

    //! @brief calculate all the bond lengths
    //! @return a vector containing bond lengths
    storage::Map< storage::VectorND< 3, size_t>, storage::Vector< double> > ConformationInterface::GetBondLengthsByType() const
    {
      storage::Map< storage::VectorND< 3, size_t>, storage::Vector< double> > bond_types_to_length;

      // compute the bond lengths
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr( GetAtomsIterator());
        itr.NotAtEnd();
        ++itr
      )
      {
        const storage::Vector< BondConformational> &connected_atoms( itr->GetBonds());

        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_connected( connected_atoms.Begin()),
            itr_connected_end( connected_atoms.End());
          itr_connected != itr_connected_end;
          ++itr_connected
        )
        {
          if( &*itr > &itr_connected->GetTargetAtom())
          {
            continue;
          }
          storage::VectorND< 3, size_t> bond_id_info;
          bond_id_info( 0) = itr->GetAtomType().GetIndex();
          bond_id_info( 1) = itr_connected->GetBondType()->WithoutAromaticOrder()->WithoutIsometry();
          bond_id_info( 2) = itr_connected->GetTargetAtom().GetAtomType().GetIndex();

          if( bond_id_info( 0) > bond_id_info( 2))
          {
            std::swap( bond_id_info( 0), bond_id_info( 2));
          }

          const double bond_length( linal::Distance( itr->GetPosition(), itr_connected->GetTargetAtom().GetPosition()));

          bond_types_to_length[ bond_id_info].PushBack( bond_length);
        }
      }
      return bond_types_to_length;
    }

    //! @brief calculate all the dihedral angles
    //! @return a vector containing dihedral angles, in ascending order of atom positions
    storage::Vector< double> ConformationInterface::GetDihedralAngles() const
    {
      // begin a vector of dihedral angles
      storage::Vector< double> dihedral_angles;

      // compute the angle A-B-C-D, where
      // B is the position of the first atom in the bond
      // C is the position of the second atom in the bond
      // A is the positions of the set of atoms that the first atom in the bond is connected to, except C
      // D is the positions of the set of atoms that the second atom in the bond is connected to, except B

      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atom_b( GetAtomsIterator());
        itr_atom_b.NotAtEnd();
        ++itr_atom_b
      )
      {
        const AtomConformationalInterface &atom_b( *itr_atom_b);
        const storage::Vector< BondConformational> &connected_atoms_b( atom_b.GetBonds());
        if( connected_atoms_b.GetSize() < 2)
        {
          continue;
        }
        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_atom_c( connected_atoms_b.Begin()), itr_atom_c_end( connected_atoms_b.End());
          itr_atom_c != itr_atom_c_end;
          ++itr_atom_c
        )
        {
          const AtomConformationalInterface &atom_c( itr_atom_c->GetTargetAtom());
          // ensure that atom c is before atom b; this orders the dihedral angles so that we don't see each
          // dihedral angle twice
          if( &atom_c > &atom_b)
          {
            continue;
          }

          const storage::Vector< BondConformational> &connected_atoms_c( atom_c.GetBonds());
          if( connected_atoms_c.GetSize() == size_t( 1))
          {
            continue;
          }
          for
          (
            storage::Vector< BondConformational>::const_iterator itr_atom_a( connected_atoms_b.Begin());
            itr_atom_a != itr_atom_c_end;
            ++itr_atom_a
          )
          {
            if( itr_atom_a == itr_atom_c)
            {
              continue;
            }
            for
            (
              storage::Vector< BondConformational>::const_iterator
              itr_atom_d( connected_atoms_c.Begin()), itr_atom_d_end( connected_atoms_c.End());
              itr_atom_d != itr_atom_d_end;
              ++itr_atom_d
            )
            {
              if( &atom_b == &itr_atom_d->GetTargetAtom())
              {
                continue;
              }

              dihedral_angles.PushBack
              (
                linal::Dihedral
                (
                  itr_atom_a->GetTargetAtom().GetPosition(),
                  atom_b.GetPosition(),
                  atom_c.GetPosition(),
                  itr_atom_d->GetTargetAtom().GetPosition()
                )
              );
            }
          }
        }
      }
      return dihedral_angles;
    }

    //! @brief calculate all the dihedral angles
    //! @return a vector containing dihedral angles, in ascending order of atom position
    storage::Map< storage::VectorND< 7, size_t>, storage::Vector< double> >
      ConformationInterface::GetDihedralAnglesByType() const
    {
      // begin a vector of dihedral angles
      storage::Map< storage::VectorND< 7, size_t>, storage::Vector< double> > dihedral_angles;

      // compute the angle A-B-C-D, where
      // B is the position of the first atom in the bond
      // C is the position of the second atom in the bond
      // A is the positions of the set of atoms that the first atom in the bond is connected to, except C
      // D is the positions of the set of atoms that the second atom in the bond is connected to, except B

      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atom_b( GetAtomsIterator());
        itr_atom_b.NotAtEnd();
        ++itr_atom_b
      )
      {
        const AtomConformationalInterface &atom_b( *itr_atom_b);
        const storage::Vector< BondConformational> &connected_atoms_b( atom_b.GetBonds());
        if( connected_atoms_b.GetSize() < 2)
        {
          continue;
        }
        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_atom_c( connected_atoms_b.Begin()), itr_atom_c_end( connected_atoms_b.End());
          itr_atom_c != itr_atom_c_end;
          ++itr_atom_c
        )
        {
          const AtomConformationalInterface &atom_c( itr_atom_c->GetTargetAtom());
          // ensure that atom c is before atom b; this orders the dihedral angles so that we don't see each
          // dihedral angle twice
          if( &atom_c > &atom_b)
          {
            continue;
          }

          const storage::Vector< BondConformational> &connected_atoms_c( atom_c.GetBonds());
          if( connected_atoms_c.GetSize() < 2)
          {
            continue;
          }
          for
          (
            storage::Vector< BondConformational>::const_iterator itr_atom_a( connected_atoms_b.Begin());
            itr_atom_a != itr_atom_c_end;
            ++itr_atom_a
          )
          {
            if( itr_atom_a == itr_atom_c)
            {
              continue;
            }
            for
            (
              storage::Vector< BondConformational>::const_iterator
              itr_atom_d( connected_atoms_c.Begin()), itr_atom_d_end( connected_atoms_c.End());
              itr_atom_d != itr_atom_d_end;
              ++itr_atom_d
            )
            {
              if( &atom_b == &itr_atom_d->GetTargetAtom())
              {
                continue;
              }

              const double dihedral_angle
              (
                math::Absolute
                (
                  linal::Dihedral
                  (
                    itr_atom_a->GetTargetAtom().GetPosition(),
                    itr_atom_b->GetPosition(),
                    itr_atom_c->GetTargetAtom().GetPosition(),
                    itr_atom_d->GetTargetAtom().GetPosition()
                  )
                )
              );

              storage::VectorND< 7, size_t> dihedral_type;
              dihedral_type( 0) = itr_atom_a->GetTargetAtom().GetAtomType().GetIndex();
              dihedral_type( 1) = itr_atom_a->GetBondType()->WithoutIsometry()->WithoutAromaticOrder();
              dihedral_type( 2) = itr_atom_b->GetAtomType().GetIndex();
              dihedral_type( 3) = itr_atom_c->GetBondType()->WithoutIsometry()->WithoutAromaticOrder();
              dihedral_type( 4) = itr_atom_c->GetTargetAtom().GetAtomType().GetIndex();
              dihedral_type( 5) = itr_atom_d->GetBondType()->WithoutIsometry()->WithoutAromaticOrder();
              dihedral_type( 6) = itr_atom_d->GetTargetAtom().GetAtomType().GetIndex();

              // order the dihedral type
              bool need_to_reverse( false);
              for( size_t i_forward( 0), i_reverse( 6); i_forward < i_reverse; ++i_forward, --i_reverse)
              {
                if( dihedral_type( i_forward) != dihedral_type( i_reverse))
                {
                  need_to_reverse = dihedral_type( i_forward) > dihedral_type( i_reverse);
                  break;
                }
              }
              if( need_to_reverse)
              {
                for( size_t i_forward( 0), i_reverse( 6); i_forward < i_reverse; ++i_forward, --i_reverse)
                {
                  std::swap( dihedral_type( i_forward), dihedral_type( i_reverse));
                }
              }
              dihedral_angles[ dihedral_type].PushBack( dihedral_angle);
            }
          }
        }
      }
      return dihedral_angles;
    }

    //! @brief compute deviations from planarity for all amide bonds
    //! @return a vector containing the indices of the N and C bonded amide atoms and
    //! the corresponding deviation from nonplanarity
    //! @return a vector of triplets where the vector is indexed by the amide bond and the triplet contains
    //! atom indices corresponding to the lower and higher amide bond atom indices, respectively, and
    //! the magnitude of the deviation
    storage::Vector< storage::Triplet< size_t, size_t, double> > ConformationInterface::GetAmideBondNonPlanarity() const
        {
          storage::Vector< storage::Triplet< size_t, size_t, double> > deviations;
      const AtomTypes &types( GetAtomTypes());
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atom( GetAtomsIterator());
        itr_atom.NotAtEnd();
        ++itr_atom
      )
      {
        const AtomConformationalInterface &atom_a( *itr_atom);
        if( !( atom_a.GetAtomType() == types.O_Tr2Tr2TrPi || atom_a.GetAtomType() == types.S_Tr2Tr2TrPi))
        {
          continue;
        }
        for
        (
          storage::Vector< BondConformational>::const_iterator
            itr_bonds( atom_a.GetBonds().Begin()), itr_bonds_end( atom_a.GetBonds().End());
          itr_bonds != itr_bonds_end;
          ++itr_bonds
        )
        {
          const AtomConformationalInterface &target_atom( itr_bonds->GetTargetAtom());
          if
          (
            target_atom.GetAtomType()->GetHybridOrbitalType() == GetHybridOrbitalTypes().e_SP3 ||
            target_atom.GetAtomType()->GetNumberBonds() != 3
          )
          {
            continue;
          }

          for
          (
            storage::Vector< BondConformational>::const_iterator
              itr_bonds_target( target_atom.GetBonds().Begin()), itr_bonds_target_end( target_atom.GetBonds().End());
            itr_bonds_target != itr_bonds_target_end;
            ++itr_bonds_target
          )
          {
            if( itr_bonds_target->GetBondType()->GetConjugation() != ConstitutionalBondTypeData::e_Amide)
            {
              continue;
            }
            const AtomConformationalInterface &target_atom_b( itr_bonds_target->GetTargetAtom());

            bool ignore_bond( false);
            storage::Vector< linal::Vector3D> target_atoms_c_coords;
            for
            (
              storage::Vector< BondConformational>::const_iterator
                itr_bonds_target_b( target_atom_b.GetBonds().Begin()), itr_bonds_target_b_end( target_atom_b.GetBonds().End());
              itr_bonds_target_b != itr_bonds_target_b_end;
              ++itr_bonds_target_b
            )
            {
              const AtomConformationalInterface &target_atom_c( itr_bonds_target_b->GetTargetAtom());
              if
              (
                target_atom_c.GetAtomType()->GetHybridOrbitalType() == GetHybridOrbitalTypes().e_SP3 ||
                target_atom_c.GetAtomType()->GetNumberBonds() == 4
              )
              {
                continue;
              }
              if( itr_bonds_target_b->GetBondType()->IsBondInRing())
              {
                continue;
              }
              if( &target_atom_c == &target_atom)
              {
                continue;

              }
              double dihedral_angle
              (
                math::Absolute
                (
                  math::Angle::Degree
                  (
                    linal::Dihedral
                    (
                      atom_a.GetPosition(),
                      target_atom.GetPosition(),
                      target_atom_b.GetPosition(),
                      target_atom_c.GetPosition()
                    )
                  )
                )
              );

              // track indices of atoms in amide bond
              size_t c_index( this->GetAtomIndex(target_atom));
              size_t n_index( this->GetAtomIndex(target_atom_b));

              // place low index first
              c_index < n_index ?
                deviations.PushBack( storage::Triplet< size_t, size_t, double>( c_index, n_index, std::min( dihedral_angle, 180.0 - dihedral_angle))) :
                deviations.PushBack( storage::Triplet< size_t, size_t, double>( n_index, c_index, std::min( dihedral_angle, 180.0 - dihedral_angle)));
            }
          }
        }
      }
      return deviations;
    }

    //! @brief Sum the deviation of amide bond planarity over the molecule
    //! @param AMIDE_DEVIATIONS from GetAmideBondNonPlanarity()
    //! @param TOLERANCES (in degrees) acceptable deformations of the amide
    //! @param PENALTIES (in BCL::Conf units) for amide deviations exceeding the tolerance;
    //! must be same size as TOLERANCES
    //! @return the per-amide bond penalty
    storage::Vector< double> ConformationInterface::GetPerAmideBondNonPlanarityPenalty
    (
      storage::Vector< storage::Triplet< size_t, size_t, double> > &AMIDE_DEVIATIONS,
      storage::Vector< double> &TOLERANCES,
      storage::Vector< double> &PENALTIES
    ) const
    {
      // require same sizes of tolerances and penalties
      BCL_Assert
      (
        TOLERANCES.GetSize() == PENALTIES.GetSize(),
        "Amide bond planarity tolerances and associated penalties are different sized vectors!"
        );

      // get the individual deviations from planarity
      storage::Vector< double> amide_penalties( AMIDE_DEVIATIONS.GetSize(), double( 0.0));

      // estimate penalty for each deviation
      for
      (
          size_t i( 0), n_deviations( AMIDE_DEVIATIONS.GetSize());
          i < n_deviations;
          ++i
      )
      {
        for
        (
            size_t p( 0), n_penalties( PENALTIES.GetSize());
            p < n_penalties;
            ++p
        )
        {
          // raw deviation from tolerance to which penalty will be applied
          double deviation_from_tolerance( std::max( 0.0, AMIDE_DEVIATIONS( i).Third() - TOLERANCES( p)));
          amide_penalties( i) += deviation_from_tolerance * PENALTIES( p);
        }
      }
      return amide_penalties;
    }

    //! @brief Sum the deviation of amide bond planarity over the molecule
    //! @param TOLERANCE tolerance (in degrees) for non-planarity, amide bonds with deviation smaller than tolerance
    //!        will contribute 0 to the sum
    //! @return deviation of amide bond planarity summed over the molecule
    double ConformationInterface::GetTotalAmideBondNonPlanarity( const double &TOLERANCE) const
    {
      // get the individual deviations from planarity
      storage::Vector< storage::Triplet< size_t, size_t, double>> amide_deviations( GetAmideBondNonPlanarity());

      // sum them up
      double total_deviation( 0.0);
      for( size_t i( 0), sz( amide_deviations.GetSize()); i < sz; ++i)
      {
        total_deviation += std::max( 0.0, amide_deviations( i).Third() - TOLERANCE);
      }

      // amide bond out-of-tolerance deviation across whole molecule (in degrees)
      return total_deviation;
    }

    //! @brief check to see if amide bonds are planer
    //! @param TOLERANCE tolerance (in degrees) for non-planarity
    //! @return true if bonds are planer else return false
    bool ConformationInterface::AreAmideBondsPlaner( const double &TOLERANCE) const
    {
      return GetTotalAmideBondNonPlanarity( TOLERANCE) <= 0.0;
    }

    //! @brief check to see if aromatic rings are planer
    //! @return true if aromatic rings are planar, if not return false
    bool ConformationInterface::AreAromaticRingsPlaner() const
    {
      FragmentEnsemble rings( FragmentSplitUnbridgedRings( true, 15)( *this));

      static const std::string s_bond_str( " -=#?");

      for
      (
        storage::List< FragmentComplete>::const_iterator itr_rings( rings.Begin()), itr_rings_end( rings.End());
          itr_rings != itr_rings_end;
        ++itr_rings
      )
      {
        if( itr_rings->GetSize() == size_t( 3))
        {
          continue;
        }

        #ifdef BCL_ShowAromaticRingPlanarityDetail
        std::string smile_str;
        if( itr_rings->GetNumberAtoms() == itr_rings->GetNumberBonds())
        {
          FragmentComplete com( *itr_rings);
          com.Canonicalize();
          std::ostringstream smile_out;
          util::SiPtr< const AtomConformationalInterface> first_atom( *com.GetAtomsIterator());
          util::SiPtr< const AtomConformationalInterface> this_atom( first_atom);
          ConfigurationalBondType cbt( this_atom->GetBonds()( 0).GetBondType());
          util::SiPtr< const AtomConformationalInterface> next_atom( this_atom->GetBonds()( 0).GetTargetAtom());
          do
          {
            smile_out << this_atom->GetAtomType().GetName()
                      << s_bond_str[ std::min( size_t( cbt->GetNumberOfElectrons() / 2), size_t( 4))];
            util::SiPtr< const AtomConformationalInterface> n_next_atom( next_atom);
            size_t next_atom_bond_index( &( next_atom->GetBonds()( 0).GetTargetAtom()) == &*this_atom ? 1 : 0);
            cbt = next_atom->GetBonds()( next_atom_bond_index).GetBondType();
            next_atom = next_atom->GetBonds()( next_atom_bond_index).GetTargetAtom();
            this_atom = n_next_atom;
          } while( this_atom != first_atom);
          smile_str = smile_out.str();
        }
        else
        {
          FragmentComplete com( *itr_rings);
          com.Canonicalize();
          smile_str += com.GetAtomTypesString() + " " + com.GetBondTypesString();
        }
        #endif

        // as a first pass, check the dihedrals. All dihedrals should at 0/360 degrees (symmetric) for small rings, and may
        // be 180 degrees for fused rings
        auto dihedrals( itr_rings->GetDihedralAngles());
        ConformationComparisonByDihedralBins ccdb;
        for
        (
          auto itr_dihedrals( dihedrals.Begin()), itr_dihedrals_end( dihedrals.End());
          itr_dihedrals != itr_dihedrals_end;
          ++itr_dihedrals
        )
        {
          int bin( ccdb.DetermineDihedralKey( *itr_dihedrals, false, false));
          if( bin != 12 && bin != 6)
          {
            return false;
          }
        }

        util::SiPtrVector< const linal::Vector3D> coordinates( itr_rings->GetAtomCoordinates());
        math::RunningAverageSD< linal::Vector3D> centroid_ave;
        size_t nrow( 0), n_atoms( coordinates.GetSize());
        linal::Matrix< double> coord_xyz( coordinates.GetSize(), 3);

        for
        (
          util::SiPtrVector< const linal::Vector3D>::const_iterator itr( coordinates.Begin()), itr_end( coordinates.End());
          itr != itr_end;
          ++itr, ++nrow
        )
        {
          // get coordinates from this
          coord_xyz( nrow, 0) = ( *itr)->X();
          coord_xyz( nrow, 1) = ( *itr)->Y();
          coord_xyz( nrow, 2) = ( *itr)->Z();

          // add it to the centroid
          centroid_ave += **itr;
        }
        for( nrow = 0; nrow < n_atoms; ++nrow)
        {
          // get coordinates from this
          coord_xyz( nrow, 0) -= centroid_ave.GetAverage()( 0);
          coord_xyz( nrow, 1) -= centroid_ave.GetAverage()( 1);
          coord_xyz( nrow, 2) -= centroid_ave.GetAverage()( 2);
        }
        auto soln_chisq( math::LinearLeastSquares::SolutionAndChiSquared( coord_xyz, linal::Vector< double>( n_atoms, 1.0)));
        linal::Matrix3x3< double> eigen_v;
        linal::Vector< double> eigen_val( 3);
        linal::Matrix3x3< double> mtm( linal::MatrixTransposeTimesMatrix( coord_xyz));
        mtm.EigenVectorsSymmetric( eigen_v, eigen_val);
        eigen_v.Transpose();

        // calculate the sum of the residuals.
        // a residual is the experimental value minus the predicted value squared
        // the predicted value is the row of the matrix times the solutions
        double mx1( 0.0), mx2( 0.0), mx3( 0.0);
        for( size_t index( 0); index < n_atoms; ++index)
        {
          // calculate the value given by the linear-least squares fit to the data
          double linear_estimate1( std::inner_product( eigen_v[ 0], eigen_v[ 0] + 3, coord_xyz[ index], 0.0));
          double linear_estimate2( std::inner_product( eigen_v[ 1], eigen_v[ 1] + 3, coord_xyz[ index], 0.0));
          double linear_estimate3( std::inner_product( eigen_v[ 2], eigen_v[ 2] + 3, coord_xyz[ index], 0.0));

          // add the difference between the estimate and the actual value, squared, to the residual sum
          mx1 = std::max( mx1, math::Absolute( linear_estimate1));
          mx2 = std::max( mx2, math::Absolute( linear_estimate2));
          mx3 = std::max( mx3, math::Absolute( linear_estimate3));
        }
        double offset( 0.05 * ( itr_rings->GetNumberBonds() - itr_rings->GetNumberAtoms()));
        if( std::min( mx1, std::min( mx2, mx3)) >= ( 0.1 + offset))
        {
          #ifdef BCL_ShowAromaticRingPlanarityDetail
          BCL_MessageStd( "Ring with types is non-aromatic: " + smile_str);
          #endif
          return false;
        }
        #ifdef BCL_ShowAromaticRingPlanarityDetail
        else
        {
          BCL_MessageStd( "Ring with types is aromatic: " + smile_str);
        }
        #endif
      }
      return true;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief move the molecule
    //! @param TRANSLATION amount to change x,y,z coordinates
    void ConformationInterface::Translate( const linal::Vector3D &TRANSLATION)
    {
      // add TRANSLATION to all the positions in the molecule
      for( iterate::Generic< AtomConformationalInterface> itr( GetAtomsIteratorNonConst()); itr.NotAtEnd(); ++itr)
      {
        itr->SetPosition( itr->GetPosition() + TRANSLATION);
      }
    }

    //! @brief transform the object by a given TransformationMatrix3D
    //! @param TRANSFORMATION_MATRIX_3D TransformationMatrix3D to be applied
    void ConformationInterface::Transform( const math::TransformationMatrix3D &TRANSFORMATION_MATRIX_3D)
    {
      // transform each position using TRANSFORMATION_MATRIX_3D
      for( iterate::Generic< AtomConformationalInterface> itr( GetAtomsIteratorNonConst()); itr.NotAtEnd(); ++itr)
      {
        // get the current position
        linal::Vector3D current_position( itr->GetPosition());
        // transform the current position using the matrix
        current_position.Transform( TRANSFORMATION_MATRIX_3D);
        // set the new positions
        itr->SetPosition( current_position);
      }
    }

    //! @brief rotate the object by a given RotationMatrix3D
    //! @param ROTATION_MATRIX_3D RotationMatrix3D to be applied
    void ConformationInterface::Rotate( const math::RotationMatrix3D &ROTATION_MATRIX_3D)
    {
      // rotate each position using ROTATION_MATRIX_3D
      for( iterate::Generic< AtomConformationalInterface> itr( GetAtomsIteratorNonConst()); itr.NotAtEnd(); ++itr)
      {
        // get the current position
        linal::Vector3D current_position( itr->GetPosition());
        // rotate the current position using the matrix
        current_position.Rotate( ROTATION_MATRIX_3D);
        // overwrite the existing atoms with new atoms containing the new positions
        itr->SetPosition( current_position);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &ConformationInterface::WriteMDL( std::ostream &OSTREAM) const
    {
      return
        sdf::MdlHandler::WriteToSDF
        (
          OSTREAM,
          GetName(),
          GetAtomInfo(),
          GetBondInfo(),
          GetStoredProperties().GetMDLProperties()
        );
    }

  } // namespace chemistry
} // namespace bcl
