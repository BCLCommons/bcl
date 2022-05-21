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
#include "chemistry/bcl_chemistry_hydrogens_handler.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_si_ptr_vector.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new HydrogensHandler
    HydrogensHandler *HydrogensHandler::Clone() const
    {
      return new HydrogensHandler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &HydrogensHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief handles hydrogens for atoms complete
    //! @param ATOMS the fragment to add hydrogens to
    //! @param H_PREF hydrogen preference
    //! @param ADD_H whether to add hydrogen to each atom, if H_PREF = e_MdlSaturate
    void HydrogensHandler::HandleHydrogenPref
    (
      AtomVector< AtomComplete> &ATOMS,
      const sdf::HydrogenHandlingPref &H_PREF,
      const storage::Vector< size_t> &ADD_H
    )
    {
      switch( H_PREF)
      {
        case sdf::e_Remove:
          Remove( ATOMS);
          break;
        case sdf::e_Saturate:
          Saturate( ATOMS);
          break;
        case sdf::e_SaturatePartial:
          SaturatePartial( ATOMS, ADD_H);
          break;
        case sdf::e_Maintain:
        default:
          break;
      }
    }

    //! @brief handles hydrogens for atoms complete
    //! @param ATOMS the fragment to add hydrogens to
    void HydrogensHandler::Saturate( AtomVector< AtomComplete> &ATOMS)
    {
      SaturatePartial( ATOMS, storage::Vector< size_t>( ATOMS.GetSize(), size_t( 1)));
    }

    //! @brief removes hydrogens for atoms complete
    //! @param ATOMS the fragment from which hydrogens will be removed
    void HydrogensHandler::Remove( AtomVector< AtomComplete> &ATOMS)
    {
      storage::Vector< size_t> non_h_indices;
      const size_t original_number_atoms( ATOMS.GetSize());
      non_h_indices.AllocateMemory( original_number_atoms);
      for( size_t atom_id( 0); atom_id < original_number_atoms; ++atom_id)
      {
        if( ATOMS( atom_id).GetElementType() != GetElementTypes().e_Hydrogen)
        {
          non_h_indices.PushBack( atom_id);
        }
      }
      ATOMS.Reorder( non_h_indices);
    }

    //! @brief handles hydrogens for atoms complete
    //! @param ATOMS the fragment to add hydrogens to
    void HydrogensHandler::SaturatePartial
    (
      AtomVector< AtomComplete> &ATOMS,
      const storage::Vector< size_t> &ADD_H
    )
    {
      // create vector for the new hydrogen ATOMS and their bonds
      storage::Vector< storage::Pair< size_t, linal::Vector3D> > atom_id_to_hydrogen_position;

      // record # atoms in ATOMS
      const size_t number_atoms( ATOMS.GetSize());

      size_t atom_id( 0);
      // iterate over the valence ATOMS
      for
      (
        AtomVector< AtomComplete>::const_iterator itr_atom( ATOMS.Begin()), itr_atom_end( ATOMS.End());
        itr_atom != itr_atom_end;
        ++itr_atom, ++atom_id
      )
      {
        // ignore undefined atom types
        if( !itr_atom->GetAtomType()->IsGasteigerAtomType())
        {
          continue;
        }

        if( !ADD_H( atom_id))
        {
          continue;
        }

        // set any unknown coordinates
        SetUndefinedHCoordinates( *itr_atom, ATOMS);

        // get the coordinates for each hydrogen
        storage::Vector< linal::Vector3D> hydrogen_positions( DetermineCoordinates( *itr_atom));

        // add hydrogens until the number of bonds on the atom == the # of bonds for the type
        for
        (
          storage::Vector< linal::Vector3D>::const_iterator
            itr_h( hydrogen_positions.Begin()),
            itr_h_end( hydrogen_positions.End());
          itr_h != itr_h_end;
          ++itr_h
        )
        {
          atom_id_to_hydrogen_position.PushBack( storage::Pair< size_t, linal::Vector3D>( atom_id, *itr_h));
        }
      }

      // create initializers for each of the hydrogens
      storage::Vector< sdf::AtomInfo> hydrogen_component;
      storage::Vector< sdf::BondInfo> bonds_to_h;
      hydrogen_component.AllocateMemory( atom_id_to_hydrogen_position.GetSize());
      bonds_to_h.AllocateMemory( atom_id_to_hydrogen_position.GetSize());
      size_t h_id( 0);
      for
      (
        storage::Vector< storage::Pair< size_t, linal::Vector3D> >::const_iterator
          itr( atom_id_to_hydrogen_position.Begin()), itr_end( atom_id_to_hydrogen_position.End());
        itr != itr_end;
        ++itr, ++h_id
      )
      {
        hydrogen_component.PushBack( sdf::AtomInfo( GetAtomTypes().H_S, e_NonChiral, itr->Second(), false));
        bonds_to_h.PushBack
        (
          sdf::BondInfo( itr->First(), h_id + number_atoms, GetConfigurationalBondTypes().e_NonConjugatedSingleBond)
        );
      }

      // create an atom vector with just the hydrogens
      const AtomVector< AtomComplete> hydrogens( hydrogen_component, storage::Vector< sdf::BondInfo>());

      ATOMS.AddAtomsWithConnectivity( hydrogens, bonds_to_h);
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &HydrogensHandler::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &HydrogensHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

    //! @brief determine the coordinates of all missing hydrogens for a given atom
    //! @param ATOM an AtomComplete that should have both type and position defined
    //! @return a storage vector containing the positions of all missing hydrogens
    storage::Vector< linal::Vector3D> HydrogensHandler::DetermineCoordinates( const AtomComplete &ATOM)
    {
      // check whether the atom has the same # of bonds that the atom type has

      // get a const reference to the bonds of this atom
      const storage::Vector< BondConformational> &bonds( ATOM.GetBonds());
      const size_t nr_known_bonds( bonds.GetSize());
      size_t number_single_bonds( ATOM.GetNumberofValenceBondsWithOrder( 1));

      const size_t nr_missing_bonds( number_single_bonds);

      if( nr_missing_bonds == 0)
      {
        return storage::Vector< linal::Vector3D>();
      }

      storage::Vector< linal::Vector3D> positions;
      positions.AllocateMemory( nr_missing_bonds);

      // store the expected bond length, which is the covalent radius of H plus the covalent radius of ATOM
      const double bond_length
      (
        GetElementTypes().e_Hydrogen->GetProperty( ElementTypeData::e_CovalentRadius) +
        ATOM.GetElementType()->GetProperty( ElementTypeData::e_CovalentRadius)
      );

      // get the hybrid orbital for the atom type
      HybridOrbitalType orbital_type( ATOM.GetAtomType()->GetHybridOrbitalType());

      // if the orbital type is unhybridized (this really shouldn't ever happen), then approximate it by considering it SP3
      if( orbital_type == GetHybridOrbitalTypes().e_Unhybridized)
      {
        orbital_type = GetHybridOrbitalTypes().e_SP3;
      }

      // get the position of this atom
      const linal::Vector3D &position( ATOM.GetPosition());

      if( !position.IsDefined())
      {
        return storage::Vector< linal::Vector3D>( nr_missing_bonds, linal::Vector3D( util::GetUndefined< double>()));
      }

      // 0-1 known bonds, just make nr_missing_bonds, with the right angle for the hybridization
      if( nr_known_bonds == 0)
      {
        positions = GetIdealizedGeometry( orbital_type, bond_length);

        // translate the points such that they are relative to position
        for
        (
          storage::Vector< linal::Vector3D>::iterator itr( positions.Begin()), itr_end( positions.End());
          itr != itr_end;
          ++itr
        )
        {
          *itr += position;
        }
      }
      else if( nr_known_bonds == 1)
      {
        if( orbital_type == GetHybridOrbitalTypes().e_SP)
        {
          positions.PushBack( linal::CoordinatesLinear( position, bonds( 0).GetTargetAtom().GetPosition(), bond_length));
        }
        else if( orbital_type == GetHybridOrbitalTypes().e_SP2)
        {
          const storage::Vector< BondConformational> &bonds_neighbor( bonds( 0).GetTargetAtom().GetBonds());
          storage::Vector< linal::Vector3D> neighboring_positions;

          for
          (
            storage::Vector< BondConformational>::const_iterator
              itr_neighbors( bonds_neighbor.Begin()), itr_neighbors_end( bonds_neighbor.End());
            itr_neighbors != itr_neighbors_end;
            ++itr_neighbors
          )
          {
            if( &itr_neighbors->GetTargetAtom() != &ATOM)
            {
              neighboring_positions.PushBack( itr_neighbors->GetTargetAtom().GetPosition());
              break;
            }
          }

          if( !neighboring_positions.IsEmpty())
          {
            // get the anchor points for bond rotation, i.e. indices of atoms that make the bond of interest
            const linal::Vector3D &start_point( bonds( 0).GetTargetAtom().GetPosition());

            // get the unit vector from the start point's other neighbor to this
            positions.PushBack( position - linal::UnitVector( start_point, neighboring_positions( 0)) * bond_length);
          }
          else
          {
            // coordinate 1
            positions.PushBack
            (
              linal::CoordinatesAngle
              (
                position,
                bonds( 0).GetTargetAtom().GetPosition(),
                linal::Vector3D( bond_length),
                bond_length,
                120.0 / 180.0 * math::g_Pi,
                120.0 / 180.0 * math::g_Pi,
                false,
                linal::Vector3D( 0.0)
              )
            );
          }
          // coordinate 2
          positions.PushBack
          (
            linal::CoordinatesTrigonal
            (
              position,
              bonds( 0).GetTargetAtom().GetPosition(),
              positions.LastElement(),
              bond_length
            )
          );
        }
        else // SP3
        {
          // coordinate 1
          positions.PushBack
          (
            linal::CoordinatesAngle
            (
              position,
              bonds( 0).GetTargetAtom().GetPosition(),
              linal::Vector3D( bond_length),
              bond_length,
              109.5 / 180.0 * math::g_Pi,
              109.5 / 180.0 * math::g_Pi,
              false,
              linal::Vector3D( 0.0)
            )
          );

          // helper coordinates
          linal::Vector3D foot_point
          (
            linal::CoordinatesTrigonal
            (
              position,
              bonds( 0).GetTargetAtom().GetPosition(),
              positions.LastElement(),
              bond_length * std::cos( 54.75 / 180 * math::g_Pi)
            )
          );

          linal::Vector3D offset
          (
            bond_length * std::sin( 54.75 / 180 * math::g_Pi) *
            linal::CrossProduct
            (
              bonds( 0).GetTargetAtom().GetPosition() - position,
              positions.LastElement() - position
            ).Normalize()
          );

          // coordinate 2
          positions.PushBack( foot_point + offset);

          // coordinate 3
          positions.PushBack( foot_point - offset);
        }
      }
      else if( nr_known_bonds == 2) // 2 known bonds
      {
        if( orbital_type == GetHybridOrbitalTypes().e_SP2) // trigonal
        {
          positions.PushBack
          (
            linal::CoordinatesTrigonal
            (
              position,
              bonds( 0).GetTargetAtom().GetPosition(),
              bonds( 1).GetTargetAtom().GetPosition(),
              bond_length
            )
          );
        }
        else // tetrahedral
        {
          // helper coordinates
          linal::Vector3D foot_point
          (
            linal::CoordinatesTrigonal
            (
              position,
              bonds( 0).GetTargetAtom().GetPosition(),
              bonds( 1).GetTargetAtom().GetPosition(),
              bond_length * std::cos( 54.75 / 180 * math::g_Pi)
            )
          );

          linal::Vector3D offset
          (
            bond_length * std::sin( 54.75 / 180 * math::g_Pi) *
            linal::CrossProduct
            (
              bonds( 0).GetTargetAtom().GetPosition() - position,
              bonds( 1).GetTargetAtom().GetPosition() - position
            ).Normalize()
          );

          // coordinate 1
          positions.PushBack( foot_point + offset);

          // coordinate 2
          positions.PushBack( foot_point - offset);
        }
      }
      else if( nr_known_bonds == 3) // tetrahedral geometry, only missing one bond
      {
        positions.PushBack
        (
          linal::CoordinatesTetrahedral
          (
            position,
            bonds( 0).GetTargetAtom().GetPosition(),
            bonds( 1).GetTargetAtom().GetPosition(),
            bonds( 2).GetTargetAtom().GetPosition(),
            bond_length
          )
        );
      }

      positions.Resize( nr_missing_bonds);
      return positions;
    }

    //! @brief get the idealized geometry for a particular hybrid orbital type and bond length
    //! @param ORBITAL_TYPE the type of orbital associated with the geometry
    //! @param BOND_LENGTH the length of the bonds in the geometry
    storage::Vector< linal::Vector3D> HydrogensHandler::GetIdealizedGeometry
    (
      const HybridOrbitalType &ORBITAL_TYPE,
      const double &BOND_LENGTH
    )
    {
      storage::Vector< linal::Vector3D> positions;

      if( ORBITAL_TYPE == GetHybridOrbitalTypes().e_SP3)
      {
        positions.Resize( 4);
        // put the H's at the vertices of a tetrahedron
        // Thus, the bond lengths should have coordinates of bond length / sqrt( 3) is
        const double bond_length_norm( BOND_LENGTH / std::sqrt( 3.0));
        positions( 0) = linal::Vector3D( bond_length_norm, bond_length_norm, bond_length_norm);
        positions( 1) = linal::Vector3D( bond_length_norm, -bond_length_norm, -bond_length_norm);
        positions( 2) = linal::Vector3D( -bond_length_norm, bond_length_norm, -bond_length_norm);
        positions( 3) = linal::Vector3D( -bond_length_norm, -bond_length_norm, bond_length_norm);
      }
      else if( ORBITAL_TYPE == GetHybridOrbitalTypes().e_SP2)
      {
        positions.Resize( 3);
        // put the vertices at the corners of an equalateral triangle centered at position
        positions( 0) = linal::Vector3D( BOND_LENGTH, 0.0, 0.0);
        positions( 1) = linal::Vector3D( -BOND_LENGTH / 2.0, BOND_LENGTH * std::sqrt( 3.0) / 2.0, 0.0);
        positions( 2) = linal::Vector3D( -BOND_LENGTH / 2.0, -BOND_LENGTH * std::sqrt( 3.0) / 2.0, 0.0);
      }
      else if( ORBITAL_TYPE == GetHybridOrbitalTypes().e_SP)
      {
        positions.Resize( 2);
        positions( 0) = linal::Vector3D( BOND_LENGTH, 0.0, 0.0);
        positions( 1) = linal::Vector3D( -BOND_LENGTH, 0.0, 0.0);
      }

      return positions;
    }

    //! @brief Set any H coordinates that are undefined
    void HydrogensHandler::SetUndefinedHCoordinates( const AtomComplete &ATOM, AtomVector< AtomComplete> &VECT)
    {
      // get the position of this atom
      const linal::Vector3D &position( ATOM.GetPosition());

      if( !position.IsDefined())
      {
        return;
      }

      // get a const reference to the bonds of this atom
      const storage::Vector< BondConformational> &bonds( ATOM.GetBonds());

      // count number of bonds to H where the H have undefined coordinates
      storage::Vector< size_t> atom_indices_undefined_hydrogen_positions;
      storage::Vector< double> bond_lengths_undefined;

      // split bonds according to whether they point to well-defined positions, hydrogen or not
      for( auto itr( bonds.Begin()), itr_end( bonds.End()); itr != itr_end; ++itr)
      {
        if( !itr->GetTargetAtom().GetPosition().IsDefined())
        {
          atom_indices_undefined_hydrogen_positions.PushBack( VECT.GetAtomIndex( itr->GetTargetAtom()));
          bond_lengths_undefined.PushBack
          (
            BondLengths::GetBondLength
            (
              ATOM.GetAtomType(),
              itr->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromatic),
              itr->GetTargetAtom().GetAtomType()
            )
          );
        }
      }

      if( atom_indices_undefined_hydrogen_positions.IsEmpty())
      {
        return;
      }

      util::SiPtrVector< const linal::Vector3D> defined_neighbor_positions;
      util::SiPtrVector< const AtomConformationalInterface> defined_neighbors;
      for( auto itr( bonds.Begin()), itr_end( bonds.End()); itr != itr_end; ++itr)
      {
        if( itr->GetTargetAtom().GetPosition().IsDefined())
        {
          defined_neighbor_positions.PushBack( itr->GetTargetAtom().GetPosition());
          defined_neighbors.PushBack( itr->GetTargetAtom());
        }
      }

      // get the hybrid orbital for the atom type
      HybridOrbitalType orbital_type( ATOM.GetAtomType()->GetHybridOrbitalType());

      // if the orbital type is unhybridized (this really shouldn't ever happen), then approximate it by considering it SP3
      if( orbital_type == GetHybridOrbitalTypes().e_Unhybridized)
      {
        orbital_type = GetHybridOrbitalTypes().e_SP3;
      }

      const size_t nr_known_bonds( defined_neighbor_positions.GetSize());
      const size_t n_valences( ATOM.GetNumberValenceBonds());
      const size_t nr_missing( atom_indices_undefined_hydrogen_positions.GetSize() + n_valences);
      for( size_t i( 0); i < n_valences; ++i)
      {
        bond_lengths_undefined.PushBack
        (
          BondLengths::GetBondLength( ATOM.GetAtomType(), 1, GetAtomTypes().H_S)
        );
      }
      storage::Vector< linal::Vector3D> positions;
      positions.AllocateMemory( nr_missing);
      // 0 known bonds, just make nr_missing_bonds, with the right angle for the hybridization
      if( nr_known_bonds == 0)
      {
        positions = GetIdealizedGeometry( orbital_type, 1.0);

        // translate the points such that they are relative to position
        auto itr_und( bond_lengths_undefined.Begin());
        for
        (
          storage::Vector< linal::Vector3D>::iterator itr( positions.Begin()), itr_end( positions.End());
          itr != itr_end;
          ++itr, ++itr_und
        )
        {
          *itr *= *itr_und;
          *itr += position;
        }
      }
      else if( nr_known_bonds == 1)
      {
        if( orbital_type == GetHybridOrbitalTypes().e_SP)
        {
          positions.PushBack
          (
            linal::CoordinatesLinear( position, *defined_neighbor_positions( 0), bond_lengths_undefined( 0))
          );
        }
        else if( orbital_type == GetHybridOrbitalTypes().e_SP2)
        {
          const storage::Vector< BondConformational> &bonds_neighbor( defined_neighbors( 0)->GetBonds());
          storage::Vector< linal::Vector3D> neighboring_positions;

          for
          (
            storage::Vector< BondConformational>::const_iterator
              itr_neighbors( bonds_neighbor.Begin()), itr_neighbors_end( bonds_neighbor.End());
            itr_neighbors != itr_neighbors_end;
            ++itr_neighbors
          )
          {
            if( &itr_neighbors->GetTargetAtom() != &ATOM && itr_neighbors->GetTargetAtom().GetPosition().IsDefined())
            {
              neighboring_positions.PushBack( itr_neighbors->GetTargetAtom().GetPosition());
              break;
            }
          }

//          if( !neighboring_positions.IsEmpty())
//          {
//            // get the anchor points for bond rotation, i.e. indices of atoms that make the bond of interest
//            const linal::Vector3D &start_point( *defined_neighbor_positions( 0));
//
//            // get the unit vector from the start point's other neighbor to this
//            positions.PushBack( position - linal::UnitVector( start_point, neighboring_positions( 0)) * bond_lengths_undefined( 0));
//          }
//          else
//          {
            // coordinate 1
            positions.PushBack
            (
              linal::CoordinatesAngle
              (
                position,
                *defined_neighbor_positions( 0),
                linal::Vector3D( bond_lengths_undefined( 0)),
                bond_lengths_undefined( 0),
                120.0 / 180.0 * math::g_Pi,
                120.0 / 180.0 * math::g_Pi,
                false,
                linal::Vector3D( 0.0)
              )
            );
//          }
          if( bond_lengths_undefined.GetSize() > size_t( 1))
          {
            // coordinate 2
            positions.PushBack
            (
              linal::CoordinatesTrigonal
              (
                position,
                *defined_neighbor_positions( 0),
                positions.LastElement(),
                bond_lengths_undefined( 1)
              )
            );
          }
        }
        else // SP3
        {
          // coordinate 1
          positions.PushBack
          (
            linal::CoordinatesAngle
            (
              position,
              *defined_neighbor_positions( 0),
              linal::Vector3D( bond_lengths_undefined( 0)),
              bond_lengths_undefined( 0),
              109.5 / 180.0 * math::g_Pi,
              109.5 / 180.0 * math::g_Pi,
              false,
              linal::Vector3D( 0.0)
            )
          );

          linal::Vector3D offset
          (
            std::sin( 54.75 / 180 * math::g_Pi) *
            linal::CrossProduct
            (
              *defined_neighbor_positions( 0) - position,
              positions.LastElement() - position
            ).Normalize()
          );

          if( nr_missing >= size_t( 2))
          {
            // helper coordinates
            linal::Vector3D foot_point_a
            (
              linal::CoordinatesTrigonal
              (
                position,
                *defined_neighbor_positions( 0),
                positions.LastElement(),
                bond_lengths_undefined( 1) * std::cos( 54.75 / 180 * math::g_Pi)
              )
            );

            // coordinate 2
            positions.PushBack( foot_point_a + bond_lengths_undefined( 1) * offset);
          }

          if( nr_missing == size_t( 3))
          {
            linal::Vector3D foot_point_b
            (
              linal::CoordinatesTrigonal
              (
                position,
                *defined_neighbor_positions( 0),
                positions( 0),
                bond_lengths_undefined( 2) * std::cos( 54.75 / 180 * math::g_Pi)
              )
            );

            // coordinate 3
            positions.PushBack( foot_point_b - bond_lengths_undefined( 2) * offset);
          }
        }
      }
      else if( nr_known_bonds == 2) // 2 known bonds
      {
        if( orbital_type == GetHybridOrbitalTypes().e_SP2) // trigonal
        {
          positions.PushBack
          (
            linal::CoordinatesTrigonal
            (
              position,
              *defined_neighbor_positions( 0),
              *defined_neighbor_positions( 1),
              bond_lengths_undefined( 0)
            )
          );
        }
        else // tetrahedral
        {
          // helper coordinates
          linal::Vector3D foot_point_a
          (
            linal::CoordinatesTrigonal
            (
              position,
              *defined_neighbor_positions( 0),
              *defined_neighbor_positions( 1),
              bond_lengths_undefined( 0) * std::cos( 54.75 / 180 * math::g_Pi)
            )
          );

          linal::Vector3D offset
          (
            std::sin( 54.75 / 180 * math::g_Pi) *
            linal::CrossProduct
            (
              *defined_neighbor_positions( 0) - position,
              *defined_neighbor_positions( 1) - position
            ).Normalize()
          );

          // coordinate 1
          positions.PushBack( foot_point_a + bond_lengths_undefined( 0) * offset);

          if( nr_missing == size_t( 2))
          {
            linal::Vector3D foot_point_b
            (
              linal::CoordinatesTrigonal
              (
                position,
                *defined_neighbor_positions( 0),
                *defined_neighbor_positions( 1),
                bond_lengths_undefined( 1) * std::cos( 54.75 / 180 * math::g_Pi)
              )
            );
            // coordinate 2
            positions.PushBack( foot_point_b - bond_lengths_undefined( 1) * offset);
          }
        }
      }
      else if( nr_known_bonds == 3) // tetrahedral geometry, only missing one bond
      {
        positions.PushBack
        (
          linal::CoordinatesTetrahedral
          (
            position,
            *defined_neighbor_positions( 0),
            *defined_neighbor_positions( 1),
            *defined_neighbor_positions( 2),
            bond_lengths_undefined( 0)
          )
        );
      }

      positions.Resize( atom_indices_undefined_hydrogen_positions.GetSize());
      auto itr_indices( atom_indices_undefined_hydrogen_positions.Begin());
      for
      (
        storage::Vector< linal::Vector3D>::const_iterator itr( positions.Begin()), itr_end( positions.End());
        itr != itr_end;
        ++itr, ++itr_indices
      )
      {
        VECT( *itr_indices).SetPosition( *itr);
      }
      return;
    }

    //! @brief Set any H coordinates that are undefined
    void HydrogensHandler::UpdateHCoordinates( const AtomComplete &ATOM, AtomVector< AtomComplete> &VECT)
    {
      size_t n_h( 0);
      for( auto itr( ATOM.GetBonds().Begin()), itr_end( ATOM.GetBonds().End()); itr != itr_end; ++itr)
      {
        if( itr->GetTargetAtom().GetAtomType() == GetAtomTypes().H_S)
        {
          VECT( VECT.GetAtomIndex( itr->GetTargetAtom())).SetPosition( linal::Vector3D( util::GetUndefinedDouble()));
          ++n_h;
        }
      }
      if( n_h)
      {
        SetUndefinedHCoordinates( ATOM, VECT);
      }
    }

  } // namespace chemistry
} // namespace bcl
