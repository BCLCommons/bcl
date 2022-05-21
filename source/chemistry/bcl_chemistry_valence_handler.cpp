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
#include "chemistry/bcl_chemistry_valence_handler.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_3d_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new ValenceHandler
    ValenceHandler *ValenceHandler::Clone() const
    {
      return new ValenceHandler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ValenceHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief determine the coordinates of all missing hydrogens for a given atom
    //! @param ATOM an AtomComplete that should have both type and position defined
    //! @return a storage vector containing the positions of all missing hydrogens
    storage::Vector< linal::Vector3D> ValenceHandler::DetermineCoordinates( const AtomConformationalInterface &ATOM)
    {
      // check whether the atom has the same # of bonds that the atom type has

      // get a const reference to the bonds of this atom
      const storage::Vector< BondConformational> &bonds( ATOM.GetBonds());
      const size_t nr_known_bonds( bonds.GetSize());

      const size_t nr_missing_bonds( ATOM.GetNumberValenceBonds());

      if( nr_missing_bonds == util::GetUndefined< size_t>())
      {
        return storage::Vector< linal::Vector3D>();
      }

      if( nr_missing_bonds == 0)
      {
        return storage::Vector< linal::Vector3D>();
      }

      storage::Vector< linal::Vector3D> positions;
      positions.AllocateMemory( nr_missing_bonds);

      // store the expected bond length, which is the covalent radius of H plus the covalent radius of ATOM
      const double bond_length( 1.0);

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

      return positions;
    }

    //! @brief get the idealized geometry for a particular hybrid orbital type and bond length
    //! @param ORBITAL_TYPE the type of orbital associated with the geometry
    //! @param BOND_LENGTH the length of the bonds in the geometry
    storage::Vector< linal::Vector3D> ValenceHandler::GetIdealizedGeometry
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

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ValenceHandler::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ValenceHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
