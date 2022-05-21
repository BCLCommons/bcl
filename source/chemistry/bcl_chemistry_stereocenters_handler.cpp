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
#include "chemistry/bcl_chemistry_stereocenters_handler.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_valence_handler.h"
#include "linal/bcl_linal_matrix3x3.h"
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
    //! @return pointer to new StereocentersHandler
    StereocentersHandler *StereocentersHandler::Clone() const
    {
      return new StereocentersHandler( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &StereocentersHandler::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    namespace
    {
      //! @brief Returns the substituents of a specified base atom in ordered by priority
      //! @param BASE Atom that the sorted atoms are all connected to; need not be a stereocenter
      //! @return vector of pairs, ordered by priorities (which may be equal)
      storage::Vector< linal::Vector3D> GetSubstituentsUnitVectorsOrderedByPriorityForRingIsomerism
      (
        const AtomConformationalInterface &BASE
      )
      {
        const size_t n_ring_bonds
        (
          BASE.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1))
        );
        storage::Vector< linal::Vector3D> unit_vecs( 2);
        if( n_ring_bonds == size_t( 2))
        {
          if( BASE.GetNumberValenceBonds())
          {
            // trivial case, we have a hydrogen, no need for sorting substituents
            for
            (
              storage::Vector< BondConformational>::const_iterator
                itr_bonded( BASE.GetBonds().Begin()),
                itr_bonded_end( BASE.GetBonds().End());
              itr_bonded != itr_bonded_end;
              ++itr_bonded
            )
            {
              if( itr_bonded->GetBondType()->IsBondInRing())
              {
                continue;
              }
              unit_vecs( 0) = linal::UnitVector( BASE.GetPosition(), itr_bonded->GetTargetAtom().GetPosition());
              break;
            }
            unit_vecs( 1) = linal::UnitVector( BASE.GetPosition(), ValenceHandler::DetermineCoordinates( BASE)( 0));
            return unit_vecs;
          }
          // just get the two chain bonds and order them properly
          size_t bond_nr( 0);

          storage::Vector< SubstituentConformational> subs( 2);
          for
          (
            storage::Vector< BondConformational>::const_iterator
              itr_bonded( BASE.GetBonds().Begin()),
              itr_bonded_end( BASE.GetBonds().End());
            itr_bonded != itr_bonded_end;
            ++itr_bonded
          )
          {
            if( itr_bonded->GetBondType()->IsBondInRing())
            {
              continue;
            }

            subs( bond_nr++) = SubstituentConformational( itr_bonded->GetTargetAtom());
          }
          if( subs( 0) < subs( 1))
          {
            unit_vecs( 0) = linal::UnitVector( BASE.GetPosition(), subs( 0).GetRootAtom()->GetPosition());
            unit_vecs( 1) = linal::UnitVector( BASE.GetPosition(), subs( 1).GetRootAtom()->GetPosition());
          }
          else
          {
            unit_vecs( 0) = linal::UnitVector( BASE.GetPosition(), subs( 1).GetRootAtom()->GetPosition());
            unit_vecs( 1) = linal::UnitVector( BASE.GetPosition(), subs( 0).GetRootAtom()->GetPosition());
          }
          return unit_vecs;
        }
        else if( n_ring_bonds >= size_t( 3))
        {
          storage::Vector< SubstituentConformational> holders;
          holders.AllocateMemory( BASE.GetBonds().GetSize());
          util::SiPtr< const AtomConformationalInterface> chain_atom;
          for
          (
            storage::Vector< BondConformational>::const_iterator
              itr_bonded( BASE.GetBonds().Begin()),
              itr_bonded_end( BASE.GetBonds().End());
            itr_bonded != itr_bonded_end;
            ++itr_bonded
          )
          {
            if( itr_bonded->GetBondType()->IsBondInRing())
            {
              holders.PushBack( SubstituentConformational( itr_bonded->GetTargetAtom()));
            }
            else
            {
              chain_atom = itr_bonded->GetTargetAtom();
            }
          }

          // sort the vector of substituents
          holders.Sort( std::less< SubstituentConformational>());
          util::SiPtr< const SubstituentConformational> last_unique_substituent_itr;

          util::SiPtr< const AtomConformationalInterface> unique_atom;

          if( holders( 0) < holders( 1))
          {
            unique_atom = holders( 0).GetRootAtom();
          }
          else if( holders( 0) < holders( 2))
          {
            unique_atom = holders( 2).GetRootAtom();
          }
          else
          {
            // all ring directions identical
            BCL_Exit( "Unknown ring issue", -1);
          }
          unit_vecs( 0) = linal::UnitVector( BASE.GetPosition(), unique_atom->GetPosition());
          if( chain_atom.IsDefined())
          {
            unit_vecs( 1) = linal::UnitVector( BASE.GetPosition(), chain_atom->GetPosition());
          }
          else
          {
            unit_vecs( 1) = linal::UnitVector( BASE.GetPosition(), ValenceHandler::DetermineCoordinates( BASE)( 0));
          }
        }
        else
        {
          return unit_vecs;
        }
        return unit_vecs;
      }
    }

    //! @brief Recalculate stereocenters given a conformation
    //! @param ATOMS atoms for which stereocenter values are determined
    //! @return storage vector with stereocenter values referenced by atom index
    linal::Vector< float> StereocentersHandler::CalculateFromConformation
    (
      const iterate::Generic< const AtomConformationalInterface> &MOLECULE
    )
    {
      linal::Vector< float> stereocenters_vector( MOLECULE.GetSize(), 0.0);
      linal::Vector< float>::iterator itr_stereo( stereocenters_vector.Begin());

      // Iterate through all atoms in MOLECULE
      storage::Vector< size_t> unknown_ring_chirality;
      for
      (
        iterate::Generic< const AtomConformationalInterface> itr_atoms( MOLECULE);
        itr_atoms.NotAtEnd();
        ++itr_atoms, ++itr_stereo
      )
      {
        *itr_stereo = CalculateFromConformation( *itr_atoms);
        if( DescriptorToChirality( *itr_stereo) == e_UnknownRingChirality)
        {
          unknown_ring_chirality.PushBack( itr_atoms.GetPosition());
        }
      }
      if( unknown_ring_chirality.GetSize() == size_t( 1))
      {
        stereocenters_vector( unknown_ring_chirality( 0)) = 0.0;
      }
      else if( !unknown_ring_chirality.IsEmpty())
      {
        iterate::Generic< const AtomConformationalInterface> mol_two( MOLECULE);
        ++mol_two;
        const size_t n_bytes_atom( ( const char *)( &*mol_two) - ( const char *)( &*MOLECULE));
        BCL_Assert
        (
          ( ( const char *)( MOLECULE.End().operator->()) - ( const char *)( MOLECULE.operator->())) / n_bytes_atom
          == MOLECULE.GetSize(),
          "Generic iterator did not reference an atom vector!"
        );
        storage::Vector< size_t> ring_system_id( MOLECULE.GetSize(), util::GetUndefined< size_t>());
        storage::List< size_t> queue;
        size_t n_ring_systems( 0);
        const char *molbegin( ( const char *)( &*MOLECULE));
        for
        (
          auto itr_unk( unknown_ring_chirality.Begin()), itr_unk_end( unknown_ring_chirality.End());
          itr_unk != itr_unk_end;
          ++itr_unk
        )
        {
          if( util::IsDefined( ring_system_id( *itr_unk)))
          {
            continue;
          }
          queue.PushBack( *itr_unk);
          while( !queue.IsEmpty())
          {
            storage::List< size_t> new_queue;
            for( auto itrq( queue.Begin()), itrq_end( queue.End()); itrq != itrq_end; ++itrq)
            {
              if( util::IsDefined( ring_system_id( *itrq)))
              {
                continue;
              }
              ring_system_id( *itrq) = n_ring_systems;
              iterate::Generic< const AtomConformationalInterface> itr_mol( MOLECULE);
              itr_mol.GotoPosition( *itrq);
              for
              (
                auto itr_bond( itr_mol->GetBonds().Begin()), itr_bond_end( itr_mol->GetBonds().End());
                itr_bond != itr_bond_end;
                ++itr_bond
              )
              {
                if( !itr_bond->GetBondType()->IsBondInRing())
                {
                  continue;
                }
                const size_t atom_index( ( ( const char *)( &itr_bond->GetTargetAtom()) - molbegin) / n_bytes_atom);
                if( util::IsDefined( ring_system_id( atom_index)))
                {
                  continue;
                }
                new_queue.PushBack( atom_index);
              }
            }
            std::swap( queue.InternalData(), new_queue.InternalData());
          }

          ++n_ring_systems;
        }
        storage::Vector< size_t> n_emergent_stereo( n_ring_systems, size_t( 0));
        for
        (
          auto itr_unk( unknown_ring_chirality.Begin()), itr_unk_end( unknown_ring_chirality.End());
          itr_unk != itr_unk_end;
          ++itr_unk
        )
        {
          ++n_emergent_stereo( ring_system_id( *itr_unk));
        }
        size_t n_ring_chirality( 0);
        storage::Vector< storage::Vector< size_t> > ring_ids_chiral( n_ring_systems);
        for
        (
          auto itr_unk( unknown_ring_chirality.Begin()), itr_unk_end( unknown_ring_chirality.End());
          itr_unk != itr_unk_end;
          ++itr_unk
        )
        {
          if( n_emergent_stereo( ring_system_id( *itr_unk)) == size_t( 1))
          {
            stereocenters_vector( *itr_unk) = double( 0.0);
          }
          else
          {
            ++n_ring_chirality;
            ring_ids_chiral( ring_system_id( *itr_unk)).PushBack( *itr_unk);
          }
        }
        for( size_t ring_sys( 0); ring_sys < n_ring_systems; ++ring_sys)
        {
          if( ring_ids_chiral( ring_sys).IsEmpty())
          {
            continue;
          }
          if( ring_ids_chiral( ring_sys).GetSize() == size_t( 2))
          {
            const size_t id_a( ring_ids_chiral( ring_sys)( 0));
            const size_t id_b( ring_ids_chiral( ring_sys)( 1));

            util::SiPtr< const linal::Vector3D> a_pos, b_pos;
            iterate::Generic< const AtomConformationalInterface> a_itr( MOLECULE), b_itr( MOLECULE);
            a_itr.GotoPosition( id_a);
            a_pos = util::ToSiPtr( a_itr->GetPosition());
            b_itr.GotoPosition( id_b);
            b_pos = util::ToSiPtr( b_itr->GetPosition());
            storage::Vector< linal::Vector3D> substituents_a( GetSubstituentsUnitVectorsOrderedByPriorityForRingIsomerism( *a_itr));
            storage::Vector< linal::Vector3D> substituents_b( GetSubstituentsUnitVectorsOrderedByPriorityForRingIsomerism( *b_itr));
            if( linal::ProjAngleCosinus( substituents_a( 0), substituents_b( 0)) > linal::ProjAngleCosinus( substituents_a( 1), substituents_b( 0)))
            {
              stereocenters_vector( id_a) = stereocenters_vector( id_b) = 4.0;
            }
            else
            {
              stereocenters_vector( id_a) = stereocenters_vector( id_b) = 5.0;
            }
          }
          else if( ring_ids_chiral( ring_sys).GetSize() == size_t( 3))
          {
            const size_t id_a( ring_ids_chiral( ring_sys)( 0));
            const size_t id_b( ring_ids_chiral( ring_sys)( 1));
            const size_t id_c( ring_ids_chiral( ring_sys)( 2));

            util::SiPtr< const linal::Vector3D> a_pos, b_pos, c_pos;
            iterate::Generic< const AtomConformationalInterface> a_itr( MOLECULE), b_itr( MOLECULE), c_itr( MOLECULE);
            a_itr.GotoPosition( id_a);
            a_pos = util::ToSiPtr( a_itr->GetPosition());
            b_itr.GotoPosition( id_b);
            b_pos = util::ToSiPtr( b_itr->GetPosition());
            c_itr.GotoPosition( id_c);
            c_pos = util::ToSiPtr( c_itr->GetPosition());
            storage::Vector< linal::Vector3D> substituents_a( GetSubstituentsUnitVectorsOrderedByPriorityForRingIsomerism( *a_itr));
            storage::Vector< linal::Vector3D> substituents_b( GetSubstituentsUnitVectorsOrderedByPriorityForRingIsomerism( *b_itr));
            storage::Vector< linal::Vector3D> substituents_c( GetSubstituentsUnitVectorsOrderedByPriorityForRingIsomerism( *c_itr));
            if
            (
              linal::ProjAngleCosinus( substituents_a( 0), substituents_b( 0))
              > linal::ProjAngleCosinus( substituents_a( 1), substituents_b( 0))
              && linal::ProjAngleCosinus( substituents_a( 0), substituents_c( 0))
                > linal::ProjAngleCosinus( substituents_a( 1), substituents_c( 0))
            )
            {
              stereocenters_vector( id_a) = stereocenters_vector( id_b) = stereocenters_vector( id_c) = 4.0;
            }
          }
        }
      }

      return stereocenters_vector;
    }

    //! @brief Recalculate stereocenters given a conformation
    //! @param MOLECULE SmallMolecule for which stereocenter values of for all atoms are determined
    //! @return storage vector with stereocenter values referenced by atom index
    float StereocentersHandler::CalculateFromConformation
    (
      const AtomConformationalInterface &ATOM
    )
    {
      // Check if this atom is a possible stereocenter
      // Current definition of stereocenter is only includes standard R/S chirality (e.g. 4 substituents)
      if( ATOM.GetAtomType()->GetNumberBonds() != size_t( 4))
      {
        return 0.0;
      }

      // Get pointers to all atoms connected with the current atom.
      const storage::Vector< BondConformational> &bonds( ATOM.GetBonds());

      // 3 explicitly bonded atoms -> 1 implicit h, so a stereocenter is possible
      // < 3 explicitly bonded atoms -> more than 1 implicit h -> no stereocenter
      if( bonds.GetSize() < size_t( 3))
      {
        // more than one implicit h -> no stereocenter
        return 0.0;
      }

      const size_t n_cov_h( ATOM.GetNumberCovalentlyBoundHydrogens());
      if( n_cov_h + ATOM.GetNumberValenceBonds() >= size_t( 2))
      {
        // more than 1 hydrogen / valence, not a stereocenter, so continue
        return 0.0;
      }

      // get the atoms sorted by priority
      const util::SiPtrVector< const AtomConformationalInterface>
        atoms_by_priority( GetUniqueConnectedSubstituentsByPriority( ATOM));

      // if the unique connected substituents is not the same size as the number of bonded atoms, then there is no
      // chirality
      if( atoms_by_priority.GetSize() == bonds.GetSize())
      {
        // Create a matrix that will hold the x,y, and z coordinates for the three atoms of highest priority.
        linal::Matrix3x3< float> xyz_coordinates( 0.0); // make a matrix of size 3 X 3
        const linal::Vector3D &root_position( ATOM.GetPosition());

        // put the positions (relative to the root atom) of the 3 highest priority substituents into a matrix
        // the sign of the determinant of the matrix will give us R and S
        for( size_t index( 0), size( 3); index < size; ++index)
        {
          // get the atom position out of the vector, which has been sorted by priority
          const linal::Vector3D &position( atoms_by_priority( index)->GetPosition());

          // put the positions of this atom relative to the root atom into the rows of the matrix
          xyz_coordinates( index, 0) = position( 0) - root_position( 0);
          xyz_coordinates( index, 1) = position( 1) - root_position( 1);
          xyz_coordinates( index, 2) = position( 2) - root_position( 2);
        }

        // Calculate the determinant of a 3 X 3 matrix that has rows sorted in descending order of priority.
        // Opposite orders will have opposite signs. The sign is assigned to a value of 1 and returned as the
        // value of the stereocenter.
        const float determinant( xyz_coordinates.Determinant());
        if( determinant != 0.0)
        {
          return math::Sign( float( 1.0), determinant);
        }
        // possibly planar chirality or 2D structure
        return 2.0;
      }

      size_t n_nonaro_ring( 0), n_aro_ring( 0), n_chain( ATOM.GetNumberValenceBonds());
      for( auto itr( atoms_by_priority.Begin()), itr_end( atoms_by_priority.End()); itr != itr_end; ++itr)
      {
        auto bond_type( ATOM.FindBondTo( **itr)->GetBondType());
        if( bond_type->IsBondInRing())
        {
          if( bond_type->GetConjugation() == ConstitutionalBondTypeData::e_Aromatic)
          {
            ++n_aro_ring;
          }
          else
          {
            ++n_nonaro_ring;
          }
        }
        else
        {
          ++n_chain;
        }
      }

      // potentially ring chiral, need to see if any other atoms in the ring are also potentially ring-chiral
      return n_chain == size_t( 2) && n_nonaro_ring == size_t( 1) && ATOM.GetElementType() != GetElementTypes().e_Nitrogen ? 3.0 : 0.0;
    }

    //! @brief convert a descriptor value to the chirality type
    //! @param CHIRALITY_DESCRIPTOR a value returned by CalculateFromConformation
    //! @return the chemistry::Chirality value
    Chirality StereocentersHandler::DescriptorToChirality( const float &VALUE)
    {
      if( VALUE == 0.0) // non-chiral
      {
        return e_NonChiral;
      }
      else if( VALUE == 1.0) // R
      {
        return e_RChirality;
      }
      else if( VALUE == -1.0) // S
      {
        return e_SChirality;
      }
      else if( VALUE == 2.0) // Chiral centers, unspecified type
      {
        return e_Chiral;
      }
      else if( VALUE == 3.0)
      {
        return e_UnknownRingChirality;
      }
      else if( VALUE == 4.0)
      {
        return e_CisRingChirality;
      }
      //else if( VALUE == 5.0)
      return e_TransRingChirality;
    }

    //! @brief Add R/S chirality information to the configuration, given a conformation
    //! @param MOLECULE Conformation upon which to add chirality information
    void StereocentersHandler::AddChiralityFromConformation( AtomVector< AtomComplete> &MOLECULE)
    {
      // determine the chirality of all the atoms
      linal::Vector< float> stereochem
      (
        CalculateFromConformation
        (
          iterate::Generic< const AtomConformationalInterface>( MOLECULE.Begin(), MOLECULE.End())
        )
      );

      // get an iterator to the stereochem vector
      linal::Vector< float>::const_iterator itr_stereo( stereochem.Begin());

      // add the chirality info to each atom with configuration
      for
      (
        AtomVector< AtomComplete>::iterator itr_atom( MOLECULE.Begin()), itr_atom_end( MOLECULE.End());
        itr_atom != itr_atom_end;
        ++itr_stereo, ++itr_atom
      )
      {
        itr_atom->SetChirality( DescriptorToChirality( *itr_stereo));
      }
    }

    //! @brief Returns the substituents of a specified base atom in ordered by priority
    //! @param BASE Atom that the sorted atoms are all connected to; need not be a stereocenter
    //! @return vector of pairs, ordered by priorities (which may be equal)
    storage::Vector< SubstituentConformational> StereocentersHandler::GetSubstituentsOrderedByPriority
    (
      const AtomConformationalInterface &BASE
    )
    {
      storage::Vector< SubstituentConformational> holders;
      holders.AllocateMemory( BASE.GetBonds().GetSize());
      for
      (
        storage::Vector< BondConformational>::const_iterator
          itr_bonded( BASE.GetBonds().Begin()),
          itr_bonded_end( BASE.GetBonds().End());
        itr_bonded != itr_bonded_end;
        ++itr_bonded
      )
      {
        holders.PushBack( SubstituentConformational( itr_bonded->GetTargetAtom()));
      }

      // sort the vector of substituents
      holders.Sort( std::less< SubstituentConformational>());
      return holders;
    }

    //! @brief Add R/S chirality information to the configuration, given a conformation
    //! @param MOLECULE Conformation upon which to add chirality information
    void StereocentersHandler::UpdateChiralityFromConformation( AtomComplete &ATOM)
    {
      ATOM.SetChirality( DescriptorToChirality( CalculateFromConformation( ATOM)));
    }

    //! @brief Returns the substituents of a specified base atom in ordered by priority, skipping any substituents that are equal
    //! @param BASE Atom that the sorted atoms are all connected to; need not be a stereocenter
    //! @return ShPtrVector of Atoms that are connected to BASE sorted by priority
    util::SiPtrVector< const AtomConformationalInterface> StereocentersHandler::GetUniqueConnectedSubstituentsByPriority
    (
      const AtomConformationalInterface &BASE
    )
    {
      storage::Vector< SubstituentConformational> substituents( GetSubstituentsOrderedByPriority( BASE));
      if( substituents.IsEmpty())
      {
        return util::SiPtrVector< const AtomConformationalInterface>();
      }

      storage::Vector< SubstituentConformational>::const_iterator last_unique_substituent_itr( substituents.Begin());

      util::SiPtrVector< const AtomConformationalInterface> unique_atoms;
      unique_atoms.AllocateMemory( substituents.GetSize());
      // add the first unique substituent to the unique atoms
      unique_atoms.PushBack( last_unique_substituent_itr->GetRootAtom());

      // make an iterator to walk ahead of last_unique_substituent_itr
      storage::Vector< SubstituentConformational>::const_iterator itr( last_unique_substituent_itr);
      ++itr;
      for( storage::Vector< SubstituentConformational>::const_iterator itr_end( substituents.End()); itr != itr_end; ++itr)
      {
        if( *last_unique_substituent_itr < *itr) // is *itr a unique substituent?
        {
          // *itr is a unique substituent
          // update the last unique substituent itr to itr, which is unique
          last_unique_substituent_itr = itr;
          unique_atoms.PushBack( last_unique_substituent_itr->GetRootAtom());
        }
      }

      return unique_atoms;
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
    std::istream &StereocentersHandler::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &StereocentersHandler::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
