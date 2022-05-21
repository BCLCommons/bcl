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
#include "chemistry/bcl_chemistry_mutate_chirality.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_mutate_bond_angles.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "graph/bcl_graph_connectivity.h"
#include "io/bcl_io_ofstream.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"

// external includes - sorted alphabetically

//#define VALIDATE_BOND_ANGLES_DO_NOT_CHANGE
//#define VALIDATE_DIHEDRAL_ANGLES

namespace bcl
{
  namespace chemistry
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateChirality::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateChirality())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor taking the member variable parameters
    //! @param FRAGMENT fragment which the mutate mutates
    //! @param UNKOWN_CHIRAL_ONLY boolean to specify whether to change unknown chiral only
    //! @param MOLECULE molecule of interest
    MutateChirality::MutateChirality
    (
      const FragmentComplete &MOLECULE,
      bool UNKOWN_CHIRAL_ONLY
    ) :
      m_Graph
      (
        ConformationGraphConverter
        (
          ConformationGraphConverter::e_CIPPriorityHighToLow,
          ConfigurationalBondTypeData::e_IsInRing
        )( MOLECULE)
      ),
      m_GraphBondOrder
      (
        ConformationGraphConverter
        (
          ConformationGraphConverter::e_CIPPriorityHighToLow,
          ConfigurationalBondTypeData::e_BondOrder
        )( MOLECULE)
      ),
      m_ChangeUnknownChiralOnly( false),
      m_BondsToDownstreamAtoms( MOLECULE.GetSize()),
      m_OnlyHasIsometry( MOLECULE.GetSize(), size_t( 0)),
      m_BondInfo( MOLECULE.GetBondInfo())
    {
      size_t atom_index( 0);

      // add the chirality info to each atom with configuration
      bool had_unknown( false);
      for
      (
        auto atom_itr( MOLECULE.GetAtomsIterator());
        atom_itr.NotAtEnd();
        ++atom_index, ++atom_itr
      )
      {
        if( atom_itr->GetChirality() != e_NonChiral)
        {
          if( atom_itr->GetChirality() == e_UnknownChirality)
          {
            m_UnknownChiralCenters.PushBack( atom_index);
            had_unknown = true;
          }
          m_ChiralCenters.PushBack( atom_index);
        }
        else if
        (
          atom_itr->CountNonValenceBondsWithProperty
          (
            ConfigurationalBondTypeData::e_IsIsometric,
            size_t( 1)
          ) == size_t( 1)
        )
        {
          m_OnlyHasIsometry( atom_index) = size_t( 1);
        }
        for
        (
          auto itr_bnd( atom_itr->GetBonds().Begin()), itr_bnd_end( atom_itr->GetBonds().End());
          itr_bnd != itr_bnd_end;
          ++itr_bnd
        )
        {
          const size_t atom_b_index( MOLECULE.GetAtomIndex( itr_bnd->GetTargetAtom()));
          if( itr_bnd->GetBondType()->IsBondInRing())
          {
            storage::Vector< size_t> tmp( size_t( 1), atom_b_index);
            for
            (
              auto itr_nxt( m_Graph.GetNeighborIndices( atom_b_index).Begin()),
                   itr_nxt_end( m_Graph.GetNeighborIndices( atom_b_index).End());
              itr_nxt != itr_nxt_end;
              ++itr_nxt
            )
            {
              if( *itr_nxt != atom_index)
              {
                tmp.PushBack( *itr_nxt);
              }
            }
            m_BondsToDownstreamAtoms( atom_index)[ atom_b_index] = tmp;
            continue;
          }
          util::ShPtr< storage::Vector< size_t> > connected_vertices_sp
          (
            graph::Connectivity::BreadthFirstSearchDirectedEdge( m_Graph, atom_b_index, atom_index)
          );
          m_BondsToDownstreamAtoms( atom_index)[ atom_b_index] = *connected_vertices_sp;
        }
      }
      if( had_unknown && UNKOWN_CHIRAL_ONLY)
      {
        m_ChangeUnknownChiralOnly = true;
      }
    }

    //! @brief Clone function
    //! @return pointer to new MutateChirality
    MutateChirality *MutateChirality::Clone() const
    {
      return new MutateChirality( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateChirality::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns atom indices that are connected to different bonds
    //! @return atom indices that are connected to different rotatable bonds
    const storage::Vector< size_t> &MutateChirality::GetChiralCenters() const
    {
      return m_ChiralCenters;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking a conformation and returns a mutated conformation
    //! @param MOLECULE conformation of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< FragmentComplete> MutateChirality::operator()
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      if( m_ChiralCenters.IsEmpty())
      {
        // create new fragment from the updated coordinaes
        util::ShPtr< FragmentComplete> new_molecule
        (
          new FragmentComplete( MOLECULE)
        );

        return math::MutateResult< FragmentComplete>( new_molecule, *this);
      }
      AtomVector< AtomComplete> atom_vec( MOLECULE.GetAtomVector());
      const storage::Vector< size_t> &vec( m_ChangeUnknownChiralOnly ? m_UnknownChiralCenters : m_ChiralCenters);
      for
      (
        auto chiral_center_itr( m_ChiralCenters.Begin()), chiral_center_itr_end( m_ChiralCenters.End());
        chiral_center_itr != chiral_center_itr_end;
        ++chiral_center_itr
      )
      {
        if( random::GetGlobalRandom().Boolean())
        {
          MutateSpecificChirality( atom_vec, *chiral_center_itr);
        }
      }
      util::ShPtr< FragmentComplete> new_molecule( new FragmentComplete( atom_vec, MOLECULE.GetName()));
      return math::MutateResult< FragmentComplete>( new_molecule, *this);
    }

    //! @brief get the downstream atoms of a particular atom or bond
    //! @param ATOMA the first atom in the bond
    //! @param ATOMB the second atom in the bond
    //! @return everything that is reachable after going one way through bond A->B. Not valid for ring bonds
    const storage::Vector< size_t> &MutateChirality::GetDownstreamAtoms( const size_t &ATOMA, const size_t &ATOMB) const
    {
      return m_BondsToDownstreamAtoms( ATOMA).GetValue( ATOMB);
    }

    //! @brief return a vector of atoms that are connected to a bond of interest
    //! @param ATOM_INDICES atom indices of the center bond around which connected bonds have to be found
    //! @param GRAPH constitution graph of the molecule of interest
    //! @param return a vector of atoms that are connected to a bond of interest
    FragmentComplete MutateChirality::MutateSpecificChirality
    (
      const FragmentComplete &MOLECULE,
      const size_t &ATOM_INDEX,
      const bool &ALTER_CHIRALITY_LABEL,
      const bool &FORCE_ALTER_LABEL
    ) const
    {
      AtomVector< AtomComplete> atom_vec( MOLECULE.GetAtomVector());
      bool success( this->MutateSpecificChirality( atom_vec, ATOM_INDEX, ALTER_CHIRALITY_LABEL, FORCE_ALTER_LABEL));
      FragmentComplete fc( atom_vec, MOLECULE.GetName());

#ifdef VALIDATE_BOND_ANGLES_DO_NOT_CHANGE
      auto ba( MOLECULE.GetBondAngles());
      auto ba2( fc.GetBondAngles());
      if( !math::EqualWithinTolerance( ba, ba2, 0.01))
      {
        BCL_Debug( linal::Vector< double>( ba.Begin(), ba.End()));
        BCL_Debug( linal::Vector< double>( ba2.Begin(), ba2.End()));
        BCL_Debug(ATOM_INDEX);
        BCL_Debug( linal::Vector< double>( ba.Begin(), ba.End()) - linal::Vector< double>(  ba2.Begin(), ba2.End()));
      }
#endif
#ifdef VALIDATE_DIHEDRAL_ANGLES
      PriorityDihedralAngles pda;
      auto saved_dihedral( pda( MOLECULE).First());
      auto new_dihedral( pda( fc).First());
      if( ( linal::Vector< double>( new_dihedral) - linal::Vector< double>( saved_dihedral)).Norm() > 0.01)
      {
        bool found_something_legit( false);

        auto moldih( MOLECULE.GetDihedralAngles());
        auto itr_moldih( moldih.Begin());
        for
        (
          iterate::Generic< const AtomConformationalInterface> itr_atom_b( fc.GetAtomsIterator());
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
                double dihedral_abcd
                (
                  linal::Dihedral
                  (
                    itr_atom_a->GetTargetAtom().GetPosition(),
                    atom_b.GetPosition(),
                    atom_c.GetPosition(),
                    itr_atom_d->GetTargetAtom().GetPosition()
                  )
                );
                double diff( math::Angle::Degree( math::Absolute( dihedral_abcd - *itr_moldih)));
                diff = std::min( diff, 360.0 - diff);
                if( diff > 10.0)
                {
                  found_something_legit = true;
                  BCL_MessageStd
                  (
                    "Chirality Change: " + util::Format()( diff)
                    + " on "
                    + util::Format()( fc.GetAtomIndex( itr_atom_a->GetTargetAtom())) + " - "
                    + util::Format()( fc.GetAtomIndex( atom_b)) + " - "
                    + util::Format()( fc.GetAtomIndex( atom_c)) + " - "
                    + util::Format()( fc.GetAtomIndex( itr_atom_d->GetTargetAtom()))
                    + " Central atom: " + util::Format()( ATOM_INDEX)
                  );
                }
                ++itr_moldih;
              }
            }
          }
        }

//      if( found_something_legit)
//      {
//        BCL_Debug( m_FragmentData->GetCentralAtomIndex());
//        BCL_Debug( m_ConnectedAtoms.GetSize());
//        BCL_Debug( m_FragmentData->GetAttachedAtomIndicesInverse().GetSize());
//        BCL_Debug( ( linal::Vector< double>( new_dihedral) - linal::Vector< double>( saved_dihedral)).Norm());
//        BCL_Debug( ( linal::Vector< double>( new_dihedral) - linal::Vector< double>( saved_dihedral)));
//        BCL_Debug( m_FragmentData->GetStandardizedCoordinateMatrix( MOLECULE));
//        BCL_Debug( m_FragmentData->GetStandardizedCoordinateMatrix( *fc));
//        BCL_Debug( chosen_matrix);
//        BCL_Debug( m_FragmentData->GetFreeEnergy( MOLECULE));
//        BCL_Debug( m_FragmentData->GetFreeEnergy( *fc));
//        BCL_Debug( m_FragmentData->GetBondAngleFreeEnergies()( angle_to_use));
//        BCL_Debug( m_FragmentData->GetAttachedAtomIndicesInverse());
//      }
      }
#endif
      return fc;
    }

    //! @brief Change the chiral, pseudochiral, or isometric-bond-adjacent identity of a given atom
    //! @param MOLECULE the molecule of interest
    //! @param ATOM_INDEX the atom whose chirality or E/Z bond identity should be updated
    //! @param ALTER_CHIRALITY_LABEL whether to swap/update the chirality or isometry label
    bool MutateChirality::MutateSpecificChirality
    (
      AtomVector< AtomComplete> &MOLECULE,
      const size_t &ATOM_INDEX,
      const bool &ALTER_CHIRALITY_LABEL,
      const bool &FORCE_ALTER_LABEL
    ) const
    {
      const storage::Vector< size_t> &neighbor_ids( m_Graph.GetNeighborIndices( ATOM_INDEX));

      // test for whether we are just flipping about a double bond, which is relatively trivial
      if( m_OnlyHasIsometry( ATOM_INDEX))
      {
        bool cannot_do_it( MOLECULE( ATOM_INDEX).CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, 1));
        if( cannot_do_it && !FORCE_ALTER_LABEL)
        {
          BCL_MessageCrt( "Cannot change ring isometry on " + util::Format()( ATOM_INDEX));
          return false;
        }
        // find the double-bonded neighbor, rotate it's dihedral by 180 degrees. Done
        for
        (
          storage::Vector< size_t>::const_iterator itr( neighbor_ids.Begin()), itr_end( neighbor_ids.End());
            itr != itr_end;
          ++itr
        )
        {
          if( m_GraphBondOrder.GetEdgeData( ATOM_INDEX, *itr) == 2)
          {
            {
              MutateBondAngles::RotateBond
              (
                MOLECULE,
                m_BondsToDownstreamAtoms( ATOM_INDEX).GetValue( *itr),
                math::g_Pi,
                ATOM_INDEX
              );
            }
            auto bond_type_original( MOLECULE( ATOM_INDEX).GetBondTypeTo( MOLECULE( *itr)));
            if( bond_type_original->GetIsometry() == e_EIsometry)
            {
              MOLECULE( ATOM_INDEX).SetBondTypeTo( MOLECULE( *itr), bond_type_original->WithIsometry( e_ZIsometry));
            }
            else if( bond_type_original->GetIsometry() == e_ZIsometry)
            {
              MOLECULE( ATOM_INDEX).SetBondTypeTo( MOLECULE( *itr), bond_type_original->WithIsometry( e_EIsometry));
            }
            else if( bond_type_original->GetIsometry() == e_UnknownIsometry)
            {
              BondIsometryHandler::AddIsometryInformation( MOLECULE, false);
            }
            return true;
          }
        }
        BCL_Exit( "Atom or bond order changed prior to mutate for isometry", -1);
      }

      storage::Vector< size_t> not_in_rings;
      storage::Vector< size_t> invariant;
      for
      (
        storage::Vector< size_t>::const_iterator itr( neighbor_ids.Begin()), itr_end( neighbor_ids.End());
          itr != itr_end;
        ++itr
      )
      {
        if( m_Graph.GetEdgeData( ATOM_INDEX, *itr) == 1 || not_in_rings.GetSize() == size_t( 2))
        {
          invariant.PushBack( *itr);
        }
        else
        {
          not_in_rings.PushBack( *itr);
        }
      }
      if( not_in_rings.GetSize() >= size_t( 2) && invariant.GetSize() == size_t( 1))
      {
        invariant.PushBack( not_in_rings.LastElement());
        not_in_rings.PopBack();
      }

      // 3 atoms in ring
      if( invariant.GetSize() == size_t( 3))
      {
        FragmentComplete old_mol( MOLECULE, "");
        // position of chiral center
        const linal::Vector3D &chiral_center_pos( MOLECULE( ATOM_INDEX).GetPosition());

        // get a normal vector to the invariant vertices
        linal::Vector3D normal
        (
          linal::CrossProduct
          (
            MOLECULE( invariant( 0)).GetPosition(),
            MOLECULE( invariant( 1)).GetPosition(),
            MOLECULE( invariant( 2)).GetPosition()
          ).Normalized()
        );

        // do a householder-flip on the two variable (not-in-ring) atom positions

        // compute the rotation, which is just I-2vV^T.
        // to make this faster, just compute I-OuterProduct(Sqrt(2)*v,Sqrt(2)*v)
        linal::Matrix3x3< double> rot;
        rot( 0, 0) = rot( 1, 1) = rot( 2, 2) = 1.0;
        linal::Vector3D normalsq2( normal);
        normalsq2 *= math::Sqrt( 2.0);
        linal::AddOuterProductToMatrix( rot, normalsq2, -normalsq2);

        // in this case we want to hold the point on the plane made up of the invariant vectors contant, so translate
        // everything there first. Likewise, we need to get the closest point on the plane to the chiral atoms position

        // average the positions of the invariants to get a point on the plane
        linal::Vector3D ave_pos
        (
          (
            MOLECULE( invariant( 0)).GetPosition()
            + MOLECULE( invariant( 1)).GetPosition()
            + MOLECULE( invariant( 2)).GetPosition()
          ) / 3.0
        );

        // get the distance from the invariant ave position and the chiral atom position
        double signed_dist_to_plane( linal::ScalarProduct( normal, ( ave_pos - chiral_center_pos)));
        linal::Vector3D invariant_point( chiral_center_pos + signed_dist_to_plane * normal);

        math::TransformationMatrix3D householder_flip;
        householder_flip( -invariant_point);
        householder_flip( math::RotationMatrix3D( rot));
        householder_flip( invariant_point);
        linal::Vector3D moved( chiral_center_pos);
        moved.Transform( householder_flip);
        double signed_dist_to_plane2( linal::ScalarProduct( normal, ( ave_pos - moved)));

        if( not_in_rings.GetSize())
        {
          // this gives us the updated position of the given atom, however, we can't just continue using the same transformation
          // across the atoms downstream from this one because it'd flip all their chiralities and dihedrals as well. Rather,
          // we compute a chirality / dihedral preserving transformation from the old position to the flipped position
          math::TransformationMatrix3D transform_matrix_a
          (
            coord::LineSegment3D( invariant_point, moved),
            coord::LineSegment3D( invariant_point, chiral_center_pos)
          );

          for( size_t nir_pos( 0), nir_sz( not_in_rings.GetSize()); nir_pos < nir_sz; ++nir_pos)
          {
            const storage::Vector< size_t> &connected_vertices_a
            (
              m_BondsToDownstreamAtoms( ATOM_INDEX).GetValue( not_in_rings( nir_pos))
            );
            for
            (
              storage::Vector< size_t>::const_iterator itr_a( connected_vertices_a.Begin()), itr_a_end( connected_vertices_a.End());
              itr_a != itr_a_end;
              ++itr_a
            )
            {
              linal::Vector3D cur_coords( MOLECULE( *itr_a).GetPosition());
              MOLECULE( *itr_a).SetPosition( cur_coords.Transform( transform_matrix_a));
            }
          }
        }
        MOLECULE( ATOM_INDEX).SetPosition( moved);
      }
      else if( not_in_rings.GetSize() == 0)
      {
        if( !FORCE_ALTER_LABEL)
        {
          BCL_MessageVrb
          (
            "Cannot change chirality type if atom has 4 bonds in a ring. Atom: "
            + util::Format()( ATOM_INDEX) + " Invariants: " + util::Format()( invariant.GetSize())
          );
          return false;
        }
      }
      else
      {
        // position of chiral center
        const linal::Vector3D &chiral_center_pos( MOLECULE( ATOM_INDEX).GetPosition());

        // get a normal vector to the invariant vertices
        linal::Vector3D normal
        (
          linal::CrossProduct
          (
            MOLECULE( ATOM_INDEX).GetPosition(),
            MOLECULE( invariant( 0)).GetPosition(),
            MOLECULE( invariant( 1)).GetPosition()
          ).Normalized()
        );

        // do a householder-flip on the two variable (not-in-ring) atom positions

        // compute the rotation, which is just I-2vV^T.
        // to make this faster, just compute I-OuterProduct(Sqrt(2)*v,Sqrt(2)*v)
        linal::Matrix3x3< double> rot;
        rot( 0, 0) = rot( 1, 1) = rot( 2, 2) = 1.0;
        normal *= math::Sqrt( 2.0);
        linal::AddOuterProductToMatrix( rot, normal, -normal);

        math::TransformationMatrix3D householder_flip;
        householder_flip( -chiral_center_pos);
        householder_flip( math::RotationMatrix3D( rot));
        householder_flip( chiral_center_pos);

        for( size_t nir_pos( 0), nir_sz( not_in_rings.GetSize()); nir_pos < nir_sz; ++nir_pos)
        {
          linal::Vector3D pos_a( MOLECULE( not_in_rings( nir_pos)).GetPosition());
          linal::Vector3D flipped_pos_a( pos_a);
          flipped_pos_a.Transform( householder_flip);

          // this gives us the updated position of the given atom, however, we can't just continue using the same transformation
          // across the atoms downstream from this one because it'd flip all their chiralities and dihedrals as well. Rather,
          // we compute a chirality / dihedral preserving transformation from the old position to the flipped position
          math::TransformationMatrix3D transform_matrix_a
          (
            coord::LineSegment3D( chiral_center_pos, flipped_pos_a),
            coord::LineSegment3D( chiral_center_pos, pos_a)
          );

          const storage::Vector< size_t> &connected_vertices_a
          (
            m_BondsToDownstreamAtoms( ATOM_INDEX).GetValue( not_in_rings( nir_pos))
          );
          for
          (
            storage::Vector< size_t>::const_iterator itr_a( connected_vertices_a.Begin()), itr_a_end( connected_vertices_a.End());
            itr_a != itr_a_end;
            ++itr_a
          )
          {
            linal::Vector3D cur_coords( MOLECULE( *itr_a).GetPosition());
            MOLECULE( *itr_a).SetPosition( cur_coords.Transform( transform_matrix_a));
          }
        }
      }

      if( ALTER_CHIRALITY_LABEL)
      {
        if( MOLECULE( ATOM_INDEX).GetChirality() == e_RChirality)
        {
          MOLECULE( ATOM_INDEX).SetChirality( e_SChirality);
        }
        else if( MOLECULE( ATOM_INDEX).GetChirality() == e_SChirality)
        {
          MOLECULE( ATOM_INDEX).SetChirality( e_RChirality);
        }
        else if( MOLECULE( ATOM_INDEX).GetChirality() == e_CisRingChirality)
        {
          MOLECULE( ATOM_INDEX).SetChirality( e_TransRingChirality);
        }
        else if( MOLECULE( ATOM_INDEX).GetChirality() == e_TransRingChirality)
        {
          MOLECULE( ATOM_INDEX).SetChirality( e_CisRingChirality);
        }
        else
        {
          // add stereocenter information
          StereocentersHandler::UpdateChiralityFromConformation( MOLECULE( ATOM_INDEX));
        }
      }
      return true;
    }

    //! @brief Copy the chirality of a different molecule onto this molecule
    //! This function assumes that MOLECULE currently has up-to-date chirality information
    //! Atoms with ambiguous or unknown chirality from TEMPLATE are ignored
    //! @return true if chirality was applied successfully
    bool MutateChirality::ApplyChirality
    (
      AtomVector< AtomComplete> &MOLECULE,
      const AtomVector< AtomComplete> &TEMPLATE,
      const bool &FORCE_ALTER_LABEL
    ) const
    {
      // add the chirality info to each atom with configuration
      BCL_Assert( TEMPLATE.GetSize() == MOLECULE.GetSize(), "Cannot apply chirality to molecules with different size");
      auto itr_template( TEMPLATE.Begin());
      size_t atom_index( 0);
      bool all_good( true);
      for
      (
        auto atom_itr( MOLECULE.Begin()), atom_itr_end( MOLECULE.End());
        atom_itr != atom_itr_end;
        ++atom_itr, ++itr_template, ++atom_index
      )
      {
        if
        (
          itr_template->GetChirality() == atom_itr->GetChirality()
          || itr_template->GetChirality() == e_UnknownChirality
          || itr_template->GetChirality() == e_UnknownRingChirality
          || itr_template->GetChirality() == e_NonChiral
        )
        {
          continue;
        }
        if( !all_good && !FORCE_ALTER_LABEL)
        {
          break;
        }
        all_good = MutateSpecificChirality( MOLECULE, atom_index, true, FORCE_ALTER_LABEL) && all_good;
      }
      return all_good;
    }

    //! @brief Copy the double bond isometries of a different molecule onto this molecule
    //! This function assumes that MOLECULE currently has up-to-date double bond isometry information
    //! @return true if isometry was applied successfully
    bool MutateChirality::ApplyDoubleBondIsometry
    (
      AtomVector< AtomComplete> &MOLECULE,
      const AtomVector< AtomComplete> &TEMPLATE,
      const bool &FORCE_ALTER_LABEL
    ) const
    {
      // add the chirality info to each atom with configuration
      BCL_Assert( TEMPLATE.GetSize() == MOLECULE.GetSize(), "Cannot apply chirality to molecules with different size");
      auto itr_template( TEMPLATE.Begin());
      size_t atom_index( 0);
      bool all_good( true);
      for
      (
        auto atom_itr( MOLECULE.Begin()), atom_itr_end( MOLECULE.End());
        atom_itr != atom_itr_end;
        ++atom_itr, ++itr_template, ++atom_index
      )
      {
        if( m_OnlyHasIsometry( atom_index))
        {
          auto itr_bnd( atom_itr->GetBonds().Begin());
          for
          (
            auto itr_bnd_end( atom_itr->GetBonds().End());
            itr_bnd != itr_bnd_end;
            ++itr_bnd
          )
          {
            if
            (
              itr_bnd->GetBondType()->GetNumberOfElectrons() == size_t( 4)
              && itr_bnd->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_IsIsometric)
            )
            {
              break;
            }
          }
          if( itr_bnd != atom_itr->GetBonds().End())
          {
            size_t bnd_atom_idx( MOLECULE.GetAtomIndex( itr_bnd->GetTargetAtom()));
            if( itr_bnd->GetBondType() != itr_template->GetBondTypeTo( TEMPLATE( bnd_atom_idx)))
            {
              all_good = MutateSpecificChirality( MOLECULE, atom_index, true, FORCE_ALTER_LABEL) && all_good;
              if( !FORCE_ALTER_LABEL && !all_good)
              {
                break;
              }
            }
          }
        }
      }
      return all_good;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateChirality::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ChiralCenters, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateChirality::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ChiralCenters, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
