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
#include "chemistry/bcl_chemistry_mutate_bond_angles.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_dihedral_angles.h"
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_conformation_comparison_by_dihedral_bins.h"
#include "chemistry/bcl_chemistry_constitution_graph_converter.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_ring_fragment_map.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "chemistry/bcl_chemistry_valence_handler.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_histogram.h"
#include "quality/bcl_quality_rmsd.h"

using bcl::chemistry::HydrogensHandler;

// external includes - sorted alphabetically

//#define DEBUG_CHIRALITY_MUTATIONS
#define BCL_PROFILE_MutateBondAngles
#ifdef BCL_PROFILE_MutateBondAngles
#include "util/bcl_util_stopwatch.h"
#endif

namespace bcl
{
  namespace chemistry
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> MutateBondAngles::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateBondAngles())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateBondAngles::MutateBondAngles() :
      m_FragmentData(),
      m_ProbabilityDistribution(),
      m_ConnectedAtoms(),
      m_Unbiased( false)
    {
    }

    //! @brief constructor taking the member variable parameters
    //! @param FRAGMENT fragment which the mutate mutates
    //! @param MOLECULE molecule of interest
    MutateBondAngles::MutateBondAngles
    (
      const util::ShPtr< BondAngleAssignment> &FRAGMENT,
      const FragmentComplete &MOLECULE,
      const MutateChirality &MUTATE_CHIRALITY,
      bool CHANGE_CHIRALITY,
      bool UNBIASED,
      bool SAMPLE_LENGTHS
    ) :
      m_FragmentData( FRAGMENT),
      m_ChangeChirality( CHANGE_CHIRALITY),
      m_Unbiased( UNBIASED),
      m_MutateChirality( MUTATE_CHIRALITY),
      m_SampleBondLengths( SAMPLE_LENGTHS)
    {
      // for each of the non-ring bonds get atoms connected to it. Iterate over list of non-ring bonds
      const AtomConformationalInterface &atom( MOLECULE.GetAtomVector()( FRAGMENT->GetCentralAtomIndex()));
      const storage::Vector< BondConformational> &bonds( atom.GetBonds());
      m_NumberRingBonds = atom.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1));
      size_t bond_index( 0);
      for( auto itr( bonds.Begin()), itr_end( bonds.End()); itr != itr_end; ++itr, ++bond_index)
      {
        // for each non rotatable bond get connected atoms
        m_ConnectedAtoms.PushBack();
        if( !itr->GetBondType()->IsBondInRing())
        {
          m_ConnectedAtoms.LastElement()
              = m_MutateChirality->GetDownstreamAtoms( FRAGMENT->GetCentralAtomIndex(), MOLECULE.GetAtomIndex( itr->GetTargetAtom()));
        }
      }

      // get counts for all the rotamers to construct a probability histogram
      storage::Map< double, double> histogram_vector;

      // get counts for all the rotamers to construct a probability histogram
      size_t rotamer_index( 0);
      for
      (
        auto itr_count( m_FragmentData->GetBondAngleCounts().Begin()),
             itr_count_end( m_FragmentData->GetBondAngleCounts().End());
        itr_count != itr_count_end;
        ++itr_count, ++rotamer_index
      )
      {
        if( !m_Unbiased)
        {
          histogram_vector[ double( rotamer_index)] = itr_count->Second().GetWeight() + m_FragmentData->GetPseudocount();
        }
        else
        {
          histogram_vector[ double( rotamer_index)] = 1.0;
        }
      }
      m_ProbabilityDistribution = util::ShPtr< random::Histogram1DDistribution>
      (
        new random::Histogram1DDistribution( math::Histogram( histogram_vector))
      );
      // get complement of rings that exist in different conformations
      m_Chirality = atom.GetChirality();

      iterate::Generic< const AtomConformationalInterface> itr_atom_b( MOLECULE.GetAtomsIterator());
      itr_atom_b.GotoPosition( m_FragmentData->GetCentralAtomIndex());
      const AtomConformationalInterface &atom_b( *itr_atom_b);
      const storage::Vector< BondConformational> &connected_atoms_b( itr_atom_b->GetBonds());
      for
      (
        storage::Vector< BondConformational>::const_iterator
          itr_atom_c( connected_atoms_b.Begin()), itr_atom_c_end( connected_atoms_b.End());
        itr_atom_c != itr_atom_c_end;
        ++itr_atom_c
      )
      {
        // skip ring bonds
        if( itr_atom_c->GetBondType()->IsBondInRing())
        {
          continue;
        }
        const AtomConformationalInterface &atom_c( itr_atom_c->GetTargetAtom());

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
            m_DihedralAngleIndices.PushBack
            (
              storage::VectorND< 4, size_t>
              (
                MOLECULE.GetAtomIndex( itr_atom_a->GetTargetAtom()),
                m_FragmentData->GetCentralAtomIndex(),
                MOLECULE.GetAtomIndex( itr_atom_c->GetTargetAtom()),
                MOLECULE.GetAtomIndex( itr_atom_d->GetTargetAtom())
              )
            );
          }
        }
      }

    }

    //! @brief Clone function
    //! @return pointer to new MutateBondAngles
    MutateBondAngles *MutateBondAngles::Clone() const
    {
      return new MutateBondAngles( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateBondAngles::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateBondAngles::GetScheme() const
    {
      static const std::string s_gen( "MutateBondAngle");
      return s_gen;
    }

    //! @brief returns the fragment that the mutate uses
    //! @return the fragment that the mutate uses
    const util::ShPtr< BondAngleAssignment> &MutateBondAngles::GetFragmentData() const
    {
      return m_FragmentData;
    }

    //! @brief returns probability distribution for rotamers of this particular fragments
    //! @return probability distribution for rotamers of this particular fragments
    const random::Histogram1DDistribution &MutateBondAngles::GetProbability() const
    {
      return *m_ProbabilityDistribution;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking a conformation and returns a mutated conformation
    //! @param MOLECULE conformation of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< FragmentComplete> MutateBondAngles::Mutate
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      return this->Mutate( MOLECULE, false);
    }

    //! @brief operator taking a conformation and returns a mutated conformation
    //! @param MOLECULE conformation of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< FragmentComplete> MutateBondAngles::Mutate
    (
      const FragmentComplete &MOLECULE,
      bool SWAP_Z
    ) const
    {
#ifdef BCL_PROFILE_MutateBondAngles
      static util::Stopwatch s_timer_2bond( "bond_angle_sampling_for_2_bonded_atoms", util::Message::e_Standard, true, false);
      static util::Stopwatch s_timer_34bond( "bond_angle_sampling_for_3-4_bonded_atoms", util::Message::e_Standard, true, false);
      static util::Stopwatch s_timer_ringbond( "bond_angles_sampling_for_ring_bonds", util::Message::e_Standard, true, false);
      static util::Stopwatch s_timer_bondanglesetup( "bond_angle_setup", util::Message::e_Standard, true, false);
      static util::Stopwatch s_timer_bondlength( "bond_length_sampling", util::Message::e_Standard, true, false);
      static util::Stopwatch s_timer_flip( "flipping_invertable_nitrogen_centers", util::Message::e_Standard, true, false);
      s_timer_bondanglesetup.Start();
#endif

      //MOLECULE.WriteMDL( util::GetLogger());
      size_t angle_to_use
      (
        GetChooseBestOnlyFlag()
        ? m_ProbabilityDistribution->DetermineMostLikelyCase()
        : m_ProbabilityDistribution->DetermineRandomCase()
      );
      auto itr_angle( m_FragmentData->GetBondAngleCounts().Begin());
      std::advance( itr_angle, angle_to_use);
      const linal::Matrix< double> &chosen_matrix( itr_angle->First());

      AtomVector< AtomComplete> new_frag( MOLECULE.GetAtomVector());

      // sample the bond lengths occasionally
      if( !SWAP_Z)
      {
        #ifdef BCL_PROFILE_MutateBondAngles
        s_timer_bondanglesetup.Stop();
        s_timer_bondlength.Start();
        #endif
        AtomComplete &atom( new_frag( m_FragmentData->GetCentralAtomIndex()));
        size_t n_bonds( m_FragmentData->GetAttachedAtomIndicesInverse().GetSize());
        size_t n_nonringindex( 0);
        for( size_t bondn( 0); bondn < n_bonds; ++bondn)
        {
          const size_t x( m_FragmentData->GetAttachedAtomIndicesInverse()( bondn));
          if( atom.GetBonds()( x).GetBondType()->IsBondInRing())
          {
            continue;
          }
          double desired_length( itr_angle->Second().GetAverage()( n_nonringindex));
          ++n_nonringindex;
          if( !m_FragmentData->GetCanUseBondLengths())
          {
            desired_length =
                BondLengths::GetBondLength
                (
                  atom.GetAtomType(),
                  atom.GetBonds()( x).GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderOrAromatic),
                  atom.GetBonds()( x).GetTargetAtom().GetAtomType()
                );
          }
          if( !util::IsDefined( desired_length))
          {
            continue;
          }
          if( m_SampleBondLengths)
          {
            desired_length += std::max( -0.08, std::min( 0.08, random::GetGlobalRandom().RandomGaussian( 0.0, 0.03)));
          }
          const AtomConformationalInterface &atomx( atom.GetBonds()( x).GetTargetAtom());
          linal::Vector3D existing_bond_direction( atomx.GetPosition() - atom.GetPosition());
          const double norm( existing_bond_direction.Norm());
          if( norm)
          {
            if( math::Absolute( norm - desired_length) < 0.01)
            {
              continue;
            }
            existing_bond_direction *= desired_length / norm - 1.0;
          }
          else
          {
            existing_bond_direction = linal::Vector3D( desired_length, 0.0, 0.0);
          }
          // existing bond direction is now the translation to add to everything on the other side of the bond
          const storage::Vector< size_t> &vectices_reachable
          (
            m_MutateChirality->GetDownstreamAtoms( m_FragmentData->GetCentralAtomIndex(), new_frag.GetAtomIndex( atomx))
          );
          for
          (
            auto itr_downstream( vectices_reachable.Begin()), itr_downstream_end( vectices_reachable.End());
            itr_downstream != itr_downstream_end;
            ++itr_downstream
          )
          {
            new_frag( *itr_downstream).SetPosition( new_frag( *itr_downstream).GetPosition() + existing_bond_direction);
          }
        }
        #ifdef BCL_PROFILE_MutateBondAngles
        s_timer_bondlength.Stop();
        s_timer_bondanglesetup.Start();
        #endif
      }
      bool did_flip( false);
      FragmentComplete new_frag_saved;
      if( m_FragmentData->CanFlip() && !SWAP_Z && random::GetGlobalRandom().Boolean())
      {
        #ifdef BCL_PROFILE_MutateBondAngles
        s_timer_bondanglesetup.Stop();
        s_timer_flip.Start();
        #endif
        m_MutateChirality->MutateSpecificChirality
        (
          new_frag,
          m_FragmentData->GetCentralAtomIndex(),
          true
        );
        new_frag_saved = FragmentComplete( new_frag, MOLECULE.GetName());
        did_flip = true;
        #ifdef BCL_PROFILE_MutateBondAngles
        s_timer_flip.Stop();
        s_timer_bondanglesetup.Start();
        #endif
        // BCL_MessageStd( "Flipped " + util::Format()( m_FragmentData->GetCentralAtomIndex()));
      }
      AtomComplete &atom( new_frag( m_FragmentData->GetCentralAtomIndex()));

      #ifdef BCL_PROFILE_MutateBondAngles
      s_timer_bondanglesetup.Stop();
      #endif
      if( m_FragmentData->GetAttachedAtomIndicesInverse().GetSize() == size_t( 2))
      {
        #ifdef BCL_PROFILE_MutateBondAngles
        s_timer_2bond.Start();
        #endif
        const size_t atom_a_bond_index( m_FragmentData->GetAttachedAtomIndicesInverse()( 0));
        const size_t atom_b_bond_index( m_FragmentData->GetAttachedAtomIndicesInverse()( 1));
        const AtomConformationalInterface &atom_a( atom.GetBonds()( atom_a_bond_index).GetTargetAtom());
        const AtomConformationalInterface &atom_b( atom.GetBonds()( atom_b_bond_index).GetTargetAtom());
        const double current_angle
        (
          linal::ProjAngle( atom.GetPosition(), atom_a.GetPosition(), atom_b.GetPosition())
        );
        double chosen_angleb( math::Sign( std::acos( chosen_matrix( 1, 0)), current_angle));

        linal::Vector3D cross( linal::CrossProduct( atom.GetPosition(), atom_a.GetPosition(), atom_b.GetPosition()));
        math::TransformationMatrix3D transb
        (
          coord::LineSegment3D( atom.GetPosition(), atom.GetPosition() + cross),
          atom.GetPosition(),
          current_angle - chosen_angleb
        );

        for( auto itratom( m_ConnectedAtoms( atom_b_bond_index).Begin()), itratom_end( m_ConnectedAtoms( atom_b_bond_index).End()); itratom != itratom_end; ++itratom)
        {
          linal::Vector3D npos( new_frag( *itratom).GetPosition());
          npos.Transform( transb);
          new_frag( *itratom).SetPosition( npos);
        }
        #ifdef BCL_PROFILE_MutateBondAngles
        s_timer_2bond.Stop();
        #endif
      }
      else
      {
        linal::Vector3D translation( -atom.GetPosition());
        for( auto itratom( new_frag.Begin()), itratom_end( new_frag.End()); itratom != itratom_end; ++itratom)
        {
          linal::Vector3D poscopy( itratom->GetPosition());
          poscopy.Translate( translation);
          itratom->SetPosition( poscopy);
        }

        static const linal::Vector3D origin( 0.0, 0.0, 0.0), x_axis( 1.0, 0.0, 0.0);
        if( !m_NumberRingBonds)
        {
          #ifdef BCL_PROFILE_MutateBondAngles
          s_timer_34bond.Start();
          #endif
          const AtomConformationalInterface &atom_b( atom.GetBonds()( m_FragmentData->GetAttachedAtomIndicesInverse()( 0)).GetTargetAtom());
          math::TransformationMatrix3D align_bond_a
          (
            coord::LineSegment3D( origin, x_axis),
            coord::LineSegment3D( origin, linal::UnitVector( origin, atom_b.GetPosition()))
          );
          for( auto itratom( new_frag.Begin()), itratom_end( new_frag.End()); itratom != itratom_end; ++itratom)
          {
            linal::Vector3D poscopy( itratom->GetPosition());
            poscopy.Transform( align_bond_a);
            itratom->SetPosition( poscopy);
          }
          double original_dihedral_abcd( 0.0);
          if( atom_b.GetBonds().GetSize() > size_t( 1))
          {
            original_dihedral_abcd =
              linal::Dihedral
              (
                new_frag( m_ConnectedAtoms( m_FragmentData->GetAttachedAtomIndicesInverse()( 0))( 1)).GetPosition(),
                x_axis,
                origin,
                new_frag( m_ConnectedAtoms( m_FragmentData->GetAttachedAtomIndicesInverse()( 1))( 0)).GetPosition()
              );
          }

          size_t n_bonds( m_FragmentData->GetAttachedAtomIndicesInverse().GetSize());
          for( size_t bondn( 1); bondn < n_bonds; ++bondn)
          {
            const size_t x( m_FragmentData->GetAttachedAtomIndicesInverse()( bondn));
            const AtomConformationalInterface &atomx( atom.GetBonds()( x).GetTargetAtom());
            linal::Vector3D atomx_pos( atomx.GetPosition());
            double old_dihedral_val( 0.0);
            if( atomx.GetBonds().GetSize() > size_t( 1))
            {
              old_dihedral_val =
                linal::Dihedral
                (
                  x_axis,
                  origin,
                  new_frag( m_ConnectedAtoms( x)( 0)).GetPosition(),
                  new_frag( m_ConnectedAtoms( x)( 1)).GetPosition()
                );
            }
            atomx_pos.Normalize();
            math::TransformationMatrix3D align_bond_x
            (
              coord::LineSegment3D( origin, linal::Vector3D( chosen_matrix.GetRow( bondn))),
              coord::LineSegment3D( origin, atomx_pos)
            );
            linal::Vector3D c( atomx_pos);
            c.Transform( align_bond_x);

            for
            (
              auto itr_con( m_ConnectedAtoms( x).Begin()), itr_con_end( m_ConnectedAtoms( x).End());
              itr_con != itr_con_end;
              ++itr_con
            )
            {
              linal::Vector3D poscopy( new_frag( *itr_con).GetPosition());
              poscopy.Transform( align_bond_x);
              new_frag( *itr_con).SetPosition( poscopy);
            }

            if( atomx.GetBonds().GetSize() > size_t( 1))
            {
              double new_dihedral_val
              (
                linal::Dihedral
                (
                  x_axis,
                  origin,
                  new_frag( m_ConnectedAtoms( x)( 0)).GetPosition(),
                  new_frag( m_ConnectedAtoms( x)( 1)).GetPosition()
                )
              );
              RotateBond
              (
                new_frag,
                m_ConnectedAtoms( x),
                new_dihedral_val - old_dihedral_val,
                m_FragmentData->GetCentralAtomIndex()
              );
            }
          }
          if( atom_b.GetBonds().GetSize() > size_t( 1))
          {
            double new_dihedral_abcd =
              linal::Dihedral
              (
                new_frag( m_ConnectedAtoms( m_FragmentData->GetAttachedAtomIndicesInverse()( 0))( 1)).GetPosition(),
                x_axis,
                origin,
                new_frag( m_ConnectedAtoms( m_FragmentData->GetAttachedAtomIndicesInverse()( 1))( 0)).GetPosition()
              );
            RotateBond
            (
              new_frag,
              m_ConnectedAtoms( m_FragmentData->GetAttachedAtomIndicesInverse()( 0)),
              new_dihedral_abcd - original_dihedral_abcd,
              m_FragmentData->GetCentralAtomIndex()
            );
          }
          #ifdef BCL_PROFILE_MutateBondAngles
          s_timer_34bond.Stop();
          #endif
        }
        else if( m_NumberRingBonds >= size_t( 2))
        {
          #ifdef BCL_PROFILE_MutateBondAngles
          s_timer_ringbond.Start();
          #endif

          // get the ring bond coordinates
          storage::Vector< linal::Vector3D> ring_coords_template( size_t( m_NumberRingBonds + 1), linal::Vector3D( 0.0));
          storage::Vector< linal::Vector3D> ring_coords_current( size_t( m_NumberRingBonds + 1), linal::Vector3D( 0.0));
          ring_coords_template( 1) = linal::Vector3D( 1.0, 0.0, 0.0);
          ring_coords_template( 2) = linal::Vector3D( 0.0, 1.0, 0.0);
          ring_coords_current( 1) = atom.GetBonds()( m_FragmentData->GetAttachedAtomIndicesInverse()( 0)).GetTargetAtom().GetPosition().Normalized();
          ring_coords_current( 2) = atom.GetBonds()( m_FragmentData->GetAttachedAtomIndicesInverse()( 1)).GetTargetAtom().GetPosition().Normalized();
          if( m_NumberRingBonds == size_t( 3))
          {
            ring_coords_template( 3) = linal::Vector3D( 0.0, 0.0, 1.0);
            ring_coords_current( 3) = atom.GetBonds()( m_FragmentData->GetAttachedAtomIndicesInverse()( 2)).GetTargetAtom().GetPosition().Normalized();
          }
          math::TransformationMatrix3D transform_ring_coords
          (
            quality::RMSD::SuperimposeCoordinates
            (
              util::SiPtrVector< const linal::Vector3D>( ring_coords_template.Begin(), ring_coords_template.End()),
              util::SiPtrVector< const linal::Vector3D>( ring_coords_current.Begin(), ring_coords_current.End())
            )
          );
          transform_ring_coords.SetTranslation( linal::Vector3D( 0.0));

          // align the ring bonds, but maintain the origin
          for( auto itratom( new_frag.Begin()), itratom_end( new_frag.End()); itratom != itratom_end; ++itratom)
          {
            linal::Vector3D poscopy( itratom->GetPosition());
            poscopy.Transform( transform_ring_coords);
            itratom->SetPosition( poscopy);
          }

          size_t n_bonds( m_FragmentData->GetAttachedAtomIndicesInverse().GetSize());
          size_t matrix_row( 0);
          for( size_t bondn( 0); bondn < n_bonds; ++bondn)
          {
            const size_t x( m_FragmentData->GetAttachedAtomIndicesInverse()( bondn));

            if( atom.GetBonds()( x).GetBondType()->IsBondInRing())
            {
              continue;
            }
            const AtomConformationalInterface &atomx( atom.GetBonds()( x).GetTargetAtom());
            linal::Vector3D atomx_pos( atomx.GetPosition());
            double old_dihedral_val( 0.0);
            if( atomx.GetBonds().GetSize() > size_t( 1))
            {
              old_dihedral_val =
                linal::Dihedral
                (
                  atom.GetBonds()( m_FragmentData->GetAttachedAtomIndicesInverse()( 0)).GetTargetAtom().GetPosition(),
                  origin,
                  new_frag( m_ConnectedAtoms( x)( 0)).GetPosition(),
                  new_frag( m_ConnectedAtoms( x)( 1)).GetPosition()
                );
            }

            atomx_pos.Normalize();
            linal::Vector3D nrow( chosen_matrix.GetRow( matrix_row));
//            if( nrow.Z() * atomx_pos.Z() < 0.0)
//            {
//              nrow.Z() *= -1.0;
//            }
            if( SWAP_Z)
            {
              nrow.Z() *= -1.0;
            }
            math::TransformationMatrix3D align_bond_x
            (
              coord::LineSegment3D( origin, nrow),
              coord::LineSegment3D( origin, atomx_pos)
            );
            linal::Vector3D c( atomx_pos);
            c.Transform( align_bond_x);

            for
            (
              auto itr_con( m_ConnectedAtoms( x).Begin()), itr_con_end( m_ConnectedAtoms( x).End());
              itr_con != itr_con_end;
              ++itr_con
            )
            {
              linal::Vector3D poscopy( new_frag( *itr_con).GetPosition());
              poscopy.Transform( align_bond_x);
              new_frag( *itr_con).SetPosition( poscopy);
            }

            if( atomx.GetBonds().GetSize() > size_t( 1))
            {
              double new_dihedral_val
              (
                linal::Dihedral
                (
                  atom.GetBonds()( m_FragmentData->GetAttachedAtomIndicesInverse()( 0)).GetTargetAtom().GetPosition(),
                  origin,
                  new_frag( m_ConnectedAtoms( x)( 0)).GetPosition(),
                  new_frag( m_ConnectedAtoms( x)( 1)).GetPosition()
                )
              );
              RotateBond
              (
                new_frag,
                m_ConnectedAtoms( x),
                new_dihedral_val - old_dihedral_val,
                m_FragmentData->GetCentralAtomIndex()
              );
            }
            ++matrix_row;
          }
          #ifdef BCL_PROFILE_MutateBondAngles
          s_timer_ringbond.Stop();
          #endif
        }
      }
      // BCL_MessageStd( "Doing: " + util::Format()( m_FragmentData->GetCentralAtomIndex()));
      #ifdef BCL_PROFILE_MutateBondAngles
      s_timer_bondanglesetup.Start();
      #endif
      HydrogensHandler::UpdateHCoordinates( atom, new_frag);
      for( auto itr_atm( new_frag.Begin()), itr_atm_end( new_frag.End()); itr_atm != itr_atm_end; ++itr_atm)
      {
        // check that coordinates are defined
        if( !itr_atm->GetPosition().IsDefined())
        {
          #ifdef BCL_PROFILE_MutateBondAngles
          BCL_MessageStd( "MBA- bad angles");
          s_timer_bondanglesetup.Stop();
          #endif
          return math::MutateResult< FragmentComplete>();
        }
      }
      if( m_Chirality == e_RChirality || m_Chirality == e_SChirality || m_Chirality == e_UnknownChirality)
      {
        StereocentersHandler::UpdateChiralityFromConformation( new_frag( m_FragmentData->GetCentralAtomIndex()));
      }
      else if( m_Chirality == e_TransRingChirality || m_Chirality == e_CisRingChirality)
      {
        StereocentersHandler::AddChiralityFromConformation( new_frag);
      }

      util::ShPtr< FragmentComplete> fc( new FragmentComplete( new_frag, MOLECULE.GetName()));
      #ifdef BCL_PROFILE_MutateBondAngles
      s_timer_bondanglesetup.Stop();
      #endif
      if( m_Chirality != atom.GetChirality() && !m_ChangeChirality)
      {
        if( !SWAP_Z)
        {
          return Mutate( did_flip ? new_frag_saved : MOLECULE, true);
        }
        BCL_MessageStd( "Chirality messed up by bond angles " + util::Format()( m_FragmentData->GetCentralAtomIndex()) + " " + util::Format()( chosen_matrix));
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
      }

//      bool found_something_legit( true);
//      bool must_stop( false);
//      double fifteen_radian( math::Angle::Radian( 15.0));
//      double thirty_radian( math::Angle::Radian( 30.0));
//
//      const FragmentComplete &moldih( did_flip ? new_frag_saved : MOLECULE);
//      const AtomVector< AtomComplete> &moldihatms( moldih.GetAtomVector());
//      size_t n_cycles( 0);
//      while( found_something_legit && ++n_cycles < size_t( 3))
//      {
//        found_something_legit = false;
//        for
//        (
//          auto itr_di( m_DihedralAngleIndices.Begin()), itr_di_end( m_DihedralAngleIndices.End());
//          itr_di != itr_di_end;
//          ++itr_di
//        )
//        {
//          double dihedral_abcd
//          (
//            linal::Dihedral
//            (
//              new_frag( itr_di->First()).GetPosition(),
//              new_frag( itr_di->Second()).GetPosition(),
//              new_frag( itr_di->Third()).GetPosition(),
//              new_frag( itr_di->Fourth()).GetPosition()
//            ) + fifteen_radian
//          );
//
//          double dihedral_mol_abcd
//          (
//            linal::Dihedral
//            (
//              moldihatms( itr_di->First()).GetPosition(),
//              moldihatms( itr_di->Second()).GetPosition(),
//              moldihatms( itr_di->Third()).GetPosition(),
//              moldihatms( itr_di->Fourth()).GetPosition()
//            ) + fifteen_radian
//          );
//          if( rint( dihedral_abcd / thirty_radian) != rint( dihedral_mol_abcd / thirty_radian))
//          {
//            RotateBond
//            (
//              new_frag,
//              m_MutateChirality->GetDownstreamAtoms( itr_di->Second(), itr_di->Third()),
//              ( dihedral_abcd - dihedral_mol_abcd),
//              itr_di->Second()
//            );
//            found_something_legit = true;
//
//            if( util::GetMessenger().GetCurrentMessageLevel() == util::Message::e_Debug)
//            {
//              double diff( math::Angle::Degree( math::Absolute( dihedral_abcd - dihedral_mol_abcd)));
//              diff = std::min( diff, 360.0 - diff);
//              BCL_MessageStd
//              (
//                "Change: " + util::Format()( diff)
//                + " " + util::Format()( rint( dihedral_abcd / math::Angle::Radian( 30.0)))
//                + " " + util::Format()( rint( dihedral_mol_abcd / math::Angle::Radian( 30.0)))
//                + " on "
//                + util::Format()( itr_di->First()) + " - "
//                + util::Format()( itr_di->Second()) + " - "
//                + util::Format()( itr_di->Third()) + " - "
//                + util::Format()( itr_di->Fourth())
//                + " Central atom: " + util::Format()( m_FragmentData->GetCentralAtomIndex())
//              );
//            }
//          }
//        }
//      }

      if( m_FragmentData->HasAnAmideBond() && fc->GetTotalAmideBondNonPlanarity( 20.0) > MOLECULE.GetTotalAmideBondNonPlanarity( 20.0) + 0.01)
      {
        BCL_MessageDbg( "Perturbed amide bond");
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
      }

      if( util::GetMessenger().GetCurrentMessageLevel() == util::Message::e_Debug) // this is a debug section
      {
        double expected_energy( m_FragmentData->GetBondAngleFreeEnergies()( angle_to_use));
        double actual_energy( m_FragmentData->GetFreeEnergy( *fc));
        if( expected_energy != actual_energy && !SWAP_Z)
        {
          BCL_MessageStd
          (
            "Messed up energies..." + util::Format()( expected_energy) + " " + util::Format()( actual_energy)
            + " " + util::Format()( angle_to_use) + " " + util::Format()( m_FragmentData->GetCentralAtomIndex())
            + " " + util::Format()( m_NumberRingBonds)
            + " " + util::Format()( m_FragmentData->GetBondAngleFreeEnergies())
          );
          BCL_Debug( m_FragmentData->GetStandardizedCoordinateMatrix( MOLECULE));
          BCL_Debug( m_FragmentData->GetStandardizedCoordinateMatrix( *fc));
          BCL_Debug( chosen_matrix);
        }
      }

//      if( found_something_legit)
//      {
//        if( util::GetMessenger().GetCurrentMessageLevel() == util::Message::e_Debug) // this is a debug section
//        {
//          BCL_Debug( m_FragmentData->GetCentralAtomIndex());
//          BCL_Debug( m_ConnectedAtoms.GetSize());
//          BCL_Debug( m_FragmentData->GetAttachedAtomIndicesInverse().GetSize());
//          BCL_Debug( m_FragmentData->GetStandardizedCoordinateMatrix( did_flip ? new_frag_saved : MOLECULE));
//          BCL_Debug( m_FragmentData->GetStandardizedCoordinateMatrix( *fc));
//          BCL_Debug( chosen_matrix);
//          BCL_Debug( m_FragmentData->GetFreeEnergy( did_flip ? new_frag_saved : MOLECULE));
//          BCL_Debug( m_FragmentData->GetFreeEnergy( *fc));
//          BCL_Debug( m_FragmentData->GetBondAngleFreeEnergies()( angle_to_use));
//          BCL_Debug( m_FragmentData->GetAttachedAtomIndicesInverse());
//        }
//        return math::MutateResult< FragmentComplete>();
//      }
//        io::OFStream output;
//        io::File::MustOpenOFStream( output,  "mba.before.sdf", std::ios::app);
//        aligned_frag.WriteMDL( output);
//        io::File::CloseClearFStream( output);
//        io::File::MustOpenOFStream( output,  "mba.after.sdf", std::ios::app);
//        fc->WriteMDL( output);
//        io::File::CloseClearFStream( output);

      return math::MutateResult< FragmentComplete>( fc, *this);
    }

    //! @brief rotate the molecule of interest at a particular bond at a particular angle
    //! @param ATOM_INFO the atoms that needs to be updated with coordinates after rotatating about a particular bond
    //! @param CONNECTED_ATOMS the connected atoms associated with bond of interest
    //! @param NEW_ANGLE the new angle for the bond of interest
    //! @param EXISTING_ANGLE the angle for the bond of interest in the argument that is passed in
    void MutateBondAngles::RotateBond
    (
      AtomVector< AtomComplete> &ATOM_INFO,
      const storage::Vector< size_t> &CONNECTED_ATOMS,
      double ANGULAR_CHANGE,
      const size_t &ORIGINATING_ATOM
    )
    {
      // get the anchor points for bond rotation, i.e. indices of atoms that make the bond of interest
      storage::Vector< size_t>::const_iterator itr_connected_atoms( CONNECTED_ATOMS.Begin());
      const linal::Vector3D &end_point( ATOM_INFO( *itr_connected_atoms).GetPosition());
      ++itr_connected_atoms;

      math::TransformationMatrix3D transformation
      (
        coord::LineSegment3D( ATOM_INFO( ORIGINATING_ATOM).GetPosition(), end_point),
        ATOM_INFO( ORIGINATING_ATOM).GetPosition(),
        ANGULAR_CHANGE < 0.0 ? ANGULAR_CHANGE + 2.0 * math::g_Pi : ANGULAR_CHANGE
      );
//      BCL_MessageStd
//      (
//        "Rotating bond beginning with: " + util::Format()( CONNECTED_ATOMS( 0)) + " by " + util::Format()( math::Angle::Degree( ANGULAR_CHANGE)) + " second atom if present " +
//        util::Format()( CONNECTED_ATOMS.GetSize() > size_t( 1) ? CONNECTED_ATOMS( 1) : CONNECTED_ATOMS( 0))
//        + " " + util::Format()( ATOM_INFO( ORIGINATING_ATOM).GetPosition())
//        + " " + util::Format()( end_point)
//      );

      // set the coordinates of the atoms from the transformed fragment cloud
      for
      (
        storage::Vector< size_t>::const_iterator itr_connected_atoms_end( CONNECTED_ATOMS.End());
        itr_connected_atoms != itr_connected_atoms_end;
        ++itr_connected_atoms
      )
      {
        linal::Vector3D point( ATOM_INFO( *itr_connected_atoms).GetPosition());
        ATOM_INFO( *itr_connected_atoms).SetPosition( point.Transform( transformation));
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateBondAngles::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_ConnectedAtoms, ISTREAM);
      io::Serialize::Read( m_ProbabilityDistribution, ISTREAM);

      // return "ISTREAM"
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &MutateBondAngles::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ConnectedAtoms, OSTREAM, INDENT);
      io::Serialize::Write( m_ProbabilityDistribution, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
