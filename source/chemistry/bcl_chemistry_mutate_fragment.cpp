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
#include "chemistry/bcl_chemistry_mutate_fragment.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "chemistry/bcl_chemistry_bond_dihedral_angles.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
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

// external includes - sorted alphabetically

//#define DEBUG_CHIRALITY_MUTATIONS
#define BCL_PROFILE_MutateFragment
#ifdef BCL_PROFILE_MutateFragment
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
    const util::SiPtr< const util::ObjectInterface> MutateFragment::s_Instance
    (
      GetObjectInstances().AddInstance( new MutateFragment())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MutateFragment::MutateFragment() :
      m_FragmentData(),
      m_ProbabilityDistribution(),
      m_ConnectedAtoms(),
      m_FragmentMap(),
      m_Unbiased( false),
      m_CheckAmideBondPlanarity( false)
    {
    }

    //! @brief constructor taking the member variable parameters
    //! @param FRAGMENT fragment which the mutate mutates
    //! @param GRAPH constitution graph of the molecule of interest
    //! @param MOLECULE molecule of interest
    MutateFragment::MutateFragment
    (
      const util::ShPtr< RotamerDihedralBondData> &FRAGMENT,
      const graph::ConstGraph< size_t, size_t> &GRAPH,
      const FragmentComplete &MOLECULE,
      const MutateChirality &MUTATE_CHIRALITY,
      bool CHANGE_CHIRALITY,
      bool UNBIASED,
      bool SAMPLE_DIHEDRALS,
      bool CHECK_AMIDE_PLANARITY
    ) :
      m_FragmentData( FRAGMENT),
      m_UpdateRings( false),
      m_ChangeChirality( CHANGE_CHIRALITY),
      m_ChangeIsometry( CHANGE_CHIRALITY),
      m_Unbiased( UNBIASED),
      m_SampleDihedralAngles( SAMPLE_DIHEDRALS),
      m_MutateChirality( MUTATE_CHIRALITY),
      m_CheckAmideBondPlanarity( CHECK_AMIDE_PLANARITY)
    {
      if( FRAGMENT->ContainsRings() && !FRAGMENT->HasIncompleteRings())
      {
        if
        (
          !FRAGMENT->GetFragment().GetRotamerCoordinates().IsEmpty()
          && FRAGMENT->GetFragment().GetNumberDihedralChainBonds() == size_t( 0)
          && FRAGMENT->GetFragment().GetNumberChainBonds() == size_t( 0)
        )
        {
          // To prevent excessive ring swapping (which is slow and potentially reverts many changes that were made unrelated
          // to the ring due to a relatively naive ring-insertion algorithm), we only swap ring conformations if the
          // isomorophism only corresponds to the ring part
          m_UpdateRings = true;
        }
      }
      m_DihedralIndicesMutated = m_FragmentData->GetCenterBondIsomorphism();
      if( !m_UpdateRings && FRAGMENT->ContainsRings())
      {
        // need to exclude ring dihedral angles
        m_DihedralIndicesMutated.Reorder( m_FragmentData->GetNonRingBonds());
      }

      // for each of the non-ring bonds get atoms connected to it. Iterate over list of non-ring bonds
      for
      (
        auto itr( m_FragmentData->GetNonRingBonds().Begin()), itr_end( m_FragmentData->GetNonRingBonds().End());
        itr != itr_end;
        ++itr
      )
      {
        // for each non rotatable bond get connected atoms
        m_ConnectedAtoms[ *itr] = GetConnectedAtoms( m_FragmentData->GetCenterBonds()( *itr), GRAPH);
      }

      // get counts for all the rotamers to construct a probability histogram
      const storage::Vector< double> &rotamer_counts( m_FragmentData->GetRotamerCounts());
      storage::Map< double, double> histogram_vector;

      // get counts for all the rotamers to construct a probability histogram
      size_t rotamer_index( 0);
      for
      (
        storage::Vector< double>::const_iterator itr_count( rotamer_counts.Begin()), itr_count_end( rotamer_counts.End());
        itr_count != itr_count_end;
        ++itr_count, ++rotamer_index
      )
      {
        if( !m_Unbiased)
        {
          histogram_vector[ double( rotamer_index)] = *itr_count;
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
      if( m_UpdateRings)
      {
        m_FragmentMap = GetFragmentLinkWithComplement( GRAPH);
      }
      if( !m_ChangeChirality)
      {
        m_ChiralityString = MOLECULE.GetChiralityString();
      }
      if( !m_ChangeIsometry)
      {
        m_IsometryString = MOLECULE.GetDoubleBondIsometryString();
      }
      m_LastChosen.First() = m_LastChosen.Second() = util::GetUndefined< size_t>();
    }

    //! @brief Clone function
    //! @return pointer to new MutateFragment
    MutateFragment *MutateFragment::Clone() const
    {
      return new MutateFragment( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateFragment::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &MutateFragment::GetScheme() const
    {
      static const std::string s_ring( "MutateRing"), s_generic( "MutateFragment"), s_bond( "BiasedMutateDihedralBond");
      if( !m_FragmentData.IsDefined())
      {
        return s_generic;
      }
      return m_FragmentMap.IsEmpty()
             ? ( m_FragmentData->GetFragment().GetNumberDihedralChainBonds() == size_t( 1) ? s_bond : s_generic)
             : s_ring;
    }

    //! @brief get the dihedral bond indices referenced by this mutate
    const storage::Vector< size_t> &MutateFragment::GetMoleculeDihedralBondIndices() const
    {
      return m_DihedralIndicesMutated;
    }

    //! @brief returns the fragment that the mutate uses
    //! @return the fragment that the mutate uses
    const util::ShPtr< RotamerDihedralBondData> &MutateFragment::GetFragmentData() const
    {
      return m_FragmentData;
    }

    //! @brief returns atom indices that are connected to different bonds
    //! @return atom indices that are connected to different rotatable bonds
    const storage::Map< size_t, storage::Vector< size_t> > &MutateFragment::GetConnectedAtoms() const
    {
      return m_ConnectedAtoms;
    }

    //! @brief returns probability distribution for rotamers of this particular fragments
    //! @return probability distribution for rotamers of this particular fragments
    const random::Histogram1DDistribution &MutateFragment::GetProbability() const
    {
      return *m_ProbabilityDistribution;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator taking a conformation and returns a mutated conformation
    //! @param MOLECULE conformation of interest
    //! @return MutateResult that results from mutating to the argument
    math::MutateResult< FragmentComplete> MutateFragment::Mutate
    (
      const FragmentComplete &MOLECULE
    ) const
    {
      // the following debug block takes substantial time but is useful periodically to validate that we aren't uninentionally
      // changing chirality or isometry anywhere
#ifdef DEBUG_CHIRALITY_MUTATIONS
      // add isometry info, which may be different for a subgraph
      if( !m_ChangeChirality && !m_ChangeIsometry)
      {
        AtomVector< AtomComplete> atomv( MOLECULE.GetAtomInfo(), m_MutateChirality->GetBondInfo());
        BondIsometryHandler::AddIsometryInformation( atomv, true);
        StereocentersHandler::AddChiralityFromConformation( atomv);
        FragmentComplete fc( atomv, "");
        if( fc.GetChiralityString() != m_ChiralityString)
        {
          BCL_MessageCrt( "Alert: Chirality unintentionally changed: " + fc.GetChiralityString() + " from " + m_ChiralityString);
          io::OFStream output;
          io::File::MustOpenOFStream( output, "chirality_changed.sdf", std::ios::app);
          MOLECULE.WriteMDL( output);
          io::File::CloseClearFStream( output);
        }
        if( fc.GetDoubleBondIsometryString() != m_IsometryString)
        {
          BCL_MessageCrt( "Alert: Isometry unintentionally changed: " + fc.GetDoubleBondIsometryString() + " from " + m_IsometryString);
          io::OFStream output;
          io::File::MustOpenOFStream( output, "chirality_changed.sdf", std::ios::app);
          MOLECULE.WriteMDL( output);
          io::File::CloseClearFStream( output);
        }
      }
#endif

      #ifdef BCL_PROFILE_MutateFragment
      static util::Stopwatch s_timer_ring( "mutate_ring", util::Message::e_Standard, true, false);
      static util::Stopwatch s_timer_chain( "mutate_chain", util::Message::e_Standard, true, false);
      #endif

      double chosen_val( random::GetGlobalRandom().Double());

      // get a random iterator to get rotatable bonds related to a random isomorphism
      auto rotamer_bond_itr( m_FragmentData->GetRotamerBonds().Begin());
      auto itr_weights( m_FragmentData->GetIsomorphismWeights().Begin());
      while( rotamer_bond_itr != m_FragmentData->GetRotamerBonds().End() && chosen_val > *itr_weights)
      {
        chosen_val -= *itr_weights;
        ++itr_weights;
        ++rotamer_bond_itr;
      }
      if( rotamer_bond_itr == m_FragmentData->GetRotamerBonds().End())
      {
        --rotamer_bond_itr;
      }
      const size_t rotamer_chosen( std::distance( m_FragmentData->GetRotamerBonds().Begin(), rotamer_bond_itr));

      // determine which rotamer to use depending on the probability distribution
      size_t rotamer_to_use( GetChooseBestOnlyFlag() ? m_ProbabilityDistribution->DetermineMostLikelyCase() : m_ProbabilityDistribution->DetermineRandomCase());

      // get rotamer angle information for the particular rotamer of the fragment
      const linal::Vector< double> &rotamer_angles( m_FragmentData->GetRotamers()( rotamer_to_use));

      storage::Vector< sdf::AtomInfo> atom_info;

      // check whether fragment needs to be updated with new coordinates depending on rotamer chosen. If the argument
      // provideed already contains the rings in the same conformation then no need to swap coordintates. If rotamer
      // contains rings in differnet conformation as compared to argument then need to swap coordinates.
#ifdef DEBUG_CHIRALITY_MUTATIONS
      BCL_Debug( MOLECULE.HasBadGeometry());
#endif
      // if ring has different ring conformations then only there may be a possibilty that whole ring will have to be
      // swapped
      if( m_UpdateRings)
      {
#ifdef BCL_PROFILE_MutateFragment
        s_timer_ring.Start();
#endif
        if( m_DisallowedPairs.Contains( storage::Pair< size_t, size_t>( rotamer_chosen, rotamer_to_use)))
        {
          // already tried this combination of rotamer and isomer...and it leads to chirality being wrong in a non-correctable
          // way, so skip this move
          #ifdef BCL_PROFILE_MutateFragment
          s_timer_ring.Stop();
          #endif
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }
        storage::Pair< size_t, size_t> this_chosen( rotamer_chosen, rotamer_to_use);

        auto itr_iso( m_FragmentData->GetIsomorphisms().Begin());
        std::advance( itr_iso, rotamer_chosen);
        const storage::Vector< size_t> &isomorphism( *itr_iso);

        // get atom info of the fragment that this class mutates and the molecule of interest
        storage::Vector< sdf::AtomInfo> fragment_atoms( m_FragmentData->GetFragment().GetFragment().GetAtomInfo());

        // get coordinates for the rotamer which will be used
        const storage::Vector< linal::Vector3D> &rotamer_coordinates( m_FragmentData->GetFragment().GetRotamerCoordinates()( rotamer_to_use));

        // create a copy molecule of atom info and set undefined coordinates
        storage::Vector< sdf::AtomInfo> updated_ring_atominfo( MOLECULE.GetAtomInfo());
        for
        (
          storage::Vector< sdf::AtomInfo>::iterator
            itr_atominfo( updated_ring_atominfo.Begin()), itr_atominfo_end( updated_ring_atominfo.End());
          itr_atominfo != itr_atominfo_end;
          ++itr_atominfo
        )
        {
          itr_atominfo->SetCoordinates
          (
            linal::Vector3D
            (
              util::GetUndefined< double>(), util::GetUndefined< double>(), util::GetUndefined< double>()
            )
          );
        }

        // set the coordinates fragment to current rotamer and also set coordinates in molecule to current rotamer
        size_t count( 0);
        for
        (
          storage::Vector< size_t>::const_iterator
            itr_iso( isomorphism.Begin()), itr_iso_end( isomorphism.End());
          itr_iso != itr_iso_end;
          ++itr_iso, ++count
        )
        {
          updated_ring_atominfo( *itr_iso).SetCoordinates( rotamer_coordinates( count));
          fragment_atoms( count).SetCoordinates( rotamer_coordinates( count));
        }

        AtomVector< AtomComplete> new_fragment( fragment_atoms, m_FragmentData->GetFragment().GetFragment().GetBondInfo());

        // set the fragment as base fragment and grow the rest of the molecule on it
        for
        (
          storage::Map< size_t, storage::Vector< RingFragmentMap> >::const_iterator
            itr_mapping( m_FragmentMap( rotamer_chosen).Begin()), itr_mapping_end( m_FragmentMap( rotamer_chosen).End());
          itr_mapping != itr_mapping_end;
          ++itr_mapping
        )
        {
          MappingFragment( updated_ring_atominfo, isomorphism, MOLECULE, new_fragment, *itr_mapping);
        }

        // updating the dihedrals must be done after building the molecule up since there are some undefined coordinates
        // in updated_ring_atominfo until the end of the previous loop (since that loop re-assembles the molecule around the ring)
        for
        (
          storage::Map< size_t, storage::Vector< RingFragmentMap> >::const_iterator
            itr_mapping( m_FragmentMap( rotamer_chosen).Begin()), itr_mapping_end( m_FragmentMap( rotamer_chosen).End());
          itr_mapping != itr_mapping_end;
          ++itr_mapping
        )
        {
          FixDihedral( updated_ring_atominfo, MOLECULE, *itr_mapping);
        }

        // since ring has been updated, make sure the non-ring dihedral angle for the updated molecule are same as the
        // molecule passed in. This is for consistency. These will be changed later to correspond to specific angle for the rotamer.
        BCL_Assert( m_ConnectedAtoms.IsEmpty(), "Expected to not have dihedral bonds to rotate");

        AtomVector< AtomComplete> atoms( updated_ring_atominfo, m_MutateChirality->GetBondInfo());

        // add isometry info, which may be different for a subgraph
        BondIsometryHandler::AddIsometryInformation( atoms, true);

        // add stereocenter information
        StereocentersHandler::AddChiralityFromConformation( atoms);

        bool altered_chirality( false);
        if( !m_ChangeChirality || !m_ChangeIsometry)
        {
          FragmentComplete temporary_mol( atoms, MOLECULE.GetName());
          std::string gen_chirality( temporary_mol.GetChiralityString());
          std::string gen_dbi( temporary_mol.GetDoubleBondIsometryString());
          if( !m_ChangeChirality && m_ChiralityString != gen_chirality)
          {
            BCL_MessageVrb( "Trying to fix chirality or double bond isometry that was wrong: " + m_ChiralityString + " " + gen_chirality + " " + m_IsometryString + " " + gen_dbi);
            for( size_t atom_index( 0), n_atoms( MOLECULE.GetSize()); atom_index != n_atoms; ++atom_index)
            {
              const AtomConformationalInterface &mol_atom( MOLECULE.GetAtomVector()( atom_index));
              if( mol_atom.GetChirality() != e_NonChiral)
              {
                if( temporary_mol.GetAtomVector()( atom_index).GetChirality() == mol_atom.GetChirality())
                {
                  continue;
                }
                if( atoms( atom_index).CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, 1) >= 4)
                {
                  BCL_MessageVrb( "Cannot fix as it's on atom: " + util::Format()( atom_index) + " which has 4 bonds in rings)");
                  m_DisallowedPairs.InsertElement( storage::Pair< size_t, size_t>( rotamer_chosen, rotamer_to_use));
                  return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
                }
                BCL_MessageVrb( "Trying to fix chirality at atom index (4-bonds): " + util::Format()( atom_index));
                temporary_mol = m_MutateChirality->MutateSpecificChirality( temporary_mol, atom_index);
                altered_chirality = true;
                if
                (
                  mol_atom.GetChirality() == e_TransRingChirality || mol_atom.GetChirality() == e_UnknownRingChirality
                  || mol_atom.GetChirality() == e_CisRingChirality
                )
                {
                  atoms = temporary_mol.GetAtomVector();
                  StereocentersHandler::AddChiralityFromConformation( atoms);
                  temporary_mol = FragmentComplete( atoms, temporary_mol.GetName());
                }
              }
            }
            if( m_ChiralityString != temporary_mol.GetChiralityString())
            {
              BCL_MessageStd
              (
                "Altered chirality in unexpected way. Should have been: " + m_ChiralityString
                + " but was " + temporary_mol.GetChiralityString()
              );
              return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
            }
          }
          if( !m_ChangeIsometry && m_IsometryString != gen_dbi)
          {
            BCL_MessageVrb( "Trying to fix isometry to " + m_IsometryString + " from " + gen_dbi);
            for( size_t atom_index( 0), n_atoms( MOLECULE.GetSize()); atom_index != n_atoms; ++atom_index)
            {
              if( atoms( atom_index).GetAtomType()->GetNumberBonds() == size_t( 3) && !m_ChangeIsometry)
              {
                const AtomConformationalInterface &atom( MOLECULE.GetAtomVector()( atom_index));
                const AtomConformationalInterface &tmp_atom( atoms( atom_index));
                // skip non-isometric bonds
                if( atom.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsIsometric, 1) != 1)
                {
                  continue;
                }
                // skip this atom if it is part of a double bond on the wrong side of the ring
                if( atom.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, 1) >= 1)
                {
                  BCL_MessageVrb( "Cannot fix as it's on atom: " + util::Format()( atom_index) + " which has bonds in ring)");
                  m_DisallowedPairs.InsertElement( storage::Pair< size_t, size_t>( rotamer_chosen, rotamer_to_use));
                  return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
                }
                // skip cases with the same isometry
                if
                (
                  tmp_atom.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_Isometry, 1)
                  == atom.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_Isometry, 1)
                )
                {
                  continue;
                }
                BCL_MessageStd( "Trying to fix isometry at atom index: " + util::Format()( atom_index));
                temporary_mol = m_MutateChirality->MutateSpecificChirality( temporary_mol, atom_index);
                altered_chirality = true;
              }
            }
          }

          // if chirality was fixed, it may be necessary to re-adjust the dihedral angles again
          if( altered_chirality)
          {
            updated_ring_atominfo = temporary_mol.GetAtomInfo();
            // updating the dihedrals must be done after building the molecule up since there are some undefined coordinates
            // in updated_ring_atominfo until the end of the previous loop (since that loop re-assembles the molecule around the ring)
            for
            (
              storage::Map< size_t, storage::Vector< RingFragmentMap> >::const_iterator
                itr_mapping( m_FragmentMap( rotamer_chosen).Begin()), itr_mapping_end( m_FragmentMap( rotamer_chosen).End());
              itr_mapping != itr_mapping_end;
              ++itr_mapping
            )
            {
              FixDihedral( updated_ring_atominfo, MOLECULE, *itr_mapping);
            }
            atoms = AtomVector< AtomComplete>( updated_ring_atominfo, m_MutateChirality->GetBondInfo());
          }
          else
          {
            atoms = temporary_mol.GetAtomVector();
          }

          // add stereocenter information
          #ifdef DEBUG_CHIRALITY_MUTATIONS
          FragmentComplete bmol( atoms, MOLECULE.GetName());
          if
          (
            bmol.GetChiralityString() != m_ChiralityString
          )
          {
            BCL_MessageStd( "Evidently Failed to fix chirality! " + bmol.GetChiralityString());
          }
          // add isometry info, which may be different for a subgraph
          BondIsometryHandler::AddIsometryInformation( atoms, true);

          // add stereocenter information
          StereocentersHandler::AddChiralityFromConformation( atoms);

          FragmentComplete tmol( atoms, MOLECULE.GetName());
          if
          (
            tmol.GetChiralityString() != m_ChiralityString
          )
          {
            BCL_MessageStd( "Really Failed to fix chirality! " + tmol.GetChiralityString());
          }
          if( tmol.GetDoubleBondIsometryString() != m_IsometryString)
          {
            BCL_MessageStd( " Failed to fix isometry! " + tmol.GetDoubleBondIsometryString());
          }
          BCL_Debug( tmol.HasBadGeometry());
          #endif
        }

        util::ShPtr< FragmentComplete> new_molecule( new FragmentComplete( atoms, MOLECULE.GetName()));

        #ifdef BCL_PROFILE_MutateFragment
          s_timer_ring.Stop();
        #endif
        // return the mutated argument
        // check for severely non-planar amide bonds
        if( m_ConnectedAtoms.GetSize() == 0)
        {
          const double total_non_planarity_final
          (
            m_CheckAmideBondPlanarity
            ? new_molecule->GetTotalAmideBondNonPlanarity( double( 25))
            : 0.0
          );

          if( total_non_planarity_final)
          {
            // ensure that we haven't made amide bonds even more non-planar (this can happen because fragments can split at
            // amide-bond boundaries or be used at amide-bond boundaries)
            const double total_non_planarity_original( MOLECULE.GetTotalAmideBondNonPlanarity( double( 25)));
            if( total_non_planarity_original + 0.01 < total_non_planarity_final)
            {
              BCL_MessageStd( "Skipping due to nonplanarity of amide bonds Multiple Bond");
              // return the mutated argument
              return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
            }
          }
          return math::MutateResult< FragmentComplete>( new_molecule, *this);
        }
        atom_info = new_molecule->GetAtomInfo();
      }
      else
      {
        atom_info = MOLECULE.GetAtomInfo();
      }

      #ifdef BCL_PROFILE_MutateFragment
      s_timer_chain.Start();
      #endif
      // get dihedral angle associated with rotatable bonds associated with rotamer_bond_itr
      storage::Vector< double> existing_angles( GetDihedralAngles( atom_info, *rotamer_bond_itr));

      for
      (
        storage::Map< size_t, storage::Vector< size_t> >::const_iterator
          itr_bonds( m_ConnectedAtoms.Begin()), itr_bonds_end( m_ConnectedAtoms.End());
        itr_bonds != itr_bonds_end;
        ++itr_bonds
      )
      {
        if
        (
          !m_ChangeIsometry
          && MOLECULE.GetAtomVector()( itr_bonds->second( 0)).GetBondTypeTo
          (
            MOLECULE.GetAtomVector()( itr_bonds->second( 1))
          )->GetBondData( ConfigurationalBondTypeData::e_IsIsometric) == size_t( 1)
        )
        {
          continue;
        }
        const double effective_std
        (
          BondDihedralAngles::GetEstimatedStdForDihedralBondAngleBin
          (
            MOLECULE.GetAtomVector()( itr_bonds->second( 0)),
            MOLECULE.GetAtomVector()( itr_bonds->second( 1))
          )
        );
        double new_angle( rotamer_angles( itr_bonds->first));
        RotateBond
        (
          atom_info,
          itr_bonds->second,
          random::GetGlobalRandom().RandomGaussian( new_angle, effective_std),
          existing_angles( itr_bonds->first)
        );
      }

      // create new fragment from the updated coordinaes
      util::ShPtr< FragmentComplete> new_molecule
      (
        new FragmentComplete( AtomVector< AtomComplete>( atom_info, m_MutateChirality->GetBondInfo()), MOLECULE.GetName())
      );
      #ifdef DEBUG_CHIRALITY_MUTATIONS
      BCL_Debug( new_molecule->HasBadGeometry());
      #endif
      #ifdef BCL_PROFILE_MutateFragment
      s_timer_chain.Stop();
      #endif
      const double total_non_planarity_final
      (
        m_CheckAmideBondPlanarity
        ? new_molecule->GetTotalAmideBondNonPlanarity( double( 20))
        : 0.0
      );

      // This commented out block was included when the algorithm was benchmarked, but at that time a bug in AreAmideBondsPlaner
      // caused it to only return true for terminal nitrogen groups (e.g. non-hydrogenated, terminal nitrogen) two bonds
      // away from an oxygen double bond; and was not checking planarity at all. Therefore, this block needs to be tested
      // with the bug-fixed version of the code before using this feature in production setting
      if( total_non_planarity_final)
      {
        // ensure that we haven't made amide bonds even more non-planar (this can happen because fragments can split at
        // amide-bond boundaries or be used at amide-bond boundaries)
        const double total_non_planarity_original( MOLECULE.GetTotalAmideBondNonPlanarity( double( 20)));
        if( total_non_planarity_original + 0.01 < total_non_planarity_final)
        {
          BCL_MessageDbg
          (
            "Skipping due to nonplanarity of amide bonds Single Bond "
            + util::Format()( total_non_planarity_original) + " degrees originally, now: " + util::Format()( total_non_planarity_final)
          );
          // return the mutated argument
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }
      }
      return math::MutateResult< FragmentComplete>( new_molecule, *this);
    }

    //! @brief return a vector of atoms that are connected to a bond of interest
    //! @param ATOM_INDICES atom indices of the center bond around which connected bonds have to be found
    //! @param GRAPH constitution graph of the molecule of interest
    //! @param return a vector of atoms that are connected to a bond of interest
    storage::Vector< size_t> MutateFragment::GetConnectedAtoms
    (
      const graph::UndirectedEdge< ConfigurationalBondType> &ATOM_INDICES,
      const graph::ConstGraph< size_t, size_t> &GRAPH
    )
    {
      // get the atom indices
      size_t second_atom_index( ATOM_INDICES.GetIndexLow());
      size_t third_atom_index( ATOM_INDICES.GetIndexHigh());

      // find the part of graph that contains second_atom_index but does not contain third_atom_index
      storage::Vector< size_t> connected_vertices_b
      (
        *graph::Connectivity::GetVerticesReachableFromDirectedEdge( GRAPH, second_atom_index, third_atom_index)
      );

      // if number of atoms connected to second_atom_index is less than half molecule size then this is the smallest
      // set of atoms that need to be rotated. Otherwise find atoms connected to third atom index
      if( connected_vertices_b.GetSize() < GRAPH.GetSize() / 2)
      {
        size_t removable_element( connected_vertices_b.Find( second_atom_index));
        if( removable_element < connected_vertices_b.GetSize())
        {
          connected_vertices_b.RemoveElements( removable_element, 1);
        }
        storage::Set< size_t> set_vert( connected_vertices_b.Begin(), connected_vertices_b.End());
        storage::Vector< size_t> connected_vertices( set_vert.Begin(), set_vert.End());
        connected_vertices.InsertElements( 0, third_atom_index, 1);
        connected_vertices.InsertElements( 1, second_atom_index, 1);
        return connected_vertices;
      }

      storage::Vector< size_t> connected_vertices_c
      (
        *graph::Connectivity::GetVerticesReachableFromDirectedEdge( GRAPH, third_atom_index, second_atom_index)
      );
      size_t removable_element( connected_vertices_b.Find( third_atom_index));
      if( removable_element < connected_vertices_b.GetSize())
      {
        connected_vertices_b.RemoveElements( removable_element, 1);
      }

      storage::Set< size_t> set_vert( connected_vertices_c.Begin(), connected_vertices_c.End());
      storage::Vector< size_t> connected_vertices( set_vert.Begin(), set_vert.End());
      connected_vertices.InsertElements( 0, second_atom_index, 1);
      connected_vertices.InsertElements( 1, third_atom_index, 1);
      return connected_vertices;
    }

    //! @brief rotate the molecule of interest at a particular bond at a particular angle
    //! @param ATOM_INFO the atoms that needs to be updated with coordinates after rotatating about a particular bond
    //! @param CONNECTED_ATOMS the connected atoms associated with bond of interest
    //! @param NEW_ANGLE the new angle for the bond of interest
    //! @param EXISTING_ANGLE the angle for the bond of interest in the argument that is passed in
    void MutateFragment::RotateBond
    (
      storage::Vector< sdf::AtomInfo> &ATOM_INFO,
      const storage::Vector< size_t> &CONNECTED_ATOMS,
      double NEW_ANGLE,
      double EXISTING_ANGLE
    )
    {
      #ifdef BCL_PROFILE_MutateFragment
      static util::Stopwatch s_timer_create_transform( "create_transformation_matrix", util::Message::e_Standard, true, false);
      static util::Stopwatch s_timer_rotate( "do_rotate_bond", util::Message::e_Standard, true, false);
      s_timer_create_transform.Start();
      #endif
      double degree_rotation( EXISTING_ANGLE - NEW_ANGLE);

      // get the anchor points for bond rotation, i.e. indices of atoms that make the bond of interest
      storage::Vector< size_t>::const_iterator itr_connected_atoms( CONNECTED_ATOMS.Begin());
      const linal::Vector3D &start_point( ATOM_INFO( *itr_connected_atoms).GetCoordinates());
      ++itr_connected_atoms;
      const linal::Vector3D &end_point( ATOM_INFO( *itr_connected_atoms).GetCoordinates());
      ++itr_connected_atoms;

      math::TransformationMatrix3D transformation
      (
        coord::LineSegment3D( start_point, end_point),
        start_point,
        math::Angle::Radian( degree_rotation)
      );
      #ifdef BCL_PROFILE_MutateFragment
      s_timer_create_transform.Stop();
      s_timer_rotate.Start();
      #endif
      // set the coordinates of the atoms from the transformed fragment cloud
      for
      (
        storage::Vector< size_t>::const_iterator itr_connected_atoms_end( CONNECTED_ATOMS.End());
        itr_connected_atoms != itr_connected_atoms_end;
        ++itr_connected_atoms
      )
      {
        linal::Vector3D point( ATOM_INFO( *itr_connected_atoms).GetCoordinates());
        ATOM_INFO( *itr_connected_atoms).SetCoordinates( point.Transform( transformation));
      }
      #ifdef BCL_PROFILE_MutateFragment
      s_timer_rotate.Stop();
      #endif
    }

    //! @brief get a vector containing angle for each rotatable bond in the molecule
    //! @param ATOM_INFO the atom info of molecule of interest
    //! @param DIHEDRAL_BONDS atom indices making up the different rotatable bonds in the molecule
    //! @return vector containing angle for each rotatable bond in the molecule
    storage::Vector< double> MutateFragment::GetDihedralAngles
    (
      const storage::Vector< sdf::AtomInfo> &ATOM_INFO,
      const storage::Vector< storage::VectorND< 4, size_t> > &DIHEDRAL_BONDS
    )
    {
      // create a vector to store the angle of interest
      storage::Vector< double> bond_angles;
      bond_angles.AllocateMemory( DIHEDRAL_BONDS.GetSize());
      for
      (
        storage::Vector< storage::VectorND< 4, size_t> >::const_iterator
          itr( DIHEDRAL_BONDS.Begin()), itr_end( DIHEDRAL_BONDS.End());
        itr != itr_end;
        ++itr
      )
      {
        // push back into the vector
        bond_angles.PushBack( GetDihedralAngle( ATOM_INFO, *itr));
      }
      return bond_angles;
    }

    //! @brief get a vector containing angle for each rotatable bond in the molecule
    //! @param ATOM_INFO the atom info of molecule of interest
    //! @param DIHEDRAL_BONDS atom indices making up the different rotatable bonds in the molecule
    //! @return vector containing angle for each rotatable bond in the molecule
    double MutateFragment::GetDihedralAngle
    (
      const storage::Vector< sdf::AtomInfo> &ATOM_INFO,
      const storage::VectorND< 4, size_t> &DIHEDRAL_BOND
    )
    {
      // get the angle measure for the angle being iterated through
      double existing_angle
      (
        math::Angle::Degree
        (
          linal::Dihedral
          (
            ATOM_INFO( DIHEDRAL_BOND.First()).GetCoordinates(),
            ATOM_INFO( DIHEDRAL_BOND.Second()).GetCoordinates(),
            ATOM_INFO( DIHEDRAL_BOND.Third()).GetCoordinates(),
            ATOM_INFO( DIHEDRAL_BOND.Fourth()).GetCoordinates()
          )
        )
      );

      // take care of the wrap around
      if( existing_angle < 15.0)
      {
        existing_angle = 360 + existing_angle;
      }
      return existing_angle;
    }

    //! @brief get complement of fragment in the molecule and information about bonds that can connect fragment to other parts of the molecule
    //! @param GRAPH the constitution graph of the molecule on which this mutates operates
    //! @return
    storage::Vector< storage::Map< size_t, storage::Vector< RingFragmentMap> > > MutateFragment::GetFragmentLinkWithComplement
    (
      const graph::ConstGraph< size_t, size_t> &GRAPH
    ) const
    {
      storage::Vector< storage::Map< size_t, storage::Vector< RingFragmentMap> > > link_complement;
      link_complement.AllocateMemory( m_FragmentData->GetIsomorphisms().GetSize());

      for
      (
        auto itr_iso( m_FragmentData->GetIsomorphisms().Begin()), itr_iso_end( m_FragmentData->GetIsomorphisms().End());
        itr_iso != itr_iso_end;
        ++itr_iso
      )
      {
        graph::Subgraph< size_t, size_t> fragment_graph
        (
          util::OwnPtr< const graph::ConstGraph< size_t, size_t> >( &GRAPH, false),
          *itr_iso
        );
        storage::List< storage::Pair< size_t, size_t> > adjacent_edges( fragment_graph.GetAdjacentEdgeIndices());
        storage::Map< size_t, storage::Vector< RingFragmentMap> > fragment_connections;
        for
        (
          storage::List< storage::Pair< size_t, size_t> >::const_iterator
            itr_adj( adjacent_edges.Begin()), itr_adj_end( adjacent_edges.End());
          itr_adj != itr_adj_end;
          ++itr_adj
        )
        {
          //BCL_Assert( GRAPH.GetEdgeData( itr_adj->First(), itr_adj->Second()) == 0, "Incomplete ring mapped...bug!");
          //if( GRAPH.GetEdgeData( itr_adj->First(), itr_adj->Second()) == 0)
          {
            fragment_connections[ itr_adj->First()].PushBack
            (
              RingFragmentMap
              (
                itr_adj->Second(),
                false,
                *graph::Connectivity::GetVerticesReachableFromDirectedEdge( GRAPH, itr_adj->Second(), itr_adj->First()),
                GetConnectedAtoms
                (
                  graph::UndirectedEdge< ConfigurationalBondType>
                  (
                    itr_adj->Second(),
                    itr_adj->First(),
                    ConfigurationalBondType()
                  ),
                  GRAPH
                )
              )
            );
          }
        }
        link_complement.PushBack( fragment_connections);
      }

      return link_complement;
    }

    //! @brief maps and grows a part of molecule on fragment that this mutate works on
    //! @param ATOM_INFO atom info of the new conformation
    //! @param ISOMORPHISM isomorphism between fragment and molecule of interest
    //! @param MOLECULE_ATOM_INFO atom info of conformation that is given
    //! @param FRAGMENT_ATOM_INFO atom info of fragment that this mutate uses
    //! @param MAPPING_INFORMATION connectivity between fragment and rest of the molecule
    void MutateFragment::MappingFragment
    (
      storage::Vector< sdf::AtomInfo> &ATOM_INFO,
      const storage::Vector< size_t> &ISOMORPHISM,
      const FragmentComplete &MOLECULE,
      const AtomVector< AtomComplete> &FRAGMENT,
      const std::pair< size_t, storage::Vector< RingFragmentMap> > &MAPPING_INFORMATION
    ) const
    {
      const size_t fragment_atom_id( MAPPING_INFORMATION.first);
      const size_t fragment_atom_isomorphism_id( ISOMORPHISM.Find( fragment_atom_id));

      const AtomComplete &fragment_atom( FRAGMENT( fragment_atom_isomorphism_id));
      storage::Vector< linal::Vector3D> fragment_valence( ValenceHandler::DetermineCoordinates( fragment_atom));
      if( fragment_valence.GetSize() > size_t( 1))
      {
        fragment_valence.Shuffle();
      }
      const linal::Vector3D &molecule_connector_coordinates( MOLECULE.GetAtomVector()( fragment_atom_id).GetPosition());
      const linal::Vector3D &fragment_connector_coordinates( fragment_atom.GetPosition());

      storage::Vector< linal::Vector3D>::const_iterator itr_val( fragment_valence.Begin());
      for
      (
        storage::Vector< RingFragmentMap>::const_iterator
          itr( MAPPING_INFORMATION.second.Begin()), itr_end( MAPPING_INFORMATION.second.End());
        itr != itr_end;
        ++itr, ++itr_val
      )
      {
        const size_t molecule_atom_id( itr->GetMoleculeAtom());

        const linal::Vector3D molecule_complement_coordinates
        (
          MOLECULE.GetAtomVector()( molecule_atom_id).GetPosition()
        );

        linal::Vector3D fragment_complement_coordinates;

        BCL_Assert( itr_val != fragment_valence.End(), "Reached end of fragment valence vector prematurely");
        const double bond_length( linal::Distance( molecule_complement_coordinates, molecule_connector_coordinates));
        fragment_complement_coordinates = ( *itr_val - fragment_connector_coordinates) * bond_length + fragment_connector_coordinates;
        //BCL_Debug( ( *itr_val - fragment_connector_coordinates).Norm());
        const storage::Vector< size_t> &non_overlapping_atoms( itr->GetConnectedVertices());

        // now we find the transformation matrix so as to superimpose defined vertices of molecule and fragment
        math::TransformationMatrix3D transform
        (
          coord::LineSegment3D( fragment_complement_coordinates, fragment_connector_coordinates),
          coord::LineSegment3D( molecule_complement_coordinates, molecule_connector_coordinates)
        );

        for
        (
            storage::Vector< size_t>::const_iterator
            itr_non_overlapping( non_overlapping_atoms.Begin()), itr_non_overlapping_end( non_overlapping_atoms.End());
            itr_non_overlapping != itr_non_overlapping_end;
            ++itr_non_overlapping
        )
        {
          linal::Vector3D point
          (
            MOLECULE.GetAtomVector()( *itr_non_overlapping).GetPosition()
          );
          ATOM_INFO( *itr_non_overlapping).SetCoordinates( point.Transform( transform));
        }
      }
    }
    //! @brief maps and grows a part of molecule on fragment that this mutate works on
    //! @param ATOM_INFO atom info of the new conformation
    //! @param ISOMORPHISM isomorphism between fragment and molecule of interest
    //! @param MOLECULE_ATOM_INFO atom info of conformation that is given
    //! @param FRAGMENT_ATOM_INFO atom info of fragment that this mutate uses
    //! @param MAPPING_INFORMATION connectivity between fragment and rest of the molecule
    void MutateFragment::FixDihedral
    (
      storage::Vector< sdf::AtomInfo> &ATOM_INFO,
      const FragmentComplete &MOLECULE,
      const std::pair< size_t, storage::Vector< RingFragmentMap> > &MAPPING_INFORMATION
    ) const
    {
      const size_t fragment_atom_id( MAPPING_INFORMATION.first);

      size_t dihedral_atom_index_a( util::GetUndefined< size_t>());
      for
      (
        auto itr_bnd( MOLECULE.GetAtomVector()( fragment_atom_id).GetBonds().Begin()),
             itr_bnd_end( MOLECULE.GetAtomVector()( fragment_atom_id).GetBonds().End());
        itr_bnd != itr_bnd_end;
        ++itr_bnd
      )
      {
        if( itr_bnd->GetBondType()->IsBondInRing())
        {
          dihedral_atom_index_a = MOLECULE.GetAtomIndex( itr_bnd->GetTargetAtom());
          break;
        }
      }

      for
      (
        storage::Vector< RingFragmentMap>::const_iterator
          itr( MAPPING_INFORMATION.second.Begin()), itr_end( MAPPING_INFORMATION.second.End());
        itr != itr_end;
        ++itr
      )
      {
        const size_t molecule_atom_id( itr->GetMoleculeAtom());

        if
        (
          util::IsDefined( dihedral_atom_index_a)
          && MOLECULE.GetAtomVector()( molecule_atom_id).GetBonds().GetSize() > size_t( 1)
        )
        {
          double dihedral( util::GetUndefined< double>());
          size_t dihedral_atom_index_d( util::GetUndefined< size_t>());
          for
          (
            auto itr_bnd( MOLECULE.GetAtomVector()( molecule_atom_id).GetBonds().Begin()),
                 itr_bnd_end( MOLECULE.GetAtomVector()( molecule_atom_id).GetBonds().End());
            itr_bnd != itr_bnd_end;
            ++itr_bnd
          )
          {
            if( MOLECULE.GetAtomIndex( itr_bnd->GetTargetAtom()) != fragment_atom_id)
            {
              dihedral_atom_index_d = MOLECULE.GetAtomIndex( itr_bnd->GetTargetAtom());
              dihedral =
                linal::Dihedral
                (
                  MOLECULE.GetAtomVector()( dihedral_atom_index_a).GetPosition(),
                  MOLECULE.GetAtomVector()( fragment_atom_id).GetPosition(),
                  MOLECULE.GetAtomVector()( molecule_atom_id).GetPosition(),
                  MOLECULE.GetAtomVector()( dihedral_atom_index_d).GetPosition()
                );
              break;
            }
          }

          if( util::IsDefined( dihedral))
          {
            double new_dihedral
            (
              linal::Dihedral
              (
                ATOM_INFO( dihedral_atom_index_a).GetCoordinates(),
                ATOM_INFO( fragment_atom_id).GetCoordinates(),
                ATOM_INFO( molecule_atom_id).GetCoordinates(),
                ATOM_INFO( dihedral_atom_index_d).GetCoordinates()
              )
            );
            if( !math::EqualWithinAbsoluteTolerance( dihedral, new_dihedral, 0.01))
            {
              RotateBond( ATOM_INFO, itr->GetConnectedVerticesForRotation(), math::Angle::Degree( dihedral), math::Angle::Degree( new_dihedral));
            }
          }
        }
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &MutateFragment::Read( std::istream &ISTREAM)
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
    std::ostream &MutateFragment::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ConnectedAtoms, OSTREAM, INDENT);
      io::Serialize::Write( m_ProbabilityDistribution, OSTREAM, INDENT);

      // return "OSTREAM"
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
