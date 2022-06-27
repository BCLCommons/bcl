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
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_constitution_set.h"
#include "chemistry/bcl_chemistry_fragment_constitution_shared.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_fragment_split_ecfp_fragments.h"
#include "chemistry/bcl_chemistry_fragment_split_rigid.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "descriptor/bcl_descriptor_molecule_druglike.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_conformational.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "descriptor/bcl_descriptor_molecule_atom_environment_map.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serialize.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // initialize static
    sched::Mutex &MoleculeDruglike::GetMutex()
    {
      static sched::Mutex s_Mutex;
      return s_Mutex;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

      //! @brief default constructor
    MoleculeDruglike::MoleculeDruglike() :
      m_MinWeight( 10.0),
      m_MaxWeight( 600.0),
      m_MaxRingSize( 10),
      m_MinLogP( -1.0),
      m_MaxLogP( 5.0),
      m_MaxTPSA( 140.0),
      m_MaxHBDA( 12),
      m_MaxNRotBonds( 10),
      m_MaxTotalBondEnergy( 0.10),
      m_MaxNF( 5),
      m_MaxNCl( 2),
      m_MaxNBr( 1),
      m_MaxNI( 1),
      m_MaxNHalogens( 5),
      m_MaxRingHalogens( 3),
      m_MaxNonAromaticClBrI( 0),
      m_MaxComplexity( 4.0),
      m_EnforceHitlike( false),
      m_EnforceDruglikeRings( false),
      m_EnforceDruglikeECFP( false)
    {
    }

    //! @brief hitlike constructor
    MoleculeDruglike::MoleculeDruglike
    (
      const bool &ENFORCE_HITLIKE
    ) :
      m_MinWeight( 10.0),
      m_MaxWeight( 600.0),
      m_MaxRingSize( 10),
      m_MinLogP( -1.0),
      m_MaxLogP( 5.0),
      m_MaxTPSA( 140.0),
      m_MaxHBDA( 12),
      m_MaxNRotBonds( 10),
      m_MaxTotalBondEnergy( 0.10),
      m_MaxNF( 5),
      m_MaxNCl( 2),
      m_MaxNBr( 1),
      m_MaxNI( 1),
      m_MaxNHalogens( 5),
      m_MaxRingHalogens( 3),
      m_MaxNonAromaticClBrI( 0),
      m_MaxComplexity( 4.0),
      m_EnforceHitlike( ENFORCE_HITLIKE),
      m_EnforceDruglikeRings( false),
      m_EnforceDruglikeECFP( false)
    {
    }

      //! @brief constructor
    MoleculeDruglike::MoleculeDruglike
      (
        const float &MIN_WEIGHT,
        const float &MAX_WEIGHT,
        const size_t &MAX_RING_SIZE,
        const float &MIN_LOGP,
        const float &MAX_LOGP,
        const float &MAX_TPSA,
        const size_t &MAX_HBDA,
        const size_t &MAX_N_ROT_BONDS,
        const float &MAX_BOND_PROPENSITY,
        const size_t &MAX_N_F,
        const size_t &MAX_N_CL,
        const size_t &MAX_N_BR,
        const size_t &MAX_N_I,
        const size_t &MAX_TOTAL_HALOGENS,
        const size_t &MAX_RING_HALOGENS,
        const size_t &MAX_NON_AROMATIC_CLBRI,
        const float &MAX_COMPLEXITY,
        const bool &ENFORCE_HITLIKE,
        const bool &ENFORCE_DRUGLIKE_RINGS,
        const bool &ENFORCE_DRUGLIKE_ECFP
      ) :
      m_MinWeight( MIN_WEIGHT),
      m_MaxWeight( MAX_WEIGHT),
      m_MaxRingSize( MAX_RING_SIZE),
      m_MinLogP( MIN_LOGP),
      m_MaxLogP( MAX_LOGP),
      m_MaxTPSA( MAX_TPSA),
      m_MaxHBDA( MAX_HBDA),
      m_MaxNRotBonds( MAX_N_ROT_BONDS),
      m_MaxTotalBondEnergy( MAX_BOND_PROPENSITY),
      m_MaxNF( MAX_N_F),
      m_MaxNCl( MAX_N_CL),
      m_MaxNBr( MAX_N_BR),
      m_MaxNI( MAX_N_I),
      m_MaxNHalogens( MAX_TOTAL_HALOGENS),
      m_MaxRingHalogens( MAX_RING_HALOGENS),
      m_MaxNonAromaticClBrI( MAX_NON_AROMATIC_CLBRI),
      m_MaxComplexity( MAX_COMPLEXITY),
      m_EnforceHitlike( ENFORCE_HITLIKE),
      m_EnforceDruglikeRings( ENFORCE_DRUGLIKE_RINGS),
      m_EnforceDruglikeECFP( ENFORCE_DRUGLIKE_ECFP)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeDruglike
    MoleculeDruglike *MoleculeDruglike::Clone() const
    {
      return new MoleculeDruglike( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeDruglike::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeDruglike::GetAlias() const
    {
      static const std::string s_name( "IsConstitutionDruglike"), s_hitlike_name( "IsConstitutionDruglikeAndHitlike");
      return m_EnforceHitlike ? s_hitlike_name : s_name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeDruglike::Calculate( linal::VectorReference< float> &STORAGE)
    {
      // Get our molecule
      util::SiPtr< const chemistry::ConformationInterface> current_mol( this->GetCurrentObject());
      chemistry::FragmentComplete molecule( *current_mol);

      // if enforce hit-like is on then do that first before bothering with the rest
      if( m_EnforceHitlike)
      {
        // count the number of aromatic rings
        size_t n_aro_rings( GetCheminfoProperties().calc_NAromaticRings->SumOverObject( molecule)( 0));

        // for molecules with only 1 aromatic ring
        if( n_aro_rings == size_t( 1))
        {
          if
          (
              GetCheminfoProperties().calc_MolWeight->SumOverObject( molecule)( 0) > float( 217.0) &&
              GetCheminfoProperties().calc_MolWeight->SumOverObject( molecule)( 0) < float( 478.0) &&
              GetCheminfoProperties().calc_TopologicalPolarSurfaceArea->SumOverObject( molecule)( 0) > float( 57.0)
          )
          {
            // continue to druglike evaluation
          }
          else
          {
            STORAGE( 0) = 0;
            return;
          }
        }
        // for molecules with more than 1 aromatic ring
        else if( n_aro_rings > size_t( 1))
        {
          if
          (
              GetCheminfoProperties().calc_MolWeight->SumOverObject( molecule)( 0) < float( 520.0) &&
              GetCheminfoProperties().calc_TopologicalPolarSurfaceArea->SumOverObject( molecule)( 0) > float( 40.0) &&
              std::max
              (
                bool( GetCheminfoProperties().calc_NRotBond->SumOverObject( molecule)( 0) < size_t( 3)),
                bool( GetCheminfoProperties().calc_HbondDonor->SumOverObject( molecule)( 0) < size_t( 4))
              )
          )
          {
            // continue to druglike evaluation
          }
          else
          {
            STORAGE( 0) = 0;
            return;
          }
        }
        // for molecules with no aromatic rings
        else
        {
          STORAGE( 0) = 0;
          return;
        }
      }

      // check for druglike rings against ring database explicitly
      if( m_EnforceDruglikeRings && !ContainsOnlyDruglikeRings( molecule)) //!< slow and generally deprecated if we filter for druglike fragments
      {
        STORAGE( 0) = 0;
        return;
      }

      // check for druglike ecfp fragments
      GetMutex().Lock();
      static CheminfoProperty atom_environment;
      if( !command::CommandState::IsInStaticInitialization())
      {
        atom_environment = CheminfoProperty( "MoleculeAtomEnvironmentMap");
        GetMutex().Unlock();
        if( m_EnforceDruglikeECFP && !atom_environment->SumOverObject( molecule)( 0))
        {
          STORAGE( 0) = 0;
          return;
        }
      }
      else
      {
        std::string atom_environment_map_file
        (
          chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") +
          (
              "atom_environments/chembl.h.atom_environment_hashmap.two_bonds.txt.gz"
          )
        );
        static MoleculeAtomEnvironmentMap aem( atom_environment_map_file, false, size_t( 2));
        GetMutex().Unlock();
        if( m_EnforceDruglikeECFP && !aem.SumOverObject( molecule)( 0))
        {
          STORAGE( 0) = 0;
          return;
        }
      }

      // create bond propensity score descriptor object
      static CheminfoProperty bonde( "MoleculeTotalBondEnergy");

    ////////////////////
    //     DEBUG      //
    ////////////////////

      // uncomment for debug
//      BCL_Debug(molecule.GetNumberAtoms());
//      BCL_Debug(GetCheminfoProperties().calc_MolWeight->SumOverObject( molecule)( 0));
//      BCL_Debug(GetCheminfoProperties().calc_XLogP->SumOverObject( molecule)( 0));
//      BCL_Debug(GetCheminfoProperties().calc_TopologicalPolarSurfaceArea->SumOverObject( molecule)( 0));
//      BCL_Debug(GetCheminfoProperties().calc_HbondAcceptor->SumOverObject( molecule)( 0) + GetCheminfoProperties().calc_HbondDonor->SumOverObject( molecule)( 0));
//      BCL_Debug(GetCheminfoProperties().calc_NRotBond->SumOverObject( molecule)( 0));
//      BCL_Debug(bonde->SumOverObject( molecule)( 3));
//      BCL_Debug(GetCheminfoProperties().calc_IsF->SumOverObject( molecule)( 0));
//      BCL_Debug(GetCheminfoProperties().calc_IsCl->SumOverObject( molecule)( 0));
//      BCL_Debug(GetCheminfoProperties().calc_IsBr->SumOverObject( molecule)( 0));
//      BCL_Debug(GetCheminfoProperties().calc_IsI->SumOverObject( molecule)( 0));
//      BCL_Debug(GetCheminfoProperties().calc_IsHalogen->SumOverObject( molecule)( 0));
//      BCL_Debug(ContainsBadNitrogenLinkages( molecule));
//      BCL_Debug(CountRingHalogens( molecule));
//      BCL_Debug(CountNonAroClBrI( molecule));
//      BCL_Debug(ContainsReactiveAlkene( molecule));
//      BCL_Debug(ContainsAromaticAmine( molecule));
//      BCL_Debug(MaxRingFragmentSize( molecule));
//      BCL_Debug(!HasGoodRingAtomHalogenRatio( molecule));
//      BCL_Debug(GetCheminfoProperties().calc_MolComplexity->SumOverObject( molecule)( 0));

      // check for violations in any of the standard druglikeness criteria
      if
      (
          molecule.GetNumberAtoms() < size_t( 2) ||
          GetCheminfoProperties().calc_MolWeight->SumOverObject( molecule)( 0) > m_MaxWeight  ||
          GetCheminfoProperties().calc_MolWeight->SumOverObject( molecule)( 0) < m_MinWeight  ||
          GetCheminfoProperties().calc_XLogP->SumOverObject( molecule)( 0) < m_MinLogP ||
          GetCheminfoProperties().calc_XLogP->SumOverObject( molecule)( 0) > m_MaxLogP ||
          GetCheminfoProperties().calc_TopologicalPolarSurfaceArea->SumOverObject( molecule)( 0) > m_MaxTPSA ||
          GetCheminfoProperties().calc_HbondAcceptor->SumOverObject( molecule)( 0) + GetCheminfoProperties().calc_HbondDonor->SumOverObject( molecule)( 0) > m_MaxHBDA ||
          GetCheminfoProperties().calc_NRotBond->SumOverObject( molecule)( 0) > m_MaxNRotBonds ||
          bonde->SumOverObject( molecule)( 3) > m_MaxTotalBondEnergy || // finds worst bond propensity score of all the bonds in the molecule
          GetCheminfoProperties().calc_IsF->SumOverObject( molecule)( 0) > m_MaxNF ||
          GetCheminfoProperties().calc_IsCl->SumOverObject( molecule)( 0) > m_MaxNCl  ||
          GetCheminfoProperties().calc_IsBr->SumOverObject( molecule)( 0) > m_MaxNBr  ||
          GetCheminfoProperties().calc_IsI->SumOverObject( molecule)( 0) > m_MaxNI ||
          GetCheminfoProperties().calc_IsHalogen->SumOverObject( molecule)( 0) > m_MaxNHalogens ||
          ContainsBadNitrogenLinkages( molecule) ||
          CountRingHalogens( molecule) > size_t( m_MaxRingHalogens) ||
          CountNonAroClBrI( molecule) > size_t( m_MaxNonAromaticClBrI) ||
//          ContainsReactiveAlkene( molecule) ||
//          ContainsAromaticAmine( molecule) ||
          MaxRingFragmentSize( molecule) > m_MaxRingSize ||
          !HasGoodRingAtomHalogenRatio( molecule) ||
          GetCheminfoProperties().calc_MolComplexity->SumOverObject( molecule)( 0) > m_MaxComplexity
      )
        // one of the above violations occurred
      {
        BCL_MessageVrb( "Failed at least 1 criteria from druglikeness set");
        STORAGE( 0) = 0;
      }
      // no violations
      else
      {
        STORAGE( 0) = 1;
      }
    } // Recalculate

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief calculate the max number of halogen atoms on a ring in a molecule
    //! @param MOLECULE the molecule of interest
    size_t MoleculeDruglike::CountRingHalogens( const chemistry::FragmentComplete &MOLECULE) const
    {
      // generate a molecule splitter
      chemistry::FragmentSplitRigid splitter;
      chemistry::FragmentEnsemble rigid_fragments( splitter( MOLECULE));

      // see if fragment is a ring
      size_t max_n_halogens( 0);
      for
      (
          auto frag_itr( rigid_fragments.Begin()), frag_itr_end( rigid_fragments.End());
          frag_itr != frag_itr_end;
          ++frag_itr
      )
      {
        if( GetCheminfoProperties().calc_NRings->SumOverObject( *frag_itr)( 0))
        {
          // if so, count number of halogens
          size_t n_halogens( GetCheminfoProperties().calc_IsHalogen->SumOverObject( *frag_itr)( 0));
          if( n_halogens > max_n_halogens)
          {
            max_n_halogens = n_halogens;
          }
        }
      }
      return max_n_halogens;
    }

    //! @brief determines if the ratio of ring heavy atoms to substituted halogens is druglike
    //! @param MOLECULE the molecule of interest
    bool MoleculeDruglike::HasGoodRingAtomHalogenRatio( const chemistry::FragmentComplete &MOLECULE) const
    {
      // yet another attempt to get make ring substitutions more reasonable
      // generate a molecule splitter
      chemistry::FragmentSplitRigid rigid_splitter;
      chemistry::FragmentSplitRings ring_splitter( true, size_t( 3));
      chemistry::FragmentEnsemble rigid_fragments( rigid_splitter( MOLECULE));

      // see if fragment is a ring
      for
      (
          auto rigid_frag_itr( rigid_fragments.Begin()), rigid_frag_itr_end( rigid_fragments.End());
          rigid_frag_itr != rigid_frag_itr_end;
          ++rigid_frag_itr
      )
      {
        // check if rigid fragment contains ring
        if( GetCheminfoProperties().calc_NRings->SumOverObject( *rigid_frag_itr)( 0))
        {
          // if so, count number of halogens attached to ring atoms
//          size_t n_halogens( GetCheminfoProperties().calc_IsHalogen->SumOverObject( *rigid_frag_itr)( 0));
          size_t n_halogens( 0);
          for
          (
              auto atom_itr( rigid_frag_itr->GetAtomVector().Begin()), atom_itr_end( rigid_frag_itr->GetAtomVector().End());
              atom_itr != atom_itr_end;
              ++atom_itr
          )
          {
            // checking if atom is in ring
            if
            (
                atom_itr->CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsInRing, size_t( 1))
            )
            {
              // if so, see if bonded to a halogen
              for
              (
                 auto bond_itr( atom_itr->GetBonds().Begin()), bond_itr_end( atom_itr->GetBonds().End());
                  bond_itr != bond_itr_end;
                  ++bond_itr
              )
              {
                if
                (
                    bond_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Fluorine ||
                    bond_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Chlorine ||
                    bond_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Bromine ||
                    bond_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Iodine
                )
                {
                  n_halogens += size_t( 1);
                  break;
                }
              }
            }
          }

          // no point in continuing this fragment if no halogens present
          if( n_halogens)
          {
            // count number of atoms in ring
            chemistry::FragmentEnsemble ring( ring_splitter( *rigid_frag_itr));
            ring.RemoveH();
            if( ring.GetSize())
            {
              // do not allow rings to simultaneously have 6 or fewer atoms and 3 or more halogens
              size_t n_ring_atoms( ring.GetMolecules().FirstElement().GetSize());
              if( n_ring_atoms < size_t( 7) && n_halogens > size_t( 2))
              {
                return false;
              }
            }
            else
            {
              continue;
            }
          }
        }
      }
      return true;
    }

    //! @brief calculate the max number of atoms in a ring in a molecule
    //! @param MOLECULE the molecule of interest
    size_t MoleculeDruglike::MaxRingFragmentSize( const chemistry::FragmentComplete &MOLECULE) const
    {
      // generate a molecule splitter
      chemistry::FragmentSplitRings splitter( true, 3);
      chemistry::FragmentEnsemble ring_fragments( splitter( MOLECULE));

      // see if fragment is a ring
      size_t max_n_ring_atoms( 0);
      for
      (
          auto frag_itr( ring_fragments.Begin()), frag_itr_end( ring_fragments.End());
          frag_itr != frag_itr_end;
          ++frag_itr
      )
      {
        // if so, count number of halogens
        size_t n_ring_atoms( frag_itr->GetNumberAtoms() - frag_itr->GetNumberHydrogens());
        if( n_ring_atoms > max_n_ring_atoms)
        {
          max_n_ring_atoms = n_ring_atoms;
        }
      }
      return max_n_ring_atoms;
    }

    //! @brief calculate the total number of non-aromatic ring Cl, Br, I halogen atoms
    //! @param MOLECULE the molecule of interest
    size_t MoleculeDruglike::CountNonAroClBrI( const chemistry::FragmentComplete &MOLECULE) const
    {
      // generate a molecule splitter
      chemistry::FragmentSplitRigid splitter;
      chemistry::FragmentEnsemble rigid_fragments( splitter( MOLECULE));

      // see if fragment is an aromatic ring
      size_t n_halogens( 0);
      for
      (
          auto frag_itr( rigid_fragments.Begin()), frag_itr_end( rigid_fragments.End());
          frag_itr != frag_itr_end;
          ++frag_itr
      )
      {
        // find non-aromatic rings
        if( !GetCheminfoProperties().calc_NAromaticRings->SumOverObject( *frag_itr)( 0))
        {
          // count number of cl br i halogens
          n_halogens += GetCheminfoProperties().calc_IsCl->SumOverObject( *frag_itr)( 0);
          n_halogens += GetCheminfoProperties().calc_IsBr->SumOverObject( *frag_itr)( 0);
          n_halogens += GetCheminfoProperties().calc_IsI->SumOverObject( *frag_itr)( 0);
        }
      }
      return n_halogens;
    }

    //! @brief determine whether molecule rings are drug-like
    //! @param MOLECULE the molecule of interest
    bool MoleculeDruglike::ContainsOnlyDruglikeRings( const chemistry::FragmentComplete &MOLECULE) const
    {
      // read in ring database
      static std::string ring_database( chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "")
        + "ring_libraries/drug_ring_database.sdf.gz");
      BCL_Assert( !ring_database.empty(), "could not find drug ring database!");
      io::IFStream input;
      io::File::MustOpenIFStream( input, ring_database);

      // remove hydrogen atoms to make faster
      static chemistry::FragmentEnsemble ensemble( input, sdf::e_Remove);
      io::File::CloseClearFStream( input);

      // collect rings in object
      static util::ShPtr< chemistry::ConstitutionSet> rings( new chemistry::ConstitutionSet());
      for
      (
          auto itr_ensemble( ensemble.Begin()), itr_ensemble_end( ensemble.End());
          itr_ensemble != itr_ensemble_end;
          ++itr_ensemble
      )
      {
        rings->Insert( chemistry::FragmentConstitutionShared( *itr_ensemble));
      }
      BCL_Assert( rings->GetSize(), "Ring library is empty!");

      // split out rings
      static chemistry::FragmentSplitRings ring_splitter( true, 3);
      chemistry::FragmentEnsemble split_rings( ring_splitter( MOLECULE));

      // make sure all rings are found in the ring dataset
      for
      (
          chemistry::FragmentEnsemble::iterator split_rings_itr( split_rings.Begin()), split_rings_itr_end( split_rings.End());
          split_rings_itr != split_rings_itr_end;
          ++split_rings_itr
      )
      {
        // remove hydrogen atoms to be consistent with how we read in ring database
        split_rings_itr->RemoveH();
        if( rings->Find( chemistry::FragmentConstitutionShared( *split_rings_itr)) == rings->End())
        {
          // we have a ring not in the database; we fail the filter
          return false;
        }
      }
      return true;
    }

    //! @brief determine if undesirable nitrogen covalent bonds exist in a molecule
    //! @param MOLECULE the molecule of interest
    bool MoleculeDruglike::ContainsBadNitrogenLinkages( const chemistry::FragmentComplete &MOLECULE) const
    {
      // iterate over all atoms in molecules
      for
      (
          auto atoms_itr( MOLECULE.GetAtomVector().Begin()), atoms_itr_end( MOLECULE.GetAtomVector().End());
          atoms_itr != atoms_itr_end;
          ++atoms_itr
      )
      {
        // if nitrogen, find bonded atoms
        if( atoms_itr->GetElementType() == chemistry::GetElementTypes().e_Nitrogen)
        {
          for
          (
              auto bonds_itr( atoms_itr->GetBonds().Begin()), bonds_itr_end( atoms_itr->GetBonds().End());
              bonds_itr != bonds_itr_end;
              ++bonds_itr
          )
          {
            // check to see if bonded to some atom type that is not okay outside of a ring
            if
            (
              // not in a ring
              !bonds_itr->GetTargetAtom().CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsInRing, size_t( 1))
              &&
              (
                // is bonded to O outside of a ring
                bonds_itr->GetTargetAtom().GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Oxygen ||
                // is double bonded to anything outside of a ring
                bonds_itr->GetBondType() == chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond ||
                // is bonded to a S and is part of an amide bond
                (
                  atoms_itr->CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsAmide, size_t( 1))
                  && bonds_itr->GetTargetAtom().GetAtomType()->GetElementType() == chemistry::GetElementTypes().e_Sulfur
                )
              )
            )
            {
              return true;
            }
          }
        }
        // if sulfur, find bond types
        else if( atoms_itr->GetElementType() == chemistry::GetElementTypes().e_Sulfur)
        {
          // go over bonds
          for
          (
              auto bond_itr( atoms_itr->GetBonds().Begin()), bond_itr_end( atoms_itr->GetBonds().End());
              bond_itr != bond_itr_end;
              ++bond_itr
          )
          {
            // if double bond, only allow oxygen as partner
            if
            (
                bond_itr->GetBondType() == chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond &&
                bond_itr->GetTargetAtom().GetElementType() != chemistry::GetElementTypes().e_Oxygen
            )
            {
              return true;
            }
            // if single bond, do not allow oxygen as partner
            if
            (
                bond_itr->GetBondType() != chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond &&
                bond_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Oxygen
            )
            {
              return true;
            }
          }
        }
      }
      return false;
    }

    //! @brief determine if molecule contains aromatic amine
    //! @param MOLECULE the molecule of interest
    bool MoleculeDruglike::ContainsAromaticAmine( const chemistry::FragmentComplete &MOLECULE) const
    {
      // generate a molecule splitter
      chemistry::FragmentSplitRigid splitter;
      chemistry::FragmentEnsemble rigid_fragments( splitter( MOLECULE));

      // iterate over fragments from splitter
      for
      (
          auto frag_itr( rigid_fragments.Begin()), frag_itr_end( rigid_fragments.End());
          frag_itr != frag_itr_end;
          ++frag_itr
      )
      {
        // see if fragment is an aromatic ring
        if( GetCheminfoProperties().calc_NAromaticRings->SumOverObject( *frag_itr)( 0))
        {
          // if so, iterate over atoms
          for
          (
              auto atom_itr( frag_itr->GetAtomVector().Begin()), atom_itr_end( frag_itr->GetAtomVector().End());
              atom_itr != atom_itr_end;
              ++atom_itr
          )
          {
            // only care if it is a nitrogen and it is not in a ring
            if
            (
                atom_itr->GetElementType() == chemistry::GetElementTypes().e_Nitrogen &&
                !atom_itr->CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsInRing, size_t( 1))
            )
            {
              // go over bonds from this terrible nitrogen
              size_t n_h( 0);
              bool aro_friend( false);
              for
              (
                  auto bond_itr( atom_itr->GetBonds().Begin()), bond_itr_end( atom_itr->GetBonds().End());
                  bond_itr != bond_itr_end;
                  ++bond_itr
              )
              {
                // if one of the bonds is to an aromatic
                if( bond_itr->GetTargetAtom().CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsAromatic, size_t( 1)))
                {
                  aro_friend = true;
                }
                else if( bond_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
                {
                  ++n_h;
                }
              }
              // check to see if it is an aromatic amine
              if( aro_friend && n_h >= size_t( 2))
              {
                return true;
              }
            }
            // I guess also care if it is an oxygen and not in a ring
            else if
            (
                atom_itr->GetElementType() == chemistry::GetElementTypes().e_Oxygen &&
                !atom_itr->CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsInRing, size_t( 1))
            )
            {

              // go over bonds from this terrible oxygen
              size_t n_h( 0);
              bool aro_friend( false);
              for
              (
                  auto bond_itr( atom_itr->GetBonds().Begin()), bond_itr_end( atom_itr->GetBonds().End());
                  bond_itr != bond_itr_end;
                  ++bond_itr
              )
              {
                // if one of the bonds is to an aromatic
                if( bond_itr->GetTargetAtom().CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsAromatic, size_t( 1)))
                {
                  aro_friend = true;
                }
                else if( bond_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
                {
                  ++n_h;
                }
              }
              // check to see if it is an aromatic alcohol
              if( aro_friend && n_h >= size_t( 1))
              {
                return true;
              }
            }
          }
        }
      }
      return false;
    }

    //! @brief determine if molecules contains a reactive alkene
    //! note that this will need to be disabled if we choose to include certain covalent warheads, e.g. acrylamide
    //! @param MOLECULE the molecule of interest
    bool MoleculeDruglike::ContainsReactiveAlkene( const chemistry::FragmentComplete &MOLECULE) const
    {
      // reactive alkenes can be good if we want to form covalent bonds in a controlled manner
      // however, it is probably usually bad
      // here, we consider any conjugated double bond C = C to be reactive if it is not in a ring

      // go over the bondinfo vector and find atom index pairs where there is a conjugated double bond
      size_t bondinfo_index( 0), bondinfo_sz( MOLECULE.GetAtomVector().GetBondInfo().GetSize());
      for
      (
          auto bondinfo_itr( MOLECULE.GetAtomVector().GetBondInfo().Begin());
          bondinfo_index < bondinfo_sz;
          ++bondinfo_itr, ++bondinfo_index
      )
      {
        // skip it if it is not a conjugated double bond
        if( std::max( bondinfo_itr->GetAtomIndexHigh(), bondinfo_itr->GetAtomIndexLow()) < MOLECULE.GetSize())
        {
          if( bondinfo_itr->GetConfigurationalBondType() != chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond)
          {
            continue;
          }
          // if the two atom indices are carbon atoms then check if either are in a ring
          if
          (
              MOLECULE.GetAtomVector()( bondinfo_itr->GetAtomIndexHigh()).GetElementType() == chemistry::GetElementTypes().e_Carbon &&
              MOLECULE.GetAtomVector()( bondinfo_itr->GetAtomIndexLow()).GetElementType() == chemistry::GetElementTypes().e_Carbon
          )
          {
            // if they are not in a ring then fail
            if
            (
                !MOLECULE.GetAtomVector()( bondinfo_itr->GetAtomIndexHigh()).CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsInRing, size_t( 1)) ||
                !MOLECULE.GetAtomVector()( bondinfo_itr->GetAtomIndexLow()).CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsInRing, size_t( 1))
            )
            {
              return true;
            }
          }
        }
      }
      return false;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeDruglike::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Returns 1 if the molecule is druglike and 0 otherwise");
      parameters.AddInitializer
      (
        "min_weight",
        "minimum molecular weight considered to be druglike",
        io::Serialization::GetAgent( &m_MinWeight),
        "10.0"
      );
      parameters.AddInitializer
      (
        "max_weight",
        "maximum molecular weight considered to be druglike",
        io::Serialization::GetAgent( &m_MaxWeight),
        "600.0"
      );
      parameters.AddInitializer
      (
        "max_ring_size",
        "maximum number of heavy atoms in a ring considered to be druglike",
        io::Serialization::GetAgent( &m_MaxRingSize),
        "10"
      );
      parameters.AddInitializer
      (
        "min_log_p",
        "minimum logP considered to be druglike",
        io::Serialization::GetAgent( &m_MinLogP),
        "-1.0"
      );
      parameters.AddInitializer
      (
        "max_log_p",
        "maximum logP considered to be druglike",
        io::Serialization::GetAgent( &m_MaxLogP),
        "5.0"
      );
      parameters.AddInitializer
      (
        "max_tpsa",
        "maximum topological polar surface area considered to be druglike",
        io::Serialization::GetAgent( &m_MaxTPSA),
        "140"
      );
      parameters.AddInitializer
      (
        "max_hbda",
        "maximum number of hydrogen bond donors + acceptors considered to be druglike",
        io::Serialization::GetAgent( &m_MaxHBDA),
        "12"
      );
      parameters.AddInitializer
      (
        "max_n_rot_bonds",
        "maximum number of rotatable bonds considered to be druglike",
        io::Serialization::GetAgent( &m_MaxNRotBonds),
        "10"
      );
      parameters.AddInitializer
      (
        "max_bond_score",
        "maximum bond propensity score considered to be druglike",
        io::Serialization::GetAgent( &m_MaxTotalBondEnergy),
        "0.10"
      );
      parameters.AddInitializer
      (
        "max_complexity",
        "maximum molecule complexity considered to be druglike",
        io::Serialization::GetAgent( &m_MaxComplexity),
        "4.0"
      );
      parameters.AddInitializer
      (
        "max_n_f",
        "maximum number of fluorine atoms considered to be druglike",
        io::Serialization::GetAgent( &m_MaxNF),
        "5"
      );
      parameters.AddInitializer
      (
        "max_n_cl",
        "maximum number of chlorine atoms considered to be druglike",
        io::Serialization::GetAgent( &m_MaxNCl),
        "2"
      );
      parameters.AddInitializer
      (
        "max_n_br",
        "maximum number of bromine atoms considered to be druglike",
        io::Serialization::GetAgent( &m_MaxNBr),
        "1"
      );
      parameters.AddInitializer
      (
        "max_n_i",
        "maximum number of iodine atoms considered to be druglike",
        io::Serialization::GetAgent( &m_MaxNI),
        "1"
      );
      parameters.AddInitializer
      (
        "max_total_halogens",
        "maximum total number of halogen atoms considered to be druglike",
        io::Serialization::GetAgent( &m_MaxNHalogens),
        "5"
      );
      parameters.AddInitializer
      (
        "max_n_ring_halogens",
        "maximum total number of halogen atoms that can be on an aromatic ring system",
        io::Serialization::GetAgent( &m_MaxRingHalogens),
        "3"
      );
      parameters.AddInitializer
      (
        "max_n_non_aromatic_cl_br_i",
        "maximum total number of Cl, Br, and I atoms that can be on a non-aromatic system",
        io::Serialization::GetAgent( &m_MaxNonAromaticClBrI),
        "0"
      );
      parameters.AddOptionalInitializer
      (
        "enforce_hitlike",
        "require that the molecule also fullfil hit-like criteria;"
        "criteria based on a simple decision tree trained with the data provided"
        "in Bickerton et al., 2012, Quantifying the chemical beauty of drugs, 10.1038/nchem.1243."
        "Specifically, it addresses the following question posed to medicinal chemists:"
        "'If this compound came out of a primary screen, would you perform hit optimization on it?'",
        io::Serialization::GetAgent( &m_EnforceHitlike)
      );
      parameters.AddOptionalInitializer
      (
        "enforce_druglike_rings",
        "require that the molecule only contain rings found in the druglike ring database",
        io::Serialization::GetAgent( &m_EnforceDruglikeRings)
      );
      parameters.AddOptionalInitializer
      (
        "enforce_druglike_fragments",
        "require that the molecule only be composed of ECFP fragments found in the ECFP fragment database",
        io::Serialization::GetAgent( &m_EnforceDruglikeECFP)
      );

      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
