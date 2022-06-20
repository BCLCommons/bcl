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
#include "chemistry/bcl_chemistry_fragment_mutate_halogenate.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "find/bcl_find_collector_interface.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "random/bcl_random_uniform_distribution.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentMutateHalogenate::s_Instance
    (
      util::Enumerated< FragmentMutateInterface>::AddInstance( new FragmentMutateHalogenate())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentMutateHalogenate::FragmentMutateHalogenate() :
      m_AllowedHalogens( storage::Vector< AtomType>()),
      m_AllowedHalogensString( "F Cl Br I"),
      m_Reversible( false),
      m_RestrictReversibility( false),
      m_DisableAromaticRingReq( false)
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief druglikeness constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    FragmentMutateHalogenate::FragmentMutateHalogenate
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const bool &CORINA_CONFS
    ) :
      m_AllowedHalogens( storage::Vector< AtomType>()),
      m_AllowedHalogensString( "F Cl Br I"),
      m_Reversible( false),
      m_RestrictReversibility( false),
      m_DisableAromaticRingReq( false)
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief local mutate constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    FragmentMutateHalogenate::FragmentMutateHalogenate
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const bool &CORINA_CONFS
    ) :
      m_AllowedHalogens( storage::Vector< AtomType>()),
      m_AllowedHalogensString( "F Cl Br I"),
      m_Reversible( false),
      m_RestrictReversibility( false),
      m_DisableAromaticRingReq( false)
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_ScaffoldFragment = SCAFFOLD_FRAGMENT;
      m_MutableFragments = MUTABLE_FRAGMENTS;
      m_MutableAtomIndices = MUTABLE_ATOM_INDICES;
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief local mutate pose-sensitive constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    //! @param MDL property label containing path to protein binding pocket PDB file
    //! @param PROPERTY_SCORER property that will be used to score interactions with protein pocket
    //! @param RESOLVE_CLASHES if true, resolve clashes with specified protein pocket after mutatation
    //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
    FragmentMutateHalogenate::FragmentMutateHalogenate
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const std::string &MDL,
      const descriptor::CheminfoProperty &PROPERTY_SCORER,
      const bool &RESOLVE_CLASHES,
      const storage::Vector< float> &BFACTORS,
      const bool &CORINA_CONFS
    ) :
      m_AllowedHalogens( storage::Vector< AtomType>()),
      m_AllowedHalogensString( "F Cl Br I"),
      m_Reversible( false),
      m_RestrictReversibility( false),
      m_DisableAromaticRingReq( false)
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_ScaffoldFragment = SCAFFOLD_FRAGMENT;
      m_MutableFragments = MUTABLE_FRAGMENTS;
      m_MutableAtomIndices = MUTABLE_ATOM_INDICES;
      m_MDL = MDL;
      m_PropertyScorer = PROPERTY_SCORER;
      m_ResolveClashes = RESOLVE_CLASHES;
      m_BFactors = BFACTORS;
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief local clash resolver constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    //! @param MDL property label containing path to protein binding pocket PDB file
    //! @param PROPERTY_SCORER property that will be used to score interactions with protein pocket
    //! @param RESOLVE_CLASHES if true, resolve clashes with specified protein pocket after mutatation
    //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
    FragmentMutateHalogenate::FragmentMutateHalogenate
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const std::string &MDL,
      const bool &RESOLVE_CLASHES,
      const storage::Vector< float> &BFACTORS,
      const bool &CORINA_CONFS
    ) :
      m_AllowedHalogens( storage::Vector< AtomType>()),
      m_AllowedHalogensString( "F Cl Br I"),
      m_Reversible( false),
      m_RestrictReversibility( false),
      m_DisableAromaticRingReq( false)
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_ScaffoldFragment = SCAFFOLD_FRAGMENT;
      m_MutableFragments = MUTABLE_FRAGMENTS;
      m_MutableAtomIndices = MUTABLE_ATOM_INDICES;
      m_MDL = MDL;
      m_ResolveClashes = RESOLVE_CLASHES;
      m_BFactors = BFACTORS;
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief clone constructor
    FragmentMutateHalogenate *FragmentMutateHalogenate::Clone() const
    {
      return new FragmentMutateHalogenate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentMutateHalogenate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentMutateHalogenate::GetAlias() const
    {
      static const std::string s_name( "Halogenate");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @return MutateResult with Constitution after the mutate
    math::MutateResult< FragmentComplete> FragmentMutateHalogenate::operator()( const FragmentComplete &FRAGMENT) const
    {
      BCL_MessageStd( "Halogenate!");

      // for cleaning and optimizing the new molecule conformer
      FragmentMapConformer cleaner
      (
        m_DrugLikenessType,
        m_MDL,
        FRAGMENT.GetMDLProperty( m_MDL),
        m_PropertyScorer,
        m_ResolveClashes,
        m_BFactors,
        m_Corina
      );

      for( size_t i( 0); i < m_NumberMaxAttempts; ++i)
      {
        // pick an atom
        util::SiPtr< const AtomConformationalInterface> picked_atom;
        if( m_MutableAtomIndices.GetSize() || m_MutableElements.GetSize() || m_MutableFragments.GetSize())
        {
          picked_atom = this->PickAtom( FRAGMENT, false);
        }
        else
        {
          picked_atom = this->PickAtom( FRAGMENT, true);
        }

        // if atom is hydrogen atom, grab the atom to which it is connected
        if( picked_atom->GetElementType() == GetElementTypes().e_Hydrogen)
        {
          if( !picked_atom->GetBonds().GetSize())
          {
            continue;
          }
          picked_atom = util::SiPtr< const AtomConformationalInterface>( picked_atom->GetBonds().Begin()->GetTargetAtom());
        }

        // must pass criteria
        bool pass( false);
        m_DisableAromaticRingReq ?
            pass = bool( picked_atom->GetElementType() == GetElementTypes().e_Carbon) : // only add halogens to carbon
            pass = bool( picked_atom->GetElementType() == GetElementTypes().e_Carbon && // only add halogens to carbon AND
            picked_atom->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1)) && // must be in a ring AND
            picked_atom->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAromatic, size_t( 1))); // must be aromatic

        // check if the hybridization could be aromatic
        if( pass)
        {
          // get index of picked atom
          size_t picked_atom_index( FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom));

          // manage reversibility
          float reversible_rand( 0.0);
          if( m_Reversible)
          {
            reversible_rand = random::GetGlobalRandom().Random< float>( 0.0, 1.0);
          }
          // TODO: bool
          // check if substitution position is accessible
          //            if( !IsRingSubstitutionDirected( FRAGMENT, picked_atom))
          //            {
          //              // if aromatic and not correctly substituted, skip mutate
          //              return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
          //            }

          // count number of allowed halogens and probabilities
          size_t n_allowed_halogens( m_AllowedHalogens.GetSize());
          storage::Vector< float> halogen_probs( n_allowed_halogens, 1.0 / n_allowed_halogens);

          // randomly pick a halogen
          float rand( random::GetGlobalRandom().Random< float>( 0.0, 1.0));
          float cumulative_weighted_sum( 0.0);
          size_t rand_halogen_i( 0);
          for
          (
              ; // we want the halogen index to persist after the loop
              rand_halogen_i < n_allowed_halogens;
              ++rand_halogen_i
          )
          {
            cumulative_weighted_sum += halogen_probs( rand_halogen_i);
            if( rand < cumulative_weighted_sum)
            {
              break;
            }
          }

          // non-reversible
          if( reversible_rand < 0.5)
          {
            // make a fragment object from the chosen halogen
            // i know it seems more intuitive to just do SetAtomType on a bonded hydrogen atom (original algorithm),
            // but using that approach it periodically gets hung up really bad when many mutations are
            // being performed. my hope is that refocusing it to a heavy-atom centric approach will
            // solve that issue, but i know this is ugly and looks like massive overkill.
            sdf::AtomInfo new_halogen( m_AllowedHalogens( rand_halogen_i), e_NonChiral);
            storage::Vector< sdf::AtomInfo> new_halogen_atominfo( size_t( 1), new_halogen);
            storage::Vector< sdf::BondInfo> empty_bonds;
            AtomVector< AtomComplete> new_halogen_v( new_halogen_atominfo, empty_bonds);
            FragmentComplete new_halogen_frag( new_halogen_v, "");

            // open a valence on the original molecule
            storage::Triplet< FragmentComplete, size_t, size_t> pair_a( OpenValence( FRAGMENT, picked_atom_index, m_OVShuffleH, m_OVReverse));

            // link fragments
            storage::Pair< bool, FragmentComplete> new_fragment
            (
              MergeFragmentComplete::MergeFragments
              (
                pair_a.First(),
                new_halogen_frag,
                GetConfigurationalBondTypes().e_ConjugatedSingleBond,
                storage::Pair< size_t, size_t>( pair_a.Second(), size_t( 0))
              )
            );

            // done
            if( new_fragment.First())
            {
              AtomVector< AtomComplete> atom_vector( new_fragment.Second().GetAtomVector());
              HydrogensHandler::Remove( atom_vector);

              util::ShPtr< FragmentComplete> new_mol_ptr
              (
                m_ScaffoldFragment.GetSize()
                ? cleaner.Clean( atom_vector, m_ScaffoldFragment, m_DrugLikenessType)
                    : cleaner.Clean( atom_vector, FRAGMENT, m_DrugLikenessType)
              );
              return math::MutateResult< FragmentComplete>( new_mol_ptr, *this);
            }
          }
          // reversible
          else
          {
            // get a modifiable version of the atom vector
            AtomVector< AtomComplete> frag_atom_v( FRAGMENT.GetAtomVector());
            size_t removal_index( util::GetUndefinedSize_t());

            // iterate over all bonded partners to find a halogen
            for
            (
                auto bond_itr( frag_atom_v( picked_atom_index).GetBonds().Begin()),
                bond_itr_end( frag_atom_v( picked_atom_index).GetBonds().End());
                bond_itr != bond_itr_end;
                ++bond_itr
            )
            {
              if( m_RestrictReversibility)
              {
                if( bond_itr->GetTargetAtom().GetElementType() == m_AllowedHalogens( rand_halogen_i)->GetElementType())
                {
                  removal_index = frag_atom_v.GetAtomIndex( bond_itr->GetTargetAtom());
                  break;
                }
              }
              else
              {
                if
                (
                    bond_itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Fluorine ||
                    bond_itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Chlorine ||
                    bond_itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Bromine ||
                    bond_itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Iodine
                )
                {
                  removal_index = frag_atom_v.GetAtomIndex( bond_itr->GetTargetAtom());
                  break;
                }
              }
            }
            // remove identified halogen
            if( util::IsDefined( removal_index))
            {
              storage::Vector< size_t> atom_indices( storage::CreateIndexVector( frag_atom_v.GetSize()));
              atom_indices.RemoveElements( removal_index, size_t( 1));
              frag_atom_v.Reorder( atom_indices);

              // make new mol
              HydrogensHandler::Remove( frag_atom_v);
              util::ShPtr< FragmentComplete> new_mol_ptr
              (
                m_ScaffoldFragment.GetSize()
                ? cleaner.Clean( frag_atom_v, m_ScaffoldFragment, m_DrugLikenessType)
                    : cleaner.Clean( frag_atom_v, FRAGMENT, m_DrugLikenessType)
              );
              return math::MutateResult< FragmentComplete>( new_mol_ptr, *this);
            }
            return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
          } // end reversibility handling
        } // end hybridization check for aromaticity
      } // end redo attempts if it fails
      // if we never get it right then just return null
      return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set the fragment mutable atom indices
    void FragmentMutateHalogenate::SetAllowedHalogens( const storage::Vector< AtomType> &ALLOWED_HALOGENS)
    {
      m_AllowedHalogens = ALLOWED_HALOGENS;
    }

    //! @brief set reversibility
    void FragmentMutateHalogenate::SetReversibility( const bool REVERSIBLE)
    {
      m_Reversible = REVERSIBLE;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    io::Serializer FragmentMutateHalogenate::GetSerializer() const
    {
      io::Serializer parameters( FragmentMutateInterface::GetSerializer());
      parameters.SetClassDescription
      (
        "Extends a molecule either by breaking a bond and inserting a motif before reconnecting, "
        "or by linking to a new ring system"
      );

      parameters.AddInitializer
      (
        "allowed_halogens",
        "halogens that are accessible via this mutate",
        io::Serialization::GetAgent( &m_AllowedHalogensString),
        "F Cl Br I"
      );

      parameters.AddInitializer
      (
        "reversible",
        "sets the probability of removing a halogen to 50%; "
        "default false indicates that halogens can only be added",
        io::Serialization::GetAgent( &m_Reversible),
        "false"
      );

      parameters.AddInitializer
      (
        "restrict_reversibility",
        "if enabled with the 'reversible' flag then only halogens specified with 'allowed_halogens' "
        "flag are eligible for removal; maybe useful when wanting to design halogen sites without "
        "perturbing existing halogens",
        io::Serialization::GetAgent( &m_RestrictReversibility),
        "false"
      );

      parameters.AddInitializer
      (
        "disable_aromatic_ring_req",
        "halogenate by default only adds halogens to aromatic rings; this behavior can be overridden "
        "by setting this flag to true",
        io::Serialization::GetAgent( &m_DisableAromaticRingReq),
        "false"
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentMutateHalogenate::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      // static initialization check
      if( command::CommandState::IsInStaticInitialization())
      {
        return true;
      }

      // call RISH function of the base class
      if( !FragmentMutateInterface::ReadInitializerSuccessHook( LABEL, ERROR_STREAM))
      {
        return false;
      }

      // read in allowed halogens
      if( m_AllowedHalogensString.size())
      {
        // parse input
        const storage::Vector< std::string> allowed_halogens
        (
          util::SplitString( util::TrimString( m_AllowedHalogensString), " \t\n\r,")
        );

        // stupid check to add only the correct halogens
        for( size_t h_i( 0), h_sz( allowed_halogens.GetSize()); h_i < h_sz; ++h_i)
        {
          // Fluorine
          if( allowed_halogens( h_i) == "F")
          {
            m_AllowedHalogens.PushBack( GetAtomTypes().F_SP2P2P2);
          }
          // Chlorine
          if( allowed_halogens( h_i) == "Cl")
          {
            m_AllowedHalogens.PushBack( GetAtomTypes().Cl_SP2P2P2);
          }
          // Bromine
          if( allowed_halogens( h_i) == "Br")
          {
            m_AllowedHalogens.PushBack( GetAtomTypes().Br_SP2P2P2);
          }
          // Iodine
          if( allowed_halogens( h_i) == "I")
          {
            m_AllowedHalogens.PushBack( GetAtomTypes().I_SP2P2P2);
          }
        }
      }

      // done
      return true;
    }

  } // namespace chemistry
} // namespace bcl
