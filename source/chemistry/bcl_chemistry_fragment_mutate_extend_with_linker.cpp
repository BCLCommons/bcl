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
#include "chemistry/bcl_chemistry_fragment_mutate_extend_with_linker.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_make_conformers.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_mutate_remove_bond.h"
#include "chemistry/bcl_chemistry_fragment_split_isolate.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"
#include "chemistry/bcl_chemistry_mutate_bond_angles.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "command/bcl_command_command_state.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "find/bcl_find_collector_interface.h"
#include "graph/bcl_graph_edge_cover_ring_perception.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "sdf/bcl_sdf_bond_info.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentMutateExtendWithLinker::s_Instance
    (
      util::Enumerated< FragmentMutateInterface>::AddInstance( new FragmentMutateExtendWithLinker())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentMutateExtendWithLinker::FragmentMutateExtendWithLinker() :
        m_Rings( util::ShPtr< FragmentEnsemble>()),
        m_RingsFilename( std::string()),
        m_ExtendWithinProb( 0.50),
        m_FragmentMinSize( 0),
        m_AllowFragmentDuplication( false),
        m_MethoxyLinkProb( 0.20),
        m_MethoxyOToAProb( 0.50),
        m_EthoxyLinkProb( 0.20),
        m_EthoxyOToAProb( 0.50),
        m_AmideLinkProb( 0.20),
        m_AmideNToAProb( 0.50),
        m_EsterLinkProb( 0.20),
        m_EsterOToAProb( 0.5),
        m_SingleElementLinkProb( 0.20),
        m_B( 0.1),
        m_C( 0.1),
        m_O( 0.1),
        m_N( 0.1),
        m_P( 0.1),
        m_S( 0.1),
        m_Se( 0.1),
        m_DirectLinkProb( 0.20),
        m_AlkylLinkProb( 0.20),
        m_ThreeC( 0.0),
        m_Alkyne( 0.5),
        m_RingLinkProb( 0.20)
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    FragmentMutateExtendWithLinker::FragmentMutateExtendWithLinker
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const bool &CORINA_CONFS
    ) :
      m_Rings( util::ShPtr< FragmentEnsemble>()),
      m_RingsFilename( std::string()),
      m_ExtendWithinProb( 0.50),
      m_FragmentMinSize( 0),
      m_AllowFragmentDuplication( false),
      m_MethoxyLinkProb( 0.20),
      m_MethoxyOToAProb( 0.50),
      m_EthoxyLinkProb( 0.20),
      m_EthoxyOToAProb( 0.50),
      m_AmideLinkProb( 0.20),
      m_AmideNToAProb( 0.50),
      m_EsterLinkProb( 0.20),
      m_EsterOToAProb( 0.5),
      m_SingleElementLinkProb( 0.20),
      m_B( 0.1),
      m_C( 0.1),
      m_O( 0.1),
      m_N( 0.1),
      m_P( 0.1),
      m_S( 0.1),
      m_Se( 0.1),
      m_DirectLinkProb( 0.20),
      m_AlkylLinkProb( 0.20),
      m_ThreeC( 0.0),
      m_Alkyne( 0.5),
      m_RingLinkProb( 0.20)
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
    FragmentMutateExtendWithLinker::FragmentMutateExtendWithLinker
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
      m_Rings( util::ShPtr< FragmentEnsemble>()),
      m_RingsFilename( std::string()),
      m_ExtendWithinProb( 0.50),
      m_FragmentMinSize( 0),
      m_AllowFragmentDuplication( false),
      m_MethoxyLinkProb( 0.20),
      m_MethoxyOToAProb( 0.50),
      m_EthoxyLinkProb( 0.20),
      m_EthoxyOToAProb( 0.50),
      m_AmideLinkProb( 0.20),
      m_AmideNToAProb( 0.50),
      m_EsterLinkProb( 0.20),
      m_EsterOToAProb( 0.5),
      m_SingleElementLinkProb( 0.20),
      m_B( 0.1),
      m_C( 0.1),
      m_O( 0.1),
      m_N( 0.1),
      m_P( 0.1),
      m_S( 0.1),
      m_Se( 0.1),
      m_DirectLinkProb( 0.20),
      m_AlkylLinkProb( 0.20),
      m_ThreeC( 0.0),
      m_Alkyne( 0.5),
      m_RingLinkProb( 0.20)
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

    //! @brief local mutate pose-sensitive constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    //! @param MDL property label containing path to protein binding pocket PDB file
    //! @param PROPERTY_SCORER property that will be used to score interactions with protein pocket
    //! @param RESOLVE_CLASHES if true, resolve clashes with specified protein pocket after mutatation
    //! @param BFACTORS vector of values indicating per-residue flexibility (higher values are more flexible)
    FragmentMutateExtendWithLinker::FragmentMutateExtendWithLinker
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
      m_Rings( util::ShPtr< FragmentEnsemble>()),
      m_RingsFilename( std::string()),
      m_ExtendWithinProb( 0.50),
      m_FragmentMinSize( 0),
      m_AllowFragmentDuplication( false),
      m_MethoxyLinkProb( 0.20),
      m_MethoxyOToAProb( 0.50),
      m_EthoxyLinkProb( 0.20),
      m_EthoxyOToAProb( 0.50),
      m_EsterLinkProb( 0.20),
      m_EsterOToAProb( 0.5),
      m_AmideLinkProb( 0.20),
      m_AmideNToAProb( 0.50),
      m_SingleElementLinkProb( 0.20),
      m_B( 0.1),
      m_C( 0.1),
      m_O( 0.1),
      m_N( 0.1),
      m_P( 0.1),
      m_S( 0.1),
      m_Se( 0.1),
      m_DirectLinkProb( 0.20),
      m_AlkylLinkProb( 0.20),
      m_ThreeC( 0.0),
      m_Alkyne( 0.5),
      m_RingLinkProb( 0.20)
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
    FragmentMutateExtendWithLinker *FragmentMutateExtendWithLinker::Clone() const
    {
      return new FragmentMutateExtendWithLinker( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentMutateExtendWithLinker::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentMutateExtendWithLinker::GetAlias() const
    {
      static const std::string s_name( "ExtendWithLinker");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @return MutateResult with Constitution after the mutate
    math::MutateResult< FragmentComplete> FragmentMutateExtendWithLinker::operator()( const FragmentComplete &FRAGMENT) const
    {
      BCL_MessageStd( "ExtendWithLinker!");

      // initialize necessities
      FragmentComplete new_mol;
      util::SiPtr< const AtomConformationalInterface> picked_atom, second_atom;
      sdf::BondInfo bond_to_remove;
      size_t second_atom_index( util::GetUndefinedSize_t());
      size_t picked_atom_index( util::GetUndefinedSize_t()), tries( 0);

      // either break the molecule or extend a ring
      float decision( random::GetGlobalRandom().Random< float>( 0.0, 1.0));

      // break the molecule and extend internally with linker
      if( decision < m_ExtendWithinProb)
      {
        BCL_MessageStd( "ExtendWithLinker - Insert linker within fragment");

        // does not work if we only have one atom
        if( FRAGMENT.GetSize() == size_t( 1))
        {
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }

        // try a few times
        for( ; tries < m_NumberMaxAttempts; ++tries)
        {
          // pick an atom
          if( m_MutableAtomIndices.GetSize() || m_MutableElements.GetSize() || m_MutableFragments.GetSize())
          {
            picked_atom = this->PickAtom( FRAGMENT, false);
          }
          else
          {
            picked_atom = this->PickAtom( FRAGMENT, true);
          }

          // if a hydrogen atom is selected then just grab its bonded atom
          if( picked_atom->GetElementType() == GetElementTypes().e_Hydrogen)
          {
            // must have bonded atom
            if( !picked_atom->GetBonds().GetSize())
            {
              continue;
            }
            picked_atom = util::SiPtr< const AtomConformationalInterface>( picked_atom->GetBonds().Begin()->GetTargetAtom());
          }
          picked_atom_index = FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom);

          // pick from paired atoms
          if( m_PairedAtomIndices.GetSize())
          {
            // second atom
            second_atom_index = m_PairedAtomIndices( random::GetGlobalRandom().Random< size_t>( 0, m_PairedAtomIndices.GetSize() - 1));
            second_atom = FRAGMENT.GetAtomVector()( second_atom_index);

            // bond between the two chosen atoms
            for
            (
                auto bond_itr( picked_atom->GetBonds().Begin()), bond_itr_end( picked_atom->GetBonds().End());
                bond_itr != bond_itr_end;
                ++bond_itr
            )
            {
              if( FRAGMENT.GetAtomVector().GetAtomIndex( bond_itr->GetTargetAtom()) == second_atom_index)
              {
                bond_to_remove = sdf::BondInfo( picked_atom_index, second_atom_index, bond_itr->GetBondType());
              }
            }
          }

          // find a bond attached to our selected atom to break
          else
          {
            storage::Vector< size_t> bonds;
            for( size_t i( 0), sz( picked_atom->GetBonds().GetSize()); i < sz; ++i)
            {
              bonds.PushBack( i);
            }
            bonds.Shuffle();
            auto picked_atom_bonds( picked_atom->GetBonds());
            for
            (
                auto bond_itr( bonds.Begin()), bond_itr_end( bonds.End());
                bond_itr != bond_itr_end;
                ++bond_itr
            )
            {
              // bonded atom must not be hydrogen atom, nor can both atoms be in the same ring (I do not want to unfold rings with this)
              if( picked_atom_bonds( *bond_itr).GetTargetAtom().GetElementType() == GetElementTypes().e_Hydrogen)
              {
                continue;
              }
              // check if both atoms in a ring
              else if
              (
                  // both atoms can be in a ring only if they are in separate rings
                  picked_atom_bonds( *bond_itr).GetTargetAtom().CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1)) &&
                  picked_atom->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1))
              )
              {
                // if so, make sure not in same ring
                if
                (
                    !IfAtomsInSameRing
                    (
                      picked_atom_index,
                      FRAGMENT.GetAtomVector().GetAtomIndex( picked_atom_bonds( *bond_itr).GetTargetAtom()),
                      FRAGMENT
                    )
                )
                {
                  // assign second atom
                  second_atom_index = FRAGMENT.GetAtomVector().GetAtomIndex( picked_atom_bonds( *bond_itr).GetTargetAtom());
                  second_atom = FRAGMENT.GetAtomVector()( second_atom_index);

                  // save bond to remove
                  bond_to_remove = sdf::BondInfo( picked_atom_index, second_atom_index, picked_atom_bonds( *bond_itr).GetBondType());

                  // time to exit second atom picking stage
                  break;
                }
              }
              // both atoms are not in a ring
              else
              {
                // assign second atom
                second_atom_index = FRAGMENT.GetAtomVector().GetAtomIndex( picked_atom_bonds( *bond_itr).GetTargetAtom());
                second_atom = FRAGMENT.GetAtomVector()( second_atom_index);

                // save bond to remove
                bond_to_remove = sdf::BondInfo( picked_atom_index, second_atom_index, picked_atom_bonds( *bond_itr).GetBondType());

                // time to exit second atom picking stage
                break;
              }
            }
          }
          // done trying to find atom pairs
          if( util::IsDefined( second_atom_index) && util::IsDefined( picked_atom_index))
          {
            break;
          }
        }

        // if we never managed to get an atom pair then bail
        if( !util::IsDefined( second_atom_index) || !util::IsDefined( picked_atom_index))
        {
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }

        // break the bond between the two picked atoms
        FragmentMutateRemoveBond remove_bond;
        AtomVector< AtomComplete> atom_vector( remove_bond.RemoveBond( FRAGMENT, bond_to_remove));
        FragmentComplete complex_fragment( atom_vector, "");

        // isolate the two new fragments
        FragmentSplitIsolate isolater( m_FragmentMinSize);
        FragmentEnsemble isolated_frags( isolater( complex_fragment));

        // we need at least one valid fragment
        if( !isolated_frags.GetSize())
        {
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }

        if
        (
            isolated_frags.GetSize() > size_t( 1) || // must have two fragments
            m_AllowFragmentDuplication // alternatively, must allow duplication of the one fragment I do have
        )
        {
          // fun part - glue the fragments back together with a linker
          // find the open valence atom indices
          size_t atom_index( 0);
          size_t isolated_frag_one_atom_index( util::GetUndefinedSize_t()), isolated_frag_two_atom_index( util::GetUndefinedSize_t());
          for
          (
              auto frag_one_itr( isolated_frags.GetMolecules().FirstElement().GetAtomsIterator());
              frag_one_itr.NotAtEnd();
              ++frag_one_itr, ++atom_index
          )
          {
            if( frag_one_itr->GetNumberValenceBonds())
            {
              isolated_frag_one_atom_index = atom_index;
            }
          }
          atom_index = 0;
          for
          (
              auto frag_two_itr( isolated_frags.GetMolecules().LastElement().GetAtomsIterator());
              frag_two_itr.NotAtEnd();
              ++frag_two_itr, ++atom_index
          )
          {
            if( frag_two_itr->GetNumberValenceBonds())
            {
              isolated_frag_two_atom_index = atom_index;
            }
          }

          // choose a link method
          std::string method( this->ChooseLinkMethod( true));

          // amide link
          if( method == "AmideLink")
          {
            new_mol = AmideLink
                (
                  isolated_frags.GetMolecules().FirstElement(),
                  isolated_frags.GetMolecules().LastElement(),
                  isolated_frag_one_atom_index,
                  isolated_frag_two_atom_index,
                  size_t( 1)
                );
          }
          // ester link
          else if( method == "EsterLink")
          {
            new_mol = EsterLink
                (
                  isolated_frags.GetMolecules().FirstElement(),
                  isolated_frags.GetMolecules().LastElement(),
                  isolated_frag_one_atom_index,
                  isolated_frag_two_atom_index,
                  size_t( 1)
                );
          }
          // methoxy link
          else if( method == "MethoxyLink")
          {
            new_mol = MethoxyLink
                (
                  isolated_frags.GetMolecules().FirstElement(),
                  isolated_frags.GetMolecules().LastElement(),
                  isolated_frag_one_atom_index,
                  isolated_frag_two_atom_index,
                  size_t( 1)
                );
          }
          else if( method == "EthoxyLink")
          {
            new_mol = EthoxyLink
                (
                  isolated_frags.GetMolecules().FirstElement(),
                  isolated_frags.GetMolecules().LastElement(),
                  isolated_frag_one_atom_index,
                  isolated_frag_two_atom_index,
                  size_t( 1)
                );
          }
          // single element link
          else if( method == "SingleElementLink")
          {
            std::string element_type( this->ChooseLinkElement());
            new_mol = SingleElementLink
                (
                  isolated_frags.GetMolecules().FirstElement(),
                  isolated_frags.GetMolecules().LastElement(),
                  isolated_frag_one_atom_index,
                  isolated_frag_two_atom_index,
                  element_type
                );
          }
          // alkyl link
          else if( method == "AlkylLink")
          {
            // determine alkyl linker size
            float repeat_prob( random::GetGlobalRandom().Random< float>( 0.0, 1.0));
            size_t n_repeats;
            repeat_prob < m_ThreeC ?
                n_repeats = size_t( 3) :
                n_repeats = size_t( 2);

            // perform alkyl link
            new_mol = AlkylLink
                (
                  isolated_frags.GetMolecules().FirstElement(),
                  isolated_frags.GetMolecules().LastElement(),
                  isolated_frag_one_atom_index,
                  isolated_frag_two_atom_index,
                  n_repeats
                );
          }
          // ring link
          else if( method == "RingLink")
          {
            new_mol = RingLink
                (
                  isolated_frags.GetMolecules().FirstElement(),
                  isolated_frags.GetMolecules().LastElement(),
                  isolated_frag_one_atom_index,
                  isolated_frag_two_atom_index
                );
          }
          // return empty
          else
          {
            return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
          }
        }
      }
      // extend from ring
      else
      {
        BCL_MessageStd( "ExtendWithLinker - Extend linker from fragment");

        for( ; tries < m_NumberMaxAttempts; ++tries)
        {
          // pick an atom
          if( m_MutableAtomIndices.GetSize() || m_MutableElements.GetSize() || m_MutableFragments.GetSize())
          {
            picked_atom = this->PickAtom( FRAGMENT, false);
          }
          else
          {
            picked_atom = this->PickAtom( FRAGMENT, true);
          }

          // if a hydrogen atom is selected then just grab its bonded atom
          if( picked_atom->GetElementType() == GetElementTypes().e_Hydrogen)
          {
            // must have bonded atom
            if( !picked_atom->GetBonds().GetSize())
            {
              continue;
            }
            picked_atom = util::SiPtr< const AtomConformationalInterface>( picked_atom->GetBonds().Begin()->GetTargetAtom());
          }

          picked_atom_index = FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom);
          if( util::IsDefined( picked_atom_index))
          {
            break;
          }
        }

        // need to pick a legitimate atom
        if( !util::IsDefined( picked_atom_index))
        {
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }

        // link to ring
        std::string element_type( this->ChooseLinkElement());
        new_mol = RingExtension
            (
              FRAGMENT,
              picked_atom_index,
              size_t
              (
                random::GetGlobalRandom().Random< float>( 0.0, 1.0) < m_ThreeC ?
                    size_t( 3) :
                    size_t( 2)
              ), // we're basically just rolling with linker size 2 right now
              element_type
            );
      }

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

      // Remove hydrogen atoms before clean to allow proper bondtype selection
      AtomVector< AtomComplete> not_empty( new_mol.GetAtomVector());
      HydrogensHandler::Remove( not_empty);

      // Check for valid atom types
      util::ShPtr< FragmentComplete> new_mol_ptr
      (
        m_ScaffoldFragment.GetSize()
        ? cleaner.Clean( not_empty, m_ScaffoldFragment, m_DrugLikenessType)
            : cleaner.Clean( not_empty, FRAGMENT, m_DrugLikenessType)
      );

      if( !new_mol_ptr.IsDefined() || new_mol_ptr->HasNonGasteigerAtomTypes())
      {
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
      }

      // return the new molecule
      return math::MutateResult< FragmentComplete>( new_mol_ptr, *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set the absolute probability of extending the fragment internally
    void FragmentMutateExtendWithLinker::SetExtendWithinProb( const float EXTEND_WITHIN_PROB)
    {
      m_ExtendWithinProb = EXTEND_WITHIN_PROB;
    }

    //! @brief set the size of the smallest fragment to be is.olated
    void FragmentMutateExtendWithLinker::SetMinFragmentSize( const size_t MIN_FRAGMENT_SIZE)
    {
      m_FragmentMinSize = MIN_FRAGMENT_SIZE;
    }

    //! @brief set whether to allow duplicates if only one fragment meets min size requirement
    void FragmentMutateExtendWithLinker::SetAllowDuplicates( const bool ALLOW_DUPLICATES)
    {
      m_AllowFragmentDuplication = ALLOW_DUPLICATES;
    }

    //! @brief set the relative probability of extending with an amide linker
    void FragmentMutateExtendWithLinker::SetAmideLinkProb( const float AMIDE_LINK_PROB)
    {
      m_AmideLinkProb = AMIDE_LINK_PROB;
    }

    //! @brief set the probability of orienting the amide N toward fragment A with the amide linker
    void FragmentMutateExtendWithLinker::SetAmideNToAProb( const float AMIDE_N_TO_A_PROB)
    {
      m_AmideNToAProb = AMIDE_N_TO_A_PROB;
    }

    //! @brief set the relative probability of extending with an ester linker
    void FragmentMutateExtendWithLinker::SetEsterLinkProb( const float ESTER_LINK_PROB)
    {
      m_EsterLinkProb = ESTER_LINK_PROB;
    }

    //! @brief set the probability of orienting the ester O toward fragment A with the ester linker
    void FragmentMutateExtendWithLinker::SetEsterOToAProb( const float ESTER_O_TO_A_PROB)
    {
      m_EsterOToAProb = ESTER_O_TO_A_PROB;
    }

    //! @brief set the relative probability of extending with a methoxy linker
    void FragmentMutateExtendWithLinker::SetMethoxyLinkProb( const float METHOXY_LINK_PROB)
    {
      m_MethoxyLinkProb = METHOXY_LINK_PROB;
    }

    //! @brief set the probability of orienting the methoxy O toward fragment A with the methoxy linker
    void FragmentMutateExtendWithLinker::SetMethoxyOToAProb( const float METHOXY_O_TO_A_PROB)
    {
      m_MethoxyOToAProb = METHOXY_O_TO_A_PROB;
    }

    //! @brief set the relative probability of extending with an ethoxy linker
    void FragmentMutateExtendWithLinker::SetEthoxyLinkProb( const float ETHOXY_LINK_PROB)
    {
      m_EthoxyLinkProb = ETHOXY_LINK_PROB;
    }

    //! @brief set the probability of orienting the ethoxy O toward fragment A with the ethoxy linker
    void FragmentMutateExtendWithLinker::SetEthoxyOToAProb( const float ETHOXY_O_TO_A_PROB)
    {
      m_EthoxyOToAProb = ETHOXY_O_TO_A_PROB;
    }

    //! @brief set the relative probability of extending with a single element linker
    void FragmentMutateExtendWithLinker::SetSingleElementLinkProb( const float SINGLE_ELEMENT_LINK_PROB)
    {
      m_SingleElementLinkProb = SINGLE_ELEMENT_LINK_PROB;
    }

    //! @brief set the relative probability of linking with a single B
    void FragmentMutateExtendWithLinker::SetBProb( const float B_PROB)
    {
      m_B = B_PROB;
    }

    //! @brief set the relative probability of linking with a single C
    void FragmentMutateExtendWithLinker::SetCProb( const float C_PROB)
    {
      m_C = C_PROB;
    }

    //! @brief set the relative probability of linking with a single O
    void FragmentMutateExtendWithLinker::SetOProb( const float O_PROB)
    {
      m_O = O_PROB;
    }

    //! @brief set the relative probability of linking with a single N
    void FragmentMutateExtendWithLinker::SetNProb( const float N_PROB)
    {
      m_N = N_PROB;
    }

    //! @brief set the relative probability of linking with a single P
    void FragmentMutateExtendWithLinker::SetPProb( const float P_PROB)
    {
      m_P = P_PROB;
    }

    //! @brief set the relative probability of linking with a single S
    void FragmentMutateExtendWithLinker::SetSProb( const float S_PROB)
    {
      m_S = S_PROB;
    }

    //! @brief set the relative probability of linking with a single Se
    void FragmentMutateExtendWithLinker::SetSeProb( const float SE_PROB)
    {
      m_Se = SE_PROB;
    }

    //! @brief set the relative probability of extending with a direct link
    void FragmentMutateExtendWithLinker::SetDirectLinkProb( const float DIRECT_LINK_PROB)
    {
      m_DirectLinkProb = DIRECT_LINK_PROB;
    }

    //! @brief set the relative probability of extending with an alkyl linker
    void FragmentMutateExtendWithLinker::SetAlkylLinkProb( const float ALKYL_LINK_PROB)
    {
      m_AlkylLinkProb = ALKYL_LINK_PROB;
    }

    //! @brief set the relative probability that the alkyl linker will contain a triple bond
    void FragmentMutateExtendWithLinker::SetAlkyneProb( const float ALKYNE_PROB)
    {
      m_Alkyne = ALKYNE_PROB;
    }

    //! @brief set the relative probability of extending with a ring
    void FragmentMutateExtendWithLinker::SetRingLinkProb( const float RING_LINK_PROB)
    {
      m_RingLinkProb = RING_LINK_PROB;
    }

    //! @brief directly link two fragments via one of the linker strategies below
    //! @param FRAGMENT_A first molecule
    //! @param FRAGMENT_B second molecule
    //! @param LINK_INDEX_A link indices in first molecule
    //! @param LINK_INDEX_B link indices in second molecule
    //! @param RING_LIBRARY optional ring library if not provided through map
    //! @param RING_LINK_INDEX_A optional specific index by which to connect fragment A to the ring
    //! @param RING_LINK_INDEX_B optional specific index by which to connect fragment B to the ring
     FragmentComplete FragmentMutateExtendWithLinker::RingLink
     (
       const FragmentComplete &FRAGMENT_A,
       const FragmentComplete &FRAGMENT_B,
       const size_t&LINK_INDEX_A,
       const size_t&LINK_INDEX_B,
       const FragmentEnsemble &RING_LIBRARY,
       const size_t RING_LINK_INDEX_A,
       const size_t RING_LINK_INDEX_B
     ) const
     {
       // info on mutate
       BCL_MessageStd( "ExtendWithLinker - RingLink - FragmentEnsemble");

       // get a random ring from mapped filename
       s_Mutex.Lock(); // lock for file reading
       const FragmentEnsemble &rings
       (
         // prefer the ring library passed to the function directly
         RING_LIBRARY.GetSize() ?
             RING_LIBRARY :
             *( s_RingLibraries.Find( m_RingsFilename)->second)
       );
       s_Mutex.Unlock();

       // return null if there are no rings
       if( !rings.GetSize())
       {
         BCL_MessageStd( "There are no molecules in the provided ring library. Returning NULL...");
         return FragmentComplete();
       }

       // choose a ring at random from library
       FragmentEnsemble::const_iterator ring_itr( rings.Begin());
       size_t rand_pos( random::GetGlobalRandom().Random< size_t>( size_t( 0), rings.GetSize() - 1));
       std::advance( ring_itr, rand_pos);

       // run ringlink with the chosen ring
       return RingLink( FRAGMENT_A, FRAGMENT_B, LINK_INDEX_A, LINK_INDEX_B, *ring_itr, RING_LINK_INDEX_A, RING_LINK_INDEX_B);
     }

     //! @brief link two fragments via a ring
     //! @param FRAGMENT_A first molecule
     //! @param FRAGMENT_B second molecule
     //! @param LINK_INDEX_A link indices in first molecule
     //! @param LINK_INDEX_B link indices in second molecule
     //! @param RING the ring with which to bridge the two fragments
     //! @param RING_LINK_INDEX_A optional specific index by which to connect fragment A to the ring
     //! @param RING_LINK_INDEX_B optional specific index by which to connect fragment B to the ring
     //! @return the newly generated molecules
     FragmentComplete FragmentMutateExtendWithLinker::RingLink
     (
       const FragmentComplete &FRAGMENT_A,
       const FragmentComplete &FRAGMENT_B,
       const size_t &LINK_INDEX_A,
       const size_t &LINK_INDEX_B,
       const FragmentComplete &RING,
       const size_t RING_LINK_INDEX_A,
       const size_t RING_LINK_INDEX_B
     ) const
     {
       // info on mutate
       BCL_MessageStd( "ExtendWithLinker - RingLink - FragmentComplete");

       // get atominfo and bondinfo of initial fragments and ring
       storage::Vector< sdf::AtomInfo> atominfo_a( FRAGMENT_A.GetAtomVector().GetAtomInfo());
       storage::Vector< sdf::AtomInfo> atominfo_ring( RING.GetAtomVector().GetAtomInfo());

       storage::Vector< sdf::BondInfo> bondinfo_a( FRAGMENT_A.GetAtomVector().GetBondInfo());
       storage::Vector< sdf::BondInfo> bondinfo_ring( RING.GetAtomVector().GetBondInfo());

       // add ring info to fragment_a
       size_t atominfo_a_old_sz( atominfo_a.GetSize()), atominfo_ring_sz( atominfo_ring.GetSize());
       size_t bondinfo_ring_sz( bondinfo_ring.GetSize());
       for( size_t atominfo_ring_index( 0); atominfo_ring_index < atominfo_ring_sz; ++atominfo_ring_index)
       {
         atominfo_a.PushBack( atominfo_ring( atominfo_ring_index));
       }
       for( size_t bondinfo_ring_index( 0); bondinfo_ring_index < bondinfo_ring_sz; ++bondinfo_ring_index)
       {
         bondinfo_a.PushBack
         (
           sdf::BondInfo
           (
             bondinfo_ring( bondinfo_ring_index).GetAtomIndexLow() + atominfo_a_old_sz,
             bondinfo_ring( bondinfo_ring_index).GetAtomIndexHigh() + atominfo_a_old_sz,
             bondinfo_ring( bondinfo_ring_index).GetConfigurationalBondType()
           )
         );
       }

       // pick two atoms from the ring to connect to fragments A and B
       size_t rand_ring_pos_one( util::IsDefined( RING_LINK_INDEX_A) ? RING_LINK_INDEX_A : RunPickAtom( RING, true, true));
       rand_ring_pos_one += atominfo_a_old_sz;
       size_t rand_ring_pos_two( util::IsDefined( RING_LINK_INDEX_B) ? RING_LINK_INDEX_B : RunPickAtom( RING, true, true));
       rand_ring_pos_two += atominfo_a_old_sz;

       size_t ring_atom_selection_tries( 1);
       while( rand_ring_pos_one == rand_ring_pos_two && ring_atom_selection_tries < m_NumberMaxAttempts)
       {
         util::IsDefined( RING_LINK_INDEX_B) ?
             RING_LINK_INDEX_B :
             rand_ring_pos_two = RunPickAtom( RING, true, true);
         rand_ring_pos_two += atominfo_a_old_sz;
       }

       // return null if there are undefined atom indices
       if( !util::IsDefined( rand_ring_pos_one) || !util::IsDefined( rand_ring_pos_two))
       {
         BCL_MessageStd( "FragmentMutateExtendWithLinker unable to identify two valid ring atom indices with which to link new ring from ring library");
         BCL_MessageStd( "Atom index 1: " + util::Format()( rand_ring_pos_one));
         BCL_MessageStd( "Atom index 2: " + util::Format()( rand_ring_pos_two));
       }

       // add bond to link fragment A with ring
       sdf::BondInfo bond_a_ring
       (
         sdf::BondInfo
         (
           LINK_INDEX_A,
           rand_ring_pos_one,
           GetConfigurationalBondTypes().e_NonConjugatedSingleBond
         )
       );
       bondinfo_a.PushBack( bond_a_ring);

       // create new fragment for first half of molecule
       AtomVector< AtomComplete> a_and_ring( atominfo_a, bondinfo_a);
       FragmentComplete frag_a_and_ring( a_and_ring, "");

       // if there is no second molecule for attachment just return after the initial connection
       if( !FRAGMENT_B.GetSize())
       {
         BCL_MessageVrb( "Ending RingLink after addition of ring to the first fragment. No second fragment available. "
             "If this is undesired, check input. This behavior is expected under certain conditions.");
         return frag_a_and_ring;
       }

       // add fragment B to linked fragmentA_ring
       storage::Pair< bool, FragmentComplete> new_fragment
       (
         MergeFragmentComplete::MergeFragments
         (
           frag_a_and_ring,
           FRAGMENT_B,
           GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
           storage::Pair< size_t, size_t>( rand_ring_pos_two, LINK_INDEX_B)
         )
       );

       if( new_fragment.First())
       {
         return new_fragment.Second();
       }
       return FragmentComplete();

     }

     //! @brief directly link two fragments
     //! @param FRAGMENT_A first molecule
     //! @param FRAGMENT_B second molecule
     //! @param LINK_INDEX_A link indices in first molecule
     //! @param LINK_INDEX_B link indices in second molecule
     //! @return the newly generated molecules
     FragmentComplete FragmentMutateExtendWithLinker::DirectLink
     (
       const FragmentComplete &FRAGMENT_A,
       const FragmentComplete &FRAGMENT_B,
       const size_t&LINK_INDEX_A,
       const size_t&LINK_INDEX_B
     ) const
     {
       BCL_MessageStd( "ExtendWithLinker - DirectLink");

       // you need two fragments for this to be useful in any way
       BCL_Assert
       (
         FRAGMENT_A.GetSize() && FRAGMENT_B.GetSize(),
         "Attempting to directly connect two fragments when at least one fragment contains no atoms! Check input."
       );

       // cannot connect hydrogen atoms to anything
       if( FRAGMENT_A.GetAtomVector()( LINK_INDEX_A).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
       {
         BCL_MessageStd( "Whoops! Specified atom index in A is a hydrogen atom! Ignoring index...");
         return FragmentComplete();
       }
       if( FRAGMENT_B.GetAtomVector()( LINK_INDEX_B).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
       {
         BCL_MessageStd( "Whoops! Specified atom index in B is a hydrogen atom! Ignoring index...");
         return FragmentComplete();
       }
       storage::Triplet< FragmentComplete, size_t, size_t> pair_a( OpenValence( FRAGMENT_A, LINK_INDEX_A, m_OVShuffleH, m_OVReverse));
       storage::Triplet< FragmentComplete, size_t, size_t> pair_b( OpenValence( FRAGMENT_B, LINK_INDEX_B, m_OVShuffleH, m_OVReverse));

       // link fragments
       storage::Pair< bool, FragmentComplete> new_fragment
       (
         MergeFragmentComplete::MergeFragments
         (
           pair_a.First(),
           pair_b.First(),
           GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
           storage::Pair< size_t, size_t>( pair_a.Second(), pair_b.Second())
         )
       );

       // add new molecule to the ensemble
       BCL_MessageCrt( "Done directly linking fragments!");
       if( new_fragment.First())
       {
         return new_fragment.Second();
       }
       return FragmentComplete();
     }

     //! @brief link two fragments via a single element
     //! @param FRAGMENT_A first molecule
     //! @param FRAGMENT_B second molecule
     //! @param LINK_INDEX_A link indices in first molecule
     //! @param LINK_INDEX_B link indices in second molecule
     //! @param ELEMENT_TYPE element type serving as a link
     //! @return the newly generated molecules
     FragmentComplete FragmentMutateExtendWithLinker::SingleElementLink
     (
       const FragmentComplete &FRAGMENT_A,
       const FragmentComplete &FRAGMENT_B,
       const size_t&LINK_INDEX_A,
       const size_t&LINK_INDEX_B,
       const std::string &ELEMENT_TYPE
     ) const
     {
       BCL_MessageStd( "ExtendWithLinker - SingleElementLink");

       // you need at least the starting fragment
       BCL_Assert
       (
         FRAGMENT_A.GetSize(),
         "Input fragment contains no atoms!"
       );

       // open a valence on fragment A
       storage::Triplet< FragmentComplete, size_t, size_t> pair_a( OpenValence( FRAGMENT_A, LINK_INDEX_A, m_OVShuffleH, m_OVReverse));
       AtomVector< AtomComplete> atom_v_a( pair_a.First().GetAtomVector());

       // cannot connect hydrogen atoms to anything
       if( atom_v_a( pair_a.Second()).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
       {
         return FragmentComplete();
       }

       // get link atom
       // boron
       sdf::AtomInfo link_atom;

       if( ELEMENT_TYPE == "B")
       {
         link_atom = sdf::AtomInfo( GetAtomTypes().B_TeTeTe, e_NonChiral);
       }
       // carbon
       else if( ELEMENT_TYPE == "C")
       {
         link_atom = sdf::AtomInfo( GetAtomTypes().C_TeTeTeTe, e_NonChiral);
       }
       // oxygen
       else if( ELEMENT_TYPE == "O")
       {
         link_atom = sdf::AtomInfo( GetAtomTypes().O_Te2Te2TeTe, e_NonChiral);
       }
       // nitrogen
       else if( ELEMENT_TYPE == "N")
       {
         link_atom = sdf::AtomInfo( GetAtomTypes().N_Te2TeTeTe, e_NonChiral);
       }
       // phosphorous
       else if( ELEMENT_TYPE == "P")
       {
         link_atom = sdf::AtomInfo( GetAtomTypes().P_Te2TeTeTe, e_NonChiral);
       }
       // sulfur
       else if( ELEMENT_TYPE == "S")
       {
         link_atom = sdf::AtomInfo( GetAtomTypes().S_Te2Te2TeTe, e_NonChiral);
       }
       // selenium
       else if( ELEMENT_TYPE == "Se")
       {
         link_atom = sdf::AtomInfo( GetAtomTypes().Se_Te2Te2TeTe, e_NonChiral);
       }
       // oops
       else
       {
         BCL_MessageStd( "Invalid element type specified for link atom");
         return FragmentComplete();
       }

       // let's make a fragment from this little atom
       // if we just append the atominfo to e.g. the molecule_A atominfo with a
       // bond of arbitrary length then when we merge the remaining fragment B
       // and clean, the irrelevant geometry of the molecule_A+link_atom pseudomolecule
       // actually influences the directionality of the bond, which is bad. we can
       // avoid that by merging them as two fragments.
       // but, i know this looks stupid
       storage::Vector< sdf::BondInfo> empty_bonds;
       storage::Vector< sdf::AtomInfo> link_atominfo( size_t( 1), link_atom);
       AtomVector< AtomComplete> link_atom_v( link_atominfo, empty_bonds);
       FragmentComplete link_atom_frag( link_atom_v, "");

       // merge link atom with molecule A
       storage::Pair< bool, FragmentComplete> atom_v_a_mod_frag_pair
       (
         MergeFragmentComplete::MergeFragments
         (
           pair_a.First(),
           link_atom_frag,
           GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
           storage::Pair< size_t, size_t>( pair_a.Second(), size_t( 0))
         )
       );

       if( atom_v_a_mod_frag_pair.First())
       {
         // why am i making this a separate object? seems stupid
         FragmentComplete atom_v_a_mod_frag( atom_v_a_mod_frag_pair.Second());

         // if the second fragment is empty then we can just return the first fragment with the new extension
         if( !FRAGMENT_B.GetSize())
         {
           return atom_v_a_mod_frag;
         }

         // open a valence on fragment b
         storage::Triplet< FragmentComplete, size_t, size_t> pair_b( OpenValence( FRAGMENT_B, LINK_INDEX_B, m_OVShuffleH, m_OVReverse));
         AtomVector< AtomComplete> atom_v_b( pair_b.First().GetAtomVector());

         // cannot connect hydrogen atoms to anything
         if( atom_v_b( pair_b.Second()).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
         {
           return FragmentComplete();
         }

         // link fragment B to the link atom
         storage::Pair< bool, FragmentComplete> new_fragment
         (
           MergeFragmentComplete::MergeFragments
           (
             atom_v_a_mod_frag,
             pair_b.First(),
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
             storage::Pair< size_t, size_t>( atom_v_a_mod_frag.GetSize() - 1, pair_b.Second())
           )
         );

         // done
         if( new_fragment.First())
         {
           return new_fragment.Second();
         }
       }
       return FragmentComplete();
     }

     //! @brief link two fragments with an alkyl chain
     //! @param FRAGMENT_A first molecule
     //! @param FRAGMENT_B second molecule
     //! @param LINK_INDEX_A link indices in first molecule
     //! @param LINK_INDEX_B link indices in second molecule
     //! @param REPEATS the number of linker repeats
     //! @return the newly generated molecules
     FragmentComplete FragmentMutateExtendWithLinker::AlkylLink
     (
       const FragmentComplete &FRAGMENT_A,
       const FragmentComplete &FRAGMENT_B,
       const size_t&LINK_INDEX_A,
       const size_t&LINK_INDEX_B,
       const size_t &REPEATS
     ) const
     {
       BCL_MessageStd( "ExtendWithLinker - AlkylLink - " + util::Format()( REPEATS) + " Repeats");

       // repeats must be at least 2 or else this is pointless and we should
       // just do a single element link
       if( REPEATS < 2)
       {
         return FragmentComplete();
       }

       // you need at least the starting fragment
       BCL_Assert
       (
         FRAGMENT_A.GetSize(),
         "Input fragment contains no atoms!"
       );

       // cannot connect hydrogen atoms to anything
       if( FRAGMENT_A.GetAtomVector()( LINK_INDEX_A).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
       {
         return FragmentComplete();
       }

       // get atominfo and bondinfo of initial fragments and ring
       storage::Vector< sdf::AtomInfo> atominfo_a( FRAGMENT_A.GetAtomVector().GetAtomInfo());
       storage::Vector< sdf::BondInfo> bondinfo_a( FRAGMENT_A.GetAtomVector().GetBondInfo());

       // Get base link atom type
       storage::Vector< sdf::AtomInfo> alkyl_chain_atominfo;
       storage::Vector< sdf::BondInfo> alkyl_chain_bondinfo;
       sdf::AtomInfo link_atom( GetAtomTypes().C_TrTrTrPi, e_NonChiral);

       // make chain
       for( size_t repeat_index( 0); repeat_index < REPEATS; ++repeat_index)
       {
         alkyl_chain_atominfo.PushBack( link_atom);
       }
       for( size_t linker_index( 1); linker_index < REPEATS; ++linker_index)
       {
         //         if( REPEATS >= size_t( 2) && frag_a_contains_ring && frag_b_contains_ring && random::GetGlobalRandom().Random< float>( 0.0, 1.0) < 0.50)
         if( random::GetGlobalRandom().Random< float>( 0.0, 1.0) < m_Alkyne)
         {
           BCL_MessageStd( "Extend with triple bond!");
           alkyl_chain_bondinfo.PushBack
           (
             sdf::BondInfo
             (
               linker_index - 1,
               linker_index,
               GetConfigurationalBondTypes().e_ConjugatedTripleBond
             )
           );
         }
         // standard alkyl linker
         else
         {
           BCL_MessageStd( "Extend with single bond!");
           alkyl_chain_bondinfo.PushBack
           (
             sdf::BondInfo
             (
               linker_index - 1,
               linker_index,
               GetConfigurationalBondTypes().e_NonConjugatedSingleBond
             )
           );
         }
       }

       // add the alkyl chain to the first fragment
       size_t atominfo_a_old_sz( atominfo_a.GetSize()), atominfo_chain_sz( alkyl_chain_atominfo.GetSize());
       size_t bondinfo_chain_sz( alkyl_chain_bondinfo.GetSize());
       for( size_t atominfo_chain_index( 0); atominfo_chain_index < atominfo_chain_sz; ++atominfo_chain_index)
       {
         atominfo_a.PushBack( alkyl_chain_atominfo( atominfo_chain_index));
       }
       for( size_t bondinfo_ring_index( 0); bondinfo_ring_index < bondinfo_chain_sz; ++bondinfo_ring_index)
       {
         bondinfo_a.PushBack
         (
           sdf::BondInfo
           (
             alkyl_chain_bondinfo( bondinfo_ring_index).GetAtomIndexLow() + atominfo_a_old_sz,
             alkyl_chain_bondinfo( bondinfo_ring_index).GetAtomIndexHigh() + atominfo_a_old_sz,
             alkyl_chain_bondinfo( bondinfo_ring_index).GetConfigurationalBondType()
           )
         );
       }

       // add bond to link fragment A with ring
       sdf::BondInfo bond_a_chain
       (
         sdf::BondInfo
         (
           LINK_INDEX_A,
           atominfo_a_old_sz,
           GetConfigurationalBondTypes().e_NonConjugatedSingleBond
         )
       );
       bondinfo_a.PushBack( bond_a_chain);

       // create new fragment for first half of molecule
       AtomVector< AtomComplete> a_and_ring( atominfo_a, bondinfo_a);
       FragmentComplete frag_a_and_ring( a_and_ring, "");

       // if the second fragment is empty then we can just return the first fragment with the new extension
       if( !FRAGMENT_B.GetSize())
       {
         return frag_a_and_ring;
       }

       // cannot connect hydrogen atoms to anything
       if( FRAGMENT_B.GetAtomVector()( LINK_INDEX_B).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
       {
         return FragmentComplete();
       }

       // add fragment B to linked fragmentA_ring
       storage::Pair< bool, FragmentComplete> new_fragment
       (
         MergeFragmentComplete::MergeFragments
         (
           frag_a_and_ring,
           FRAGMENT_B,
           GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
           storage::Pair< size_t, size_t>( frag_a_and_ring.GetSize() - 1, LINK_INDEX_B)
         )
       );

       if( new_fragment.First())
       {
         return new_fragment.Second();
       }
       return FragmentComplete();
     }

     //! @brief link two fragments with a methoxy repeat
     //! @param FRAGMENT_A first molecule
     //! @param FRAGMENT_B second molecule
     //! @param LINK_INDEX_A link indices in first molecule
     //! @param LINK_INDEX_B link indices in second molecule
     //! @param REPEATS the number of linker repeats
     //! @return the newly generated molecules
     FragmentComplete FragmentMutateExtendWithLinker::MethoxyLink
     (
       const FragmentComplete &FRAGMENT_A,
       const FragmentComplete &FRAGMENT_B,
       const size_t&LINK_INDEX_A,
       const size_t&LINK_INDEX_B,
       const size_t &REPEATS
     ) const
     {
       BCL_MessageStd( "ExtendWithLinker - MethoxyLink");

       // you need at least the starting fragment
       BCL_Assert
       (
         FRAGMENT_A.GetSize(),
         "Input fragment contains no atoms!"
       );

       // open a valence on A
       storage::Triplet< FragmentComplete, size_t, size_t> pair_a( OpenValence( FRAGMENT_A, LINK_INDEX_A, m_OVShuffleH, m_OVReverse));

       storage::Triplet< FragmentComplete, size_t, size_t> pair_b;
       if( FRAGMENT_B.GetSize())
       {
         // open a valence on B
         pair_b = OpenValence( FRAGMENT_B, LINK_INDEX_B, m_OVShuffleH, m_OVReverse);
       }

       // get atom vector to modify
       AtomVector< AtomComplete> atom_v_a( pair_a.First().GetAtomVector());

       // cannot connect hydrogen atoms to anything
       if( atom_v_a( pair_a.Second()).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
       {
         return FragmentComplete();
       }

       // this will be our final molecule
       storage::Pair< bool, FragmentComplete> complete_fragment;

       // get atominfo and bondinfo vectors for fragment A
       storage::Vector< sdf::AtomInfo> atom_info_a( atom_v_a.GetAtomInfo());
       storage::Vector< sdf::BondInfo> bond_info_a( atom_v_a.GetBondInfo());

       // create the methoxy atoms
       sdf::AtomInfo c( GetAtomTypes().C_TeTeTeTe, e_NonChiral);
       sdf::AtomInfo o( GetAtomTypes().O_Te2Te2TeTe, e_NonChiral);

       // set reasonable coordinates on amide atoms
       c.SetCoordinates( linal::Vector3D( -0.037, 0.053, 0.018));
       o.SetCoordinates( linal::Vector3D( 0.875, 1.124, -0.141));

       // determine directionality
       float rand_direction( random::GetGlobalRandom().Random< float>( 0.0, 1.0));

       // add amide N to fragment A, add so that new atoms are in the back of the atom vector
       if( rand_direction < m_MethoxyOToAProb)
       {
         // fill the atom info vector
         atom_info_a.PushBack( o);
         atom_info_a.PushBack( c);

         // connect the methoxy O and C
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             size_t( atom_v_a.GetSize()), // methoxy O
             size_t( atom_v_a.GetSize() + 1), // methoxy C
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond
           )
         );

         // connect the methoxy O with the fragment_a link atom
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             pair_a.Second(), // fragment_a link atom
             size_t( atom_v_a.GetSize()), // methoxy O
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond
           )
         );
         AtomVector< AtomComplete> fragment_a_methoxy_v( atom_info_a, bond_info_a);
         FragmentComplete fragment_a_methoxy( fragment_a_methoxy_v, "");

         // return the current fragment with the methoxy extension if there are no atoms in B
         if( !FRAGMENT_B.GetSize())
         {
           return fragment_a_methoxy;
         }

         // cannot connect hydrogen atoms to anything
         if( FRAGMENT_B.GetAtomVector()( LINK_INDEX_B).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
         {
           return FragmentComplete();
         }

         // link the amide C with fragment B
         complete_fragment = storage::Pair< bool, FragmentComplete>
         (
           MergeFragmentComplete::MergeFragments
           (
             fragment_a_methoxy,
             pair_b.First(),
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
             storage::Pair< size_t, size_t>( fragment_a_methoxy.GetSize() - 1, pair_b.Second())
           )
         );
       }
       // add methoxy C to fragment A, retain directionality such that new atoms are added to back of vector
       else
       {
         // fill the atom info vector
         atom_info_a.PushBack( c);
         atom_info_a.PushBack( o);

         // connect the methoxy C and O
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             size_t( atom_v_a.GetSize()), // methoxy C
             size_t( atom_v_a.GetSize() + 1), // methoxy O
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond
           )
         );

         // connect the methoxy C with the fragment_a link atom
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             pair_a.Second(), // fragment_a link atom
             size_t( atom_v_a.GetSize()), // methoxy C
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond
           )
         );
         AtomVector< AtomComplete> fragment_a_methoxy_v( atom_info_a, bond_info_a);
         FragmentComplete fragment_a_methoxy( fragment_a_methoxy_v, "");

         // return the current fragment with the methoxy extension if there are no atoms in B
         if( !FRAGMENT_B.GetSize())
         {
           return fragment_a_methoxy;
         }

         // cannot connect hydrogen atoms to anything
         if( FRAGMENT_B.GetAtomVector()( LINK_INDEX_B).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
         {
           return FragmentComplete();
         }

         // link the amide N with fragment B
         complete_fragment = storage::Pair< bool, FragmentComplete>
         (
           MergeFragmentComplete::MergeFragments
           (
             fragment_a_methoxy,
             pair_b.First(),
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
             storage::Pair< size_t, size_t>( fragment_a_methoxy.GetSize() - 1, pair_b.Second())
           )
         );
       }
       // return new combined fragment if successful
       if( complete_fragment.First())
       {
         return complete_fragment.Second();
       }
       return FragmentComplete();
     }

     //! @brief link two fragments with an ethoxy repeat
     //! @param FRAGMENT_A first molecule
     //! @param FRAGMENT_B second molecule
     //! @param LINK_INDEX_A link indices in first molecule
     //! @param LINK_INDEX_B link indices in second molecule
     //! @param REPEATS the number of linker repeats
     //! @return the newly generated molecules
     FragmentComplete FragmentMutateExtendWithLinker::EthoxyLink
     (
       const FragmentComplete &FRAGMENT_A,
       const FragmentComplete &FRAGMENT_B,
       const size_t&LINK_INDEX_A,
       const size_t&LINK_INDEX_B,
       const size_t &REPEATS
     ) const
     {
       BCL_MessageStd( "ExtendWithLinker - EthoxyLink");

       // you need at least the starting fragment
       BCL_Assert
       (
         FRAGMENT_A.GetSize(),
         "Input fragment contains no atoms!"
       );

       // open a valence on A
       storage::Triplet< FragmentComplete, size_t, size_t> pair_a( OpenValence( FRAGMENT_A, LINK_INDEX_A, m_OVShuffleH, m_OVReverse));

       storage::Triplet< FragmentComplete, size_t, size_t> pair_b;
       if( FRAGMENT_B.GetSize())
       {
         // open a valence on B
         pair_b = OpenValence( FRAGMENT_B, LINK_INDEX_B, m_OVShuffleH, m_OVReverse);
       }

       // get atom vector to modify
       AtomVector< AtomComplete> atom_v_a( pair_a.First().GetAtomVector());

       // cannot connect hydrogen atoms to anything
       if( atom_v_a( pair_a.Second()).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
       {
         return FragmentComplete();
       }

       // this will be our final molecule
       storage::Pair< bool, FragmentComplete> complete_fragment;

       // get atominfo and bondinfo vectors for each fragment
       storage::Vector< sdf::AtomInfo> atom_info_a( atom_v_a.GetAtomInfo());
       storage::Vector< sdf::BondInfo> bond_info_a( atom_v_a.GetBondInfo());

       // create the ethoxy atoms
       sdf::AtomInfo c1( GetAtomTypes().C_TeTeTeTe, e_NonChiral);
       sdf::AtomInfo c2( GetAtomTypes().C_TeTeTeTe, e_NonChiral);
       sdf::AtomInfo o( GetAtomTypes().O_Te2Te2TeTe, e_NonChiral);

       // set reasonable coordinates on amide atoms
       c1.SetCoordinates(linal::Vector3D( 0.564, -1.263, -0.549));
       c2.SetCoordinates(linal::Vector3D( 0.058, 0.005, 0.106));
       o.SetCoordinates( linal::Vector3D( 1.047, 1.016, -0.018));

       // determine directionality
       float rand_direction( random::GetGlobalRandom().Random< float>( 0.0, 1.0));

       // add ethoxy O to fragment A, add so that new atoms are in the back of the atom vector
       if( rand_direction < m_EthoxyOToAProb)
       {
         // fill the atom info vector
         atom_info_a.PushBack( o);
         atom_info_a.PushBack( c2);
         atom_info_a.PushBack( c1);

         // connect the ethoxy c1 and c2 atoms
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             size_t( atom_v_a.GetSize() + 2), // methoxy c1
             size_t( atom_v_a.GetSize() + 1), // methoxy c2
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond
           )
         );
         // connect the methoxy c2 and o
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             size_t( atom_v_a.GetSize() + 1), // methoxy c2
             size_t( atom_v_a.GetSize()), // methoxy o
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond
           )
         );

         // connect the methoxy o with the fragment_a link atom
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             pair_a.Second(), // fragment_a link atom
             size_t( atom_v_a.GetSize()), // methoxy o
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond
           )
         );
         AtomVector< AtomComplete> fragment_a_ethoxy_v( atom_info_a, bond_info_a);
         FragmentComplete fragment_a_ethoxy( fragment_a_ethoxy_v, "");

         // return the current fragment with the ethoxy extension if there are no atoms in B
         if( !FRAGMENT_B.GetSize())
         {
           return fragment_a_ethoxy;
         }

         // cannot connect hydrogen atoms to anything
         if( FRAGMENT_B.GetAtomVector()( LINK_INDEX_B).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
         {
           return FragmentComplete();
         }

         // link the ethoxy C with fragment B
         complete_fragment = storage::Pair< bool, FragmentComplete>
         (
           MergeFragmentComplete::MergeFragments
           (
             fragment_a_ethoxy,
             pair_b.First(),
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
             storage::Pair< size_t, size_t>( fragment_a_ethoxy.GetSize() - 3, pair_b.Second())
           )
         );
       }
       // add ethoxy C to fragment A, retain directionality such that new atoms are added to back of vector
       else
       {
         // fill the atom info vector
         atom_info_a.PushBack( c1);
         atom_info_a.PushBack( c2);
         atom_info_a.PushBack( o);

         // connect the methoxy o and c2
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             size_t( atom_v_a.GetSize() + 2), // methoxy o
             size_t( atom_v_a.GetSize() + 1), // methoxy c2
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond
           )
         );
         // connect the methoxy carbon atoms
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             size_t( atom_v_a.GetSize() + 1), // methoxy c2
             size_t( atom_v_a.GetSize()), // methoxy c1
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond
           )
         );

         // connect the methoxy with the fragment_a link atom
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             pair_a.Second(), // fragment_a link atom
             size_t( atom_v_a.GetSize()), // methoxy c1
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond
           )
         );
         AtomVector< AtomComplete> fragment_a_ethoxy_v( atom_info_a, bond_info_a);
         FragmentComplete fragment_a_ethoxy( fragment_a_ethoxy_v, "");

         // return the current fragment with the ethoxy extension if there are no atoms in B
         if( !FRAGMENT_B.GetSize())
         {
           return fragment_a_ethoxy;
         }

         // cannot connect hydrogen atoms to anything
         if( FRAGMENT_B.GetAtomVector()( LINK_INDEX_B).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
         {
           return FragmentComplete();
         }

         // link the ethoxy O with fragment B
         complete_fragment = storage::Pair< bool, FragmentComplete>
         (
           MergeFragmentComplete::MergeFragments
           (
             fragment_a_ethoxy,
             pair_b.First(),
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
             storage::Pair< size_t, size_t>( fragment_a_ethoxy.GetSize() - 3, pair_b.Second())
           )
         );
       }

       // return new combined fragment if successful
       if( complete_fragment.First())
       {
         return complete_fragment.Second();
       }
       return FragmentComplete();
     }

     //! @brief link two fragments with an amide repeat
     //! @param FRAGMENT_A first molecule
     //! @param FRAGMENT_B second molecule
     //! @param LINK_INDEX_A link indices in first molecule
     //! @param LINK_INDEX_B link indices in second molecule
     //! @param REPEATS the number of linker repeats
     //! @return the newly generated molecules
     FragmentComplete FragmentMutateExtendWithLinker::AmideLink
     (
       const FragmentComplete &FRAGMENT_A,
       const FragmentComplete &FRAGMENT_B,
       const size_t&LINK_INDEX_A,
       const size_t&LINK_INDEX_B,
       const size_t &REPEATS
     ) const
     {
       BCL_MessageStd( "ExtendWithLinker - AmideLink");

       // you need at least the starting fragment
       BCL_Assert
       (
         FRAGMENT_A.GetSize(),
         "Input fragment contains no atoms!"
       );

       // open a valence on A
       storage::Triplet< FragmentComplete, size_t, size_t> pair_a( OpenValence( FRAGMENT_A, LINK_INDEX_A, m_OVShuffleH, m_OVReverse));

       storage::Triplet< FragmentComplete, size_t, size_t> pair_b;
       if( FRAGMENT_B.GetSize())
       {
         // open a valence on B
         pair_b = OpenValence( FRAGMENT_B, LINK_INDEX_B, m_OVShuffleH, m_OVReverse);
       }

       // get atom vector to modify
       AtomVector< AtomComplete> atom_v_a( pair_a.First().GetAtomVector());

       // cannot connect hydrogen atoms to anything
       if( atom_v_a( pair_a.Second()).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
       {
         return FragmentComplete();
       }

       // this will be our final molecule
       storage::Pair< bool, FragmentComplete> complete_fragment;

       // get atominfo vectors for each fragment
       storage::Vector< sdf::AtomInfo> atom_info_a( atom_v_a.GetAtomInfo());

       // get bondinfo vectors for each fragment
       storage::Vector< sdf::BondInfo> bond_info_a( atom_v_a.GetBondInfo());

       // create the amide atoms
       sdf::AtomInfo n( GetAtomTypes().N_Te2TeTeTe, e_NonChiral);
       sdf::AtomInfo c( GetAtomTypes().C_TrTrTrPi, e_NonChiral);
       sdf::AtomInfo o( GetAtomTypes().O_Tr2Tr2TrPi, e_NonChiral);

       // set reasonable coordinates on amide atoms
       n.SetCoordinates(linal::Vector3D( -1.1906, 1.8666, 0.0178));
       c.SetCoordinates(linal::Vector3D( -0.0144, 1.2086, 0.0087));
       o.SetCoordinates(linal::Vector3D( 0.0021, -0.0041, 0.0020));

       // determine directionality
       float rand_direction( random::GetGlobalRandom().Random< float>( 0.0, 1.0));

       // add amide N to fragment A, add so that new atoms are in the back of the atom vector
       if( rand_direction < m_AmideNToAProb)
       {
         // fill the atom info vector
         atom_info_a.PushBack( n);
         atom_info_a.PushBack( o);
         atom_info_a.PushBack( c);

         // connect the amide N and C
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             size_t( atom_v_a.GetSize()), // amide N
             size_t( atom_v_a.GetSize() + 2), // amide C
             GetConfigurationalBondTypes().e_AmideSingleBond
           )
         );
         // connect the amide C and O
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             size_t( atom_v_a.GetSize() + 1), // amide O
             size_t( atom_v_a.GetSize() + 2), // amide C
             GetConfigurationalBondTypes().e_ConjugatedDoubleBond
           )
         );

         // connect the amide N with the fragment_a link atom
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             pair_a.Second(), // fragment_a link atom
             size_t( atom_v_a.GetSize()), // amide N
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond
           )
         );
         AtomVector< AtomComplete> fragment_a_amide_v( atom_info_a, bond_info_a);
         FragmentComplete fragment_a_amide( fragment_a_amide_v, "");

         // return the current fragment with the amide extension if there are no atoms in B
         if( !FRAGMENT_B.GetSize())
         {
           return fragment_a_amide;
         }

         // cannot connect hydrogen atoms to anything
         if( FRAGMENT_B.GetAtomVector()( LINK_INDEX_B).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
         {
           return FragmentComplete();
         }

         // link the amide C with fragment B
         complete_fragment = storage::Pair< bool, FragmentComplete>
         (
           MergeFragmentComplete::MergeFragments
           (
             fragment_a_amide,
             pair_b.First(),
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
             storage::Pair< size_t, size_t>( fragment_a_amide.GetSize() - 1, pair_b.Second())
           )
         );
       }
       // add amide C to fragment A, retain directionality such that new atoms are added to back of vector
       else
       {
         // fill the atom info vector
         atom_info_a.PushBack( c);
         atom_info_a.PushBack( o);
         atom_info_a.PushBack( n);

         // connect the amide C and N
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             size_t( atom_v_a.GetSize()), // amide C
             size_t( atom_v_a.GetSize() + 2), // amide N
             GetConfigurationalBondTypes().e_AmideSingleBond
           )
         );
         // connect the amide C and O
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             size_t( atom_v_a.GetSize()), // amide C
             size_t( atom_v_a.GetSize() + 1), // amide O
             GetConfigurationalBondTypes().e_ConjugatedDoubleBond
           )
         );

         // connect the amide C with the fragment_a link atom
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             pair_a.Second(), // fragment_a link atom
             size_t( atom_v_a.GetSize()), // amide C
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond
           )
         );
         AtomVector< AtomComplete> fragment_a_amide_v( atom_info_a, bond_info_a);
         FragmentComplete fragment_a_amide( fragment_a_amide_v, "");

         // return the current fragment with the ethoxy extension if there are no atoms in B
         if( !FRAGMENT_B.GetSize())
         {
           return fragment_a_amide;
         }

         // cannot connect hydrogen atoms to anything
         if( FRAGMENT_B.GetAtomVector()( LINK_INDEX_B).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
         {
           return FragmentComplete();
         }

         // link the amide N with fragment B
         complete_fragment = storage::Pair< bool, FragmentComplete>
         (
           MergeFragmentComplete::MergeFragments
           (
             fragment_a_amide,
             pair_b.First(),
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
             storage::Pair< size_t, size_t>( fragment_a_amide.GetSize() - 1, pair_b.Second())
           )
         );
       }

       // return new combined fragment if successful
       if( complete_fragment.First())
       {
         return complete_fragment.Second();
       }
       return FragmentComplete();
     }

     //! @brief link two fragments with an ester repeat
     //! @param FRAGMENT_A first molecule
     //! @param FRAGMENT_B second molecule
     //! @param LINK_INDEX_A link indices in first molecule
     //! @param LINK_INDEX_B link indices in second molecule
     //! @param REPEATS the number of linker repeats
     //! @return the newly generated molecules
     FragmentComplete FragmentMutateExtendWithLinker::EsterLink
     (
       const FragmentComplete &FRAGMENT_A,
       const FragmentComplete &FRAGMENT_B,
       const size_t&LINK_INDEX_A,
       const size_t&LINK_INDEX_B,
       const size_t &REPEATS
     ) const
     {
       BCL_MessageStd( "ExtendWithLinker - EsterLink");

       // you need at least the starting fragment
       BCL_Assert
       (
         FRAGMENT_A.GetSize(),
         "Input fragment contains no atoms!"
       );

       // open a valence on A
       storage::Triplet< FragmentComplete, size_t, size_t> pair_a( OpenValence( FRAGMENT_A, LINK_INDEX_A, m_OVShuffleH, m_OVReverse));

       storage::Triplet< FragmentComplete, size_t, size_t> pair_b;
       if( FRAGMENT_B.GetSize())
       {
         // open a valence on B
         pair_b = OpenValence( FRAGMENT_B, LINK_INDEX_B, m_OVShuffleH, m_OVReverse);
       }

       // get atom vector to modify
       AtomVector< AtomComplete> atom_v_a( pair_a.First().GetAtomVector());

       // cannot connect hydrogen atoms to anything
       if( atom_v_a( pair_a.Second()).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
       {
         return FragmentComplete();
       }

       // this will be our final molecule
       storage::Pair< bool, FragmentComplete> complete_fragment;

       // get atominfo vectors for each fragment
       storage::Vector< sdf::AtomInfo> atom_info_a( atom_v_a.GetAtomInfo());

       // get bondinfo vectors for each fragment
       storage::Vector< sdf::BondInfo> bond_info_a( atom_v_a.GetBondInfo());

       // create the ester atoms
       sdf::AtomInfo os( GetAtomTypes().O_Te2Te2TeTe, e_NonChiral);
       sdf::AtomInfo c( GetAtomTypes().C_TrTrTrPi, e_NonChiral);
       sdf::AtomInfo od( GetAtomTypes().O_Tr2Tr2TrPi, e_NonChiral);

       // set reasonable coordinates on ester atoms
       // -1.186,   1.859,   0.018
       // -0.014,   1.204,   0.009
       //  0.002,  -0.004,   0.002
       os.SetCoordinates(linal::Vector3D( -1.186, 1.859, 0.018));
       c.SetCoordinates(linal::Vector3D( -0.014, 1.204, 0.009));
       od.SetCoordinates(linal::Vector3D( 0.002, -0.004, 0.002));

       // determine directionality
       float rand_direction( random::GetGlobalRandom().Random< float>( 0.0, 1.0));

       // add ester O to fragment A, add so that new atoms are in the back of the atom vector
       if( rand_direction < m_EsterOToAProb)
       {
         // fill the atom info vector
         atom_info_a.PushBack( os);
         atom_info_a.PushBack( od);
         atom_info_a.PushBack( c);

         // connect the ester O and C
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             size_t( atom_v_a.GetSize()), // ester Os
             size_t( atom_v_a.GetSize() + 2), // ester C
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond
           )
         );
         // connect the ester C and O
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             size_t( atom_v_a.GetSize() + 1), // ester Od
             size_t( atom_v_a.GetSize() + 2), // ester C
             GetConfigurationalBondTypes().e_ConjugatedDoubleBond
           )
         );

         // connect the ester Os with the fragment_a link atom
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             pair_a.Second(), // fragment_a link atom
             size_t( atom_v_a.GetSize()), // ester Os
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond
           )
         );
         AtomVector< AtomComplete> fragment_a_ester_v( atom_info_a, bond_info_a);
         FragmentComplete fragment_a_ester( fragment_a_ester_v, "");

         // return the current fragment with the ester extension if there are no atoms in B
         if( !FRAGMENT_B.GetSize())
         {
           return fragment_a_ester;
         }

         // cannot connect hydrogen atoms to anything
         if( FRAGMENT_B.GetAtomVector()( LINK_INDEX_B).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
         {
           return FragmentComplete();
         }

         // link the ester C with fragment B
         complete_fragment = storage::Pair< bool, FragmentComplete>
         (
           MergeFragmentComplete::MergeFragments
           (
             fragment_a_ester,
             pair_b.First(),
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
             storage::Pair< size_t, size_t>( fragment_a_ester.GetSize() - 1, pair_b.Second())
           )
         );
       }
       // add ester C to fragment A, retain directionality such that new atoms are added to back of vector
       else
       {
         // fill the atom info vector
         atom_info_a.PushBack( c);
         atom_info_a.PushBack( od);
         atom_info_a.PushBack( os);

         // connect the ester C and N
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             size_t( atom_v_a.GetSize()), // ester C
             size_t( atom_v_a.GetSize() + 2), // ester Os
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond
           )
         );
         // connect the ester C and Od
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             size_t( atom_v_a.GetSize()), // ester C
             size_t( atom_v_a.GetSize() + 1), // ester Od
             GetConfigurationalBondTypes().e_ConjugatedDoubleBond
           )
         );

         // connect the ester C with the fragment_a link atom
         bond_info_a.PushBack
         (
           sdf::BondInfo
           (
             pair_a.Second(), // fragment_a link atom
             size_t( atom_v_a.GetSize()), // ester C
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond
           )
         );
         AtomVector< AtomComplete> fragment_a_ester_v( atom_info_a, bond_info_a);
         FragmentComplete fragment_a_ester( fragment_a_ester_v, "");

         // return the current fragment with the ethoxy extension if there are no atoms in B
         if( !FRAGMENT_B.GetSize())
         {
           return fragment_a_ester;
         }

         // cannot connect hydrogen atoms to anything
         if( FRAGMENT_B.GetAtomVector()( LINK_INDEX_B).GetAtomType()->GetElementType() == GetElementTypes().e_Hydrogen)
         {
           return FragmentComplete();
         }

         // link the ester O with fragment B
         complete_fragment = storage::Pair< bool, FragmentComplete>
         (
           MergeFragmentComplete::MergeFragments
           (
             fragment_a_ester,
             pair_b.First(),
             GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
             storage::Pair< size_t, size_t>( fragment_a_ester.GetSize() - 1, pair_b.Second())
           )
         );
       }

       // return new combined fragment if successful
       if( complete_fragment.First())
       {
         return complete_fragment.Second();
       }
       return FragmentComplete();
     }

     //! @brief link two fragments with an amide repeat
     //! @param CONNECTION whether to connect the amide to molecule A via C or N or both
     //! @param REPEATS the number of linker repeats
     //! @return the newly generated molecules
     // TODO DEPRECATE this function; currently in use in apps/internal/chemistry/bcl_app_link_fragments.cpp
     FragmentEnsemble FragmentMutateExtendWithLinker::GenerateAmideLinker
     (
       const std::string &CONNECTION,
       const size_t &REPEATS
     ) const
     {
       // Begin
       FragmentEnsemble amide;

       // Read in amide linker file
       io::IFStream file;
       io::File::MustOpenIFStream
       (
         file,
         RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "combichem_linkers/amide.sdf.gz"
       );
       amide.ReadMoreFromMdl( file, sdf::e_Maintain);

       io::File::CloseClearFStream( file);

       // amide linker repeat length of 1
       if( REPEATS == size_t( 1))
       {
         return amide;
       }
       else
       {
         return amide;
       }
       return FragmentEnsemble();
     }

     //! @brief extends a new ring from a current ring with one of the above linkers
     //! @param FRAGMENT_A first molecule
     //! @param LINK_INDEX_A link indices in first molecule
     //! @param REPEATS the number of linker repeats
     //! @param ELEMENT_TYPE element type serving as a link of SingleElementLinker
     //! @return the newly generated molecules
     FragmentComplete FragmentMutateExtendWithLinker::RingExtension
     (
       const FragmentComplete &FRAGMENT_A,
       const size_t &LINK_INDEX_A,
       const size_t &REPEATS,
       const std::string &ELEMENT_TYPE
     ) const
     {
       BCL_MessageStd( "ExtendWithLinker - RingExtension");
       BCL_Assert( m_RingsFilename.size(), "No ring library filename specified! Cannot perform 'ExtendWithLinker - RingExtension'");

       // get a ring
       s_Mutex.Lock();
       const FragmentEnsemble &rings( *( s_RingLibraries.Find( m_RingsFilename)->second));
       s_Mutex.Unlock();
       BCL_Assert( rings.GetSize(), "No molecules in the provided ring library file! Cannot perform 'ExtendWithLinker - RingExtension'");

       // advance to random ring
       size_t rand_pos( random::GetGlobalRandom().Random< size_t>( size_t( 0), rings.GetSize() - 1));
       auto ring_itr( rings.Begin());
       std::advance( ring_itr, rand_pos);

       // get a random ring atom
//       size_t rand_ring_atom( random::GetGlobalRandom().Random< size_t>( 0, ring_itr->GetSize() - 1));
       size_t rand_ring_atom( RunPickAtom( *ring_itr, true, true));

       // link the ring to the starting fragment ring
       FragmentComplete new_mol;
       std::string method( this->ChooseLinkMethod( false));
       if( method == "AmideLink")
       {
         new_mol = AmideLink
             (
               FRAGMENT_A,
               *ring_itr,
               LINK_INDEX_A,
               rand_ring_atom,
               size_t( 1)
             );
       }
       else if( method == "EsterLink")
       {
         new_mol = EsterLink
             (
               FRAGMENT_A,
               *ring_itr,
               LINK_INDEX_A,
               rand_ring_atom,
               size_t( 1)
             );
       }
       else if( method == "MethoxyLink")
       {
         new_mol = MethoxyLink
             (
               FRAGMENT_A,
               *ring_itr,
               LINK_INDEX_A,
               rand_ring_atom,
               size_t( 1)
             );
       }
       else if( method == "EthoxyLink")
       {
         new_mol = EthoxyLink
             (
               FRAGMENT_A,
               *ring_itr,
               LINK_INDEX_A,
               rand_ring_atom,
               size_t( 1)
             );
       }
       else if( method == "DirectLink")
       {
         new_mol = DirectLink
             (
               FRAGMENT_A,
               *ring_itr,
               LINK_INDEX_A,
               rand_ring_atom
             );
       }
       else if( method == "SingleElementLink")
       {
         new_mol = SingleElementLink
             (
               FRAGMENT_A,
               *ring_itr,
               LINK_INDEX_A,
               rand_ring_atom,
               ELEMENT_TYPE
             );
       }
       else if( method == "AlkylLink")
       {
         new_mol = AlkylLink
             (
               FRAGMENT_A,
               *ring_itr,
               LINK_INDEX_A,
               rand_ring_atom,
               REPEATS
             );
       }

       // end
       if( new_mol.GetSize())
       {
         return new_mol;
       }

       return FragmentComplete();
     }

     //! @brief selects an element to use with the SingleElementLink function
     //! @return string corresponding to chosen element
     std::string FragmentMutateExtendWithLinker::ChooseLinkMethod( const bool &EXTEND_WITHIN) const
     {
       // determine what type of linkage to perform
       float sum_probs( 0.0);

       // different probabilities depending on whether or not we are extending the fragment internally
       EXTEND_WITHIN ?
       sum_probs = m_AmideLinkProb + m_EsterLinkProb + m_MethoxyLinkProb + m_EthoxyLinkProb + m_SingleElementLinkProb + m_AlkylLinkProb + m_RingLinkProb :
       sum_probs = m_AmideLinkProb + m_EsterLinkProb + m_MethoxyLinkProb + m_EthoxyLinkProb + m_SingleElementLinkProb + m_AlkylLinkProb + m_DirectLinkProb;

       // if everything is zero then return an empty string
       if( !sum_probs)
       {
         return std::string();
       }

       // options
       storage::Vector< std::string> options;
       options.PushBack( "AmideLink");
       options.PushBack( "EsterLink");
       options.PushBack( "MethoxyLink");
       options.PushBack( "EthoxyLink");
       options.PushBack( "SingleElementLink");
       options.PushBack( "AlkylLink");
       EXTEND_WITHIN ?
         options.PushBack( "RingLink") :
         options.PushBack( "DirectLink");

       // probs
       storage::Vector< float> probs;
       probs.PushBack( m_AmideLinkProb);
       probs.PushBack( m_EsterLinkProb);
       probs.PushBack( m_MethoxyLinkProb);
       probs.PushBack( m_EthoxyLinkProb);
       probs.PushBack( m_SingleElementLinkProb);
       probs.PushBack( m_AlkylLinkProb);
       EXTEND_WITHIN ?
         probs.PushBack( m_RingLinkProb) :
         probs.PushBack( m_DirectLinkProb);

       // obtain a link method
       float rand( random::GetGlobalRandom().Random< float>( 0, sum_probs));
       float cumulative_weighted_sum( 0.0);
       size_t i( 0);
       for( ; i < options.GetSize(); ++i)
       {
         cumulative_weighted_sum += probs( i);
         if( rand < cumulative_weighted_sum)
         {
           break;
         }
       }

       // return the link method string
       return options( i);
     }

     //! @brief selects an element to use with the SingleElementLink function
     //! @return string corresponding to chosen element
     std::string FragmentMutateExtendWithLinker::ChooseLinkElement() const
     {
       // determine what type of SingleElementLink element to use
       float sum_probs( m_B + m_C + m_O + m_N + m_P + m_S + m_Se);

       // if everything is zero then return an empty string
       if( !sum_probs)
       {
         return std::string();
       }

       // options
       storage::Vector< std::string> options;
       options.PushBack( "B");
       options.PushBack( "C");
       options.PushBack( "O");
       options.PushBack( "N");
       options.PushBack( "P");
       options.PushBack( "S");
       options.PushBack( "Se");

       // probs
       storage::Vector< float> probs;
       probs.PushBack( m_B);
       probs.PushBack( m_C);
       probs.PushBack( m_O);
       probs.PushBack( m_N);
       probs.PushBack( m_P);
       probs.PushBack( m_S);
       probs.PushBack( m_Se);

       // obtain a link method
       float rand( random::GetGlobalRandom().Random< float>( 0, sum_probs));
       float cumulative_weighted_sum( 0.0);
       size_t i( 0);
       for( ; i < options.GetSize(); ++i)
       {
         cumulative_weighted_sum += probs( i);
         if( rand < cumulative_weighted_sum)
         {
           break;
         }
       }
       // return the link method string
       return options( i);
     }

     //! @brief checks if two atoms in rings are in the same ring or different rings
     //! @param ATOM_INDEX_A first atom
     //! @param ATOM_INDEX_B second atom
     //! @param PARENT_MOL the fragment containing both atoms
     //! @return true if the atoms are in the same ring, false otherwise
     bool FragmentMutateExtendWithLinker::IfAtomsInSameRing
     (
       const size_t &ATOM_INDEX_A,
       const size_t &ATOM_INDEX_B,
       const FragmentComplete &PARENT_MOL
     ) const
     {
       // create a ring splitter
       FragmentSplitRings ringsplit( size_t( 2), true);

       // create a graph of our molecule of interest
       ConformationGraphConverter::t_AtomGraph mol_graph( ConformationGraphConverter::CreateGraphWithAtoms( PARENT_MOL));

       // split out the rings
       auto ring_systems( ringsplit.GetComponentVertices( PARENT_MOL, mol_graph));

       // find out which ring contains atom a
       size_t ring_index( 0);
       bool this_ring_has_a( false);
       for
       (
           auto ring_itr( ring_systems.Begin()), ring_itr_end( ring_systems.End());
           ring_itr != ring_itr_end;
           ++ring_itr, ++ring_index
       )
       {
         for
         (
             auto atom_itr( ring_itr->Begin()), atom_itr_end( ring_itr->End());
             atom_itr != atom_itr_end;
             ++atom_itr
         )
         {
           if( *atom_itr == ATOM_INDEX_A)
           {
             this_ring_has_a = true;
           }
         }
         if( this_ring_has_a)
         {
           break;
         }
       }

       // if atom a is not found in any rings then just be done
       if( !this_ring_has_a)
       {
         return false;
       }

       // now search only the thing that has atom a for atom b
       auto ring_itr( ring_systems.Begin());
       std::advance( ring_itr, ring_index);
       for
       (
           auto atom_itr( ring_itr->Begin()), atom_itr_end( ring_itr->End());
           atom_itr != atom_itr_end;
           ++atom_itr
       )
       {
         if( *atom_itr == ATOM_INDEX_B)
         {
           return true;
         }
       }
       return false;
     }

  //////////////////////
  // helper functions //
  //////////////////////

    io::Serializer FragmentMutateExtendWithLinker::GetSerializer() const
    {
      io::Serializer parameters( FragmentMutateInterface::GetSerializer());
      parameters.SetClassDescription
      (
        "Extends a molecule either by breaking a bond and inserting a motif before reconnecting, "
        "or by linking to a new ring system"
      );

      parameters.AddInitializer
      (
        "ring_library",
        "path to the ring library",
        io::Serialization::GetAgent( &m_RingsFilename),
        command::CommandState::IsInStaticInitialization() ?
            "" :
            RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "ring_libraries/drug_ring_database.simple.aro.sdf.gz"
      );

      parameters.AddInitializer
      (
        "extend_within_prob",
        "probability of extending the molecule from within two bonds. "
        "this grows the molecule from within by fragmenting the molecule "
        "into two pieces, inserting a linker, and re-attaching the original "
        "fragments to the new linker. the remainder of the probability goes "
        "toward extending the molecule outward using a linker and/or ring.",
        io::Serialization::GetAgent( &m_ExtendWithinProb),
        "0.50"
      );

      parameters.AddInitializer
      (
        "min_fragment_size",
        "if extending the fragment internally, the fragment will be discarded unless it contains "
        "at least 'min_fragment_size' atoms + 1. If neither fragment contains a valid size, returns "
        "null and tries the mutate again up until the max attempts are reached. If only one "
        "fragment is of a valid size, behavior is determined by 'allow_fragment_duplication' flag.",
        io::Serialization::GetAgent( &m_FragmentMinSize),
        "0"
      );

      parameters.AddInitializer
      (
        "allow_fragment_duplication",
        "If only one of two fragments is of a valid fragment size as specified via "
        "'min_fragment_size' during a within-fragment extension, enabling this flag allows "
        "the linker to connect a duplicate copy of the single valid fragment. Otherwise, "
        "returns null",
        io::Serialization::GetAgent( &m_AllowFragmentDuplication),
        "false"
      );

      parameters.AddInitializer
      (
        "paired_atoms",
        "atom indices (0-indexed) of a second pool of atoms that can be mutated; "
        "this flag has no effect if there are no valid 'mutable_atoms' for "
        "the chosen mutate; defaults to all atoms in the molecule, which may "
        "have different implications for different mutates; for example, in "
        "RemoveBond the 'paired_atoms' flag can be specified along with "
        "'mutable_atoms' to delete a specific bond between two atoms",
        io::Serialization::GetAgent( &m_PairedAtoms),
        ""
      );

      parameters.AddInitializer
      (
        "amide_link_prob",
        "relative probability of linking fragments with an amide.",
        io::Serialization::GetAgent( &m_AmideLinkProb),
        "0.20"
      );

      parameters.AddInitializer
      (
        "amide_n_attach_prob",
        "probability of AmideLink connecting two fragments such that fragment A "
        "is attached via the amide nitrogen. the remainder of the probability is the likelihood "
        "that fragment A will be attached at the amide carbon.",
        io::Serialization::GetAgent( &m_AmideNToAProb),
        "0.50"
      );

      parameters.AddInitializer
      (
        "ester_link_prob",
        "relative probability of linking fragments with an ester.",
        io::Serialization::GetAgent( &m_EsterLinkProb),
        "0.20"
      );

      parameters.AddInitializer
      (
        "ester_o_attach_prob",
        "probability of EsterLink connecting two fragments such that fragment A "
        "is attached via the acidic ester single-bond oxygen. the remainder of the probability is the likelihood "
        "that fragment A will be attached at the ester carbon.",
        io::Serialization::GetAgent( &m_EsterOToAProb),
        "0.50"
      );

      parameters.AddInitializer
      (
        "methoxy_link_prob",
        "relative probability of linking fragments with a methoxy.",
        io::Serialization::GetAgent( &m_MethoxyLinkProb),
        "0.20"
      );

      parameters.AddInitializer
      (
        "methoxy_o_attach_prob",
        "probability of MethoxyLink connecting two fragments such that fragment A "
        "is attached via the methoxy oxygen. the remainder of the probability is the likelihood "
        "that fragment A will be attached at the methoxy carbon.",
        io::Serialization::GetAgent( &m_MethoxyOToAProb),
        "0.50"
      );

      parameters.AddInitializer
      (
        "ethoxy_link_prob",
        "relative probability of linking fragments with an ethoxy.",
        io::Serialization::GetAgent( &m_EthoxyLinkProb),
        "0.20"
      );

      parameters.AddInitializer
      (
        "ethoxy_o_attach_prob",
        "probability of EthoxyLink connecting two fragments such that fragment A "
        "is attached via the ethoxy oxygen. the remainder of the probability is the likelihood "
        "that fragment A will be attached at the ethoxy carbon.",
        io::Serialization::GetAgent( &m_EthoxyOToAProb),
        "0.50"
      );

      parameters.AddInitializer
      (
        "single_element_link_prob",
        "relative probability of linking fragments with a single element",
        io::Serialization::GetAgent( &m_SingleElementLinkProb),
        "0.20"
      );

      parameters.AddInitializer
      (
        "direct_link_prob",
        "relative probability of linking fragments directly; "
        "note that this does not apply to within-fragment linking",
        io::Serialization::GetAgent( &m_DirectLinkProb),
        "0.20"
      );

      parameters.AddInitializer
      (
        "alkyl_link_prob",
        "relative probability of linking fragments with a carbon chain; "
        "note that it is called 'alkyl' linker, but in actuality it can "
        "also link via alkynes.",
        io::Serialization::GetAgent( &m_AlkylLinkProb),
        "0.20"
      );

      parameters.AddInitializer
      (
        "ring_link_prob",
        "relative probability of linking fragments with a new ring; "
        "no effect if not extending within a fragment",
        io::Serialization::GetAgent( &m_RingLinkProb),
        "0.20"
      );

      parameters.AddInitializer
      (
        "B_prob",
        "relative probability of the SingleElementLinker linking with Boron; "
        "no effect if SingleElementLink is not called",
        io::Serialization::GetAgent( &m_B),
        "0.1"
      );

      parameters.AddInitializer
      (
        "C_prob",
        "relative probability of the SingleElementLinker linking with Carbon; "
        "no effect if SingleElementLink is not called",
        io::Serialization::GetAgent( &m_C),
        "0.1"
      );

      parameters.AddInitializer
      (
        "O_prob",
        "relative probability of the SingleElementLinker linking with Oxygen; "
        "no effect if SingleElementLink is not called",
        io::Serialization::GetAgent( &m_O),
        "0.1"
      );

      parameters.AddInitializer
      (
        "N_prob",
        "relative probability of the SingleElementLinker linking with Nitrogen; "
        "no effect if SingleElementLink is not called",
        io::Serialization::GetAgent( &m_N),
        "0.1"
      );

      parameters.AddInitializer
      (
        "P_prob",
        "relative probability of the SingleElementLinker linking with Phosphorous; "
        "no effect if SingleElementLink is not called",
        io::Serialization::GetAgent( &m_P),
        "0.1"
      );

      parameters.AddInitializer
      (
        "S_prob",
        "relative probability of the SingleElementLinker linking with Sulfur; "
        "no effect if SingleElementLink is not called",
        io::Serialization::GetAgent( &m_S),
        "0.1"
      );

      parameters.AddInitializer
      (
        "Se_prob",
        "relative probability of the SingleElementLinker linking with Selenium; "
        "no effect if SingleElementLink is not called",
        io::Serialization::GetAgent( &m_Se),
        "0.1"
      );

//      parameters.AddInitializer
//      (
//        "alkyl_three_C_prob",
//        "relative probability of linking with 3 carbon atoms instead of 2; "
//        "no effect if AlkylLinker is not called",
//        io::Serialization::GetAgent( &m_ThreeC),
//        "0.0"
//      );

      parameters.AddInitializer
      (
        "alkyne_prob",
        "probability of linking via an alkyne instead of an alkane; "
        "no effect if AlkylLinker is not called",
        io::Serialization::GetAgent( &m_Alkyne),
        "0.5"
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentMutateExtendWithLinker::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
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

      // read in ring library
      if( m_RingsFilename.size())
      {
        s_Mutex.Lock();
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_RingsFilename);
        FragmentEnsemble ensemble( input, sdf::e_Remove);
        m_Rings = util::CloneToShPtr( ensemble);
        io::File::CloseClearFStream( input);
        s_Mutex.Unlock();
        BCL_Assert( m_Rings->GetSize(), "Ring library is empty!");
        s_Mutex.Lock();
        s_RingLibraries.Insert( std::make_pair( m_RingsFilename, m_Rings));
        s_Mutex.Unlock();
      }

      // read in paired atom indices
      if( m_PairedAtoms.size())
      {
        m_PairedAtomIndices = util::SplitStringToNumerical< size_t>( m_PairedAtoms);
      }

      // done
      return true;
    }

  } // namespace chemistry
} // namespace bcl
