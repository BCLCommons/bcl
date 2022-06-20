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
#include "chemistry/bcl_chemistry_fragment_mutate_fluorinate.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "command/bcl_command_command_state.h"
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
    const util::SiPtr< const util::ObjectInterface> FragmentMutateFluorinate::s_Instance
    (
      util::Enumerated< FragmentMutateInterface>::AddInstance( new FragmentMutateFluorinate())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentMutateFluorinate::FragmentMutateFluorinate() :
      m_Reversible( false)
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief druglikeness constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    FragmentMutateFluorinate::FragmentMutateFluorinate
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const bool &CORINA_CONFS
    ) :
      m_Reversible( false)
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief full constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    FragmentMutateFluorinate::FragmentMutateFluorinate
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const bool &CORINA_CONFS
    ) :
      m_Reversible( false)
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
    FragmentMutateFluorinate::FragmentMutateFluorinate
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
      m_Reversible( false)
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
    FragmentMutateFluorinate::FragmentMutateFluorinate
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
      m_Reversible( false)
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
    FragmentMutateFluorinate *FragmentMutateFluorinate::Clone() const
    {
      return new FragmentMutateFluorinate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentMutateFluorinate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentMutateFluorinate::GetAlias() const
    {
      static const std::string s_name( "Fluorinate");
      return s_name;
    }

    //! @brief returns whether fluorinate is reversible
    //! @return reversibility of the fluorinate mutate
    const bool FragmentMutateFluorinate::GetReversibility() const
    {
      return m_Reversible;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @return MutateResult with Constitution after the mutate
    math::MutateResult< FragmentComplete> FragmentMutateFluorinate::operator()( const FragmentComplete &FRAGMENT) const
    {
      BCL_MessageStd( "Fluorinate!");
      AtomVector< AtomComplete> atom_vector( FRAGMENT.GetAtomVector());

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

      // make equal probability to add or remove fluorine atoms
      float rand( 0.0);
      if( m_Reversible)
      {
        rand = random::GetGlobalRandom().Random< float>( 0.0, 1.0);
      }

      // add fluorine atoms
      if( rand < 0.50)
      {
        BCL_MessageStd( "Adding fluorine atoms!");
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

          size_t picked_atom_index( FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom));

          // TODO: consider making a flag for this
          // make sure not aromatic to avoid double-counting with Halogenate
//          if( picked_atom->GetElementType() == GetElementTypes().e_Hydrogen)
//          {
//            util::SiPtr< const AtomConformationalInterface> bonded_atom( picked_atom->GetBonds().Begin()->GetTargetAtom());
//            size_t n_aro_neigh( bonded_atom->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAromatic, 1));
//            if( n_aro_neigh)
//            {
//              return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
//            }
//          }
//          else
//          {
//            size_t n_aro_neigh( picked_atom->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAromatic, 1));
//            if( n_aro_neigh)
//            {
//              return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
//            }
//          }

          // find a carbon or nitrogen with some hydrogen atoms
          size_t n_added_f( 0), n_removed_h( 0);
          if( picked_atom->GetElementType() == GetElementTypes().e_Carbon) // || picked_atom->GetElementType() == GetElementTypes().e_Nitrogen)
          {
            for( auto bonds_itr( atom_vector( picked_atom_index).GetBonds().Begin()), bonds_itr_end( atom_vector( picked_atom_index).GetBonds().End()); bonds_itr != bonds_itr_end; ++bonds_itr)
            {
              if( bonds_itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Hydrogen)
              {
                // convert to fluorine
                atom_vector( atom_vector.GetAtomIndex( bonds_itr->GetTargetAtom())).SetAtomType( GetAtomTypes().F_SP2P2P2);
                ++n_added_f, ++n_removed_h;
              }
              else if( bonds_itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Fluorine && m_IncludeExistingFInMinCount)
              {
                ++n_added_f;
              }
            }

            // check min f req
            if( n_added_f < m_MinF || n_removed_h < m_MinH || n_removed_h > m_MaxH)
            {
              atom_vector = FRAGMENT.GetAtomVector();
              continue;
            }

            HydrogensHandler::Remove( atom_vector);
            util::ShPtr< FragmentComplete> new_mol_ptr
            (
              m_ScaffoldFragment.GetSize() ?
              cleaner.Clean( atom_vector, m_ScaffoldFragment, m_DrugLikenessType) :
              cleaner.Clean( atom_vector, FRAGMENT, m_DrugLikenessType)
            );
            if( !new_mol_ptr.IsDefined())
            {
              continue;
            }
            return math::MutateResult< FragmentComplete>( new_mol_ptr, *this);
          }
        }
      }
      // remove fluorine atoms
      else
      {
        BCL_MessageStd( "Removing fluorine atoms!");
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

          size_t picked_atom_index( FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom));

          // find a carbon or nitrogen with some fluorine atoms
          size_t n_removed_f( 0), n_added_h( 0);
          if( picked_atom->GetElementType() == GetElementTypes().e_Carbon) // || picked_atom->GetElementType() == GetElementTypes().e_Nitrogen)
          {
            for( auto bonds_itr( atom_vector( picked_atom_index).GetBonds().Begin()), bonds_itr_end( atom_vector( picked_atom_index).GetBonds().End()); bonds_itr != bonds_itr_end; ++bonds_itr)
            {
              if( bonds_itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Fluorine)
              {
                // convert to hydrogen
                atom_vector( atom_vector.GetAtomIndex( bonds_itr->GetTargetAtom())).SetAtomType( GetAtomTypes().H_S);
                ++n_removed_f, ++n_added_h;
              }
            }

            // check min f req
            if( n_removed_f < m_MinF || n_added_h < m_MinH || n_added_h > m_MaxH)
            {
              atom_vector = FRAGMENT.GetAtomVector();
              continue;
            }

            HydrogensHandler::Remove( atom_vector);
            util::ShPtr< FragmentComplete> new_mol_ptr
            (
              m_ScaffoldFragment.GetSize() ?
              cleaner.Clean( atom_vector, m_ScaffoldFragment, m_DrugLikenessType) :
              cleaner.Clean( atom_vector, FRAGMENT, m_DrugLikenessType)
            );
            if( !new_mol_ptr.IsDefined())
            {
              continue;
            }
            return math::MutateResult< FragmentComplete>( new_mol_ptr, *this);
          }
        }
      }
      // return null if we do not identify any suitable carbon or nitrogen atoms
      return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set reversibility
    void FragmentMutateFluorinate::SetReverisibility( const bool REVERSIBLE)
    {
      m_Reversible = REVERSIBLE;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    io::Serializer FragmentMutateFluorinate::GetSerializer() const
    {
      io::Serializer parameters( FragmentMutateInterface::GetSerializer());
      parameters.SetClassDescription
      (
        "Add fluorine atoms to molecules; distinct from "
        "halogenate because patterns of fluorine placement "
        "in organic molecules tend to differ from patterns "
        "of bulkier halogens"
      );

      parameters.AddInitializer
      (
        "reversible",
        "sets the probability of removing fluorine atoms to 50%; "
        "default false indicates that fluorine can only be added",
        io::Serialization::GetAgent( &m_Reversible),
        "false"
      );

      parameters.AddInitializer
      (
        "n_min_f",
        "sets the minimum number of fluorine atoms that must be added or removed from a single heavy atom for a successful mutate; "
        "if this many fluorines fail to be added or removed then the mutate will try again up to n_max_mutates "
        "number of times; if failure occurs beyond the maximum number of attempts then this mutate returns null; "
        "this option may be useful to guide mutable atom selection; "
        "TODO: expand this to optionally cover multiple mutable atoms instead of just a single selected mutable atom",
        io::Serialization::GetAgent( &m_MinF),
        "0"
      );

      parameters.AddInitializer
      (
        "include_existing_f_in_min_count",
        "when counting the number of fluorine atoms added, "
        "include in that count fluorine atoms already attached to the chosen heavy atom",
        io::Serialization::GetAgent( &m_IncludeExistingFInMinCount),
        "false"
      );

      parameters.AddInitializer
      (
        "n_min_h_sub",
        "sets the minimum number of hydrogen atoms that must be replaced on a single heavy atom for a successful mutate; "
        "TODO: expand this to optionally cover multiple mutable atoms instead of just a single selected mutable atom",
        io::Serialization::GetAgent( &m_MinH),
        "0"
      );

      parameters.AddInitializer
      (
        "n_max_h_sub",
        "sets the maximum number of hydrogen atoms that can be replaced on a single heavy atom for a successful mutate; "
        "TODO: expand this to optionally cover multiple mutable atoms instead of just a single selected mutable atom",
        io::Serialization::GetAgent( &m_MaxH),
        "6"
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentMutateFluorinate::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
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

      // done
      return true;
    }

  } // namespace chemistry
} // namespace bcl
