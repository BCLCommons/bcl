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
#include "chemistry/bcl_chemistry_fragment_mutate_alchemy.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_mutate_remove_bond.h"
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "find/bcl_find_collector_interface.h"
#include "graph/bcl_graph_common_subgraph_isomorphism_base.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "random/bcl_random_uniform_distribution.h"
#include "util/bcl_util_string_functions.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentMutateAlchemy::s_Instance
    (
      util::Enumerated< FragmentMutateInterface>::AddInstance( new FragmentMutateAlchemy())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentMutateAlchemy::FragmentMutateAlchemy() :
        m_AllowedElements( storage::Vector< ElementType>()),
        m_AllowedElementsString( "H C O N S"),
        m_FormalCharge( util::GetUndefined< float>()),
        m_Chirality( Chirality::e_UnknownChirality),
        m_RestrictToBondedH( false)
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief druglikeness constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    FragmentMutateAlchemy::FragmentMutateAlchemy
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const bool &CORINA_CONFS
    ) :
      m_AllowedElements( storage::Vector< ElementType>()),
      m_AllowedElementsString( "H C O N S"),
      m_FormalCharge( util::GetUndefined< float>()),
      m_Chirality( Chirality::e_UnknownChirality),
      m_RestrictToBondedH( false)
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief secondary constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    FragmentMutateAlchemy::FragmentMutateAlchemy
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const bool &CORINA_CONFS
    ) :
      m_AllowedElements( storage::Vector< ElementType>()),
      m_AllowedElementsString( "H C O N S"),
      m_FormalCharge( util::GetUndefined< float>()),
      m_Chirality( Chirality::e_UnknownChirality),
      m_RestrictToBondedH( false)
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_ScaffoldFragment = SCAFFOLD_FRAGMENT;
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief local mutate constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    FragmentMutateAlchemy::FragmentMutateAlchemy
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const bool &CORINA_CONFS
    ) :
      m_AllowedElements( storage::Vector< ElementType>()),
      m_AllowedElementsString( "H C O N S"),
      m_FormalCharge( util::GetUndefined< float>()),
      m_Chirality( Chirality::e_UnknownChirality),
      m_RestrictToBondedH( false)
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
    FragmentMutateAlchemy::FragmentMutateAlchemy
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
      m_AllowedElements( storage::Vector< ElementType>()),
      m_AllowedElementsString( "H C O N S"),
      m_FormalCharge( util::GetUndefined< float>()),
      m_Chirality( Chirality::e_UnknownChirality),
      m_RestrictToBondedH( false)
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
    FragmentMutateAlchemy::FragmentMutateAlchemy
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
      m_AllowedElements( storage::Vector< ElementType>()),
      m_AllowedElementsString( "H C O N S"),
      m_FormalCharge( util::GetUndefined< float>()),
      m_Chirality( Chirality::e_UnknownChirality),
      m_RestrictToBondedH( false)
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
    FragmentMutateAlchemy *FragmentMutateAlchemy::Clone() const
    {
      return new FragmentMutateAlchemy( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentMutateAlchemy::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentMutateAlchemy::GetAlias() const
    {
      static const std::string s_name( "Alchemy");
      return s_name;
    }

    //! @brief returns the element type chosen during the mutate
    //! @return the element type chosen during the mutate; if the
    //! mutate has not yet been run, this will return an undefined
    //! element type object, which is different than the element
    //! type notated 'X' for undefined.
    const ElementType &FragmentMutateAlchemy::GetChosenElementType() const
    {
      return m_ChosenElementType;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @return MutateResult with Constitution after the mutate
    math::MutateResult< FragmentComplete> FragmentMutateAlchemy::operator()( const FragmentComplete &FRAGMENT) const
    {
      BCL_MessageStd( "FragmentMutateAlchemy!");

      // try this a few times
      for( size_t counter( 0); counter < m_NumberMaxAttempts; ++counter)
      {

        // pick random atom to transform
        util::SiPtr< const AtomConformationalInterface> picked_atom;
        if( m_MutableAtomIndices.GetSize() || m_MutableElements.GetSize() || m_MutableFragments.GetSize())
        {
          picked_atom = this->PickAtom( FRAGMENT, false);
        }
        else
        {
          picked_atom = this->PickAtom( FRAGMENT, true);
        }

        // if restricted to hydrogen atoms get the new picked atom
        if( m_RestrictToBondedH && picked_atom->GetElementType() != GetElementTypes().e_Hydrogen)
        {
          // loop over bonds and find the hydrogen atom indices
          storage::Vector< size_t> h_indices;
          for
          (
              auto bond_itr( picked_atom->GetBonds().Begin()),
              bond_itr_end( picked_atom->GetBonds().End());
              bond_itr != bond_itr_end;
              ++bond_itr
          )
          {
            if( bond_itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Hydrogen)
            {
              h_indices.PushBack( FRAGMENT.GetAtomVector().GetAtomIndex( bond_itr->GetTargetAtom()));
            }
          }

          // require some hydrogen atoms that can be mutated
          if( !h_indices.GetSize())
          {
            continue;
          }
          // set new picked atom
          else if( h_indices.GetSize() > size_t( 1))
          {
            h_indices.Shuffle();
          }
          picked_atom = util::SiPtr< const AtomConformationalInterface>( FRAGMENT.GetAtomVector()( h_indices( 0)));
        }

        AtomType new_atom_type( picked_atom->GetAtomType());
        ElementType starting_ele( picked_atom->GetElementType());

        // TODO: consider wrapping these bits of logic regarding ring substitution in a serializable bool for all mutates
        // if hydrogen atom, check if in aromatic ring
//        if( picked_atom->GetElementType() == GetElementTypes().e_Hydrogen)
//        {
//          // check if bonded to something in an aromatic ring
//          util::SiPtr< const AtomConformationalInterface> bonded_atom( picked_atom->GetBonds().Begin()->GetTargetAtom());
//          for
//          (
//              auto bond_itr( bonded_atom->GetBonds().Begin()), bond_itr_end( bonded_atom->GetBonds().End());
//              bond_itr != bond_itr_end;
//              ++bond_itr
//          )
//          {
//            // if bonded to something in aromatic ring, then we need to obey ring activation/deactivation rules
//            if( bond_itr->GetTargetAtom().CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAromatic, size_t( 1)))
//            {
//              if( !IsRingSubstitutionDirected( FRAGMENT, picked_atom))
//              {
//                // if aromatic and not correctly substituted, skip mutate
//                return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
//              }
//            }
//          }
//        }

        // count number of allowed elements and probabilities
        size_t n_possible_ele_types( m_AllowedElements.GetSize());
        storage::Vector< float> ele_swap_types_probs( n_possible_ele_types, 1.0 / n_possible_ele_types);

        // obtain a swap type
        float rand( random::GetGlobalRandom().Random< float>( 0, 1));
        float cumulative_weighted_sum( 0.0);
        size_t rand_ele_index( 0);
        for( ; rand_ele_index < n_possible_ele_types; ++rand_ele_index)
        {
          cumulative_weighted_sum += ele_swap_types_probs( rand_ele_index);
          if( rand < cumulative_weighted_sum)
          {
            break;
          }
        }

        // remember the chosen element type
        SetChosenElement( m_AllowedElements( rand_ele_index));

        // is our atom in an aromatic ring?
        bool picked_atom_aromatic( picked_atom->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAromatic, 1));
        if( GetChosenElementType() == GetElementTypes().e_Hydrogen)
        {
          new_atom_type = GetAtomTypes().H_S;
        }
        else if( GetChosenElementType() == GetElementTypes().e_Undefined)
        {
          new_atom_type = GetAtomTypes().e_Undefined;
        }
        else
        {
          // find number of bonds being made by starting atom
          size_t n_h( picked_atom->GetNumberCovalentlyBoundHydrogens());
          size_t n_nonh_bonds( picked_atom->GetAtomType()->GetNumberBonds() - n_h);
          for( size_t n_current_bonds( n_nonh_bonds); n_current_bonds <= 6; ++n_current_bonds)
          {
            size_t n_nonh_e( picked_atom->GetAtomType()->GetNumberElectronsInBonds() - n_h + ( n_current_bonds - n_nonh_bonds));
            BCL_MessageDbg( "n_h: " + util::Format()( n_h));
            BCL_MessageDbg( "n_nonh_bonds: " + util::Format()( n_current_bonds));
            BCL_MessageDbg( "n_e_bonds: " + util::Format()( picked_atom->GetAtomType()->GetNumberElectronsInBonds()));
            BCL_MessageDbg( "n_nonh_e: " + util::Format()( n_nonh_e));
            PossibleAtomTypesForAtom available_atom_types
            (
              GetChosenElementType(),
              n_nonh_e,
              n_current_bonds,
              util::IsDefined( m_FormalCharge) ? m_FormalCharge : picked_atom->GetCharge(),
              picked_atom_aromatic
            );
            // find something with the requested formal charge
            if( available_atom_types.GetNumberPossibleTypes() && util::IsDefined( m_FormalCharge))
            {
              if( available_atom_types.GetAlternateTypeWithCharge( m_FormalCharge).IsDefined())
              {
                new_atom_type = available_atom_types.GetAlternateTypeWithCharge( m_FormalCharge);
                break;
              }
            }
            // prefer original charge
            else if( available_atom_types.GetNumberPossibleTypes() && available_atom_types.GetMostStableType()->IsGasteigerAtomType())
            {
              if( available_atom_types.GetAlternateTypeWithCharge( picked_atom->GetCharge()).IsDefined())
              {
                new_atom_type = available_atom_types.GetAlternateTypeWithCharge( picked_atom->GetCharge());
                break;
              }
              new_atom_type = available_atom_types.GetMostStableType();
              break;
            }
          }
        }

        // Change the element type of our picked atom
        size_t picked_atom_index( FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom));
        AtomVector< AtomComplete> new_atom_vector( FRAGMENT.GetAtomVector());
        new_atom_vector( picked_atom_index).SetAtomType( new_atom_type);

        // set chirality if it is not unknown
        if( m_Chirality != Chirality::e_UnknownChirality)
        {
          new_atom_vector( picked_atom_index).SetChirality( m_Chirality);
        }

        // make sure the formal charge is consistent with the desired formal charge
        float fc( new_atom_vector( picked_atom_index).GetAtomType()->GetFormalCharge());
        if( util::IsDefined( m_FormalCharge))
        {
          if( fc != m_FormalCharge)
          {
            return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
          }
        }

        // TODO check if this can be removed; not really doing much
        // maybe make this an option
        if( fc != 0.0 && !util::IsDefined( m_FormalCharge))
        {
          new_atom_vector( picked_atom_index).SetCharge( 0.0);
        }

        // TODO same here; check if this can be removed; not really doing much
        // neutralize the new atom and reset the bond conjugation; the cleaner will re-aromatize if needed
        if
        (
            (
                fc != 0.0 &&
                picked_atom_aromatic &&
                new_atom_vector( picked_atom_index).GetElementType() != picked_atom->GetElementType()
            ) ||
            // we also want to make sure we can change bond orders to accomodate hydrogen atoms subbing in for a heavy atom
            new_atom_type == GetAtomTypes().H_S
        )
        {
          // find connected aromatic bonds and convert them to non-aromatic bond order
          size_t bond_index( 0);
          for
          (
              auto bond_itr( new_atom_vector( picked_atom_index).GetBonds().Begin()),
              bond_itr_end( new_atom_vector( picked_atom_index).GetBonds().End());
              bond_itr != bond_itr_end;
              ++bond_itr, ++bond_index
          )
          {
            if
            (
                (
                    bond_itr->GetBondType() == GetConfigurationalBondTypes().e_AromaticSingleBond ||
                    bond_itr->GetBondType() == GetConfigurationalBondTypes().e_AromaticDoubleBond
                ) ||
                (
                    new_atom_type == GetAtomTypes().H_S &&
                    bond_itr->GetBondType() != GetConfigurationalBondTypes().e_NonConjugatedSingleBond
                )
            )
            {
              size_t i( new_atom_vector.GetAtomIndex( new_atom_vector( picked_atom_index).GetBonds()( bond_index).GetTargetAtom()));
              new_atom_vector( picked_atom_index).SetBondTypeTo( new_atom_vector( i), GetConfigurationalBondTypes().e_NonConjugatedSingleBond);
            }
          }
        }

        // skip geometry cleanup if we mutated into an undefined atom
        if( GetChosenElementType() == GetElementTypes().e_Undefined)
        {
          BCL_MessageStd( "Chosen element type from allowed element type list is Undefined; skipping fragment cleaning");
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>( new FragmentComplete( new_atom_vector, "")), *this);
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
//          storage::Vector< size_t>(),
//          false,
//          false,
//          4
        );

        // remove hydrogen atoms to ease burden on the isomorphism search during cleaning
        HydrogensHandler::Remove( new_atom_vector);

        // Standardize and return
        if( m_ScaffoldFragment.GetSize())
        {
          return math::MutateResult< FragmentComplete>( cleaner.Clean( new_atom_vector, m_ScaffoldFragment, m_DrugLikenessType), *this);
        }
        else
        {
          return math::MutateResult< FragmentComplete>( cleaner.Clean( new_atom_vector, FRAGMENT, m_DrugLikenessType), *this);
        }
      }

      // if no luck, return null ptr
      return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set the fragment mutable atom indices
    void FragmentMutateAlchemy::SetAllowedElements( const storage::Vector< ElementType> &ALLOWED_ELEMENTS)
    {
      m_AllowedElements = ALLOWED_ELEMENTS;
    }

    //! @brief set the fragment mutable atom indices
    void FragmentMutateAlchemy::SetRestrictions( const bool RESTRICT_TO_BOND_H)
    {
      m_RestrictToBondedH = RESTRICT_TO_BOND_H;
    }

    //! @brief set the chosen element type to which we are mutating
    void FragmentMutateAlchemy::SetChosenElement( const ElementType &ELEMENT_TYPE) const
    {
      m_ChosenElementType = ELEMENT_TYPE;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    io::Serializer FragmentMutateAlchemy::GetSerializer() const
    {
      io::Serializer parameters( FragmentMutateInterface::GetSerializer());
      parameters.SetClassDescription
      (
        "Transforms atom/element types inside of a molecule"
      );

      parameters.AddInitializer
      (
        "allowed_elements",
        "elements that are accessible via this mutate",
        io::Serialization::GetAgent( &m_AllowedElementsString),
        "H C O N S"
      );

      parameters.AddInitializer
      (
        "set_formal_charge",
        "explicitly set the formal charge of an atom; if the charged species does not exist as a valid atom "
        "type in the context of the rest of the molecule then the mutate fails.",
        io::Serialization::GetAgent( &m_FormalCharge),
        "nan"
      );

      parameters.AddInitializer
      (
        "set_chirality",
        "explicitly set the chirality of an atom; no effect on atoms that are not at a chiral center; "
        "default is to UnknownChirality; if chirality is unknown and unspecified, will default to "
        "trying to determine based on the input geometry and new atom type",
        io::Serialization::GetAgent( &m_Chirality),
        "?"
      );

      parameters.AddInitializer
      (
        "restrict_to_bonded_h",
        "If false, this mutate does as expected and converts an allowed atom into an appropriate atom type "
        "for the list of allowed elements at a mutable atom index. If true, then only hydrogen atoms that are "
        "bonded to mutable heavy atoms (or mutable hydrogen atoms) will be mutated. The reason for this flag "
        "is that hydrogen atoms can be difficult to track. Generally, the mutates strip and re-add hydrogen "
        "atoms during the cleaning phase, which means that hydrogen atom indices change at a much easier than "
        "heavy atom indices. So, if you are using multiple mutates sequentially and you want to perturb hydrogen "
        "atoms directly with FragmentMutateAlchemy, consider doing it by specifying the bonded heavy atom index as "
        "mutable and setting this flag true.",
        io::Serialization::GetAgent( &m_RestrictToBondedH),
        "false"
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentMutateAlchemy::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
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

      // read in allowed elements
      if( m_AllowedElementsString.size())
      {
        // parse input
        const storage::Vector< std::string> allowed_elements
        (
          util::SplitString( util::TrimString( m_AllowedElementsString), " \t\n\r,")
        );

        // stupid check to add only the correct elements
        // TODO: make this more directly serializable from element types
        m_AllowedElements.Reset();
        for( size_t e_i( 0), e_sz( allowed_elements.GetSize()); e_i < e_sz; ++e_i)
        {
          // Hydrogen
          if( allowed_elements( e_i) == "H")
          {
            m_AllowedElements.PushBack( GetElementTypes().e_Hydrogen);
          }
          // Boron
          if( allowed_elements( e_i) == "B")
          {
            m_AllowedElements.PushBack( GetElementTypes().e_Boron);
          }
          // Carbon
          if( allowed_elements( e_i) == "C")
          {
            m_AllowedElements.PushBack( GetElementTypes().e_Carbon);
          }
          // Oxygen
          if( allowed_elements( e_i) == "O")
          {
            m_AllowedElements.PushBack( GetElementTypes().e_Oxygen);
          }
          // Nitrogen
          if( allowed_elements( e_i) == "N")
          {
            m_AllowedElements.PushBack( GetElementTypes().e_Nitrogen);
          }
          // Phosphorous
          if( allowed_elements( e_i) == "P")
          {
            m_AllowedElements.PushBack( GetElementTypes().e_Phosphorus);
          }
          // Sulfur
          if( allowed_elements( e_i) == "S")
          {
            m_AllowedElements.PushBack( GetElementTypes().e_Sulfur);
          }
          // Selenium
          if( allowed_elements( e_i) == "Se")
          {
            m_AllowedElements.PushBack( GetElementTypes().e_Selenium);
          }
          // Fluorine
          if( allowed_elements( e_i) == "F")
          {
            m_AllowedElements.PushBack( GetElementTypes().e_Fluorine);
          }
          // Chlorine
          if( allowed_elements( e_i) == "Cl")
          {
            m_AllowedElements.PushBack( GetElementTypes().e_Chlorine);
          }
          // Bromine
          if( allowed_elements( e_i) == "Br")
          {
            m_AllowedElements.PushBack( GetElementTypes().e_Bromine);
          }
          // Iodine
          if( allowed_elements( e_i) == "I")
          {
            m_AllowedElements.PushBack( GetElementTypes().e_Iodine);
          }
          // Undefined
          if( allowed_elements( e_i) == "X")
          {
            m_AllowedElements.PushBack( GetElementTypes().e_Undefined);
          }
        }
      }

      // done
      return true;
    }

  } // namespace chemistry
} // namespace bcl
