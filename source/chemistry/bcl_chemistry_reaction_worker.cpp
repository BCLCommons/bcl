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
#include "chemistry/bcl_chemistry_reaction_worker.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_clash_score.h"
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_angle_assignment.h"
#include "chemistry/bcl_chemistry_bond_isometry_handler.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_align_to_scaffold.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"
#include "chemistry/bcl_chemistry_mutate_bond_angles.h"
#include "chemistry/bcl_chemistry_mutate_bond_lengths.h"
#include "chemistry/bcl_chemistry_mutate_clash_resolver.h"
#include "chemistry/bcl_chemistry_mutate_dihedrals_interface.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "chemistry/bcl_chemistry_search_fragment_library_from_tree.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_comparisons.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! @param ATOM_COMPARISON how to compare atoms in structures
    //! @param BOND_COMPARISON how to compare bonds in structures
    ReactionWorker::ReactionWorker
    (
      const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON,
      const ConfigurationalBondTypeData::Data &BOND_COMPARISON
    ) :
      m_ReactantCache(),
      m_ProductCache(),
      m_Reactants( FragmentEnsemble()),
      m_Reference( FragmentComplete()),
      m_ReferenceReactantIndex( size_t()), // intentionally undefined
      m_ReactantsReactiveAtoms( storage::Vector< storage::Map< size_t, size_t> >()),
      m_ProductConformerArbitrary( true),
      m_CorrectGeometry( false),
      m_CorrectNonReferenceRingGeometry( false),
      m_AdditionalAdjacentAtoms( size_t( 0)) // intentionally 0
    {
    }

    //! @brief copy constructor; does not copy caches
    ReactionWorker::ReactionWorker
    (
      const ReactionWorker &OTHER
    ) :
      m_ReactantCache(),
      m_ProductCache(),
      m_Reactants( FragmentEnsemble()),
      m_Reference( FragmentComplete()),
      m_ReferenceReactantIndex( size_t()), // intentionally undefined
      m_ReactantsReactiveAtoms( storage::Vector< storage::Map< size_t, size_t> >()),
      m_ProductConformerArbitrary( OTHER.GetIsProductConformerArbitrary()),
      m_CorrectGeometry( OTHER.GetIsCorrectGeometry()),
      m_CorrectNonReferenceRingGeometry( OTHER.GetIsCorrectNonReferenceRingGeometry()),
      m_AdditionalAdjacentAtoms( OTHER.GetAdditionalAdjacentAtoms())
    {
    }

    //! @brief clone constructor
    //! @return a pointer to a copy of this class
    ReactionWorker *ReactionWorker::Clone() const
    {
      return new ReactionWorker( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ReactionWorker::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the reactant cache
    const storage::Map< util::SiPtr< const ReactionComplete>, storage::Vector< ReactionStructure> >
    &ReactionWorker::GetReactantCache() const
    {
      return m_ReactantCache;
    }

    //! @brief return the product cache
    const storage::Map< util::SiPtr< const ReactionComplete>, storage::Vector< ReactionStructure> >
    &ReactionWorker::GetProductCache() const
    {
      return m_ProductCache;
    }

    //! @brief get the passed reagents
    const FragmentEnsemble &ReactionWorker::GetReactants() const
    {
      return m_Reactants;
    }

    //! @brief get the target fragment used to initialize the reaction
    const FragmentComplete &ReactionWorker::GetReference() const
    {
      return m_Reference;
    }

    //! @brief get the reference target fragment reactant index
    const size_t &ReactionWorker::GetReferenceReactantIndex() const
    {
      return m_ReferenceReactantIndex;
    }

    //! @brief get the maps relating reactants atom indices to reactive atoms
    const storage::Vector< storage::Map< size_t, size_t> > &ReactionWorker::GetReactantsReactiveAtoms() const
    {
      return m_ReactantsReactiveAtoms;
    }

    //! @brief return whether or not an arbitrary 3D conformer will be generated for final product(s)
    const bool &ReactionWorker::GetIsProductConformerArbitrary() const
    {
      return m_ProductConformerArbitrary;
    }

    //! @brief return whether or not an arbitrary 3D conformer will be generated for final product(s)
    const bool &ReactionWorker::GetIsCorrectGeometry() const
    {
      return m_CorrectGeometry;
    }

    //! @brief return whether or not an arbitrary 3D conformer will be generated for final product(s)
    const bool &ReactionWorker::GetIsCorrectNonReferenceRingGeometry() const
    {
      return m_CorrectNonReferenceRingGeometry;
    }

    //! @brief return whether or not an arbitrary 3D conformer will be generated for final product(s)
    const size_t &ReactionWorker::GetAdditionalAdjacentAtoms() const
    {
      return m_AdditionalAdjacentAtoms;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief assignment operator
    ReactionWorker &ReactionWorker::operator =( const ReactionWorker &OTHER)
    {
      m_ReactantCache.Reset();
      m_ProductCache.Reset();
      return *this;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief searches reactions for those that MOLECULE could satisfy as a reactant
    //! @param MOLECULE the query molecule
    //! @param REACTION the reaction to search
    //! @param REACTANT_NUMS indices of reactants of interest (empty is all reactants)
    //! @return a vector of the reactant indices that match MOLECULE
    storage::Vector< size_t> ReactionWorker::MatchesReactants
    (
      const ConformationInterface &MOLECULE,
      const ReactionComplete &REACTION,
      const storage::Vector< size_t> &REACTANT_NUMS
    ) const
    {
      size_t n_specified_reactants( REACTANT_NUMS.GetSize());

      // local storage of matching reactants
      storage::Vector< size_t> reactant_indices;

      // Worker for matching reactants in this reaction
      const storage::Vector< ReactionStructure> rxt_structs( GetReactantStructures( REACTION));

      // Search the molecule for each reactant structure
      for( size_t r_no( 0), end_no( rxt_structs.GetSize()); r_no < end_no; ++r_no)
      {

        // Skip if the graph is not defined or if the reactant was not specified for search
        if( n_specified_reactants && REACTANT_NUMS.Find( r_no) >= n_specified_reactants)
        {
          continue;
        }

        if( rxt_structs( r_no).ContainedIn( MOLECULE))
        {
          reactant_indices.PushBack( r_no);
        }
      }
      return reactant_indices;
    }

    //! @brief searches reactions for those that MOLECULE could be a product
    //! @param MOLECULE the molecule to use as a worker
    //! @param REACTION the reaction to search
    //! @return a vector of the product indices that match MOLECULE
    storage::Vector< size_t> ReactionWorker::MatchesProducts
    (
      const ConformationInterface &MOLECULE,
      const ReactionComplete &REACTION,
      const storage::Vector< size_t> &PRODUCT_NUMS
    ) const
    {
      size_t n_specified_products( PRODUCT_NUMS.GetSize());

      // local storage of matching products
      storage::Vector< size_t> product_indices;

      // Worker for matching products in this reaction
      const storage::Vector< ReactionStructure> rxt_structs( GetProductStructures( REACTION));

      // Search the molecule for each product structure
      for( size_t prod_no( 0); prod_no < rxt_structs.GetSize(); ++prod_no)
      {

        // Skip if the graph is not defined or if the product was not specified for search
        if( n_specified_products && PRODUCT_NUMS.Find( prod_no) >= n_specified_products)
        {
          continue;
        }

        if( rxt_structs( prod_no).ContainedIn( MOLECULE))
        {
          product_indices.PushBack( prod_no);
        }
      }
      return product_indices;
    }

    //! @brief searches reactions for those that MOLECULE could satisfy as a reactant
    //! @param MOLECULE the query molecule
    //! @param REACTION the reaction to search
    //! @param INDEX the reactant index against which the molecule is queried
    //! @param REACTANT_NUMS indices of reactants of interest (empty is all reactants)
    //! @return true if the MOLECULE matches the reactant at index INDEX
    bool ReactionWorker::MatchesReactantIndex
    (
      const ConformationInterface &MOLECULE,
      const ReactionComplete &REACTION,
      const size_t &INDEX,
      const storage::Vector< size_t> &REACTANT_NUMS
    ) const
    {
      // match reactants
      storage::Vector< size_t> check_mol_rxt_pos( this->MatchesReactants( MOLECULE, REACTION, REACTANT_NUMS));

      // check index
      if( check_mol_rxt_pos.GetSize() && check_mol_rxt_pos.Find( INDEX) < check_mol_rxt_pos.GetSize())
      {
        return true;
      }
      return false;
    }

    //! @brief Executes a reaction with provided reactants, does not depend on the order reactants are given in
    //! @param REACTION the reaction to use
    //! @param MOLECULES molecules used to do the reaction
    //! @return an ensemble of products
    FragmentEnsemble ReactionWorker::React
    (
      const ReactionComplete &REACTION,
      const FragmentEnsemble &MOLECULES
    ) const // dummy test
    {

      FragmentEnsemble result;

      size_t n_mols( MOLECULES.GetSize());

      // Check that molecules were provided
      if( !n_mols)
      {
        BCL_MessageStd( "No molecules provided");
        return result;
      }

      size_t n_reactants( REACTION.GetNumberReactants());

      // If enough molecules were given to specify every reactant check that they actually match
      // and reorder them so that ExecuteReaction can handle them
      if( n_mols == n_reactants)
      {
        storage::Vector< storage::Vector< size_t> > mol_matches( n_mols);
        storage::Vector< storage::Vector< size_t> > reactant_matches( n_reactants);

        // Determine which reactants match this molecule
        size_t mol_no( 0);
        for
        (
          FragmentEnsemble::const_iterator itr_mol( MOLECULES.Begin()), itr_mol_end( MOLECULES.End());
          itr_mol != itr_mol_end;
          ++itr_mol, ++mol_no
        )
        {
          mol_matches( mol_no) = MatchesReactants( *itr_mol, REACTION);
          if( mol_matches( mol_no).IsEmpty())
          {
            return result;
          }

          for
          (
            storage::Vector< size_t>::const_iterator itr_rxt( mol_matches( mol_no).Begin()), itr_rxt_end( mol_matches( mol_no).End());
            itr_rxt != itr_rxt_end;
            ++itr_rxt
          )
          {
            reactant_matches( *itr_rxt).PushBack( mol_no);
          }
        }

        // make sure each reactant has at least one associated molecule (quick check)
        for
        ( 
          storage::Vector< storage::Vector< size_t> >::const_iterator itr_rxt_match( reactant_matches.Begin()), 
            itr_rxt_match_end( reactant_matches.End());
          itr_rxt_match != itr_rxt_match_end;
          ++itr_rxt_match
        )
        {
          if( itr_rxt_match->IsEmpty())
          {
            return result;
          }
        }

        // Reorder the molecules so that they correspond to the reactants in the given reaction
        storage::Vector< size_t> reorder_mols( n_mols, util::GetUndefined< size_t>());
        storage::Vector< size_t> picked_reactants( n_reactants, size_t( 0));
        if( ChooseValueCombination( mol_matches, 0, picked_reactants, reorder_mols))
        {
          FragmentEnsemble new_order;
          for( size_t i( 0); i < reorder_mols.GetSize(); ++i)
          {
            FragmentEnsemble::const_iterator itr_mol( MOLECULES.Begin());
            std::advance( itr_mol, reorder_mols( i));
            new_order.PushBack( *itr_mol);
          }
          return ExecuteReaction( REACTION, new_order);
        }
      }
      else if
      ( 
        n_mols == 1 
        && MatchesReactants( *MOLECULES.Begin(), REACTION).GetSize() == n_reactants
      )
      {
        return ExecuteIntramolecularReaction( REACTION, *MOLECULES.Begin());
      }

      return result; 
    }

    //! @brief tests if a molecule would be reactive with other molecules of itself under a given reaction; does not count for intramolecular reactions
    //! @param REACTION the reaction to inspect
    //! @param MOLECULE the molecule in question
    //! @return true if the molecule would react with itself, false otherwise
    bool ReactionWorker::IsSelfReactive
    (
      const ReactionComplete &REACTION,
      const FragmentComplete &MOLECULE
    ) const
    {
      const storage::Vector< ReactionStructure> &rxt_structs( GetReactantStructures( REACTION));

      size_t n_matches( 0);
      for
      (
        storage::Vector< ReactionStructure>::const_iterator itr_str( rxt_structs.Begin()), itr_str_end( rxt_structs.End());
        itr_str != itr_str_end;
        ++itr_str
      )
      {
        if( !itr_str->ContainedIn( MOLECULE))
        {
          ++n_matches;
        }
      }
      return n_matches > 1;
    }

    //! @brief set the reagents
    void ReactionWorker::SetReactants( const FragmentEnsemble &REACTANTS) const
    {
      m_Reactants = REACTANTS;
    }

    //! @brief set the target fragment used to initialize the reaction
    void ReactionWorker::SetReference( const FragmentComplete &REFERENCE) const
    {
      m_Reference = REFERENCE;
    }

    //! @brief set the reference target fragment reactant index
    void ReactionWorker::SetReferenceReactantIndex( const size_t REFERENCE_REACTANT_INDEX) const
    {
      m_ReferenceReactantIndex = REFERENCE_REACTANT_INDEX;
    }

    //! @brief set the maps relating reactants atom indices to reactive atoms
    void ReactionWorker::SetReactantsReactiveAtoms
    (
      const storage::Vector< storage::Map< size_t, size_t> > &REACTANTS_REACTIVE_ATOMS
    ) const
    {
      m_ReactantsReactiveAtoms = REACTANTS_REACTIVE_ATOMS;
    }

    //! @brief set whether or not an arbitrary 3D conformer will be generated for final product(s)
    void ReactionWorker::SetProductConformerArbitrary( const bool PRODUCT_CONFORMER_ARBITRARY)
    {
      m_ProductConformerArbitrary = PRODUCT_CONFORMER_ARBITRARY;
    }

    //! @brief set whether or not an arbitrary 3D conformer will be generated for final product(s)
    void ReactionWorker::SetCorrectGeometry( const bool CORRECT_GEOMETRY)
    {
      m_CorrectGeometry = CORRECT_GEOMETRY;
    }

    //! @brief set whether or not an arbitrary 3D conformer will be generated for final product(s)
    void ReactionWorker::SetCorrectNonReferenceRingGeometry( const bool CORRECT_NON_REF_RING_GEOMETRY)
    {
      m_CorrectNonReferenceRingGeometry = CORRECT_NON_REF_RING_GEOMETRY;
    }

    //! @brief set whether or not an arbitrary 3D conformer will be generated for final product(s)
    void ReactionWorker::SetAdditionalAdjacentAtoms( const size_t N_ADDITIONAL_ADJACENT_ATOMS)
    {
      m_AdditionalAdjacentAtoms = N_ADDITIONAL_ADJACENT_ATOMS;
    }

    //! @brief aligns products to mapped reactive atoms
    //! @param PRODUCT the product to be aligned to reactants
    //! @param PRODUCT_REACTIVE_ATOMS map of reactive atoms to atom indices
    //! @param REACTANTS the reactants to which the product is aligned
    //! @param REACTANTS_REACTIVE_ATOMS maps relating atom indices to reactive atoms
    //! @return true if the alignment succeeds
    bool ReactionWorker::AlignProductToReactants
    (
      FragmentComplete &PRODUCT,
      const storage::Map< size_t, size_t> &PRODUCT_REACTIVE_ATOMS,
      const FragmentComplete &REFERENCE,
      const storage::Map< size_t, size_t> &REACTANTS_REACTIVE_ATOMS
    ) const
    {
      // initialize alignment object
      FragmentAlignToScaffold ats
      (
        ConformationGraphConverter::AtomComparisonType::e_ElementType,
        ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness,
        size_t( 2)
      );

      // align product to reactant
      bool result
      (
        ats.AlignToScaffold
        (
          PRODUCT,
          REFERENCE,
          PRODUCT_REACTIVE_ATOMS.GetMappedValues(),
          REACTANTS_REACTIVE_ATOMS.GetKeysAsVector()
        )
      );
      // true if success
      return result;
    }

    //! @brief executes a reaction given a set of molecules which must correspond exactly (in order) to reactants
    //! @param REACTION the reaction to execute
    //! @param MOLECULES the molecules that must match the reactants of REACTION in the same order
    //! @param REFERENCE the target reagent used to initialize the reaction
    //! @param REFERENCE_REACTANT_INDEX the index of the reference reagent within the reaction
    //! @return an ensemble of product molecules
    FragmentEnsemble ReactionWorker::ExecuteReaction
    (
      const ReactionComplete &REACTION,
      const FragmentEnsemble &MOLECULES,
      const FragmentComplete &REFERENCE,
      const size_t &REFERENCE_REACTANT_INDEX
    ) const
    {
      const storage::Vector< ReactionStructure> &rxt_structs( GetReactantStructures( REACTION));
      m_ReferenceReactantIndex = REFERENCE_REACTANT_INDEX;

      // products that will be returned on successful execution (or empty until then)
      FragmentEnsemble result;

      // Reactants must be matched exactly; return if not
      size_t n_molecules( MOLECULES.GetSize());
      if( n_molecules != REACTION.GetNumberReactants())
      {
        BCL_MessageStd
        ( 
          "Reaction was given an inappropriate number of reactants; provided "
          + util::Format()( n_molecules) + ", needed " + util::Format()( REACTION.GetNumberReactants()) 
        );
        return result;
      }

      // TODO this class design is so fucked; demands refactoring
      // if so, save to member
      if( MOLECULES.GetSize())
      {
        SetReactants( MOLECULES);
      }
      if( REFERENCE.GetSize())
      {
        SetReference( REFERENCE);
      }

      // isomorphisms between each molecule and its corresponding reactant substructure;
      // outer vector tracks molecules (reactants);
      // inner vector acts as a map between .rxn file atom ranks (vector indices)
      // and the atom ranks (0-indices) of of the user-provided molecules (values)
      storage::Vector< storage::Vector< size_t> > reactant_substr_isos;
      reactant_substr_isos.AllocateMemory( n_molecules);

      // atom graphs of each molecule
      storage::Vector< ConformationGraphConverter::t_AtomGraph> mol_atom_graphs;
      mol_atom_graphs.AllocateMemory( n_molecules);

      // reaction graph resolution
      ConformationGraphConverter converter
      (
        ConformationGraphConverter::e_ElementType,
        ConfigurationalBondTypeData::e_BondOrderOrAromatic
      );

      // initialize graphs for all molecules
      storage::Vector< graph::ConstGraph< size_t, size_t> > mol_graphs;
      mol_graphs.AllocateMemory( n_molecules);

      // a vector to make molecule access easier
      util::SiPtrVector< const FragmentComplete> mol_ptrs;
      mol_ptrs.AllocateMemory( n_molecules);

      // make a ReactionStructure, atom graph, and pointer for each molecule; check that the molecule matches its reactant 
      size_t n( 0);
      for
      (
        FragmentEnsemble::const_iterator itr_mol( MOLECULES.Begin()), itr_mol_end( MOLECULES.End());
        itr_mol != itr_mol_end;
        ++itr_mol, ++n
      )
      {
        mol_graphs.PushBack( converter( *itr_mol));
        mol_atom_graphs.PushBack( ConformationGraphConverter::CreateGraphWithAtoms( *itr_mol));
        mol_ptrs.PushBack( &( *itr_mol));

        // All disparate isomorphisms that satisfy the reactant substructure
        // TODO: see if we can increase the resolution of the substructure comparisons to include amide bonds
        storage::Vector< storage::Vector< size_t> > matched_isos( rxt_structs( n).GetMatchingSubstructures( *itr_mol));

        // bail if there are no matching substructures
        if( matched_isos.IsEmpty())
        {
          std::string desc( REACTION.GetDescription());
          desc.erase( std::remove( desc.begin(), desc.end(), '\n'), desc.end());
          BCL_MessageStd
          ( 
            "Could not execute reaction \"" + desc + "\"!  Molecule given as reactant " + util::Format()( n) +
            " did not match the specified substructure in the reaction"
          );
          return result;
        } 
        reactant_substr_isos.PushBack( matched_isos( 0));
      }

      //
      // determine which atoms in each provided molecule are reactive atoms 
      //

      // mapping between reaction mapped atoms (10th column after element type in .rxn)
      // and atom rank in .rxn file; also referred to as the generic atom mapping;
      // for clarity, 'rank' is used to indicate 0-indexed atom numbers
      const storage::Map< size_t, storage::Pair< size_t, size_t> > &atom_map
      (
        REACTION.GetReactiveAtomsAllReactants()
      );

      // translate the generic atom mapping to a mapping relevant to the given molecules
      m_ReactantsReactiveAtoms = storage::Vector< storage::Map< size_t, size_t> >( n_molecules);

      for
      (
        storage::Map< size_t, storage::Pair< size_t, size_t> >::const_iterator
        itr_map( atom_map.Begin()), itr_map_end( atom_map.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        const size_t &mol_no( itr_map->second.First()); // reactant number in .rxn
        const size_t &substr_atom( itr_map->second.Second()); // atom rank in .rxn
        size_t new_atom_no( reactant_substr_isos( mol_no)( substr_atom)); // atom rank in user-provided molecule

        // the translated map will use atom ranks from the user-provided reagent molecules as keys
        // and the .rxn file atom mapping number (10th column after element type in .rxn) as values
        m_ReactantsReactiveAtoms( mol_no)[ new_atom_no] = itr_map->first;
      }

      //
      // build up sets of atoms in each reactant that should be connected together on their own (not through reactive atoms)
      //

      // substituent information
      std::list< ReactionWorker::SubstituentInfo> mol_sub_infos;

      // get substituents from each molecule
      for( size_t mol_no( 0); mol_no < n_molecules; ++mol_no)
      {
        std::list< ReactionWorker::SubstituentInfo> subs
        (
          GetSubstituents
          (
            mol_graphs( mol_no),
            mol_atom_graphs( mol_no),
            m_ReactantsReactiveAtoms( mol_no),
            reactant_substr_isos( mol_no)
          )
        );
        mol_sub_infos.splice( mol_sub_infos.end(), subs);
      }

      // Add substituents to each product
      return AddSubstituentsToProducts( REACTION, mol_sub_infos);
    }

    //! @brief executes a reaction where a single molecule matches all reactants
    //! @param REACTION the reaction to use
    //! @param MOLECULE the molecule to use as a reactant
    //! @return an ensemble of products from the reaction
    FragmentEnsemble ReactionWorker::ExecuteIntramolecularReaction
    (
      const ReactionComplete &REACTION,
      const FragmentComplete &MOLECULE
    ) const
    {
      const storage::Vector< ReactionStructure> &rxt_structs( GetReactantStructures( REACTION));

      size_t n_reactants( REACTION.GetNumberReactants());

      FragmentEnsemble result;

      // Graph of the molecule
      ConformationGraphConverter converter( ConformationGraphConverter::e_ElementType, ConfigurationalBondTypeData::e_BondOrderOrAromatic); 
      graph::ConstGraph< size_t, size_t> mol_graph( converter( MOLECULE));
      ConformationGraphConverter::t_AtomGraph mol_atom_graph( ConformationGraphConverter::CreateGraphWithAtoms( MOLECULE));
      
      storage::Vector< storage::Vector< storage::Vector< size_t> > > isomorphisms( n_reactants);

      for( size_t r( 0); r < n_reactants; ++r)
      {
        storage::Vector< storage::Vector< size_t> > iso_matches( rxt_structs( r).GetMatchingSubstructures( MOLECULE));

        if( iso_matches.IsEmpty())
        {
          BCL_MessageStd
          ( 
            "Cannot perform intramolecular reaction: reactant " + util::Format()( r) + " was missing"
          );
          return result;
        }
        
        // Store the isomorphisms
        isomorphisms( r).Append( iso_matches);
      }

      storage::Vector< size_t> solution( n_reactants, size_t( 0));
      storage::Vector< size_t> chosen_atoms( MOLECULE.GetNumberAtoms(), size_t( 0));
      if( !ChooseSubgraphCombination( isomorphisms, 0, chosen_atoms, solution))
      {
        return result;
      }

      // Consolidate the isomorphisms
      storage::Vector< size_t> all_reactant_iso;
      all_reactant_iso.AllocateMemory( chosen_atoms.GetSize());
      storage::Map< size_t, size_t> reactive_atom_map;
      for( size_t r( 0); r < n_reactants; ++r)
      {
        storage::Vector< size_t> &cur_iso( isomorphisms( r)( solution( r)));
        size_t cur_iso_size( cur_iso.GetSize());

        // Append the isomorphism to the big iso
        all_reactant_iso.Append( cur_iso);

        // Figure out a mapping to reactive atoms
        storage::Map< size_t, size_t> reactive_atoms( REACTION.GetReactiveAtomsInReactant( r));
        
        // Determine how reactive atoms map to MOLECULE
        for
        (
          storage::Map< size_t, size_t>::const_iterator itr_atom( reactive_atoms.Begin()),
            itr_atom_end( reactive_atoms.End());
          itr_atom != itr_atom_end;
          ++itr_atom
        )
        {
          // Find where the reactant atom shows up in the isomorphism
          if( itr_atom->second > cur_iso_size)
          {
            BCL_MessageStd
            ( 
              "IntramolecularReaction: Something went wrong with atom mapping.  "
              "Could not find reactive atom " + util::Format()( itr_atom->first) + " in the input molecule"
            );
            return result;
          }
          // Set the discovered atom as this reactive atom
          reactive_atom_map[ cur_iso( itr_atom->second)] = itr_atom->first;
        }
      }
      
      // Get Substituents
      std::list< SubstituentInfo> subs
      ( 
        GetSubstituents( mol_graph, mol_atom_graph, reactive_atom_map, all_reactant_iso)
      );

      FragmentEnsemble new_ensemble( AddSubstituentsToProducts( REACTION, subs)); 
      return new_ensemble;
    }

    /*
    //! TODO: Implement this
    //! @brief executes a reaction wherein a subset of the reactants may be present
    //!        in a single provided molecule.
    //! @details tries to maximally match each molecule with as many reactants as possible
    //! @return an ensemble of product molecules
    FragmentEnsemble ReactionWorker::ExecutePartialIntramolecularReaction
    (
      const ReactionComplete &REACTION,
      const FragmentEnsemble &MOLECULES
    ) const
    {
      return FragmentEnsemble();
    }
    */

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ReactionWorker::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT number of indentations
    //! @return ostream which was written to
    std::ostream &ReactionWorker::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    /*
    //! @brief get substituents with
    //! @param MOL_GRAPH a graph of the molecule labeled by ElementType and BondOrderOrAromatic
    //! @param MOL_ATOM_GRAPH a graph of the molecule with atoms
    //! @param REACTIVE_SUBSTR isomorphism of the reactive substructure; anything in this substructure
    //!        that is not mapped will be removed
    //! @param MAPPED_ATOMS reaction mapped atoms (keys) corresponding to atoms in the given molecule
    BuildSubstituents
    (
      graph::ConstGraph< size_t, size_t> MOL_GRAPH,
      const ConformationGraphConverter::t_AtomGraph &MOL_ATOM_GRAPH,
      const storage::Vector< size_t> &REACTIVE_SUBSTR
      const storage::Map< size_t, size_t> &MAPPED_ATOMS
    )
    {
      graph::ConstGraph< size_t, size_t> &mol_graph( MOL_GRAPH);

      storage::Vector< bool> keep_atoms( MOL_GRAPH.GetSize(), true);
      
      // Discard reactive substructure atoms
      for
      (
        storage::Vector< size_t>::const_iterator itr_substr( REACTIVE_SUBSTR.Begin()), itr_substr_end( REACTIVE_SUBSTR.End());
        itr_substr != itr_substr_end;
        ++itr_substr
      )
      {
        keep_atoms( *itr_substr) = false;
      }

      // Keep mapped atoms
      for
      (
        storage::Map< size_t, size_t>::const_iterator itr_map( MAPPED_ATOMS.Begin()), itr_map_end( MAPPED_ATOMS.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        keep_atoms( itr_map->second) = true;
      }

      // remove bonds between kept and discarded atoms
      for( size_t i( 0), end_i( keep_atoms.GetSize()); i < end_i; ++i)
      {
        if( keep_atoms( i))
        {
          const storage::Vector< size_t> &neighbors( mol_graph.GetNeighborIndices( i));
          for
          (
            storage::Vector< size_t>::const_iterator itr_neighbor( neighbors.Begin()), itr_neighbor_end( neighbors.End());
            itr_neighbor != itr_neighbor_end;
            ++itr_neighbor
          )
          {
          }
        }
      }
    }

    */

    //! @brief gets which atoms are substituents of a molecule that should be preserved
    //! @param MOL_GRAPH the graph of the molecule of interest
    //! @param REACTIVE_ATOM_MAP map of which vertices (keys) are which reactive atoms (values)
    //! @param REACTANT_SUBSTRUCTURE vertices that are specified as parts of any reactant
    //! @return a list of substituent information structures 
    std::list< ReactionWorker::SubstituentInfo> ReactionWorker::GetSubstituents
    (
      const graph::ConstGraph< size_t, size_t> &MOL_GRAPH,
      const ConformationGraphConverter::t_AtomGraph &MOL_ATOM_GRAPH,
      const storage::Map< size_t, size_t> &REACTIVE_ATOM_MAP,
      const storage::Vector< size_t> &REACTANT_SUBSTRUCTURE
    ) const
    {
      // copy the graph since bonds will be added/removed
      graph::ConstGraph< size_t, size_t> mol_graph( MOL_GRAPH);

      // keep everything by default
      storage::Vector< size_t> keep_indices( mol_graph.GetSize(), size_t( 1));
      
      storage::Vector< storage::Pair< size_t, size_t> > reactive_atom_connections;
      storage::Vector< size_t> bond_types;

      // flag atoms that should be discarded
      for( size_t r( 0); r < REACTANT_SUBSTRUCTURE.GetSize(); ++r)
      {
        // discard atoms that are mapped to a reactant atom
        const size_t &reactant_atom_index( REACTANT_SUBSTRUCTURE( r));
        keep_indices( reactant_atom_index) = 0;
      }

      // isolate atoms from the reactant substructure;
      // save bonds from reactive atoms to non-substructure atoms
      for( size_t r( 0); r < REACTANT_SUBSTRUCTURE.GetSize(); ++r)
      {
        const size_t &reactant_atom_index( REACTANT_SUBSTRUCTURE( r));

        // determine if the atom is a reactive atom 
        size_t reactive_atom( 0); 
        bool store_bonds( REACTIVE_ATOM_MAP.Has( reactant_atom_index));
        if( store_bonds)
        {
          // get the reactive atom mapping number corresponding to .rxn map
          reactive_atom = REACTIVE_ATOM_MAP.GetValue( reactant_atom_index);
        }

        // adjacent atoms
        storage::Vector< size_t> neighbors( mol_graph.GetNeighborIndices( reactant_atom_index));

        // remove bonds to neighbors, i.e. isolate the atom
        for( size_t n( 0); n < neighbors.GetSize(); ++n)
        {
          const size_t &neighbor_atom( neighbors( n));

          // store information about the bond if this atom is a reactive atom and its neighbor was not specified
          // as part of the reactant substructure
          if( store_bonds && keep_indices( neighbor_atom))
          {
            bond_types.PushBack( MOL_ATOM_GRAPH.GetEdgeData( reactant_atom_index, neighbor_atom));
            reactive_atom_connections.PushBack( storage::Pair< size_t, size_t>( reactive_atom, neighbor_atom));
          }
          mol_graph.RemoveEdge( reactant_atom_index, neighbor_atom);
        }
      }

      //
      // collect connected components not present in the reactant substruct
      //

      storage::List< storage::Vector< size_t> > components( graph::Connectivity::GetComponents( mol_graph));

      // substituent info
      std::list< SubstituentInfo> sub_info;

      // discard any component that is part of the reactant structure (single atom with keep_indices == 0)
      for
      ( 
        storage::List< storage::Vector< size_t> >::iterator itr_comp( components.Begin()), 
          itr_comp_end( components.End());
        itr_comp != itr_comp_end;
        ++itr_comp
      )
      {
        // skip the component if it's a discarded atom  
        if( itr_comp->GetSize() == 1 && !keep_indices( ( *itr_comp)( 0)))
        {
          continue;
        }
        
        // if we get here then keep these components
        storage::Vector< size_t> &component( *itr_comp);
        
        // Information for this component/substituent
        SubstituentInfo component_info;
        storage::Map< size_t, storage::Vector< storage::Pair< size_t, size_t> > > &bond_map( component_info.m_Connections);

        // find bonds to reactive atoms in this component
        for( size_t i( 0); i < component.GetSize(); ++i)
        {
          // check the connection list for any bonds to this atom
          for( size_t c( 0); c < reactive_atom_connections.GetSize(); ++c)
          {
            const size_t &target_atom( reactive_atom_connections( c).Second());
            if( target_atom == component( i))
            {
              // there is a reactive atom connected to this target atom;
              // save the bond as part of this component's bond map
              const size_t &reactive_atom( reactive_atom_connections( c).First());
              bond_map[ reactive_atom].PushBack
              ( 
                storage::Pair< size_t, size_t>
                (
                  i,              // atom index after extracting this component substructure from the parent
                  bond_types( c)  // the bond type to the reactive atom
                )                
              ); 
            }
          }
        }
        if( !component_info.m_Connections.IsEmpty())
        {
          component_info.m_Atoms = ConformationGraphConverter::CreateAtomsFromGraph( MOL_ATOM_GRAPH.GetSubgraph( component));
          sub_info.push_back( component_info);
        }
      }
      return sub_info;
    }

    //! @brief adds substituents to products in a reaction
    //! @param REACTION the reaction
    //! @param SUBSTITUENTS the substituents to add
    //! @return an ensemble of products that have been substituted at the reactive atoms
    FragmentEnsemble ReactionWorker::AddSubstituentsToProducts
    (
      const ReactionComplete &REACTION,
      const std::list< SubstituentInfo> &SUBSTITUENTS
    ) const
    {
      FragmentEnsemble substituted_mols;

      // reference to products
      const storage::Vector< FragmentComplete> &products( REACTION.GetProducts());
      const storage::Vector< FragmentComplete> &reactants( REACTION.GetReactants());

      // Copy this so it can be updated
      std::list< ReactionWorker::SubstituentInfo> mol_sub_infos( SUBSTITUENTS);

      // do this for each product
      for( size_t p_no( 0); p_no < products.GetSize(); ++p_no)
      {

        // Get the reactive atoms in this product
        storage::Map< size_t, size_t> product_reactive_atoms( REACTION.GetReactiveAtomsInProduct( p_no));

        // a new atom vector derived from this product
        const AtomVector< AtomComplete> &product_atoms( products( p_no).GetAtomVector());

        // in the legacy code this atom vector was just connected to subs, standardized, and pushed back
        AtomVector< AtomComplete> parent_atoms( product_atoms); // do not remove hydrogen atoms until affter the product is assembled
        storage::Vector< size_t> mobile_atoms, connecting_atoms, secondary_mobile;
        std::string mobile_atoms_str;
        FragmentMapConformer cleaner
        (
          "",                                 // druglikeness type
          "",                                 // receptor MDL property
          "",                                 // receptor filename
          descriptor::CheminfoProperty(),     // affinity net
          false,                              // resolve protein-ligand clashes
          storage::Vector< float>(),          // b-factors
          false,                              // corina confs
          storage::Vector< size_t>(),         // moveable indices
          false,                              // best aligned conf
          m_CorrectGeometry,                  // fix geometry
          m_AdditionalAdjacentAtoms           // adjacent atoms
        );

        BCL_MessageVrb( "Product conformer arbitrary? " + util::Format()( m_ProductConformerArbitrary ? "true" : "false"));
        if( !m_ProductConformerArbitrary)
        {
          // give the product realistic 3D coordinates
          BCL_MessageVrb( "Cleaning product atoms");
          AtomVector< AtomComplete> temp_vec( cleaner.CleanAtoms( parent_atoms, "None", true, true));
          FragmentComplete clean_mol( temp_vec, "");
          BCL_MessageVrb( "Cleaning 3D coordinates of product");
          clean_mol = cleaner.Clean3DCoords( clean_mol);
          if( !clean_mol.GetSize())
          {
            continue;
          }

          // align the product
          BCL_MessageVrb( "Aligning product to reference reactant");
          bool alignment_success
          (
            AlignProductToReactants
            (
              clean_mol,
              product_reactive_atoms,
              m_Reference,
              m_ReactantsReactiveAtoms( m_ReferenceReactantIndex)
            )
          );

          // set the product atoms to be mobile
          parent_atoms = clean_mol.GetAtomVector();
          if( !alignment_success)
          {
            for
            (
                auto atom_itr( parent_atoms.Begin()), atom_itr_end( parent_atoms.End());
                atom_itr != atom_itr_end;
                ++atom_itr
            )
            {
              mobile_atoms.PushBack( parent_atoms.GetAtomIndex( *atom_itr));
              mobile_atoms_str.append( util::Format()( parent_atoms.GetAtomIndex( *atom_itr)) + " ");
            }
          }
        }
        // Iterate through available substituents and determine if they share any reactive atoms with this product
        // if so, merge them into the product and update the substituent
        size_t n_parent_atoms( parent_atoms.GetSize());
        storage::Vector< storage::Vector< size_t>> updated_substituent_indices( mol_sub_infos.size());
        size_t subs_i( 0);
        for
        (
            std::list< ReactionWorker::SubstituentInfo>::iterator
            itr_sub( mol_sub_infos.begin()), itr_sub_end( mol_sub_infos.end());
            itr_sub != itr_sub_end;
            ++itr_sub, ++subs_i
        )
        {
          // iterate through reactive atom bonds and determine if any need to be formed between
          // this substituent and the product
          storage::Vector< sdf::BondInfo> form_bonds;
          for
          (
              storage::Map< size_t, storage::Vector< storage::Pair< size_t, size_t> > >::const_iterator
              itr_conn( itr_sub->m_Connections.Begin()), itr_conn_end( itr_sub->m_Connections.End());
              itr_conn != itr_conn_end;
              ++itr_conn
          )
          {
            // component shares a reactive atom with the product; store the bond information
            if( product_reactive_atoms.Has( itr_conn->first))
            {
              const storage::Vector< storage::Pair< size_t, size_t> > &connections( itr_conn->second);
              for( size_t b_no( 0); b_no < connections.GetSize(); ++b_no)
              {
                // add directly connected atom to list of mobile atoms
                connecting_atoms.PushBack( connections( b_no).First() + n_parent_atoms);
                mobile_atoms_str.append( util::Format()( connections( b_no).First() + n_parent_atoms) + " ");

                form_bonds.PushBack
                ( 
                  sdf::BondInfo
                  (
                    product_reactive_atoms.GetValue( itr_conn->first), // product atom
                    connections( b_no).First() + n_parent_atoms,       // sub atom after appending to the parent
                    ConfigurationalBondType( connections( b_no).Second()) // bond type to form
                  )
                );
              }
            }
          }

          // If there are bonds to form, merge this molecule into the parent
          if( !form_bonds.IsEmpty())
          {
            size_t n_atoms( itr_sub->m_Atoms.GetSize());
            parent_atoms.AddAtomsWithConnectivity
            (
              itr_sub->m_Atoms,
              form_bonds
            );

            // get new atom indices for this substituent
            for( size_t i( 0); i < itr_sub->m_Atoms.GetSize(); ++i)
            {
              updated_substituent_indices( subs_i).PushBack(  i + n_parent_atoms);
            }

            // update the substituent to reflect the merge 
            itr_sub->m_Atoms = parent_atoms;
            for
            (
                storage::Map< size_t, storage::Vector< storage::Pair< size_t, size_t> > >::iterator itr_conn( itr_sub->m_Connections.Begin()),
                itr_conn_end( itr_sub->m_Connections.End());
                itr_conn != itr_conn_end;
                ++itr_conn
            )
            {
              storage::Vector< storage::Pair< size_t, size_t> > &connections( itr_conn->second);
              for( size_t b_no( 0); b_no < connections.GetSize(); ++b_no)
              {
                connections( b_no).First() += n_parent_atoms;
              }
            }
            n_parent_atoms += n_atoms; 
          }
        }

        // standardize atoms and get a new fragment
        AtomsCompleteStandardizer std
        (
          parent_atoms,               // modify on parent atom vector
          "",                         // unimportant
          m_ProductConformerArbitrary // if arbitrary then force recalculation; otherwise do not force
        );
        FragmentComplete new_frag( parent_atoms, "");

        // more extensive and slower if we care about the pose
        if( !m_ProductConformerArbitrary)
        {
          // first fix bond lengths
          MutateBondLengths bond_length_mutate( 5, 0.01, 0.01, mobile_atoms);
          bond_length_mutate.Initialize( new_frag);

          // track which atoms were the ones with bad bond lengths
          const storage::Vector< double> &bad_bonds( bond_length_mutate.GetInfo().m_BondLengthsDiff);
          storage::Set< size_t> atoms_in_bad_bonds;
          for
          (
              size_t i( 0); i < bad_bonds.GetSize(); ++i
          )
          {
            if( bad_bonds( i))
            {
              atoms_in_bad_bonds.InsertElement( bond_length_mutate.GetInfo().m_BondInfo( i).GetAtomIndexLow());
              atoms_in_bad_bonds.InsertElement( bond_length_mutate.GetInfo().m_BondInfo( i).GetAtomIndexHigh());
            }
          }

          // run bond length mutate
          BCL_MessageVrb( "Fix product bond lengths");
          auto fixed_frag( bond_length_mutate( new_frag));
          if( fixed_frag.GetArgument().IsDefined())
          {
            // contains at least one bad bond length
            new_frag = FragmentComplete( *( fixed_frag.GetArgument()));
          }
          new_frag.SaturateWithH();

          // then add in extra atoms directly connected to product base fragment
          // for clash resolution
          storage::Set< size_t> nbr_indices;

          // loop over the atoms historically involved in bad bonds
          storage::Vector< size_t> bad_atoms( atoms_in_bad_bonds.Begin(), atoms_in_bad_bonds.End());
          for( size_t bi( 0), bsz( bad_atoms.GetSize()); bi < bsz; ++bi)
          {
            // loop over each of the substituent fragments
            for( size_t si( 0), ssz( updated_substituent_indices.GetSize()); si < ssz; ++si)
            {
              // check if bad bond atoms are in substituent
              if( updated_substituent_indices( si).Find( bad_atoms( bi)) < updated_substituent_indices( si).GetSize())
              {
                secondary_mobile.InsertElements( secondary_mobile.GetSize(), updated_substituent_indices( si));
              }
            }
          }

          // add neighbor atoms to atoms historically involved in bad bonds
          nbr_indices.InsertElements( secondary_mobile.Begin(), secondary_mobile.End());
          for( size_t i( 0); i < secondary_mobile.GetSize(); ++i)
          {
            // get indices of bonded atoms and add to set
            for
            (
                auto bond_itr( new_frag.GetAtomVector()( secondary_mobile( i)).GetBonds().Begin()),
                bond_itr_end( new_frag.GetAtomVector()( secondary_mobile( i)).GetBonds().End());
                bond_itr != bond_itr_end;
                ++bond_itr
            )
            {
              nbr_indices.InsertElement( new_frag.GetAtomVector().GetAtomIndex( bond_itr->GetTargetAtom()));
            }
          }
          secondary_mobile = storage::Vector< size_t>( nbr_indices.Begin(), nbr_indices.End());
          mobile_atoms.InsertElements( mobile_atoms.GetSize(), secondary_mobile);

          // adding hydrogen atoms now will not change indices since they are added to end
          new_frag.RemoveH(); // remove for standardization
          AtomVector< AtomComplete> mol_atom_v( new_frag.GetAtomVector());
          AtomsCompleteStandardizer standardizer( mol_atom_v, "", true);
          standardizer.SetConjugationOfBondTypes( mol_atom_v);
          BondIsometryHandler::AddIsometryInformation( mol_atom_v, true);
          StereocentersHandler::AddChiralityFromConformation( mol_atom_v);
          new_frag = FragmentComplete( mol_atom_v, "");
          new_frag.SaturateWithH(); // add back before clash resolver

          // then resolve clashes
          MutateClashResolver clash_resolver( 100, mobile_atoms);
          static RotamerLibraryFile rotamer_library;
          util::ShPtr< SearchFragmentLibraryFromTree> rotlib_searcher;
          rotlib_searcher = util::ShPtr< SearchFragmentLibraryFromTree>( new SearchFragmentLibraryFromTree( rotamer_library));
          util::ShPtrVector< BondAngleAssignment> bond_angle_assignments
          (
            rotlib_searcher->GetBondAngleAssignments( new_frag, true)
          );

          util::ShPtrList< MutateDihedralsInterface> mutate_bond_angles_refinement;
          MutateChirality chirality_changer( new_frag);
          for( auto itr_bas( bond_angle_assignments.Begin()), itr_bas_end( bond_angle_assignments.End()); itr_bas != itr_bas_end; ++itr_bas)
          {
            mutate_bond_angles_refinement.PushBack
            (
              util::ShPtr< MutateBondAngles>
              (
                new MutateBondAngles
                (
                  *itr_bas,
                  new_frag,
                  chirality_changer,
                  false,
                  true,
                  true
                )
              )
            );
          }
          util::ShPtr< AtomClashScore> clash_score_sp( new AtomClashScore( true));
          clash_resolver.Setup
          (
            new_frag,
            util::ShPtrVector< MutateBondAngles>( mutate_bond_angles_refinement.Begin(), mutate_bond_angles_refinement.End()),
            clash_score_sp
          );
          BCL_MessageVrb( "Resolve product clashes");
          auto fixed_mol( clash_resolver( new_frag));

          // if we have to resolve clashes, then afterward also check the dihedrals
          if( fixed_mol.GetArgument().IsDefined())
          {
            // restandardize again
            new_frag = *( fixed_mol.GetArgument());
            new_frag.RemoveH(); // we do not reference atom indices explicitly anymore and this will help with standardization and iso
            AtomVector< AtomComplete> v( new_frag.GetAtomVector());
            AtomsCompleteStandardizer std_v( v, "", true);
            std_v.SetConjugationOfBondTypes( v);
            BondIsometryHandler::AddIsometryInformation( v, true);
            StereocentersHandler::AddChiralityFromConformation( v);
            new_frag = FragmentComplete( v, "");

            // sample by parts just the mobile atom indices
            BCL_MessageVrb( "Sample product conformations");
            static SampleConformations sample_confs_single
            (
              rotamer_library,        // rotamer library file
              "",                     // conformation comparer type
              0.0,                    // conformational comparer tolerance
              100,                    // number of conformations
              2000,                   // number of iterations
              false,                  // no change chirality
              0.0,                    // random dihedral change weight
              true,                   // generate 3d
              0.1,                    // clash tolerance
              true                    // clustering
            );

            // map atoms of interest to new mobile
            ConformationGraphConverter graph_maker
            (
              ConformationGraphConverter::AtomComparisonType::e_ElementType,
              ConfigurationalBondTypeData::Data::e_BondOrderAmideOrAromaticWithRingness,
              false
            );

            // set up the isomorphisms with pointers to the graphs
            graph::CommonSubgraphIsomorphism< size_t, size_t> isomorphism( graph::CommonSubgraphIsomorphismBase::SolutionType::e_Connected);
            isomorphism.SetGraphA( graph_maker( new_frag));
            isomorphism.SetGraphB( graph_maker( m_Reference));

            // get the isomorphism - allocate array based on maximum possible atom matches without considering connectivity
            isomorphism.FindIsomorphism( isomorphism.EstimateUpperBounds());

            // get graph of the new mol
            graph::Subgraph< size_t, size_t> subgraph
            (
              util::OwnPtr< const graph::ConstGraph< size_t, size_t> >( &isomorphism.GetGraphA(), false),
              isomorphism.GetIsomorphism().GetKeysAsVector()
            );

            // allow sampling on all atoms not in the reference substructure
            auto complement_subgraph( subgraph.GetComplement());
            mobile_atoms = complement_subgraph.GetVertexIndices();

            // also allow sampling on the bad geometry atoms
            if( m_CorrectGeometry)
            {
              storage::Vector< size_t> bad_geo_atoms( new_frag.GetAtomsWithBadGeometry());
              mobile_atoms.InsertElements( mobile_atoms.GetSize(), bad_geo_atoms);
              std::string mobile_atoms_str;
              for
              (
                  auto itr( bad_geo_atoms.Begin()), itr_end( bad_geo_atoms.End());
                  itr != itr_end;
                  ++itr
              )
              {
                mobile_atoms_str.append( util::Format()( *itr) + " ");
              }
              BCL_MessageStd( "The following atoms with bad geometry will be explicitly added for conformer sampling: " + util::Format()( mobile_atoms_str));
            }

            // add atoms adjacent to mobile atoms
            if( m_AdditionalAdjacentAtoms)
            {
              storage::Set< size_t> adjacent_indices( cleaner.MapSubgraphAdjacentAtoms( complement_subgraph, m_AdditionalAdjacentAtoms));
              storage::Vector< size_t> adjacent_indices_v( adjacent_indices.Begin(), adjacent_indices.End());
              mobile_atoms.InsertElements( mobile_atoms.GetSize(), adjacent_indices_v);
              std::string mobile_atoms_str;
              for
              (
                  auto itr( adjacent_indices_v.Begin()), itr_end( adjacent_indices_v.End());
                  itr != itr_end;
                  ++itr
              )
              {
                mobile_atoms_str.append( util::Format()( *itr) + " ");
              }
              BCL_MessageStd( "The following extended adjacent atoms will be explicitly added for conformer sampling: " + util::Format()( mobile_atoms_str));
            }

            // add all rings not in subgraph
            if( m_CorrectNonReferenceRingGeometry)
            {
              ConformationGraphConverter::t_AtomGraph molecule_graph( ConformationGraphConverter::CreateGraphWithAtoms( new_frag));
              FragmentSplitRings splitter( true, size_t( 2));
              storage::List< storage::Vector< size_t> > ring_components( splitter.GetComponentVertices( new_frag, molecule_graph));
              storage::Vector< storage::Vector< size_t>> ring_components_vec( ring_components.Begin(), ring_components.End());
              storage::Set< size_t> important_ring_atoms( cleaner.MapSubgraphRingAtoms( complement_subgraph, ring_components));
              storage::Vector< size_t> important_ring_atoms_v( important_ring_atoms.Begin(), important_ring_atoms.End());
              mobile_atoms.InsertElements( mobile_atoms.GetSize(), important_ring_atoms_v);
              std::string mobile_atoms_str;
              for
              (
                  auto itr( important_ring_atoms_v.Begin()), itr_end( important_ring_atoms_v.End());
                  itr != itr_end;
                  ++itr
              )
              {
                mobile_atoms_str.append( util::Format()( *itr) + " ");
              }
              BCL_MessageStd( "The following ring atoms will be explicitly added for conformer sampling: " + util::Format()( mobile_atoms_str));
            }

            // sample conformations
            storage::Set< size_t> mobile_atoms_set( mobile_atoms.Begin(), mobile_atoms.End());
            sample_confs_single.SetSampleByPartsIndices( storage::Vector< size_t>( mobile_atoms_set.Begin(), mobile_atoms_set.End()));
            new_frag.SaturateWithH(); // add back for our final time before conformer sampling
            auto confs( sample_confs_single( new_frag));

            std::string mobile_atoms_str;
            for
            (
                auto itr( mobile_atoms_set.Begin()), itr_end( mobile_atoms_set.End());
                itr != itr_end;
                ++itr
            )
            {
              mobile_atoms_str.append( util::Format()( *itr) + " ");
            }
            BCL_MessageStd( "The following atoms will be added to SampleByParts indices: " + util::Format()( mobile_atoms_str));

            // filter the ensemble to get the best conformer with the non-mobile atoms at starting position
            // TODO add construction/get-set/seralization control over these options
            math::Comparisons< float>::Comparison compare_type( "less_equal");
            auto filtered_confs
            (
              sample_confs_single.FilterConformerEnsemble
              (
                confs.First(),
                new_frag,
                "SymmetryRealSpaceRMSD",
                compare_type,
                1.0,
                isomorphism.GetIsomorphism().GetKeysAsVector()
              )
            );

            // save best conf score (should be 0-index from filter)
            if( filtered_confs.GetSize())
            {
              // assign properties to output molecule
              new_frag = filtered_confs.GetMolecules().FirstElement();
            }
            new_frag.StoreProperty( "Reaction", REACTION.GetDescription());
          }

          // re-align to initial reagent one last time
          AlignProductToReactants
          (
            new_frag,
            storage::Map< size_t, size_t>(),
            m_Reference,
            storage::Map< size_t, size_t>()
          );
        }
        // legacy: similar to Alex's original implementation
        else
        {
          std.SetConjugationOfBondTypes( parent_atoms);
          BondIsometryHandler::AddIsometryInformation( parent_atoms, true);
          StereocentersHandler::AddChiralityFromConformation( parent_atoms);
          new_frag = FragmentComplete( parent_atoms, "");
        }
        substituted_mols.PushBack( new_frag);
      }
      return substituted_mols;
    }

    //! @brief choose a combination of numbers from a set, given constrains on which numbers 
    //!        can be chosen for each solution slot
    //! @param ALLOWED_VALUES values that are allowed at each solution slot
    //! @param SOLN_NO which solution slot should be chosen 
    //! @param CHOSEN_VALUES whether a value has been used already (index = value, 0 or 1)
    //! @param SOLUTION the combination of values that satisfies the constraints
    //! @return true if a combination of values could be found for this and subsequent slots 
    bool ReactionWorker::ChooseValueCombination
    (
      const storage::Vector< storage::Vector< size_t> > &ALLOWED_VALUES,
      const size_t &SOLN_NO,
      storage::Vector< size_t> &CHOSEN_VALUES,
      storage::Vector< size_t> &SOLUTION
    ) const
    {
      // if we have run out of slots, return true; nothing to do
      if( SOLN_NO >= SOLUTION.GetSize())
      {
        return true;
      }

      const storage::Vector< size_t> &allowed_values( ALLOWED_VALUES( SOLN_NO));
      size_t n_allowed_values( allowed_values.GetSize());

      // iterate through available values and pick each one until we find a satisfactory choice
      for( size_t i( 0); i < n_allowed_values; ++i)
      {
        const size_t &chosen_value( allowed_values( i));

        // if this chosen_value hasn't been picked before, recurse to the next slot to see
        // if it can be chosen from the remaining values 
        if( !CHOSEN_VALUES( chosen_value))
        {
          CHOSEN_VALUES( chosen_value) = 1;

          // values could be picked that satisfied the constraints of all subsequent slots
          // use this value for this  
          if( ChooseValueCombination( ALLOWED_VALUES, SOLN_NO + 1, CHOSEN_VALUES, SOLUTION))
          {
            SOLUTION( SOLN_NO) = chosen_value;
            return true;
          }

          // satisfactory solution wasn't found, so clear this chosen_value
          CHOSEN_VALUES( chosen_value) = 0;
        }
      }

      // Nothing could be found
      return false;
    }

    //! @brief choose a combination of subgraphs (of a supergraph) which do not overlap, given 
    //!        constraints on which subgraphs can be chosen for each solution slot (recursive algorithm) 
    //! @param ALLOWED_SUBGRAPHS subgraphs (i.e. lists of indices) that are allowed for each solution slot
    //! @param SOLN_NO which solution slot should be chosen during the current iteration
    //! @param COVERED_VERTICES which vertices are already covered by chosen subgraphs
    //! @param SOLUTION which subgraph should be used for each solution slot (used to return info) 
    //! @return true if non-overlapping sugraphs for current and subsequent solution slots could be found
    bool ReactionWorker::ChooseSubgraphCombination
    (
      const storage::Vector< storage::Vector< storage::Vector< size_t> > > &ALLOWED_SUBGRAPHS,
      const size_t &SOLN_NO,
      storage::Vector< size_t> &COVERED_VERTICES,
      storage::Vector< size_t> &SOLUTION
    ) const
    {
      // the number of slots (i.e. solution slots) that must be picked
      size_t n_slots( ALLOWED_SUBGRAPHS.GetSize());

      // Made it to the end of the recursion, return true
      if( SOLN_NO >= n_slots)
      {
        return true;
      }

      // Iterate through allowed subgraphs for this slot
      const storage::Vector< storage::Vector< size_t> > &subgraphs( ALLOWED_SUBGRAPHS( SOLN_NO));
      size_t subgraphs_size( subgraphs.GetSize());
      for( size_t i( 0); i < subgraphs_size; ++i)
      {
        const storage::Vector< size_t> &cur_subgraph( subgraphs( i));
        size_t cur_subgraph_size( cur_subgraph.GetSize());

        // check if vertices in this subgraph have been covered already
        bool use_subgraph( true);
        for( size_t j( 0); j < cur_subgraph_size; ++j)
        {
          // If any of the atom indices are used, skip this subgraph
          if( COVERED_VERTICES( cur_subgraph( j)) == 1)
          {
            use_subgraph = false;
            break;
          }
        }

        // Shouldn't use the subgraph
        if( !use_subgraph)
        {
          continue;
        }
        
        // mark the vertices in this subgraph as covered 
        for( size_t j( 0); j < cur_subgraph_size; ++j)
        {
          COVERED_VERTICES( cur_subgraph( j)) = 1;
        }

        // recurse to the next solution slot; checks if picking this subgraph will
        // allow subsequent slots to be picked 
        if( ChooseSubgraphCombination( ALLOWED_SUBGRAPHS, SOLN_NO + 1, COVERED_VERTICES, SOLUTION))
        {
          // A solution could be found, so use this subgraph for this slot
          SOLUTION( SOLN_NO) = i;
          return true;
        }

        // picking this subgraph would prevent subsequent satisfying subsequent slots
        // unmark the vertices and try a different subgraph 
        for( size_t j( 0); j < cur_subgraph_size; ++j)
        {
          COVERED_VERTICES( cur_subgraph( j)) = 0;
        }
      }

      // solution couldn't be found
      return false;
    }

    //! @brief caches the reactants and products of a reaction if necessary
    //! @param REACTION the reaction to cache
    void ReactionWorker::Cache( const ReactionComplete &REACTION) const
    {
      util::SiPtr< const ReactionComplete> rxn_ptr( &REACTION);
      if( !m_ReactantCache.Has( rxn_ptr))
      {
        storage::Vector< FragmentComplete> reactants( REACTION.GetReactants());
        storage::Vector< ReactionStructure> &reactant_structs( m_ReactantCache[ rxn_ptr]);
        reactant_structs.AllocateMemory( reactants.GetSize());
        for
        (
          storage::Vector< FragmentComplete>::const_iterator itr_mol( reactants.Begin()), itr_mol_end( reactants.End());
          itr_mol != itr_mol_end;
          ++itr_mol
        )
        {
          // add a default RXNProperty if one does not exist already
          m_ReactantCache[ rxn_ptr].PushBack( ReactionStructure( *itr_mol));
        }
        
        storage::Vector< FragmentComplete> products( REACTION.GetProducts());
        storage::Vector< ReactionStructure> &product_structs( m_ProductCache[ rxn_ptr]);
        product_structs.AllocateMemory( products.GetSize());
        for
        (
          storage::Vector< FragmentComplete>::const_iterator itr_mol( products.Begin()), itr_mol_end( products.End());
          itr_mol != itr_mol_end;
          ++itr_mol
        )
        {
          product_structs.PushBack( ReactionStructure( *itr_mol));
        }
      }
    }

    //! @brief gets the reactant ReactionStructures from the reaction; retrieves from cache if available, caches first if not
    const storage::Vector< ReactionStructure> &ReactionWorker::GetReactantStructures
    ( 
      const ReactionComplete &REACTION
    ) const
    {
      Cache( REACTION);
      util::SiPtr< const ReactionComplete> rxn_ptr( &REACTION);
      return m_ReactantCache.GetValue( rxn_ptr);
    }

    //! @brief gets the product ReactionStructures from the reaction; retrieves from cache if available, caches first if not
    const storage::Vector< ReactionStructure> &ReactionWorker::GetProductStructures
    ( 
      const ReactionComplete &REACTION
    ) const
    {
      Cache( REACTION);
      util::SiPtr< const ReactionComplete> rxn_ptr( &REACTION);
      return m_ProductCache.GetValue( rxn_ptr);
    }

    // create a single instance of this class
    const util::SiPtr< const util::ObjectInterface> ReactionWorker::s_Instance
    (
      GetObjectInstances().AddInstance( new ReactionWorker())
    );

  } // namespace chemistry
} // namespace bcl
