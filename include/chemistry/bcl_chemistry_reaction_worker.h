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

#ifndef BCL_CHEMISTRY_REACTION_WORKER_H_
#define BCL_CHEMISTRY_REACTION_WORKER_H_

// include the namespace header
#include "bcl_chemistry.h"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_reaction_complete.h"
#include "bcl_chemistry_reaction_structure.h"
#include "graph/bcl_graph_const_graph.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_sh_ptr_vector.h"
#include "util/bcl_util_si_ptr.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ReactionWorker
    //! @brief this class uses reactions to manipulate molecules.  this includes executing reactions and  
    //!
    //! @see @link example_chemistry_reaction_worker.cpp @endlink
    //! @author geanesar, brownbp1
    //! @date Jan 28, 2015
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ReactionWorker :
      public util::ObjectInterface
    {
    protected:

    ////////////////////
    // helper classes //
    ////////////////////
      
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class SubstituentInfo
      //! @brief a simple container class used within ReactionWorker.  It keeps track of which atoms are attached to 
      //!        reactive atoms from a reaction
      //!
      //! @author geanesar
      //! @date Jan 28, 2015
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      struct SubstituentInfo 
      {
        //! atoms in this component
        AtomVector< AtomComplete> m_Atoms;

        //! the atom types of the mapped atoms in the molecule this substituent came from
        storage::Map< size_t, AtomType> m_MappedAtomTypes;

        //! a map from reactive atoms (keys) to pairs of connecting atoms (values)
        //! stored as atom index (first) and bond type (second)
        storage::Map< size_t, storage::Vector< storage::Pair< size_t, size_t> > > m_Connections;

        //! @brief default constructor
        SubstituentInfo() :
          m_Atoms(),
          m_Connections()
        {
        }
      };

    //////////
    // data //
    //////////

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! a cache of reactants from the reactions
      mutable storage::Map< util::SiPtr< const ReactionComplete>, storage::Vector< ReactionStructure> > m_ReactantCache;

      //! a cache of products from the reactions
      mutable storage::Map< util::SiPtr< const ReactionComplete>, storage::Vector< ReactionStructure> > m_ProductCache;

      //! the passed reagents
      mutable FragmentEnsemble m_Reactants;

      //! the target fragment used to initialize the reaction
      mutable FragmentComplete m_Reference;

      //! the reference target fragment reactant index
      mutable size_t m_ReferenceReactantIndex;

      //! maps relating reactants atom indices to reactive atoms
      mutable storage::Vector< storage::Map< size_t, size_t> > m_ReactantsReactiveAtoms;

      //! true if the final product(s) conformation does not need to be constrained in any way (i.e. is of arbitrary 3D coordinates)
      bool m_ProductConformerArbitrary;

      //! if 3D conformer matters, fix atoms with bad geometry even if they are in reference structure
      bool m_CorrectGeometry;

      //! if 3D conformer matters, add all ring atoms from non-reference scaffolds to mobile selection
      bool m_CorrectNonReferenceRingGeometry;

      //! if 3D conformer matters, fix atoms this many bonds out from any other mobile atom
      size_t m_AdditionalAdjacentAtoms;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor
      //! @param ATOM_COMPARISON how to compare atoms in structures
      //! @param BOND_COMPARISON how to compare bonds in structures
      ReactionWorker
      (
        const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON = 
          ConformationGraphConverter::e_ElementType,
        const ConfigurationalBondTypeData::Data &BOND_COMPARISON = 
          ConfigurationalBondTypeData::e_BondOrderOrAromatic
      );

      //! @brief copy constructor; does not copy cache
      ReactionWorker( const ReactionWorker &OTHER);

      //! @brief clone constructor
      //! @return a pointer to a copy of this class
      ReactionWorker *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the reactant cache
      const storage::Map< util::SiPtr< const ReactionComplete>, storage::Vector< ReactionStructure> >
      &GetReactantCache() const;

      //! @brief return the product cache
      const storage::Map< util::SiPtr< const ReactionComplete>, storage::Vector< ReactionStructure> >
      &GetProductCache() const;

      //! @brief get the passed reagents
      const FragmentEnsemble &GetReactants() const;

      //! @brief get the target fragment used to initialize the reaction
      const FragmentComplete &GetReference() const;

      //! @brief get the reference target fragment reactant index
      const size_t &GetReferenceReactantIndex() const;

      //! @brief get the maps relating reactants atom indices to reactive atoms
      const storage::Vector< storage::Map< size_t, size_t> > &GetReactantsReactiveAtoms() const;

      //! @brief return whether or not an arbitrary 3D conformer will be generated for final product(s)
      const bool &GetIsProductConformerArbitrary() const;

      //! @brief return whether or not an arbitrary 3D conformer will be generated for final product(s)
      const bool &GetIsCorrectGeometry() const;

      //! @brief return whether or not an arbitrary 3D conformer will be generated for final product(s)
      const bool &GetIsCorrectNonReferenceRingGeometry() const;

      //! @brief return whether or not an arbitrary 3D conformer will be generated for final product(s)
      const size_t &GetAdditionalAdjacentAtoms() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief assignment operator
      ReactionWorker &operator =( const ReactionWorker &OTHER);

    ////////////////
    // operations //
    ////////////////

      //! @brief searches reactions for those that MOLECULE could satisfy as a reactant
      //! @param MOLECULE the query molecule
      //! @param REACTION the reaction to search
      //! @param REACTANT_NUMS indices of reactants of interest (empty is all reactants)
      //! @return a vector of the reactant indices that match MOLECULE
      storage::Vector< size_t> MatchesReactants
      (
        const ConformationInterface &MOLECULE,
        const ReactionComplete &REACTION,
        const storage::Vector< size_t> &REACTANT_NUMS = storage::Vector< size_t>()
      ) const;

      //! @brief searches reactions for those that MOLECULE could satisfy as a product
      //! @param MOLECULE the query molecule
      //! @param REACTION the reaction to search
      //! @param PRODUCT_NUMS indices of products of interest (empty is all products)
      //! @return a vector of the product indices that match MOLECULE
      storage::Vector< size_t> MatchesProducts
      (
        const ConformationInterface &MOLECULE,
        const ReactionComplete &REACTION,
        const storage::Vector< size_t> &PRODUCT_NUMS = storage::Vector< size_t>()
      ) const;

      //! @brief searches reactions for those that MOLECULE could satisfy as a reactant
      //! @param MOLECULE the query molecule
      //! @param REACTION the reaction to search
      //! @param INDEX the reactant index against which the molecule is queried
      //! @param REACTANT_NUMS indices of reactants of interest (empty is all reactants)
      //! @return true if the MOLECULE matches the reactant at index INDEX
      bool MatchesReactantIndex
      (
        const ConformationInterface &MOLECULE,
        const ReactionComplete &REACTION,
        const size_t &INDEX,
        const storage::Vector< size_t> &REACTANT_NUMS = storage::Vector< size_t>()
      ) const;

      //! @brief tests if a molecule would be reactive with other molecules of itself under a given reaction; does not count for intramolecular reactions
      //! @param REACTION the reaction to inspect
      //! @param MOLECULE the molecule in question
      //! @return true if the molecule would react with itself, false otherwise
      bool IsSelfReactive
      (
        const ReactionComplete &REACTION,
        const FragmentComplete &MOLECULE
      ) const;
    
      //! @brief Executes a reaction with provided reactants, does not depend on the order reactants are given in
      //! @param REACTION the reaction to use
      //! @param MOLECULES molecules used to do the reaction
      //! @return an ensemble of products
      FragmentEnsemble React
      (
        const ReactionComplete &REACTION,
        const FragmentEnsemble &MOLECULES
      ) const;

      //! @brief executes a reaction given a set of molecules which must correspond exactly (in order) to reactants
      //! @param REACTION the reaction to execute
      //! @param MOLECULES the molecules that must match the reactants of REACTION in the same order
      //! @param REFERENCE the target reagent used to initialize the reaction
      //! @param REFERENCE_REACTANT_INDEX the index of the reference reagent within the reaction
      //! @return an ensemble of product molecules
      FragmentEnsemble ExecuteReaction
      (
        const ReactionComplete &REACTION,
        const FragmentEnsemble &MOLECULES,
        const FragmentComplete &REFERENCE = FragmentComplete(),
        const size_t &REFERENCE_REACTANT_INDEX = size_t()
      ) const;
      
      //! @brief executes a reaction where a single molecule matches all reactants
      //! @param REACTION the reaction to use
      //! @param MOLECULE the molecule to use as a reactant
      //! @return an ensemble of products from the reaction
      FragmentEnsemble ExecuteIntramolecularReaction
      (
        const ReactionComplete &REACTION,
        const FragmentComplete &MOLECULE
      ) const;

      /*
      //! TODO: Implement this
      //! @brief executes a reaction wherein a subset of the reactants may be present
      //!        in a single provided molecule.
      //! @details tries to maximally match each molecule with as many reactants as possible
      //! @return an ensemble of product molecules
      FragmentEnsemble ExecutePartialIntramolecularReaction
      (
        const ReactionComplete &REACTION,
        const FragmentEnsemble &MOLECULES
      ) const;
      */

    private:

      //! @brief set the reagents
      void SetReactants( const FragmentEnsemble &REACTANTS) const;

      //! @brief set the target fragment used to initialize the reaction
      void SetReference( const FragmentComplete &REFERENCE) const;

      //! @brief set the reference target fragment reactant index
      void SetReferenceReactantIndex( const size_t REFERENCE_REACTANT_INDEX) const;

      //! @brief set the maps relating reactants atom indices to reactive atoms
      void SetReactantsReactiveAtoms
      (
        const storage::Vector< storage::Map< size_t, size_t> > &REACTANTS_REACTIVE_ATOMS
      ) const;

    public:

      //! @brief set whether or not an arbitrary 3D conformer will be generated for final product(s)
      void SetProductConformerArbitrary( const bool PRODUCT_CONFORMER_ARBITRARY);

      //! @brief set whether or not an arbitrary 3D conformer will be generated for final product(s)
      void SetCorrectGeometry( const bool CORRECT_GEOMETRY);

      //! @brief set whether or not an arbitrary 3D conformer will be generated for final product(s)
      void SetCorrectNonReferenceRingGeometry( const bool CORRECT_NON_REF_RING_GEOMETRY);

      //! @brief set whether or not an arbitrary 3D conformer will be generated for final product(s)
      void SetAdditionalAdjacentAtoms( const size_t N_ADDITIONAL_ADJACENT_ATOMS);

    private:

      //! @brief aligns products to mapped reactive atoms
      //! @param PRODUCT the product to be aligned to reactants
      //! @param PRODUCT_REACTIVE_ATOMS map of reactive atoms to atom indices
      //! @param REACTANTS the reactants to which the product is aligned
      //! @param REACTANTS_REACTIVE_ATOMS maps relating atom indices to reactive atoms
      //! @return true if the alignment succeeds
      bool AlignProductToReactants
      (
        FragmentComplete &PRODUCT,
        const storage::Map< size_t, size_t> &PRODUCT_REACTIVE_ATOMS,
        const FragmentComplete &REFERENCE,
        const storage::Map< size_t, size_t> &REACTANTS_REACTIVE_ATOMS
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream
      //! @param INDENT number of indentations
      //! @return ostream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief gets which atoms are substituents of a molecule that should be preserved
      //! @param MOL_GRAPH the graph of the molecule of interest
      //! @param REACTIVE_ATOM_MAP map of which graph vertices (keys) are which reactive atoms (values)
      //! @param SPECIFIED_IN_REACTANT vertices that are specified as parts of any reactant
      //! @return a list of connected components of the graph that should be used
      std::list< SubstituentInfo> GetSubstituents
      (
        const graph::ConstGraph< size_t, size_t> &MOL_GRAPH,
        const ConformationGraphConverter::t_AtomGraph &MOL_ATOM_GRAPH,
        const storage::Map< size_t, size_t> &REACTIVE_ATOM_MAP,
        const storage::Vector< size_t> &SPECIFIED_IN_REACTANT
      ) const;

      //! @brief adds substituents to products in a reaction
      //! @param REACTION the reaction
      //! @param SUBSTITUENTS the substituents to add
      //! @return an ensemble of products that have been substituted at the reactive atoms
      FragmentEnsemble AddSubstituentsToProducts
      (
        const ReactionComplete &REACTION,
        const std::list< SubstituentInfo> &SUBSTITUENTS
      ) const;
    
    protected:

      //! @brief choose a combination of numbers from a set, given constrains on which numbers 
      //!        can be chosen for each solution slot
      //! @param ALLOWED_VALUES values that are allowed at each solution slot
      //! @param SOLN_NO which solution slot should be chosen 
      //! @param CHOSEN_VALUES values that have already been selected
      //! @param SOLUTION the combination of values that satisfies the constraints
      //! @return true if a combination of values could be found for this and subsequent slots 
      bool ChooseValueCombination
      (
        const storage::Vector< storage::Vector< size_t> > &ALLOWED_VALUES,
        const size_t &SOLN_NO,
        storage::Vector< size_t> &CHOSEN_VALUES,
        storage::Vector< size_t> &SOLUTION
      ) const;

      //! @brief choose a combination of subgraphs (of a supergraph) which do not overlap, given 
      //!        constraints on which subgraphs can be chosen for each solution slot (recursive algorithm) 
      //! @param ALLOWED_SUBGRAPHS subgraphs (i.e. lists of indices) that are allowed for each solution slot
      //! @param SOLN_NO which solution slot should be chosen during the current iteration
      //! @param COVERED_VERTICES which vertices are already covered by chosen subgraphs
      //! @param SOLUTION which subgraph should be used for each solution slot (used to return info) 
      //! @return true if non-overlapping sugraphs for current and subsequent solution slots could be found
      bool ChooseSubgraphCombination
      (
        const storage::Vector< storage::Vector< storage::Vector< size_t> > > &ALLOWED_SUBGRAPHS,
        const size_t &SOLN_NO,
        storage::Vector< size_t> &COVERED_VERTICES,
        storage::Vector< size_t> &SOLUTION
      ) const;

      //! @brief caches the reactants and products of a reaction if necessary
      //! @param REACTION the reaction to cache
      void Cache( const ReactionComplete &REACTION) const;

      //! @brief gets the reactant structures from the reaction; retrieves from cache if available, stores them if not
      const storage::Vector< ReactionStructure> &GetReactantStructures( const ReactionComplete &REACTION) const;

      //! @brief gets the product structures from the reaction; retrieves from cache if available, stores them if not
      const storage::Vector< ReactionStructure> &GetProductStructures( const ReactionComplete &REACTION) const;

    }; // class ReactionWorker

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_REACTION_WORKER_H_
