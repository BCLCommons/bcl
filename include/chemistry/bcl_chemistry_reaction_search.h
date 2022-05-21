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

#ifndef BCL_CHEMISTRY_REACTION_SEARCH_H_
#define BCL_CHEMISTRY_REACTION_SEARCH_H_

// include the namespace header
#include "bcl_chemistry.h"

// includes from bcl - sorted alphabetically
#include "bcl_chemistry_conformation_graph_converter.h"
#include "bcl_chemistry_fragment_ensemble.h"
#include "bcl_chemistry_reaction_ensemble.h"
#include "bcl_chemistry_reaction_worker.h"
#include "bcl_chemistry_rotamer_library_file.h"
#include "bcl_chemistry_search_fragment_library_from_tree.h"
#include "graph/bcl_graph_const_graph.h"
#include "sched/bcl_sched_mutex.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_si_ptr.h"
#include "util/bcl_util_si_ptr_list.h"
#include "util/bcl_util_si_ptr_vector.h"

namespace bcl
{
  namespace chemistry
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ReactionSearch
    //! @brief stores reactions and associated available reactants.  stored in such a way that looking up
    //!        reactants that match a reaction or reactions that match a molecule can be done rapidly.
    //!        internally stored as a graph of associations between molecular substructures 
    //!
    //! @see @link example_chemistry_reaction_search.cpp @endlink
    //! @author geanesar
    //! @date Feb 18, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ReactionSearch :
      public util::ObjectInterface 
    {
    public:

      // typedefs for convenience. *_p = SiPtr<*>, *_sp = ShPtr<*>, *_pv = SiPtrVector<*>
      typedef util::ShPtr< ReactionStructure> ReactionStructure_sp;
      typedef util::SiPtr< const ReactionStructure> ReactionStructure_p;
      typedef util::SiPtrVector< const ReactionStructure> ReactionStructure_pv;

      typedef util::ShPtr< ReactionComplete> ReactionComplete_sp;
      typedef util::SiPtr< const ReactionComplete> ReactionComplete_p;
      typedef util::SiPtrVector< const ReactionComplete> ReactionComplete_pv;

      typedef util::SiPtr< FragmentComplete> FragmentComplete_sp;
      typedef util::SiPtr< const FragmentComplete> FragmentComplete_p;
      typedef util::SiPtrVector< const FragmentComplete> FragmentComplete_pv;

    private:

    ////////////////////
    // helper classes //
    ////////////////////
      
    //////////
    // data //
    //////////

      //! a singleton for this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      // mutex for initialization
      mutable sched::Mutex m_Mutex; 

      //! whether the class needs to be initialized
      mutable bool m_Initialized;

      //! the filename of where to get reactants from
      std::string m_ReactantFilename;

      //! the filename of where to get reactions from
      std::string m_ReactionDirectory;

      //! pool of molecules available as reactants in reactions
      mutable util::ShPtr< storage::Vector< FragmentComplete> > m_Molecules;

      //! reactions that can be associated with reactants
      mutable util::ShPtr< storage::Vector< ReactionComplete> > m_Reactions;

      //! reaction structures that will be searched
      mutable util::ShPtr< storage::Vector< ReactionStructure> > m_ReactionStructures;

      //! the reaction tree; relates substructures to the reactions and molecules that match them,
      //! and how the substructures relate to each other 
      mutable graph::ConstGraph< ReactionStructure_p, int> m_ReactionStructureTree;

      //! the vertices to start any structure search at (have 0 neighbors)
      mutable storage::Vector< size_t> m_StartSearchVertices;

      //! a reverse map of reaction structures to their reactions and which reactants they match
      mutable storage::Map
      <
        ReactionStructure_p,
        storage::Map< ReactionComplete_p, storage::Set< size_t> >
      > m_AvailableReactions;

      //! a map from reactions (keys) to available reactants (values) for each reactant number
      mutable storage::Map
      < 
        ReactionComplete_p,
        storage::Vector< FragmentComplete_pv>
      > m_AvailableReactants;

      //! whether to remove reactions which have zero reactants
      bool m_RemoveDeficientRXNs;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor using existing reactants and reactions
      //! @param REACTANTS the reactant molecules that should be available for reactions
      //! @param REACTIONS the reactions that are available
      ReactionSearch
      (
        const FragmentEnsemble &REACTANTS = FragmentEnsemble(),
        const ReactionEnsemble &REACTIONS = ReactionEnsemble(),
        const bool &REMOVE_DEFICIENT_RXNS = true
      );

      //! @brief constructor to read reactants and reactions from a file
      //! @param REACTANT_FILENAME name of file containing reactants
      //! @param REACTION_DIRECTORY name of file containing reactions
      ReactionSearch
      (
        const std::string &REACTANT_FILENAME,
        const std::string &REACTION_DIRECTORY
      );

      //! @brief constructor using existing reactants and reading new reactions from file
      //! @param REACTANTS the reactant molecules that should be available for reactions
      //! @param REACTIONS the reactions that are available
      ReactionSearch
      (
        const FragmentEnsemble &REACTANTS,
        const std::string &REACTION_DIRECTORY,
        const bool &REMOVE_DEFICIENT_RXNS = true
      );

      //! @brief copy constructor
      ReactionSearch( const ReactionSearch &OTHER);

      //! @brief clone constructor
      //! @return a pointer to a copy of this class
      ReactionSearch *Clone() const;
    
    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief clear internal data
      void Reset();

      //! @brief build the reaction structure tree
      void Initialize() const;

      //! @brief test whether this class is empty
      //! @return true if the class has no reactions or molecules in it
      bool IsEmpty() const;

      //! @brief get internal data
      const graph::ConstGraph< ReactionStructure_p, int> &GetStructureTree() const
      {
        return m_ReactionStructureTree;
      }

      //! @brief get a vector of reactions that are stored in this class
      //! @return a vector of reactions
      const util::ShPtr< storage::Vector< ReactionComplete> > &GetReactions() const;

      //! @brief get a vector of reactants that are stored in this class
      //! @return a vector of molecules
      const util::ShPtr< storage::Vector< FragmentComplete> > &GetMolecules() const;

      //! @brief set an ensemble to become the new reactants in the catalog
      //! @param MOLECULES the molecules to become the catalog
      void SetMolecules( const util::ShPtr< storage::Vector< FragmentComplete> > &MOLECULES);

      //! @brief adds a molecule to the reaction catalog
      //! @param MOLECULE the molecule to add to the catalog
      void AddMolecule( const FragmentComplete &MOLECULE);

      //! @brief finds reactions that a molecule can participate in
      //! @param MOLECULE the query molecule
      //! @return a list of ReactionStructureNodes containing matching reactions 
      util::SiPtrVector< const ReactionComplete> FindReactions( const ConformationInterface &MOLECULE) const;

      //! @brief removes reactions that are incapable of providing molecules for every reactant
      ReactionComplete_pv RemoveReactantDeficientReactions() const;
      
      //! @brief gets a list of available reactants for a reaction
      //! @param REACTION a pointer to a reaction of interest
      //! @return a vector of pointers to molecules which are available for each reactant in the given reaction
      //! outer vector is reactant index in the reaction and inner vector are fragments that match the reactant
      //! for that index
      const storage::Vector< FragmentComplete_pv> &GetAvailableReactants
      (
        const ReactionComplete_p &REACTION
      ) const;

      //! @brief Choose a random reaction that a molecule will participate in
      //! @param MOL the molecule to match
      //! @return a pointer to a ReactionComplete that can use MOL as a reactant
      ReactionComplete_p ChooseRandomAvailableRxn
      ( 
        const FragmentComplete &MOL
      ) const;

      //! @brief select random reactants for each position of a given reaction
      //! @param RXN_P pointer to the reaction of interest
      //! @return returns pointer vector of reactants where the vector position corresponds to
      //! the position of the reactant in the reaction
      FragmentComplete_pv ChooseRandomReactants
      ( 
        const ReactionComplete_p &RXN_P
      ) const;

      //! @brief choose a random reaction that a molecule can participate in, and other reaction partners
      //! @param MOLECULE the molecule to search for
      //! @return pointer to ReactionComplete and pointer vector to partner reactants
      storage::Pair< ReactionComplete_p, FragmentComplete_pv> ChooseRandomRxnAndReactants
      ( 
        const FragmentComplete &MOLECULE
      ) const;

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    private:

    //////////////////////
    // helper functions //
    //////////////////////
      
      //! @brief locate a ReactionStructure
      size_t FindReactionStructure( const ReactionStructure &STRUCT) const;

      //! @brief create the structure tree from reactions
      void MakeStructureTree() const;

      //! @brief update vertices that any searches will begin from 
      void DetermineSearchVertices() const;

      //! @brief adds a molecule to a matching set of reactions
      //! @param MOLECULE the molecule to add
      void AddMoleculeToMatchedReactions( const FragmentComplete &MOLECULE) const;
      
      //! @brief update the available reactants mapping (which reactants can be used for each reaction)
      void DetermineAvailableReactants() const;

      //! @brief read RXN files from a directory to create reactions
      //! @param DIRECTORY the directory to use.
      //! @return a shptr to a vector of ReactionCompletes built from RXN files
      static util::ShPtr< storage::Vector< ReactionComplete> > ReadReactions( const std::string &DIRECTORY);

      //! @brief read the reactions and reactants from files
      void ReadFiles( const bool &FORCE = false) const;
      
      //! @brief retrieve the available reactant map
      //! @return the available reactant map
      const storage::Map
      < 
        ReactionComplete_p,
        storage::Vector< FragmentComplete_pv>
      > &GetAvailableReactantsMap() const;

      //! @brief depth-first search of the reaction trees for a query
      //! @param QUERY the query to search
      //! @return vector of ptrs to matching ReactionStructure objects
      util::SiPtrVector< const ReactionStructure> FindMatchingStructures
      ( 
        const ConformationInterface &QUERY
      ) const;

      //! @brief helper function for FindMatchingStructures
      //! @param QUERY the query to search for
      //! @param CURRENT_NODE index of the current structure node in the graph
      //! @param CHECKED_NODES a list of nodes which have already been visited
      //! @param MATCHING_STRUCTS a pointer vector populated with reaction structures matching QUERY
      void FindMatchingStructuresRecurse
      ( 
        const ConformationInterface &QUERY,
        const size_t &CURRENT_NODE,
        storage::Vector< size_t> &CHECKED_NODES,
        ReactionStructure_pv &MATCHING_STRUCTS
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

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

    }; // class ReactionSearch

  } // namespace chemistry
} // namespace bcl

#endif // BCL_CHEMISTRY_REACTION_SEARCH_H_
