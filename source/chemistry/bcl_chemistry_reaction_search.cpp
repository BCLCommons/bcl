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
#include "chemistry/bcl_chemistry_reaction_search.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "graph/bcl_graph_edge_cover_ring_perception.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_rxn_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //
    // Reaction Catalog
    //
    
    typedef ReactionSearch::ReactionStructure_sp ReactionStructure_sp;
    typedef ReactionSearch::ReactionStructure_p ReactionStructure_p;
    typedef ReactionSearch::ReactionStructure_pv ReactionStructure_pv;

    typedef ReactionSearch::ReactionComplete_sp ReactionComplete_sp;
    typedef ReactionSearch::ReactionComplete_p ReactionComplete_p;
    typedef ReactionSearch::ReactionComplete_pv ReactionComplete_pv;

    typedef ReactionSearch::FragmentComplete_sp FragmentComplete_sp;
    typedef ReactionSearch::FragmentComplete_p FragmentComplete_p;
    typedef ReactionSearch::FragmentComplete_pv FragmentComplete_pv;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! @param REACTANTS the reactant molecules that should be available for reactions
    //! @param REACTIONS the reactions that are available
    ReactionSearch::ReactionSearch
    (
      const FragmentEnsemble &REACTANTS,
      const ReactionEnsemble &REACTIONS,
      const bool &REMOVE_DEFICIENT_RXNS
    ) :
      m_Initialized( false),
      m_ReactantFilename(),
      m_ReactionDirectory(),
      m_Molecules( new storage::Vector< FragmentComplete>( REACTANTS.Begin(), REACTANTS.End())),
      m_Reactions( new storage::Vector< ReactionComplete>( REACTIONS.Begin(), REACTIONS.End())),
      m_ReactionStructures( new storage::Vector< ReactionStructure>()),
      m_ReactionStructureTree(),
      m_StartSearchVertices(),
      m_AvailableReactants(),
      m_RemoveDeficientRXNs( REMOVE_DEFICIENT_RXNS)
    {
    }

    //! @brief constructor
    //! @param REACTANT_FILENAME name of file containing reactants
    //! @param REACTION_DIRECTORY name of file containing reactions
    ReactionSearch::ReactionSearch
    (
      const std::string &REACTANT_FILENAME,
      const std::string &REACTION_DIRECTORY
    ) :
      m_Initialized( false),
      m_ReactantFilename( REACTANT_FILENAME),
      m_ReactionDirectory( REACTION_DIRECTORY),
      m_Molecules( new storage::Vector< FragmentComplete>()),
      m_Reactions( new storage::Vector< ReactionComplete>()),
      m_ReactionStructures( new storage::Vector< ReactionStructure>()),
      m_ReactionStructureTree(),
      m_StartSearchVertices(),
      m_AvailableReactions(),
      m_AvailableReactants() 
    {
    }

    //! @brief constructor using existing reactants and reading new reactions from file
    //! @param REACTANTS the reactant molecules that should be available for reactions
    //! @param REACTIONS the reactions that are available
    ReactionSearch::ReactionSearch
    (
      const FragmentEnsemble &REACTANTS,
      const std::string &REACTION_DIRECTORY,
      const bool &REMOVE_DEFICIENT_RXNS
    ) :
      m_Initialized( false),
      m_ReactantFilename( ),
      m_ReactionDirectory( REACTION_DIRECTORY),
      m_Molecules( new storage::Vector< FragmentComplete>( REACTANTS.Begin(), REACTANTS.End())),
      m_Reactions( new storage::Vector< ReactionComplete>()),
      m_ReactionStructures( new storage::Vector< ReactionStructure>()),
      m_ReactionStructureTree(),
      m_StartSearchVertices(),
      m_AvailableReactions(),
      m_AvailableReactants(),
      m_RemoveDeficientRXNs( REMOVE_DEFICIENT_RXNS)
    {
    }

    void ReactionSearch::Reset()
    {
      m_Mutex.Lock();
      m_Initialized = false;
      m_Molecules = util::ShPtr< storage::Vector< FragmentComplete> >();
      m_Reactions = util::ShPtr< storage::Vector< ReactionComplete> >();
      m_ReactionStructures = util::ShPtr< storage::Vector< ReactionStructure> >();
      m_ReactionStructureTree = graph::ConstGraph< ReactionStructure_p, int>();
      m_StartSearchVertices = storage::Vector< size_t>();
      m_AvailableReactions = storage::Map
                             <
                               ReactionStructure_p,
                               storage::Map< ReactionComplete_p, storage::Set< size_t> >
                             >();
      m_AvailableReactants = storage::Map
                             < 
                               ReactionComplete_p, 
                               storage::Vector< FragmentComplete_pv> 
                             >();
      m_Mutex.Unlock();
    }

    //! @brief copy constructor
    ReactionSearch::ReactionSearch( const ReactionSearch &OTHER) :
      m_ReactantFilename( OTHER.m_ReactantFilename),
      m_ReactionDirectory( OTHER.m_ReactionDirectory)
    {

      //BCL_Assert( false, "Bailout");
      //OTHER.Initialize();
      
      // lock this to copy the molecule, reaction, and reaction structure sh pointers
      OTHER.m_Mutex.Lock();

      m_Molecules = OTHER.m_Molecules;
      m_Reactions = OTHER.m_Reactions;
      m_ReactionStructures = OTHER.m_ReactionStructures;

      OTHER.m_Mutex.Unlock();
      
      if( ( m_Molecules.IsDefined() && m_Reactions.IsDefined() && m_ReactionStructures.IsDefined()))
      {

        m_ReactionStructureTree = OTHER.m_ReactionStructureTree;
        m_StartSearchVertices = OTHER.m_StartSearchVertices;
        m_AvailableReactions = OTHER.m_AvailableReactions;
        m_AvailableReactants = OTHER.m_AvailableReactants;
        m_Initialized = OTHER.m_Initialized;

      }
      else
      {
        Reset();
      }
    }

    //! @brief clone constructor
    //! @return a pointer to a copy of this class
    ReactionSearch *ReactionSearch::Clone() const
    {
      return new ReactionSearch( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ReactionSearch::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief test whether this class is empty
    //! @return true if the class has no reactions or no molecules in it
    bool ReactionSearch::IsEmpty() const
    {
      Initialize();
      return ( m_Reactions.IsDefined() && m_Molecules.IsDefined())
        ? m_Reactions->IsEmpty() || m_Molecules->IsEmpty()
        : true;
    }

    //! @brief get a vector of reactions that are stored in this class
    //! @return a vector of reactions
    const util::ShPtr< storage::Vector< ReactionComplete> > &ReactionSearch::GetReactions() const
    {
      return m_Reactions;
    }

    //! @brief get a vector of reactions that are stored in this class
    //! @return a vector of reactions
    const util::ShPtr< storage::Vector< FragmentComplete> > &ReactionSearch::GetMolecules() const
    {
      Initialize();
      return m_Molecules;
    }

    //! @brief set an ensemble to become the new reactants in the catalog
    //! @param MOLECULES the molecules to become the catalog
    void ReactionSearch::SetMolecules( const util::ShPtr< storage::Vector< FragmentComplete> > &MOLECULES)
    {
      m_Molecules = MOLECULES;
    }

    //! @brief adds a molecule to the reaction catalog
    //! @param MOLECULE the molecule to add to the catalog
    void ReactionSearch::AddMolecule( const FragmentComplete &MOLECULE)
    {
      m_Molecules->PushBack( MOLECULE);
    }

    //! @brief finds reactions that a molecule can participate in
    //! @param MOLECULE the query molecule
    //! @return a list of ReactionStructureNodes containing matching reactions 
    util::SiPtrVector< const ReactionComplete> ReactionSearch::FindReactions
    ( 
      const ConformationInterface &MOLECULE
    ) const
    {
      Initialize();

      ReactionStructure_pv found_reactions;
      if( MOLECULE.GetNumberAtoms() == 0)
      {
        return found_reactions; 
      }

      ReactionStructure_pv found_structs( FindMatchingStructures( MOLECULE));

      util::SiPtrVector< const ReactionComplete> rxn_ptrs;

      // Find reactions associated with the matched structures
      for( size_t str_idx( 0), end_idx( found_structs.GetSize()); str_idx < end_idx; ++str_idx)
      {
        const ReactionStructure_p &str_ptr( found_structs( str_idx));
        if( m_AvailableReactions.Has( str_ptr))
        {
          const storage::Map< ReactionComplete_p, storage::Set< size_t> > &map( m_AvailableReactions.GetValue( str_ptr));
          size_t total( 0);
          for
          (
            storage::Map< ReactionComplete_p, storage::Set< size_t> >::const_iterator itr_map( map.Begin()), 
              itr_map_end( map.End());
            itr_map != itr_map_end;
            ++itr_map, ++total
          )
          {
            rxn_ptrs.PushBack( itr_map->first);
          }
          BCL_MessageDbg( "  Found " + util::Format()( total) + " reactions associated with a structure");
        }
        else
        {
          BCL_MessageStd( "A structure wasn't found in the reactions map, something is wrong...");
        }
      }

      return rxn_ptrs;
    }

    //! @brief gets a list of available reactants for a reaction
    //! @param REACTION a pointer to a reaction of interest
    //! @return a vector of pointers to molecules which are available for each reactant in the given reaction;
    //! outer vector is reactant index in the reaction and inner vector are fragments that match the reactant
    //! for that index
    const storage::Vector< FragmentComplete_pv> &ReactionSearch::GetAvailableReactants
    (
      const ReactionComplete_p &REACTION
    ) const
    {
      static storage::Vector< FragmentComplete_pv> empty_v;

      BCL_MessageVrb( "Initialize reactions");
      Initialize();

      if( REACTION.IsDefined() && m_AvailableReactants.Has( REACTION))
      {
        return m_AvailableReactants.GetValue( REACTION);
      }
      return empty_v; 
    }

    ReactionComplete_p ReactionSearch::ChooseRandomAvailableRxn
    ( 
      const FragmentComplete &MOL
    ) const
    {
      Initialize();
      ReactionComplete_p res;

      // no atoms in input struct, bail
      if( !MOL.GetNumberAtoms())
      {
        BCL_MessageVrb( "ReactionSearch No atoms in molecule");
        return res;
      }

      util::SiPtrVector< const ReactionComplete> available_rxns( FindReactions( MOL));
      if( available_rxns.IsEmpty())
      {
        BCL_MessageVrb( "ReactionSearch: ChooseRandomAvailableReaction couldn't find any available reactions");
        return res;
      }

      // randomly pick a reaction
      size_t picked_rxn( 0);
      if( available_rxns.GetSize() > 1)
      {
        picked_rxn = random::GetGlobalRandom().Random< size_t>( 0, available_rxns.GetSize() - 1);
      }

      // return a random reaction that MOL can participate in
      res = available_rxns( picked_rxn);
      return res;
    }

    //! @brief select random reactants for each position of a given reaction
    //! @param RXN_P pointer to the reaction of interest
    //! @return returns pointer vector of reactants where the vector position corresponds to
    //! the position of the reactant in the reaction
    FragmentComplete_pv ReactionSearch::ChooseRandomReactants
    ( 
      const ReactionComplete_p &RXN_P
    ) const
    {
      FragmentComplete_pv res;

      if( !RXN_P.IsDefined())
      {
        return res;
      }

      size_t n_rxt_pos( RXN_P->GetNumberReactants());
      res.AllocateMemory( n_rxt_pos);

      storage::Vector< FragmentComplete_pv> avail_reacts
      (
        GetAvailableReactants( *RXN_P) 
      );

      if( avail_reacts.GetSize() != n_rxt_pos)
      {
        return res;
      }

      for( size_t i( 0), l( avail_reacts.GetSize()); i < l; ++i)
      {
        size_t rxt_size( avail_reacts( i).GetSize());
        if( !rxt_size)
        {
          break;
        }

        size_t picked_rxt( rxt_size > 1 ? random::GetGlobalRandom().Random< size_t>( 0, rxt_size - 1) : 0);
        res.PushBack( avail_reacts( i)( picked_rxt));
      }

      // if
      if( res.GetSize() != n_rxt_pos)
      {
        res.Reset();
      }

      return res;
    }

    storage::Pair< ReactionComplete_p, FragmentComplete_pv> 
      ReactionSearch::ChooseRandomRxnAndReactants
      ( 
        const FragmentComplete &MOLECULE
      ) const
    {
      storage::Pair< ReactionComplete_p, FragmentComplete_pv> res; 

      res.First() = ChooseRandomAvailableRxn( MOLECULE);
      ReactionComplete_p &rxn_p( res.First()); // reaction pointer

      // if molecule is empty then rxn_p will be undefined, so this also acts as a # atoms check
      if( !rxn_p.IsDefined())
      {
        // reaction didn't exist
        return res;
      }

      res.Second() = ChooseRandomReactants( rxn_p); // random reactants for each slot in the reaction
      FragmentComplete_pv &rxts_p( res.Second()); // pointers to reactant structures

      if( rxts_p.IsEmpty())
      {
        res.First() = ReactionComplete_p();
      }

      return res;
    }

    //! @brief detects reactions that are incapable of providing molecules for every reactant
    util::SiPtrVector< const ReactionComplete> ReactionSearch::RemoveReactantDeficientReactions() const
    {
      BCL_MessageStd( "Removing deficient reactions");
      util::SiPtrVector< const ReactionComplete> deficient_rxns;

      for
      (
        storage::Map< ReactionComplete_p, storage::Vector< FragmentComplete_pv> >::iterator
          itr_map( m_AvailableReactants.Begin()), itr_map_end( m_AvailableReactants.End());
        itr_map != itr_map_end;
        ++itr_map
      )
      {
        for( size_t i( 0), end_i( itr_map->second.GetSize()); i < end_i; ++i)
        {
          if( itr_map->second( i).IsEmpty())
          {
            deficient_rxns.PushBack( itr_map->first);
            break;
          }
        }
      }

      if( !deficient_rxns.IsEmpty())
      {

        // remove deficient reactions from the 
        for
        ( 
          storage::Map< ReactionStructure_p, storage::Map< ReactionComplete_p, storage::Set< size_t> > >::iterator
            itr_map( m_AvailableReactions.Begin()), itr_map_end( m_AvailableReactions.End());
          itr_map != itr_map_end;
          ++itr_map
        )
        {
          for( size_t r( 0), end_r( deficient_rxns.GetSize()); r < end_r; ++r)
          {
            // remove deficient reactions from the available reactions map
            itr_map->second.Erase( deficient_rxns( r));
          }
        }

        // remove empty reactions from the available_reactions map
        for( size_t r( 0), end_r( deficient_rxns.GetSize()); r < end_r; ++r)
        {

          m_AvailableReactions.Erase( deficient_rxns( r));
        }
      }

      return deficient_rxns;
    }
      
    //! @brief Build the structure tree
    void ReactionSearch::Initialize() const
    {
      // quick check to see if this instance is initialized already; don't do anything if it is
      // first check determines if the mutex needs to be locked at all
      m_Mutex.Lock();
      if( !m_Initialized)
      {
        BCL_MessageStd( "Initializing reaction search database");

        // read any necessary files
        if( !m_Molecules.IsDefined() || !m_Reactions.IsDefined() || m_Molecules->IsEmpty() || m_Reactions->IsEmpty())
        {
          ReadFiles();
        }

        // make the structure tree from reactant substructures present in m_Reactions
        BCL_MessageStd( "Make structure tree");
        MakeStructureTree();

        m_Initialized = true;

        // prune empty reactions if necessary
        if( m_RemoveDeficientRXNs)
        {
          util::SiPtrVector< const ReactionComplete> removed( RemoveReactantDeficientReactions());
          BCL_MessageVrb
          ( 
            "Removed " + util::Format()( removed.GetSize()) + " reactant-deficient reactions from the search tree"
          );
        }

        // Determine which molecules match reaction substructures
        DetermineAvailableReactants();

      }
      m_Mutex.Unlock();
    }

    size_t ReactionSearch::FindReactionStructure( const ReactionStructure &STRUCT) const
    {
      size_t index( util::GetUndefined< size_t>());
      if( m_ReactionStructures.IsDefined())
      {
        storage::Vector< ReactionStructure> &rxn_structs( *m_ReactionStructures);
        size_t end( rxn_structs.GetSize());
        for( index = 0; index < end; ++index)
        {
          if( STRUCT == rxn_structs( index))
          {
            break;
          }
        }
      }
      return index;
    }

    void ReactionSearch::MakeStructureTree() const
    {
      BCL_MessageStd( "Constructing the search tree");
      if( !( m_Reactions.IsDefined() && m_Molecules.IsDefined()))
      {
        BCL_MessageStd( "No defined reactions or molecules for MakeStructureTree");
        return;
      }
      m_ReactionStructures = util::ShPtr< storage::Vector< ReactionStructure> >( new storage::Vector< ReactionStructure>());

      storage::Map< size_t, storage::Set< size_t> > rxn_edges; // smaller structure to larger structure

      storage::Vector< ReactionStructure> structures; // temporarily holds reactant structures
      storage::Vector< storage::Map< ReactionComplete_p, storage::Set< size_t> > > struct_rxn_map; // maps reactant structures to reactions and their indices

      // Make one node for each reactant structure in each reaction
      for
      (
        storage::Vector< ReactionComplete>::const_iterator itr_rxn( m_Reactions->Begin()), itr_rxn_end( m_Reactions->End());
        itr_rxn != itr_rxn_end;
        ++itr_rxn
      )
      {
        const storage::Vector< FragmentComplete> &reactants( itr_rxn->GetReactants());
        ReactionComplete_p rxn_pointer( &( *itr_rxn));

        size_t mol_no( 0);
        for
        (
          storage::Vector< FragmentComplete>::const_iterator itr_mol( reactants.Begin()), itr_mol_end( reactants.End());
          itr_mol != itr_mol_end;
          ++itr_mol, ++mol_no
        )
        {
          FragmentComplete mol( *itr_mol);
          if( !mol.GetNumberAtoms())
          {
            BCL_MessageStd( "Reaction " + itr_rxn->GetDescription() + " reactant " + util::Format()( mol_no) + " was empty");
            continue;
          }

          ReactionStructure rxt_struct( mol);

          if( rxt_struct.GetSize() != mol.GetNumberAtoms())
          {
            BCL_MessageStd( "Reactant number " + util::Format()( mol_no) + " in reaction: " + itr_rxn->GetDescription() + " could not make a proper reactant structure");
            continue;
          }
          
          // determine if this structure is repeated
          bool should_insert( true);
          for( size_t i( 0), end( structures.GetSize()); i < end; ++i)
          {
            if( rxt_struct == structures( i))
            {
              struct_rxn_map( i)[ rxn_pointer].Insert( mol_no); // associate this structure with the reaction and the reactant index
              should_insert = false;
              break;
            }
          }

          // if the structure wasn't found, save it
          if( should_insert)
          {
            structures.PushBack( rxt_struct);
            struct_rxn_map.PushBack( storage::Map< ReactionComplete_p, storage::Set< size_t> >());
            struct_rxn_map.LastElement()[ rxn_pointer].Insert( mol_no);
          }
        }
      }
      *m_ReactionStructures = structures;
      storage::Vector< ReactionStructure_p> graph_vertices;
      graph_vertices.AllocateMemory( m_ReactionStructures->GetSize());

      // vertex_order pairs are <size of structure, index of structure>
      std::vector< std::pair< size_t, size_t> > vertex_order;
      vertex_order.reserve( graph_vertices.GetSize());

      // copy the structure->reaction maps
      for( size_t i( 0), end( m_ReactionStructures->GetSize()); i < end; ++i)
      {
        BCL_MessageStd( "ReactionStructure " + util::Format()( i));
        ReactionStructure_p str_ptr( &( *m_ReactionStructures)( i));
        if( str_ptr.IsDefined())
        {
          graph_vertices.PushBack( str_ptr); // add a pointer to the structure
          vertex_order.push_back( std::pair< size_t, size_t>( str_ptr->GetSize(), graph_vertices.GetSize() - 1));
          m_AvailableReactions[ str_ptr] = struct_rxn_map( i);
        }
      }

      // determine vertex order
      std::sort( vertex_order.begin(), vertex_order.end());

      // make a graph with no connections first, and leave it undirected for ring detection
      storage::Vector< graph::UndirectedEdge< int> > empty_edges;
      m_ReactionStructureTree = graph::ConstGraph< ReactionStructure_p, int>( graph_vertices, empty_edges, int( 0));

      // connect smaller structures to larger ones
      for( size_t i( 0), end( vertex_order.size()); i < end; ++i)
      {
        size_t smaller_vertex( vertex_order[ i].second);
        const ReactionStructure_p &smaller( m_ReactionStructureTree.GetVertexData( smaller_vertex));

        // compare the smallest vertex against all larger vertices
        for( size_t j( i + 1); j < end; ++j)
        {
          size_t larger_vertex( vertex_order[ j].second);
          const ReactionStructure_p &larger( m_ReactionStructureTree.GetVertexData( larger_vertex));
          if( larger->Contains( *smaller))
          {
            m_ReactionStructureTree.AddEdge( smaller_vertex, larger_vertex, int( 1));
          }
        }
      }

      // since this should be an acyclic graph, remove any cycles by removing edges between smallest and largest vertices
      for
      ( 
        storage::List< graph::Ring> rings( graph::EdgeCoverRingPerception( m_ReactionStructureTree).GetRings());
        !rings.IsEmpty();
        rings = graph::EdgeCoverRingPerception( m_ReactionStructureTree).GetRings()
      )
      {
        BCL_MessageVrb( "  Found " + util::Format()( rings.GetSize()) + " cycles in the structure graph");

        // Break each ring
        for
        (
          storage::List< graph::Ring>::const_iterator itr_ring( rings.Begin()), itr_ring_end( rings.End());
          itr_ring != itr_ring_end;
          ++itr_ring
        )
        {
          const storage::Vector< size_t> &ring_vertices( itr_ring->GetVertices());
          
          // find the smallest vertex
          size_t smallest( 0);
          size_t smallest_size( m_ReactionStructureTree.GetVertexData( ring_vertices( smallest))->GetSize());
          for( size_t i( 1), end( ring_vertices.GetSize()); i < end; ++i)
          {
            if( m_ReactionStructureTree.GetVertexData( ring_vertices( i))->GetSize() < smallest_size)
            {
              smallest = i;
              smallest_size = m_ReactionStructureTree.GetVertexData( ring_vertices( i))->GetSize();
            }
          }

          // find the largest neighbor of the smallest vertex
          std::pair< size_t, size_t> to_check;
          if( smallest == 0)
          {
            to_check.first = ring_vertices.GetSize() - 1;
          }
          else
          {
            to_check.first = smallest - 1;
          }
          if( smallest == ring_vertices.GetSize() - 1)
          {
            to_check.second = 0;
          }
          else
          {
            to_check.second = smallest + 1;
          }

          // remove connection between smallest and largest vertices
          if( m_ReactionStructureTree.GetVertexData( ring_vertices( to_check.first))->GetSize() > m_ReactionStructureTree.GetVertexData( ring_vertices( to_check.second))->GetSize())
          {
            m_ReactionStructureTree.RemoveEdge( ring_vertices( smallest), ring_vertices( to_check.first));
          }
          else
          {
            m_ReactionStructureTree.RemoveEdge( ring_vertices( smallest), ring_vertices( to_check.second));
          }

        }
      }

      // Set the graph to directed, and update the edges so that smaller structures point toward larger ones with positive (+1) connectivity
      // and larger structures point toward smaller one with negative (-1) connectivity
      m_ReactionStructureTree.SetDirected();
      for( size_t v( 0), end_v( m_ReactionStructureTree.GetSize()); v < end_v; ++v)
      {
        const ReactionStructure_p &v_data( m_ReactionStructureTree.GetVertexData( v));
        // get the neighbors of this vertex
        storage::Vector< size_t> neighbors( m_ReactionStructureTree.GetNeighborIndices( v));
        for( size_t n( 0), end_n( neighbors.GetSize()); n < end_n; ++n)
        {
          // only consider vertices with higher indices than the current one
          if( neighbors( n) < v)
          {
            continue;
          }

          const ReactionStructure_p &neighbor_data( m_ReactionStructureTree.GetVertexData( neighbors( n)));
          m_ReactionStructureTree.RemoveEdge( v, neighbors( n));
          m_ReactionStructureTree.RemoveEdge( neighbors( n), v);
          if( v_data->GetSize() <= neighbor_data->GetSize())
          {
            BCL_MessageDbg( "Directing " + util::Format()( v) + " to " + util::Format()( neighbors( n)));
            m_ReactionStructureTree.AddEdge( v, neighbors( n), int( 1));
            m_ReactionStructureTree.AddEdge( neighbors( n), v, int( -1));
          }
          else
          {
            BCL_MessageDbg( "Directing " + util::Format()( neighbors( n)) + " to " + util::Format()( v));
            m_ReactionStructureTree.AddEdge( v, neighbors( n), int( -1));
            m_ReactionStructureTree.AddEdge( neighbors( n), v, int( 1));
          }
        }
      }
      BCL_MessageStd( "Done setting up tree");
      DetermineSearchVertices();
    }

    //! @brief update vertices that any searches will begin from 
    void ReactionSearch::DetermineSearchVertices() const
    {
      m_StartSearchVertices.Reset();

      // find any structures that ONLY have outgoing (larger) structures attached to them
      //std::stringstream msg;
      //msg << "Search vertices: ";
      for( size_t v( 0), end_v( m_ReactionStructureTree.GetSize()); v < end_v; ++v)
      {
        const storage::Vector< size_t> &neighbors( m_ReactionStructureTree.GetNeighborIndices( v));
        bool keep_vertex( true);
        for( size_t n( 0), end_n( neighbors.GetSize()); n < end_n && keep_vertex; ++n)
        {
          if( m_ReactionStructureTree.GetEdgeData( v, neighbors( n)) < 0)
          {
            keep_vertex = false;
          }
        }

        if( keep_vertex)
        {
          m_StartSearchVertices.PushBack( v);
        }
      }
      BCL_MessageStd( "There are " + util::Format()( m_StartSearchVertices.GetSize()) + " search starting points");
    }

    //! @brief update the available reactants mapping (which reactants can be used for each reaction)
    void ReactionSearch::DetermineAvailableReactants() const
    {
      if( !( m_Molecules.IsDefined()))
      {
        BCL_MessageStd( "DetermineAvailableReactants - molecules undefined!");
        return;
      }

      BCL_MessageStd( "Associating reactants with reactions");
      for( size_t mol_idx( 0), last_idx( m_Molecules->GetSize()); mol_idx < last_idx; ++mol_idx)
      {
        if( !( mol_idx % 1000))
        {
          util::GetLogger().LogStatus( "Associated " + util::Format()( mol_idx) + " molecules");
        }
        BCL_MessageVrb( "Searching molecule " + util::Format()( mol_idx));
        storage::Vector< size_t> checked_nodes( m_ReactionStructureTree.GetSize(), size_t( 0));

        ReactionStructure_pv matching_structs( FindMatchingStructures( ( *m_Molecules)( mol_idx)));

        // Collect reactions
        for( size_t str_idx( 0), end_idx( matching_structs.GetSize()); str_idx < end_idx; ++str_idx)
        {
          const ReactionStructure_p &str( matching_structs( str_idx));
          if( m_AvailableReactions.Has( str))
          {
            // associate this molecule with this reaction
            const storage::Map< ReactionComplete_p, storage::Set< size_t> > &map( m_AvailableReactions.GetValue( str));
            for
            (
              storage::Map< ReactionComplete_p, storage::Set< size_t> >::const_iterator itr_map( map.Begin()), 
                itr_map_end( map.End());
              itr_map != itr_map_end;
              ++itr_map
            )
            {
              for
              (
                storage::Set< size_t>::const_iterator itr_set( itr_map->second.Begin()), itr_set_end( itr_map->second.End());
                itr_set != itr_set_end;
                ++itr_set
              )
              {
                storage::Vector< FragmentComplete_pv> &rxt_vector( m_AvailableReactants[ itr_map->first]);
                if( rxt_vector.IsEmpty())
                {
                  rxt_vector = storage::Vector< FragmentComplete_pv>( itr_map->first->GetNumberReactants());
                }
                rxt_vector( *itr_set).PushBack( FragmentComplete_p( ( *m_Molecules)( mol_idx)));
              }
            }
          }
        }
      }
      BCL_MessageStd( "Done associating reactants with reactions");
    }

    //! @brief read RXN files from a directory to create reactions
    //! @param DIRECTORY the directory to use.
    //! @return a shptr to a vector of ReactionCompletes built from RXN files
    util::ShPtr< storage::Vector< ReactionComplete> > ReactionSearch::ReadReactions( const std::string &DIRECTORY)
    {
      util::ShPtr< storage::Vector< ReactionComplete> > reactions;
      io::Directory rxn_dir( DIRECTORY);
      if( !rxn_dir.DoesExist())
      {
        BCL_MessageCrt( "Warning: reaction directory \"" + DIRECTORY + "\" does not exist");
        return reactions;
      }
      
      storage::List< io::DirectoryEntry> files( rxn_dir.ListEntries());
      BCL_MessageStd( "Reading from directory \"" + DIRECTORY + "\"");
      io::IFStream input;
      for
      (
        storage::List< io::DirectoryEntry>::const_iterator itr_entry( files.Begin()), itr_entry_end( files.End());
        itr_entry != itr_entry_end;
        ++itr_entry
      )
      {
        if( !itr_entry->DoesExist() || itr_entry->GetType() != io::Directory::e_File)
        {
          BCL_MessageStd( "\"" + itr_entry->GetName() + "\" is not a file, skipping");
          continue;
        }
        else
        {
          BCL_MessageStd( "  Reading \"" + itr_entry->GetName() + "\"");
        }
        io::File::MustOpenIFStream( input, itr_entry->GetFullName());

        // TODO modularize so sdf namespace isn't used directly
        sdf::RXNHandler rxn_handler( input);
        if( rxn_handler.IsValid())
        {
          if( !reactions.IsDefined())
          {
            reactions = util::ShPtr< storage::Vector< ReactionComplete> >( new storage::Vector< ReactionComplete>());
          }
          reactions->PushBack( sdf::RXNFactory::MakeReactionComplete( rxn_handler));
        }
        else
        {
          BCL_MessageStd( "File \"" + itr_entry->GetName() + "\" did not contain RXN data");
        }
        io::File::CloseClearFStream( input);
      }
      return reactions;
    }

    //! @brief read the reactions and reactants from files
    void ReactionSearch::ReadFiles( const bool &FORCE) const
    {
      io::IFStream input;

      // Read in reactions
      if( ( FORCE || !m_Reactions.IsDefined() || m_Reactions->IsEmpty()) && !m_ReactionDirectory.empty())
      {
        io::Directory rxn_dir( m_ReactionDirectory);
        if( !rxn_dir.DoesExist())
        {
          BCL_MessageStd( "Reaction directory \"" + m_ReactionDirectory + "\" does not exist");
          m_Initialized = true;
          m_Reactions.Reset();
        }
        else
        {
          //BCL_Assert
          //( 
          //  !m_ReactionDirectory.empty(),
          //  "No reactions were provided, and no file containing reactions was given for the reaction catalog."
          //);
          // mark the class as dirty since we are reading in new data
          m_Initialized = false;

          BCL_MessageStd( "Reading reactions from " + util::Format()( m_ReactionDirectory));
          m_Reactions = util::ShPtr< storage::Vector< ReactionComplete> >( ReadReactions( m_ReactionDirectory));
          BCL_MessageStd( "  done reading reactions, read " + util::Format()( m_Reactions.IsDefined() ? m_Reactions->GetSize() : 0));
        }
      }

      // Read in molecules
      if( ( FORCE || !m_Molecules.IsDefined() || m_Molecules->IsEmpty()) && !m_ReactantFilename.empty())
      {
        io::DirectoryEntry rxt_file( m_ReactantFilename);
        if( !rxt_file.DoesExist() || rxt_file.GetType() != io::Directory::e_File)
        {
          BCL_MessageStd( "Reactant file \"" + m_ReactantFilename + "\" cannot be opened for reading");
          m_Initialized = true;
          m_Molecules.Reset();
        }
        else
        {
          // mark the class as dirty since we are reading in new data
          m_Initialized = false;

          BCL_MessageStd( "Reading reactants from " + util::Format()( m_ReactantFilename));
          io::File::MustOpenIFStream( input, m_ReactantFilename);
          FragmentEnsemble new_mols( input, sdf::e_Saturate);
          m_Molecules = util::ShPtr< storage::Vector< FragmentComplete> >( new storage::Vector< FragmentComplete>( new_mols.Begin(), new_mols.End()));
          io::File::CloseClearFStream( input);
          BCL_MessageStd( "  done reading reactants, read " + util::Format()( m_Molecules->GetSize()));
        }
      }
    }

    //! @brief retrieve the available reactant map
    //! @return the available reactant map
    const storage::Map
    < 
      ReactionComplete_p, 
      storage::Vector< FragmentComplete_pv>
    > &ReactionSearch::GetAvailableReactantsMap() const
    {
      return m_AvailableReactants;
    }

    ReactionStructure_pv ReactionSearch::FindMatchingStructures
    ( 
      const ConformationInterface &QUERY
    ) const
    {
      ReactionStructure_pv found_structs;

      size_t n_nodes( m_ReactionStructureTree.GetSize());

      storage::Vector< size_t> seen_vertices( n_nodes, size_t( 0));
      
      // Begin a search at each of the smallest substructures
      for( size_t i( 0), end_i( m_StartSearchVertices.GetSize()); i < end_i; ++i)
      {
        const size_t &start_index( m_StartSearchVertices( i));
        
        FindMatchingStructuresRecurse
        (
          QUERY,
          start_index,
          seen_vertices,
          found_structs
        );
      }
      return found_structs;
    }

    void ReactionSearch::FindMatchingStructuresRecurse
    ( 
      const ConformationInterface &QUERY,
      const size_t &CURRENT_NODE,
      storage::Vector< size_t> &CHECKED_NODES,
      ReactionStructure_pv &MATCHING_STRUCTS
    ) const
    {
      if( !CHECKED_NODES( CURRENT_NODE))
      {
        CHECKED_NODES( CURRENT_NODE) = 1;
        const util::SiPtr< const ReactionStructure> &node_info( m_ReactionStructureTree.GetVertexData( CURRENT_NODE));
        BCL_MessageDbg( "    Marking node " + util::Format()( CURRENT_NODE) + " as checked");
        if( node_info->ContainedIn( QUERY))
        {
          MATCHING_STRUCTS.PushBack( node_info);

          // do a depth-first search
          const storage::Vector< size_t> &neighbors( m_ReactionStructureTree.GetNeighborIndices( CURRENT_NODE));
          for( size_t n( 0), end_n( neighbors.GetSize()); n < end_n; ++n)
          {
            if( m_ReactionStructureTree.GetEdgeData( CURRENT_NODE, neighbors( n)) > 0)
            {
              FindMatchingStructuresRecurse( QUERY, neighbors( n), CHECKED_NODES, MATCHING_STRUCTS);
            }
          }
        }
      }
    }

  ///////////////
  // operators //
  ///////////////

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ReactionSearch::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT number of indentations
    //! @return ostream which was written to
    std::ostream &ReactionSearch::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////
    
    // create a single instance of this class
    const util::SiPtr< const util::ObjectInterface> ReactionSearch::s_Instance
    (
      GetObjectInstances().AddInstance( new ReactionSearch())
    );
  } // namespace chemistry
} // namespace bcl
