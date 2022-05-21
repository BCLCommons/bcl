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
#include "chemistry/bcl_chemistry_search_fragment_library_from_tree.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_bond_angle_assignment.h"
#include "chemistry/bcl_chemistry_rotamer_library_interface.h"
#include "sched/bcl_sched_thunk_job.h"

// external includes - sorted alphabetically

// Uncomment the next line to perform profiling of the tree search
//#define BCL_PROFILE_SearchFragmentLibraryFromTree

namespace bcl
{
  namespace chemistry
  {

    //! @brief default constructor
    SearchFragmentLibraryFromTree::SearchFragmentLibraryFromTree()
    {
    }

    //! @brief constructor
    //! @param FRAGMENT_LIBRARY pointer to ensemble of fragments
    SearchFragmentLibraryFromTree::SearchFragmentLibraryFromTree
    (
      const RotamerLibraryInterface &FRAGMENT_LIBRARY
    ) :
      m_FragmentLibrary( FRAGMENT_LIBRARY)
    {
      // load rotamer library information, except during static initialization
      if( m_FragmentLibrary->IsDefined())
      {
        LoadRotamerLibraryInformation();
      }
    }

    //! @brief Clone function
    //! @return pointer to new SearchFragmentLibraryFromTree
    SearchFragmentLibraryFromTree *SearchFragmentLibraryFromTree::Clone() const
    {
      return new SearchFragmentLibraryFromTree( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SearchFragmentLibraryFromTree::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief load the rotamer library information. It is necessary to call this function before the first call to
    //!        any of the operations. It should not, however, be called in the constructor, because this class may be
    //!        constructed with a RotamerLibraryInterface during static initialization, and where that happens, the
    //!        rotamer library information does not exist
    void SearchFragmentLibraryFromTree::LoadRotamerLibraryInformation() const
    {
      // bool, to circumvent mutex. This bool will only be set once all rotamer library information has been read
      static sched::Mutex s_mutex;
      s_mutex.Lock();
      if( m_SubstructureTree.GetSize())
      {
        s_mutex.Unlock();
        return;
      }
      // get the substructure tree of constiutions
      m_SubstructureTree = m_FragmentLibrary->RetrieveRotamerSubstructureTree();

      // set size of m_ConstitutionGraphs and m_ConstitutionsSeen
      size_t number_constitutions( m_SubstructureTree.GetSize());
      m_ConstitutionGraphs = storage::Vector< graph::ConstGraph< size_t, size_t> >( number_constitutions);
      m_ConstitutionsSeen = std::string( number_constitutions, '0');
      m_Muteces.Resize( number_constitutions);

      const storage::Vector< graph::ConstGraph< size_t, size_t> > root_constitutions_graph( m_FragmentLibrary->GetRootConstitutions());

      size_t node_index( 0);

      for
      (
        storage::Vector< graph::ConstGraph< size_t, size_t> >::const_iterator
          itr( root_constitutions_graph.Begin()), itr_end( root_constitutions_graph.End());
          itr != itr_end;
        ++itr, ++node_index
      )
      {
        m_ConstitutionGraphs( node_index) = *itr;
        m_ConstitutionsSeen[ node_index] = '1';
      }

      linal::Vector< size_t> root_nodes( linal::FillVector( root_constitutions_graph.GetSize(), size_t( 0), size_t( 1)));
      m_RootNodes = storage::List< size_t>( root_nodes.Begin(), root_nodes.End());

      m_ConfigurationMapping = m_FragmentLibrary->RetrieveConfigurationMapping();

      s_mutex.Unlock();
    }

  ////////////////
  // operations //
  ////////////////

    namespace
    {
      //! @brief fast check for whether a molecule has any dihedral edges in chains
      bool HasChainDihedralEdges( const ConformationInterface &MOL)
      {
        for( auto itr( MOL.GetAtomsIterator()); itr.NotAtEnd(); ++itr)
        {
          if( itr->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 0)) >= 2)
          {
            return itr.GetSize() > size_t( 3);
          }
        }
        return false;
      }
    }

    //! @brief find fragments that have isomorphism with molecule
    //! @param MOLECULE molecule of interest whose fragments need to be searched from the given library of fragments
    //! @return vector of SmallMoleculeFragmentIsomorphism
    util::ShPtrVector< SmallMoleculeFragmentIsomorphism> SearchFragmentLibraryFromTree::FindFragmentsOfMolecule
    (
      const ConformationInterface &MOLECULE
    ) const
    {
      #ifdef BCL_PROFILE_SearchFragmentLibraryFromTree
      static util::Stopwatch s_load_rot( "LoadRotamerLibraryInformation", util::Time( 1, 0), util::Message::e_Standard, true, false);
      s_load_rot.Start();
      #endif

      // load the rotamer library information, if necessary
      LoadRotamerLibraryInformation();

      #ifdef BCL_PROFILE_SearchFragmentLibraryFromTree
      s_load_rot.Stop();
      #endif
      // graph maker to create molecule graph
      ConformationGraphConverter graph_maker( ConformationGraphConverter::e_AtomType, ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness);
      graph::ConstGraph< size_t, size_t> molecule_graph( graph_maker( MOLECULE));
      ConformationGraphConverter graph_maker_nh( ConformationGraphConverter::e_AtomTypeAndDistinguishHydrogens, ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness);
      graph::ConstGraph< size_t, size_t> molecule_graph_nh( graph_maker_nh( MOLECULE));

      #ifdef BCL_PROFILE_SearchFragmentLibraryFromTree
      static util::Stopwatch s_iso( "Search fragments", util::Time( 1, 0), util::Message::e_Standard, true, false);
      s_iso.Start();
      #endif

      storage::Vector< util::ShPtrList< SmallMoleculeFragmentIsomorphism> > isomorphism_object( MOLECULE.GetNumberAtoms() + 1);
      storage::List< size_t> current_nodes( m_RootNodes);
      storage::List< size_t> retrieve_configurations;

      #ifdef BCL_PROFILE_SearchFragmentLibraryFromTree
      static util::Stopwatch s_vertcon( "GetVerticesConnectedFrom", util::Time( 1, 0), util::Message::e_Standard, true, false);
      static util::Stopwatch s_check( "CheckSubStructure", util::Time( 1, 0), util::Message::e_Standard, true, false);
      static util::Stopwatch s_root( "Checking root nodes", util::Time( 1, 0), util::Message::e_Standard, true, false);
      #endif

      std::string have_seen( m_ConstitutionGraphs.GetSize(), '0');
      #ifdef BCL_PROFILE_SearchFragmentLibraryFromTree
      size_t current_itr( 0);
      #endif
      while( !current_nodes.IsEmpty())
      {
        storage::List< size_t> new_cur_nodes;
        #ifdef BCL_PROFILE_SearchFragmentLibraryFromTree
        if( !current_itr)
        {
          s_root.Start();
          ++current_itr;
        }
        #endif
        for
        (
          storage::List< size_t>::const_iterator
            itr_cur( current_nodes.Begin()), itr_cur_end( current_nodes.End());
          itr_cur != itr_cur_end;
          ++itr_cur
        )
        {
          if( have_seen[ *itr_cur] == '0')
          {
            have_seen[ *itr_cur] = '1';
            m_Muteces( *itr_cur).Lock();
            if( m_ConstitutionsSeen[ *itr_cur] == '0')
            {
              m_ConstitutionGraphs( *itr_cur) =
                m_FragmentLibrary->RetrieveConstitution( storage::Vector< size_t>::Create( *itr_cur))( 0);
              m_ConstitutionsSeen[ *itr_cur] = '1';
            }
            m_Muteces( *itr_cur).Unlock();

            #ifdef BCL_PROFILE_SearchFragmentLibraryFromTree
            s_check.Start();
            #endif
            if( CheckSubStructure( molecule_graph, m_ConstitutionGraphs( *itr_cur)))
            {
              #ifdef BCL_PROFILE_SearchFragmentLibraryFromTree
              s_check.Stop();
              #endif
              const storage::Vector< size_t> &neighboring_nodes( m_SubstructureTree.GetNeighborIndices( *itr_cur));
              new_cur_nodes.InsertElements( new_cur_nodes.End(), neighboring_nodes.Begin(), neighboring_nodes.End());
              retrieve_configurations.PushBack( *itr_cur);
            }
            else
            {
              have_seen[ *itr_cur] = 'X';
              #ifdef BCL_PROFILE_SearchFragmentLibraryFromTree
              s_check.Stop();
              s_vertcon.Start();
              #endif

              // initialize seen_vertices with the first vertex
              storage::Vector< size_t> seen_vertices_queue( 1, *itr_cur);

              size_t vertex_queue_position( 0); // index of the active vertex in the breadth-first-search

              // So long as there are vertices left in the queue whose connections haven't been examined, this loop will continue
              // unless all vertices are put into the graph.
              for( ; vertex_queue_position < seen_vertices_queue.GetSize(); ++vertex_queue_position)
              {
                const size_t current_vertex( seen_vertices_queue( vertex_queue_position));
                // target row is a reference to the edges reachable from the current vertex
                const storage::Vector< size_t> &target_row( m_SubstructureTree.GetNeighborIndices( current_vertex));
                for( size_t i( 0), number_seen( target_row.GetSize()); i < number_seen; ++i)
                {
                  const size_t new_vertex( target_row( i));
                  if( have_seen[ new_vertex] == '0') // found a vertex in the target list of vertex seen_vertices_queue(vertex_queue_position)
                  {
                    seen_vertices_queue.PushBack( new_vertex);
                    have_seen[ new_vertex] = 'X';
                  }
                }
              }
              #ifdef BCL_PROFILE_SearchFragmentLibraryFromTree
              s_vertcon.Stop();
              #endif
            }
          }
        }
        #ifdef BCL_PROFILE_SearchFragmentLibraryFromTree
        if( current_itr == size_t( 1))
        {
          s_root.Stop();
          ++current_itr;
        }
        #endif
        current_nodes = new_cur_nodes;
      }

      storage::List< size_t> configuration_indices;
      for
      (
        storage::List< size_t>::const_iterator
          itr( retrieve_configurations.Begin()), itr_end( retrieve_configurations.End());
        itr != itr_end;
        ++itr
      )
      {
        const storage::Set< size_t> &configurations( m_ConfigurationMapping( *itr));
        configuration_indices.PushBack( *configurations.Begin());
      }
      #ifdef BCL_PROFILE_SearchFragmentLibraryFromTree
      s_iso.Stop();
      static util::Stopwatch s_load( "Load configurations", util::Time( 1, 0), util::Message::e_Standard, true, false);
      s_load.Start();
      #endif

      FragmentEnsemble ensemble( m_FragmentLibrary->RetrieveAssociatedConfigurations( storage::Set< size_t>( configuration_indices.Begin(), configuration_indices.End())));
      #ifdef BCL_PROFILE_SearchFragmentLibraryFromTree
      s_load.Stop();
      static util::Stopwatch s_create( "CreateFragmentIsomorphismObject", util::Time( 1, 0), util::Message::e_Standard, true, false);
      s_create.Start();
      #endif
      for
      (
        storage::List< FragmentComplete>::const_iterator
          itr_configurations( ensemble.Begin()), itr_configurations_end( ensemble.End());
        itr_configurations != itr_configurations_end;
        ++itr_configurations
      )
      {
        if( !itr_configurations->GetNumberAtoms())
        {
          FragmentComplete molecule( CreateMoleculeFromRotlib( *itr_configurations));
          molecule.StoreProperties( *itr_configurations);
          bool has_chain( HasChainDihedralEdges( molecule));
          if( itr_configurations->GetNumberAtoms() > size_t( 4) && has_chain)
          {
            isomorphism_object( molecule.GetNumberAtoms()).PushBack
            (
              CreateFragmentIsomorphismObject( MOLECULE, molecule_graph_nh, graph_maker_nh( molecule), molecule)
            );
          }
          else
          {
            isomorphism_object( molecule.GetNumberAtoms()).PushBack
            (
              CreateFragmentIsomorphismObject( MOLECULE, molecule_graph, graph_maker( molecule), molecule)
            );
          }
        }
        else
        {
          bool has_chain( HasChainDihedralEdges( *itr_configurations));
          if( itr_configurations->GetNumberAtoms() > size_t( 4) && has_chain)
          {
            isomorphism_object( itr_configurations->GetNumberAtoms()).PushBack
                (
                  CreateFragmentIsomorphismObject( MOLECULE, molecule_graph_nh, graph_maker_nh( *itr_configurations), *itr_configurations)
                );
          }
          else
          {
            isomorphism_object( itr_configurations->GetNumberAtoms()).PushBack
                (
                  CreateFragmentIsomorphismObject( MOLECULE, molecule_graph, graph_maker( *itr_configurations), *itr_configurations)
                );
          }
          if( !isomorphism_object( itr_configurations->GetNumberAtoms()).LastElement().IsDefined())
          {
            isomorphism_object( itr_configurations->GetNumberAtoms()).PopBack();
          }
        }
      }
      // create a vector of SmallMoleculeFragmentIsomorphism objects for each fragment which will be returned
      util::ShPtrVector< SmallMoleculeFragmentIsomorphism> fragments;
      fragments.AllocateMemory( ensemble.GetSize());
      for( auto itr( isomorphism_object.ReverseBegin()), itr_end( isomorphism_object.ReverseEnd()); itr != itr_end; ++itr)
      {
        for( auto itr_iso( itr->Begin()), itr_iso_end( itr->End()); itr_iso != itr_iso_end; ++itr_iso)
        {
          fragments.PushBack( *itr_iso);
        }
      }
      #ifdef BCL_PROFILE_SearchFragmentLibraryFromTree
      s_create.Stop();
      #endif
      return fragments;
    }

    //! @brief checks whether fragment graph is part of molecule graph
    //! @param MOLECULE_GRAPH molecule graph for which fragments need to be checked
    //! @param FRAGMENT_GRAPH fragment graph which needs to be checked if contained in molecule
    bool SearchFragmentLibraryFromTree::CheckSubStructure
    (
      const graph::ConstGraph< size_t, size_t> &MOLECULE_GRAPH,
      const graph::ConstGraph< size_t, size_t> &FRAGMENT_GRAPH
    )
    {
      if( FRAGMENT_GRAPH.GetSize() < 2)
      {
        return false;
      }
      // check that the # of vertices and edges are the same
      if
      (
        MOLECULE_GRAPH.GetSize() < FRAGMENT_GRAPH.GetSize()
        || MOLECULE_GRAPH.NumEdges() < FRAGMENT_GRAPH.NumEdges()
        || !graph::CSISubstructure::IsContainedIn( FRAGMENT_GRAPH.GetVertexTypeCountMap(), MOLECULE_GRAPH.GetVertexTypeCountMap())
        || !graph::CSISubstructure::IsContainedIn( FRAGMENT_GRAPH.GetEdgeTypeCountMap(), MOLECULE_GRAPH.GetEdgeTypeCountMap())
      )
      {
        return false;
      }

      storage::List< graph::SubgraphIsomorphism< size_t, size_t> >::iterator simple_csi_itr( AcquireIsomorphism());
      // object that will compute the fragment isomorphism
      graph::SubgraphIsomorphism< size_t, size_t> &simple_csi( *simple_csi_itr);

      // set the smaller graph of the substructure to the molecule graph
      simple_csi.SetGraphExternalOwnership( MOLECULE_GRAPH);

      // set the smaller graph to
      simple_csi.SetSubgraphExternalOwnership( FRAGMENT_GRAPH);

      // there enough edges and vertices of the right type in the molecule's graph
      // only select isomorphisms with the  whole fragment and insert in map for storing SmallMoleculeFragmentIsomorphism
      const bool has_isomorphism( simple_csi.FindIsomorphism());
      ReleaseIsomorphism( simple_csi_itr);
      return has_isomorphism;
    }

    namespace
    {
      //! @brief get the vertices involved in dihedral edges
      storage::Vector< size_t> GetVerticesInDihedralEdges( const FragmentComplete &FRAGMENT)
      {
        auto dihedral_angle_edges( PriorityDihedralAngles()( FRAGMENT).Second());
        storage::Vector< size_t> is_in_dihedral_edge( FRAGMENT.GetSize(), size_t( 0));
        for( auto itr( dihedral_angle_edges.Begin()), itr_end( dihedral_angle_edges.End()); itr != itr_end; ++itr)
        {
          is_in_dihedral_edge( itr->First()) = size_t( 1);
          is_in_dihedral_edge( itr->Second()) = size_t( 1);
          is_in_dihedral_edge( itr->Third()) = size_t( 1);
          is_in_dihedral_edge( itr->Fourth()) = size_t( 1);
        }
        storage::Vector< size_t> dihedral_indices;
        dihedral_indices.AllocateMemory( FRAGMENT.GetSize());
        for( size_t i( 0), sz( FRAGMENT.GetSize()); i < sz; ++i)
        {
          if( is_in_dihedral_edge( i))
          {
            dihedral_indices.PushBack( i);
          }
        }
        return dihedral_indices;
      }
    }

    //! @brief creates smallmolecule fragment isomorphism object
    //! @param MOLECULE_GRAPH molecule graph for which fragments need to be checked
    //! @param FRAGMENT_GRAPH fragment graph which needs to be checked if contained in molecule
    //! @param FRAGMENT fragment which is part of molecule
    //! @return return a small molecule fragment isomorphism object
    util::ShPtr< SmallMoleculeFragmentIsomorphism> SearchFragmentLibraryFromTree::CreateFragmentIsomorphismObject
    (
      const FragmentComplete &MOLECULE,
      const graph::ConstGraph< size_t, size_t> &MOLECULE_GRAPH,
      const graph::ConstGraph< size_t, size_t> &FRAGMENT_GRAPH,
      const FragmentComplete &FRAGMENT
    )
    {
      // get the vertices that matter for dihedral angles. For symmetric positions (mainly hydrogens)
      storage::Vector< size_t> dihedral_indices( GetVerticesInDihedralEdges( FRAGMENT));
      storage::List< graph::SubgraphIsomorphism< size_t, size_t> >::iterator simple_csi_itr( AcquireIsomorphism());
      // object that will compute the fragment isomorphism
      graph::SubgraphIsomorphism< size_t, size_t> &simple_csi( *simple_csi_itr);

      // set the smaller graph of the substructure to the molecule graph
      simple_csi.SetGraphExternalOwnership( MOLECULE_GRAPH);

      // set the smaller graph to
      simple_csi.SetSubgraphExternalOwnership( FRAGMENT_GRAPH);

      // there enough edges and vertices of the right type in the molecule's graph
      // only select isomorphisms with the  whole fragment and insert in map for storing SmallMoleculeFragmentIsomorphism
      bool found( simple_csi.FindAllIsomorphisms());
      simple_csi.PruneIsomorphismsToThoseUniqueInField( dihedral_indices);
      util::ShPtr< SmallMoleculeFragmentIsomorphism> iso
      (
        new SmallMoleculeFragmentIsomorphism( MOLECULE, FRAGMENT, simple_csi)
      );
      ReleaseIsomorphism( simple_csi_itr);
      if( !found)
      {
        return util::ShPtr< SmallMoleculeFragmentIsomorphism>();
      }
      return iso;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SearchFragmentLibraryFromTree::Read( std::istream &ISTREAM)
    {
      // read members
      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SearchFragmentLibraryFromTree::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief acquire mutex for acquiring an isomorphism object
    //! @return mutex used in AcquireIsomorphism and ReleaseIsomorphism
    sched::Mutex &SearchFragmentLibraryFromTree::GetIsomorphismMutex()
    {
      static sched::Mutex s_mutex;
      return s_mutex;
    }

    //! @brief acquire an isomorphism object
    //! @return iterator to the isomorphism object
    storage::List< graph::SubgraphIsomorphism< size_t, size_t> >::iterator SearchFragmentLibraryFromTree::AcquireIsomorphism()
    {
      GetIsomorphismMutex().Lock();

      storage::List< graph::SubgraphIsomorphism< size_t, size_t> > &avail( GetAvailableIsomorphisms());
      storage::List< graph::SubgraphIsomorphism< size_t, size_t> > &alloc( GetAllocatedIsomorphisms());
      // all isomorphisms have been allocated
      if( avail.IsEmpty())
      {
        avail.PushBack( graph::SubgraphIsomorphism< size_t, size_t>());
      }
      // splice the last hidden vector set from the available list onto the allocated list
      storage::List< graph::SubgraphIsomorphism< size_t, size_t> >::iterator avail_itr( avail.Begin());
      alloc.InternalData().splice( alloc.Begin(), avail.InternalData(), avail_itr);

      // save the iterator to the now-allocated array set
      avail_itr = alloc.Begin();
      GetIsomorphismMutex().Unlock();
      return avail_itr;
    }

    //! @brief release a given isomorphism object
    //! @param ITR iterator to the isomorphism object
    void SearchFragmentLibraryFromTree::ReleaseIsomorphism
    (
      const storage::List< graph::SubgraphIsomorphism< size_t, size_t> >::iterator &ITR
    )
    {
      GetIsomorphismMutex().Lock();
      // splice the iterator from the allocated back onto available
      GetAvailableIsomorphisms().InternalData().splice( GetAvailableIsomorphisms().Begin(), GetAllocatedIsomorphisms().InternalData(), ITR);
      GetIsomorphismMutex().Unlock();
    }

    //! @brief get the list of isomorphisms that have been allocated
    storage::List< graph::SubgraphIsomorphism< size_t, size_t> > &SearchFragmentLibraryFromTree::GetAllocatedIsomorphisms()
    {
      static storage::List< graph::SubgraphIsomorphism< size_t, size_t> > s_allocated;
      return s_allocated;
    }

    //! @brief get the list of isomorphisms that are available
    storage::List< graph::SubgraphIsomorphism< size_t, size_t> > &SearchFragmentLibraryFromTree::GetAvailableIsomorphisms()
    {
      static storage::List< graph::SubgraphIsomorphism< size_t, size_t> > s_available;
      return s_available;
    }

    //! @brief create molecule from empty molecule containing atom vector and bond vector information in properties
    //! @param molecule which contains atomvector and bondvector information in properties
    //! @return molecule with all the atom type, bond type and connectivity information
    FragmentComplete SearchFragmentLibraryFromTree::CreateMoleculeFromRotlib
    (
      FragmentComplete MOLECULE
    )
    {
      linal::Vector< double> atoms_string( MOLECULE.GetStoredProperties().GetMDLPropertyAsVector( "AtomVector"));
      size_t number_atoms( atoms_string.GetSize());
      storage::Vector< sdf::AtomInfo> atom_info;
      atom_info.AllocateMemory( number_atoms);

      linal::Vector< double> bond_string( MOLECULE.GetStoredProperties().GetMDLPropertyAsVector( "BondVector"));
      size_t number_bonds( bond_string.GetSize() / size_t( 3));
      storage::Vector< sdf::BondInfo> bond_info;
      bond_info.AllocateMemory( number_bonds);
      for( size_t entry_index( 0); entry_index < number_atoms; entry_index++)
      {
        size_t string_position( entry_index);
        atom_info.PushBack( sdf::AtomInfo( AtomType( atoms_string( string_position)), e_UnknownChirality));
      }
      for( size_t entry_index( 0); entry_index < number_bonds; entry_index++)
      {
        size_t string_position( entry_index * 3);
        bond_info.PushBack( sdf::BondInfo( bond_string( string_position + 0), bond_string( string_position + 1), ConfigurationalBondType( bond_string( string_position + 2))));
      }
      return FragmentComplete( AtomVector< AtomComplete>( atom_info, bond_info), MOLECULE.GetName());
    }

    //! @brief search bond angle library for a given atom
    typename RotamerLibraryInterface::t_BondAngleMap::const_iterator
    SearchFragmentLibraryFromTree::SearchForBondAngles
    (
      const AtomConformationalInterface &ATOM,
      const ConformationGraphConverter::AtomComparisonType &ATOM_COMPARISON
    ) const
    {
      typename RotamerLibraryInterface::t_BondAngleMapKey key;
      key.First() = ATOM_COMPARISON;
      key.Second() = ATOM.GetAtomType().GetIndex();

      for
      (
        auto itr_bond( ATOM.GetBonds().Begin()), itr_bond_end( ATOM.GetBonds().End());
        itr_bond != itr_bond_end;
        ++itr_bond
      )
      {
        if( itr_bond->GetTargetAtom().GetAtomType() == GetAtomTypes().H_S)
        {
          continue;
        }
        key.Third().PushBack();
        key.Third().LastElement().First() =
          itr_bond->GetBondType()->GetBondData
          (
            ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness
          );
        key.Third().LastElement().Second() = ConformationGraphConverter::ConvertAtomTypeData
                                             (
                                               itr_bond->GetTargetAtom().GetAtomType(),
                                               ATOM_COMPARISON
                                             );
      }
      key.Third().Sort( std::greater< storage::Pair< size_t, size_t> >());

      return m_BondAngleMap.Find( key);
    }

    //! @brief find bond angle information objects for a given molecule
    //! @param MOLECULE molecule of interest whose fragments need to be searched from the given library of fragments
    //! @return vector of SmallMoleculeFragmentIsomorphism
    util::ShPtrVector< BondAngleAssignment> SearchFragmentLibraryFromTree::GetBondAngleAssignments
    (
      const ConformationInterface &MOLECULE,
      const bool &ENFORCE_CHIRALITY
    ) const
    {
      static sched::Mutex s_mutex;
      s_mutex.Lock();
      if( m_BondAngleMap.IsEmpty())
      {
        m_BondAngleMap = m_FragmentLibrary->GetBondAngleMap();
      }
      s_mutex.Unlock();
      util::ShPtrVector< BondAngleAssignment> bas;
      bas.AllocateMemory( MOLECULE.GetSize());

      // Key <- atom-type <> Vector< bond type, atom type>
      // Value <- Matrix with unit-vector coordinates of all atoms after the first.
      // The first atom is always moved to 1 0 0, second atom is moved such that it is 0 in
      for( iterate::Generic< const AtomConformationalInterface> itr_atm( MOLECULE.GetAtomsIterator()); itr_atm.NotAtEnd(); ++itr_atm)
      {
        if
        (
          itr_atm->GetBonds().GetSize() <= size_t( 1)
          || itr_atm->GetBonds().GetSize() > size_t( 4)
          || !itr_atm->GetAtomType()->GetHybridOrbitalType().IsDefined()
          || itr_atm->GetAtomType()->GetHybridOrbitalType() == GetHybridOrbitalTypes().e_Unhybridized
        )
        {
          continue;
        }
        const size_t n_ring( itr_atm->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1)));
        size_t n_effbonds( itr_atm->GetBonds().GetSize() - itr_atm->GetNumberCovalentlyBoundHydrogens());
        if( n_ring == n_effbonds || n_effbonds <= size_t( 1))
        {
          continue;
        }
        // during testing, we will limit this to the simple case of two bonds
        storage::Vector< ConformationGraphConverter::AtomComparisonTypeEnum> searches
        (
          storage::Vector< ConformationGraphConverter::AtomComparisonTypeEnum>::Create
          (
            ConformationGraphConverter::e_AtomType,
            ConformationGraphConverter::e_ElementType,
            ConformationGraphConverter::e_Identity
          )
        );
        bool found( false);
        for( auto itr_searches( searches.Begin()), itr_searches_end( searches.End()); itr_searches != itr_searches_end; ++itr_searches)
        {
          auto itr_map( SearchForBondAngles( *itr_atm, *itr_searches));
          if( itr_map != m_BondAngleMap.End())
          {
            bas.PushBack
            (
              util::ShPtr< BondAngleAssignment>
              (
                new BondAngleAssignment
                (
                  itr_map->second,
                  itr_map->first,
                  MOLECULE,
                  itr_atm.GetPosition(),
                  !ENFORCE_CHIRALITY,
                  *itr_searches
                )
              )
            );
            found = true;
            break;
          }
        }
        if( !found)
        {
          std::ostringstream outs;
          outs << itr_atm->GetAtomType().GetName() << ' ';
          for
          (
            auto itr_bond( itr_atm->GetBonds().Begin()), itr_bond_end( itr_atm->GetBonds().End());
            itr_bond != itr_bond_end;
            ++itr_bond
          )
          {
            if( itr_bond->GetTargetAtom().GetAtomType() == GetAtomTypes().H_S)
            {
              continue;
            }
            outs
              << itr_bond->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness)
              << ' ';
          }
          BCL_MessageStd( "Could not find a bond angle set for " + outs.str() + " even ignoring neighboring atom types!");
        }
      }
      return bas;
    }

  } // namespace chemistry

} // namespace bcl
