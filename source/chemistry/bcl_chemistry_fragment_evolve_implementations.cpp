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
#include "chemistry/bcl_chemistry_fragment_evolve_implementations.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_collector_valence.h"
#include "chemistry/bcl_chemistry_conformation_set.h"
#include "chemistry/bcl_chemistry_fragment_grow.h"
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"
#include "chemistry/bcl_chemistry_pick_atom_random.h"
#include "chemistry/bcl_chemistry_pick_fragment_random.h"
#include "chemistry/bcl_chemistry_stereocenters_handler.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "io/bcl_io_file.h"
#include "sdf/bcl_sdf_fragment_factory.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    //! @param EVOLVE_TYPE the evolution type to use (defines operator() behavior)
    //! @param CHECK_BAD_BONDS whether to ensure no bad bonds a formed
    FragmentEvolveImplementations::FragmentEvolveImplementations( const EvolveType &EVOLVE_TYPE, const bool &ALLOW_BAD_BONDS) :
      util::SerializableInterface(),
      m_EvolveType( EVOLVE_TYPE),
      m_FragmentPool(),
      m_AllowBadBonds( ALLOW_BAD_BONDS)
    {
    }

    //! @brief Clone function
    //! @return pointer to new FragmentEvolveImplementations
    FragmentEvolveImplementations *FragmentEvolveImplementations::Clone() const
    {
      return new FragmentEvolveImplementations( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &FragmentEvolveImplementations::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns class name of the object when used in a dynamic context
    //! @return the class name
    const std::string &FragmentEvolveImplementations::GetAlias() const
    {
      return GetEvolveTypeString( m_EvolveType);
    }

    //! @brief gets the type of evolution operations that will be performed
    //! @return the evolution type
    const FragmentEvolveImplementations::EvolveType &FragmentEvolveImplementations::GetEvolveType() const
    {
      return m_EvolveType;
    }

    //! @brief gets the string associated with each evolve type
    //! @param EVOLVE_TYPE the evolve type to get the string for
    //! @return a string describing the evolve type
    const std::string &FragmentEvolveImplementations::GetEvolveTypeString( const FragmentEvolveImplementations::EvolveType &EVOLVE_TYPE)
    {
      static const std::string s_clone_alias( "Clone");
      static const std::string s_add_alias( "AddFragment");
      static const std::string s_del_alias( "DelFragment");
      static const std::string s_combine_alias( "Combine");
      static const std::string s_any_alias( "Invalid");

      switch( EVOLVE_TYPE)
      {
        case e_Clone:
          return s_clone_alias;
          break;
        case e_FragAdd:
          return s_add_alias;
          break;
        case e_FragDel:
          return s_del_alias;
          break;
        case e_Combine:
          return s_combine_alias;
          break;
        default:
          return s_any_alias;
          break;
      }
      return s_any_alias;
    }

    //! @brief gets the fragment list that is in use
    //! @return a ShPtr to a FragmentEnsemble that holds the fragments
    const util::ShPtr< FragmentEnsemble> &FragmentEvolveImplementations::GetFragmentList() const
    {
      return m_FragmentPool;
    }

    //! @brief remove a fragment from a molecule
    //! @param MOLECULE the molecule to remove a fragment from
    //! @return a new molecule with a removed fragment
    void FragmentEvolveImplementations::SetFragmentList( const FragmentEnsemble &FRAGMENT_LIST)
    {
      if( FRAGMENT_LIST.GetSize() == 0)
      {
        BCL_MessageStd( "fragment list is empty.  not replacing current fragment list");
      }
      else
      {
        m_FragmentPool = util::CloneToShPtr< FragmentEnsemble>( FRAGMENT_LIST);
      }
    }

    //! @brief sets the fragment list for this class
    //! @param FRAGMENT_LIST the fragment list to use
    void FragmentEvolveImplementations::SetFragmentList( const util::ShPtr< FragmentEnsemble> &FRAGMENT_LIST)
    {
      m_FragmentPool = FRAGMENT_LIST;
    }

    //! @brief sets the fragment filename to a new file 
    //! @param FILENAME the file to read from
    void FragmentEvolveImplementations::SetFragmentList( const std::string &FILENAME)
    {
      m_FragmentFilename = FILENAME;
      m_FragmentPool.Reset();
//      BCL_Assert( ReadFragmentsFromFile( FILENAME), "Could not read fragments from file");
    }

    //! @brief sets whether to check for bad bonds
    //! @param ALLOW_BAD_BONDS if false, adding fragments will check for bad bonds
    void FragmentEvolveImplementations::SetAllowBadBonds( const bool &ALLOW_BAD_BONDS)
    {
      m_AllowBadBonds = ALLOW_BAD_BONDS;
    }

    //! @brief gets the number of molecules needed by operator()
    //! @return the number needed
    size_t FragmentEvolveImplementations::NumRequiredMols() const
    {
      if( m_EvolveType == e_Combine)
      {
        return 2;
      }
      return 1;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief run the operation specified by the class
    //! @param MEMBERS an ensemble of molecules to use
    //! @return a new ensemble of mutated molecules
    util::ShPtrVector< FragmentComplete> FragmentEvolveImplementations::MakeMolecules
    (
      const util::SiPtrVector< const FragmentComplete> &MEMBERS
    ) const
    {
      if( MEMBERS.GetSize() < NumRequiredMols())
      {
        return util::ShPtrVector< FragmentComplete>();
      }
      const FragmentComplete &first_mol( *MEMBERS( 0));
      
      BCL_MessageVrb( GetAlias());
      switch( m_EvolveType)
      {
        case e_Clone:
          {
            util::ShPtrVector< FragmentComplete> frag_vector;
            frag_vector.PushBack( util::ShPtr< FragmentComplete>( new FragmentComplete( first_mol)));
            return frag_vector;
          }
          break;
        case e_FragAdd:
          if( !m_FragmentPool.IsDefined())
          {
            BCL_Assert( ReadFragmentsFromFile(), "Could not read fragments from file");
          }
          return MutateAdd( first_mol, *m_FragmentPool);
          break;
        case e_FragDel:
          return MutateDel( first_mol);
          break;
        case e_Combine:
          {
            const FragmentComplete &second_mol( *MEMBERS( 1));
            return Combine( first_mol, second_mol);
            break;
          }
        case s_End:
          BCL_Assert( false, "invalid evolution type");
          break;
      }
      return util::ShPtrVector< FragmentComplete>();
    }

    //! @brief combines sections of two molecules
    //! @param FIRST_MOLECULE the first molecule
    //! @param SECOND_MOLECULE the second molecule
    //! @param MAX_FRAG_SIZE the maximum size of fragments that should be connected together
    //! @return a new fragment consisting of pieces of the inputs
    util::ShPtrVector< FragmentComplete> FragmentEvolveImplementations::Combine
    (
      const FragmentComplete &FIRST_MOLECULE,
      const FragmentComplete &SECOND_MOLECULE,
      const size_t MAX_FRAG_SIZE
    ) const
    {
      if( FIRST_MOLECULE.GetNumberAtoms() == 0 || SECOND_MOLECULE.GetNumberAtoms() == 0)
      {
        return util::ShPtrVector< FragmentComplete>();
      }

      // Remove hydrogens
      FragmentComplete first_molecule( FIRST_MOLECULE);
      first_molecule.RemoveH();
      FragmentComplete second_molecule( SECOND_MOLECULE);
      second_molecule.RemoveH();

      ConformationGraphConverter graph_maker
      (
        ConformationGraphConverter::e_AtomTypeAndChirality,
        ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness
      );

      // Make a size_t graph of the molecules
      graph::ConstGraph< size_t, size_t> first_mol_graph( graph_maker( first_molecule));
      graph::ConstGraph< size_t, size_t> second_mol_graph( graph_maker( second_molecule));
      
      // Break a random single bond in each molecule and retrieve the fragments
      storage::List< storage::Vector< size_t> > first_mol_frags
      (
        FragmentsFromRandomBondBreakage
        (
          first_molecule, first_mol_graph
        )
      );

      storage::List< storage::Vector< size_t> > second_mol_frags
      (
        FragmentsFromRandomBondBreakage
        (
          second_molecule, second_mol_graph
        )
      );
      
      size_t first_mol_frags_size( first_mol_frags.GetSize());
      size_t second_mol_frags_size( second_mol_frags.GetSize());

      // there should be only two fragments from each molecule, otherwise combining pieces is senseless
      if( first_mol_frags_size != 2 || second_mol_frags_size != 2)
      {
        BCL_MessageVrb
        (
          "First mol had " + util::Format()( first_mol_frags_size) + " fragments, second had "
          + util::Format()( second_mol_frags_size)
        );
        return util::ShPtrVector< FragmentComplete>();
      }

      // If everything succeeded, make atom graphs of the molecules
      ConformationGraphConverter::t_AtomGraph first_mol_atom_graph( graph_maker.CreateGraphWithAtoms( first_molecule));
      ConformationGraphConverter::t_AtomGraph second_mol_atom_graph( graph_maker.CreateGraphWithAtoms( second_molecule));

      // OwnPtrs to graphs; to be used in subgraph construction
      util::OwnPtr< graph::ConstGraph< size_t, size_t> > optr_first_graph( &first_mol_graph, false);
      util::OwnPtr< graph::ConstGraph< size_t, size_t> > optr_second_graph( &second_mol_graph, false);

      // The FragmentEnsemble that will be returned
      util::ShPtrVector< FragmentComplete> new_molecules;

      ///////////////
      // Determine which bond was broken in each molecule
      ///////////////

      graph::Subgraph< size_t, size_t> first_subgraph( optr_first_graph, *first_mol_frags.Begin());
      storage::List< storage::Pair< size_t, size_t> > first_adj_edges( first_subgraph.GetAdjacentEdgeIndices());
      if( first_adj_edges.GetSize() != 1)
      {
        return util::ShPtrVector< FragmentComplete>();
      }
      storage::Pair< size_t, size_t> &first_broken_bond( *first_adj_edges.Begin());

      graph::Subgraph< size_t, size_t> second_subgraph( optr_second_graph, *second_mol_frags.Begin());
      storage::List< storage::Pair< size_t, size_t> > second_adj_edges( second_subgraph.GetAdjacentEdgeIndices());
      if( second_adj_edges.GetSize() != 1)
      {
        return util::ShPtrVector< FragmentComplete>();
      }
      storage::Pair< size_t, size_t> &second_broken_bond( *second_adj_edges.Begin());

      // Combine each fragment with each other fragment
      for
      (
        storage::List< storage::Vector< size_t> >::const_iterator itr_first_frag( first_mol_frags.Begin()),
          itr_first_frag_end( first_mol_frags.End());
        itr_first_frag != itr_first_frag_end;
        ++itr_first_frag
      )
      {

        if( itr_first_frag->GetSize() > MAX_FRAG_SIZE)
        {
          continue;
        }

        // get which atom from the first fragment's subgraph was in the broken bond
        size_t first_atom
        (
          itr_first_frag->Find( first_broken_bond.First()) < itr_first_frag->GetSize()
          ? first_broken_bond.First()
          : first_broken_bond.Second()
        );

        for
        (
          storage::List< storage::Vector< size_t> >::const_iterator itr_second_frag( second_mol_frags.Begin()),
            itr_second_frag_end( second_mol_frags.End());
            itr_second_frag != itr_second_frag_end;
          ++itr_second_frag
        )
        {

          if( itr_second_frag->GetSize() > MAX_FRAG_SIZE)
          {
            continue;
          }

          // get which atom from the second fragment's subgraph was in the broken bond
          size_t second_atom
          (
            itr_second_frag->Find( second_broken_bond.First()) < itr_second_frag->GetSize()
            ? second_broken_bond.First()
            : second_broken_bond.Second()
          );

          storage::Vector< size_t> first_mapping;
          first_mapping.AllocateMemory( first_mol_graph.GetSize());
          storage::Vector< size_t> second_mapping;
          second_mapping.AllocateMemory( second_mol_graph.GetSize());

          util::SiPtr< storage::Vector< size_t> > sip_first_map( &first_mapping);
          util::SiPtr< storage::Vector< size_t> > sip_second_map( &second_mapping);

          // Get the atom graphs and the mappings for each fragment
          ConformationGraphConverter::t_AtomGraph first_frag_graph
          (
            first_mol_atom_graph.GetSubgraph( *itr_first_frag, sip_first_map)
          );

          ConformationGraphConverter::t_AtomGraph second_frag_graph
          (
            second_mol_atom_graph.GetSubgraph( *itr_second_frag, sip_second_map)
          );

          // Make FragmentCompletes
          FragmentComplete first_frag( graph_maker.CreateAtomsFromGraph( first_frag_graph), "");
          FragmentComplete second_frag( graph_maker.CreateAtomsFromGraph( second_frag_graph), "");

          // Determine if the bond should be conjugated or not
          const bool is_conjugated
          (
            first_frag.GetAtomVector()( first_mapping( first_atom)).GetAtomType()->IsConjugated() &&
            second_frag.GetAtomVector()( second_mapping( second_atom)).GetAtomType()->IsConjugated()
          );

          storage::Pair< bool, FragmentComplete> merged
          (
            MergeFragmentComplete::MergeFragments
            (
              first_frag,
              second_frag,
              is_conjugated
                ? GetConfigurationalBondTypes().e_ConjugatedSingleBond
                : GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
              storage::Pair< size_t, size_t>
              (
                first_mapping( first_atom),
                second_mapping( second_atom)
              )
            )
          );

          if( merged.First())
          {
            // If we are checking for bad bonds, do so
            if( !m_AllowBadBonds)
            {
              bool bad_bonds_result( FragmentEvolveBase::IsConstitutionDruglike( merged.Second()));
              if( bad_bonds_result)
              {
                BCL_MessageVrb( GetAlias() + ": bad bonds were formed, rejecting one molecule");
                continue;
              }
            }
            util::ShPtr< FragmentComplete> mol_3D( FragmentEvolveBase::FinalizeMolecule( merged.Second()));
            if( mol_3D.IsDefined() && mol_3D->GetNumberAtoms())
            {
              new_molecules.PushBack( mol_3D);
            }
          }
        }
      }
      return new_molecules;
    }

    //! @brief add a fragment to a molecule
    //! @param MOLECULE the molecule to add a fragment to
    //! @return a new molecule with an added fragment
    util::ShPtrVector< FragmentComplete> FragmentEvolveImplementations::MutateAdd
    (
      const FragmentComplete &MOLECULE,
      const FragmentEnsemble &FRAGMENTS
    ) const
    {
      // If we don't have a list of fragments to add, return the an empty molecule
      if( FRAGMENTS.GetSize() == 0 || MOLECULE.GetNumberAtoms() == 0)
      {
        return util::ShPtrVector< FragmentComplete>();
      }

      // A modifiable form of the molecule (for removing hydrogens, if needed)
      FragmentComplete molecule_mutable( MOLECULE);

      // If there are already open valences, use them.  Otherwise remove the hydrogens to open valences up
      if( CollectorValence().Collect( molecule_mutable).GetSize() == 0)
      {
        molecule_mutable.RemoveH();
      }

      if( CollectorValence().Collect( molecule_mutable).GetSize() == 0)
      {
        return util::ShPtrVector< FragmentComplete>();
      }

      // Create the fragment grower - adds fragments to open atoms
      util::ShPtr< FragmentEnsemble> null_ptr;
      FragmentGrow fragment_grower
      (
        null_ptr,
        CollectorValence(),
        PickAtomRandom(),
        PickFragmentRandom()
      );

      // Add fragment and return the result
      util::ShPtr< FragmentComplete> result
      (
        fragment_grower.AddFragmentFromList( molecule_mutable, FRAGMENTS)
      );

      if( !result.IsDefined())
      {
        BCL_MessageVrb( "frag add: could not add a fragment");
        return util::ShPtrVector< FragmentComplete>();
      }

      result->RemoveH();

      // If we are checking for bad bonds, do so
      if( !m_AllowBadBonds)
      {
        bool bad_bonds_result( FragmentEvolveBase::IsConstitutionDruglike( *result));
        if( bad_bonds_result)
        {
          BCL_MessageVrb( GetAlias() + ": contained bad bonds, rejecting");
          return util::ShPtrVector< FragmentComplete>();
        }
      }

      util::ShPtrVector< FragmentComplete> new_mol;
      util::ShPtr< FragmentComplete> finalized_mol( FragmentEvolveBase::FinalizeMolecule( *result));
      if( finalized_mol.IsDefined() && finalized_mol->GetNumberAtoms())
      {
        new_mol.PushBack( finalized_mol);
      }
      return new_mol;
    }

    //! @brief remove a fragment from a molecule
    //! @param MOLECULE the molecule to remove a fragment from
    //! @return a new molecule with a removed fragment
    util::ShPtrVector< FragmentComplete> FragmentEvolveImplementations::MutateDel( const FragmentComplete &MOLECULE) const
    {
      FragmentComplete molecule_mutable( MOLECULE);

      molecule_mutable.RemoveH();

      if( !molecule_mutable.GetNumberAtoms())
      {
        return util::ShPtrVector< FragmentComplete>();
      }

      ConformationGraphConverter graph_maker
      (
        ConformationGraphConverter::e_ElementType,
        ConfigurationalBondTypeData::e_BondOrderOrAromaticWithRingness
      );

      graph::ConstGraph< size_t, size_t> mol_graph( graph_maker( molecule_mutable));

      storage::List< storage::Vector< size_t> > fragments
      (
        FragmentsFromRandomBondBreakage( molecule_mutable, mol_graph)
      );

      if( fragments.GetSize() != 2)
      {
        BCL_MessageVrb( "could not delete pieces of the molecule");
        return util::ShPtrVector< FragmentComplete>();
      }

      // atom graphs
      ConformationGraphConverter::t_AtomGraph mol_atom_graph
      (
        graph_maker.CreateGraphWithAtoms( molecule_mutable)
      );

      util::ShPtrVector< FragmentComplete> new_molecules;

      // Iterate through the fragments, make FragmentCompletes, and add them to new_fragments
      for
      (
        storage::List< storage::Vector< size_t> >::const_iterator itr_frag( fragments.Begin()), itr_frag_end( fragments.End());
        itr_frag != itr_frag_end;
        ++itr_frag
      )
      {
        ConformationGraphConverter::t_AtomGraph frag_graph( mol_atom_graph.GetSubgraph( *itr_frag));

        FragmentComplete new_frag( graph_maker.CreateAtomsFromGraph( frag_graph), "");

        util::ShPtr< FragmentComplete> new_molecule( FragmentEvolveBase::FinalizeMolecule( new_frag));

        if( new_molecule.IsDefined() && new_molecule->GetNumberAtoms())
        {
          new_molecules.PushBack( new_molecule);
        }
      }
      return new_molecules;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief determines what fragments would result from breaking a bond in a graph
    //! @param MOLECULE_GRAPH the graph that will have its bond broken
    //! @param FROM one vertex that makes up the bond to break
    //! @param TO the other vertex
    //! @return a list of vectors of indices which correspond to connected components of the graph
    storage::List< storage::Vector< size_t> > FragmentEvolveImplementations::CollectFragmentsFromBondBreakage
    (
      graph::ConstGraph< size_t, size_t> &MOLECULE_GRAPH,
      const size_t &FROM,
      const size_t &TO
    ) const
    {
      if( FROM >= MOLECULE_GRAPH.GetSize() || TO >= MOLECULE_GRAPH.GetSize() || FROM == TO)
      {
        return storage::List< storage::Vector< size_t> >();
      }

      // Save the bond info
      size_t bond_info( MOLECULE_GRAPH.GetEdgeData( FROM, TO));

      // Break the bond
      MOLECULE_GRAPH.RemoveEdge( FROM, TO);

      // Get the pieces of the graph
      storage::List< storage::Vector< size_t> > components( graph::Connectivity::GetComponents( MOLECULE_GRAPH));

      // Restore the bond
      MOLECULE_GRAPH.AddEdge( FROM, TO, bond_info);

      return components;
    }

    //! @brief determines what fragments would result from breaking a bond in a graph
    //! @param MOLECULE the molecule that will have a bond broken
    //! @param MOLECULE_GRAPH the graph MOLECULE
    //! @return a list of vectors of indices which correspond to connected components of the graph
    storage::List< storage::Vector< size_t> > FragmentEvolveImplementations::FragmentsFromRandomBondBreakage
    (
      const FragmentComplete &MOLECULE,
      graph::ConstGraph< size_t, size_t> &MOLECULE_GRAPH,
      const size_t &EDGE_TYPE
    ) const
    {
      // Make sure everything matches
      if( MOLECULE_GRAPH.GetSize() == 0 || MOLECULE_GRAPH.GetSize() != MOLECULE.GetNumberAtoms())
      {
        return storage::List< storage::Vector< size_t> >();
      }

      // Get a list of bonds of the molecule
      storage::Vector< sdf::BondInfo> bonds( MOLECULE.GetBondInfo());

      // Determine which ones can be broken; don't break ring bonds
      storage::Vector< size_t> available_bonds;
      available_bonds.AllocateMemory( bonds.GetSize());

      for( size_t pos( 0), end( bonds.GetSize()); pos < end; ++pos)
      {
        // Check to make sure the edge isn't in a ring, and it matches the edge type
        if
        (
          !bonds( pos).GetConstitutionalBondType()->IsBondInRing()
          && MOLECULE_GRAPH.GetEdgeData( bonds( pos).GetAtomIndexLow(), bonds( pos).GetAtomIndexHigh()) == EDGE_TYPE
        )
        {
          available_bonds.PushBack( pos);
        }
      }

      if( !available_bonds.GetSize())
      {
        return storage::List< storage::Vector< size_t> >();
      }

      size_t which_bond( random::GetGlobalRandom().Random< size_t>( available_bonds.GetSize() - 1));

      return CollectFragmentsFromBondBreakage
          (
            MOLECULE_GRAPH,
            bonds( available_bonds( which_bond)).GetAtomIndexLow(),
            bonds( available_bonds( which_bond)).GetAtomIndexHigh()
          );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief reads in fragments from a file
    //! @param FILENAME the file to read fragments from
    //! @return true if the fragments were read
    bool FragmentEvolveImplementations::ReadFragmentsFromFile() const
    {
      BCL_MessageStd( "Reading fragments from " + m_FragmentFilename);

      util::ShPtr< FragmentEnsemble> new_ensemble( new FragmentEnsemble);
      io::IFStream input;
      io::File::MustOpenIFStream( input, m_FragmentFilename);
      new_ensemble->ReadMoreFromMdl( input);
      io::File::CloseClearFStream( input);

      if( !new_ensemble->GetSize())
      {
        BCL_MessageStd( m_FragmentFilename + " did not contain any fragments");
        return false;
      }

      m_FragmentPool = new_ensemble;
      return true;

    }

    //! @brief called after TryRead successfully reads a serializer containing only initializer info (no data variables)
    //! @param SERIALIZER the serializer object with initialization information
    //! @param ERR_STREAM stream to write out errors to
    //! @return true, unless there were new errors
    bool FragmentEvolveImplementations::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &SERIALIZER,
      std::ostream &ERR_STREAM
    )
    {
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer FragmentEvolveImplementations::GetSerializer() const
    {
      io::Serializer member_data;
      switch( m_EvolveType)
      {
        case e_Clone:
          member_data.SetClassDescription( "makes a copy of a molecule");
          break;
        case e_FragAdd:
          member_data.SetClassDescription( "adds a fragment to a molecule");
          member_data.AddInitializer
          (
            "file",
            "filename containing fragments to add",
            io::Serialization::GetAgentWithCheck( &m_FragmentFilename, command::ParameterCheckFileExistence())
          );
          break;
        case e_FragDel:
          member_data.SetClassDescription( "deletes a fragment from a molecule");
          break;
        case e_Combine:
          member_data.SetClassDescription( "combines two molecules together using a single bond");
          break;
        default:
          BCL_Assert( m_EvolveType < s_End, "unknown operation type");
          break;
      }
      if( m_EvolveType != e_Clone && m_EvolveType != e_FragDel)
      {
        member_data.AddInitializer
        (
          "allow bad bonds",
          "whether bad bonds (two heteroatoms of the same type) should be allowed to form",
          io::Serialization::GetAgent( &m_AllowBadBonds),
          "true"
        );
      }
      return member_data;
    }

  } // namespace chemistry
} // namespace bcl
