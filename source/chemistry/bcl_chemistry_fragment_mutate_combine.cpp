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
#include "chemistry/bcl_chemistry_fragment_mutate_combine.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_fragment_evolve_base.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "find/bcl_find_collector_interface.h"
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_subgraph.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
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
    const util::SiPtr< const util::ObjectInterface> FragmentMutateCombine::s_Instance
    (
      util::Enumerated< FragmentMutateInterface>::AddInstance( new FragmentMutateCombine())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentMutateCombine::FragmentMutateCombine() :
      m_FragmentPool( util::ShPtr< FragmentEnsemble>()),
      m_FragmentFilename( std::string()),
      m_MaxFragmentSize( util::GetUndefinedSize_t())
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief construct with a pool of external fragments for fragment grow
    //! @param FRAGMENT_POOL external fragments to add to base fragment
    FragmentMutateCombine::FragmentMutateCombine
    (
      const util::ShPtr< FragmentEnsemble> &FRAGMENT_POOL,
      const bool &CORINA_CONFS
    ) :
      m_FragmentPool( util::ShPtr< FragmentEnsemble>()),
      m_FragmentFilename( std::string()),
      m_MaxFragmentSize( util::GetUndefinedSize_t())
    {
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief druglikeness constructor
    //! @param FRAGMENT_POOL external fragments to add to base fragment
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    FragmentMutateCombine::FragmentMutateCombine
    (
      const util::ShPtr< FragmentEnsemble> &FRAGMENT_POOL,
      const std::string &DRUG_LIKENESS_TYPE,
      const bool &CORINA_CONFS
    ) :
      m_FragmentPool( FRAGMENT_POOL),
      m_FragmentFilename( std::string()),
      m_MaxFragmentSize( util::GetUndefinedSize_t())
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief local mutate constructor
    //! @param FRAGMENT_POOL external fragments to add to base fragment
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    FragmentMutateCombine::FragmentMutateCombine
    (
      const util::ShPtr< FragmentEnsemble> &FRAGMENT_POOL,
      const std::string &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const bool &CORINA_CONFS
    ) :
      m_FragmentPool( FRAGMENT_POOL),
      m_FragmentFilename( std::string()),
      m_MaxFragmentSize( util::GetUndefinedSize_t())
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
    FragmentMutateCombine::FragmentMutateCombine
    (
      const util::ShPtr< FragmentEnsemble> &FRAGMENT_POOL,
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
      m_FragmentPool( FRAGMENT_POOL),
      m_FragmentFilename( std::string()),
      m_MaxFragmentSize( util::GetUndefinedSize_t())
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
    FragmentMutateCombine::FragmentMutateCombine
    (
      const util::ShPtr< FragmentEnsemble> &FRAGMENT_POOL,
      const std::string &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const std::string &MDL,
      const bool &RESOLVE_CLASHES,
      const storage::Vector< float> &BFACTORS,
      const bool &CORINA_CONFS
    ) :
      m_FragmentPool( FRAGMENT_POOL),
      m_FragmentFilename( std::string()),
      m_MaxFragmentSize( util::GetUndefinedSize_t())
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
    FragmentMutateCombine *FragmentMutateCombine::Clone() const
    {
      return new FragmentMutateCombine( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentMutateCombine::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentMutateCombine::GetAlias() const
    {
      static const std::string s_name( "Combine");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @return MutateResult with Constitution after the mutate
    math::MutateResult< FragmentComplete> FragmentMutateCombine::operator()( const FragmentComplete &FRAGMENT) const
    {
      // mutate label
      BCL_MessageStd( "Combine!");

      // redo the whole thing n-max times; increment can also be made in an inner while-loop during atom index selection
      size_t try_index( 0);
      for( ; try_index < m_NumberMaxAttempts; ++try_index)
      {
        // select random medchem fragment
        iterate::Generic< const FragmentComplete> itr_gen( m_FragmentPool->Begin(), m_FragmentPool->End());
        itr_gen.GotoRandomPosition();
        FragmentComplete medchem_frag( *itr_gen);

        if( !FRAGMENT.GetNumberAtoms() || !medchem_frag.GetNumberAtoms())
        {
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }

        // Remove hydrogens
        // TODO modernize valence handling so that we can track atoms
        FragmentComplete first_molecule( FRAGMENT);
        first_molecule.RemoveH();
        FragmentComplete second_molecule( medchem_frag);
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
          FragmentEvolveBase::FragmentsFromRandomBondBreakage
          (
            first_molecule, first_mol_graph
          )
        );

        storage::List< storage::Vector< size_t> > second_mol_frags
        (
          FragmentEvolveBase::FragmentsFromRandomBondBreakage
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
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }

        // If everything succeeded, make atom graphs of the molecules
        ConformationGraphConverter::t_AtomGraph first_mol_atom_graph( graph_maker.CreateGraphWithAtoms( first_molecule));
        ConformationGraphConverter::t_AtomGraph second_mol_atom_graph( graph_maker.CreateGraphWithAtoms( second_molecule));

        // OwnPtrs to graphs; to be used in subgraph construction
        util::OwnPtr< graph::ConstGraph< size_t, size_t> > optr_first_graph( &first_mol_graph, false);
        util::OwnPtr< graph::ConstGraph< size_t, size_t> > optr_second_graph( &second_mol_graph, false);

        // The FragmentEnsemble that will be returned
        storage::Vector< FragmentComplete> new_molecules;

        ///////////////
        // Determine which bond was broken in each molecule
        ///////////////

        graph::Subgraph< size_t, size_t> first_subgraph( optr_first_graph, *first_mol_frags.Begin());
        storage::List< storage::Pair< size_t, size_t> > first_adj_edges( first_subgraph.GetAdjacentEdgeIndices());
        if( first_adj_edges.GetSize() != 1)
        {
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }
        storage::Pair< size_t, size_t> &first_broken_bond( *first_adj_edges.Begin());

        graph::Subgraph< size_t, size_t> second_subgraph( optr_second_graph, *second_mol_frags.Begin());
        storage::List< storage::Pair< size_t, size_t> > second_adj_edges( second_subgraph.GetAdjacentEdgeIndices());
        if( second_adj_edges.GetSize() != 1)
        {
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
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

          if( itr_first_frag->GetSize() > m_MaxFragmentSize)
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

            if( itr_second_frag->GetSize() > m_MaxFragmentSize)
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
              new_molecules.PushBack( merged.Second());
            }
          }
        }

        // exit if nothing
        if( !new_molecules.GetSize())
        {
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
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

        // clean and output
        // TODO: make some type of selector to choose which fragment to take (e.g., Largest, Smallest, Closest to some feature space, etc)
        iterate::Generic< const FragmentComplete> mol_itr( new_molecules.Begin(), new_molecules.End());
        mol_itr.GotoRandomPosition();
        AtomVector< AtomComplete> atoms( mol_itr->GetAtomVector());

        // Remove hydrogen atoms to allow bond type adjustment
        HydrogensHandler::Remove( atoms);
        if( m_ScaffoldFragment.GetSize())
        {
          return math::MutateResult< FragmentComplete>( cleaner.Clean( atoms, m_ScaffoldFragment, m_DrugLikenessType), *this);
        }
        else
        {
          return math::MutateResult< FragmentComplete>( cleaner.Clean( atoms, FRAGMENT, m_DrugLikenessType), *this);
        }
      }
      return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set medchem fragment library from filename
    void FragmentMutateCombine::SetFragmentLibraryFromFilename( const std::string &FRAGMENTS_FILENAME)
    {
      s_Mutex.Lock();
      io::IFStream input;
      io::File::MustOpenIFStream( input, FRAGMENTS_FILENAME);
      FragmentEnsemble medchem_groups;
      medchem_groups.ReadMoreFromMdl( input);
      m_FragmentPool = util::CloneToShPtr( medchem_groups);
      io::File::CloseClearFStream( input);
      s_Mutex.Unlock();
    }

  //////////////////////
  // helper functions //
  //////////////////////

    io::Serializer FragmentMutateCombine::GetSerializer() const
    {
      io::Serializer parameters( FragmentMutateInterface::GetSerializer());
      parameters.SetClassDescription
      (
        "Combines a fragment from a target molecule with a fragment from an external library. "
        "WARNING - This mutate is created from a legacy chemical perturbation function. "
        "This mutate does NOT obey atom selection rules, and it was modified from "
        "its original version to fit the FragmentMutateInterface framework for "
        "compatibility purposes. "
      );

      parameters.AddInitializer
      (
        "fragment_library",
        "path to the fragment library",
        io::Serialization::GetAgent( &m_FragmentFilename),
        ""
      );

      parameters.AddInitializer
      (
        "max_fragment_size",
        "maximum allowed size of fragments used for recombination.",
        io::Serialization::GetAgent( &m_MaxFragmentSize),
        "50"
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentMutateCombine::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
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

      // read in ring library filename
      if( m_FragmentFilename.size())
      {
        s_Mutex.Lock();
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_FragmentFilename);
        FragmentEnsemble medchem_groups;
        medchem_groups.ReadMoreFromMdl( input);
        m_FragmentPool = util::CloneToShPtr( medchem_groups);
        io::File::CloseClearFStream( input);
        s_Mutex.Unlock();
      }

      // done
      return true;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentMutateCombine::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_FragmentPool, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT number of indentations
    //! @return ostream which was written to
    std::ostream &FragmentMutateCombine::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_FragmentPool, OSTREAM, INDENT) << '\n';
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
