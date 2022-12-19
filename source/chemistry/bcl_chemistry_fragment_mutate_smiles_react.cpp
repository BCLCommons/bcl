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
#include "chemistry/bcl_chemistry_fragment_mutate_smiles_react.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_merge_fragment_complete.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "find/bcl_find_collector_interface.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "random/bcl_random_uniform_distribution.h"
#include "smiles/bcl_smiles_rdkit_smiles_parser.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically
#include <algorithm>
#include <fstream>

using bcl::chemistry::MergeFragmentComplete;

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentMutateSmilesReact::s_Instance
    (
      util::Enumerated< FragmentMutateInterface>::AddInstance( new FragmentMutateSmilesReact())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentMutateSmilesReact::FragmentMutateSmilesReact() :
      m_ReactionFilename( std::string())
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief construct with a pool of external fragments for fragment grow
    //! @param FRAGMENT_POOL external fragments to add to base fragment
    FragmentMutateSmilesReact::FragmentMutateSmilesReact
    (
      const bool &CORINA_CONFS
    ) :
      m_ReactionFilename( std::string())
    {
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief druglikeness constructor
    //! @param FRAGMENT_POOL external fragments to add to base fragment
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    FragmentMutateSmilesReact::FragmentMutateSmilesReact
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const bool &CORINA_CONFS
    ) :
      m_ReactionFilename( std::string())
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
    FragmentMutateSmilesReact::FragmentMutateSmilesReact
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const bool &CORINA_CONFS
    ) :
      m_ReactionFilename( std::string())
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
    FragmentMutateSmilesReact::FragmentMutateSmilesReact
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
      m_ReactionFilename( std::string())
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
    FragmentMutateSmilesReact::FragmentMutateSmilesReact
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
      m_ReactionFilename( std::string())
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
    FragmentMutateSmilesReact *FragmentMutateSmilesReact::Clone() const
    {
      return new FragmentMutateSmilesReact( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentMutateSmilesReact::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentMutateSmilesReact::GetAlias() const
    {
      static const std::string s_name( "SmilesReact");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @return MutateResult with Constitution after the mutate
    math::MutateResult< FragmentComplete> FragmentMutateSmilesReact::operator()( const FragmentComplete &FRAGMENT) const
    {
      // mutate label
      BCL_MessageStd( "SmilesReact!");
      io::OFStream debug_out;

      // try N times to have a successful mutate
      size_t try_index( 0);
      for( ; try_index < m_NumberMaxAttempts; ++try_index)
      {

        // sanity check
        if( !m_InitializedReactions || !m_InitializedReagents)
        {
          BCL_MessageStd("FragmentMutateSmilesReact::operator() not initialized! Returning null...");
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }

        // choose a reaction and reagents
        std::string rxn_id( GetRandomReactionID());
//        BCL_Debug( rxn_id);
        std::string rxn_line( GetReactionFromID( rxn_id));
//        BCL_Debug( rxn_line);
        SmirksReactor reactor( rxn_line);
        size_t start_mol_rxn_pos( GetRandomReactionPosition( m_AllowedRxnPositionIndices));

        // identify the dummy atom(s) based on the starting mol position in the rxn
        BCL_Debug( start_mol_rxn_pos);
        BCL_MessageStd( "Parse start dummy atoms.");
        storage::Vector< ElementType> start_mol_dummy_elements( reactor.ParseDummyAtom( reactor.m_SmirksReagents( start_mol_rxn_pos)));
        BCL_MessageStd( "Parse start dummy atom bonds");
        storage::Map< ElementType, ConfigurationalBondType> start_mol_dummy_bondtypes( reactor.ParseDummyAtomBonds( reactor.m_SmirksReagents( start_mol_rxn_pos)));

        // remove dummy atom and track attached index (map this index from the element type of the removed dummy atom)
        storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> stripped_mol
        (
          RemoveDummyElement( FRAGMENT, start_mol_dummy_elements)
        );
        io::File::MustOpenOFStream( debug_out, "stripped_mol.sdf");
        stripped_mol.First().WriteMDL( debug_out);
        io::File::CloseClearFStream( debug_out);

        //!!
        // TODO Need to know what the bond type is that will connect the two fragments
        //!!

        // create fragment for the other reagent
        storage::Vector< size_t> rxn_pos;
        for( size_t i( 0); i < reactor.m_NumberReagents; ++i)
        {
          if( i != start_mol_rxn_pos)
          {
            rxn_pos.PushBack( i);
          }
        }

        // remove dummy atom(s) and track attached indices (map from element type)
        storage::Vector< storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>>> stripped_reagents;
        for( size_t i( 0); i < rxn_pos.GetSize(); ++i)
        {
          BCL_Debug( rxn_pos( i));
          BCL_MessageStd( "Parse reagent dummy atoms.");
          storage::Vector< ElementType> dummy_elements( reactor.ParseDummyAtom( reactor.m_SmirksReagents( rxn_pos( i))));
          BCL_MessageStd( "Parse reagent dummy atom bonds");
          storage::Map< ElementType, ConfigurationalBondType> dummy_bondtypes( reactor.ParseDummyAtomBonds( reactor.m_SmirksReagents( rxn_pos( i))));
          start_mol_dummy_bondtypes.InsertElements( dummy_bondtypes.Begin(), dummy_bondtypes.End());
          BCL_MessageStd( "GetRandomReagent!");
          FragmentComplete reagent( GetRandomReagent( rxn_id, rxn_pos( i)));
          storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> stripped_reagent
          (
            RemoveDummyElement( reagent, dummy_elements)
          );
          io::File::MustOpenOFStream( debug_out, "stripped_reagent." + util::Format()( i) + ".sdf");
          stripped_reagent.First().WriteMDL( debug_out);
          io::File::CloseClearFStream( debug_out);
          stripped_reagents.PushBack( stripped_reagent);
        }
        stripped_reagents.PushBack( stripped_mol);

        // perform our reaction(s)
//        BCL_MessageStd(" React reagent fragments!");
        // TODO but what if the reagents do not react together, but they both react with the start mol?
//        storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> reagents_product( ReactFragments( stripped_reagents, start_mol_dummy_bondtypes));
        storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> product( ReactFragments( stripped_reagents, start_mol_dummy_bondtypes));
//        BCL_MessageStd(" React reagent product with start molecule!");
//        storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> product( ReactFragments( stripped_mol, reagents_product, start_mol_dummy_bondtypes));
//        storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> product( ReactFragments( stripped_mol, stripped_reagents( 0), start_mol_dummy_bondtypes));

        // DEBUG out
        io::File::MustOpenOFStream( debug_out, "debug_out.sdf");
        product.First().WriteMDL( debug_out);
        io::File::CloseClearFStream( debug_out);
        io::File::MustOpenOFStream( debug_out, "debug_out.smiles");
        product.First().WriteSMILES( debug_out);
        io::File::CloseClearFStream( debug_out);

        // for cleaning and optimizing the new molecule conformer
        BCL_MessageStd( "Build my cleaner");
        FragmentMapConformer cleaner
        (
          m_DrugLikenessType,
          m_MDL,
          FRAGMENT.GetMDLProperty( m_MDL),
          m_PropertyScorer,
          m_ResolveClashes,
          m_BFactors,
          m_Corina,
          storage::Vector< size_t>(), // moveable indices (empty)
          false, // choose best aligned conf (false)
          true, // fix geometry (true)
          1, // adjacent nbrs (1)
          true, // map subgraph rings (true)
          false, // skip confgen (false)
          true // geo opt MMFF94s (false)
        );

        // clean and output
        BCL_MessageStd( "Getting product atom v");
        AtomVector< AtomComplete> atoms( product.First().GetAtomVector());

        // Remove hydrogen atoms to allow bond type adjustment
        HydrogensHandler::Remove( atoms);
        if( m_ScaffoldFragment.GetSize())
        {
          return math::MutateResult< FragmentComplete>( cleaner.Clean( atoms, m_ScaffoldFragment, m_DrugLikenessType), *this);
        }
        else
        {
          return math::MutateResult< FragmentComplete>( cleaner.Clean( atoms, stripped_mol.First(), m_DrugLikenessType), *this);
        }
      }
      // return null on failure
      return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set medchem fragment library from filename
    void FragmentMutateSmilesReact::SetFragmentLibraryFromFilename( const std::string &FRAGMENTS_FILENAME)
    {
      s_Mutex.Lock();
      io::IFStream input;
      io::File::MustOpenIFStream( input, FRAGMENTS_FILENAME);

      // TODO
      // do something

      io::File::CloseClearFStream( input);
      s_Mutex.Unlock();
    }

    //! @brief get a random reagent matching the reaction id
    //! @return a reagent molecule
    // TODO: (1) make a variant of this for getting a fragment probabilistically based on difference in substructure Tanimoto similarity
    // TODO: (2) add option to make 3D conformer of reagent (add hydrogens, generate conformer, remove hydrogens, return)
    FragmentComplete FragmentMutateSmilesReact::GetRandomReagent( const std::string &REACTION_ID, const size_t REACTION_POS) const
    {
      BCL_Debug( REACTION_ID);
      BCL_Debug( REACTION_POS);
      std::vector< SmilesReactionComponent> reagent_smiles
      (
        m_AssociatedReactions.find( std::pair< std::string, size_t>( REACTION_ID, REACTION_POS))->second
      );
      BCL_Debug( reagent_smiles.size());
      const size_t rand_pos( random::GetGlobalRandom().Random<size_t>( 0, reagent_smiles.size() - 1));
      BCL_Debug( rand_pos);
//      for( size_t i( 0); i < reagent_smiles.size(); ++i)
//      {
//        BCL_Debug( reagent_smiles[ i].m_ReagentSmiles);
//      }
      BCL_Debug( reagent_smiles[ rand_pos].m_ReagentSmiles);
      return smiles::RdkitSmilesParser::ConvertSMILESToMOL
          (
            reagent_smiles[ rand_pos].m_ReagentSmiles,
            false,
            true,
            true,
            100,
            "MMFF94s",
            10.0,
            true
          );
    }


    //! @brief get a random reaction
    //! @return a reaction ID string
    std::string FragmentMutateSmilesReact::GetRandomReactionID() const
    {
      size_t pos( random::GetGlobalRandom().Random<size_t>( m_ReactionIDs.GetSize() - 1));
      return m_ReactionIDs( pos);
    }


    //! @brief get a random reaction
    //! @return a reaction ID string
    std::string FragmentMutateSmilesReact::GetReactionFromID( const std::string &RXN_ID) const
    {
      return m_Reactions.find( RXN_ID)->second;
    }


    //! @brief get random allowed reaction position from user-specified options
    size_t FragmentMutateSmilesReact::GetRandomReactionPosition( const storage::Vector< size_t> &RXN_POSITIONS) const
    {
      size_t pos( random::GetGlobalRandom().Random<size_t>( RXN_POSITIONS.GetSize() - 1));
      return RXN_POSITIONS( pos);
    }


    //! @brief generate a 3D conformer of a molecule
    void FragmentMutateSmilesReact::Generate3DConformer( FragmentComplete &MOLECULE) const
    {
      static RotamerLibraryFile rotamer_library;
      static SampleConformations sample_confs
      (
        rotamer_library,  // rotamer library file
        "SymmetryRMSD",   // conformation comparer type
        0.25,             // conformational comparer tolerance
        1,                // number of conformations
        1000,             // number of iterations
        false,            // change chirality?
        0.0,              // random dihedral change weight
        false,            // generate 3d?
        0.1,              // clash tolerance
        true              // cluster?
      );
      auto confs( sample_confs( MOLECULE).First());
      MOLECULE = confs.GetMolecules().FirstElement();
    }


    //! @brief combine two fragments through their pseudoreaction scheme
    //! @details attach two fragments using their mutually matched dummy atom element types
    //! and mapped attachment atom indices. After attaching the two fragments, it is
    //! possible that the new product molecule contains additional mutually matched
    //! dummy elements within a single fragment. Therefore, after combining the two
    //! fragments, this function will perform intramolecular pseudoreactions until
    //! there are no more mutually matched dummy atoms
    //! @return a product molecule with an associated map between un-reacted dummy atom
    //! element types and the attachment indices using the new product molecule indexing
    storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> FragmentMutateSmilesReact::ReactFragments
    (
      const storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> &REAGENT_A,
      const storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> &REAGENT_B,
      const storage::Map< ElementType, ConfigurationalBondType> &BONDS
    ) const
    {

      io::OFStream io_reagent;

      // the final molecule will contain all dummy atom attachment sites minus the ones used to connect here
      BCL_Debug( REAGENT_A);
      BCL_Debug( REAGENT_B);
      BCL_Debug( BONDS);
      storage::Set< ElementType> uniq_keys( REAGENT_A.Second().GetKeys());
      storage::Vector< ElementType> b_keys( REAGENT_B.Second().GetKeysAsVector());
      BCL_Debug( REAGENT_A.Second().GetKeys().GetSize());
      BCL_Debug( b_keys.GetSize());
      for
      (
          auto key_itr( b_keys.Begin()),
          key_itr_end( b_keys.End());
          key_itr != key_itr_end;
          ++key_itr
      )
      {
        BCL_Debug( key_itr->GetName());
        ElementType ele( GetElementTypes().ElementTypeLookup( ( *key_itr)->GetChemicalSymbol()));
        uniq_keys.Insert( ele);
      }
      BCL_Debug( uniq_keys);

      // we can only attach at mutually matched dummy atom sites, so always get the smaller map
      BCL_MessageStd( "Do some key stuff");
      storage::Vector< ElementType> keys
      (
        REAGENT_A.Second().GetSize() < REAGENT_B.Second().GetSize() ?
        REAGENT_A.Second().GetKeysAsVector() :
        REAGENT_B.Second().GetKeysAsVector()
      );
      BCL_Debug( keys);

      // our final molecule will begin from our first reagent
//      storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> product( REAGENT_A);
      storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> product;
      storage::Map< ElementType, storage::Pair< size_t, size_t>> product_react_atoms; // this is a pair now for intramolecular attachment

      // for all potential dummy elements
      for( size_t e_i( 0), e_sz( keys.GetSize()); e_i < e_sz; ++e_i)
      {
        // must mutually match or else we cannot connect
        BCL_Debug( e_i);
        if( !REAGENT_A.Second().Count( keys( e_i)) || !REAGENT_B.Second().Count( keys( e_i)))
        {
          BCL_MessageStd( "not mutually matched");
          continue; // TODO add a second condition for if we have 2 of the same dummy atoms in the same molecule
        }

        // convenience
        const FragmentComplete &a_frag( REAGENT_A.First());
        const FragmentComplete &b_frag( REAGENT_B.First());
        io::File::MustOpenOFStream( io_reagent, "reagent.a.check.sdf");
        a_frag.WriteMDL( io_reagent);
        io::File::CloseClearFStream( io_reagent);
        io::File::MustOpenOFStream( io_reagent, "reagent.b.check.sdf");
        b_frag.WriteMDL( io_reagent);
        io::File::CloseClearFStream( io_reagent);
        const size_t a_connect( REAGENT_A.Second().Find( keys( e_i))->second);
        const size_t b_connect( REAGENT_B.Second().Find( keys( e_i))->second);
        BCL_Debug( a_connect);
        BCL_Debug( b_connect);
        BCL_Debug( BONDS.Find( keys( e_i))->second);

        // join our fragments via the mapped dummy atoms;
        // OpenValence step can be skipped because we have done an equivalent step in RemoveDummyElement
        BCL_MessageStd( "Merge fragment");
        storage::Pair< bool, FragmentComplete> new_fragment
        (
          MergeFragmentComplete::MergeFragments
          (
            a_frag, // these atoms occur first in sequence
            b_frag, // these atoms are appended to the list of atoms from a_frag
            BONDS.Find( keys( e_i))->second,
            storage::Pair< size_t, size_t>( a_connect, b_connect)
          )
        );
        BCL_MessageStd( "Tried the merge");
        BCL_Debug( new_fragment.First());
        io::File::MustOpenOFStream( io_reagent, "reagent.post_merge.check.sdf");
        new_fragment.Second().WriteMDL( io_reagent);
        io::File::CloseClearFStream( io_reagent);

        // do not update
        if( !new_fragment.First()) // TODO should this be an Assert?
        {
          BCL_MessageStd( "[WARNING] FragmentMutateSmilesReact::ReactFragments merge fragments failed! Skipping reaction step!");
          continue;
        }

        // update the remaining dummy atom attachment sites, the attachment indices, and finally the molecule
        BCL_MessageStd( "Doing post-merge stuff");
        for
        (
            auto key_itr( uniq_keys.Begin()),
            key_itr_end( uniq_keys.End());
            key_itr != key_itr_end;
            ++key_itr
        )
        {
          // add mapped attachment sites, excluding the key we just used to attach the two reagents
          // A atom indices do not need adjustment
          // B atom indices must be offset by the size of reagent A fragment
          BCL_Debug( *key_itr);
          BCL_Debug( keys( e_i));
          if( *key_itr != keys( e_i))
          {
            product_react_atoms.Insert
            (
              storage::Pair< ElementType, storage::Pair< size_t, size_t>>
              (
                  *key_itr,
                  storage::Pair< size_t, size_t>
                  (
                      REAGENT_A.Second().Count( *key_itr) ? REAGENT_A.Second().Find( *key_itr)->second : util::GetUndefinedSize_t(),
                      REAGENT_B.Second().Count( *key_itr) ? REAGENT_B.Second().Find( *key_itr)->second + a_frag.GetSize() : util::GetUndefinedSize_t()
                  )
              )
            );
          }
        }
        product.First() = new_fragment.Second();
        break; // after joining the two fragments, any additional reactions are intramolecular
      }

      // return null if there is no product
      if( !product.First().GetSize())
      {
        return storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>>();
      }

      // perform intramolecular reactions
      BCL_MessageStd( "Intramolecular reaction stuff");
      storage::Map< ElementType, size_t> product_dummy_sites;
      storage::Vector< sdf::AtomInfo> product_atominfo( product.First().GetAtomInfo());
      storage::Vector< sdf::BondInfo> product_bondinfo( product.First().GetBondInfo());
      for
      (
          auto itr( product_react_atoms.Begin()), itr_end( product_react_atoms.End());
          itr != itr_end;
          ++itr
      )
      {
        // check that there are defined atom indices for both pair values;
        // if this is the case, then this is an intramolecular reaction and we
        // need to connect the atoms with a bond and discard the (not save) the
        // dummy atom attachment sites
        if( util::IsDefined( itr->second.First()) && util::IsDefined( itr->second.Second()))
        {
          // create a bond between those atoms
          product_bondinfo.PushBack
          (
            sdf::BondInfo
            (
              itr->second.First(),
              itr->second.Second(),
              BONDS.Find( itr->first)->second
            )
          );
        }

        // add unpaired dummy atom attachment sites to keys list for the product molecule
        else if( util::IsDefined( itr->second.First()))
        {
          product_dummy_sites.Insert( storage::Pair< ElementType, size_t>( itr->first, itr->second.First()));
        }

        else if( util::IsDefined( itr->second.Second()))
        {
          product_dummy_sites.Insert( storage::Pair< ElementType, size_t>( itr->first, itr->second.Second()));
        }
      }

      // TODO 3D?
      // create the final molecule
      BCL_MessageStd( "Finalizing recursive reagent product thing");
      AtomVector< AtomComplete> product_atoms( product_atominfo, product_bondinfo);
      AtomsCompleteStandardizer standardizer( product_atoms, "", true);
      standardizer.SetConjugationOfBondTypes( product_atoms);
      FragmentComplete final_product( product_atoms, "");
      final_product.StandardizeBondLengths();
      BCL_MessageStd( "Returning pairwise reagent product thing");
      return storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>>( final_product, product_dummy_sites);
    }

    //! @brief recursively call ReactFragments to combine reagents into a single product
    //! @return a product molecule with an associated map between un-reacted dummy atom
    //! element types and the attachment indices using the new product molecule indexing
    storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> FragmentMutateSmilesReact::ReactFragments
    (
      const storage::Vector< storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>>> &REAGENTS,
      const storage::Map< ElementType, ConfigurationalBondType> &BONDS
    ) const
//    {
//      // copy
//      size_t n_reagents( REAGENTS.GetSize());
////      BCL_Debug( n_reagents);
////      BCL_Debug( REAGENTS);
//      storage::Vector< storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>>> reagents( REAGENTS);
//
//      // react until only one molecule remaining
////      BCL_MessageStd("Pre while-loop");
//      size_t blah( 0);
//      while( n_reagents > size_t( 1))
//      {
////        BCL_MessageStd("In while-loop");
////        BCL_MessageStd( "n_reagents counter at: " + util::Format()( n_reagents));
//        storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> product
//        (
//          // react the last two reagents in the vector
//          ReactFragments( reagents( n_reagents - 1), reagents( n_reagents - 2), BONDS)
//        );
//
//        io::OFStream debug_out;
//        io::File::MustOpenOFStream( debug_out, "reagent_rxn_product." + util::Format()( blah) + ".sdf");
//        product.First().WriteMDL( debug_out);
//        io::File::CloseClearFStream( debug_out);
//
//        // kick the two reagents off the back of the vector
//        reagents.PopBack();
//        reagents.PopBack();
//
//        // add the product to the end of the vector
//        reagents.PushBack( product);
//        BCL_Debug( reagents);
//
//        // increment counter
//        n_reagents = reagents.GetSize();
//        ++blah;
//      }
////      BCL_MessageStd("Post (or skip) while-loop");
////      BCL_Debug( reagents( 0));
//      return reagents( 0);
//    }
    {
      int n_reagents( REAGENTS.GetSize());
      BCL_Debug( n_reagents);
      storage::Vector< storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>>> reagents( REAGENTS);
      io::OFStream debug_out;

      // need a pair of reagents for each reaction
      for( int i( 0); i < n_reagents; ++i)
      {
        BCL_MessageStd( "Start of for-loop i: " + util::Format()( i));
        for( int j( i + 1); j < n_reagents; ++j)
        {
          BCL_MessageStd( "Start of for-loop j: " + util::Format()( j));

          // invalid if there are not at least 2 reagents
          if( n_reagents < int( 2))
          {
            BCL_MessageStd( "Too few reagents");
            break;
          }

          // perform reaction
          io::File::MustOpenOFStream
          (
            debug_out,
            "reaction.r_" + util::Format()( i) + "_" + util::Format()( j) + ".i_" + util::Format()( i) + ".sdf",
            std::ios::app
          );
          reagents( i).First().WriteMDL( debug_out);
          io::File::CloseClearFStream( debug_out);

          io::File::MustOpenOFStream
          (
            debug_out,
            "reaction.r_" + util::Format()( i) + "_" + util::Format()( j) + ".j_" + util::Format()( j) + ".sdf",
            std::ios::app
          );
          reagents( j).First().WriteMDL( debug_out);
          io::File::CloseClearFStream( debug_out);

          storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> product
          (
            ReactFragments( reagents( i), reagents( j), BONDS)
          );

          // skip next part if no valid product formed
          if( !product.First().GetSize())
          {
            BCL_MessageStd( "Invalid product formed");
            continue;
          }

          // remove reagents from reagent pool
          reagents.RemoveElements( i, 1);
          j < i ?
            reagents.RemoveElements( j, 1) : // j-index unchanged
            reagents.RemoveElements( j - 1, 1); // j-index shifted down from removal of i

          // add new product to reagent pool
          reagents.PushBack( product);

          // reset our progression through the for-loop
          j = int( -1);
          BCL_MessageStd( "Reset j: " + util::Format()( j));
          n_reagents = reagents.GetSize();
          i = int( -1);
          BCL_MessageStd( "Reset i: " + util::Format()( i));
          break;
        }
      }
      return reagents( 0);
    }

    //! @brief remove the dummy element from the input molecule
    //! @return the new molecule with all dummy elements removed, as well
    //! as a map between the dummy element type and the atom that used to connect
    //! to that dummy element
    storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> FragmentMutateSmilesReact::RemoveDummyElement
    (
      const FragmentComplete &MOLECULE,
      const storage::Vector< ElementType> &ELEMENTS
    ) const
    {
      // does not work if there are no elements
      BCL_Assert( ELEMENTS.GetSize(), "[ERROR] FragmentMutateSmilesReact::RemoveDummyElement no elements to remove!");

      // we want to throw an error if multiple dummy elements of the same type are in our molecule
      // because then the reaction is ambiguous
      storage::Map< ElementType, size_t> element_count;

      // track the atoms to keep and the bonded atoms to the dummy atom
      storage::Vector< size_t> non_dummy_atoms;
      storage::Map< ElementType, size_t> bonded_atoms;

      // track the atom indices of the dummy atoms to be removed
      AtomVector< AtomComplete> atoms( MOLECULE.GetAtomVector());
      for( size_t atom_i( 0), n_atoms( MOLECULE.GetSize()); atom_i < n_atoms; ++atom_i)
      {
        // check if the current atom is of the element type marked for dummy atoms
        const ElementType &atom_e( atoms( atom_i).GetElementType());
        const size_t e_i( ELEMENTS.Find( atom_e));

        // if so, collect the atoms for removal
        if( e_i < ELEMENTS.GetSize())
        {
          // get the index of the atom bonded to the dummy atom; TODO assume just 1
          const auto &bonded_atom( atoms( atom_i).GetBonds().Begin()->GetTargetAtom());
          size_t bonded_atom_i( atoms.GetAtomIndex( bonded_atom));

          // account for index change when dummy atom removed
          atom_i < bonded_atom_i ?
            bonded_atoms.Insert( storage::Pair< ElementType, size_t>( ELEMENTS( e_i), bonded_atom_i - 1) ):
            bonded_atoms.Insert( storage::Pair< ElementType, size_t>( ELEMENTS( e_i), bonded_atom_i) );

          // increment the count
          if( element_count.Count( ELEMENTS( e_i)))
          {
            element_count.Find( ELEMENTS( e_i))->second += size_t( 1);
          }
          else
          {
            element_count.Insert( storage::Pair< ElementType, size_t>( ELEMENTS( e_i), 1));
          }
        }
        // save these so that we can reorder the atom vector to only have these
        else
        {
          non_dummy_atoms.PushBack( atom_i);
        }
      }

      // check for ambiguous dummy atoms
      for( size_t e_i( 0), e_sz( ELEMENTS.GetSize()); e_i < e_sz; ++e_i)
      {
          BCL_Assert
          (
            element_count.Find( ELEMENTS( e_i))->second < size_t( 2), // TODO second assert that each dummy element is found?
            "[ERROR] FragmentMutateSmilesReact: multiple dummy atoms of the same element type identified; "
            "ambiguous reaction cannot be performed."
          ); // TODO edge case - unambiguous if we are doing the reaction on a single fragment
      }

      // remove the dummy atoms
      atoms.Reorder( non_dummy_atoms);

      AtomsCompleteStandardizer standardizer( atoms, "", true);
      standardizer.SetConjugationOfBondTypes( atoms);

      // make new fragment
      FragmentComplete fragment( atoms, MOLECULE.GetName());
      fragment.StandardizeBondLengths();
      return storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>>( fragment, bonded_atoms);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the data members corresponding to reagents files
    //! @return true if initialization successful; false otherwise
    bool FragmentMutateSmilesReact::InitializeReagents()
    {
      // early stop
      if( m_ReagentsFileContents.empty())
      {
        m_InitializedReagents = false;
        return false;
      }

      // associate our reaction IDs with every reagent
      storage::Vector< std::string> lines( util::SplitString( m_ReagentsFileContents, "\n") );
      for
      (
          auto line_itr( lines.Begin()), line_itr_end( lines.End());
          line_itr != line_itr_end;
          ++line_itr
      )
      {
        // SMILES reagent_id reaction_id reaction_pos
        const storage::Vector< std::string> line( util::SplitString( *line_itr, " \t") );
//        BCL_Debug( line);
        SmilesReactionComponent rxn_components
        (
          line( 0),                                                      // Reagent SMILES
          util::ConvertStringToNumericalValue< size_t>( line( 1)),       // Reagent ID
          line( 2),                                                      // Reaction ID
          util::ConvertStringToNumericalValue< size_t>( line( 3)) - 1    // Reaction position 0-indexed
        );
//        BCL_Debug( util::ConvertStringToNumericalValue< size_t>( line( 1)));
//        BCL_Debug( util::ConvertStringToNumericalValue< size_t>( line( 3)) - 1);

        // if we already have a key and we need to update the value
        const std::pair< std::string, size_t> key_pair( rxn_components.m_ReactionID, rxn_components.m_ReactionPosition);
//        BCL_Debug( key_pair.first);
//        BCL_Debug( key_pair.second);
        if( m_AssociatedReactions.count( key_pair))
        {
          m_AssociatedReactions.find( key_pair)->second.push_back( rxn_components);
        }
        // adding a new key with the first element in the value vector
        else
        {
          std::vector< SmilesReactionComponent> smiles_v( 1, rxn_components);
          m_AssociatedReactions.insert( std::pair< std::pair< std::string, size_t>, std::vector< SmilesReactionComponent> >( key_pair, smiles_v ) );
        }

//        BCL_Debug( m_AssociatedReactions.find( key_pair)->second.back().m_ReagentSmiles);

      }
//      for( auto debug_itr( m_AssociatedReactions.begin()), debug_itr_end( m_AssociatedReactions.end()); debug_itr != debug_itr_end; ++debug_itr)
//      {
//        BCL_Debug( debug_itr->first.first);
//        BCL_Debug( debug_itr->first.second);
//        BCL_Debug( debug_itr->second.size());
//      }

      // unsuccessful if there are no reagents
      if( m_AssociatedReactions.empty())
      {
        BCL_MessageStd("FragmentMutateSmilesReact::InitializeReagents no reagents!");
        m_InitializedReagents = false;
        return false;
      }

      // ended successfully
      m_InitializedReagents = true;
      return true;
    }

    //! @brief Set the data members corresponding to reactions files
    //! @return true if initialization successful; false otherwise
    bool FragmentMutateSmilesReact::InitializeReactions()
    {
      // early stop
      if( m_ReactionFileContents.empty())
      {
        m_InitializedReactions = false;
        return false;
      }

      // associate our reaction IDs with every reagent
      storage::Vector< std::string> lines( util::SplitString( m_ReactionFileContents, "\n") );
      size_t line_index( 0);
      for
      (
          auto line_itr( lines.Begin()), line_itr_end( lines.End());
          line_itr != line_itr_end;
          ++line_itr, ++line_index
      )
      {
        // reactions are in the first column; save unique
        const storage::Vector< std::string> line( util::SplitString( *line_itr, " \t") );
        const std::string &rxn_id( line( 0));
        if( !m_Reactions.count( rxn_id))
        {
          m_ReactionIDs.PushBack( rxn_id);
          m_Reactions.insert( std::pair< std::string, std::string>( rxn_id, lines( line_index)));
        }
      }

      // unsuccessful if there are no reactions
      if( m_Reactions.empty())
      {
        BCL_MessageStd("FragmentMutateSmilesReact::InitializeReactions no reactions!");
        m_InitializedReactions = false;
        return false;
      }
      m_InitializedReactions = true;
      return true;
    }


    io::Serializer FragmentMutateSmilesReact::GetSerializer() const
    {
      io::Serializer parameters( FragmentMutateInterface::GetSerializer());
      parameters.SetClassDescription
      (
        "Chemically perturb a molecule based on a SMILES/SMIRKS/SMARTS reaction file."
      );

      parameters.AddInitializer
      (
        "reactions_filename",
        "file containing allowed reactions",
        io::Serialization::GetAgent( &m_ReactionFilename),
        ""
      );

      parameters.AddInitializer
      (
        "reagents_filename",
        "file containing allowed reagents; "
        "file should be organized in the following 4-column format: SMILES reagent_id reaction_id reaction_position; "
        "no headers should be provided in the input file",
        io::Serialization::GetAgent( &m_ReagentsFilename),
        ""
      );

      parameters.AddInitializer
      (
        "allowed_reaction_positions",
        "indices (0-indexed) of reaction positions that the starting molecule is allowed to take",
        io::Serialization::GetAgent( &m_AllowedRxnPositions),
        ""
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentMutateSmilesReact::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
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

      // read in reaction files
      if( m_ReactionFilename.size())
      {
        // read from input stream
        std::ifstream ifstream( m_ReactionFilename);
        std::stringstream rxn_file_stream;
        rxn_file_stream << ifstream.rdbuf();
        m_ReactionFileContents = rxn_file_stream.str();
      }

      // read in ring library filename
      if( m_ReagentsFilename.size())
      {
        // read from input stream
        std::ifstream ifstream( m_ReagentsFilename);
        std::stringstream rgt_file_stream;
        rgt_file_stream << ifstream.rdbuf();
        m_ReagentsFileContents = rgt_file_stream.str();
      }

      // verify that we can initialize reactions and reagents
      if( !InitializeReagents() || !InitializeReactions())
      {
        return false;
      }

      // read in allowed positions in the reaction mechanism for the starting molecule
      if( m_AllowedRxnPositions.size())
      {
        m_AllowedRxnPositionIndices.Reset();
        m_AllowedRxnPositionIndices = util::SplitStringToNumerical< size_t>( m_AllowedRxnPositions);
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
    std::istream &FragmentMutateSmilesReact::Read( std::istream &ISTREAM)
    {
//      io::Serialize::Read( m_FragmentPool, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT number of indentations
    //! @return ostream which was written to
    std::ostream &FragmentMutateSmilesReact::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
//      io::Serialize::Write( m_FragmentPool, OSTREAM, INDENT) << '\n';
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
