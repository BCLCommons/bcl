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
        std::string rxn_line( GetReactionFromID( rxn_id));
        SmirksReactor reactor( rxn_line);
        size_t start_mol_rxn_pos( GetRandomReactionPosition( m_AllowedRxnPositionIndices));

        // identify the dummy atom(s) based on the starting mol position in the rxn
        storage::Vector< ElementType> start_mol_dummy_elements( reactor.ParseDummyAtom( reactor.m_SmirksReagents( start_mol_rxn_pos)));
        storage::Map< ElementType, ConfigurationalBondType> start_mol_dummy_bondtypes( reactor.ParseDummyAtomBonds( reactor.m_SmirksReagents( start_mol_rxn_pos)));

        // remove dummy atom and track attached index (map this index from the element type of the removed dummy atom)
        storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> stripped_mol
        (
          RemoveDummyElement( FRAGMENT, start_mol_dummy_elements)
        );

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
          storage::Vector< ElementType> dummy_elements( reactor.ParseDummyAtom( reactor.m_SmirksReagents( rxn_pos( i))));
          storage::Map< ElementType, ConfigurationalBondType> dummy_bondtypes( reactor.ParseDummyAtomBonds( reactor.m_SmirksReagents( rxn_pos( i))));
          start_mol_dummy_bondtypes.InsertElements( dummy_bondtypes.Begin(), dummy_bondtypes.End());
          FragmentComplete reagent( GetRandomReagent( rxn_id, rxn_pos( i)));
          storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> stripped_reagent
          (
            RemoveDummyElement( reagent, dummy_elements)
          );
          stripped_reagents.PushBack( stripped_reagent);
        }

        // perform our reaction(s)
        storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> reagents_product( ReactFragments( stripped_reagents, start_mol_dummy_bondtypes));
        storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> product( ReactFragments( stripped_mol, reagents_product, start_mol_dummy_bondtypes));

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
        AtomVector< AtomComplete> atoms( product.First().GetAtomVector());

        // Remove hydrogen atoms to allow bond type adjustment
        HydrogensHandler::Remove( atoms);
        if( m_ScaffoldFragment.GetSize())
        {
          return math::MutateResult< FragmentComplete>( cleaner.Clean( atoms, m_ScaffoldFragment, m_DrugLikenessType), *this);
        }
        else
        {
          return math::MutateResult< FragmentComplete>( cleaner.Clean( atoms, product.First(), m_DrugLikenessType), *this);
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
      std::vector< SmilesReactionComponent> reagent_smiles
      (
        m_AssociatedReactions.find( std::pair< std::string, size_t>( REACTION_ID, REACTION_POS))->second
      );
      const size_t rand_pos( random::GetGlobalRandom().Random<size_t>( 0, reagent_smiles.size()));

      // note that hydrogen atoms are not added if they are not present in SMILES
      return smiles::RdkitSmilesParser::ConvertSMILESToMOL( reagent_smiles[ rand_pos].m_ReagentSmiles);
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
      // the final molecule will contain all dummy atom attachment sites minus the ones used to connect here
      storage::Set< ElementType> uniq_keys( REAGENT_A.Second().GetKeys());
      for
      (
          auto key_itr( REAGENT_B.Second().GetKeys().Begin()),
          key_itr_end( REAGENT_B.Second().GetKeys().End());
          key_itr != key_itr_end;
          ++key_itr
      )
      {
        uniq_keys.Insert( *key_itr);
      }

      // we can only attach at mutually matched dummy atom sites, so always get the smaller map
      storage::Vector< ElementType> keys
      (
        REAGENT_A.Second().GetSize() < REAGENT_B.Second().GetSize() ?
        REAGENT_A.Second().GetKeysAsVector() :
        REAGENT_B.Second().GetKeysAsVector()
      );

      // our final molecule will begin from our first reagent
      storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> product( REAGENT_A);
      storage::Map< ElementType, storage::Pair< size_t, size_t>> product_react_atoms; // this is a pair now for intramolecular attachment

      // for all potential dummy elements
      for( size_t e_i( 0), e_sz( keys.GetSize()); e_i < e_sz; ++e_i)
      {
        // must mutually match or else we cannot connect
        if( !REAGENT_A.Second().Count( keys( e_i)) || !REAGENT_B.Second().Count( keys( e_i)))
        {
          continue; // TODO add a second condition for if we have 2 of the same dummy atoms in the same molecule
        }

        // convenience
        const FragmentComplete &a_frag( REAGENT_A.First());
        const FragmentComplete &b_frag( REAGENT_B.First());
        const size_t a_connect( REAGENT_A.Second().Find( keys( e_i))->second);
        const size_t b_connect( REAGENT_B.Second().Find( keys( e_i))->second);

        // join our fragments via the mapped dummy atoms;
        // OpenValence step can be skipped because we have done an equivalent step in RemoveDummyElement
        storage::Pair< bool, FragmentComplete> new_fragment
        (
          MergeFragmentComplete::MergeFragments
          (
            a_frag, // these atoms occur first in sequence
            b_frag, // these atoms are appended to the list of atoms from a_frag
            BONDS.Find( keys( e_i))->second,
//            GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
            storage::Pair< size_t, size_t>( a_connect, b_connect)
          )
        );

        // do not update
        if( !new_fragment.First()) // TODO should this be an Assert?
        {
          BCL_MessageStd( "[WARNING] FragmentMutateSmilesReact::ReactFragments merge fragments failed! Skipping reaction step!");
          continue;
        }

        // update the remaining dummy atom attachment sites, the attachment indices, and finally the molecule
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

      // perform intramolecular reactions
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
//              GetConfigurationalBondTypes().e_NonConjugatedSingleBond
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
      AtomVector< AtomComplete> product_atoms( product_atominfo, product_bondinfo);
      AtomsCompleteStandardizer standardizer( product_atoms, "", true);
      standardizer.SetConjugationOfBondTypes( product_atoms);
      FragmentComplete final_product( product_atoms, "");
      final_product.StandardizeBondLengths();
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
    {
      // copy
      size_t n_reagents( REAGENTS.GetSize());
      storage::Vector< storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>>> reagents( REAGENTS);

      // react until only one molecule remaining
      while( n_reagents > size_t( 1))
      {
        storage::Pair< FragmentComplete, storage::Map< ElementType, size_t>> product
        (
          // react the last two reagents in the vector
          ReactFragments( reagents( n_reagents - 1), reagents( n_reagents - 2), BONDS)
        );

        // kick the two reagents off the back of the vector
        reagents.PopBack();
        reagents.PopBack();

        // add the product to the end of the vector
        reagents.PushBack( product);

        // increment counter
        n_reagents = reagents.GetSize();
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
        const storage::Vector< std::string> line( util::SplitString( *line_itr, " ") );
        SmilesReactionComponent rxn_components
        (
          line( 0),                                                 // Reagent SMILES
          util::ConvertStringToNumericalValue< size_t>( line( 1)),  // Reagent ID
          line( 2),                                                 // Reaction ID
          util::ConvertStringToNumericalValue< size_t>( line( 3))   // Reaction position
        );

        // if we already have a key and we need to update the value
        const std::pair< std::string, size_t> key_pair( rxn_components.m_ReactionID, rxn_components.m_ReactionPosition);
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
      }

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
        const storage::Vector< std::string> line( util::SplitString( *line_itr, " ") );
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
