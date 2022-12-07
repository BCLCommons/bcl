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
#include "chemistry/bcl_chemistry_fragment_mutate_connect.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_clash_score.h"
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_bond_angle_assignment.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_flex_field.h"
#include "chemistry/bcl_chemistry_fragment_make_conformers.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_split_largest_component.h"
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
#include "chemistry/bcl_chemistry_mutate_bond_lengths.h"
#include "chemistry/bcl_chemistry_mutate_chirality.h"
#include "chemistry/bcl_chemistry_mutate_clash_resolver.h"
#include "chemistry/bcl_chemistry_mutate_dihedrals_interface.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "chemistry/bcl_chemistry_sample_conformations.h"
#include "chemistry/bcl_chemistry_search_fragment_library_from_tree.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "command/bcl_command_command_state.h"
#include "find/bcl_find_collector_interface.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "quality/bcl_quality_rmsd.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentMutateConnect::s_Instance
    (
      util::Enumerated< FragmentMutateInterface>::AddInstance( new FragmentMutateConnect())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentMutateConnect::FragmentMutateConnect() :
//        m_RotamerLibrarySearcher( util::ShPtr< SearchFragmentLibraryFromTree>()),
        m_RotamerLibrarySearcher( new SearchFragmentLibraryFromTree()),
        m_Linkers( storage::Vector< std::string>()),
        m_LinkersString( std::string()),
        m_TerminalFragment( FragmentComplete()),
        m_TerminalFragmentFilename( std::string()),
        m_TerminalFragmentMutableAtomIndices( storage::Vector< size_t>()),
        m_TerminalFragmentMutableAtoms( std::string())
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief druglikeness constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    FragmentMutateConnect::FragmentMutateConnect
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const bool &CORINA_CONFS
    ) :
//    m_RotamerLibrarySearcher( util::ShPtr< SearchFragmentLibraryFromTree>()),
      m_RotamerLibrarySearcher( new SearchFragmentLibraryFromTree()),
      m_Linkers( storage::Vector< std::string>()),
      m_LinkersString( std::string()),
      m_TerminalFragment( FragmentComplete()),
      m_TerminalFragmentFilename( std::string()),
      m_TerminalFragmentMutableAtomIndices( storage::Vector< size_t>()),
      m_TerminalFragmentMutableAtoms( std::string())
    {
      m_DrugLikenessType = DRUG_LIKENESS_TYPE;
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief full constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    FragmentMutateConnect::FragmentMutateConnect
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const bool &CORINA_CONFS
    ) :
//    m_RotamerLibrarySearcher( util::ShPtr< SearchFragmentLibraryFromTree>()),
      m_RotamerLibrarySearcher( new SearchFragmentLibraryFromTree()),
      m_Linkers( storage::Vector< std::string>()),
      m_LinkersString( std::string()),
      m_TerminalFragment( FragmentComplete()),
      m_TerminalFragmentFilename( std::string()),
      m_TerminalFragmentMutableAtomIndices( storage::Vector< size_t>()),
      m_TerminalFragmentMutableAtoms( std::string())
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
    FragmentMutateConnect::FragmentMutateConnect
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
//    m_RotamerLibrarySearcher( util::ShPtr< SearchFragmentLibraryFromTree>()),
      m_RotamerLibrarySearcher( new SearchFragmentLibraryFromTree()),
      m_Linkers( storage::Vector< std::string>()),
      m_LinkersString( std::string()),
      m_TerminalFragment( FragmentComplete()),
      m_TerminalFragmentFilename( std::string()),
      m_TerminalFragmentMutableAtomIndices( storage::Vector< size_t>()),
      m_TerminalFragmentMutableAtoms( std::string())
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
    FragmentMutateConnect::FragmentMutateConnect
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
//    m_RotamerLibrarySearcher( util::ShPtr< SearchFragmentLibraryFromTree>()),
      m_RotamerLibrarySearcher( new SearchFragmentLibraryFromTree()),
      m_Linkers( storage::Vector< std::string>()),
      m_LinkersString( std::string()),
      m_TerminalFragment( FragmentComplete()),
      m_TerminalFragmentFilename( std::string()),
      m_TerminalFragmentMutableAtomIndices( storage::Vector< size_t>()),
      m_TerminalFragmentMutableAtoms( std::string())
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
    FragmentMutateConnect *FragmentMutateConnect::Clone() const
    {
      return new FragmentMutateConnect( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentMutateConnect::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentMutateConnect::GetAlias() const
    {
      static const std::string s_name( "Connect");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @return MutateResult with Constitution after the mutate
    math::MutateResult< FragmentComplete> FragmentMutateConnect::operator()( const FragmentComplete &FRAGMENT) const
    {
      BCL_MessageStd( "Connect!");

      // pick the first atom
      util::SiPtr< const AtomConformationalInterface> first_atom, second_atom;
      size_t first_atom_index( util::GetUndefinedSize_t()), second_atom_index( util::GetUndefinedSize_t());

      // try at least once to get an atom, but bail if we fail too much
      // TODO: consider refactoring to make this atom selection code in each of these derived classes into
      // a virtual or override function in the base class. most of the time we do the same thing, but then in special
      // cases we can re-define it in the derived class (like here for second atom index).
      size_t failures( 0);
      do
      {
        // bail
        if( failures > m_NumberMaxAttempts)
        {
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }

        // pick random atom to transform
        if( m_MutableAtomIndices.GetSize() || m_MutableElements.GetSize() || m_MutableFragments.GetSize())
        {
          first_atom = this->PickAtom( FRAGMENT, false);
        }
        else
        {
          first_atom = this->PickAtom( FRAGMENT, true);
        }

        // increment count
        ++failures;

      } while( first_atom->GetElementType() == GetElementTypes().e_Hydrogen);

      // first atom index
      first_atom_index = size_t( FRAGMENT.GetAtomVector().GetAtomIndex( *first_atom));

      // do it again but this time for the terminal fragment
      failures = size_t( 0);
      do
      {
        // bail
        if( failures > m_NumberMaxAttempts)
        {
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }

        if( m_TerminalFragmentMutableAtomIndices.GetSize() || m_TerminalFragmentMutableElements.GetSize() || m_TerminalFragmentMutableFragments.GetSize())
        {
          second_atom = util::SiPtr< const AtomConformationalInterface>
          (
            FragmentTrackMutableAtoms::GetAtomFromMutable
            (
              m_TerminalFragment,
              false,
              m_TerminalFragmentMutableAtomIndices,
              m_TerminalFragmentMutableElements,
              m_TerminalFragmentMutableFragments,
              m_TerminalFragmentFixedAtomindices,
              m_TerminalFragmentFixedElements,
              m_TerminalFragmentFixedFragments,
              m_TerminalFragmentComplementMutableFragments,
              m_TerminalFragmentComplementFixedFragments,
              ConformationGraphConverter
              (
                m_TerminalFragmentMutableAtomComparisonType,
                m_TerminalFragmentMutableBondComparisonType
              ),
              ConformationGraphConverter
              (
                m_TerminalFragmentFixedAtomComparisonType,
                m_TerminalFragmentFixedBondComparisonType
              )
            )
          );
        }
        else
        {
          second_atom = util::SiPtr< const AtomConformationalInterface>
          (
            (
                FragmentTrackMutableAtoms::GetAtomFromMutable
                (
                  m_TerminalFragment,
                  true
                )
            )
          );
        }

        // increment count
        ++failures;

        // set second atom index
        second_atom_index = size_t( m_TerminalFragment.GetAtomVector().GetAtomIndex( *second_atom));

      } while( second_atom->GetElementType() == GetElementTypes().e_Hydrogen || second_atom_index == first_atom_index);

      // catch my fuckups from this stupid couple of do-while loops
      if( !util::IsDefined( first_atom_index) || !util::IsDefined( second_atom_index))
      {
        BCL_MessageStd( "FragmentMutateConnect failed to define indices for the chosen atoms! Returning NULL...");
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
      }

      // initialize sample conformations object
      // TODO make static member data and/or consider moving to base class
      // ideally would like to only have 1 SampleConfs object and 1 rotamer library
      // for entirety of FragmentMutateInterface shit, maybe give them to MapConformer at construction,
      // reduce memory footprint
      static RotamerLibraryFile rotamer_library;
      static SampleConformations sample_confs
      (
        rotamer_library,  // rotamer library file
        "",               // conformation comparer type
        0.0,              // conformational comparer tolerance
        250,              // number of conformations
        2000,             // number of iterations
        false,            // no change chirality
        0.0,              // random dihedral change weight
        false,            // generate 3d?
        0.1,              // clash tolerance
        true              // cluster?
      );
      InitializeSampleConfsGlobal( sample_confs);

      // for cleaning and optimizing the new molecule conformer
      FragmentMapConformer cleaner
      (
        m_DrugLikenessType,
        m_MDL,
        FRAGMENT.GetMDLProperty( m_MDL),
        m_PropertyScorer,
        m_ResolveClashes,
        m_BFactors,
        m_Corina,
        storage::Vector< size_t>(),
        false,
        false,
        size_t( 1)
      );

      io::OFStream output_ens;
      if( m_OutputEnsembleFilename.size())
      {
        io::File::MustOpenOFStream( output_ens, m_OutputEnsembleFilename);
      }

      // now begin the fun part and grow the linker
      FragmentComplete frag_a( FRAGMENT), frag_b( m_TerminalFragment);

      // we need to obtain the rings and indices that we will use so that we can control ring attach points
      storage::Vector< storage::Triplet< FragmentComplete, size_t, size_t> > rings_and_indices( ChooseRings( m_Linkers));

      auto fwd_linked_fragment( ExtendHalfLinkerForward( frag_a, m_Linkers, rings_and_indices, first_atom_index, cleaner, sample_confs));
      if( !fwd_linked_fragment.First().GetSize())
      {
        BCL_MessageStd
        (
          "FragmentMutateConnect could not generate conformers of the input fragment coupled with "
          "the half-linker. Returning NULL..."
        );
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
      }

      // grow the other half of the linker from fragment_b
      auto rev_linked_fragment( ExtendHalfLinkerReverse( frag_b, m_Linkers, rings_and_indices, second_atom_index, cleaner, sample_confs));
      if( !rev_linked_fragment.First().GetSize())
      {
        BCL_MessageStd
        (
          "FragmentMutateConnect could not generate conformers of the terminal fragment coupled with "
          "the half-linker. Returning NULL..."
        );
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
      }

      // connect the linker at the split point
      storage::Pair< FragmentEnsemble, storage::Vector< size_t> > confs
      (
        JoinHalfExtendedFragments
        (
          fwd_linked_fragment.First(),
          fwd_linked_fragment.Second(),
          fwd_linked_fragment.Third(),
          rev_linked_fragment.First(),
          rev_linked_fragment.Second(),
          rev_linked_fragment.Third(),
          rotamer_library
        )
      );

      if( !confs.First().GetSize())
      {
        BCL_MessageStd
        (
          "FragmentMutateConnect could not generate any valid conformers of the fully linked "
          "molecule. Consider changing your distance cutoff criteria and/or linker length. "
          "Returning NULL..."
        );
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
      }

      // make local conformers for each solution
      FragmentEnsemble final_confs;
      FragmentComplete supermol( CreateSuperMolecule( FRAGMENT, m_TerminalFragment));
      for
      (
          auto conf_itr( confs.First().Begin()), conf_itr_end( confs.First().End());
          conf_itr != conf_itr_end;
          ++conf_itr
      )
      {
        // set sampling to linker only
        InitializeSampleConfsLocal( sample_confs);
        conf_itr->GetStoredPropertiesNonConst().SetMDLProperty( "SampleByParts", confs.Second());
        storage::Pair< FragmentEnsemble, FragmentEnsemble> conf_confs( sample_confs( *conf_itr));
        final_confs.Append( conf_confs.First());
      }

      // get common indices; start by looping over all atom indices
      storage::Vector< size_t> common_atom_indices;
      for( size_t i( 0); i < final_confs.GetMolecules().FirstElement().GetSize(); ++i)
      {
        // if the atom index is NOT in the linker, keep it
        if( confs.Second().Find( i) >= confs.Second().GetSize())
        {
          common_atom_indices.PushBack( i);
        }
      }

      // get the first/best
      storage::Pair< double, size_t> best_rmsd( 999.0, util::GetUndefinedSize_t());
      size_t conf_index( 0);
      for
      (
          auto conf_itr( final_confs.Begin()), conf_itr_end( final_confs.End());
          conf_itr != conf_itr_end;
          ++conf_itr, ++conf_index
      )
      {
        const double rmsd( ComputeRMSDToSuperMolecule( *conf_itr, common_atom_indices, supermol));
        if( rmsd < best_rmsd.First())
        {
          best_rmsd.First() = rmsd;
          best_rmsd.Second() = conf_index;
        }
        if( rmsd < m_StartingFragmentRMSDCutoff && m_OutputEnsembleFilename.size())
        {
          s_Mutex.Lock();
          conf_itr->WriteMDL( output_ens);
          s_Mutex.Unlock();
        }
      }
      io::File::CloseClearFStream( output_ens);

      // get best
      BCL_MessageStd( "Saving best linked molecule with RMSD of " + util::Format()( best_rmsd.First()) + " Angstroms to starting fragments");
      storage::Vector< FragmentComplete> final_confs_v( final_confs.Begin(), final_confs.End());
      FragmentComplete new_mol( final_confs_v( best_rmsd.Second()));

      // standardize and return
      AtomVector< AtomComplete> atoms( new_mol.GetAtomVector());
      HydrogensHandler::Remove( atoms);
      atoms = AtomVector< AtomComplete>( cleaner.CleanAtoms( atoms, m_DrugLikenessType));
      if( !atoms.GetSize())
      {
        return math::MutateResult< FragmentComplete>();
      }
      return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>( util::CloneToShPtr( FragmentComplete( atoms, ""))), *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set options for the global conformational ensembles made of each fragment + half linker
    void FragmentMutateConnect::InitializeSampleConfsGlobal( SampleConformations &SAMPLER) const
    {
      SAMPLER.SetSamplingPreferences( true, true, true, false);
      SAMPLER.SetConfComparer( m_ConfComparerGlobal);
      SAMPLER.SetTolerance( m_ConfToleranceGlobal);
      SAMPLER.SetMaxNumConfs( m_NMaxConfsGlobal);
      SAMPLER.SetMaxIterations( m_NMaxItersGlobal);
      SAMPLER.SetIfCluster( m_ClusterGlobal);
      SAMPLER.SetClashResolution( m_ClashResolutionGlobal);
      SAMPLER.SetMaxClashScore( m_ClashScoreGlobal);
    }

    //! @brief set options for the local conformational ensembles made of the final joined molecule
    void FragmentMutateConnect::InitializeSampleConfsLocal( SampleConformations &SAMPLER) const
    {
      SAMPLER.SetSamplingPreferences( false, false, true, false);
      SAMPLER.SetConfComparer( m_ConfComparerLocal);
      SAMPLER.SetTolerance( m_ConfToleranceLocal);
      SAMPLER.SetMaxNumConfs( m_NMaxConfsLocal);
      SAMPLER.SetMaxIterations( m_NMaxItersLocal);
      SAMPLER.SetIfCluster( m_ClusterLocal);
      SAMPLER.SetClashResolution( m_ClashResolutionLocal);
      SAMPLER.SetMaxClashScore( m_ClashScoreLocal);
    }

    //! @brief return an ensemble of conformers for the starting fragment extended with half the linker
    //! @param FRAGMENT_A molecule from which extension will start
    //! @param LINKER_COMPOSITION vector of strings indicating type of linker to add
    //! @param LINK_INDEX_A atom in FRAGMENT_A at which to begin extension
    //! @return an ensemble of conformers for the first half of the extended linker (for odd numbers, half+1),
    //! the index of the link point to the next half of the molecule (in hydrogenated numbering), and
    //! the samplebyparts indices
    storage::Triplet< FragmentEnsemble, size_t, storage::Vector< size_t> > FragmentMutateConnect::ExtendHalfLinkerForward
    (
      const FragmentComplete FRAGMENT_A,
      const storage::Vector< std::string> &LINKER_COMPOSITION,
      const storage::Vector< storage::Triplet< FragmentComplete, size_t, size_t> > &RINGS,
      const size_t LINK_INDEX_A,
      const FragmentMapConformer &CLEANER,
      SampleConformations &CONFORMATOR
    ) const
    {
      BCL_MessageStd( "========== BEGIN ExtendHalfLinkerForward ==========");

      // get molecule size
      FragmentComplete fragment_a( FRAGMENT_A);
      size_t final_connection_index( util::GetUndefinedSize_t()), n_preceeding_h_atoms( util::GetUndefinedSize_t());

      // get number of linkers to add
      size_t n_linker_fragments
      (
        LINKER_COMPOSITION.GetSize() % size_t( 2) ?
            LINKER_COMPOSITION.GetSize() / size_t( 2) + 1 :
            LINKER_COMPOSITION.GetSize() / size_t( 2)
      );

      // add half of the linker to our starting fragment
      for( size_t i( 0); i < n_linker_fragments; ++i)
      {
        BCL_MessageStd("Linker index: " + util::Format()( i));
        BCL_MessageStd( "Linked atom index: " + util::Format()( !i ? LINK_INDEX_A : fragment_a.GetSize() - 1));
        BCL_MessageStd(
          "Current atom: " + util::Format()(
              fragment_a.GetAtomVector()( !i ? LINK_INDEX_A : fragment_a.GetSize() - 1).GetAtomType().GetName()
          )
        );
        if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Verbose))
        {
          for
          (
              auto bond_itr( fragment_a.GetAtomVector()( !i ? LINK_INDEX_A : fragment_a.GetSize() - 1).GetBonds().Begin()),
              bond_itr_end( fragment_a.GetAtomVector()( !i ? LINK_INDEX_A : fragment_a.GetSize() - 1).GetBonds().End());
              bond_itr != bond_itr_end;
              ++bond_itr
          )
          {
            BCL_MessageStd( "Bonded neighbor: " + util::Format()( bond_itr->GetTargetAtom().GetAtomType()));
          }
        }

        // the connection index is not necessarily the final index in RingLink
        size_t connection_index( !i ? LINK_INDEX_A : fragment_a.GetSize() - 1);

        // TODO: fix this
        if( RINGS( i).First().GetSize() && i == n_linker_fragments - 1)
        {
          final_connection_index = RINGS( i).Third() + fragment_a.GetSize();
        }

        // update fragment
        fragment_a = FragmentComplete
        (
          AddLinkFragment
          (
            fragment_a,                // base fragment
            connection_index,
            LINKER_COMPOSITION( i),    // current linker string
            RINGS( i).First(),
            RINGS( i).Second()
          )
        );
      }

      if( util::IsDefined( final_connection_index))
      {
        n_preceeding_h_atoms = fragment_a.GetNumberPreceedingHydrogenAtoms( final_connection_index);
      }

      // now we want to make a clean, hydrogenated version of the extended molecule to generate a high quality ensemble
      fragment_a.RemoveH();
      CLEANER.ResetMoveableIndices();
      math::MutateResult< FragmentComplete> fragment_a_clean( CLEANER.Clean( fragment_a.GetAtomVector(), FRAGMENT_A, m_DrugLikenessType), *this);

      // if we cannot clean it, return null
      if( !fragment_a_clean.GetArgument().IsDefined())
      {
        BCL_MessageStd
        (
          "Unable to create valid 3D conformer of the input fragment coupled to the half-linker. "
          "Returning NULL..."
        );
        return storage::Triplet< FragmentEnsemble, size_t, storage::Vector< size_t> >();
      }

      // obtain the mobile indices
      FragmentComplete fragment_a_h( *( fragment_a_clean.GetArgument()));
      storage::Vector< size_t> sample_by_parts_indices( util::SplitStringToNumerical< size_t>( fragment_a_h.GetMDLProperty( "SampleByParts")));

      // generate conformers
      FragmentEnsemble confs( CONFORMATOR( fragment_a_h).First());

      BCL_MessageStd( "========== END ExtendHalfLinkerForward ==========");

      // end
      return storage::Triplet< FragmentEnsemble, size_t, storage::Vector< size_t> >
          (
            confs,
            util::IsDefined( final_connection_index) ?
              final_connection_index - n_preceeding_h_atoms :
              sample_by_parts_indices( sample_by_parts_indices.GetSize() - 1),
            sample_by_parts_indices
          );
    }

    //! @brief return an ensemble of conformers for the terminal fragment extended with half the linker
    //! @param FRAGMENT_B terminal fragment from which extension will start
    //! @param LINKER_COMPOSITION vector of strings indicating type of linker to add
    //! @param LINK_INDEX_B atom in FRAGMENT_B at which to begin extension
    //! @return an ensemble of conformers for the terminal half of the extended linker (for odd numbers, half-1)
    //! the index of the link point to the other half of the molecule (in hydrogenated numbering),
    //! and the samplebyparts indices
    storage::Triplet< FragmentEnsemble, size_t, storage::Vector< size_t> > FragmentMutateConnect::ExtendHalfLinkerReverse
    (
      const FragmentComplete &FRAGMENT_B,
      const storage::Vector< std::string> &LINKER_COMPOSITION,
      const storage::Vector< storage::Triplet< FragmentComplete, size_t, size_t> > &RINGS,
      const size_t LINK_INDEX_B,
      const FragmentMapConformer &CLEANER,
      SampleConformations &CONFORMATOR
    ) const
    {
      BCL_MessageStd( "========== BEGIN ExtendHalfLinkerReverse ==========");

      // get molecule size
      FragmentComplete fragment_b( FRAGMENT_B);
      size_t final_connection_index( util::GetUndefinedSize_t()), n_preceeding_h_atoms( util::GetUndefinedSize_t());
      size_t linker_sz( LINKER_COMPOSITION.GetSize());

      // get number of linkers to add
      size_t n_linker_fragments( linker_sz / size_t( 2));

      // add half of the linker to our starting fragment
      for( size_t i( linker_sz - 1); i >= size_t( linker_sz - n_linker_fragments); --i)
      {
        BCL_MessageStd("Linker index: " + util::Format()( i));
        BCL_MessageStd( "Linked atom index: " + util::Format()( i == linker_sz - 1 ? LINK_INDEX_B : fragment_b.GetSize() - 1));
        BCL_MessageStd(
          "Current atom: " + util::Format()(
            fragment_b.GetAtomVector()
            ( !i ? LINK_INDEX_B : fragment_b.GetSize() - 1).GetAtomType().GetName()
          )
        );
        if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Verbose))
        {
          for
          (
              auto bond_itr( fragment_b.GetAtomVector()( i == linker_sz - 1 ? LINK_INDEX_B : fragment_b.GetSize() - 1).GetBonds().Begin()),
              bond_itr_end( fragment_b.GetAtomVector()( i == linker_sz - 1 ? LINK_INDEX_B : fragment_b.GetSize() - 1).GetBonds().End());
              bond_itr != bond_itr_end;
              ++bond_itr
          )
          {
            BCL_MessageStd( "Bonded neighbor: " + util::Format()( bond_itr->GetTargetAtom().GetAtomType()));
          }
        }

        // the connection index is not necessarily the final index in RingLink
        size_t connection_index( i == linker_sz - 1 ? LINK_INDEX_B : fragment_b.GetSize() - 1);

        // TODO: fix this
        if( RINGS( i).First().GetSize() && i == n_linker_fragments - 1)
        {
          final_connection_index = RINGS( i).Third() + fragment_b.GetSize();
        }

        fragment_b = FragmentComplete
        (
          AddLinkFragment
          (
            fragment_b,
            connection_index,
            LINKER_COMPOSITION( i),
            RINGS( i).First(),
            RINGS( i).Second()
          )
        );
      }

      // use this to obtain the correct atom index after adjusting hydrogen atoms
      if( util::IsDefined( final_connection_index))
      {
        n_preceeding_h_atoms = fragment_b.GetNumberPreceedingHydrogenAtoms( final_connection_index);
      }

      // now we want to make a clean, hydrogenated version of the extended molecule to generate a high quality ensemble
      fragment_b.RemoveH();
      CLEANER.ResetMoveableIndices();
      math::MutateResult< FragmentComplete> fragment_b_clean( CLEANER.Clean( fragment_b.GetAtomVector(), FRAGMENT_B, m_DrugLikenessType), *this);

      // if we cannot clean it, return null
      if( !fragment_b_clean.GetArgument().IsDefined())
      {
        BCL_MessageStd
        (
          "Unable to create valid 3D conformer of the terminal fragment coupled to the half-linker. "
          "Returning NULL..."
        );
        return storage::Triplet< FragmentEnsemble, size_t, storage::Vector< size_t> >();
      }

      // obtain the mobile indices
      FragmentComplete fragment_b_h( *( fragment_b_clean.GetArgument()));
      storage::Vector< size_t> sample_by_parts_indices( util::SplitStringToNumerical< size_t>( fragment_b_h.GetMDLProperty( "SampleByParts")));

      // generate conformers
      FragmentEnsemble confs( CONFORMATOR( fragment_b_h).First());

      BCL_MessageStd( "========== END ExtendHalfLinkerReverse ==========");

      // end
      return storage::Triplet< FragmentEnsemble, size_t, storage::Vector< size_t> >
          (
            confs,
            util::IsDefined( final_connection_index) ?
              final_connection_index - n_preceeding_h_atoms :
              sample_by_parts_indices( sample_by_parts_indices.GetSize() - 1),
            sample_by_parts_indices
          );
    }

    //! @brief propose solutions for the fully linked molecule given two extended halves
    //! @param FWD_ENS fragment that was extended with half of the linker in the forward direction
    //! @param FWD_LINK_ATOM the atom at which this fragment will be joined
    //! @param FWD_SBP_INDICES the atoms that compose conformationally flexible dihedrals
    //! @param REV_ENS fragment that was extended with half of the linker in the reverse direction
    //! @param REV_LINK_ATOM the atom at which this fragment will be joined
    //! @param REV_SBP_INDICES the atoms that compose conformationally flexible dihedrals
    //! @return a conformational ensemble of the fully linked molecule and the final linker indices
    storage::Pair< FragmentEnsemble, storage::Vector< size_t> > FragmentMutateConnect::JoinHalfExtendedFragments
    (
      const FragmentEnsemble &FWD_ENS,
      const size_t FWD_LINK_ATOM,
      const storage::Vector< size_t> &FWD_SBP_INDICES,
      const FragmentEnsemble &REV_ENS,
      const size_t REV_LINK_ATOM,
      const storage::Vector< size_t> &REV_SBP_INDICES,
      const RotamerLibraryFile &ROTLIB
    ) const
    {
      BCL_MessageStd( "========== BEGIN JoinHalfExtendedFragments ==========");

      BCL_MessageStd("FWD_LINK_ATOM: " + util::Format()( FWD_LINK_ATOM));
      BCL_MessageStd("REV_LINK_ATOM: " + util::Format()( REV_LINK_ATOM));

      // final conformer ensemble
      FragmentEnsemble final_conformer_ensemble;
      storage::Vector< size_t> linker_atom_indices, linker_atom_indices_plus_connect;
      storage::Vector< size_t> final_linker_atom_indices;

      // try all pairs
      for
      (
          auto a_itr( FWD_ENS.Begin()), a_itr_end( FWD_ENS.End());
          a_itr != a_itr_end;
          ++a_itr
      )
      {
        for
        (
            auto b_itr( REV_ENS.Begin()), b_itr_end( REV_ENS.End());
            b_itr != b_itr_end;
            ++b_itr
        )
        {
          // filter stupid distances
          if( !DistanceFilter( *a_itr, FWD_LINK_ATOM, *b_itr, REV_LINK_ATOM, m_JoinDistanceCutoffMin, m_JoinDistanceCutoffMax))
          {
            continue;
          }

          // open valence to enable linkage
          storage::Triplet< FragmentComplete, size_t, size_t> pair_a( OpenValence( *a_itr, FWD_LINK_ATOM, m_OVShuffleH, m_OVReverse));
          storage::Triplet< FragmentComplete, size_t, size_t> pair_b( OpenValence( *b_itr, REV_LINK_ATOM, m_OVShuffleH, m_OVReverse));

          // get atominfo and bondinfo, which we will append to shortly
          storage::Vector< sdf::AtomInfo> atominfo_a( pair_a.First().GetAtomVector().GetAtomInfo());
          storage::Vector< sdf::BondInfo> bondinfo_a( pair_a.First().GetAtomVector().GetBondInfo());
          storage::Vector< sdf::BondInfo> b_itr_bondinfo( pair_b.First().GetAtomVector().GetBondInfo());

          // start by adding the atominfos of the second fragment into the first fragment
          atominfo_a.Append( pair_b.First().GetAtomVector().GetAtomInfo());

          // do the same for the bondinfo, but update the indices to match the new molecule size
          for( size_t b_itr_bondinfo_i( 0); b_itr_bondinfo_i < b_itr_bondinfo.GetSize(); ++b_itr_bondinfo_i)
          {
            bondinfo_a.PushBack
            (
              sdf::BondInfo
              (
                b_itr_bondinfo( b_itr_bondinfo_i).GetAtomIndexLow() + pair_a.First().GetSize(),
                b_itr_bondinfo( b_itr_bondinfo_i).GetAtomIndexHigh() + pair_a.First().GetSize(),
                b_itr_bondinfo( b_itr_bondinfo_i).GetConfigurationalBondType()
              )
            );
          }

          // now join the two halves at the cutpoint with a single bond
          bondinfo_a.PushBack
          (
            sdf::BondInfo
            (
              pair_a.Second(),
              pair_b.Second() + pair_a.First().GetSize(),
              GetConfigurationalBondTypes().e_NonConjugatedSingleBond
            )
          );

          // create new molecule
          AtomVector< AtomComplete> new_mol_v( atominfo_a, bondinfo_a);
          FragmentComplete new_mol( new_mol_v, "");

          // add fragment_a linker indices
          if( !linker_atom_indices.GetSize())
          {
            linker_atom_indices = storage::Vector< size_t>( FWD_SBP_INDICES.Begin(), FWD_SBP_INDICES.End());
            linker_atom_indices_plus_connect = linker_atom_indices;

            // remove the atom in fragment_a to which the linker connects from the available linker indices
            linker_atom_indices.RemoveElements( 0);

            // collapse the mobile sbp atom indices from the reverse linker into the
            // total collection of linker atom indices; offset assigns correct atom index;
            // start from 1 so that we do not include the atom in fragment_b to which the reverse linker is connected
            for( size_t i( 0); i < REV_SBP_INDICES.GetSize(); ++i)
            {
              // i realize I could do this less stupidly, but i do not want the first rev index to be at the end of the vector;
              // when we fix bond lengths, we want to include the connection point
              linker_atom_indices_plus_connect.PushBack( REV_SBP_INDICES( i) + pair_a.First().GetSize());
              if( i)
              {
                // otherwise, we just want the linker
                linker_atom_indices.PushBack( REV_SBP_INDICES( i) + pair_a.First().GetSize());
              }
            }

            // get the hydrogen atoms bonded to the linker heavy atoms, too
            storage::Vector< size_t> linker_h_atom_indices;
            for( size_t l_i( 0), l_i_sz( linker_atom_indices.GetSize()); l_i < l_i_sz; ++l_i)
            {
              for
              (
                  auto bond_itr( new_mol.GetAtomVector()( linker_atom_indices( l_i)).GetBonds().Begin()),
                  bond_itr_end( new_mol.GetAtomVector()( linker_atom_indices( l_i)).GetBonds().End());
                  bond_itr != bond_itr_end;
                  ++bond_itr
              )
              {
                if( bond_itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Hydrogen)
                {
                  linker_h_atom_indices.PushBack( new_mol.GetAtomVector().GetAtomIndex( bond_itr->GetTargetAtom()));
                }
              }
            }
            for( size_t i( 0); i < linker_h_atom_indices.GetSize(); ++i)
            {
              linker_atom_indices.PushBack( linker_h_atom_indices( i));
            }
          }

          // correct bad bond lengths
//          CorrectBadBondLengths( new_mol, linker_atom_indices_plus_connect);
          CorrectBadBondLengths( new_mol, storage::Vector< size_t>());

          // try to massage the linker with mutate clash resolver
//          ResolveClashes(new_mol, linker_atom_indices, ROTLIB);

          // save
          final_conformer_ensemble.PushBack( new_mol);

        } // molecule a conformers
      } // molecule b conformers

      BCL_MessageStd( "========== END JoinHalfExtendedFragments ==========");
      return storage::Pair< FragmentEnsemble, storage::Vector< size_t> >( final_conformer_ensemble, linker_atom_indices);
    }

    //! @brief select rings from the library and identify attachment atoms
    //! @params LINKER_COMPONENTS the components to check for rings
    //! @return a collection of rings and corresponding attachment indices
    storage::Vector< storage::Triplet< FragmentComplete, size_t, size_t> > FragmentMutateConnect::ChooseRings
    (
      const storage::Vector< std::string> &LINKER_COMPONENTS
    ) const
    {
      // initialize output
      storage::Vector< storage::Triplet< FragmentComplete, size_t, size_t> > rings_and_indices( LINKER_COMPONENTS.GetSize());

      // check which linker components are rings
      size_t i( 0);
      for
      (
          auto itr( LINKER_COMPONENTS.Begin()), itr_end( LINKER_COMPONENTS.End());
          itr != itr_end;
          ++itr, ++i
      )
      {
        // increment undefined values for non-ring components
        if( *itr != "ring")
        {
          rings_and_indices( i) = storage::Triplet< FragmentComplete, size_t, size_t>();
          continue;
        }

        // if it is a ring then choose from our library
        FragmentEnsemble::const_iterator ring_itr( m_Rings->Begin());
        size_t rand_pos( random::GetGlobalRandom().Random< size_t>( size_t( 0), m_Rings->GetSize() - 1));
        std::advance( ring_itr, rand_pos);

        // and then choose the link atoms in the ring
        size_t i_a( RunPickAtom( *ring_itr, true, true));
        size_t i_b( RunPickAtom( *ring_itr, true, true));
        size_t i_tries( 1);
        while( i_a == i_b && i_tries < m_NumberMaxAttempts)
        {
          ++i_tries;
          i_b = RunPickAtom( *ring_itr, true, true);
        }
        BCL_MessageStd( "Chosen ring atom 0 for ring in linker index " + util::Format()( i) + ": " + util::Format()( i_a));
        BCL_MessageStd( "Chosen ring atom 1 for ring in linker index " + util::Format()( i) + ": " + util::Format()( i_b));

        // save
        rings_and_indices( i) = storage::Triplet< FragmentComplete, size_t, size_t>( *ring_itr, i_a, i_b);
      }
      return rings_and_indices;
    }

    //! @brief modify input molecule by adding a link fragment
    //! @param FRAGMENT the fragment to be modified
    //! @param INDEX the index of FRAGMENT to which the linker fragment will be appended
    //! @param LINK_TYPE the string indicating the link type
    FragmentComplete FragmentMutateConnect::AddLinkFragment
    (
      const FragmentComplete FRAGMENT,
      const size_t INDEX,
      const std::string &LINK_TYPE,
      const FragmentComplete &RING,
      const size_t ATTACH_INDEX
    ) const
    {
      FragmentComplete fragment( FRAGMENT);

      // TODO: collapse the single element linker conditionals into a single conditional that takes an element type;
      // may require refactor in ExtendWithLinker SingleElementLink function
      // extend by a single carbon atom
      if( LINK_TYPE == "B")
      {
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.SingleElementLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t(), // if connecting to empty fragment, this does not matter
                "B"
              )
            );
      }
      // extend by a single carbon atom
      else if( LINK_TYPE == "C")
      {
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.SingleElementLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t(), // if connecting to empty fragment, this does not matter
                "C"
              )
            );
      }
      // extend by a single oxygen atom
      else if( LINK_TYPE == "O")
      {
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.SingleElementLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t(), // if connecting to empty fragment, this does not matter
                "O"
              )
            );
      }
      // extend by a single nitrogen atom
      else if( LINK_TYPE == "N")
      {
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.SingleElementLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t(), // if connecting to empty fragment, this does not matter
                "N"
              )
            );
      }
      // extend by a single phosphorous atom
      else if( LINK_TYPE == "P")
      {
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.SingleElementLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t(), // if connecting to empty fragment, this does not matter
                "P"
              )
            );
      }
      // extend by a single sulfur atom
      else if( LINK_TYPE == "S")
      {
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.SingleElementLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t(), // if connecting to empty fragment, this does not matter
                "S"
              )
            );
      }
      // extend by a single selenium atom
      else if( LINK_TYPE == "Se")
      {
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.SingleElementLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t(), // if connecting to empty fragment, this does not matter
                "Se"
              )
            );
      }
      // extend by multiple carbon atoms simultaneously
      else if( LINK_TYPE == "alkyl")
      {
        m_ExtendWithLinker.SetAlkyneProb( 0);
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.AlkylLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t(), // if connecting to empty fragment, this does not matter
                size_t( 2)
              )
            );
      }
      // extend by an n-amide
      else if( LINK_TYPE == "n_amide")
      {
        m_ExtendWithLinker.SetAmideNToAProb( 1.0);
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.AmideLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t(), // if connecting to empty fragment, this does not matter
                size_t( 1)
              )
            );
      }
      // extend by an c-amide
      else if( LINK_TYPE == "c_amide")
      {
        m_ExtendWithLinker.SetAmideNToAProb( 0.0);
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.AmideLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t(), // if connecting to empty fragment, this does not matter
                size_t( 1)
              )
            );
      }
      // extend by an o-ester
      else if( LINK_TYPE == "o_ester")
      {
        m_ExtendWithLinker.SetEsterOToAProb( 1.0);
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.EsterLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t(), // if connecting to empty fragment, this does not matter
                size_t( 1)
              )
            );
      }
      // extend by an c-ester
      else if( LINK_TYPE == "c_ester")
      {
        m_ExtendWithLinker.SetEsterOToAProb( 0.0);
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.EsterLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t(), // if connecting to empty fragment, this does not matter
                size_t( 1)
              )
            );
      }
      // extend by a methoxy with the o linked to fragment A
      else if( LINK_TYPE == "oa_methoxy")
      {
        m_ExtendWithLinker.SetMethoxyOToAProb( 1.0);
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.MethoxyLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t(), // if connecting to empty fragment, this does not matter
                size_t( 1)
              )
            );
      }
      // extend by a methoxy with the o linked to fragment B
      else if( LINK_TYPE == "ob_methoxy")
      {
        m_ExtendWithLinker.SetMethoxyOToAProb( 0.0);
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.MethoxyLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t(), // if connecting to empty fragment, this does not matter
                size_t( 1)
              )
            );
      }
      // extend by a ethoxy with the o linked to fragment A
      else if( LINK_TYPE == "oa_ethoxy")
      {
        m_ExtendWithLinker.SetEthoxyOToAProb( 1.0);
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.EthoxyLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t(), // if connecting to empty fragment, this does not matter
                size_t( 1)
              )
            );
      }
      // extend by a ethoxy with the o linked to fragment B
      else if( LINK_TYPE == "ob_ethoxy")
      {
        m_ExtendWithLinker.SetEthoxyOToAProb( 0.0);
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.EthoxyLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t(), // if connecting to empty fragment, this does not matter
                size_t( 1)
              )
            );
      }
      // ring link
      else if( LINK_TYPE == "ring")
      {
        fragment = FragmentComplete
            (
              m_ExtendWithLinker.RingLink
              (
                fragment,
                FragmentComplete(),
                INDEX,
                size_t( 1),
                RING,
                ATTACH_INDEX
              )
            );
      }
      return fragment;
    }

    //! @brief create supermolecule of the two fragments without a linker
    FragmentComplete FragmentMutateConnect::CreateSuperMolecule
    (
      const FragmentComplete &FRAGMENT_A,
      const FragmentComplete &FRAGMENT_B
    ) const
    {
      // get some references to the atom vectors for convenience
      AtomVector< AtomComplete> atom_v_a( FRAGMENT_A.GetAtomVector());
      const AtomVector< AtomComplete> &atom_v_b( FRAGMENT_B.GetAtomVector());
      storage::Vector< sdf::BondInfo> empty_bonds( 0);
      atom_v_a.AddAtomsWithConnectivity( atom_v_b, empty_bonds);
      return FragmentComplete( atom_v_a, "");
    }

    //! @brief compute RMSD of linked molecule to the unlinked supermolecule of the two fragments
    double FragmentMutateConnect::ComputeRMSDToSuperMolecule
    (
      const FragmentComplete &MOLECULE,
      const storage::Vector< size_t> &COMMON_INDICES,
      const FragmentComplete &SUPERMOLECULE
    ) const
    {

      AtomVector< AtomComplete> mol_atom_v( MOLECULE.GetAtomVector());
      mol_atom_v.Reorder( COMMON_INDICES);
      FragmentComplete molecule( mol_atom_v, "");

      // get the heavy atom coordinates for both molecules
      util::SiPtrVector< const linal::Vector3D> mol_coords( molecule.GetHeavyAtomCoordinates());
      util::SiPtrVector< const linal::Vector3D> supermol_coords( SUPERMOLECULE.GetHeavyAtomCoordinates());

      // measure RMSD
      return quality::RMSD::RealSpaceRMSD( mol_coords, supermol_coords);
    }

    //! @brief filter out linker pairs whose link atoms are too far apart
    //! @param FRAGMENT_A first fragment to be linked
    //! @param DIST_ATOM_A atom at which FRAGMENT_A will be linked to FRAGMENT_B
    //! @param FRAGMENT_B second fragment to be linked
    //! @param DIST_ATOM_B atom at which FRAGMENT_B will be linked to FRAGMENT_A
    //! @param DIST_THRESHOLD maximum allowed distance between DIST_ATOM_A and DIST_ATOM_B
    //! @return true if equal or below DIST_THRESHOLD, false otherwise
    bool FragmentMutateConnect::DistanceFilter
    (
      const FragmentComplete &FRAGMENT_A,
      const size_t DIST_ATOM_A,
      const FragmentComplete &FRAGMENT_B,
      const size_t DIST_ATOM_B,
      const float MIN_DIST_THRESHOLD, // TODO consider making a function of CovalentBondRadius
      const float MAX_DIST_THRESHOLD // TODO consider making a function of CovalentBondRadius
    ) const
    {
      // get positions of atoms in molecules
      const linal::Vector3D &atom_a( FRAGMENT_A.GetAtomVector()( DIST_ATOM_A).GetPosition());
      const linal::Vector3D &atom_b( FRAGMENT_B.GetAtomVector()( DIST_ATOM_B).GetPosition());

      // check if distance is above threshold
      const double distance( linal::Distance( atom_a, atom_b));
      if( distance > MAX_DIST_THRESHOLD || distance < MIN_DIST_THRESHOLD)
      {
        return false;
      }

      // below threshold
      return true;
    }

    void FragmentMutateConnect::ResolveClashes
    (
      FragmentComplete &MOLECULE,
      const storage::Vector< size_t> &LINKER_ATOM_INDICES,
      const RotamerLibraryFile &ROTLIB
    ) const
    {

      storage::Set< size_t> flexible_atoms( LINKER_ATOM_INDICES.Begin(), LINKER_ATOM_INDICES.End());
      MOLECULE.GetStoredPropertiesNonConst().SetMDLProperty( "SampleByParts", flexible_atoms);
      MutateClashResolver geometry_fixer( 10);

      util::ShPtr< SearchFragmentLibraryFromTree> rotlib_searcher;
      rotlib_searcher = util::ShPtr< SearchFragmentLibraryFromTree>( new SearchFragmentLibraryFromTree( ROTLIB));
      util::ShPtrVector< BondAngleAssignment> bond_angle_assignments
      (
        rotlib_searcher->GetBondAngleAssignments( MOLECULE, true)
      );

      util::ShPtrList< MutateDihedralsInterface> mutate_bond_angles_refinement;
      MutateChirality chirality_changer( MOLECULE);
      for( auto itr_bas( bond_angle_assignments.Begin()), itr_bas_end( bond_angle_assignments.End()); itr_bas != itr_bas_end; ++itr_bas)
      {
        mutate_bond_angles_refinement.PushBack
        (
          util::ShPtr< MutateBondAngles>
          (
            new MutateBondAngles
            (
              *itr_bas,
              MOLECULE,
              chirality_changer,
              false,
              false,
              true
            )
          )
        );
      }

      util::ShPtr< AtomClashScore> clash_score_sp( new AtomClashScore( true));
      geometry_fixer.Setup
      (
        MOLECULE,
        util::ShPtrVector< MutateBondAngles>( mutate_bond_angles_refinement.Begin(), mutate_bond_angles_refinement.End()),
        clash_score_sp
      );

      auto fixed_mol( geometry_fixer( MOLECULE));
      MOLECULE = *( fixed_mol.GetArgument());
    }

    void FragmentMutateConnect::CorrectBadBondLengths
    (
      FragmentComplete &MOLECULE,
      const storage::Vector< size_t> &LINKER_ATOM_INDICES
    ) const
    {
      MutateBondLengths bond_length_mutate( 100, 0.01, 0.01, LINKER_ATOM_INDICES);
      bond_length_mutate.Initialize( MOLECULE);

      // run bond length mutate
      auto fixed_frag( bond_length_mutate( MOLECULE));
      if( fixed_frag.GetArgument().IsDefined())
      {
        // contains at least one bad bond length
        MOLECULE = FragmentComplete( *( fixed_frag.GetArgument()));
      }
      MOLECULE.SaturateWithH();
    }

  //////////////////////
  // helper functions //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    io::Serializer FragmentMutateConnect::GetSerializer() const
    {
      io::Serializer parameters( FragmentMutateInterface::GetSerializer());
      parameters.SetClassDescription
      (
        "Connects two disconnected fragments using linkers available in FragmentExtendWithLinker.\n"
        "The connection algorithm proceeds in two primary pieces:\n(1) The user-specified linker string "
        "is interpreted and attached to each of the two fragments being connected in halves, such that "
        "both fragments contain half of the total linker; conformational ensembles are generated for "
        "these half linkers using the 'global' settings.\n(2) Conformer pairs of the two ensembles are "
        "chosen based on terminal atom distance/angle criteria, linked, and then refined through 'local' "
        "conformer ensemble generation of just the linker region. The final structure is chosen as a "
        "compromise between the best conformer by ConfScore and the nearest conformer by RMSD to the "
        "starting fragments. 'Global' and 'local' refer to the SetSamplingPreferences options in the "
        "SampleConformations class, which dictate whether sampling will include dihedrals, rings, bond "
        "angles, bond lengths, and/or chirality."
      );

      parameters.AddInitializer
      (
        "terminal_fragment",
        "molecule to which input fragment will be connected",
        io::Serialization::GetAgent( &m_TerminalFragmentFilename),
        ""
      );

      parameters.AddInitializer
      (
        "join_min_distance_cutoff",
        "if the distance between two half-linker terminal heavy atoms is less than "
        "this value then that pair is rejected",
        io::Serialization::GetAgent( &m_JoinDistanceCutoffMin),
        "1.0"
      );

      parameters.AddInitializer
      (
        "join_max_distance_cutoff",
        "if the distance between two half-linker terminal heavy atoms is greater than "
        "this value then that pair is rejected",
        io::Serialization::GetAgent( &m_JoinDistanceCutoffMax),
        "2.5"
      );

      parameters.AddInitializer
      (
        "output_ensemble",
        "if a name is provided, write a conformational ensemble of the linked molecule; no effect "
        "if an output filename is not provided",
        io::Serialization::GetAgent( &m_OutputEnsembleFilename),
        ""
      );

      parameters.AddInitializer
      (
        "rmsd_cutoff",
        "starting fragment rmsd cutoff below which conformers can be added to the optional "
        "output ensemble",
        io::Serialization::GetAgent( &m_StartingFragmentRMSDCutoff),
        "0.75"
      );

      parameters.AddInitializer
      (
        "global_conformer_comparer",
        "method to compare conformers with",
        io::Serialization::GetAgent( &m_ConfComparerGlobal),
        "SymmetryRMSD"
      );

      parameters.AddInitializer
      (
        "global_conformer_tolerance",
        "amount of tolerance allowed between two conformers",
        io::Serialization::GetAgent( &m_ConfToleranceGlobal),
        "0.25"
      );

      parameters.AddInitializer
      (
        "global_max_conformations",
        "maximum number conformations to generate per molecule",
        io::Serialization::GetAgent( &m_NMaxConfsGlobal),
        "250"
      );

      parameters.AddInitializer
      (
        "global_conformer_max_iterations",
        "number of iterations",
        io::Serialization::GetAgent( &m_NMaxItersGlobal),
        "2000"
      );

      parameters.AddInitializer
      (
        "global_conformer_cluster",
        "perform leader-follower clustering during conformer generation",
        io::Serialization::GetAgent( &m_ClusterGlobal),
        "true"
      );

      parameters.AddInitializer
      (
        "global_max_cr_iterations",
        "max number of clash resolution iterations performed (as fraction of # dihedrals & bond angles)."
        "More iterations will produced less clashed molecules but may "
        "take much longer. Often it is faster to just make more molecules and discard the excessively "
        "clashed ones if they can't be resolved in a small fraction of moves (e.g. 10% of the total dihedrals + bond angles",
        io::Serialization::GetAgent( &m_ClashResolutionGlobal),
        "0.1"
      );

      parameters.AddInitializer
      (
        "global_max_avg_clash",
        "max average clash allowed for conformations =sum(vdw_overlap for all atoms in molecule) / atoms in molecule."
        " larger values will yield potentially more clashed conformers for some molecules",
        io::Serialization::GetAgent( &m_ClashScoreGlobal),
        "0.1"
      );

      parameters.AddInitializer
      (
        "local_conformer_comparer",
        "method to compare conformers with",
        io::Serialization::GetAgent( &m_ConfComparerLocal),
        ""
      );

      parameters.AddInitializer
      (
        "local_conformer_tolerance",
        "amount of tolerance allowed between two conformers",
        io::Serialization::GetAgent( &m_ConfToleranceLocal),
        "0.125"
      );

      parameters.AddInitializer
      (
        "local_max_conformations",
        "maximum number conformations to generate per molecule",
        io::Serialization::GetAgent( &m_NMaxConfsLocal),
        "50"
      );

      parameters.AddInitializer
      (
        "local_conformer_max_iterations",
        "number of iterations",
        io::Serialization::GetAgent( &m_NMaxItersLocal),
        "500"
      );

      parameters.AddInitializer
      (
        "local_conformer_cluster",
        "perform leader-follower clustering during conformer generation",
        io::Serialization::GetAgent( &m_ClusterLocal),
        "false"
      );

      parameters.AddInitializer
      (
        "local_max_cr_iterations",
        "max number of clash resolution iterations performed (as fraction of # dihedrals & bond angles)."
        "More iterations will produced less clashed molecules but may "
        "take much longer. Often it is faster to just make more molecules and discard the excessively "
        "clashed ones if they can't be resolved in a small fraction of moves (e.g. 10% of the total dihedrals + bond angles",
        io::Serialization::GetAgent( &m_ClashResolutionLocal),
        "0.5"
      );

      parameters.AddInitializer
      (
        "local_max_avg_clash",
        "max average clash allowed for conformations =sum(vdw_overlap for all atoms in molecule) / atoms in molecule."
        " larger values will yield potentially more clashed conformers for some molecules",
        io::Serialization::GetAgent( &m_ClashScoreLocal),
        "0.1"
      );

      parameters.AddInitializer
      (
        "ring_library",
        "path to the ring library",
        io::Serialization::GetAgent( &m_RingsFilename),
        command::CommandState::IsInStaticInitialization() ?
            "" :
            RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "ring_libraries/drug_ring_database.simple.aro.sdf.gz"
      );

      parameters.AddInitializer
      (
        "terminal_fragment_mutable_atoms",
        "atom indices (0-indexed) that can be mutated",
        io::Serialization::GetAgent( &m_TerminalFragmentMutableAtoms),
        ""
      );

      parameters.AddInitializer
      (
        "terminal_fragment_mutable_elements",
        "element types that can be mutated",
        io::Serialization::GetAgent( &m_TerminalFragmentMutableElementsString),
        ""
      );

      parameters.AddInitializer
      (
        "terminal_fragment_mutable_fragments",
        "allow atoms that match these fragments in substructure comparisons to be "
        "mutated; can be inverted to set mutable atoms to the complement subgraph "
        "isomorphism atoms via 'complement_mutable_fragments'",
        io::Serialization::GetAgent( &m_TerminalFragmentMutableFragmentsFilename),
        ""
      );

      parameters.AddInitializer
      (
        "terminal_fragment_complement_mutable_fragments",
        "invert the subgraph isomorphisms between the molecule of interest and "
        "the mutable fragments such that the derived mutable atoms are the "
        "non-common atoms",
        io::Serialization::GetAgent( &m_TerminalFragmentComplementMutableFragments),
        "false"
      );

      parameters.AddInitializer
      (
        "terminal_fragment_mutable_atom_comparison",
        "atom data that are compared to determine whether atoms are equivalent "
        "during mutable fragment substructure comparisons",
        io::Serialization::GetAgent( &m_TerminalFragmentMutableAtomComparisonType),
        "ElementType"
      );

      parameters.AddInitializer
      (
        "terminal_fragment_mutable_bond_comparison",
        "bond data that are compared to determine whether bonds are equivalent "
        "during mutable fragment substructure comparisons",
        io::Serialization::GetAgent( &m_TerminalFragmentMutableBondComparisonType),
        "BondOrderAmideOrAromaticWithRingness"
      );

      parameters.AddInitializer
      (
        "terminal_fragment_fixed_atoms",
        "atom indices (0-indexed) that will be fixed after mutable atoms are assigned",
        io::Serialization::GetAgent( &m_TerminalFragmentFixedAtoms),
        ""
      );

      parameters.AddInitializer
      (
        "terminal_fragment_fixed_elements",
        "element types that will be fixed after mutable atoms are assigned",
        io::Serialization::GetAgent( &m_TerminalFragmentFixedElementsString),
        ""
      );

      parameters.AddInitializer
      (
        "terminal_fragment_fixed_fragments",
        "fix atoms that match these fragments in substructure comparisons; "
        "can be inverted to set fixed atoms to the complement subgraph "
        "isomorphism atoms via 'complement_fixed_fragments'",
        io::Serialization::GetAgent( &m_TerminalFragmentFixedFragmentsFilename),
        ""
      );

      parameters.AddInitializer
      (
        "terminal_fragment_complement_fixed_fragments",
        "invert the subgraph isomorphisms between the molecule of interest and "
        "the fixed fragments such that the derived mutable atoms are the "
        "non-common atoms",
        io::Serialization::GetAgent( &m_TerminalFragmentComplementFixedFragments),
        "false"
      );

      parameters.AddInitializer
      (
        "terminal_fragment_fixed_atom_comparison",
        "atom data that are compared to determine whether atoms are equivalent "
        "during fixed fragment substructure comparisons",
        io::Serialization::GetAgent( &m_TerminalFragmentFixedAtomComparisonType),
        "ElementType"
      );

      parameters.AddInitializer
      (
        "terminal_fragment_fixed_bond_comparison",
        "bond data that are compared to determine whether bonds are equivalent "
        "during fixed fragment substructure comparisons",
        io::Serialization::GetAgent( &m_TerminalFragmentFixedBondComparisonType),
        "BondOrderAmideOrAromaticWithRingness"
      );

      parameters.AddInitializer
      (
        "linker_composition",
        "specify the linker types that will be used to connect the specified fragments; "
        "the linker subunits will be added sequentially as they are given on the command-line; "
        "the full linker is cut in half by naively counting the subunits specified in the string, "
        "such that odd numbers will favor an extra subunit on the input fragment; "
        "the first half-linker will be connected to the input fragment and the second half-linker "
        "will be inverted and connected to the 'terminal_fragment' to yield the complete linker in "
        "the correct sequence as specified here; closure happens internally in the linker."
        "Linker types are as follows: "
        "B - ; "
        "C - ; "
        "O - ; "
        "N - ; "
        "P - ; "
        "S - ; "
        "Se - ; "
        "n_amide - ; "
        "c_amide - ; "
        "o_ester - ; "
        "c_ester - ; "
        "oa_methoxy - ; "
        "ob_methoxy - ; "
        "oa_ethoxy - ; "
        "ob_ethoxy - ; "
        "ring - ; ",
        io::Serialization::GetAgent( &m_LinkersString),
        "N C C O C C O C C O C C O C C c_ester"
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentMutateConnect::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
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

      // read in terminal fragment filename
      if( m_TerminalFragmentFilename.size())
      {
        // lock to avoid multiple thread read/write issues
        s_Mutex.Lock();
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_TerminalFragmentFilename);
        FragmentEnsemble terminal_mol;

        // shuffle the ensemble
        terminal_mol.ReadMoreFromMdl( input);
        terminal_mol.Shuffle();

        // choose random
        m_TerminalFragment = terminal_mol.GetMolecules().FirstElement();

        // close and unlock
        io::File::CloseClearFStream( input);
        s_Mutex.Unlock();
      }
      else
      {
        // mandatory
        BCL_MessageStd( "No molecule provided for the terminal fragment. Check terminal fragment file.");
        return false;
      }

      // read in ring library
      if( m_RingsFilename.size())
      {
        s_Mutex.Lock();
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_RingsFilename);
        FragmentEnsemble ensemble( input, sdf::e_Remove);
        m_Rings = util::CloneToShPtr( ensemble);
        io::File::CloseClearFStream( input);
        s_Mutex.Unlock();
        BCL_Assert( m_Rings->GetSize(), "Ring library is empty!");
        s_Mutex.Lock();
        s_RingLibraries.Insert( std::make_pair( m_RingsFilename, m_Rings));
        s_Mutex.Unlock();
      }

      // read in reference filename
      s_Mutex.Lock();
      if( m_TerminalFragmentMutableFragmentsFilename.size())
      {
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_TerminalFragmentMutableFragmentsFilename);
        m_TerminalFragmentMutableFragments.ReadMoreFromMdl( input);
        io::File::CloseClearFStream( input);
      }
      s_Mutex.Unlock();

      // read in mutable atom indices
      if( m_TerminalFragmentMutableAtoms.size())
      {
        m_TerminalFragmentMutableAtomIndices.Reset();
        m_TerminalFragmentMutableAtomIndices = util::SplitStringToNumerical< size_t>( m_TerminalFragmentMutableAtoms);
      }

      // read in allowed mutable elements when choosing mutable atom
      if( m_TerminalFragmentMutableElementsString.size())
      {
        // parse input
        const storage::Vector< std::string> mutable_elements
        (
          util::SplitString( util::TrimString( m_TerminalFragmentMutableElementsString), " \t\n\r,")
        );

        // stupid check to add only the correct elements
        // TODO: this should be directly serializable from element types to make concise
        m_TerminalFragmentMutableElements.Reset();
        for( size_t e_i( 0), e_sz( mutable_elements.GetSize()); e_i < e_sz; ++e_i)
        {
          // Hydrogen
          if( mutable_elements( e_i) == "H")
          {
            m_TerminalFragmentMutableElements.PushBack( GetElementTypes().e_Hydrogen);
          }
          // Boron
          if( mutable_elements( e_i) == "B")
          {
            m_TerminalFragmentMutableElements.PushBack( GetElementTypes().e_Boron);
          }
          // Carbon
          if( mutable_elements( e_i) == "C")
          {
            m_TerminalFragmentMutableElements.PushBack( GetElementTypes().e_Carbon);
          }
          // Oxygen
          if( mutable_elements( e_i) == "O")
          {
            m_TerminalFragmentMutableElements.PushBack( GetElementTypes().e_Oxygen);
          }
          // Nitrogen
          if( mutable_elements( e_i) == "N")
          {
            m_TerminalFragmentMutableElements.PushBack( GetElementTypes().e_Nitrogen);
          }
          // Potassium
          if( mutable_elements( e_i) == "P")
          {
            m_TerminalFragmentMutableElements.PushBack( GetElementTypes().e_Phosphorus);
          }
          // Sulfur
          if( mutable_elements( e_i) == "S")
          {
            m_TerminalFragmentMutableElements.PushBack( GetElementTypes().e_Sulfur);
          }
          // Selenium
          if( mutable_elements( e_i) == "Se")
          {
            m_TerminalFragmentMutableElements.PushBack( GetElementTypes().e_Selenium);
          }
          // Fluorine
          if( mutable_elements( e_i) == "F")
          {
            m_TerminalFragmentMutableElements.PushBack( GetElementTypes().e_Fluorine);
          }
          // Chlorine
          if( mutable_elements( e_i) == "Cl")
          {
            m_TerminalFragmentMutableElements.PushBack( GetElementTypes().e_Chlorine);
          }
          // Bromine
          if( mutable_elements( e_i) == "Br")
          {
            m_TerminalFragmentMutableElements.PushBack( GetElementTypes().e_Bromine);
          }
          // Iodine
          if( mutable_elements( e_i) == "I")
          {
            m_TerminalFragmentMutableElements.PushBack( GetElementTypes().e_Iodine);
          }
          // Undefined
          if( mutable_elements( e_i) == "X")
          {
            m_TerminalFragmentMutableElements.PushBack( GetElementTypes().e_Undefined);
          }
        }
      }

      // read in reference filename
      s_Mutex.Lock();
      if( m_TerminalFragmentFixedFragmentsFilename.size())
      {
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_TerminalFragmentFixedFragmentsFilename);
        m_TerminalFragmentFixedFragments.ReadMoreFromMdl( input);
        io::File::CloseClearFStream( input);
      }
      s_Mutex.Unlock();

      // read in mutable atom indices
      if( m_TerminalFragmentFixedAtoms.size())
      {
        m_TerminalFragmentFixedAtomindices.Reset();
        m_TerminalFragmentFixedAtomindices = util::SplitStringToNumerical< size_t>( m_TerminalFragmentFixedAtoms);
      }

      // read in allowed mutable elements when choosing mutable atom
      if( m_TerminalFragmentFixedElementsString.size())
      {
        // parse input
        const storage::Vector< std::string> fixed_elements
        (
          util::SplitString( util::TrimString( m_TerminalFragmentFixedElementsString), " \t\n\r,")
        );

        // check to add only the correct elements
        m_TerminalFragmentFixedElements.Reset();
        for( size_t e_i( 0), e_sz( fixed_elements.GetSize()); e_i < e_sz; ++e_i)
        {
          // Hydrogen
          if( fixed_elements( e_i) == "H")
          {
            m_TerminalFragmentFixedElements.PushBack( GetElementTypes().e_Hydrogen);
          }
          // Boron
          if( fixed_elements( e_i) == "B")
          {
            m_TerminalFragmentFixedElements.PushBack( GetElementTypes().e_Boron);
          }
          // Carbon
          if( fixed_elements( e_i) == "C")
          {
            m_TerminalFragmentFixedElements.PushBack( GetElementTypes().e_Carbon);
          }
          // Oxygen
          if( fixed_elements( e_i) == "O")
          {
            m_TerminalFragmentFixedElements.PushBack( GetElementTypes().e_Oxygen);
          }
          // Nitrogen
          if( fixed_elements( e_i) == "N")
          {
            m_TerminalFragmentFixedElements.PushBack( GetElementTypes().e_Nitrogen);
          }
          // Potassium
          if( fixed_elements( e_i) == "P")
          {
            m_TerminalFragmentFixedElements.PushBack( GetElementTypes().e_Phosphorus);
          }
          // Sulfur
          if( fixed_elements( e_i) == "S")
          {
            m_TerminalFragmentFixedElements.PushBack( GetElementTypes().e_Sulfur);
          }
          // Selenium
          if( fixed_elements( e_i) == "Se")
          {
            m_TerminalFragmentFixedElements.PushBack( GetElementTypes().e_Selenium);
          }
          // Fluorine
          if( fixed_elements( e_i) == "F")
          {
            m_TerminalFragmentFixedElements.PushBack( GetElementTypes().e_Fluorine);
          }
          // Chlorine
          if( fixed_elements( e_i) == "Cl")
          {
            m_TerminalFragmentFixedElements.PushBack( GetElementTypes().e_Chlorine);
          }
          // Bromine
          if( fixed_elements( e_i) == "Br")
          {
            m_TerminalFragmentFixedElements.PushBack( GetElementTypes().e_Bromine);
          }
          // Iodine
          if( fixed_elements( e_i) == "I")
          {
            m_TerminalFragmentFixedElements.PushBack( GetElementTypes().e_Iodine);
          }
          // Undefined
          if( fixed_elements( e_i) == "X")
          {
            m_TerminalFragmentFixedElements.PushBack( GetElementTypes().e_Undefined);
          }
        }
      }

      // read in mutable atom indices
      if( m_TerminalFragmentMutableAtoms.size())
      {
        m_TerminalFragmentMutableAtomIndices = util::SplitStringToNumerical< size_t>( m_TerminalFragmentMutableAtoms);
      }

      // read in linker types
      if( m_LinkersString.size())
      {
        // parse input
        m_Linkers = storage::Vector< std::string>
        (
          util::SplitString( util::TrimString( m_LinkersString), " \t\n\r,")
        );
      }

      // done
      return true;
    }

  } // namespace chemistry
} // namespace bcl
