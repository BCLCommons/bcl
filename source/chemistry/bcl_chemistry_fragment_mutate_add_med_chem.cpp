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
#include "chemistry/bcl_chemistry_fragment_mutate_add_med_chem.h"

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
    const util::SiPtr< const util::ObjectInterface> FragmentMutateAddMedChem::s_Instance
    (
      util::Enumerated< FragmentMutateInterface>::AddInstance( new FragmentMutateAddMedChem())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentMutateAddMedChem::FragmentMutateAddMedChem() :
      m_FragmentPool( util::ShPtr< FragmentEnsemble>()),
      m_MedChemFilename( std::string()),
      m_RestrictAdditionsToAroRings( false)
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief construct with a pool of external fragments for fragment grow
    //! @param FRAGMENT_POOL external fragments to add to base fragment
    FragmentMutateAddMedChem::FragmentMutateAddMedChem
    (
      const util::ShPtr< FragmentEnsemble> &FRAGMENT_POOL,
      const bool &CORINA_CONFS
    ) :
      m_FragmentPool( util::ShPtr< FragmentEnsemble>()),
      m_MedChemFilename( std::string()),
      m_RestrictAdditionsToAroRings( false)
    {
      m_Corina = CORINA_CONFS;

      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief druglikeness constructor
    //! @param FRAGMENT_POOL external fragments to add to base fragment
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    FragmentMutateAddMedChem::FragmentMutateAddMedChem
    (
      const util::ShPtr< FragmentEnsemble> &FRAGMENT_POOL,
      const std::string &DRUG_LIKENESS_TYPE,
      const bool &CORINA_CONFS
    ) :
      m_FragmentPool( FRAGMENT_POOL),
      m_MedChemFilename( std::string()),
      m_RestrictAdditionsToAroRings( false)
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
    FragmentMutateAddMedChem::FragmentMutateAddMedChem
    (
      const util::ShPtr< FragmentEnsemble> &FRAGMENT_POOL,
      const std::string &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const bool &CORINA_CONFS
    ) :
      m_FragmentPool( FRAGMENT_POOL),
      m_MedChemFilename( std::string()),
      m_RestrictAdditionsToAroRings( false)
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
    FragmentMutateAddMedChem::FragmentMutateAddMedChem
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
      m_MedChemFilename( std::string()),
      m_RestrictAdditionsToAroRings( false)
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
    FragmentMutateAddMedChem::FragmentMutateAddMedChem
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
      m_MedChemFilename( std::string()),
      m_RestrictAdditionsToAroRings( false)
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
    FragmentMutateAddMedChem *FragmentMutateAddMedChem::Clone() const
    {
      return new FragmentMutateAddMedChem( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentMutateAddMedChem::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentMutateAddMedChem::GetAlias() const
    {
      static const std::string s_name( "AddMedChem");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @return MutateResult with Constitution after the mutate
    math::MutateResult< FragmentComplete> FragmentMutateAddMedChem::operator()( const FragmentComplete &FRAGMENT) const
    {
      // mutate label
      BCL_MessageStd( "AddMedChem!");

      // redo the whole thing n-max times; increment can also be made in an inner while-loop during atom index selection
      size_t try_index( 0);
      for( ; try_index < m_NumberMaxAttempts; ++try_index)
      {

        // select random medchem fragment
        iterate::Generic< const FragmentComplete> itr_gen( m_FragmentPool->Begin(), m_FragmentPool->End());
        itr_gen.GotoRandomPosition();
        FragmentComplete medchem_frag( *itr_gen);

        // fragment pool marks "reactive" atom as the undefined atom
        // get the undefined atom index
        AtomVector< AtomComplete> medchem_atom_v( medchem_frag.GetAtomVector());
        size_t undefined_index( util::GetUndefinedSize_t());
        for( size_t i( 0), end_i( medchem_atom_v.GetSize()); i < end_i; ++i)
        {
          if( medchem_atom_v( i).GetElementType() == GetElementTypes().e_Undefined)
          {
            undefined_index = i;
          }
        }
        if( undefined_index == util::GetUndefinedSize_t())
        {
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }

        // get the element to which the reactive atom is bonded
        auto &atom_bonded_to_undefined( medchem_atom_v( undefined_index).GetBonds().Begin()->GetTargetAtom());
        size_t atom_bonded_to_undefined_index( medchem_atom_v.GetAtomIndex( atom_bonded_to_undefined));

        // now remove undefined atom and update tracking of the bonded atom's index
        if( undefined_index < atom_bonded_to_undefined_index)
        {
          atom_bonded_to_undefined_index -= size_t( 1);
        }
        storage::Vector< size_t> defined_indices;
        for( size_t i( 0), end_i( medchem_atom_v.GetSize()); i < end_i; ++i)
        {
          if( medchem_atom_v( i).GetElementType() != GetElementTypes().e_Undefined)
          {
            defined_indices.PushBack( i);
          }
        }

        // create new medchem fragment that does not contain the undefined atom
        medchem_atom_v.Reorder( defined_indices);

        // clean up the medchem fragment
        AtomsCompleteStandardizer standardizer(medchem_atom_v,"",true);
        standardizer.SetConjugationOfBondTypes(medchem_atom_v);

        // make new fragment
        FragmentComplete new_medchem_frag( medchem_atom_v, "");
        new_medchem_frag.StandardizeBondLengths();

        // now pick a random heavy atom from the base fragment
        util::SiPtr< const AtomConformationalInterface> picked_atom;
        if( m_MutableAtomIndices.GetSize() || m_MutableElements.GetSize() || m_MutableFragments.GetSize())
        {
          picked_atom = this->PickAtom( FRAGMENT, false);
        }
        else
        {
          picked_atom = this->PickAtom( FRAGMENT, true);
        }

        // if the chosen atom is undefined then just grab a bonded atom
        // this is biased to the lower index bonded atom, but should not generally matter
        size_t undefined_base_index( util::GetUndefinedSize_t());
        if( picked_atom->GetElementType() == GetElementTypes().e_Undefined)
        {
          undefined_base_index = FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom);
          picked_atom = util::SiPtr< const AtomConformationalInterface>( picked_atom->GetBonds().Begin()->GetTargetAtom());
        }

        // restrict medchem additions to aromatic rings if desired to see if it is actually aromatic
        size_t picked_atom_index( FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom));
        bool aromatic( false);
        if
        (
            m_RestrictAdditionsToAroRings &&
            picked_atom->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1))
        )
        {
          for
          (
              auto bonds_itr( FRAGMENT.GetAtomVector()( picked_atom_index).GetBonds().Begin()),
              bonds_itr_end( FRAGMENT.GetAtomVector()( picked_atom_index).GetBonds().End());
              bonds_itr != bonds_itr_end;
              ++bonds_itr
          )
          {
            if( bonds_itr->GetBondType()->GetConjugation() == ConstitutionalBondTypeData::e_Aromatic)
            {
              aromatic = true;
              break;
            }
          }
        }

        // make sure it is a heavy atom
        while( picked_atom->GetElementType() == GetElementTypes().e_Hydrogen && try_index < m_NumberMaxAttempts)
        {
          // pick random atom to transform
          if( m_MutableAtomIndices.GetSize() || m_MutableElements.GetSize() || m_MutableFragments.GetSize())
          {
            picked_atom = this->PickAtom( FRAGMENT, false);
          }
          else
          {
            picked_atom = this->PickAtom( FRAGMENT, true);
          }

          if( picked_atom->GetElementType() == GetElementTypes().e_Undefined)
          {
            undefined_base_index = FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom);
            picked_atom = util::SiPtr< const AtomConformationalInterface>( picked_atom->GetBonds().Begin()->GetTargetAtom());
          }

          // restrict medchem additions to aromatic rings if desired to see if it is actually aromatic
          picked_atom_index = FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom);
          aromatic = false;
          if
          (
              m_RestrictAdditionsToAroRings &&
              picked_atom->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1))
          )
          {
            for
            (
                auto bonds_itr( FRAGMENT.GetAtomVector()( picked_atom_index).GetBonds().Begin()),
                bonds_itr_end( FRAGMENT.GetAtomVector()( picked_atom_index).GetBonds().End());
                bonds_itr != bonds_itr_end;
                ++bonds_itr
            )
            {
              if( bonds_itr->GetBondType()->GetConjugation() == ConstitutionalBondTypeData::e_Aromatic)
              {
                aromatic = true;
                break;
              }
            }
          }
          // update
          ++try_index;
        }

        // enforce heavy atom requirement
        if( picked_atom->GetElementType() == GetElementTypes().e_Hydrogen)
        {
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }

        // enforce aromatic ring requirement
        if
        (
            m_RestrictAdditionsToAroRings &&
            !aromatic
        )
        {
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }

        // move forward with the chosen index
        picked_atom_index = FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom);

        // adjust the picked atom index for when we remove the undefined atom from the base fragment
        if( undefined_base_index < picked_atom_index)
        {
          picked_atom_index -= size_t( 1);
        }

        // get an atom vector of starting molecule that only contains the defined atoms
        storage::Vector< size_t> keep_indices;
        for( size_t i( 0); i < FRAGMENT.GetSize(); ++i)
        {
          if( i == undefined_base_index)
          {
            continue;
          }
          keep_indices.PushBack( i);
        }
        AtomVector< AtomComplete> all_defined_atom_v( FRAGMENT.GetAtomVector());
        all_defined_atom_v.Reorder( keep_indices);
        AtomsCompleteStandardizer standardizer_2( all_defined_atom_v, "", true);
        standardizer.SetConjugationOfBondTypes( all_defined_atom_v);

        // make new fragment
        FragmentComplete fragment( all_defined_atom_v, FRAGMENT.GetName());
        fragment.StandardizeBondLengths();

        // TODO: control with bool
        //      // if aromatic, check if substitution will be directed correctly
        //      size_t n_aro_neigh( picked_atom->CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsAromatic, 1));
        //      if( n_aro_neigh)
        //      {
        //        if( !IsRingSubstitutionDirected( FRAGMENT, picked_atom))
        //        {
        //          // if aromatic and not correctly substituted, skip mutate
        //          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        //        }
        //      }

        // open a valence for the addition
        storage::Triplet< FragmentComplete, size_t, size_t> pair_a( OpenValence( fragment, picked_atom_index, m_OVShuffleH, m_OVReverse));

        // join the medchem fragment to our base molecule
        storage::Pair< bool, FragmentComplete> new_fragment
        (
          MergeFragmentComplete::MergeFragments
          (
            util::IsDefined( undefined_base_index) ? fragment : pair_a.First(),
                new_medchem_frag,
                GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
                util::IsDefined( undefined_base_index) ?
                    storage::Pair< size_t, size_t>( picked_atom_index, atom_bonded_to_undefined_index) :
                    storage::Pair< size_t, size_t>( pair_a.Second(), atom_bonded_to_undefined_index)
          )
        );

        // try again
        if( !new_fragment.First())
        {
          continue;
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
        AtomVector< AtomComplete> atoms( new_fragment.Second().GetAtomVector());

        // Remove hydrogen atoms to allow bond type adjustment
        HydrogensHandler::Remove( atoms);
        if( m_ScaffoldFragment.GetSize())
        {
          return math::MutateResult< FragmentComplete>( cleaner.Clean( atoms, m_ScaffoldFragment, m_DrugLikenessType), *this);
        }
        else
        {
          return math::MutateResult< FragmentComplete>( cleaner.Clean( atoms, fragment, m_DrugLikenessType), *this);
        }
      }
      return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set medchem fragment library from filename
    void FragmentMutateAddMedChem::SetFragmentLibraryFromFilename( const std::string &FRAGMENTS_FILENAME)
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

    io::Serializer FragmentMutateAddMedChem::GetSerializer() const
    {
      io::Serializer parameters( FragmentMutateInterface::GetSerializer());
      parameters.SetClassDescription
      (
        "Appends a classic medicinal chemistry functional group to the current molecule"
      );

      parameters.AddInitializer
      (
        "medchem_library",
        "path to the medchem library",
        io::Serialization::GetAgent( &m_MedChemFilename),
        RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "medchem_fragments/bcl_buildfrag_0.sdf.gz"
      );

      parameters.AddInitializer
      (
        "restrict_additions_to_aromatic_rings",
        "only permit medchem additions if the base atom is part of an aromatic ring",
        io::Serialization::GetAgent( &m_RestrictAdditionsToAroRings),
        "false"
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentMutateAddMedChem::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
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
      if( m_MedChemFilename.size())
      {
        s_Mutex.Lock();
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_MedChemFilename);
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
    std::istream &FragmentMutateAddMedChem::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_FragmentPool, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT number of indentations
    //! @return ostream which was written to
    std::ostream &FragmentMutateAddMedChem::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_FragmentPool, OSTREAM, INDENT) << '\n';
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
