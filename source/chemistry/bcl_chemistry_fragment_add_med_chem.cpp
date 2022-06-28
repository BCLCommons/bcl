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
#include "chemistry/bcl_chemistry_fragment_add_med_chem.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_element_types.h"
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
    const util::SiPtr< const util::ObjectInterface> FragmentAddMedChem::s_Instance
    (
      util::Enumerated< FragmentMutateInterface>::AddInstance( new FragmentAddMedChem())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentAddMedChem::FragmentAddMedChem() :
      m_FragmentPool( util::ShPtr< FragmentEnsemble>()),
      m_MedChemFilename( std::string()),
      m_RestrictAdditionsToAroRings( false)
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief construct with a pool of external fragments for fragment grow
    //! @param FRAGMENT_POOL external fragments to add to base fragment
    FragmentAddMedChem::FragmentAddMedChem
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
    FragmentAddMedChem::FragmentAddMedChem
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
    FragmentAddMedChem::FragmentAddMedChem
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
    FragmentAddMedChem::FragmentAddMedChem
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
    FragmentAddMedChem::FragmentAddMedChem
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
    FragmentAddMedChem *FragmentAddMedChem::Clone() const
    {
      return new FragmentAddMedChem( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentAddMedChem::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentAddMedChem::GetAlias() const
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
    math::MutateResult< FragmentComplete> FragmentAddMedChem::operator()( const FragmentComplete &FRAGMENT) const
    {
      // mutate label
      BCL_MessageStd( "AddMedChem!");

      // this will cause issues so it's banned
      if
      (
         ( m_TargetMoleculeLinkElementType.empty() && m_EnableDummyAtom ) ||
         ( !m_TargetMoleculeLinkElementType.empty() && !m_EnableDummyAtom )
      )
      {
        BCL_MessageStd
        (
          "\n"
          "Invalid combination of target molecule atom selection options. "
          "To obtain pseudo-reaction-style control over the reaction without "
          "specifying specific atom indices, the following options are available: \n"
          "1. Set 'mutable_elements=X', do not specify any other atom selectors. \n"
          "2. Choose an element type that is unique in the target molecule, such as "
          "Rb, set 'mutable_elements=Rb', 'target_molecule_link_element=Rb', and "
          "'enable_target_dummy_atom=true'. \n"
          "For either option 1 or 2, you are free to change the link/dummy element type "
          "of the medchem fragments via the 'medchem_fragment_link_element' flag depending "
          "on how your library is constructed."
          "\n"
        );
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
      }
      else if
      (
          !m_TargetMoleculeLinkElementType.empty() &&
          m_EnableDummyAtom &&
          m_MutableElements.Find( GetElementTypes().ElementTypeLookup( m_TargetMoleculeLinkElementType)) >= m_MutableElements.GetSize()
      )
      {
        BCL_MessageStd
        (
          "\n"
          "A custom target molecule link element type was specified and enabled, but "
          "the input molecule does not contain any elements of the desired type. "
          "Alternatively, the specified element type is mismatched with the allowed mutable "
          "element types."
          "\n"
        );
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
      }

      // get connecting element types
      const ElementType medchem_fragment_link_element
      (
        m_MedChemFragmentLinkElementType.empty() ? GetElementTypes().e_Undefined : GetElementTypes().ElementTypeLookup( m_MedChemFragmentLinkElementType)
      );
      const ElementType target_molecule_link_element
      (
        m_TargetMoleculeLinkElementType.empty() ? GetElementTypes().e_Undefined : GetElementTypes().ElementTypeLookup( m_TargetMoleculeLinkElementType)
      );

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
          if( medchem_atom_v( i).GetElementType() == medchem_fragment_link_element)
          {
            undefined_index = i;
            break;
          }
        }
        if( undefined_index == util::GetUndefinedSize_t())
        {
          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
        }

        // dummy atoms are allowed only one bonded partner, which must be a heavy atom
        BCL_Assert
        (
          medchem_atom_v( undefined_index).GetBonds().GetSize() == size_t( 1),
          "Encountered a medchem fragment from the library whose target link atom "
          "contains more than one bond! This is not allowed. Atoms designating "
          "pseudo-reactions must be bonded to only one heavy atom. Exiting..."
        );

        BCL_Assert
        (
          medchem_atom_v( undefined_index).GetBonds().Begin()->GetTargetAtom().GetElementType() != GetElementTypes().e_Hydrogen,
          "Encountered a medchem fragment from the library whose target link atom "
          "is bonded to a hdyrogen atom. This is not allowed. Atoms designating "
          "pseudo-reactions must be bonded to only one heavy atom. Exiting..."
        );

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
          if
          (
              // remove undefined elements anyway, even if we use a separate link element type
              medchem_atom_v( i).GetElementType() != GetElementTypes().e_Undefined &&
              i != undefined_index
          )
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
        if
        (
            picked_atom->GetElementType() == GetElementTypes().e_Undefined ||
            ( picked_atom->GetElementType() == target_molecule_link_element && m_EnableDummyAtom )
        )
        {
          // this is the dummy atom
          undefined_base_index = FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom);

          // dummy atoms are allowed only one bonded partner, which must be a heavy atom
          BCL_Assert
          (
            FRAGMENT.GetAtomVector()( undefined_base_index).GetBonds().GetSize() == size_t( 1),
            "The user atom selection specified a target atom dummy atom for linking "
            "that contains more than one bond! This is not allowed. Atoms designating "
            "pseudo-reactions must be bonded to only one heavy atom. Exiting..."
          );

          BCL_Assert
          (
            FRAGMENT.GetAtomVector()( undefined_base_index).GetBonds().Begin()->GetTargetAtom().GetElementType() != GetElementTypes().e_Hydrogen,
            "The user atom selection specified a target atom dummy atom for linking "
            "that is bonded to a hydrogen atom! This is not allowed. Atoms designating "
            "pseudo-reactions must be bonded to only one heavy atom. Exiting..."
          );

          // reassign picked atom to the atom bonded to the dummy atom
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

          if
          (
              picked_atom->GetElementType() == GetElementTypes().e_Undefined ||
              ( picked_atom->GetElementType() == target_molecule_link_element && m_EnableDummyAtom )
          )
          {
            // this is the dummy atom
            undefined_base_index = FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom);

            // dummy atoms are allowed only one bonded partner, which must be a heavy atom
            BCL_Assert
            (
              FRAGMENT.GetAtomVector()( undefined_base_index).GetBonds().GetSize() == size_t( 1),
              "The user atom selection specified a target atom dummy atom for linking "
              "that contains more than one bond! This is not allowed. Atoms designating "
              "pseudo-reactions must be bonded to only one heavy atom. Exiting..."
            );

            BCL_Assert
            (
              FRAGMENT.GetAtomVector()( undefined_base_index).GetBonds().Begin()->GetTargetAtom().GetElementType() != GetElementTypes().e_Hydrogen,
              "The user atom selection specified a target atom dummy atom for linking "
              "that is bonded to a hydrogen atom! This is not allowed. Atoms designating "
              "pseudo-reactions must be bonded to only one heavy atom. Exiting..."
            );

            // reassign picked atom to the atom bonded to the dummy atom
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
        standardizer_2.SetConjugationOfBondTypes( all_defined_atom_v);

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
    void FragmentAddMedChem::SetFragmentLibraryFromFilename( const std::string &FRAGMENTS_FILENAME)
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

    io::Serializer FragmentAddMedChem::GetSerializer() const
    {
      io::Serializer parameters( FragmentMutateInterface::GetSerializer());
      parameters.SetClassDescription
      (
        "Appends a classic medicinal chemistry functional group to the current molecule. "
        "By default, fragments passed with the 'medchem_library' flag will be appended "
        "to the input "
        "molecule; "
        "To obtain pseudo-reaction-style control over the reaction without "
          "specifying specific atom indices, the following options are available: "
          "1. Set 'mutable_elements=X', do not specify any other atom selectors. "
          "2. Choose an element type that is unique in the target molecule, such as "
          "Rb, set 'mutable_elements=Rb', 'target_molecule_link_element=Rb', and "
          "'enable_target_dummy_atom=true'. "
          "For either option 1 or 2, you are free to change the link/dummy element type "
          "of the medchem fragments via the 'medchem_fragment_link_element' flag depending "
          "on how your library is constructed."
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

      parameters.AddInitializer
      (
        "medchem_fragment_link_element",
        "alternative link element type for the medchem fragments; if unspecified, defaults to "
        "the undefined element type (X). "
        "Dummy/linker atoms are required to have only one bond to the target atom of interest. ",
        io::Serialization::GetAgent( &m_MedChemFragmentLinkElementType),
        "X"
      );

      parameters.AddInitializer
      (
        "target_molecule_link_element",
        "alternative link element type for the input molecules; "
        "if you are not using an undefined element (specific 'X' in SDF) to mark the attachment site "
        "by specifying 'mutable_elements=X', then use this flag to change the element type; "
        "requires that 'enable_target_dummy_atoms' is set; "
        "be careful that this is applied appropriately with the mutable_elements atom selector "
        "Dummy/linker atoms are required to have only one bond to the target atom of interest. ",
        io::Serialization::GetAgent( &m_TargetMoleculeLinkElementType),
        ""
      );

      parameters.AddInitializer
      (
        "enable_target_dummy_atom",
        "allows users to specify dummy element types other than undefined (X) for directed "
        "pseudo-reaction-style attachment of medchem fragments; "
        "by default if 'mutable_elements' is set to X and no other atom selectors are specified "
        "then only X elements will 'react' with the link element type in the medchem library "
        "fragments (default is also X). ",
        io::Serialization::GetAgent( &m_EnableDummyAtom),
        "0"
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentAddMedChem::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
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
    std::istream &FragmentAddMedChem::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_FragmentPool, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @param INDENT number of indentations
    //! @return ostream which was written to
    std::ostream &FragmentAddMedChem::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_FragmentPool, OSTREAM, INDENT) << '\n';
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
