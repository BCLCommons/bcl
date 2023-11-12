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
#include "chemistry/bcl_chemistry_fragment_mutate_interface.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_split_largest_component.h"
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "find/bcl_find_collector_interface.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "random/bcl_random_uniform_distribution.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////
  // data //
  //////////

    // static initialization
    storage::Map< std::string, util::ShPtr< FragmentEnsemble> > FragmentMutateInterface::s_RingLibraries = storage::Map< std::string, util::ShPtr< FragmentEnsemble> >();
    storage::Map< std::string, util::ShPtr< FragmentEnsemble> > FragmentMutateInterface::s_MedChemFragmentLibraries = storage::Map< std::string, util::ShPtr< FragmentEnsemble> >();
    sched::Mutex FragmentMutateInterface::s_Mutex = sched::Mutex();

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    // please no constructor

  /////////////////
  // data access //
  /////////////////

//    //! @brief returns all of the mutable atoms from all atom selection methods
//    //! @return the mutable atoms
//    const storage::Vector< size_t> &FragmentMutateInterface::GetAllMutableAtomIndices() const
//    {
//      storage::Set< size_t> unique_atom_indices( m_MutableAtomIndices.Begin(), m_MutableAtomIndices.End());
//      unique_atom_indices.InsertElements( m_MutableElementsAtomindices.Begin(), m_MutableElementsAtomindices.End());
//      unique_atom_indices.InsertElements( m_MutableFragmentsAtomindices.Begin(), m_MutableFragmentsAtomindices.End());
//      unique_atom_indices.EraseKeys( m_FixedAtomindices.Begin(), m_FixedAtomindices.End());
//      unique_atom_indices.EraseKeys( m_FixedElementsAtomindices.Begin(), m_FixedElementsAtomindices.End());
//      unique_atom_indices.EraseKeys( m_FixedFragmentsAtomIndices.Begin(), m_FixedFragmentsAtomIndices.End());
//      const storage::Vector< size_t> index_v( storage::Vector< size_t>( unique_atom_indices.Begin(), unique_atom_indices.End()));
//      return index_v;
//    }

    //! @brief returns individually specified mutable atoms
    //! @return the mutable atoms
    const storage::Vector< size_t> &FragmentMutateInterface::GetMutableAtomIndices() const
    {
      return m_MutableAtomIndices;
    }

    //! @brief returns mutable elements
    //! @return the mutable elements
    const storage::Vector< ElementType> &FragmentMutateInterface::GetMutableElements() const
    {
      return m_MutableElements;
    }

    //! @brief returns atom indices of mutable elements
    //! @return the mutable atoms
    const storage::Vector< size_t> &FragmentMutateInterface::GetMutableElementsAtomIndices() const
    {
      return m_MutableElementsAtomindices;
    }

    //! @brief returns mutable fragments
    //! @return the mutable fragments
    const FragmentEnsemble &FragmentMutateInterface::GetMutableFragments() const
    {
      return m_MutableFragments;
    }

    //! @brief returns atom indices of mutable fragments
    //! @return the mutable atoms
    const storage::Vector< size_t> &FragmentMutateInterface::GetMutableFragmentsAtomIndices() const
    {
      return m_MutableFragmentsAtomindices;
    }

    //! @brief returns passed mutable atoms
    //! @return the mutable atoms
    const storage::Vector< size_t> &FragmentMutateInterface::GetFixedAtomIndices() const
    {
      return m_FixedAtomindices;
    }

    //! @brief returns mutable elements
    //! @return the mutable elements
    const storage::Vector< ElementType> &FragmentMutateInterface::GetFixedElements() const
    {
      return m_FixedElements;
    }

    //! @brief returns atom indices of mutable elements
    //! @return the mutable atoms
    const storage::Vector< size_t> &FragmentMutateInterface::GetFixedElementsAtomIndices() const
    {
      return m_FixedElementsAtomindices;
    }

    //! @brief returns mutable fragments
    //! @return the mutable fragments
    const FragmentEnsemble &FragmentMutateInterface::GetFixedFragments() const
    {
      return m_FixedFragments;
    }

    //! @brief returns atom indices of mutable fragments
    //! @return the mutable atoms
    const storage::Vector< size_t> &FragmentMutateInterface::GetFixedFragmentsAtomIndices() const
    {
      return m_FixedFragmentsAtomIndices;
    }

    //! @brief returns individually specified mutable atoms
    //! @return the mutable atoms
    const storage::Vector< size_t> &FragmentMutateInterface::GetRingLibraryMutableAtomIndices() const
    {
      return m_RingLibraryMutableAtomindices;
    }

    //! @brief returns mutable elements
    //! @return the mutable elements
    const storage::Vector< ElementType> &FragmentMutateInterface::GetRingLibraryMutableElements() const
    {
      return m_RingLibraryMutableElements;
    }

    //! @brief returns atom indices of mutable elements
    //! @return the mutable atoms
    const storage::Vector< size_t> &FragmentMutateInterface::GetRingLibraryMutableElementsAtomIndices() const
    {
      return m_RingLibraryMutableElementsAtomindices;
    }

    //! @brief returns mutable fragments
    //! @return the mutable fragments
    const FragmentEnsemble &FragmentMutateInterface::GetRingLibraryMutableFragments() const
    {
      return m_RingLibraryMutableFragments;
    }

    //! @brief returns atom indices of mutable fragments
    //! @return the mutable atoms
    const storage::Vector< size_t> &FragmentMutateInterface::GetRingLibraryMutableFragmentsAtomIndices() const
    {
      return m_RingLibraryMutableFragmentsAtomindices;
    }

    //! @brief returns passed mutable atoms
    //! @return the mutable atoms
    const storage::Vector< size_t> &FragmentMutateInterface::GetRingLibraryFixedAtomIndices() const
    {
      return m_RingLibraryFixedAtomindices;
    }

    //! @brief returns mutable elements
    //! @return the mutable elements
    const storage::Vector< ElementType> &FragmentMutateInterface::GetRingLibraryFixedElements() const
    {
      return m_RingLibraryFixedElements;
    }

    //! @brief returns atom indices of mutable elements
    //! @return the mutable atoms
    const storage::Vector< size_t> &FragmentMutateInterface::GetRingLibraryFixedElementsAtomIndices() const
    {
      return m_RingLibraryFixedElementsAtomindices;
    }

    //! @brief returns mutable fragments
    //! @return the mutable fragments
    const FragmentEnsemble &FragmentMutateInterface::GetRingLibraryFixedFragments() const
    {
      return m_RingLibraryFixedFragments;
    }

    //! @brief returns atom indices of mutable fragments
    //! @return the mutable atoms
    const storage::Vector< size_t> &FragmentMutateInterface::GetRingLibraryFixedFragmentsAtomIndices() const
    {
      return m_RingLibraryFixedFragmentsAtomIndices;
    }

  ///////////////
  // operators //
  ///////////////

  ////////////////
  // operations //
  ////////////////

    //! @brief set druglikeness type
    void FragmentMutateInterface::SetDruglikenessType( const std::string &DRUGLIKENESS_TYPE)
    {
      m_DrugLikenessType = DRUGLIKENESS_TYPE;
    }

    //! @brief set MDL SDF property label for the receptor path for BCL structure-based scoring
    void FragmentMutateInterface::SetMDLReceptorLabel( const std::string &MDL)
    {
      m_MDL = MDL;
    }

    //! @brief set the scoring filter type
    void FragmentMutateInterface::SetPropertyScorer( const descriptor::CheminfoProperty &PROPERTY_SCORER)
    {
      m_PropertyScorer = PROPERTY_SCORER;
    }

    //! @brief set maximum number of mutate attempts
    void FragmentMutateInterface::SetNumberMaxMutates( const size_t N_MAX_MUTATES)
    {
      m_NumberMaxAttempts = N_MAX_MUTATES;
    }

    //! @brief set bool to shuffle sequence of hydrogen atom removal when opening valences
    void FragmentMutateInterface::SetOVShuffleH( const bool SHUFFLE_H)
    {
      m_OVShuffleH = SHUFFLE_H;
    }

    //! @brief set bool to reverse the order of hydrogen atom removal when opening valences
    void FragmentMutateInterface::SetOVReverseH( const bool REVERSE_H)
    {
      m_OVReverse = REVERSE_H;
    }

    //! @brief set generic mutable atoms
    void FragmentMutateInterface::SetMutableAtomIndices( const storage::Vector< size_t> &MUTABLE_ATOM_INDICES)
    {
      m_MutableAtomIndices = MUTABLE_ATOM_INDICES;
    }

    //! @brief set the mutable atoms based on mutable elements
    void FragmentMutateInterface::SetMutableElementsAtomIndices
    (
      const FragmentComplete &FRAGMENT,
      const storage::Vector< ElementType> &MUTABLE_ELEMENTS
    )
    {
      m_MutableElementsAtomindices = FragmentTrackMutableAtoms::SetMutableElements( FRAGMENT, MUTABLE_ELEMENTS);
    }

    //! @brief set the mutable atoms based on mutable fragments
    void FragmentMutateInterface::SetMutableFragmentsAtomIndices
    (
      const FragmentComplete &FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const bool COMPLEMENT
    )
    {
      // add atom indices of allowed mutable fragments
      m_MutableFragmentsAtomindices = FragmentTrackMutableAtoms::SetMutableFragments( FRAGMENT, MUTABLE_FRAGMENTS, COMPLEMENT);
    }

    //! @brief set generic mutable atoms
    void FragmentMutateInterface::SetFixedAtomIndices( const storage::Vector< size_t> &FIXED_ATOM_INDICES)
    {
      m_FixedAtomindices = FIXED_ATOM_INDICES;
    }

    //! @brief set the mutable atoms based on mutable elements
    void FragmentMutateInterface::SetFixedElementsAtomIndices
    (
      const FragmentComplete &FRAGMENT,
      const storage::Vector< ElementType> &FIXED_ELEMENTS
    )
    {
      m_FixedElementsAtomindices = FragmentTrackMutableAtoms::SetMutableElements( FRAGMENT, FIXED_ELEMENTS);
    }

    //! @brief set the mutable atoms based on mutable fragments
    void FragmentMutateInterface::SetFixedFragmentsAtomIndices
    (
      const FragmentComplete &FRAGMENT,
      const FragmentEnsemble &FIXED_FRAGMENTS,
      const bool COMPLEMENT
    )
    {
      // add atom indices of allowed mutable fragments
      m_FixedFragmentsAtomIndices = FragmentTrackMutableAtoms::SetMutableFragments( FRAGMENT, FIXED_FRAGMENTS, COMPLEMENT);
    }

    //! @brief set generic mutable atoms
    void FragmentMutateInterface::SetRingLibraryMutableAtomIndices( const storage::Vector< size_t> &MUTABLE_ATOM_INDICES)
    {
      m_RingLibraryMutableAtomindices = MUTABLE_ATOM_INDICES;
    }

    //! @brief set the mutable atoms based on mutable elements
    void FragmentMutateInterface::SetRingLibraryMutableElementsAtomIndices
    (
      const FragmentComplete &FRAGMENT,
      const storage::Vector< ElementType> &MUTABLE_ELEMENTS
    )
    {
      m_RingLibraryMutableElementsAtomindices = FragmentTrackMutableAtoms::SetMutableElements( FRAGMENT, MUTABLE_ELEMENTS);
    }

    //! @brief set the mutable atoms based on mutable fragments
    void FragmentMutateInterface::SetRingLibraryMutableFragmentsAtomIndices
    (
      const FragmentComplete &FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const bool COMPLEMENT
    )
    {
      // add atom indices of allowed mutable fragments
      m_RingLibraryMutableFragmentsAtomindices = FragmentTrackMutableAtoms::SetMutableFragments( FRAGMENT, MUTABLE_FRAGMENTS, COMPLEMENT);
    }

    //! @brief set generic mutable atoms
    void FragmentMutateInterface::SetRingLibraryFixedAtomIndices( const storage::Vector< size_t> &FIXED_ATOM_INDICES)
    {
      m_RingLibraryFixedAtomindices = FIXED_ATOM_INDICES;
    }

    //! @brief set the mutable atoms based on mutable elements
    void FragmentMutateInterface::SetRingLibraryFixedElementsAtomIndices
    (
      const FragmentComplete &FRAGMENT,
      const storage::Vector< ElementType> &FIXED_ELEMENTS
    )
    {
      m_RingLibraryFixedElementsAtomindices = FragmentTrackMutableAtoms::SetMutableElements( FRAGMENT, FIXED_ELEMENTS);
    }

    //! @brief set the mutable atoms based on mutable fragments
    void FragmentMutateInterface::SetRingLibraryFixedFragmentsAtomIndices
    (
      const FragmentComplete &FRAGMENT,
      const FragmentEnsemble &FIXED_FRAGMENTS,
      const bool COMPLEMENT
    )
    {
      // add atom indices of allowed mutable fragments
      m_RingLibraryFixedFragmentsAtomIndices = FragmentTrackMutableAtoms::SetMutableFragments( FRAGMENT, FIXED_FRAGMENTS, COMPLEMENT);
    }

    //! @brief set the scaffold fragment
    void FragmentMutateInterface::SetScaffoldFragment( const FragmentComplete &SCAFFOLD_FRAGMENT)
    {
      m_ScaffoldFragment = SCAFFOLD_FRAGMENT;
    }

    //! @brief set the properties from one molecule onto the clean molecule;
    //! useful when properties computed during mutate are desirable to keep
    void FragmentMutateInterface::SetPropertiesFromOther
    (
      FragmentComplete &TARGET_FRAGMENT,
      const FragmentComplete &REFERENCE_FRAGMENT
    ) const
    {
      // loop over properties of reference fragment
      for
      (
          auto prop_itr( REFERENCE_FRAGMENT.GetStoredProperties().Begin()),
          prop_itr_end( REFERENCE_FRAGMENT.GetStoredProperties().End());
          prop_itr != prop_itr_end;
          ++prop_itr
      )
      {
        TARGET_FRAGMENT.StoreProperty( prop_itr->first, prop_itr->second);
      }
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief remove a hydrogen atom from a target atom
    //! @param FRAGMENT the molecule of interest
    //! @param ATOM_INDEX the index of the atom in the molecule of interest
    //! @param SHUFFLE_H if true, randomly select a hydrogen atom to remove
    //! @param REVERSE_H if true and not SHUFFLE_H, begin removal with the highest index hydrogen atom
    //! @return the new molecule, the index of the desired atom, and the original index of the removed hydrogen atom
    storage::Triplet< FragmentComplete, size_t, size_t> FragmentMutateInterface::OpenValence
    (
      const FragmentComplete &FRAGMENT,
      const size_t &ATOM_INDEX,
      const bool SHUFFLE_H,
      const bool REVERSE_H
    ) const
    {
      // find a hydrogen atom attached to specified atom index
      storage::Vector< size_t> h_indices;

      // grab the largest atom index hydrogen atoms first
      if( REVERSE_H)
      {
        for
        (
            auto bond_itr( FRAGMENT.GetAtomVector()( ATOM_INDEX).GetBonds().ReverseBegin()),
            bond_itr_end( FRAGMENT.GetAtomVector()( ATOM_INDEX).GetBonds().ReverseEnd());
            bond_itr != bond_itr_end;
            ++bond_itr
        )
        {
          if( bond_itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Hydrogen)
          {
            // add hydrogen atom
            h_indices.PushBack( FRAGMENT.GetAtomVector().GetAtomIndex( bond_itr->GetTargetAtom()));

            // break after first hydrogen atom
            if( !SHUFFLE_H)
            {
              break;
            }
          }
        }
      }
      else
      {
        for
        (
            auto bond_itr( FRAGMENT.GetAtomVector()( ATOM_INDEX).GetBonds().Begin()),
            bond_itr_end( FRAGMENT.GetAtomVector()( ATOM_INDEX).GetBonds().End());
            bond_itr != bond_itr_end;
            ++bond_itr
        )
        {
          if( bond_itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Hydrogen)
          {
            h_indices.PushBack( FRAGMENT.GetAtomVector().GetAtomIndex( bond_itr->GetTargetAtom()));
            if( !SHUFFLE_H)
            {
              break;
            }
          }
        }
      }

      // if no hydrogen atoms then return input
      if( !h_indices.GetSize())
      {
        return storage::Triplet< FragmentComplete, size_t, size_t>( FRAGMENT, ATOM_INDEX, util::GetUndefinedSize_t());
      }

      // randomly select one of the hydrogen atoms to remove
      if( SHUFFLE_H)
      {
        h_indices.Shuffle();
      }
      size_t h_index( h_indices( 0));

      // adjust the reference atom index as needed
      size_t new_atom_index( ATOM_INDEX);
      if( h_index < ATOM_INDEX)
      {
        new_atom_index -= size_t( 1);
      }

      // generate new atom indices excluding the hydrogen atom
      storage::Vector< size_t> keep_indices;
      for
      (
          size_t i( 0); i < FRAGMENT.GetSize(); ++i
      )
      {
        if( i == h_index)
        {
          continue;
        }
        keep_indices.PushBack( i);
      }

      // make new molecule without hydrogen atom
      AtomVector< AtomComplete> new_atom_v( FRAGMENT.GetAtomVector());
      new_atom_v.Reorder( keep_indices);
      FragmentComplete new_frag( new_atom_v, "");

      // return new molecule with the updated atom index
      return storage::Triplet< FragmentComplete, size_t, size_t>( new_frag, new_atom_index, h_index);
    }

    //! @brief checks whether substitution at this atom is ortho, meta, or para directed
    //! @param MOLECULE the small molecule of interest
    //! @param ATOM simple pointer to the atom of interest in the molecule
    //! @return return true if the substitution is directed correctly
    bool FragmentMutateInterface::IsRingSubstitutionDirected
    (
      const FragmentComplete &MOLECULE,
      util::SiPtr< const AtomConformationalInterface> &ATOM
    ) const
    {
      // compute our chosen atom pi charge
      double atom_pi_charge
      (
        descriptor::GetCheminfoProperties().calc_PiCharge->CollectValuesOnEachElementOfObject
        (
          MOLECULE
        )( MOLECULE.GetAtomVector().GetAtomIndex( *ATOM))
      );

      // check if the pi-charge is negative on the picked atom
      // ortho/para-directing
      if( atom_pi_charge < 0.0)
      {
        return true;
      }
      // if not negative, check whether less positive than adjacent atoms
      // meta-directing
      else
      {
        size_t n_bonded_atoms( ATOM->GetBonds().GetSize());
        size_t n_h_cov_bonds( ATOM->GetNumberCovalentlyBoundHydrogens());
        size_t n_less_pos( 0);
        for
        (
            auto bond_itr( ATOM->GetBonds().Begin()), bond_itr_end( ATOM->GetBonds().End());
            bond_itr != bond_itr_end;
            ++bond_itr
        )
        {
          // compute bonded atom pi charge
          double bonded_atom_pi_charge
          (
            descriptor::GetCheminfoProperties().calc_PiCharge->CollectValuesOnEachElementOfObject( MOLECULE)
            (
                MOLECULE.GetAtomVector().GetAtomIndex( bond_itr->GetTargetAtom())
            )
          );
          if( atom_pi_charge < bonded_atom_pi_charge)
          {
            n_less_pos += size_t( 1);
          }
        }
        if( n_less_pos == n_bonded_atoms - n_h_cov_bonds)
        {
          return true;
        }
      }
      // not directed correctly
      return false;
    }

    //! @brief select an atom from the target fragment
    //! @brief MOLECULE molecule from which to choose atom
    //! @param GLOBAL randomly select an atom from all atom indices
    //! @param RING_LIBRARY randomly select an atom from a ring (previously chosen
    //! from the user-passed ring library)
    //! @return return the chosen atom; base class default chooses global random atom
    //! using the input molecule atom selection specifications (not ring library data)
    util::SiPtr< const AtomConformationalInterface> FragmentMutateInterface::PickAtom
    (
      const FragmentComplete &MOLECULE,
      const bool GLOBAL,
      const bool RING_LIBRARY
    ) const
    {
      util::SiPtr< const AtomConformationalInterface> picked_atom;
      if( RING_LIBRARY)
      {
        picked_atom = util::SiPtr< const AtomConformationalInterface>
        (
          FragmentTrackMutableAtoms::GetAtomFromMutable
          (
            MOLECULE,
            GLOBAL,
            m_RingLibraryMutableAtomindices,
            m_RingLibraryMutableElements,
            m_RingLibraryMutableFragments,
            m_RingLibraryFixedAtomindices,
            m_RingLibraryFixedElements,
            m_RingLibraryFixedFragments,
            m_RingLibraryComplementMutableFragments,
            m_RingLibraryComplementFixedFragments,
            ConformationGraphConverter
            (
              m_RingLibraryMutableAtomComparisonType,
              m_RingLibraryMutableBondComparisonType
            ),
            ConformationGraphConverter
            (
              m_RingLibraryFixedAtomComparisonType,
              m_RingLibraryFixedBondComparisonType
            )
          )
        );
      }
      else
      {
        picked_atom = util::SiPtr< const AtomConformationalInterface>
        (
          FragmentTrackMutableAtoms::GetAtomFromMutable
          (
            MOLECULE,
            GLOBAL,
            m_MutableAtomIndices,
            m_MutableElements,
            m_MutableFragments,
            m_FixedAtomindices,
            m_FixedElements,
            m_FixedFragments,
            m_ComplementMutableFragments,
            m_ComplementFixedFragments,
            ConformationGraphConverter
            (
              m_MutableAtomComparisonType,
              m_MutableBondComparisonType
            ),
            ConformationGraphConverter
            (
              m_FixedAtomComparisonType,
              m_FixedBondComparisonType
            )
          )
        );
      }
      return picked_atom;
    }

    //! @brief selects a connection atom from the ring chosen from the ring library
    //! @param FRAGMENT molecule from which atom will be chosen
    //! @param RING_LIBRARY randomly select an atom from a ring (previously chosen
    //! from the user-passed ring library)
    //! @param H_BONDED_HEAVY if a hydrogen atom is selected, use its bonded heavy atom
    //! as the chosen atom instead
    //! @return the index of the chosen atom in the parent ring
    size_t FragmentMutateInterface::RunPickAtom
    (
      const FragmentComplete &FRAGMENT,
      const bool RING_LIBRARY,
      const bool H_BONDED_HEAVY
    ) const
    {
      // initialize a pointer to our atom
      util::SiPtr< const AtomConformationalInterface> picked_atom;

      // the atom is being chosen from a ring in the ring library
      if( ( m_RingLibraryMutableAtomindices.GetSize() || m_RingLibraryMutableElements.GetSize() || m_RingLibraryMutableFragments.GetSize()) && RING_LIBRARY)
      {
        picked_atom = this->PickAtom( FRAGMENT, false, true);
      }

      // the atom is being chosen from a ring in the ring library
      else if( RING_LIBRARY)
      {
        picked_atom = this->PickAtom( FRAGMENT, true, true);
      }

      // the atom is being chose from the input fragment
      else if( m_MutableAtomIndices.GetSize() || m_MutableElements.GetSize() || m_MutableFragments.GetSize())
      {
        picked_atom = this->PickAtom( FRAGMENT, false, false);
      }

      // the atom is being chose from the input fragment
      else
      {
        picked_atom = this->PickAtom( FRAGMENT, true, false);
      }

      // if undefined then return undefined index
      if( !picked_atom.IsDefined())
      {
        BCL_MessageStd( "Unable to identify a valid mutable atom");
        return util::GetUndefinedSize_t();
      }

      // if a hydrogen atom is selected then just grab its bonded atom
      if( picked_atom->GetElementType() == GetElementTypes().e_Hydrogen && H_BONDED_HEAVY)
      {
        // must have bonded atom
        if( !picked_atom->GetBonds().GetSize())
        {
          BCL_MessageStd( "There are no bonds from this hydrogen atom");
          return util::GetUndefinedSize_t();
        }
        // Note that this does not account for certain edge cases, e.g., H2
        picked_atom = util::SiPtr< const AtomConformationalInterface>( picked_atom->GetBonds().Begin()->GetTargetAtom());
      }
      return FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom);
    }

    io::Serializer FragmentMutateInterface::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Base class for alchemical fragment mutations."
      );

      parameters.AddInitializer
      (
        "druglikeness_type",
        "the type of druglikeness filter to apply; "
        "returns input molecule if fails filter",
        io::Serialization::GetAgent( &m_DrugLikenessType),
        "None"
      );

      parameters.AddInitializer
      (
        "scaffold_mol",
        "scaffold molecule for substructure-based alignment during 3D conformer construction",
        io::Serialization::GetAgent( &m_ScaffoldFragmentFilename),
        ""
      );

      parameters.AddInitializer
      (
        "n_max_attempts",
        "max number of attempts at performing the mutate; "
        "useful to increase this count if there are a large number of mutable atoms",
        io::Serialization::GetAgent( &m_NumberMaxAttempts),
        "10"
      );

      parameters.AddInitializer
      (
        "ov_shuffle_h",
        "If true, randomly select a hydrogen atom for removal when opening valences during mutation. If false, "
        "select the lowest index hydrogen atom bonded to the target heavy atom. Several mutates, such as "
        "FragmentAddMedChem and FragmentExtendWithLinker, explicitly remove hydrogen atoms to open a valence "
        "for subsequent covalent modification. If an atom contains more than one hydrogen atom, it may be "
        "useful to select the order in which these valences will be opened. For example, if FragmentAddMedchem "
        "is being used to modify a potential stereocenter with two hydrogen atoms, then setting 'ov_shuffle_h' to "
        "false would mean that one could control the stereochemistry and selectively modify only the first position "
        "(or the second; see 'ov_reverse') or apply two mutates in succession with a determined output. In contrast, "
        "you may just be trying to grow a molecule from a methyl head group, in which case you do not want to be "
        "constantly growing from the same orientation as the lowest index hydrogen atom. In this case, setting "
        "'ov_shuffle_h' to true will allow more complete geometric sampling.",
        io::Serialization::GetAgent( &m_OVShuffleH),
        "true"
      );

      parameters.AddInitializer
      (
        "ov_reverse",
        "If false and 'ov_shuffle_h' set to false, mutates that require explicit valence opening will proceed "
        "from the lowest index hydrogen atom. If true and 'ov_shuffle_h' false, proceed from the highest index "
        "hydrogen atom. No effect if 'ov_shuffle_h' is true. See 'ov_shuffle_h' for more details.",
        io::Serialization::GetAgent( &m_OVReverse),
        "false"
      );

      parameters.AddInitializer
      (
        "corina",
        "optionally generate a 3D conformer with Corina instead of BCL::Conf"
        "only a valid option for users with the Corina program installed; "
        "useful if the 3D conformer does not need to retain pose information, needs to be "
        "generated very rapidly, and/or the dataset is going to be scored with a QSAR model "
        "trained on Corina conformers.",
        io::Serialization::GetAgent( &m_Corina),
        "false"
      );

      parameters.AddInitializer
      (
        "ring_library_mutable_atoms",
        "atom indices (0-indexed) that can be mutated; "
        "no effect in mutates that do not utilize an external ring library for design",
        io::Serialization::GetAgent( &m_RingLibraryMutableAtoms),
        ""
      );

      parameters.AddInitializer
      (
        "ring_library_mutable_elements",
        "element types that can be mutated; "
        "no effect in mutates that do not utilize an external ring library for design",
        io::Serialization::GetAgent( &m_RingLibraryMutableElementsString),
        ""
      );

      parameters.AddInitializer
      (
        "ring_library_mutable_fragments",
        "allow atoms that match these fragments in substructure comparisons to be "
        "mutated; can be inverted to set mutable atoms to the complement subgraph "
        "isomorphism atoms via 'complement_mutable_fragments'; "
        "no effect in mutates that do not utilize an external ring library for design",
        io::Serialization::GetAgent( &m_RingLibraryMutableFragmentsFilename),
        ""
      );

      parameters.AddInitializer
      (
        "ring_library_complement_mutable_fragments",
        "invert the subgraph isomorphisms between the molecule of interest and "
        "the mutable fragments such that the derived mutable atoms are the "
        "non-common atoms; "
        "no effect in mutates that do not utilize an external ring library for design",
        io::Serialization::GetAgent( &m_RingLibraryComplementMutableFragments),
        "false"
      );

      parameters.AddInitializer
      (
        "ring_library_mutable_atom_comparison",
        "atom data that are compared to determine whether atoms are equivalent "
        "during mutable fragment substructure comparisons; "
        "no effect in mutates that do not utilize an external ring library for design",
        io::Serialization::GetAgent( &m_RingLibraryMutableAtomComparisonType),
        "ElementType"
      );

      parameters.AddInitializer
      (
        "ring_library_mutable_bond_comparison",
        "bond data that are compared to determine whether bonds are equivalent "
        "during mutable fragment substructure comparisons; "
        "no effect in mutates that do not utilize an external ring library for design",
        io::Serialization::GetAgent( &m_RingLibraryMutableBondComparisonType),
        "BondOrderAmideOrAromaticWithRingness"
      );

      parameters.AddInitializer
      (
        "ring_library_fixed_atoms",
        "atom indices (0-indexed) that will be fixed after mutable atoms are assigned; "
        "no effect in mutates that do not utilize an external ring library for design",
        io::Serialization::GetAgent( &m_RingLibraryFixedAtoms),
        ""
      );

      parameters.AddInitializer
      (
        "ring_library_fixed_elements",
        "element types that will be fixed after mutable atoms are assigned; "
        "no effect in mutates that do not utilize an external ring library for design",
        io::Serialization::GetAgent( &m_RingLibraryFixedElementsString),
        ""
      );

      parameters.AddInitializer
      (
        "ring_library_fixed_fragments",
        "fix atoms that match these fragments in substructure comparisons; "
        "can be inverted to set fixed atoms to the complement subgraph "
        "isomorphism atoms via 'complement_fixed_fragments'; "
        "no effect in mutates that do not utilize an external ring library for design",
        io::Serialization::GetAgent( &m_RingLibraryFixedFragmentsFilename),
        ""
      );

      parameters.AddInitializer
      (
        "ring_library_complement_fixed_fragments",
        "invert the subgraph isomorphisms between the molecule of interest and "
        "the fixed fragments such that the derived mutable atoms are the "
        "non-common atoms; "
        "no effect in mutates that do not utilize an external ring library for design",
        io::Serialization::GetAgent( &m_RingLibraryComplementFixedFragments),
        "false"
      );

      parameters.AddInitializer
      (
        "ring_library_fixed_atom_comparison",
        "atom data that are compared to determine whether atoms are equivalent "
        "during fixed fragment substructure comparisons; "
        "no effect in mutates that do not utilize an external ring library for design",
        io::Serialization::GetAgent( &m_RingLibraryFixedAtomComparisonType),
        "ElementType"
      );

      parameters.AddInitializer
      (
        "ring_library_fixed_bond_comparison",
        "bond data that are compared to determine whether bonds are equivalent "
        "during fixed fragment substructure comparisons; "
        "no effect in mutates that do not utilize an external ring library for design",
        io::Serialization::GetAgent( &m_RingLibraryFixedBondComparisonType),
        "BondOrderAmideOrAromaticWithRingness"
      );

      parameters.AddInitializer
      (
        "mutable_atoms",
        "atom indices (0-indexed) that can be mutated",
        io::Serialization::GetAgent( &m_MutableAtoms),
        ""
      );

      parameters.AddInitializer
      (
        "mutable_elements",
        "element types that can be mutated",
        io::Serialization::GetAgent( &m_MutableElementsString),
        ""
      );

      parameters.AddInitializer
      (
        "mutable_fragments",
        "allow atoms that match these fragments in substructure comparisons to be "
        "mutated; can be inverted to set mutable atoms to the complement subgraph "
        "isomorphism atoms via 'complement_mutable_fragments'",
        io::Serialization::GetAgent( &m_MutableFragmentsFilename),
        ""
      );

      parameters.AddInitializer
      (
        "complement_mutable_fragments",
        "invert the subgraph isomorphisms between the molecule of interest and "
        "the mutable fragments such that the derived mutable atoms are the "
        "non-common atoms",
        io::Serialization::GetAgent( &m_ComplementMutableFragments),
        "false"
      );

      parameters.AddInitializer
      (
        "mutable_atom_comparison",
        "atom data that are compared to determine whether atoms are equivalent "
        "during mutable fragment substructure comparisons",
        io::Serialization::GetAgent( &m_MutableAtomComparisonType),
        "ElementType"
      );

      parameters.AddInitializer
      (
        "mutable_bond_comparison",
        "bond data that are compared to determine whether bonds are equivalent "
        "during mutable fragment substructure comparisons",
        io::Serialization::GetAgent( &m_MutableBondComparisonType),
        "BondOrderAmideOrAromaticWithRingness"
      );

      parameters.AddInitializer
      (
        "fixed_atoms",
        "atom indices (0-indexed) that will be fixed after mutable atoms are assigned",
        io::Serialization::GetAgent( &m_FixedAtoms),
        ""
      );

      parameters.AddInitializer
      (
        "fixed_elements",
        "element types that will be fixed after mutable atoms are assigned",
        io::Serialization::GetAgent( &m_FixedElementsString),
        ""
      );

      parameters.AddInitializer
      (
        "fixed_fragments",
        "fix atoms that match these fragments in substructure comparisons; "
        "can be inverted to set fixed atoms to the complement subgraph "
        "isomorphism atoms via 'complement_fixed_fragments'",
        io::Serialization::GetAgent( &m_FixedFragmentsFilename),
        ""
      );

      parameters.AddInitializer
      (
        "complement_fixed_fragments",
        "invert the subgraph isomorphisms between the molecule of interest and "
        "the fixed fragments such that the derived mutable atoms are the "
        "non-common atoms",
        io::Serialization::GetAgent( &m_ComplementFixedFragments),
        "false"
      );

      parameters.AddInitializer
      (
        "fixed_atom_comparison",
        "atom data that are compared to determine whether atoms are equivalent "
        "during fixed fragment substructure comparisons",
        io::Serialization::GetAgent( &m_FixedAtomComparisonType),
        "ElementType"
      );

      parameters.AddInitializer
      (
        "fixed_bond_comparison",
        "bond data that are compared to determine whether bonds are equivalent "
        "during fixed fragment substructure comparisons",
        io::Serialization::GetAgent( &m_FixedBondComparisonType),
        "BondOrderAmideOrAromaticWithRingness"
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentMutateInterface::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      // read in scaffold filename
      s_Mutex.Lock();
      if( m_ScaffoldFragmentFilename.size())
      {
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_ScaffoldFragmentFilename);
        FragmentEnsemble scaffold_mol;
        scaffold_mol.ReadMoreFromMdl( input);
        m_ScaffoldFragment = scaffold_mol.GetMolecules().FirstElement();
        io::File::CloseClearFStream( input);
      }
      s_Mutex.Unlock();

      //////// Input fragment atom selection options ////////

      // read in reference filename
      s_Mutex.Lock();
      if( m_MutableFragmentsFilename.size())
      {
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_MutableFragmentsFilename);
        m_MutableFragments.ReadMoreFromMdl( input);
        io::File::CloseClearFStream( input);
      }
      s_Mutex.Unlock();

      // read in mutable atom indices
      if( m_MutableAtoms.size())
      {
        m_MutableAtomIndices.Reset();
        m_MutableAtomIndices = util::SplitStringToNumerical< size_t>( m_MutableAtoms);
      }

      // read in allowed mutable elements when choosing mutable atom
      if( m_MutableElementsString.size())
      {
        // parse input
        const storage::Vector< std::string> mutable_elements
        (
          util::SplitString( util::TrimString( m_MutableElementsString), " \t\n\r,")
        );

        // stupid check to add only the correct elements
        // TODO: this should be directly serializable from element types to make concise
        m_MutableElements.Reset();
        for( size_t e_i( 0), e_sz( mutable_elements.GetSize()); e_i < e_sz; ++e_i)
        {
          // Hydrogen
          if( mutable_elements( e_i) == "H")
          {
            m_MutableElements.PushBack( GetElementTypes().e_Hydrogen);
          }
          // Boron
          if( mutable_elements( e_i) == "B")
          {
            m_MutableElements.PushBack( GetElementTypes().e_Boron);
          }
          // Carbon
          if( mutable_elements( e_i) == "C")
          {
            m_MutableElements.PushBack( GetElementTypes().e_Carbon);
          }
          // Oxygen
          if( mutable_elements( e_i) == "O")
          {
            m_MutableElements.PushBack( GetElementTypes().e_Oxygen);
          }
          // Nitrogen
          if( mutable_elements( e_i) == "N")
          {
            m_MutableElements.PushBack( GetElementTypes().e_Nitrogen);
          }
          // Potassium
          if( mutable_elements( e_i) == "P")
          {
            m_MutableElements.PushBack( GetElementTypes().e_Phosphorus);
          }
          // Sulfur
          if( mutable_elements( e_i) == "S")
          {
            m_MutableElements.PushBack( GetElementTypes().e_Sulfur);
          }
          // Selenium
          if( mutable_elements( e_i) == "Se")
          {
            m_MutableElements.PushBack( GetElementTypes().e_Selenium);
          }
          // Fluorine
          if( mutable_elements( e_i) == "F")
          {
            m_MutableElements.PushBack( GetElementTypes().e_Fluorine);
          }
          // Chlorine
          if( mutable_elements( e_i) == "Cl")
          {
            m_MutableElements.PushBack( GetElementTypes().e_Chlorine);
          }
          // Bromine
          if( mutable_elements( e_i) == "Br")
          {
            m_MutableElements.PushBack( GetElementTypes().e_Bromine);
          }
          // Iodine
          if( mutable_elements( e_i) == "I")
          {
            m_MutableElements.PushBack( GetElementTypes().e_Iodine);
          }
          // Undefined
          if( mutable_elements( e_i) == "X")
          {
            m_MutableElements.PushBack( GetElementTypes().e_Undefined);
          }
        }
      }

      // read in reference filename
      s_Mutex.Lock();
      if( m_FixedFragmentsFilename.size())
      {
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_FixedFragmentsFilename);
        m_FixedFragments.ReadMoreFromMdl( input);
        io::File::CloseClearFStream( input);
      }
      s_Mutex.Unlock();

      // read in mutable atom indices
      if( m_FixedAtoms.size())
      {
        m_FixedAtomindices.Reset();
        m_FixedAtomindices = util::SplitStringToNumerical< size_t>( m_FixedAtoms);
      }

      // read in allowed mutable elements when choosing mutable atom
      if( m_FixedElementsString.size())
      {
        // parse input
        const storage::Vector< std::string> fixed_elements
        (
          util::SplitString( util::TrimString( m_FixedElementsString), " \t\n\r,")
        );

        // check to add only the correct elements
        m_FixedElements.Reset();
        for( size_t e_i( 0), e_sz( fixed_elements.GetSize()); e_i < e_sz; ++e_i)
        {
          // Hydrogen
          if( fixed_elements( e_i) == "H")
          {
            m_FixedElements.PushBack( GetElementTypes().e_Hydrogen);
          }
          // Boron
          if( fixed_elements( e_i) == "B")
          {
            m_FixedElements.PushBack( GetElementTypes().e_Boron);
          }
          // Carbon
          if( fixed_elements( e_i) == "C")
          {
            m_FixedElements.PushBack( GetElementTypes().e_Carbon);
          }
          // Oxygen
          if( fixed_elements( e_i) == "O")
          {
            m_FixedElements.PushBack( GetElementTypes().e_Oxygen);
          }
          // Nitrogen
          if( fixed_elements( e_i) == "N")
          {
            m_FixedElements.PushBack( GetElementTypes().e_Nitrogen);
          }
          // Potassium
          if( fixed_elements( e_i) == "P")
          {
            m_FixedElements.PushBack( GetElementTypes().e_Phosphorus);
          }
          // Sulfur
          if( fixed_elements( e_i) == "S")
          {
            m_FixedElements.PushBack( GetElementTypes().e_Sulfur);
          }
          // Selenium
          if( fixed_elements( e_i) == "Se")
          {
            m_FixedElements.PushBack( GetElementTypes().e_Selenium);
          }
          // Fluorine
          if( fixed_elements( e_i) == "F")
          {
            m_FixedElements.PushBack( GetElementTypes().e_Fluorine);
          }
          // Chlorine
          if( fixed_elements( e_i) == "Cl")
          {
            m_FixedElements.PushBack( GetElementTypes().e_Chlorine);
          }
          // Bromine
          if( fixed_elements( e_i) == "Br")
          {
            m_FixedElements.PushBack( GetElementTypes().e_Bromine);
          }
          // Iodine
          if( fixed_elements( e_i) == "I")
          {
            m_FixedElements.PushBack( GetElementTypes().e_Iodine);
          }
          // Undefined
          if( fixed_elements( e_i) == "X")
          {
            m_FixedElements.PushBack( GetElementTypes().e_Undefined);
          }
        }
      }

      //////// Ring library atom selection options ////////

      // read in mutable atom indices
      if( m_RingLibraryMutableAtoms.size())
      {
        m_RingLibraryMutableAtomindices.Reset();
        m_RingLibraryMutableAtomindices = util::SplitStringToNumerical< size_t>( m_RingLibraryMutableAtoms);
      }

      // read in allowed mutable elements when choosing mutable atom
      if( m_RingLibraryMutableElementsString.size())
      {
        // parse input
        const storage::Vector< std::string> mutable_elements
        (
          util::SplitString( util::TrimString( m_RingLibraryMutableElementsString), " \t\n\r,")
        );

        // stupid check to add only the correct elements
        // TODO: this should be directly serializable from element types to make concise
        m_RingLibraryMutableElements.Reset();
        for( size_t e_i( 0), e_sz( mutable_elements.GetSize()); e_i < e_sz; ++e_i)
        {
          // Hydrogen
          if( mutable_elements( e_i) == "H")
          {
            m_RingLibraryMutableElements.PushBack( GetElementTypes().e_Hydrogen);
          }
          // Boron
          if( mutable_elements( e_i) == "B")
          {
            m_RingLibraryMutableElements.PushBack( GetElementTypes().e_Boron);
          }
          // Carbon
          if( mutable_elements( e_i) == "C")
          {
            m_RingLibraryMutableElements.PushBack( GetElementTypes().e_Carbon);
          }
          // Oxygen
          if( mutable_elements( e_i) == "O")
          {
            m_RingLibraryMutableElements.PushBack( GetElementTypes().e_Oxygen);
          }
          // Nitrogen
          if( mutable_elements( e_i) == "N")
          {
            m_RingLibraryMutableElements.PushBack( GetElementTypes().e_Nitrogen);
          }
          // Potassium
          if( mutable_elements( e_i) == "P")
          {
            m_RingLibraryMutableElements.PushBack( GetElementTypes().e_Phosphorus);
          }
          // Sulfur
          if( mutable_elements( e_i) == "S")
          {
            m_RingLibraryMutableElements.PushBack( GetElementTypes().e_Sulfur);
          }
          // Selenium
          if( mutable_elements( e_i) == "Se")
          {
            m_RingLibraryMutableElements.PushBack( GetElementTypes().e_Selenium);
          }
          // Fluorine
          if( mutable_elements( e_i) == "F")
          {
            m_RingLibraryMutableElements.PushBack( GetElementTypes().e_Fluorine);
          }
          // Chlorine
          if( mutable_elements( e_i) == "Cl")
          {
            m_RingLibraryMutableElements.PushBack( GetElementTypes().e_Chlorine);
          }
          // Bromine
          if( mutable_elements( e_i) == "Br")
          {
            m_RingLibraryMutableElements.PushBack( GetElementTypes().e_Bromine);
          }
          // Iodine
          if( mutable_elements( e_i) == "I")
          {
            m_RingLibraryMutableElements.PushBack( GetElementTypes().e_Iodine);
          }
          // Undefined
          if( mutable_elements( e_i) == "X")
          {
            m_RingLibraryMutableElements.PushBack( GetElementTypes().e_Undefined);
          }
        }
      }

      // read in reference filename
      s_Mutex.Lock();
      if( m_RingLibraryFixedFragmentsFilename.size())
      {
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_RingLibraryFixedFragmentsFilename);
        m_RingLibraryFixedFragments.ReadMoreFromMdl( input);
        io::File::CloseClearFStream( input);
      }
      s_Mutex.Unlock();

      // read in mutable atom indices
      if( m_RingLibraryFixedAtoms.size())
      {
        m_RingLibraryFixedAtomindices.Reset();
        m_RingLibraryFixedAtomindices = util::SplitStringToNumerical< size_t>( m_RingLibraryFixedAtoms);
      }

      // read in allowed mutable elements when choosing mutable atom
      if( m_RingLibraryFixedElementsString.size())
      {
        // parse input
        const storage::Vector< std::string> fixed_elements
        (
          util::SplitString( util::TrimString( m_RingLibraryFixedElementsString), " \t\n\r,")
        );

        // check to add only the correct elements
        m_RingLibraryFixedElements.Reset();
        for( size_t e_i( 0), e_sz( fixed_elements.GetSize()); e_i < e_sz; ++e_i)
        {
          // Hydrogen
          if( fixed_elements( e_i) == "H")
          {
            m_RingLibraryFixedElements.PushBack( GetElementTypes().e_Hydrogen);
          }
          // Boron
          if( fixed_elements( e_i) == "B")
          {
            m_RingLibraryFixedElements.PushBack( GetElementTypes().e_Boron);
          }
          // Carbon
          if( fixed_elements( e_i) == "C")
          {
            m_RingLibraryFixedElements.PushBack( GetElementTypes().e_Carbon);
          }
          // Oxygen
          if( fixed_elements( e_i) == "O")
          {
            m_RingLibraryFixedElements.PushBack( GetElementTypes().e_Oxygen);
          }
          // Nitrogen
          if( fixed_elements( e_i) == "N")
          {
            m_RingLibraryFixedElements.PushBack( GetElementTypes().e_Nitrogen);
          }
          // Potassium
          if( fixed_elements( e_i) == "P")
          {
            m_RingLibraryFixedElements.PushBack( GetElementTypes().e_Phosphorus);
          }
          // Sulfur
          if( fixed_elements( e_i) == "S")
          {
            m_RingLibraryFixedElements.PushBack( GetElementTypes().e_Sulfur);
          }
          // Selenium
          if( fixed_elements( e_i) == "Se")
          {
            m_RingLibraryFixedElements.PushBack( GetElementTypes().e_Selenium);
          }
          // Fluorine
          if( fixed_elements( e_i) == "F")
          {
            m_RingLibraryFixedElements.PushBack( GetElementTypes().e_Fluorine);
          }
          // Chlorine
          if( fixed_elements( e_i) == "Cl")
          {
            m_RingLibraryFixedElements.PushBack( GetElementTypes().e_Chlorine);
          }
          // Bromine
          if( fixed_elements( e_i) == "Br")
          {
            m_RingLibraryFixedElements.PushBack( GetElementTypes().e_Bromine);
          }
          // Iodine
          if( fixed_elements( e_i) == "I")
          {
            m_RingLibraryFixedElements.PushBack( GetElementTypes().e_Iodine);
          }
          // Undefined
          if( fixed_elements( e_i) == "X")
          {
            m_RingLibraryFixedElements.PushBack( GetElementTypes().e_Undefined);
          }
        }
      }

      // done
      return true;
    }

  } // namespace chemistry
} // namespace bcl
