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
#include "chemistry/bcl_chemistry_fragment_mutate_remove_bond.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_voxel_grid_atom.h"
#include "command/bcl_command_command_state.h"
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

    // add the interfaces to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentMutateRemoveBond::s_AddBondInstance
    (
      util::Enumerated< FragmentMutateInterface>::AddInstance( new FragmentMutateRemoveBond( FragmentMutateRemoveBond::BondTreatment::e_AddBond))
    );
    const util::SiPtr< const util::ObjectInterface> FragmentMutateRemoveBond::s_RemoveBondInstance
    (
      util::Enumerated< FragmentMutateInterface>::AddInstance( new FragmentMutateRemoveBond( FragmentMutateRemoveBond::BondTreatment::e_RemoveBond))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentMutateRemoveBond::FragmentMutateRemoveBond() :
      m_BondChange( FragmentMutateRemoveBond::BondTreatment::e_RemoveBond),
      m_BondType( GetConfigurationalBondTypes().e_NonConjugatedSingleBond),
      m_PairedAtomIndices( storage::Vector< size_t>()),
      m_PairedAtoms( std::string())
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief bond change constructor
    //! @param BOND_CHANGE whether to add or remove bond
    FragmentMutateRemoveBond::FragmentMutateRemoveBond( const BondTreatment &BOND_CHANGE) :
      m_BondChange( BOND_CHANGE),
      m_BondType( GetConfigurationalBondTypes().e_NonConjugatedSingleBond),
      m_PairedAtomIndices( storage::Vector< size_t>()),
      m_PairedAtoms( std::string())
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief druglikeness constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    FragmentMutateRemoveBond::FragmentMutateRemoveBond
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const bool &CORINA_CONFS
    ) :
      m_BondChange( FragmentMutateRemoveBond::BondTreatment::e_RemoveBond),
      m_BondType( GetConfigurationalBondTypes().e_NonConjugatedSingleBond),
      m_PairedAtomIndices( storage::Vector< size_t>()),
      m_PairedAtoms( std::string())
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
    FragmentMutateRemoveBond::FragmentMutateRemoveBond
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const bool &CORINA_CONFS
    ) :
      m_BondChange( FragmentMutateRemoveBond::BondTreatment::e_RemoveBond),
      m_BondType( GetConfigurationalBondTypes().e_NonConjugatedSingleBond),
      m_PairedAtomIndices( storage::Vector< size_t>()),
      m_PairedAtoms( std::string())
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
    FragmentMutateRemoveBond::FragmentMutateRemoveBond
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
      m_BondChange( FragmentMutateRemoveBond::BondTreatment::e_RemoveBond),
      m_BondType( GetConfigurationalBondTypes().e_NonConjugatedSingleBond),
      m_PairedAtomIndices( storage::Vector< size_t>()),
      m_PairedAtoms( std::string())
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
    FragmentMutateRemoveBond::FragmentMutateRemoveBond
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
      m_BondChange( FragmentMutateRemoveBond::BondTreatment::e_RemoveBond),
      m_BondType( GetConfigurationalBondTypes().e_NonConjugatedSingleBond),
      m_PairedAtomIndices( storage::Vector< size_t>()),
      m_PairedAtoms( std::string())
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
    FragmentMutateRemoveBond *FragmentMutateRemoveBond::Clone() const
    {
      return new FragmentMutateRemoveBond( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentMutateRemoveBond::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentMutateRemoveBond::GetAlias() const
    {
      static const std::string s_remove( "RemoveBond"), s_add( "AddBond");
      return m_BondChange == e_AddBond ? s_add : s_remove;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @return MutateResult with Constitution after the mutate
    math::MutateResult< FragmentComplete> FragmentMutateRemoveBond::operator()( const FragmentComplete &FRAGMENT) const
    {
      std::string name( "FragmentMutateRemoveBond!");
      if( m_BondChange == e_AddBond)
      {
        name = "FragmentAddBond!";
      }
      BCL_MessageStd( util::Format()( name));

      // try a few times to select a heavy atom
      for( size_t i( 0); i < m_NumberMaxAttempts; ++i)
      {

        // choose the first atom involved in the bond to be moved
        util::SiPtr< const AtomConformationalInterface> picked_atom;
        if( m_MutableAtomIndices.GetSize() || m_MutableElements.GetSize() || m_MutableFragments.GetSize())
        {
          picked_atom = this->PickAtom( FRAGMENT, false);
        }
        else
        {
          picked_atom = this->PickAtom( FRAGMENT, true);
        }

        // if atom is hydrogen atom, grab the atom to which it is connected
        if( picked_atom->GetElementType() == GetElementTypes().e_Hydrogen)
        {
          if( !picked_atom->GetBonds().GetSize())
          {
            continue;
          }
          picked_atom = util::SiPtr< const AtomConformationalInterface>( picked_atom->GetBonds().Begin()->GetTargetAtom());
        }

        // we now have the index of our chosen heavy atom
        size_t picked_atom_index( FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom));

        // remove bonds
        AtomVector< AtomComplete> atom_vector;
        if( m_BondChange == e_RemoveBond)
        {
          // Identify eligible bonds
          storage::Vector< size_t> bond_indices;
          size_t bond_index( 0);

          for
          (
              auto bonds_itr( FRAGMENT.GetAtomVector()( picked_atom_index).GetBonds().Begin()),
              bonds_itr_end( FRAGMENT.GetAtomVector()( picked_atom_index).GetBonds().End());
              bonds_itr != bonds_itr_end;
              ++bonds_itr, ++bond_index
          )
          {
            // restrict search to paired atoms
            if( m_PairedAtomIndices.GetSize())
            {
              m_PairedAtomIndices.GetSize();
              if( m_PairedAtomIndices.Find( FRAGMENT.GetAtomIndex( bonds_itr->GetTargetAtom())) < m_PairedAtomIndices.GetSize())
              {
                bond_indices.PushBack( bond_index);
              }
            }

            // otherwise look for any heavy atom
            else if( bonds_itr->GetTargetAtom().GetElementType() != GetElementTypes().e_Hydrogen)
            {
              bond_indices.PushBack( bond_index);
            }
          }

          // Pick a bond to remove
          if( !bond_indices.GetSize())
          {
            return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
          }
          else
          {
            bond_indices.Shuffle();
            bond_index = bond_indices( 0);
          }

          // Remove bond
          sdf::BondInfo bond_ab
          (
            picked_atom_index,
            FRAGMENT.GetAtomVector().GetAtomIndex( FRAGMENT.GetAtomVector()( picked_atom_index).GetBonds()( bond_index).GetTargetAtom()),
            FRAGMENT.GetAtomVector()( picked_atom_index).GetBonds()( bond_index).GetBondType()
          );
          atom_vector = RemoveBond( FRAGMENT, bond_ab);
        }
        else if( m_BondChange == e_AddBond)
        {
          // pick an atom from paired atom list
          size_t second_atom_index( util::GetUndefinedSize_t());
          if( m_PairedAtomIndices.GetSize())
          {
            second_atom_index = m_PairedAtomIndices( random::GetGlobalRandom().Random< size_t>( 0, m_PairedAtomIndices.GetSize() - 1));
          }
          storage::Triplet< FragmentComplete, size_t, size_t> pair_a( OpenValence( FRAGMENT, picked_atom_index, m_OVShuffleH, m_OVReverse));
          if( pair_a.Third() < second_atom_index)
          {
            second_atom_index -= size_t( 1);
          }
          storage::Triplet< FragmentComplete, size_t, size_t> pair_b( OpenValence( pair_a.First(), second_atom_index, m_OVShuffleH, m_OVReverse));
          if( pair_b.Third() < pair_a.Second())
          {
            pair_a.Second() -= size_t( 1);
          }

          // Add bond
          sdf::BondInfo bond_ab
          (
            pair_a.Second(),
            pair_b.Second(),
            m_BondType
          );
          atom_vector = AddBond( pair_b.First(), bond_ab);
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

        // Remove hydrogen atoms to allow bond type adjustment
        HydrogensHandler::Remove( atom_vector);
        if( m_ScaffoldFragment.GetSize())
        {
          return math::MutateResult< FragmentComplete>( cleaner.Clean( atom_vector, m_ScaffoldFragment, m_DrugLikenessType), *this);
        }
        else
        {
          return math::MutateResult< FragmentComplete>( cleaner.Clean( atom_vector, FRAGMENT, m_DrugLikenessType), *this);
        }
      }
      // failed all tries; return null
      return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief set the bond change type
    void FragmentMutateRemoveBond::SetBondChange( const BondTreatment &BOND_CHANGE)
    {
      m_BondChange = BOND_CHANGE;
    }

    //! @brief a function that removes a bond between two atoms
    //! @param FRAGMENT the small molecule of interest
    //! @param BOND the bond to remove
    //! @return atom vector after removal of the bond
    AtomVector< AtomComplete> FragmentMutateRemoveBond::RemoveBond( const FragmentComplete &FRAGMENT, const sdf::BondInfo &BOND) const
    {
      // remove bond
      storage::Vector< sdf::AtomInfo> atominfo( FRAGMENT.GetAtomInfo());
      storage::Vector< sdf::BondInfo> bondinfo;
      for
      (
          auto atom_itr( FRAGMENT.GetAtomVector().Begin()), atom_itr_end( FRAGMENT.GetAtomVector().End());
          atom_itr != atom_itr_end;
          ++atom_itr
      )
      {
        size_t index_a( FRAGMENT.GetAtomVector().GetAtomIndex( *atom_itr));
        size_t current_bond( 0);
        for
        (
            auto bond_itr( atom_itr->GetBonds().Begin()), bond_itr_end( atom_itr->GetBonds().End());
            bond_itr != bond_itr_end;
            ++bond_itr, ++current_bond
        )
        {
          // if the current atom is the original chosen atom, do not store bond information
          size_t index_b( FRAGMENT.GetAtomVector().GetAtomIndex( bond_itr->GetTargetAtom()));
          if( index_a > index_b)
          {
            continue;
          }
          sdf::BondInfo temp_bondinfo( index_a, index_b, bond_itr->GetBondType());
          if( !( temp_bondinfo == BOND))
          {
            bondinfo.PushBack( temp_bondinfo);
          }
        }
      }

      // Make new molecule
      return AtomVector< AtomComplete>( atominfo, bondinfo);
    }

    //! @brief a function that adds a bond between two atoms
    //! @param FRAGMENT the small molecule of interest
    //! @param BOND the bond to add
    //! @return atom vector after addition of the bond
    AtomVector< AtomComplete> FragmentMutateRemoveBond::AddBond( const FragmentComplete &FRAGMENT, const sdf::BondInfo &BOND) const
    {
      // add bond
      storage::Vector< sdf::AtomInfo> atominfo( FRAGMENT.GetAtomInfo());
      storage::Vector< sdf::BondInfo> bondinfo( FRAGMENT.GetBondInfo());
      bondinfo.PushBack( BOND);

      // Make new molecule
      return AtomVector< AtomComplete>( atominfo, bondinfo);
    }

  //////////////////////
  // helper functions //
  //////////////////////

    io::Serializer FragmentMutateRemoveBond::GetSerializer() const
    {
      io::Serializer parameters( FragmentMutateInterface::GetSerializer());
      parameters.SetClassDescription
      (
        m_BondChange == e_AddBond ?
        "Adds a bond to a molecule." :
        "Removes a bond from a molecule. If the bond removal splits the molecule into "
        "multiple non-bonded components, save only the largest component."
      );

      parameters.AddInitializer
      (
        "paired_atoms",
        "atom indices (0-indexed) of a second pool of atoms that can be mutated; "
        "this flag has no effect if there are no valid 'mutable_atoms' for "
        "the chosen mutate; defaults to all atoms in the molecule, which may "
        "have different implications for different mutates; for example, in "
        "RemoveBond the 'paired_atoms' flag can be specified along with "
        "'mutable_atoms' to delete a specific bond between two atoms",
        io::Serialization::GetAgent( &m_PairedAtoms),
        ""
      );

      parameters.AddInitializer
      (
        "bond_type",
        "bond type of the bond to be added if adding and not removing a bond",
        io::Serialization::GetAgent( &m_BondType),
        "NonConjugatedSingleBond"
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentMutateRemoveBond::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
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

      // read in paired atom indices
      if( m_PairedAtoms.size())
      {
        m_PairedAtomIndices = util::SplitStringToNumerical< size_t>( m_PairedAtoms);
      }

      // done
      return true;
    }

  } // namespace chemistry
} // namespace bcl
