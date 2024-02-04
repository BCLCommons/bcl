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
#include "chemistry/bcl_chemistry_fragment_mutate_cyclize.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_conformational_interface.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_map_conformer.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "chemistry/bcl_chemistry_fragment_track_mutable_atoms.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_base.h"
#include "graph/bcl_graph_connectivity.h"
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

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> FragmentMutateCyclize::s_Instance
    (
      util::Enumerated< FragmentMutateInterface>::AddInstance( new FragmentMutateCyclize())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FragmentMutateCyclize::FragmentMutateCyclize() :
        m_Rings( new ConstitutionSet()),
        m_RingsFilename( std::string())
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief constructor
    //! @param DRUG_LIKENESS_TYPE type of druglikeness filter to apply during clean
    //! @param SCAFFOLD_FRAGMENT fragment to which the new mutated molecule will be aligned based on substructure
    //! @param MUTABLE_FRAGMENTS non-mutable component of the current molecule
    //! @param MUTABLE_ATOM_INDICES indices of atoms that can be mutated
    FragmentMutateCyclize::FragmentMutateCyclize
    (
      const std::string &DRUG_LIKENESS_TYPE,
      const FragmentComplete &SCAFFOLD_FRAGMENT,
      const FragmentEnsemble &MUTABLE_FRAGMENTS,
      const storage::Vector< size_t> &MUTABLE_ATOM_INDICES,
      const bool &CORINA_CONFS
    ) :
      m_Rings( new ConstitutionSet()),
      m_RingsFilename( std::string())
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
    FragmentMutateCyclize::FragmentMutateCyclize
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
      m_Rings( new ConstitutionSet()),
      m_RingsFilename( std::string())
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
    FragmentMutateCyclize::FragmentMutateCyclize
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
      m_Rings( new ConstitutionSet()),
      m_RingsFilename( std::string())
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
    FragmentMutateCyclize *FragmentMutateCyclize::Clone() const
    {
      return new FragmentMutateCyclize( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FragmentMutateCyclize::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &FragmentMutateCyclize::GetAlias() const
    {
      static const std::string s_name( "Cyclize");
      return s_name;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an fragment and generating a new fragment by growing on a valence
    //! @param FRAGMENT small molecule of interest
    //! @return MutateResult with Constitution after the mutate
    math::MutateResult< FragmentComplete> FragmentMutateCyclize::operator()( const FragmentComplete &FRAGMENT) const
    {
      BCL_MessageStd( "FragmentMutateCyclize!");

      // get fragment graph
      ConformationGraphConverter graph_maker;
      graph::ConstGraph< size_t, size_t> fragment_graph( graph_maker( FRAGMENT));
      AtomVector< AtomComplete> atom_vector( FRAGMENT.GetAtomVector());

      // try a few times
      size_t picked_atom_index( util::GetUndefinedSize_t()), second_atom_index( util::GetUndefinedSize_t());
      for( size_t tries( 0); tries < m_NumberMaxAttempts; ++tries)
      {
        // pick an atom
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

        // two criteria: (1) must not be in a ring, and (2) must have a potential valence through loss of hydrogen atom
        if
        (
            picked_atom->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1)) ||
            !picked_atom->GetNumberCovalentlyBoundHydrogens()
        )
        {
          continue;
//          // if it is in a ring and is protonated, convert proton to carbon and change selected atom index to that carbon
//          if
//          (
//              picked_atom->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1)) &&
//              picked_atom->GetNumberCovalentlyBoundHydrogens()
//          )
//          {
//            for
//            (
//                auto bond_itr( picked_atom->GetBonds().Begin()), bond_itr_end( picked_atom->GetBonds().End());
//                bond_itr != bond_itr_end;
//                ++bond_itr
//            )
//            {
//              if( bond_itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Hydrogen)
//              {
//                atom_vector( atom_vector.GetAtomIndex( bond_itr->GetTargetAtom())).SetAtomType( GetAtomTypes().C_TrTrTrPi);
//
//
//              }
//            }
//          }
//
//          // TODO: if the picked atom does not have hydrogen atoms, see if a bond order can be reduced to protonate it
//          else if( !picked_atom->GetNumberCovalentlyBoundHydrogens())
//          {
//            continue;
//          }
        }

        picked_atom_index = FRAGMENT.GetAtomVector().GetAtomIndex( *picked_atom);

        // determine the distance from picked_atom vertex to all vertices
        util::ShPtr< storage::Vector< size_t> > bond_distances
        (
          graph::Connectivity::DistancesToOtherVertices( fragment_graph, picked_atom_index)
        );

        // make storage for atoms 3, 4, 5, and 6 bonds away
        storage::Vector< size_t> bond_distance_three, bond_distance_four, bond_distance_five, bond_distance_six;

        // collect atom indices corresponding to each bond distance
        size_t atom_index( 0);
        for
        (
            auto
            itr_distances( bond_distances->Begin()), itr_distances_end( bond_distances->End());
            itr_distances != itr_distances_end;
            ++itr_distances, ++atom_index
        )
        {
          // make sure the atom at each index is (1) not a hydrogen atom, and (2) bound to at least one hydrogen atom
          if
          (
              FRAGMENT.GetAtomVector()( atom_index).GetElementType() != GetElementTypes().e_Hydrogen &&
              FRAGMENT.GetAtomVector()( atom_index).GetNumberCovalentlyBoundHydrogens()
          )
          {
            // we only care about distances at 3, 4, 5, and 6 bonds
            if( *itr_distances == size_t( 3) || *itr_distances == size_t( 4) || *itr_distances == size_t( 5) || *itr_distances == size_t( 6))
            {
              auto path( graph::Connectivity::FindMinimalPath( fragment_graph, picked_atom_index, atom_index));
              size_t bond_ring_cnt( 0), bond_ring_cnt_aromatic( 0);
              for
              (
                  auto path_itr( path.Begin()), path_itr_nxt( ++path.Begin()), path_itr_end( path.End());
                  path_itr_nxt != path_itr_end;
                  ++path_itr, ++path_itr_nxt
              )
              {
                auto bond_type( FRAGMENT.GetAtomVector()( *path_itr).GetBondTypeTo( FRAGMENT.GetAtomVector()( *path_itr_nxt)));
                if( bond_type->IsBondInRing())
                {
                  // crosses two aromatic bonds max
                  if( bond_type->GetConjugation() == ConstitutionalBondTypeData::e_Aromatic)
                  {
                    if( ++bond_ring_cnt_aromatic > size_t( 2))
                    {
                      break;
                    }
                  }
                  // crosses one non-aromatic bond max
                  else if( ++bond_ring_cnt > size_t( 1))
                  {
                    break;
                  }
                }
              }
//              if( bond_ring_cnt > size_t( 1) || ( bond_ring_cnt_aromatic && bond_ring_cnt) || bond_ring_cnt_aromatic > size_t( 2))
              if( bond_ring_cnt > size_t( 1) || bond_ring_cnt_aromatic > size_t( 2))
              {
                continue;
              }
              if( *itr_distances == size_t( 3) && bond_ring_cnt_aromatic < size_t( 2)) // long-term needs a better fix
              {
                bond_distance_three.PushBack( atom_index);
              }
              else if( *itr_distances == size_t( 4) && bond_ring_cnt_aromatic < size_t( 2))
              {
                bond_distance_four.PushBack( atom_index);
              }
              else if( *itr_distances == size_t( 5) && bond_ring_cnt_aromatic < size_t( 3))
              {
                bond_distance_five.PushBack( atom_index);
              }
              else if( *itr_distances == size_t( 6) && bond_ring_cnt_aromatic < size_t( 3))
              {
                bond_distance_six.PushBack( atom_index);
              }
            }
          }
        }
        // decide which ring size to make
        float rand_prob( random::GetGlobalRandom().Random< float>( 0.0, 1.0));

        // six-membered rings are ~3x more common in drug-like molecules
        bool has_six( bond_distance_six.GetSize()),
             has_five( bond_distance_five.GetSize()),
             has_four( bond_distance_four.GetSize()),
             has_three( bond_distance_three.GetSize());
        if( rand_prob < 0.75 && has_five)
         {
          // six
          second_atom_index = bond_distance_five( random::GetGlobalRandom().Random< size_t>( 0, bond_distance_five.GetSize() - 1));
          break;
        }
        else if( rand_prob < 0.75 && has_three)
        {
          second_atom_index = bond_distance_three( random::GetGlobalRandom().Random< size_t>( 0, bond_distance_three.GetSize() - 1));
          if( !FRAGMENT.GetAtomVector()( second_atom_index).GetNumberCovalentlyBoundHydrogens())
          {
            continue;
          }
          else
          {
            // change hydrogen atoms to carbon atoms and cyclize the two new carbon atoms to make a new six-membered ring
            // second atom
            for
            (
                auto bond_itr( atom_vector( second_atom_index).GetBonds().Begin()), bond_itr_end( atom_vector( second_atom_index).GetBonds().End());
                bond_itr != bond_itr_end;
                ++bond_itr
            )
            {
              if( bond_itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Hydrogen)
              {
                atom_vector( atom_vector.GetAtomIndex( bond_itr->GetTargetAtom())).SetAtomType( GetAtomTypes().C_TrTrTrPi);
                second_atom_index = atom_vector.GetAtomIndex( bond_itr->GetTargetAtom());
                break;
              }
            }

            // picked atom
            for
            (
                auto bond_itr( atom_vector( picked_atom_index).GetBonds().Begin()), bond_itr_end( atom_vector( picked_atom_index).GetBonds().End());
                bond_itr != bond_itr_end;
                ++bond_itr
            )
            {
              if( bond_itr->GetTargetAtom().GetElementType() == GetElementTypes().e_Hydrogen)
              {
                atom_vector( atom_vector.GetAtomIndex( bond_itr->GetTargetAtom())).SetAtomType( GetAtomTypes().C_TrTrTrPi);
                picked_atom_index = atom_vector.GetAtomIndex( bond_itr->GetTargetAtom());
                break;
              }
            }
          }
          break;
        }
        else if( has_four)
        {
          // five
          second_atom_index = bond_distance_four( random::GetGlobalRandom().Random< size_t>( 0, bond_distance_four.GetSize() - 1));
          break;
        }
        else if( has_six)
        {
          // seven
          second_atom_index = bond_distance_six( random::GetGlobalRandom().Random< size_t>( 0, bond_distance_six.GetSize() - 1));
          break;
        }
      }
      if( !util::IsDefined( picked_atom_index) || !util::IsDefined( second_atom_index))
      {
        util::ShPtr< FragmentComplete> new_mol_ptr( new FragmentComplete( FRAGMENT));
        return math::MutateResult< FragmentComplete>( new_mol_ptr, *this);
      }

      // Create new bond
      storage::Vector< sdf::AtomInfo> atominfo( atom_vector.GetAtomInfo());
      storage::Vector< sdf::BondInfo> bonds( atom_vector.GetBondInfo());

      //float rand( random::GetGlobalRandom().Random< float>( 0.0, 1.0));
      if
      (
          atom_vector( picked_atom_index).GetAtomType()->GetNumberBonds() - atom_vector( picked_atom_index).GetBonds().GetSize() > 1 &&
          atom_vector( second_atom_index).GetAtomType()->GetNumberBonds() - atom_vector( second_atom_index).GetBonds().GetSize() > 1
      )
      {
        bonds.PushBack( sdf::BondInfo( picked_atom_index, second_atom_index, GetConfigurationalBondTypes().e_ConjugatedDoubleBondInRing));
      }
      else
      {
        bonds.PushBack( sdf::BondInfo( picked_atom_index, second_atom_index, GetConfigurationalBondTypes().e_ConjugatedSingleBondInRing));
      }
      AtomVector< AtomComplete> not_empty( atominfo, bonds);

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

      // Remove hydrogen atoms before clean to allow proper bondtype selection
      HydrogensHandler::Remove( not_empty);

      // Check for valid atom types
      util::ShPtr< FragmentComplete> new_mol_ptr
      (
        m_ScaffoldFragment.GetSize()
         ? cleaner.Clean( not_empty, m_ScaffoldFragment, m_DrugLikenessType)
         : cleaner.Clean( not_empty, FRAGMENT, m_DrugLikenessType)
      );

      if( !new_mol_ptr.IsDefined() || new_mol_ptr->HasNonGasteigerAtomTypes())
      {
        return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
      }

      // split out rings
      FragmentSplitRings ring_splitter( true, 4);
      FragmentEnsemble split_rings( ring_splitter( *new_mol_ptr));

      // make sure all rings are found in the ring dataset
//      for
//      (
//          FragmentEnsemble::iterator split_rings_itr( split_rings.Begin()), split_rings_itr_end( split_rings.End());
//          split_rings_itr != split_rings_itr_end;
//          ++split_rings_itr
//      )
//      {
//        if( m_Rings->Find( FragmentConstitutionShared( *split_rings_itr)) == m_Rings->End())
//        {
//          BCL_MessageStd("Rings not all in the dataset \n");
//          return math::MutateResult< FragmentComplete>( util::ShPtr< FragmentComplete>(), *this);
//        }
//      }
      return math::MutateResult< FragmentComplete>( new_mol_ptr, *this);
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    io::Serializer FragmentMutateCyclize::GetSerializer() const
    {
      io::Serializer parameters( FragmentMutateInterface::GetSerializer());
      parameters.SetClassDescription
      (
        "Forms intramolecular bonds between two non-ring atoms, "
        "or between one ring and one non-ring atom"
      );

      parameters.AddInitializer
      (
        "ring_library",
        "path to the ring library",
        io::Serialization::GetAgent( &m_RingsFilename),
        ""
      );

      return parameters;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool FragmentMutateCyclize::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
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
      if( m_RingsFilename.size())
      {
        s_Mutex.Lock();
        io::IFStream input;
        io::File::MustOpenIFStream( input, m_RingsFilename);
        FragmentEnsemble ensemble( input, sdf::e_Remove);
        io::File::CloseClearFStream( input);
        for
        (
            auto itr_ensemble( ensemble.Begin()), itr_ensemble_end( ensemble.End());
            itr_ensemble != itr_ensemble_end;
            ++itr_ensemble
        )
        {
          m_Rings->Insert( FragmentConstitutionShared( *itr_ensemble));
        }
        s_Mutex.Unlock();
        BCL_Assert( m_Rings->GetSize(), "Ring library is empty!");
      }

      // done
      return true;
    }

  } // namespace chemistry
} // namespace bcl
