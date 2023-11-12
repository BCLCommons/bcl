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

// (c)   The BCL copyright and license yields to non-BCL copyrights and licenses where indicated by code comments.
// (c) (for academic users) or bcl-support-commercial@meilerlab.org (for commercial users)
// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "descriptor/bcl_descriptor_molecule_total_toxic_fragments.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_fragment_split_rigid.h"
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "chemistry/bcl_chemistry_rotamer_library_file.h"
#include "command/bcl_command_command_state.h"
#include "descriptor/bcl_descriptor_atom_sigma_charge.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "graph/bcl_graph_common_subgraph_isomorphism.h"
#include "graph/bcl_graph_const_graph.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_ifstream.h"
#include "io/bcl_io_serialization.h"
#include "sdf/bcl_sdf_bond_info.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    // initialize static
    sched::Mutex &MoleculeTotalToxicFragments::GetMutex()
    {
      static sched::Mutex s_Mutex;
      return s_Mutex;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    MoleculeTotalToxicFragments::MoleculeTotalToxicFragments() :
        m_StabilizingFragments( chemistry::FragmentEnsemble())
    {
      this->ReadInitializerSuccessHook( util::ObjectDataLabel(), util::GetLogger());
    }

    //! @brief virtual copy constructor
    //! @return pointer to new MoleculeTotalToxicFragments
    MoleculeTotalToxicFragments *MoleculeTotalToxicFragments::Clone() const
    {
      return new MoleculeTotalToxicFragments( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &MoleculeTotalToxicFragments::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the name of the property without any parameters
    //! @return name of the property as string
    const std::string &MoleculeTotalToxicFragments::GetAlias() const
    {
      static const std::string s_label( "MoleculeTotalToxicFragments");
      return s_label;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void MoleculeTotalToxicFragments::Calculate( linal::VectorReference< float> &STORAGE)
    {
      // Get our molecule
      util::SiPtr< const chemistry::ConformationInterface> current_mol( this->GetCurrentObject());
      chemistry::FragmentComplete molecule( *current_mol);

      // Count all toxic groups
      size_t n_aromatic_amines( CountAromaticAmines( molecule));
      size_t n_aromatic_nitro( CountAromaticNitro( molecule));
      size_t n_aliphaticHalide( CountAliphaticHalide( molecule));
      size_t n_three_mem_heterocycle( CountThreeMemberHeterocycle( molecule));
      size_t n_nitroso( CountNitroso( molecule));
      size_t n_unsubstituted_heteratoms( CountUnsubstitutedHeteroAtoms( molecule));
      size_t n_azo( CountAzoType( molecule));
      size_t n_polycyclic_arom( CountPolycyclicAromatic( molecule));

      size_t final_count
      (
        n_aromatic_amines +
        n_aromatic_nitro +
        n_aliphaticHalide +
        n_three_mem_heterocycle +
        n_nitroso +
        n_unsubstituted_heteratoms +
        n_azo +
        n_polycyclic_arom
      );

      final_count ?
          STORAGE( 0) = float( 1) :
          STORAGE( 0) = float( 0);
    } // Recalculate

    //! @brief print out the result for each mutagenic fragment filter
    void MoleculeTotalToxicFragments::PrintOutput( const storage::Vector< size_t> results)const
    {
      storage::Vector< std::string> frag_names;
      frag_names.Append("n_aromatic_amines");
      frag_names.Append("n_aromatic_nitro");
      frag_names.Append("n_aliphaticHalide");
      frag_names.Append("n_three_mem_heterocycle");
      frag_names.Append("n_nitroso");
      frag_names.Append("n_unsubstituted_heteratoms");
      frag_names.Append("n_azo");
      frag_names.Append("n_polycyclic_arom");

      for( size_t i( 0); i < frag_names.GetSize(); ++i)
      {
        BCL_MessageStd(util::Format()(frag_names(i)));
        BCL_MessageStd(util::Format()(results(i)));
        BCL_MessageStd("\n");
      }
    }

    //! @brief calculate the total number of aromatic amines
    //! @param MOLECULE the molecule of interest
    size_t MoleculeTotalToxicFragments::CountAromaticAmines( const chemistry::FragmentComplete &MOLECULE) const
    {
      // evaluate on per-atom basis
      const chemistry::AtomVector< chemistry::AtomComplete> &atom_v( MOLECULE.GetAtomVector());

      /// get the index of nitrogen atom by going through the vector
      size_t num_aromatic_amines( 0);
      for( size_t i( 0), sz( atom_v.GetSize()); i < sz; ++i)
      {
        // check to see if the atom is a nitrogen with a single heavy atom substitution
        if
        (
            atom_v( i).GetElementType() == chemistry::GetElementTypes().e_Nitrogen && // is N
            CountNumConnectedHeavyAtoms( atom_v( i)) == size_t( 1) // single heavy atom substitution
        )
        {
          // iterate over bonds of nitrogen atom
          for
          (
              storage::Vector< chemistry::BondConformational>::const_iterator
              bond_itr( atom_v( i).GetBonds().Begin()), bond_itr_end( atom_v( i).GetBonds().End());
              bond_itr != bond_itr_end;
              ++bond_itr
          )
          {
            // increment aromatic amine count if bonded heavy atom is aromatic
            if
            (
                bond_itr->GetTargetAtom().CountNonValenceBondsWithProperty
                (
                  chemistry::ConfigurationalBondTypeData::e_IsAromatic,
                  size_t( 1)
                )
            )
            {
              ++num_aromatic_amines;
            }
          }
        }
      }

      return num_aromatic_amines;
    }

    //! @brief calculate the total number of aromatic nitro groups
    //! @param MOLECULE the molecule of interest
    size_t MoleculeTotalToxicFragments::CountAromaticNitro( const chemistry::FragmentComplete &MOLECULE) const
    {
      // evaluate on per-atom basis
      const chemistry::AtomVector< chemistry::AtomComplete> &atom_v( MOLECULE.GetAtomVector());

      // a counter for aromatic nitro substructure
      size_t num_aromatic_nitro( 0);
      for( size_t i( 0); i < atom_v.GetSize(); ++i)
      {
        // check to see if the selected atom meets basic criteria
        if
        (
            atom_v( i).GetElementType() == chemistry::GetElementTypes().e_Nitrogen && // is N
            CountNumConnectedHeavyAtoms( atom_v( i)) > size_t( 1) // to qualify must have 2 bonded O atoms
        )
        {
          // get the bonds off of the nitrogen atom
          storage::Vector< chemistry::BondConformational> bonds( atom_v( i).GetBonds());

          // vectors to store the indices of oxygen and aromatic carbons
          size_t num_oxygen( 0);
          size_t n_aro( 0);

          // evaluate bonded atom for each bond
          for( size_t j = 0; j < bonds.GetSize(); ++j)
          {
            // get the atom indices of any oxygen atoms bonded to our nitrogen
            if
            (
                bonds( j).GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Oxygen && // is O
                CountNumConnectedHeavyAtoms( bonds( j).GetTargetAtom()) == size_t( 1) // only O heavy bond is to N
            )
            {
              ++num_oxygen;
            }

            // determine if there are aromatic atoms attached to our nitrogen
            if
            (
                bonds( j).GetTargetAtom().CountNonValenceBondsWithProperty
                (
                  chemistry::ConfigurationalBondTypeData::e_IsAromatic,
                  size_t( 1))
            )
            {
              ++n_aro;
            }
          }

          // nitrogen bonded to an aromatic atom and two otherwise unsubstituted oxygen atoms
          if( num_oxygen == size_t( 2) && n_aro)
          {
            ++num_aromatic_nitro;
          }
        }
      }

      return num_aromatic_nitro;

    }

    //! @brief calculate the total number of aliphatic halides
    //! @param MOLECULE the molecule of interest
    //! TODO: later add options to detect more complex conditons ( 2nd/ 3rd/ stabilizing groups exclusion
    //! see Chemical structure, Salmonella mutagenicity and extent of carcinogenicity as indicators of
    //! genotoxic carcinogenesis among 222 chemicals tested in rodents by the U.S. NCI/NTP p.72
    size_t MoleculeTotalToxicFragments::CountAliphaticHalide( const chemistry::FragmentComplete &MOLECULE) const
    {
      // evaluate on per-atom basis
      size_t num_aliphatic_halide( 0);
      const chemistry::AtomVector< chemistry::AtomComplete> &atom_v( MOLECULE.GetAtomVector());
      for( size_t i = 0; i < atom_v.GetSize(); ++i)
      {
        // check if atom is a reactive halogen
        if
        (
            atom_v( i).GetElementType() == chemistry::GetElementTypes().e_Chlorine ||
            atom_v( i).GetElementType() == chemistry::GetElementTypes().e_Bromine ||
            atom_v( i).GetElementType() == chemistry::GetElementTypes().e_Iodine ||
            atom_v( i).GetElementType() == chemistry::GetElementTypes().e_Astatine
        )
        {
          // halogens can only make one bond
          //storage::Vector< chemistry::BondConformational>::const_iterator bond_itr( atom_v( i).GetBonds().Begin());
          for
          ( auto bond_itr( atom_v( i).GetBonds().Begin()),
              bond_itr_end( atom_v( i).GetBonds().End());
              bond_itr != bond_itr_end;
              ++bond_itr
          )
          {
            // if not attached to an aromatic then we'll call it an aliphatic halide for now
            if
            (
                !bond_itr->GetTargetAtom().CountNonValenceBondsWithProperty
                (
                  chemistry::ConfigurationalBondTypeData::e_IsAromatic,
                  size_t( 1)
                ) &&
                ( 2 == CountNumConnectedHeavyAtoms( bond_itr->GetTargetAtom()))
            )
            {
              ++num_aliphatic_halide;
            }

          }//end for the bond itr for loop
        }
      }

      return num_aliphatic_halide;
    }

    //! @brief calculate the total number of three-membered heterocycles
    //! @param MOLECULE the molecule of interest
    //TODO: revise to check each non-carb heavy atom and neighboring connections
    size_t MoleculeTotalToxicFragments::CountThreeMemberHeterocycle( const chemistry::FragmentComplete &MOLECULE) const
    {
      // counter for three-membered heterocycles
      size_t tmh_count( 0);
      // split out rings
      chemistry::FragmentSplitRings split_rings( true, 2);
      chemistry::FragmentEnsemble rings( split_rings( MOLECULE));

      // iterate over rings
      for
      (

          chemistry::FragmentEnsemble::const_iterator ring_itr( rings.Begin()), ring_itr_end( rings.End());
          ring_itr != ring_itr_end;
          ++ring_itr

      )
      {

        // for each ring structure , go through each molecule to search for hetero atoms
        //for each hetero atom, go through neighbor recursively to search for itself and record num bond
        //if reach self within 3 bonds return ture
        storage::Vector< bool> visited( ring_itr->GetSize(), false);
        const chemistry::AtomVector< chemistry::AtomComplete> atom_v( ring_itr->GetAtomVector());

        for
        (

            //go through each atom to find hetero atoms
            auto a_itr( atom_v.Begin()), a_itr_end( atom_v.End());
            a_itr != a_itr_end;
            ++a_itr

        )
        {

          if
          (
              a_itr -> GetElementType() != chemistry::GetElementTypes().e_Carbon &&
              a_itr -> GetElementType() != chemistry::GetElementTypes().e_Hydrogen &&
              a_itr -> GetElementType() != chemistry::GetElementTypes().e_Sulfur

          )
          {

            // search three member heterocycle from the current hetero atom
            tmh_count +=
                SearchThreeMemberCycle
                (

                  atom_v,
                  a_itr->GetBonds().Begin()->GetTargetAtom(),
                  *a_itr,
                  size_t( 1)

                );
          }

        }

      }
      return tmh_count;
    }

    //! @brief calculate the total number of nitroso groups
    //! @param MOLECULE the molecule of interest
    size_t MoleculeTotalToxicFragments::CountNitroso( const chemistry::FragmentComplete &MOLECULE) const
    {
      size_t num_nitroso( 0);
      const chemistry::AtomVector< chemistry::AtomComplete> &atom_v( MOLECULE.GetAtomVector());

      // find N atoms in the molecule
      for
      (
          chemistry::AtomVector< chemistry::AtomComplete>::const_iterator
          atom_itr( atom_v.Begin()), atom_end( atom_v.End());
          atom_itr != atom_end;
          ++atom_itr
      )
      {
        // kekulized aromatic double bonds are not nitroso groups
        if
        (
            atom_itr->GetElementType() == chemistry::GetElementTypes().e_Nitrogen &&
            !atom_itr->CountNonValenceBondsWithProperty
            (
              chemistry::ConfigurationalBondTypeData::e_IsAromatic,
              size_t( 1)
            )
        )
        {

          //should only have 2 connections around this N for nitroso
          if( atom_itr->GetBonds().GetSize() != 2)
          {
            continue;
          }

          //check if there is one oxygen bonded
          size_t num_oxygen( 0);

          for
          (
              storage::Vector< chemistry::BondConformational>::const_iterator
              bond_itr( atom_itr->GetBonds().Begin()), bond_end( atom_itr->GetBonds().End());
              bond_itr != bond_end;
              ++bond_itr
          )
          {
            //also check if the oxygen is not connected with other heavy atoms
            if
            (
                bond_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Oxygen && // is O
                bond_itr->GetBondType() == chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond && // is double bond
                CountNumConnectedHeavyAtoms( bond_itr->GetTargetAtom()) == size_t( 1) // ignore charged ring species for now
            )
            {
              ++num_oxygen;
            }

          }

          // incrememt nitroso count
          if( num_oxygen   )
          {
            ++num_nitroso;
          }
        }
      }
      return num_nitroso;
    }

    //! @brief calculate the total number of unsubstituted heteroatoms
    //! @param MOLECULE the molecule of interest
    // and figure out the db/single bond representation ( how can 2 nitrogen be counted in nitroso.sdf example?
    size_t MoleculeTotalToxicFragments::CountUnsubstitutedHeteroAtoms( const chemistry::FragmentComplete &MOLECULE) const
    {
      // initialize counter for unsubstituted heteroatom pairs
      size_t num_unsub_heteros( 0);

      // operate on per atom basis
      const chemistry::AtomVector< chemistry::AtomComplete> &atom_v( MOLECULE.GetAtomVector());
      for
      (
          chemistry::AtomVector< chemistry::AtomComplete>::const_iterator
          itr( atom_v.Begin()), itr_end( atom_v.End());
          itr != itr_end;
          ++itr
      )
      {
        //find a hetero atom
        if
        (

            itr->GetElementType() != chemistry::GetElementTypes().e_Carbon && // must be heteroatom
            itr->GetElementType() != chemistry::GetElementTypes().e_Hydrogen && // must be heavy
            CountNumConnectedHeavyAtoms( *itr) == size_t( 1) // must be the unsubstituted heteroatom

        )
        {
          // iterate over bonds of our heteroatom
          for
          (
              storage::Vector< chemistry::BondConformational>::const_iterator
              b_itr( itr->GetBonds().Begin()), b_itr_end( itr->GetBonds().End());
              b_itr != b_itr_end;
              ++b_itr
          )
          {
            // bonded heavy atom must be another heteroatom

            if
            (
                b_itr->GetTargetAtom().GetElementType() != chemistry::GetElementTypes().e_Carbon && // must be heteroatom
                b_itr->GetTargetAtom().GetElementType() != chemistry::GetElementTypes().e_Hydrogen && // must be heavy
                ( b_itr->GetBondType() == chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond ||
                    b_itr->GetBondType() == chemistry::GetConfigurationalBondTypes().e_ConjugatedSingleBond)

                    // must be single bond between the two heteroatoms

            )
            {
              //check P-O and S-O single bonds
              if
              (
                  ( itr->GetElementType() == chemistry::GetElementTypes().e_Phosphorus &&
                      b_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Oxygen) ||
                      ( itr->GetElementType() == chemistry::GetElementTypes().e_Oxygen &&
                          b_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Phosphorus) ||
                          ( itr->GetElementType() == chemistry::GetElementTypes().e_Sulfur &&
                              b_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Oxygen) ||
                              ( itr->GetElementType() == chemistry::GetElementTypes().e_Oxygen &&
                                  b_itr->GetTargetAtom().GetElementType() == chemistry::GetElementTypes().e_Sulfur)

              )
              {
                continue;
              }
              ++num_unsub_heteros;
            }
          }

        }
      }
      return num_unsub_heteros;
    }

    //! @brief calculate the total number of azo groups
    //! @param MOLECULE the molecule of interest
    size_t MoleculeTotalToxicFragments::CountAzoType( const chemistry::FragmentComplete &MOLECULE) const
    {
      // initialize azo counter
      size_t num_azo( 0);

      //store index of visited nitrogens to avoid repeated checking and counting
      storage::Vector< size_t> visited_nitrogen;
      const chemistry::AtomVector< chemistry::AtomComplete> &atom_v( MOLECULE.GetAtomVector());

      // iterate over bondinfo
      const storage::Vector< sdf::BondInfo> &bond_info( atom_v.GetBondInfo());
      for
      (
          auto bond_info_itr( bond_info.Begin()), bond_info_itr_end( bond_info.End());
          bond_info_itr != bond_info_itr_end;
          ++bond_info_itr
      )
      {
        // check for double bonded nitrogen atoms
        if
        (
            atom_v( bond_info_itr->GetAtomIndexLow()).GetElementType() == chemistry::GetElementTypes().e_Nitrogen &&
            atom_v( bond_info_itr->GetAtomIndexHigh()).GetElementType() == chemistry::GetElementTypes().e_Nitrogen &&
            bond_info_itr->GetConfigurationalBondType() == chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond &&
            //check if both atoms are not aromatic
            !atom_v( bond_info_itr->GetAtomIndexLow()).CountNonValenceBondsWithProperty
            (
              chemistry::ConfigurationalBondTypeData::e_IsAromatic,
              size_t( 1)
            ) &&
            !atom_v( bond_info_itr->GetAtomIndexHigh()).CountNonValenceBondsWithProperty
            (
              chemistry::ConfigurationalBondTypeData::e_IsAromatic,
              size_t( 1)
            )
        )
        {
          ++num_azo;
        }
      }
      return num_azo;
    }

    //! @brief calculate and return the total number of ring structures that contains polycyclic aromatic substructure
    //! @param MOLECULE the molecule of interest
    size_t MoleculeTotalToxicFragments::CountPolycyclicAromatic( const chemistry::FragmentComplete &MOLECULE) const
    {

      size_t num_polycyclic_arom( 0); //a counter for connected rings that contain polycyclic aromatic substructures

      //split the molecules and get a list of ring substructures
      chemistry::FragmentSplitRings split_rings( true, 2);
      chemistry::FragmentEnsemble rings( split_rings( MOLECULE)); // a list of FragmentComplete

      for
      (
          auto ring_itr( rings.Begin()), ring_end( rings.End());
          ring_itr != ring_end;
          ++ring_itr

      ) //iterating through each ring substructures
      {
        std::vector< bool> visited( rings.GetSize(), false); //a list of bools to record if each atom within this ring is visited

        //call the recursive helper function on each atom to check for interface atoms
        for
        (
            auto atom_itr( ring_itr->GetAtomVector().Begin()),
            atom_end( ring_itr->GetAtomVector().End());
            atom_itr != atom_end;
            ++atom_itr
        )
        {
          // a counter to count interface atoms
          size_t n_interface_atom
          (
            CountAromaticInterfaceAtoms
            (
              *ring_itr,
              *atom_itr,
              visited
            )
          );
          if( n_interface_atom >= 4)
          {
            //if more than 4 interface atoms within this ring substructure, it contains at least 3 connected sub-rings in a fused aromatic system
            ++num_polycyclic_arom;
          }
        }
      }
      return num_polycyclic_arom;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Count the number of heavy atoms connected to a given atom
    //! @param ATOM the atom of interest
    //! @return the number of heavy atoms bonded to the atom of interest
    size_t MoleculeTotalToxicFragments::CountNumConnectedHeavyAtoms( const chemistry::AtomConformationalInterface &ATOM) const
    {
      return ATOM.GetBonds().GetSize() - ATOM.GetNumberCovalentlyBoundHydrogens();
    }

    //! @brief Helper function to recursively count the number of interface atoms between 2 aromatic rings in a given connected ring substructure
    //! @param sub_ring the connected ring substructure of a molecule
    //! @param curr_atom a reference of a atom that is being visited in current call
    //! @param visited: a list of boolean keeping track of if the atoms in the structure has been visited in the recursive calls or not
    //! @return the number of interface atom in the Ring substructure

    size_t MoleculeTotalToxicFragments::CountAromaticInterfaceAtoms
    (
      const chemistry::FragmentComplete &SUB_RING,
      const chemistry::AtomConformationalInterface &CURR_ATOM,
      std::vector< bool> &VISITED
    ) const
    {
      //counter for atoms that are at interfaces between aromatic systems, indicated by number of aromatic atoms connected
      size_t num_interface_atom( 0);

      //number of aromatic atoms connected to the current atom, if >2, the current atom is on an interface
      size_t num_connected_arm( 0);

      if( VISITED[ SUB_RING.GetAtomIndex( CURR_ATOM)])
      {
        return 0;
      }
      else
      {
        VISITED[ SUB_RING.GetAtomIndex( CURR_ATOM)] = true;
      }

      for
      (
          storage::Vector< chemistry::BondConformational>::const_iterator
          b_itr( CURR_ATOM.GetBonds().Begin()), b_itr_end( CURR_ATOM.GetBonds().End());
          b_itr != b_itr_end;
          ++b_itr
      )
        // iterate through bonds of the current atom
      {
        if
        (
            //check if neighbor atoms are aromatic
            b_itr->GetTargetAtom().CountNonValenceBondsWithProperty
            (
              chemistry::ConfigurationalBondTypeData::e_IsAromatic,
              size_t( 1)
            )
        )
        {
          ++num_connected_arm;

          //recursively count interface atoms among neighbor atoms if they are aromatic
          num_interface_atom += CountAromaticInterfaceAtoms( SUB_RING, b_itr->GetTargetAtom(), VISITED);
        }
      }

      if( num_connected_arm > 2)
      {
        //after iterating through neighbor atoms, check if the current atom is on an interface
        ++num_interface_atom;
      }
      return num_interface_atom;
    }

    //! @brief Helper function to recursively visited neighbors of atoms in a subring structure to search for 3 member ring
    //! @param sub_ring the connected ring substructure of a molecule
    //! @param curr_idx, the index of current atom within the sub_ring fragment
    //! @param target_idx, the index of the original hetero atom from which we started the search
    //! @param num_bond_visited, number of bonded having stepped through since start, should be <= 3 for a 3 member ring
    //! @param visited: a list of boolean keeping track of if the atoms in the structure has been visited in the recursive calls or not
    //! @return 1 if the 3rd atom visited is the target atom (the original heteroatom), 0 otherwise

    size_t MoleculeTotalToxicFragments::SearchThreeMemberCycle
    (
      const chemistry::AtomVector< chemistry::AtomComplete> &SUB_RING,
      const chemistry::AtomConformationalInterface &CURR_ATOM,
      const chemistry::AtomConformationalInterface &ORIG_ATOM,
      const size_t NUM_BOND_VISITED
    ) const
    {

      if
      (
          SUB_RING.GetAtomIndex( CURR_ATOM) == SUB_RING.GetAtomIndex( ORIG_ATOM)
          && NUM_BOND_VISITED == 3
      )
      {
        return 1; // find the original heteroatom after stepping through 3 bonds --> a 3 member cycle with a hetero atom found
      }

      //check if more than 3 bonds has been stepped through
      if( NUM_BOND_VISITED > 3)
      {
        return 0;
      }

      // a binary value to indicate if a three member heterocycle has been found
      size_t found( 0);

      //otherwise, continue to go check each neighbor of the current atom
      for
      (
          auto b_itr( CURR_ATOM.GetBonds().Begin()),
          b_itr_end( CURR_ATOM.GetBonds().End());
          b_itr != b_itr_end;
          ++b_itr
      )
      {
        found +=
            SearchThreeMemberCycle
            (
              SUB_RING,
              b_itr->GetTargetAtom(),
              ORIG_ATOM,
              NUM_BOND_VISITED + 1
            );
      }
      return found;
    }

    //! @brief Searches for stabilizing groups nearby mutagenic substructures
    //! @return a vector (corresponding to each stabilizing group) mapping
    //! the bond distance between the stabilizing group and mutagenic group
    //! connecting atoms
    storage::Vector< storage::Map< storage::Pair< size_t, size_t>, size_t> >
    MoleculeTotalToxicFragments::FindStabilizingGroups
    (
      const chemistry::FragmentComplete &MOLECULE
    ) const
    {
      // requires defined stabilizing groups
      if( !m_StabilizingFragments.GetSize())
      {
        return storage::Vector< storage::Map< storage::Pair< size_t, size_t>, size_t> >();
      }

      // initialize graph; this resolution is probably sufficient
      chemistry::ConformationGraphConverter converter
      (
        chemistry::ConformationGraphConverter::AtomComparisonType::e_ElementType,
        chemistry::ConfigurationalBondTypeData::Data::e_BondOrderAmideOrAromaticWithRingness
      );

      // generate graphs of our input molecule
      graph::ConstGraph< size_t, size_t> mol_graph( converter( MOLECULE));
      util::OwnPtr< graph::ConstGraph< size_t, size_t> > mol_graph_ptr( &mol_graph, false);

      // initialize isomorphism and set our primary molecule as graph A
      graph::CommonSubgraphIsomorphism< size_t, size_t> common_subgraph_iso;
      common_subgraph_iso.SetGraphA( mol_graph_ptr);

      // iterate over the stabilizing fragments
      for
      (
          auto frag_itr( m_StabilizingFragments.Begin()), frag_itr_end( m_StabilizingFragments.End());
          frag_itr != frag_itr_end;
          ++frag_itr
      )
      {
        // make graph of the current stabilizing fragment
        graph::ConstGraph< size_t, size_t> frag_graph( converter( *frag_itr));
        util::OwnPtr< graph::ConstGraph< size_t, size_t> > frag_graph_ptr( &frag_graph, false);
        common_subgraph_iso.SetGraphB( frag_graph_ptr);

        // make a new version of the starting molecule
        chemistry::FragmentComplete start_mol( MOLECULE);

        // get the isomorphism
        common_subgraph_iso.FindIsomorphism( common_subgraph_iso.EstimateUpperBounds(), 1);

        graph::Subgraph< size_t, size_t> subgraph
        (
          common_subgraph_iso.GetSubgraphIsomorphismsOfGraphA().FirstElement()
        );

      }
      return storage::Vector< storage::Map< storage::Pair< size_t, size_t>, size_t> >();
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERROR_STREAM the stream to write errors to
    bool MoleculeTotalToxicFragments::ReadInitializerSuccessHook
    (
      const util::ObjectDataLabel &LABEL,
      std::ostream &ERROR_STREAM
    )
    {
      // static initialization check
      if( command::CommandState::IsInStaticInitialization())
      {
        return true;
      }

      // read in stabilizing fragments
      GetMutex().Lock();
      io::IFStream input;
      io::File::MustOpenIFStream
      (
        input,
        !m_StabilizingFragmentsFilename.empty() ?
            m_StabilizingFragmentsFilename :
            chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "stabilizing_groups/stabilizing_groups.sdf.gz"
      );
      m_StabilizingFragments.ReadMoreFromMdl( input, sdf::e_Remove);
      io::File::CloseClearFStream( input);
      GetMutex().Unlock();

      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer MoleculeTotalToxicFragments::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Returns true if the molecule contains a predicted mutagenic substructure; false otherwise");

      parameters.AddInitializer
      (
        "stabilizing_fragments",
        "path to the SDF containing stabilizing fragments to consider",
        io::Serialization::GetAgent( &m_StabilizingFragmentsFilename),
        chemistry::RotamerLibraryFile::GetRotamerFinder().FindFile( "") + "stabilizing_groups/stabilizing_groups.sdf.gz"
      );

      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
