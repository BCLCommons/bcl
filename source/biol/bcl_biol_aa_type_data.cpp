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
#include "biol/bcl_biol_aa_type_data.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_types.h"
#include "biol/bcl_biol_atom_group_types.h"
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_molecular_configuration_shared.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "linal/bcl_linal_matrix.h"
#include "math/bcl_math_sum_function.h"
#include "sdf/bcl_sdf_bond_info.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace biol
  {

    //! @brief PropertyType as string
    //! @param PROPERTY_TYPE the PropertyType
    //! @return the string for the PropertyType
    const std::string &AATypeData::GetPropertyDescriptor( const PropertyType &PROPERTY_TYPE)
    {
      static const std::string s_descriptors[] =
      {
        "AA_NaturalPrevalence",
        "AA_StericalParameter",
        "AA_Polarizability",
        "AA_Volume",
        "AA_Hydrophobicity",
        "AA_IsoelectricPoint",
        "AA_Mass",
        "AA_Charge",
        "AA_pK_EMBOSS",
        "AA_pK_DTASelect",
        "AA_pk_Solomon",
        "AA_pK_Sillero",
        "AA_pK_Rodwell",
        "AA_pK_Patrickios",
        "AA_pK_Wikipedia",
        "AA_pK_Lehninger",
        "AA_pK_Grimsely",
        "AA_pK_Bjellqvist",
        "AA_pK_ProMoST",
        "AA_pK_Bjellqvist_NTerm",
        "AA_pK_Bjellqvist_CTerm",
        "AA_pK_ProMoST_NTerm",
        "AA_pK_ProMoST_CTerm",
        "AA_pK_Carey_NTerm",
        "AA_pK_Carey_CTerm",
        "AA_HelixProbability",
        "AA_StrandProbability",
        "AA_FreeEnergyHelix",
        "AA_FreeEnergyStrand",
        "AA_FreeEnergyCoil",
        "AA_TransferFreeEnergyWhimleyWhite",
        "AA_TransferFreeEnergyEngelmanSeitzGoldman",
        "AA_TransferFreeEnergyKyteDoolittle",
        "AA_TransferFreeEnergyEisenberg",
        "AA_TransferFreeEnergyHoppWoods",
        "AA_TransferFreeEnergyGuy",
        "AA_TransferFreeEnergyJanin",
        "AA_TransferFreeEnergyPuntaMaritan1D",
        "AA_TransferFreeEnergyPuntaMaritan3D",
        "AA_FreeEnergyCore",
        "AA_FreeEnergyTransition",
        "AA_FreeEnergySolution",
        "AA_FreeEnergyCoreHelix",
        "AA_FreeEnergyTransitionHelix",
        "AA_FreeEnergySolutionHelix",
        "AA_FreeEnergyCoreStrand",
        "AA_FreeEnergyTransitionStrand",
        "AA_FreeEnergySolutionStrand",
        "AA_FreeEnergyCoreCoil",
        "AA_FreeEnergyTransitionCoil",
        "AA_FreeEnergySolutionCoil",
        "AA_FreeEnergyCorePore",
        "AA_FreeEnergyCoreMembrane",
        "AA_MembraneStrandOrientationHydrophobicity",
        "AA_FreeEnergyExtracellularBlastBB",
        "AA_FreeEnergyExtracellularTypeBB",
        "AA_SASA",
        "AA_SideChainGirth",
        "AA_SideChainPolarizability",
        "AA_TopologicalPolarSurfaceArea",
        "AA_VdwSurfaceArea",
        "AA_HAcceptors",
        "AA_HDonors",
        "AA_Aromatic",
        GetStaticClassName< PropertyType>()
      };

      return s_descriptors[ PROPERTY_TYPE];
    }

    //! @brief helper function to create a molecular configuration given the component atom types
    //! @param ATOM_TYPES atom types from the AA
    //! @param ONE_LETTER_CODE one letter code; needed to process the exceptions to the general AA connectivity rules
    //! @param THREE_LETTER_CODE three letter code, used for naming the molecule
    //! @param IS_NATURAL all natural AAs must be saturated; so this serves as a check
    //! @param C_TERMINAL whether to include the C-terminal OH group in the molecule
    //! @param N_TERMINAL whether to include the N-terminal H in the molecule
    //! @param HYDROGENATE_TERMINAL_CARBOXYL whether the terminal O is hydrogenated (ignored if !C_TERMINAL)
    //! @param HYDROGENATE_TERMINAL_AMIDE whether the terminal N is hydrogenated (ignored if !N_TERMINAL)
    chemistry::FragmentConfigurationShared CreateFragment
    (
      const storage::Set< AtomType> &ATOM_TYPES,
      const char &ONE_LETTER_CODE,
      const std::string &THREE_LETTER_CODE,
      const bool &IS_NATURAL,
      const bool &C_TERMINAL,
      const bool &N_TERMINAL,
      const bool &HYDROGENATE_TERMINAL_CARBOXYL,
      const bool &HYDROGENATE_TERMINAL_AMIDE
    )
    {
      // create a vector with all component atom types
      storage::Vector< AtomType> atom_type_vec( ATOM_TYPES.Begin(), ATOM_TYPES.End());
      storage::Vector< int> charge_vec( atom_type_vec.GetSize(), 0);
      if( N_TERMINAL && !ATOM_TYPES.IsEmpty())
      {
        atom_type_vec.PushBack( GetAtomTypes().H1);
        charge_vec.PushBack( 0);
        const size_t h_pos( atom_type_vec.Find( GetAtomTypes().H));
        if( h_pos < atom_type_vec.GetSize())
        {
          atom_type_vec.RemoveElements( h_pos, 1);
          charge_vec.RemoveElements( h_pos, 1);
          atom_type_vec.PushBack( GetAtomTypes().H2);
          charge_vec.PushBack( 0);
          if( HYDROGENATE_TERMINAL_AMIDE)
          {
            atom_type_vec.PushBack( GetAtomTypes().H3);
            charge_vec.PushBack( 0);
            charge_vec( atom_type_vec.Find( GetAtomTypes().N)) = 1;
          }
        }
        else if( HYDROGENATE_TERMINAL_AMIDE)
        {
          // terminal proline-like residue
          atom_type_vec.PushBack( GetAtomTypes().H2);
          charge_vec.PushBack( 0);
          charge_vec( atom_type_vec.Find( GetAtomTypes().N)) = 1;
        }
      }
      if( C_TERMINAL && !ATOM_TYPES.IsEmpty())
      {
        atom_type_vec.PushBack( GetAtomTypes().OXT);
        charge_vec.PushBack( 0);
        if( HYDROGENATE_TERMINAL_CARBOXYL)
        {
          atom_type_vec.PushBack( GetAtomTypes().HXT);
          charge_vec.PushBack( 0);
        }
        else
        {
          charge_vec.LastElement() = -1;
        }
      }

      // create a vector with just the IDs from each of those atom types
      const size_t n_atom_types( atom_type_vec.GetSize());

      // create a matrix to hold the bond types
      linal::Matrix< size_t> bond_type_matrix( n_atom_types, n_atom_types, size_t( 0));

      // track number of electrons in bonds for each atom
      storage::Vector< size_t> bond_counts( n_atom_types, size_t( 0));
      storage::Vector< size_t> is_side_chain( n_atom_types, size_t( 0));
      // handle mandatory connections
      for( size_t atom_type_index_a( 0); atom_type_index_a < n_atom_types; ++atom_type_index_a)
      {
        if( atom_type_vec( atom_type_index_a)->IsSideChain())
        {
          is_side_chain( atom_type_index_a) = 1;
        }
        for( size_t atom_type_index_b( atom_type_index_a + 1); atom_type_index_b < n_atom_types; ++atom_type_index_b)
        {
          if( atom_type_vec( atom_type_index_a)->GetConnections().Contains( atom_type_vec( atom_type_index_b)))
          {
            bond_type_matrix( atom_type_index_a, atom_type_index_b) = 1;
            bond_type_matrix( atom_type_index_b, atom_type_index_a) = 1;
            ++bond_counts( atom_type_index_a);
            ++bond_counts( atom_type_index_b);
          }
        }
      }

      // handle the two AAs with non-standard bonds
      if( ONE_LETTER_CODE == 'P')
      {
        // proline, add connection between N and CD
        const size_t atom_type_index_a( atom_type_vec.Find( GetAtomTypes().N));
        const size_t atom_type_index_b( atom_type_vec.Find( GetAtomTypes().CD));
        bond_type_matrix( atom_type_index_a, atom_type_index_b) = 1;
        bond_type_matrix( atom_type_index_b, atom_type_index_a) = 1;
        ++bond_counts( atom_type_index_a);
        ++bond_counts( atom_type_index_b);
      }
      else if( ONE_LETTER_CODE == 'H')
      {
        // histidine, need to set a double bond between CD2 and CG
        const size_t atom_type_index_a( atom_type_vec.Find( GetAtomTypes().CG));
        const size_t atom_type_index_b( atom_type_vec.Find( GetAtomTypes().CD2));
        bond_type_matrix( atom_type_index_a, atom_type_index_b) = 2;
        bond_type_matrix( atom_type_index_b, atom_type_index_a) = 2;
        ++bond_counts( atom_type_index_a);
        ++bond_counts( atom_type_index_b);
      }

      // add double bonds to all necessary atoms
      for( size_t atom_type_index_a( 0); atom_type_index_a < n_atom_types; ++atom_type_index_a)
      {
        // skip atoms that have 4 bonds already
        if( bond_counts( atom_type_index_a) == size_t( 4))
        {
          continue;
        }
        for( size_t atom_type_index_b( atom_type_index_a + 1); atom_type_index_b < n_atom_types; ++atom_type_index_b)
        {
          // skip unbonded atoms
          if( !bond_type_matrix( atom_type_index_a, atom_type_index_b))
          {
            continue;
          }

          // skip atoms that have 4 bonds already
          if( bond_counts( atom_type_index_b) == size_t( 4))
          {
            continue;
          }

          // insert the double bond if it would be made to this atom
          if( atom_type_vec( atom_type_index_a)->GetDoubleBondConnections().Contains( atom_type_vec( atom_type_index_b)))
          {
            bond_type_matrix( atom_type_index_a, atom_type_index_b) = 2;
            bond_type_matrix( atom_type_index_b, atom_type_index_a) = 2;
            ++bond_counts( atom_type_index_a);
            ++bond_counts( atom_type_index_b);
          }
        }
      }

      if( IS_NATURAL)
      {
        // fix charges on side chain atoms
        for( size_t atom_type_index( 0); atom_type_index < n_atom_types; ++atom_type_index)
        {
          // only handle side chain atoms
          if( !is_side_chain( atom_type_index))
          {
            continue;
          }

          // ignore elements that are never ionized in the natural amino acids
          const AtomType &atom_type( atom_type_vec( atom_type_index));
          const chemistry::ElementType &element( atom_type->GetElementType());
          size_t main_group( element->GetMainGroup());

          // handle nitrogen group
          if( main_group == chemistry::GetElementTypes().e_Nitrogen->GetMainGroup())
          {
            if( bond_counts( atom_type_index) == size_t( 4))
            {
              charge_vec( atom_type_index) = 1;
            }
          }
          else if( main_group == chemistry::GetElementTypes().e_Oxygen->GetMainGroup())
          {
            // handle oxygen group (O, S, Se)
            if( bond_counts( atom_type_index) == size_t( 1))
            {
              charge_vec( atom_type_index) = -1;
            }
          }
          // other groups are never ionized in the natural amino acids
        }
      }

      storage::Vector< sdf::BondInfo> bond_info;
      storage::Vector< sdf::AtomInfo> atom_info( n_atom_types);
      std::string biol_atom_types;
      for( size_t atom_type_index_a( 0); atom_type_index_a < n_atom_types; ++atom_type_index_a)
      {
        atom_info( atom_type_index_a)
          = sdf::AtomInfo
            (
              chemistry::GetAtomTypes().GetAtomType
              (
                atom_type_vec( atom_type_index_a)->GetElementType(),
                charge_vec( atom_type_index_a)
              ),
              chemistry::e_NonChiral
            );
        biol_atom_types += atom_type_vec( atom_type_index_a).GetName() + " ";
        for( size_t atom_type_index_b( atom_type_index_a + 1); atom_type_index_b < n_atom_types; ++atom_type_index_b)
        {
          const size_t bond_order( bond_type_matrix( atom_type_index_a, atom_type_index_b));
          if( bond_order)
          {
            bond_info.PushBack
            (
              sdf::BondInfo
              (
                atom_type_index_a,
                atom_type_index_b,
                bond_order == 2
                ? chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond
                : chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond
              )
            );
          }
        }
      }
      for( size_t atom_type_index_a( 0); atom_type_index_a < n_atom_types; ++atom_type_index_a)
      {
        if( bond_counts( atom_type_index_a) == size_t( 0))
        {
          BCL_MessageCrt
          (
            "Atom type: " + atom_type_vec( atom_type_index_a).GetName()
            + " has no connected atoms in " + THREE_LETTER_CODE
          );
        }
        if( atom_type_vec( atom_type_index_a)->GetElementType() == chemistry::GetElementTypes().e_Hydrogen)
        {
          BCL_Assert
          (
            bond_counts( atom_type_index_a) == size_t( 1),
            "H should only have 1 bond, but has " + util::Format()( bond_counts( atom_type_index_a))
            + " in " + THREE_LETTER_CODE
          );
        }
        else
        {
          BCL_Assert
          (
            bond_counts( atom_type_index_a) < size_t( 5),
            "Atom type: " + atom_type_vec( atom_type_index_a).GetName()
            + " has " + util::Format()( bond_counts( atom_type_index_a)) + " connected atoms in " + THREE_LETTER_CODE
          );
        }
      }

      // set chirality of CA and CB

      // set CA chirality
      if( ATOM_TYPES.Contains( GetAtomTypes().CA))
      {
        chemistry::ChiralityEnum ca_chirality( chemistry::e_SChirality);
        if( ONE_LETTER_CODE == 'C' || THREE_LETTER_CODE == "DAL" || THREE_LETTER_CODE == "DPN")
        {
          // Cysteine is the sole natural AA with R chirality
          ca_chirality = chemistry::e_RChirality;
        }
        atom_info( atom_type_vec.Find( GetAtomTypes().CA)).SetChirality( ca_chirality);
      }

      // handle CB chirality
      if( ATOM_TYPES.Contains( GetAtomTypes().CB))
      {
        chemistry::ChiralityEnum cb_chirality( chemistry::e_NonChiral);
        if( ONE_LETTER_CODE == 'T')
        {
          // threonine has R-chirality on its CB atom
          cb_chirality = chemistry::e_RChirality;
        }
        else if( ONE_LETTER_CODE == 'I')
        {
          // isoleucine has S-chirality on its CB atom
          cb_chirality = chemistry::e_SChirality;
        }
        atom_info( atom_type_vec.Find( GetAtomTypes().CB)).SetChirality( cb_chirality);
      }

      // proline has variable N-chirality
      if( ONE_LETTER_CODE == 'P')
      {
        atom_info( atom_type_vec.Find( GetAtomTypes().N)).SetChirality( chemistry::e_UnknownChirality);
      }

      chemistry::AtomVector< chemistry::AtomComplete> atoms( atom_info, bond_info);
      // determine atom and bond types
      chemistry::AtomsCompleteStandardizer( atoms, THREE_LETTER_CODE, false);
      const size_t atoms_old( atoms.GetSize());
      chemistry::HydrogensHandler::SaturatePartial( atoms, is_side_chain);
      const size_t atoms_new( atoms.GetSize());
      BCL_Assert
      (
        atoms_new == atoms_old || !IS_NATURAL,
        util::Format()( atoms_new - atoms_old) + " H should be added to " + THREE_LETTER_CODE
      );

      // create a vector with indices of all side chain atoms, and a separate vector with all non-side-chain-atom ids
      storage::Vector< size_t> side_chain_atom_ids, nonside_chain_atom_ids;
      for( size_t atom_type_index( 0); atom_type_index < n_atom_types; ++atom_type_index)
      {
        if( atom_type_vec( atom_type_index)->IsSideChain())
        {
          side_chain_atom_ids.PushBack( atom_type_index);
        }
        else
        {
          nonside_chain_atom_ids.PushBack( atom_type_index);
        }
      }
      // add HUnk atom type for all missing H
      for( size_t added_h( atoms_old); added_h < atoms_new; ++added_h)
      {
        biol_atom_types += "HUnk ";
        side_chain_atom_ids.PushBack( added_h);
      }
      chemistry::FragmentComplete fragment_complete
      (
        atoms,
        THREE_LETTER_CODE,
        storage::Map< std::string, std::string>::Create
        (
          std::pair< std::string, std::string>( "BiolAtomTypes", biol_atom_types)
        )
      );

      // calculate all conformation & neighbor AA independent descriptors and store them on the molecule
      static storage::Vector< descriptor::CheminfoProperty> s_descriptors
      (
        storage::Vector< descriptor::CheminfoProperty>::Create
        (
          descriptor::GetCheminfoProperties().calc_NAtoms,
          descriptor::GetCheminfoProperties().calc_NStereo,
          descriptor::GetCheminfoProperties().calc_NRotBond,
          descriptor::GetCheminfoProperties().calc_NRings,
          descriptor::GetCheminfoProperties().calc_NAromaticRings,
          descriptor::GetCheminfoProperties().calc_NConjugatedRings,
          descriptor::GetCheminfoProperties().calc_NNonConjugatedRings,
          descriptor::GetCheminfoProperties().calc_HbondAcceptor,
          descriptor::GetCheminfoProperties().calc_HbondDonor,
          descriptor::GetCheminfoProperties().calc_Polarizability,
          descriptor::GetCheminfoProperties().calc_MolWeight,
          descriptor::GetCheminfoProperties().calc_TopologicalPolarSurfaceArea,
          descriptor::GetCheminfoProperties().calc_EstCovalentSurfaceArea,
          descriptor::GetCheminfoProperties().calc_EstVdwSurfaceArea,
          descriptor::GetCheminfoProperties().calc_EstVdwSurfaceAreaCSD,
          descriptor::GetCheminfoProperties().calc_MolTotalFormalCharge
        )
      );
      // create the side chain fragment
      chemistry::AtomVector< chemistry::AtomComplete> sidechain_atoms
      (
        fragment_complete.GetAtomInfo(),
        fragment_complete.GetBondInfo()
      );
      sidechain_atoms.Reorder( side_chain_atom_ids);
      chemistry::FragmentComplete sidechain( sidechain_atoms, "");

      // create the non-side chain atoms
      chemistry::AtomVector< chemistry::AtomComplete> nonsidechain_atoms
      (
        fragment_complete.GetAtomInfo(),
        fragment_complete.GetBondInfo()
      );
      nonsidechain_atoms.Reorder( nonside_chain_atom_ids);
      chemistry::FragmentComplete backbone( nonsidechain_atoms, "");
      for
      (
        storage::Vector< descriptor::CheminfoProperty>::iterator
          itr_props( s_descriptors.Begin()), itr_props_end( s_descriptors.End());
        itr_props != itr_props_end;
        ++itr_props
      )
      {
        // store the property for each fragment complete
        fragment_complete.StoreProperty( itr_props->GetString(), ( *itr_props)->SumOverObject( fragment_complete));
        fragment_complete.StoreProperty( "SideChain" + itr_props->GetString(), ( *itr_props)->SumOverObject( sidechain));
        fragment_complete.StoreProperty( "BackBone" + itr_props->GetString(), ( *itr_props)->SumOverObject( backbone));
      }

      chemistry::FragmentConfigurationShared fragment( fragment_complete);
      if( util::GetMessenger().GetCurrentMessageLevel() == util::Message::e_Debug)
      {
        // write the AA fragment off to the screen
        fragment_complete.WriteMDL( util::GetLogger());
      }
      return fragment;
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AATypeData::s_Instance( GetObjectInstances().AddInstance( new AATypeData()));

    //! @brief construct undefined amino acid
    AATypeData::AATypeData() :
      m_ThreeLetterCode( "UND"),
      m_OneLetterCode( 'U'),
      m_IsNaturalAminoAcid( false),
      m_ParentType( "UND"),
      m_AtomTypes(),
      m_FirstSidechainAtomType()
    {
      // fill properties with undefined
      std::fill_n( m_Properties, size_t( s_NumberPropertyTypes), util::GetUndefined< double>());
    }

    //! @brief construct AATypeData with undefined properties
    //! @param THREE_LETTER_CODE Three letter code for this amino acid type
    //! @param ONE_LETTER_CODE One letter code for this amino acid type
    //! @param IS_NATURAL_AA boolean to indicate whether this amino acid type is one of 20 natural ones
    //! @param PARENT_TYPE name of parent AAType, this is the same as the AAType for natural amino acids
    AATypeData::AATypeData
    (
      const std::string &THREE_LETTER_CODE,
      const char ONE_LETTER_CODE,
      const bool IS_NATURAL_AA,
      const std::string &PARENT_TYPE
    ) :
      m_ThreeLetterCode( THREE_LETTER_CODE),
      m_OneLetterCode( ONE_LETTER_CODE),
      m_IsNaturalAminoAcid( IS_NATURAL_AA),
      m_ParentType( PARENT_TYPE),
      m_AtomTypes()
    {
      // fill properties with undefined
      std::fill_n( m_Properties, size_t( s_NumberPropertyTypes), util::GetUndefined< double>());
    }

    //! @brief construct amino acid from all its data
    AATypeData::AATypeData
    (
      const std::string &THREE_LETTER_CODE,
      const char ONE_LETTER_CODE,
      const bool IS_NATURAL_AA,
      const std::string &PARENT_TYPE,
      const storage::Set< AtomType> &ATOM_TYPES,
      const AtomType &FIRST_SIDECHAIN_ATOM_TYPE,
      const double NATURAL_PREVALENCE,
      const double STERICAL_PARAMETER,
      const double POLARIZABILITY,
      const double VOLUME,
      const double HYDROPHOBICITY,
      const double ISOELECTRIC_POINT,
      const double CHARGE,
      const double PK_EMBOSS,
      const double PK_DTASELECT,
      const double PK_SOLOMON,
      const double PK_SILLERO,
      const double PK_RODWELL,
      const double PK_PATRICKIOS,
      const double PK_WIKIPEDIA,
      const double PK_LEHNINGER,
      const double PK_GRIMSELY,
      const double PK_BJELLQVIST,
      const double PK_PROMOST,
      const double PK_BJELLQVIST_NTERM,
      const double PK_BJELLQVIST_CTERM,
      const double PK_PROMOST_NTERM,
      const double PK_PROMOST_CTERM,
      const double PK_CAREY_NTERM,
      const double PK_CAREY_CTERM,
      const double HELIX_PROBABILITY,
      const double STRAND_PROBABILITY,
      const double FREE_ENERGY_HELIX,
      const double FREE_ENERGY_STRAND,
      const double FREE_ENERGY_COIL,
      const double TRANSFER_FREE_ENERGY_WHIMLEY_WHITE,
      const double TRANSFER_FREE_ENERGY_ENGELMAN_SEITZ_GOLDMAN,
      const double TRANSFER_FREE_ENERGY_KYTE_DOOLITTLE,
      const double TRANSFER_FREE_ENERGY_EISENBERG,
      const double TRANSFER_FREE_ENERGY_HOPP_WOODS,
      const double TRANSFER_FREE_ENERGY_GUY,
      const double TRANSFER_FREE_ENERGY_JANIN,
      const double TRANSFER_FREE_ENERGY_PUNTA_MARITAN_1D,
      const double TRANSFER_FREE_ENERGY_PUNTA_MARITAN_3D,
      const double FREE_ENERGY_CORE,
      const double FREE_ENERGY_TRANSITION,
      const double FREE_ENERGY_SOLUTION,
      const double FREE_ENERGY_CORE_HELIX,
      const double FREE_ENERGY_TRANSITION_HELIX,
      const double FREE_ENERGY_SOLUTION_HELIX,
      const double FREE_ENERGY_CORE_STRAND,
      const double FREE_ENERGY_TRANSITION_STRAND,
      const double FREE_ENERGY_SOLUTION_STRAND,
      const double FREE_ENERGY_CORE_COIL,
      const double FREE_ENERGY_TRANSITION_COIL,
      const double FREE_ENERGY_SOLUTION_COIL,
      const double FREE_ENERGY_CORE_PORE,
      const double FREE_ENERGY_CORE_MEMBRANE,
      const double BETA_BARREL_HYDROPHOBICITY,
      const double FREE_ENERGY_EXT_BLAST_BB,
      const double FREE_ENERGY_EXT_TYPE_BB,
      const double SASA,
      const double SIDE_CHAIN_GIRTH,
      const double SIDE_CHAIN_POLARIZABILITY,
      const double SIDE_CHAIN_TOPOLOGICAL_POLAR_SURFACE_AREA,
      const double SIDE_CHAIN_VAN_DER_WAALS_SURFACE_AREA,
      const int    HBOND_ACCEPTORS,
      const int    HBOND_DONORS,
      const bool   IS_AROMATIC,
      const double VDW_RADIUS_CB
    ) :
      m_ThreeLetterCode( THREE_LETTER_CODE),
      m_OneLetterCode( ONE_LETTER_CODE),
      m_IsNaturalAminoAcid( IS_NATURAL_AA),
      m_ParentType( PARENT_TYPE),
      m_AtomTypes( ATOM_TYPES),
      m_FirstSidechainAtomType( FIRST_SIDECHAIN_ATOM_TYPE),
      m_DihedralAtoms()
    {
      m_Properties[ e_NaturalPrevalence]                      = NATURAL_PREVALENCE;
      m_Properties[ e_StericalParameter]                      = STERICAL_PARAMETER;
      m_Properties[ e_Polarizability]                         = POLARIZABILITY;
      m_Properties[ e_Volume]                                 = VOLUME;
      m_Properties[ e_Hydrophobicity]                         = HYDROPHOBICITY;
      m_Properties[ e_IsoelectricPoint]                       = ISOELECTRIC_POINT;
      m_Properties[ e_Mass]                                   = CalculateMass();
      m_Properties[ e_Charge]                                 = CHARGE;
      m_Properties[ e_pK_EMBOSS]                              = PK_EMBOSS;
      m_Properties[ e_pK_DTASelect]                           = PK_DTASELECT;
      m_Properties[ e_pk_Solomon]                             = PK_SOLOMON;
      m_Properties[ e_pK_Sillero]                             = PK_SILLERO;
      m_Properties[ e_pK_Rodwell]                             = PK_RODWELL;
      m_Properties[ e_pK_Patrickios]                          = PK_PATRICKIOS;
      m_Properties[ e_pK_Wikipedia]                           = PK_WIKIPEDIA;
      m_Properties[ e_pK_Lehninger]                           = PK_LEHNINGER;
      m_Properties[ e_pK_Grimsely]                            = PK_GRIMSELY;
      m_Properties[ e_pK_Bjellqvist]                          = PK_BJELLQVIST;
      m_Properties[ e_pK_ProMoST]                             = PK_PROMOST;
      m_Properties[ e_pK_Bjellqvist_NTerm]                    = PK_BJELLQVIST_NTERM;
      m_Properties[ e_pK_Bjellqvist_CTerm]                    = PK_BJELLQVIST_CTERM;
      m_Properties[ e_pK_ProMoST_NTerm]                       = PK_PROMOST_NTERM;
      m_Properties[ e_pK_ProMoST_CTerm]                       = PK_PROMOST_CTERM;
      m_Properties[ e_pK_Carey_NTerm]                         = PK_CAREY_NTERM;
      m_Properties[ e_pK_Carey_CTerm]                         = PK_CAREY_CTERM;
      m_Properties[ e_HelixProbability]                       = HELIX_PROBABILITY;
      m_Properties[ e_StrandProbability]                      = STRAND_PROBABILITY;
      m_Properties[ e_FreeEnergyHelix]                        = FREE_ENERGY_HELIX;
      m_Properties[ e_FreeEnergyStrand]                       = FREE_ENERGY_STRAND;
      m_Properties[ e_FreeEnergyCoil]                         = FREE_ENERGY_COIL;
      m_Properties[ e_TransferFreeEnergyWhimleyWhite]         = TRANSFER_FREE_ENERGY_WHIMLEY_WHITE;
      m_Properties[ e_TransferFreeEnergyEngelmanSeitzGoldman] = TRANSFER_FREE_ENERGY_ENGELMAN_SEITZ_GOLDMAN;
      m_Properties[ e_TransferFreeEnergyKyteDoolittle]        = TRANSFER_FREE_ENERGY_KYTE_DOOLITTLE;
      m_Properties[ e_TransferFreeEnergyEisenberg]            = TRANSFER_FREE_ENERGY_EISENBERG;
      m_Properties[ e_TransferFreeEnergyHoppWoods]            = TRANSFER_FREE_ENERGY_HOPP_WOODS;
      m_Properties[ e_TransferFreeEnergyGuy]                  = TRANSFER_FREE_ENERGY_GUY;
      m_Properties[ e_TransferFreeEnergyJanin]                = TRANSFER_FREE_ENERGY_JANIN;
      m_Properties[ e_TransferFreeEnergyPuntaMaritan1D]       = TRANSFER_FREE_ENERGY_PUNTA_MARITAN_1D;
      m_Properties[ e_TransferFreeEnergyPuntaMaritan3D]       = TRANSFER_FREE_ENERGY_PUNTA_MARITAN_3D;
      m_Properties[ e_FreeEnergyCore]                         = FREE_ENERGY_CORE;
      m_Properties[ e_FreeEnergyTransition]                   = FREE_ENERGY_TRANSITION;
      m_Properties[ e_FreeEnergySolution]                     = FREE_ENERGY_SOLUTION;
      m_Properties[ e_FreeEnergyCoreHelix]                    = FREE_ENERGY_CORE_HELIX;
      m_Properties[ e_FreeEnergyTransitionHelix]              = FREE_ENERGY_TRANSITION_HELIX;
      m_Properties[ e_FreeEnergySolutionHelix]                = FREE_ENERGY_SOLUTION_HELIX;
      m_Properties[ e_FreeEnergyCoreStrand]                   = FREE_ENERGY_CORE_STRAND;
      m_Properties[ e_FreeEnergyTransitionStrand]             = FREE_ENERGY_TRANSITION_STRAND;
      m_Properties[ e_FreeEnergySolutionStrand]               = FREE_ENERGY_SOLUTION_STRAND;
      m_Properties[ e_FreeEnergyCoreCoil]                     = FREE_ENERGY_CORE_COIL;
      m_Properties[ e_FreeEnergyTransitionCoil]               = FREE_ENERGY_TRANSITION_COIL;
      m_Properties[ e_FreeEnergySolutionCoil]                 = FREE_ENERGY_SOLUTION_COIL;
      m_Properties[ e_FreeEnergyCoreFacingPore]               = FREE_ENERGY_CORE_PORE;
      m_Properties[ e_FreeEnergyCoreFacingMembrane]           = FREE_ENERGY_CORE_MEMBRANE;
      m_Properties[ e_MembraneStrandOrientationHydrophobicity]= BETA_BARREL_HYDROPHOBICITY;
      m_Properties[ e_FreeEnergyExtracellularBlastBB]         = FREE_ENERGY_EXT_BLAST_BB;
      m_Properties[ e_FreeEnergyExtracellularTypeBB]          = FREE_ENERGY_EXT_TYPE_BB;
      m_Properties[ e_SASA]                                   = SASA;
      m_Properties[ e_SideChainGirth]                         = SIDE_CHAIN_GIRTH;
      m_Properties[ e_SideChainPolarizability]                = SIDE_CHAIN_POLARIZABILITY;
      m_Properties[ e_SideChainTopologicalPolarSurfaceArea]   = SIDE_CHAIN_TOPOLOGICAL_POLAR_SURFACE_AREA;
      m_Properties[ e_SideChainVanDerWaalsSurfaceArea]        = SIDE_CHAIN_VAN_DER_WAALS_SURFACE_AREA;
      m_Properties[ e_SideChainHBondAcceptors]                = HBOND_ACCEPTORS;
      m_Properties[ e_SideChainHBondDonors]                   = HBOND_DONORS;
      m_Properties[ e_Aromatic]                               = IS_AROMATIC ? 1 : 0;

      m_ChemistryTypes.Resize( GetAtomTypes().GetEnumCount());
      m_VdwRadii.Resize( m_ChemistryTypes.GetSize(), util::GetUndefinedDouble());
      // get both a fully-terminated fragment and a fully integrated fragment. Together, they have the
      // full set of atom types
      const chemistry::FragmentConfigurationShared &terminal_fragment( GetFragment( true, true));
      const chemistry::FragmentConfigurationShared &std_fragment( GetFragment( false, false));
      const storage::Vector< std::string> term_biol_types
      (
        util::SplitString( terminal_fragment.GetMDLProperty( "BiolAtomTypes"))
      );
      auto itr_biol_atom_types_str( term_biol_types.Begin());
      for( auto itr_atom( terminal_fragment.GetAtomsIterator()); itr_atom.NotAtEnd(); ++itr_atom, ++itr_biol_atom_types_str)
      {
        if( !GetAtomTypes().HaveEnumWithName( *itr_biol_atom_types_str))
        {
          continue;
        }
        AtomType atom_type( *itr_biol_atom_types_str);
        m_ChemistryTypes( atom_type.GetIndex()) = itr_atom->GetAtomType();
      }
      const storage::Vector< std::string> std_biol_types
      ( 
        util::SplitString( std_fragment.GetMDLProperty( "BiolAtomTypes"))
      );
      itr_biol_atom_types_str = std_biol_types.Begin();
      for( auto itr_atom( std_fragment.GetAtomsIterator()); itr_atom.NotAtEnd(); ++itr_atom, ++itr_biol_atom_types_str)
      {
        if( !GetAtomTypes().HaveEnumWithName( *itr_biol_atom_types_str))
        {
          continue;
        }
        AtomType atom_type( *itr_biol_atom_types_str);
        m_ChemistryTypes( atom_type.GetIndex()) = itr_atom->GetAtomType();
        m_VdwRadii( atom_type.GetIndex())
          = itr_atom->GetAtomType()->GetAtomTypeProperty( chemistry::AtomTypeData::e_VdWaalsRadiusCSD);
      }
      // backbone vdw radii; computed ignoring residues on the same SSE (except for coils) and only considering other
      // backbone and Cb atoms, also ignoring backbone O h-bonded to glycine H
      m_VdwRadii( GetAtomTypes().C)   = 1.62;
      m_VdwRadii( GetAtomTypes().CA)  = 1.68;
      m_VdwRadii( GetAtomTypes().N)   = 1.55;
      m_VdwRadii( GetAtomTypes().O)   = 0.80;
      if( VDW_RADIUS_CB)
      {
        m_VdwRadii( m_FirstSidechainAtomType) = VDW_RADIUS_CB;
      }
    }

    //! @brief virtual copy constructor
    AATypeData *AATypeData::Clone() const
    {
      return new AATypeData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AATypeData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the parent amino acid type
    //! this maps uncommon aas to their natural counterpart like Seleno-methionin to methionin
    //! @return AAType of parent for this uncommon amino acid
    const AAType &AATypeData::GetParentType() const
    {
      return *m_ParentTypePtr;
    }

    //! @brief atom types for that amino acid type without hydrogens
    //! @return a set of atom types that are part of that amino acid type without hydrogens
    storage::Set< AtomType> AATypeData::GetAllowedHeavyAtomTypes() const
    {
      // initialize set
      storage::Set< AtomType> atom_types;

      // iterate over atom types
      for
      (
        storage::Set< AtomType>::const_iterator itr( m_AtomTypes.Begin()), itr_end( m_AtomTypes.End());
        itr != itr_end; ++itr
      )
      {
        // add to set if not hydrogen
        if( ( *itr)->GetElementType() != chemistry::GetElementTypes().e_Hydrogen)
        {
          atom_types.Insert( *itr);
        }
      }

      return atom_types;
    }

    //! @brief return the pdb atom name for the given AtomType
    //! @param ATOM_TYPE atom type
    //! @return the pdb atom name with proper spacing, so that it can be written to the pdb file, empty if this amino acid type does not have this atom type
    const std::string &AATypeData::GetPDBAtomName( const AtomType &ATOM_TYPE) const
    {
      // check if the atom type is part of the amino acid type
      if( !m_AtomTypes.Contains( ATOM_TYPE))
      {
        // check if it is a terminal type
        const storage::Map< AtomType, std::string>::const_iterator term_itr
        (
          GetAtomTypes().GetTerminalExtraAtomTypes().Find( ATOM_TYPE)
        );

        // could be found, return the pdb name
        if( term_itr != GetAtomTypes().GetTerminalExtraAtomTypes().End())
        {
          return term_itr->second;
        }

        static const std::string s_undefined;
        return s_undefined;
      }

      // find entry for that atom type in map of mapped types
      const storage::Map< AtomType, std::string>::const_iterator find_itr( m_PDBAtomName.Find( ATOM_TYPE));

      // special pdb atom name defined
      if( find_itr != m_PDBAtomName.End())
      {
        return find_itr->second;
      }

      // end
      return ATOM_TYPE->AtomTypeData::GetName();
    }

    //! @brief return the associated atom type from a given atom name (IUPAC or PDB)
    //! @param ATOM_NAME the atom name defined by either the pdb or IUPAC for this aa type
    //! @return the atom type for that pdb atom name, if it is defined for that amino acid type
    AtomType AATypeData::GetAtomTypeFromAtomName( const std::string &ATOM_NAME) const
    {
      // check if AtomType can be constructed from PDB_ATOM_NAME
      const AtomType type( GetAtomTypes().TypeFromPDBAtomName( ATOM_NAME));

      // valid IUPAC name
      if( type.IsDefined())
      {
        // is pdb atom name a valid atom for this type or is atom type a terminal atom type
        if
        (
             ( m_AtomTypes.Find( type) != m_AtomTypes.End())
          || GetAtomTypes().GetTerminalExtraAtomTypes().Has( type)
        )
        {
          return type;
        }

        // wrong atom type for this aa type
        return GetAtomTypes().e_Undefined;
      }

      // iterate over pdb atom map
      for
      (
        storage::Map< AtomType, std::string>::const_iterator itr( m_PDBAtomName.Begin()), itr_end( m_PDBAtomName.End());
        itr != itr_end;
        ++itr
      )
      {
        if( itr->second.find( ATOM_NAME) != std::string::npos)
        {
          return itr->first;
        }
      }

      // could still be a pdb terminal atom type name
      const AtomType term_type( GetAtomTypes().GetTerminalAtomTypeFromName( ATOM_NAME));

      // will be defined, if it is a pdb terminal atom type, or undefined otherwise
      return term_type;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief map a given atom type to a pdb atom name - is only used in the constructor of the AtomTypes class
    //! @param ATOM_TYPE the atom type
    //! @param PDB_ATOM_NAME the atom name defined by the pdb for this aa type
    //! @return true if mapping was successful (false if the atom type was already added)
    bool AATypeData::AddAtomTypePDBAtomNameMapping( const AtomType &ATOM_TYPE, const std::string &PDB_ATOM_NAME)
    {
      return m_PDBAtomName.Insert( std::pair< AtomType, std::string>( ATOM_TYPE, PDB_ATOM_NAME)).second;
    }

    //! @brief check if given ATOM_TYPE is part of the amino acid
    //! @param ATOM_TYPE type that is checked for
    //! @return true is this amino acid type has the given ATOM_TYPE
    bool AATypeData::DoesContainAtomType( const AtomType &ATOM_TYPE) const
    {
      return m_AtomTypes.Find( ATOM_TYPE) != m_AtomTypes.End();
    }

    //! @brief get properties to be used in ANN as a vector
    //! @return properties to be used in ANN as a vector
    linal::Vector< double> AATypeData::GetPropertiesForANN() const
    {
      // initialize vector of properties
      linal::Vector< double> properties_for_ANN( s_NumberPropertyTypesForANN);

      // initialize ptr to first of property vector
      double *ptr = properties_for_ANN.Begin();

      //Get 7 AA properties
      *ptr++ = m_Properties[ e_StericalParameter]; // sterical parameter
      *ptr++ = m_Properties[ e_Polarizability];    // polarizability
      *ptr++ = m_Properties[ e_Volume];            // volume
      *ptr++ = m_Properties[ e_Hydrophobicity];    // hydrophobicity
      *ptr++ = m_Properties[ e_IsoelectricPoint];  // isoelectric point
      *ptr++ = m_Properties[ e_HelixProbability];  // helix probability
      *ptr++ = m_Properties[ e_StrandProbability]; // strand probability

      // return
      return properties_for_ANN;
    }

    //! @brief get a chemistry::FragmentConfigurationalShared representation of the AA
    //! @param C_TERMINAL whether to return a representation with the C-terminal carboxyl group
    //! @param N_TERMINAL whether to return a representation with the N-terminal amino group
    //! @return a chemistry::FragmentConfigurationalShared representation of the AA
    const chemistry::FragmentConfigurationShared &AATypeData::GetFragment
    (
      const bool &C_TERMINAL,
      const bool &N_TERMINAL
    ) const
    {
      // create static representations of this AA as a fragment
      // 1. The isolated AA (fully saturated, with terminal amino/carboxyl groups)
      // 2. The AA at the N-terminus of the molecule
      // 3. The AA at the C-terminus of the molecule
      // 4. The AA with both sides having valences for a continued protein chain
      static storage::Map< std::string, storage::VectorND< 4, chemistry::FragmentConfigurationShared> > s_map;
      storage::Map< std::string, storage::VectorND< 4, chemistry::FragmentConfigurationShared> >::const_iterator
        itr( s_map.Find( m_ThreeLetterCode));
      static float s_default_ph( 7.0);
      if( itr == s_map.End())
      {
        // choose which atoms to include based on the pH
        storage::Set< AtomType> atom_types_at_ph( m_AtomTypes);

        // check for a side chain H that is ionized at the default ph
        if( m_SideChainIonizableHydrogen.IsDefined() && s_default_ph > m_Properties[ e_pK_Rodwell])
        {
          BCL_MessageVrb( "Ionizing " + m_SideChainIonizableHydrogen.GetName() + " from default " + m_ThreeLetterCode);
          atom_types_at_ph.Erase( m_SideChainIonizableHydrogen);
        }

        // determine whether to hydrogenate the terminal carboxyl
        const bool hydrogenate_term_c( s_default_ph < m_Properties[ e_pK_Carey_CTerm]);
        const bool hydrogenate_term_n( s_default_ph < m_Properties[ e_pK_Carey_NTerm]);

        // need to create the fragments
        s_map[ m_ThreeLetterCode] =
          storage::VectorND< 4, chemistry::FragmentConfigurationShared>
          (
            CreateFragment( atom_types_at_ph, m_OneLetterCode, m_ThreeLetterCode, m_IsNaturalAminoAcid, true,  true,  hydrogenate_term_c, hydrogenate_term_n),
            CreateFragment( atom_types_at_ph, m_OneLetterCode, m_ThreeLetterCode, m_IsNaturalAminoAcid, true,  false, hydrogenate_term_c, hydrogenate_term_n),
            CreateFragment( atom_types_at_ph, m_OneLetterCode, m_ThreeLetterCode, m_IsNaturalAminoAcid, false, true,  hydrogenate_term_c, hydrogenate_term_n),
            CreateFragment( atom_types_at_ph, m_OneLetterCode, m_ThreeLetterCode, m_IsNaturalAminoAcid, false, false, hydrogenate_term_c, hydrogenate_term_n)
          );
        itr = s_map.Find( m_ThreeLetterCode);
      }
      return itr->second( C_TERMINAL ? ( N_TERMINAL ? 0 : 1) : ( N_TERMINAL ? 2 : 3));
    }

    //! @brief Get the chemistry atom type of a particular atom in this aa type
    //! @param ATOM the biol::AtomType of interest
    //! @return chemistry::AtomType of that Atom
    chemistry::AtomType AATypeData::GetChemistryAtomType( const AtomType &ATOM_TYPE) const
    {
      return ATOM_TYPE.IsDefined() ? m_ChemistryTypes( ATOM_TYPE.GetIndex()) : chemistry::AtomType();
    }

    //! @brief Get the effective VdW radius of an atom type for this AA to other AAs
    //! @param ATOM the biol::AtomType of interest
    //! @return van-der-waals radius to external AAs of the given atom type
    double AATypeData::GetVdwRadiusToOtherAA( const AtomType &ATOM_TYPE) const
    {
      return ATOM_TYPE.IsDefined() ? m_VdwRadii( ATOM_TYPE.GetIndex()) : 0.0;
    }

    //! @brief gets the structure factor for this AA Type
    //! @return the structure factor for this AA Type
    util::ShPtr< math::FunctionInterfaceSerializable< restraint::SasDataParameters, double> > AATypeData::GetStructureFactor() const
    {
      // initialize structure factors
      util::ShPtr< math::SumFunction< restraint::SasDataParameters, double> > structure_factors
      (
        new math::SumFunction< restraint::SasDataParameters, double>()
      );

      // iterate over the atoms of a given amino acid
      for
      (
        storage::Set< AtomType>::const_iterator atom_itr( m_AtomTypes.Begin()), atom_itr_end( m_AtomTypes.End());
        atom_itr != atom_itr_end; ++atom_itr
      )
      {
        // Do not add hydrogen to the structure_factor function
        if( ( *atom_itr)->GetElementType() != chemistry::GetElementTypes().e_Hydrogen)
        {
          // add the structure factor
          *structure_factors += *( GetAtomGroupTypes().GetType( GetAATypes().AATypeFromThreeLetterCode( m_ThreeLetterCode), *atom_itr));
        }
      }

      // end
      return structure_factors;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AATypeData::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_ThreeLetterCode       , ISTREAM);
      io::Serialize::Read( m_OneLetterCode         , ISTREAM);
      io::Serialize::Read( m_IsNaturalAminoAcid    , ISTREAM);
      io::Serialize::Read( m_ParentType            , ISTREAM);
      io::Serialize::Read( m_AtomTypes             , ISTREAM);
      io::Serialize::Read( m_PDBAtomName           , ISTREAM);
      io::Serialize::Read( m_FirstSidechainAtomType, ISTREAM);
      io::Serialize::Read( m_DihedralAtoms         , ISTREAM);

      BCL_MessageCrt( "AAProperty reading not implemented yet");

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream
    //! @return ostream which was written to
    std::ostream &AATypeData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ThreeLetterCode       , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_OneLetterCode         , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_IsNaturalAminoAcid    , OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_ParentType            , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AtomTypes             , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_PDBAtomName           , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FirstSidechainAtomType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DihedralAtoms         , OSTREAM, INDENT);

      BCL_MessageCrt( "AAProperty writing not implemented yet");

      // return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set the dihedral angles for this AA type
    //! @param ANGLES the dihedral angles for this AAType
    void AATypeData::SetDihedralAngles( const storage::Map< ChiAngle::ChiEnum, storage::VectorND< 4, AtomType> > &ANGLES)
    {
      m_DihedralAtoms = ANGLES;
    }

    //! @brief set the ionizable H
    //! @param ATOM_TYPE the ionizable hydrogen
    void AATypeData::SetSideChainIonizableHType( const AtomType &ATOM_TYPE)
    {
      BCL_Assert
      (
        ATOM_TYPE->GetElementType() == chemistry::GetElementTypes().e_Hydrogen,
        "Cannot set pka of non-hydrogen atom type"
      );
      m_SideChainIonizableHydrogen = ATOM_TYPE;
    }

    //! @brief get the counts of what type and how many atoms in a given amino acid
    //! @return Return the map of type and count of atoms in a given amino acid
    storage::Map< chemistry::ElementType, size_t> AATypeData::GetElementTypeCount() const
    {
      // initialize map to hold element type counts
      storage::Map< chemistry::ElementType, size_t> map;

      // iterate over the atom types
      for
      (
        storage::Set< AtomType>::const_iterator itr( m_AtomTypes.Begin()), itr_end( m_AtomTypes.End());
        itr != itr_end; ++itr
      )
      {
        // increment the counter if the map already has the element type, otherwise set it to one
        // also ensure that the atom is not part of the backbone
        const chemistry::ElementType &element_type( ( *itr)->GetElementType());

        // if this element type has been seen
        if( map.Has( element_type))
        {
          // increment the counter
          map[ element_type]++;
        }
        // first time the element has been seen
        else
        {
          // set the counter to 1
          map[ element_type] = 1;
        }
      }

      // end
      return map;
    }

    //! @brief Calculate the mass of this amino acid (not when peptide bonded!)
    //! @return the mass of the amino acid as the sum of all atom masses; non-monoisotopic
    double AATypeData::CalculateMass() const
    {
      double mass( 0.0);
      for
      (
        storage::Set< AtomType>::const_iterator itr( m_AtomTypes.Begin()), itr_end( m_AtomTypes.End());
        itr != itr_end;
        ++itr
      )
      {
        mass += ( *itr)->GetElementType()->GetProperty( chemistry::ElementTypeData::e_Mass);
      }

      return mass;
    }

    //! @brief determine the parent type of amino acid
    void AATypeData::DetermineParentType( const AATypes &AATYPES)
    {
      for
      (
        AATypes::const_iterator itr( AATYPES.Begin()), itr_end( AATYPES.End());
        itr != itr_end;
        ++itr
      )
      {
        if( m_ParentType == ( *itr)->GetThreeLetterCode())
        {
          m_ParentTypePtr = util::ToSiPtr( *itr);
          break;
        }
      }
      if( !m_ParentTypePtr.IsDefined())
      {
        m_ParentTypePtr = util::ToSiPtr( AATYPES.XXX);
      }
    }

  } // namespace biol
} // namespace bcl
