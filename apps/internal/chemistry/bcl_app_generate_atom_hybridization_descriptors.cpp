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
#include "bcl_app_generate_atom_hybridization_descriptors.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_complete.h"
#include "chemistry/bcl_chemistry_atom_vector.h"
#include "chemistry/bcl_chemistry_bond_lengths.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_fragment_split_rigid.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_flag_static.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "graph/bcl_graph_connectivity.h"
#include "io/bcl_io_file.h"
#include "iterate/bcl_iterate_generic.h"
#include "math/bcl_math_linear_least_squares.h"
#include "model/bcl_model_retrieve_dataset_subset.h"
#include "sdf/bcl_sdf_bond_info.h"
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_template_instantiations.h"
#include "storage/bcl_storage_triplet.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    const ApplicationType GenerateAtomHybridizationDescriptors::GenerateAtomHybridizationDescriptors_Instance
    (
      GetAppGroups().AddAppToGroup( new GenerateAtomHybridizationDescriptors(), GetAppGroups().e_InternalChem)
    );

  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief standard constructor
    GenerateAtomHybridizationDescriptors::GenerateAtomHybridizationDescriptors() :
      m_OutputFilenameBase
      (
        new command::FlagStatic
        (
          "output",
          "base name for output of histograms",
          command::Parameter( "output", "base name for output of histograms")
        )
      ),
      m_WriteBondLengthsInitializer
      (
        new command::FlagStatic
        (
          "write_bond_lengths_initializer",
          "Use bond lengths to write code that can be copy/pasted into the initializer for chemistry::BondLengths"
        )
      ),
      m_WriteDescriptors
      (
        new command::FlagStatic
        (
          "write_hybridization_descriptors",
          "Write hybridization descriptors and bin files"
        )
      ),
      m_WriteElementBondStatistics
      (
        new command::FlagStatic
        (
          "write_element_bond_statistics",
          "count the normalized number of element_type-element_type-bond_type occurrences"
        )
      ),
      m_SkipAromaticRingDescriptors
      (
        new command::FlagStatic
        (
          "skip_aromatic_rings",
          "Whether to skip atoms in aromatic rings when calculating descriptors"
        )
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new GenerateAtomHybridizationDescriptors
    GenerateAtomHybridizationDescriptors *GenerateAtomHybridizationDescriptors::Clone() const
    {
      return new GenerateAtomHybridizationDescriptors( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &GenerateAtomHybridizationDescriptors::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> GenerateAtomHybridizationDescriptors::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

    ////////////////////
    // common options //
    ////////////////////

      // add flags for input
      chemistry::FragmentFeed::AddFlags( *sp_cmd);

      // skip aromatic rings
      sp_cmd->AddFlag( m_SkipAromaticRingDescriptors);

      // Output filename base
      sp_cmd->AddFlag( m_OutputFilenameBase);

      // whether to write initializer for chemistry::BondLengths
      sp_cmd->AddFlag( m_WriteBondLengthsInitializer);

      // whether to collect element-element-bond statistics
      sp_cmd->AddFlag( m_WriteElementBondStatistics);

      // whether to write descriptors for hybridization
      sp_cmd->AddFlag( m_WriteDescriptors);

    ///////////////////
    // default flags //
    ///////////////////

      // default flags are unnecessary for this application, but message level and read-me are useful
      command::GetAppDefaultFlags().AddRequiredCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief returns a simple id for a bond type given as a size_t into a string
    //! @returns one of XXXXXX, Aromatic, Single, Double, or Triple depending on the bond type
    size_t GetBasicBondTypeId( const size_t &BOND_TYPE_ID)
    {
      static storage::Vector< size_t> simple_bond_types
      (
        chemistry::GetConfigurationalBondTypes().GetEnumCount(), util::GetUndefined< size_t>()
      );
      if( BOND_TYPE_ID == chemistry::GetConfigurationalBondTypes().e_Undefined.GetIndex())
      {
        return 0;
      }
      else if( !util::IsDefined( simple_bond_types( BOND_TYPE_ID)))
      {
        chemistry::ConfigurationalBondType bond_type( BOND_TYPE_ID);
        simple_bond_types( BOND_TYPE_ID) = bond_type->GetBondData( chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromatic);
      }
      return simple_bond_types( BOND_TYPE_ID);
    }

    namespace
    {
      //! @brief test two (assumed sorted) vectors for any overlapping indices
      bool HaveOverlap( const storage::Vector< size_t> &COMPONENT_A, const storage::Vector< size_t> &COMPONENT_B)
      {
        storage::Vector< size_t>::const_iterator itr_b( COMPONENT_B.Begin()), itr_b_end( COMPONENT_B.End());
        for
        (
          storage::Vector< size_t>::const_iterator itr_a( COMPONENT_A.Begin()), itr_a_end( COMPONENT_A.End());
          itr_a != itr_a_end && itr_b != itr_b_end;
          ++itr_a
        )
        {
          while( itr_b != itr_b_end)
          {
            if( *itr_a > *itr_b)
            {
              ++itr_b;
            }
            else if( *itr_a == *itr_b)
            {
              return true;
            }
            else
            {
              break;
            }
          }
        }
        return false;
      }
    }

    //! @brief generate the descriptors and bond length info as requested
    //! @param FRAG the fragment for which to calculate the descriptors and bond lengths info
    void GenerateAtomHybridizationDescriptors::Process( const chemistry::FragmentComplete &FRAG) const
    {
      static const size_t s_n_descriptors( 55);
      static const size_t s_n_results( 1);
      storage::VectorND< 2, linal::Vector< double> > result_undef
      (
        linal::Vector< double>( s_n_descriptors, double( 0)),
        linal::Vector< double>( s_n_results, double( 0))
      );
      const bool skip_aromatic_rings( m_SkipAromaticRingDescriptors->GetFlag());

      linal::Vector< size_t> n_bonds(      FRAG.GetNumberAtoms(),            size_t( 2));
      linal::Vector< size_t> n_e_in_bonds( FRAG.GetNumberAtoms(),            size_t( 0));
      linal::Vector< size_t> n_empty_orbs( FRAG.GetNumberAtoms(),            size_t( 0));
      linal::Vector< size_t> n_lone_pairs( FRAG.GetNumberAtoms(),            size_t( 0));
      linal::Vector< size_t> is_in_aromatic_ring( FRAG.GetNumberAtoms(),     size_t( 0));
      linal::Vector< size_t> is_in_non_aromatic_ring( FRAG.GetNumberAtoms(), size_t( 0));
      linal::Vector< size_t> is_in_any_ring( FRAG.GetNumberAtoms(),          size_t( 0));
      linal::Vector< size_t> is_in_chain( FRAG.GetNumberAtoms(),             size_t( 0));
      linal::Vector< size_t> max_e_in_bonds( FRAG.GetNumberAtoms(),          size_t( 0));
      linal::Vector< size_t> min_e_in_bonds( FRAG.GetNumberAtoms(),          size_t( 0));
      linal::Vector< double> atom_radius( FRAG.GetNumberAtoms(),             double( 0));
      linal::Vector< size_t> ring_size( FRAG.GetNumberAtoms(),               size_t( 0));
      linal::Vector< double> electronegativity( FRAG.GetNumberAtoms(),       double( 0));
      int molecule_charge( descriptor::GetCheminfoProperties().calc_MolTotalFormalCharge->SumOverObject( FRAG)( 0));

      storage::Vector< chemistry::AtomType> atom_types( FRAG.GetNumberAtoms());
      graph::ConstGraph< size_t, size_t> graph( chemistry::ConformationGraphConverter()( FRAG));
      size_t count( 0);
      for( iterate::Generic< const chemistry::AtomConformationalInterface> itr( FRAG.GetAtomsIterator()); itr.NotAtEnd(); ++itr, ++count)
      {
        n_bonds( count) = itr->GetAtomType()->GetNumberBonds();
        n_e_in_bonds( count) = itr->GetAtomType()->GetNumberElectronsInBonds();
        if( !itr->GetAtomType()->IsGasteigerAtomType())
        {
          return;
        }

        n_lone_pairs( count) = itr->GetAtomType()->GetNumberHybridLonePairs() + itr->GetAtomType()->GetNumberUnhybridizedLonePairs();
        if( itr->GetElementType()->GetPeriod() > size_t( 1))
        {
          if( itr->GetElementType()->GetMainGroup() > size_t( 2))
          {
            n_empty_orbs( count) = std::max( 4 - int( n_lone_pairs( count)) - int( n_e_in_bonds( count)), 0);
          }
          else
          {
            n_empty_orbs( count) = std::max( 2 - int( n_e_in_bonds( count)), 0);
          }
        }

        atom_types( count) = itr->GetAtomType();

        const size_t aromatic_bond_count
        (
          itr->CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsAromatic, size_t( 1))
        );
        const size_t ring_bond_count
        (
          itr->CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsInRing, size_t( 1))
        );
        is_in_aromatic_ring( count) = ( aromatic_bond_count > 0 ? aromatic_bond_count - 1 : 0);
        is_in_any_ring( count) = ( ring_bond_count > 0 ? ring_bond_count - 1 : 0);
        is_in_non_aromatic_ring( count) = ( aromatic_bond_count < ring_bond_count ? ( aromatic_bond_count == 1 ? ring_bond_count - 1 : 1) : 0);
        is_in_chain( count) = 3 - ring_bond_count;
        if( n_bonds( count) == n_e_in_bonds( count))
        {
          max_e_in_bonds( count) = min_e_in_bonds( count) = 2;
        }
        else
        {
          min_e_in_bonds( count) = 6;
          for
          (
            storage::Vector< chemistry::BondConformational>::const_iterator
              itr_bonds( itr->GetBonds().Begin()), itr_bonds_end( itr->GetBonds().End());
            itr_bonds != itr_bonds_end;
            ++itr_bonds
          )
          {
            if( itr_bonds->GetBondType()->GetNumberOfElectrons() > max_e_in_bonds( count))
            {
              max_e_in_bonds( count) = itr_bonds->GetBondType()->GetNumberOfElectrons();
            }
            if( itr_bonds->GetBondType()->GetNumberOfElectrons() < min_e_in_bonds( count))
            {
              min_e_in_bonds( count) = itr_bonds->GetBondType()->GetNumberOfElectrons();
            }
          }
        }
        atom_radius( count) = chemistry::BondLengths().GetAverageCovalentRadius( *itr);

        if( is_in_any_ring( count))
        {
          ring_size( count) = graph::Connectivity::LengthOfSmallestCycleWithVertex( graph, count);
        }
        electronegativity( count) = itr->GetElementType()->GetProperty( chemistry::ElementTypeData::e_ElectroNegativity);
      }

      if( m_WriteDescriptors->GetFlag())
      {
        size_t count_all( 0);

        for
        (
          iterate::Generic< const chemistry::AtomConformationalInterface> itr( FRAG.GetAtomsIterator());
          itr.NotAtEnd();
          ++itr, ++count_all
        )
        {
          if( n_bonds( count_all) == size_t( 2))
          {
          }
          else if
          (
            n_bonds( count_all) != size_t( 3)
            || ( itr->GetAtomType()->GetNumberElectronsInBonds() != size_t( 3) && itr->GetElementType()->GetPeriod() <= size_t( 2))
          )
          {
            continue;
          }
          if
          (
            itr->GetBonds().GetSize() != n_bonds( count_all)
            || ( skip_aromatic_rings && is_in_aromatic_ring( count_all))
          )
          {
            continue;
          }
          storage::VectorND< 2, linal::Vector< double> > feature_result( result_undef);
          linal::Vector< double> &feature( feature_result.First());
          linal::Vector< double> &result( feature_result.Second());
          ++count;

          size_t max_lone_pair_neighbor( 0);
          size_t min_lone_pair_neighbor( 4);
          size_t lone_pair_sums( 0);
          size_t lone_pair_products( 1);
          size_t max_n_bonds( 0);
          size_t min_n_bonds( 4);
          size_t bond_sum( 0);
          size_t max_e_in_neighbor_bonds( 0);
          size_t min_e_in_neighbor_bonds( 6);
          size_t neighbor_counts_in_aromatic_ring( 0);
          size_t neighbor_counts_in_non_aromatic_ring( 0);
          size_t neighbor_counts_in_any_ring( 0);
          size_t neighbor_counts_unsaturated( 0);
          size_t neighbor_period_3_plus( 0);
          size_t neighbor_group_7_count( 0);
          size_t neighbor_group_6_count( 0);
          size_t neighbor_group_5_count( 0);
          size_t neighbor_group_4_count( 0);
          size_t neighbor_group_3_count( 0);
          size_t neighbor_count_vsepr_4( 0);
          size_t neighbor_count_vsepr_3( 0);
          size_t neighbor_count_b( 0);
          size_t neighbor_count_c( 0);
          size_t neighbor_count_n( 0);
          size_t neighbor_count_o( 0);
          size_t neighbor_count_f( 0);
          size_t neighbor_count_si( 0);
          size_t neighbor_count_p( 0);
          size_t neighbor_count_s( 0);
          size_t neighbor_count_cl( 0);
          size_t neighbor_count_as( 0);
          size_t neighbor_count_se( 0);
          size_t neighbor_count_br( 0);

          size_t neighbor_count_c_saturated( 0);
          size_t neighbor_count_c_unsaturated( 0);
          size_t neighbor_count_n_def_saturated( 0);
          size_t neighbor_count_n_def_unsaturated( 0);
          size_t neighbor_count_n_pos_unsaturated( 0);
          size_t neighbor_count_o_def_unsaturated( 0);
          size_t neighbor_count_o_minus( 0);
          size_t neighbor_count_o_pos_unsaturated( 0);
          size_t neighbor_count_s_3_plus_e( 0);
          size_t neighbor_count_s_2_e( 0);
          int    sum_neighbor_formal_charge( 0);

          size_t neighbor_counts_unsaturated_non_te( 0);
          linal::Vector< double> bond_length( itr->GetAtomType()->GetNumberBonds(), double( 0.0));
          size_t count_bond( 0);
          storage::Vector< size_t> atom_indices( itr->GetAtomType()->GetNumberBonds());
          for
          (
            storage::Vector< chemistry::BondConformational>::const_iterator
              itr_bonds( itr->GetBonds().Begin()), itr_bonds_end( itr->GetBonds().End());
            itr_bonds != itr_bonds_end;
            ++itr_bonds, ++count_bond
          )
          {
            const size_t atom_index( FRAG.GetAtomIndex( itr_bonds->GetTargetAtom()));
            atom_indices( count_bond) = atom_index;
            if( n_lone_pairs( atom_index) > max_lone_pair_neighbor)
            {
              max_lone_pair_neighbor = n_lone_pairs( atom_index);
            }
            if( n_lone_pairs( atom_index) < min_lone_pair_neighbor)
            {
              min_lone_pair_neighbor = n_lone_pairs( atom_index);
            }
            lone_pair_sums += n_lone_pairs( atom_index);
            lone_pair_products *= n_lone_pairs( atom_index);
            bond_sum += n_bonds( atom_index);
            sum_neighbor_formal_charge += atom_types( atom_index)->GetFormalCharge();
            if( n_bonds( atom_index) > max_n_bonds)
            {
              max_n_bonds = n_bonds( atom_index);
            }
            if( n_bonds( atom_index) < min_n_bonds)
            {
              min_n_bonds = n_bonds( atom_index);
            }
            if( max_e_in_bonds( atom_index) > max_e_in_neighbor_bonds)
            {
              max_e_in_neighbor_bonds = max_e_in_bonds( atom_index);
            }
            if( min_e_in_bonds( atom_index) < min_e_in_neighbor_bonds)
            {
              min_e_in_neighbor_bonds = min_e_in_bonds( atom_index);
            }
            neighbor_counts_in_aromatic_ring += is_in_aromatic_ring( atom_index);
            neighbor_counts_in_non_aromatic_ring += is_in_non_aromatic_ring( atom_index);
            neighbor_counts_in_any_ring += is_in_any_ring( atom_index);
            bond_length( count_bond) =
              chemistry::BondLengths().GetBondLength
              (
                atom_types( count_all),
                itr_bonds->GetBondType()->GetBondData( chemistry::ConfigurationalBondTypeData::e_BondOrderOrAromatic),
                atom_types( atom_index)
              );
            if
            (
              max_e_in_bonds( atom_index) > 2
              || is_in_aromatic_ring( atom_index) > size_t( 0)
              || atom_types( atom_index)->GetElementType() == chemistry::GetElementTypes().e_Boron
            )
            {
              neighbor_counts_unsaturated += 1;
              if( n_bonds( atom_index) != size_t( 4))
              {
                neighbor_counts_unsaturated_non_te += 1;
              }
            }

            if( atom_types( atom_index)->GetElementType() == chemistry::GetElementTypes().e_Carbon)
            {
              if( n_bonds( atom_index) > 3.5)
              {
                ++neighbor_count_c_saturated;
              }
              else
              {
                ++neighbor_count_c_unsaturated;
              }
            }
            else if( atom_types( atom_index)->GetElementType() == chemistry::GetElementTypes().e_Nitrogen)
            {
              if( n_bonds( atom_index) > 3.5)
              {
                ++neighbor_count_n_def_saturated;
              }
              else if( is_in_aromatic_ring( atom_index) > 0.5 || n_e_in_bonds( atom_index) > 3.5)
              {
                ++neighbor_count_n_def_unsaturated;
              }
              else
              {
                ++neighbor_count_n_pos_unsaturated;
              }
            }
            else if( atom_types( atom_index)->GetElementType() == chemistry::GetElementTypes().e_Oxygen)
            {
              if( max_e_in_bonds( atom_index) > 2 || is_in_aromatic_ring( atom_index) > 0.5)
              {
                ++neighbor_count_o_def_unsaturated;
              }
              else if( atom_types( atom_index)->GetFormalCharge() == -1)
              {
                ++neighbor_count_o_minus;
              }
              else
              {
                ++neighbor_count_o_pos_unsaturated;
              }
            }
            else if( atom_types( atom_index)->GetElementType() == chemistry::GetElementTypes().e_Sulfur)
            {
              if( max_e_in_bonds( atom_index) > 2 || is_in_aromatic_ring( atom_index) > 0.5)
              {
                ++neighbor_count_s_3_plus_e;
              }
              else
              {
                ++neighbor_count_s_2_e;
              }
            }

            if( atom_types( atom_index)->GetElementType()->GetPeriod() > size_t( 2))
            {
              ++neighbor_period_3_plus;
            }
            if( atom_types( atom_index)->GetElementType()->GetMainGroup() <= size_t( 3))
            {
              ++neighbor_group_3_count;
              if( atom_types( atom_index)->GetElementType() == chemistry::GetElementTypes().e_Boron)
              {
                ++neighbor_count_b;
              }
            }
            else if( atom_types( atom_index)->GetElementType()->GetMainGroup() == size_t( 4))
            {
              ++neighbor_group_4_count;
              if( atom_types( atom_index)->GetElementType() == chemistry::GetElementTypes().e_Carbon)
              {
                ++neighbor_count_c;
              }
              else
              {
                ++neighbor_count_si;
              }
            }
            else if( atom_types( atom_index)->GetElementType()->GetMainGroup() == size_t( 5))
            {
              ++neighbor_group_5_count;
              if( atom_types( atom_index)->GetElementType() == chemistry::GetElementTypes().e_Nitrogen)
              {
                ++neighbor_count_n;
              }
              else if( atom_types( atom_index)->GetElementType() == chemistry::GetElementTypes().e_Phosphorus)
              {
                ++neighbor_count_p;
              }
              else
              {
                ++neighbor_count_as;
              }
            }
            else if( atom_types( atom_index)->GetElementType()->GetMainGroup() == size_t( 6))
            {
              ++neighbor_group_6_count;
              if( atom_types( atom_index)->GetElementType() == chemistry::GetElementTypes().e_Oxygen)
              {
                ++neighbor_count_o;
              }
              else if( atom_types( atom_index)->GetElementType() == chemistry::GetElementTypes().e_Sulfur)
              {
                ++neighbor_count_s;
              }
              else
              {
                ++neighbor_count_se;
              }
            }
            else
            {
              ++neighbor_group_7_count;
              if( atom_types( atom_index)->GetElementType() == chemistry::GetElementTypes().e_Fluorine)
              {
                ++neighbor_count_f;
              }
              else if( atom_types( atom_index)->GetElementType() == chemistry::GetElementTypes().e_Chlorine)
              {
                ++neighbor_count_cl;
              }
              else
              {
                ++neighbor_count_br;
              }
            }
            const size_t vsepr( n_bonds( atom_index) + n_lone_pairs( atom_index));
            if( vsepr == size_t( 3))
            {
              neighbor_count_vsepr_3 += 1;
            }
            else if( vsepr == size_t( 4))
            {
              neighbor_count_vsepr_4 += 1;
            }
          }
          size_t is_in_small_ring_medium_ring_macrocycle_or_no_cycle( 0);
          if( is_in_any_ring( count_all))
          {
            if( ring_size( count_all) >= 11)
            {
              is_in_small_ring_medium_ring_macrocycle_or_no_cycle = 2;
            }
            else if( ring_size( count_all) <= size_t( 5))
            {
              is_in_small_ring_medium_ring_macrocycle_or_no_cycle = 0;
            }
            else
            {
              is_in_small_ring_medium_ring_macrocycle_or_no_cycle = 1;
            }
          }
          else
          {
            is_in_small_ring_medium_ring_macrocycle_or_no_cycle = 3;
          }

          feature( 0 ) = molecule_charge;
          feature( 1 ) = ring_size( count_all) ? std::min( ring_size( count_all), size_t( 6)) : 7;
          feature( 2 ) = sum_neighbor_formal_charge;
          feature( 3 ) = n_lone_pairs( count_all);
          feature( 4 ) = is_in_aromatic_ring( count_all);
          feature( 5 ) = is_in_non_aromatic_ring( count_all);
          feature( 6 ) = is_in_any_ring( count_all);
          feature( 7 ) = neighbor_counts_unsaturated;
          feature( 8 ) = is_in_small_ring_medium_ring_macrocycle_or_no_cycle;
          feature( 9 ) = max_lone_pair_neighbor;
          feature( 10) = min_lone_pair_neighbor;
          feature( 11) = lone_pair_sums;
          feature( 12) = max_n_bonds;
          feature( 13) = min_n_bonds;
          feature( 14) = neighbor_period_3_plus - neighbor_count_s - neighbor_count_se;
          feature( 15) = neighbor_counts_unsaturated_non_te;
          feature( 16) = neighbor_counts_in_aromatic_ring;
          feature( 17) = neighbor_counts_in_non_aromatic_ring;
          feature( 18) = neighbor_counts_in_any_ring;
          feature( 19) = neighbor_period_3_plus;
          feature( 20) = neighbor_group_3_count;
          feature( 21) = neighbor_group_4_count;
          feature( 22) = neighbor_group_5_count;
          feature( 23) = neighbor_group_6_count;
          feature( 24) = neighbor_group_7_count;
          feature( 25) = neighbor_count_vsepr_3;
          feature( 26) = neighbor_count_vsepr_4;
          feature( 27) = neighbor_count_b;
          feature( 28) = neighbor_count_c;
          feature( 29) = neighbor_count_n;
          feature( 30) = neighbor_count_o;
          feature( 31) = neighbor_count_f;
          feature( 32) = neighbor_count_si;
          feature( 33) = neighbor_count_p;
          feature( 34) = neighbor_count_s;
          feature( 35) = neighbor_count_cl;
          feature( 36) = neighbor_count_as;
          feature( 37) = neighbor_count_se;
          feature( 38) = neighbor_count_br;
          feature( 39) = neighbor_count_c_saturated;
          feature( 40) = neighbor_count_c_unsaturated;
          feature( 41) = neighbor_count_n_def_saturated;
          feature( 42) = neighbor_count_n_def_unsaturated;
          feature( 43) = neighbor_count_n_pos_unsaturated;
          feature( 44) = neighbor_count_o_def_unsaturated;
          feature( 45) = neighbor_count_o_minus;
          feature( 46) = neighbor_count_o_pos_unsaturated;
          feature( 47) = neighbor_count_s_3_plus_e;
          feature( 48) = neighbor_count_s_2_e;
          // compute planarity
          const linal::Vector3D &position_a( itr->GetPosition());
          const linal::Vector3D &position_b( itr->GetBonds()( 0).GetTargetAtom().GetPosition());
          const linal::Vector3D &position_c( itr->GetBonds()( 1).GetTargetAtom().GetPosition());
          const double angle_bac( math::Angle::Degree( linal::ProjAngle( position_a, position_b, position_c)));
          if( n_bonds( count_all) == size_t( 3))
          {
            feature( 49) = std::max( std::max( electronegativity( atom_indices( 0)), electronegativity( atom_indices( 1))), electronegativity( atom_indices( 2)));
            feature( 50) = ( electronegativity( atom_indices( 0)) + electronegativity( atom_indices( 1)) + electronegativity( atom_indices( 2))) / 3.0;
            feature( 51) = std::min( std::min( electronegativity( atom_indices( 0)), electronegativity( atom_indices( 1))), electronegativity( atom_indices( 2)));
            feature( 52) = std::min( std::min( bond_length( 0), bond_length( 1)), bond_length( 2));
            feature( 53) = ( bond_length( 0) + bond_length( 1) + bond_length( 2)) / 3.0;
            feature( 54) = std::max( std::max( bond_length( 0), bond_length( 1)), bond_length( 2));

            const linal::Vector3D &position_d( itr->GetBonds()( 2).GetTargetAtom().GetPosition());

            double angle_cad( math::Angle::Degree( linal::ProjAngle( position_a, position_c, position_d)));
            double angle_dab( math::Angle::Degree( linal::ProjAngle( position_a, position_d, position_b)));

            const double max_angle( std::max( angle_bac, std::max( angle_cad, angle_dab)));
            double angle_sum( angle_bac + angle_cad + angle_dab);
            const double min_angle_sum( angle_sum - max_angle);

            if( 360.0 - max_angle + min_angle_sum < 361.0 && math::EqualWithinAbsoluteTolerance( max_angle, min_angle_sum, 1.0))
            {
              BCL_MessageStd
              (
                "Messed up geometry: " + util::Format()( max_angle)
                + " " + util::Format()( min_angle_sum)
                + " " + util::Format()( angle_bac)
                + " " + util::Format()( angle_cad)
                + " " + util::Format()( angle_dab)
              );
              angle_sum = 360.0 - max_angle + min_angle_sum;
            }

            result( 0) = angle_sum / 3.0;
          }
          else
          {
            feature( 49) = std::max( electronegativity( atom_indices( 0)), electronegativity( atom_indices( 1)));
            feature( 50) = ( electronegativity( atom_indices( 0)) + electronegativity( atom_indices( 1))) / 2.0;
            feature( 51) = std::min( electronegativity( atom_indices( 0)), electronegativity( atom_indices( 1)));
            feature( 52) = std::min( bond_length( 0), bond_length( 1));
            feature( 53) = ( bond_length( 0) + bond_length( 1)) / 3.0;
            feature( 54) = std::max( bond_length( 0), bond_length( 1));
            result( 0) = angle_bac;
          }
          // the most accurate separator is
          // feature( 30) = 1.0017 * ( 0.92419 * ( 148.986 + 3.02961 * feature( 7) - 0.44094 * feature( 28) - 13.2087 * feature( 32)) + 0.03048 * feature( 10) + 3.07132 * feature( 33)) - 0.0813809 * feature( 27);
          // however, this version is boosted very little with decision trees, so the resulting accuracy is not as good sa
          // the version given below
          //feature( 30) = 148.986 + 3.02961 * feature( 7) - 13.2087 * feature( 32);

          size_t node( 0);
          std::string predicted_hyb;
          if( atom_types( count_all)->GetElementType() == chemistry::GetElementTypes().e_Nitrogen)
          {
            if
            (
              ( n_bonds( count_all) == size_t( 3) && feature( 1) == size_t( 3))
              ||
              ( n_bonds( count_all) == size_t( 2) && feature( 1) < size_t( 6))
            )
            {
              node = 1;
              predicted_hyb = "Tet";
            }
            else if( feature( 7) > 1.5)
            {
              node = 2;
              predicted_hyb = "Tri";
            }
            else if( feature( 24) > 0.5)
            {
              node = 3;
              predicted_hyb = "Tet";
            }
            else if( feature( 7) > 0.5)
            {
              node = 4;
              predicted_hyb = "Tri";
            }
            else if( feature( 14) > 0.5)
            {
              node = 5;
              predicted_hyb = "Tri";
            }
            else if( feature( 28) > 2.5 || feature( 6) > 1.5)
            {
              node = 6;
              predicted_hyb = "Tet";
            }
            else if( ring_size( count_all) == size_t( 6))
            {
              if( feature( 30) > 0.5)
              {
                node = 7;
                predicted_hyb = "Tri";
              }
              else
              {
                node = 8;
                predicted_hyb = "Tet";
              }
            }
            else if( feature( 28) < 1.5)
            {
              node = 9;
              predicted_hyb = "Tri";
            }
            else
            {
              node = 10;
              predicted_hyb = "Tet";
            }
          }
          else if( atom_types( count_all)->GetElementType() == chemistry::GetElementTypes().e_Oxygen)
          {
            if
            (
              ( n_bonds( count_all) == size_t( 3) && feature( 1) == size_t( 3))
              ||
              ( n_bonds( count_all) == size_t( 2) && feature( 1) < size_t( 6))
            )
            {
              node = 1;
              predicted_hyb = "Tet";
            }
            else if( feature( 29) > 0.5 || feature( 24) > 0.5 || feature( 30) > 0.5 || feature( 48) > 0.5)
            {
              node = 2;
              predicted_hyb = "Tet";
            }
            else if( feature( 28) < 1.5)
            {
              node = 3;
              predicted_hyb = "Tri";
            }
            else if( feature( 7) > 0.5)
            {
              node = 4;
              predicted_hyb = "Tri";
            }
            else
            {
              node = 5;
              predicted_hyb = "Tet";
            }
          }
          BCL_MessageStd
          (
            atom_types( count_all)->GetElementType()->GetChemicalSymbol()
            + util::Format()( n_bonds( count_all)) + util::Format()( n_e_in_bonds( count_all)) + " " + util::Format()( node)
            + " average angle: " + util::Format()( result( 0)) + " " + util::Format()( ring_size( count_all))
            + " " + util::Format()( is_in_aromatic_ring( count_all)) + " " + predicted_hyb + " " + util::Format()( feature( 6))
          );

          if( ring_size( count_all) != size_t( 3))
          {
            const std::string type_str
            (
              atom_types( count_all)->GetElementType()->GetChemicalSymbol()
              + util::Format()( n_bonds( count_all)) + util::Format()( n_e_in_bonds( count_all))
            );
            m_AngleDescriptorsByType[ type_str].PushBack( feature_result);
          }
        }
      }

      if( m_WriteBondLengthsInitializer->GetFlag())
      {
        chemistry::ConformationGraphConverter::t_AtomGraph molgraph( chemistry::ConformationGraphConverter::CreateGraphWithAtoms( FRAG));
        chemistry::FragmentSplitRigid splitter( size_t( 1), false);
        storage::List< storage::Vector< size_t> > components( splitter.GetComponentVertices( FRAG, molgraph));
        size_t rigid_id( 0);
        storage::Vector< storage::Vector< size_t> > rigid_components( FRAG.GetSize());
        for
        (
          storage::List< storage::Vector< size_t> >::const_iterator
            itr_rigid( components.Begin()), itr_rigid_end( components.End());
          itr_rigid != itr_rigid_end;
          ++itr_rigid, ++rigid_id
        )
        {
          for
          (
            storage::Vector< size_t>::const_iterator
              itr_rigid_atom( itr_rigid->Begin()), itr_rigid_atom_end( itr_rigid->End());
            itr_rigid_atom != itr_rigid_atom_end;
            ++itr_rigid_atom
          )
          {
            rigid_components( *itr_rigid_atom).PushBack( rigid_id);
          }
        }

        static storage::Map< std::string, std::string> s_substitutions;
        if( s_substitutions.IsEmpty())
        {
          s_substitutions[ "N1"] = "NX";
          s_substitutions[ "N2"] = "NX";
          s_substitutions[ "S1" ] = "SX";
          s_substitutions[ "S2" ] = "SX";
        }

        // compute the bond lengths
        count = 0;
        for
        (
          iterate::Generic< const chemistry::AtomConformationalInterface> itr( FRAG.GetAtomsIterator());
          itr.NotAtEnd();
          ++itr, ++count
        )
        {
          const storage::Vector< chemistry::BondConformational> &connected_atoms( itr->GetBonds());
          storage::Vector< size_t> distances( *graph::Connectivity::DistancesToOtherVertices( graph, count));
          const bool is_in_ring
          (
            itr->CountNonValenceBondsWithProperty( chemistry::ConfigurationalBondTypeData::e_IsInRing, 1) >= 2
          );

          for
          (
            storage::Vector< chemistry::BondConformational>::const_iterator
              itr_connected( connected_atoms.Begin()),
              itr_connected_end( connected_atoms.End());
            itr_connected != itr_connected_end;
            ++itr_connected
          )
          {
            if( &*itr > &itr_connected->GetTargetAtom())
            {
              continue;
            }

            const size_t atom_index_b( FRAG.GetAtomIndex( itr_connected->GetTargetAtom()));

            chemistry::AtomType atom_type_a( itr->GetAtomType());
            chemistry::AtomType atom_type_b( itr_connected->GetTargetAtom().GetAtomType());
            if( atom_type_a > atom_type_b)
            {
              std::swap( atom_type_a, atom_type_b);
            }
            chemistry::ConfigurationalBondType
              bond_type( itr_connected->GetBondType()->WithoutAromaticOrder()->WithoutIsometry());

            const size_t min_ring_size
            (
              bond_type->IsBondInRing() && std::max( ring_size( count), ring_size( atom_index_b)) < s_NumberRingTypes + size_t( 1)
              ? std::max( ring_size( count), ring_size( atom_index_b))
              : size_t( 2)
            );
            const std::string base_atom_type_a_name
            (
              atom_type_a->GetElementType()->GetChemicalSymbol()
              + util::Format()( atom_type_a->GetNumberBonds())
              + util::Format()( atom_type_a->GetNumberElectronsInBonds())
            );
            const std::string base_atom_type_b_name
            (
              atom_type_b->GetElementType()->GetChemicalSymbol()
              + util::Format()( atom_type_b->GetNumberBonds())
              + util::Format()( atom_type_b->GetNumberElectronsInBonds())
            );

            const size_t bond_type_id( GetBasicBondTypeId( bond_type));

            m_AtomTypesForSimpleBondType( min_ring_size - 2)( bond_type_id).Insert( base_atom_type_a_name);
            m_AtomTypesForSimpleBondType( min_ring_size - 2)( bond_type_id).Insert( base_atom_type_b_name);

            const std::string bond_length_description( base_atom_type_a_name + ' ' + base_atom_type_b_name);
            math::RunningAverageSD< double> &statistic
            (
              m_BondLengths( min_ring_size - 2)( bond_type_id)[ bond_length_description]
            );
            statistic += linal::Distance( itr->GetPosition(), itr_connected->GetTargetAtom().GetPosition());
          }

          for( size_t desired_dist( 2), max_dist( 5); desired_dist < max_dist; ++desired_dist)
          {
            size_t count_b( 0);
            double min_dist( 4.54);
            iterate::Generic< const chemistry::AtomConformationalInterface> itr_closest;
            storage::Vector< storage::Pair< iterate::Generic< const chemistry::AtomConformationalInterface>, double> > not_the_closest;
            size_t was_closest( util::GetUndefinedSize_t());
            for
            (
              iterate::Generic< const chemistry::AtomConformationalInterface> itr_b( FRAG.GetAtomsIterator());
              itr_b.NotAtEnd();
              ++itr_b, ++count_b
            )
            {
              if( count_b == count)
              {
                continue;
              }
              if( distances( count_b) < desired_dist)
              {
                continue;
              }
              if( distances( count_b) != desired_dist && desired_dist != size_t( 4))
              {
                continue;
              }
              if( distances( count_b) != size_t( 2) && HaveOverlap( rigid_components( count), rigid_components( count_b)))
              {
                continue;
              }

              const double dist
              (
                linal::Distance( itr->GetPosition(), itr_b->GetPosition())
              );
              not_the_closest.PushBack
              (
                storage::Pair< iterate::Generic< const chemistry::AtomConformationalInterface>, double>
                (
                  itr_b,
                  dist
                )
              );
              if( dist < min_dist)
              {
                min_dist = dist;
                itr_closest = itr_b;
                was_closest = count_b;
              }
            }
            std::string atom_type_a_name( itr->GetAtomType()->GetTwoLetterCode());
            if( !util::IsDefined( was_closest))
            {
              m_VdwBondLengthsByBondDist
              (
                std::min( desired_dist - size_t( 2), m_VdwBondLengthsByBondDist.GetSize() - size_t( 1))
              )[ atom_type_a_name + " XX"].PushBack( min_dist);
            }
            else
            {
              m_TimesNotTheClosest
              (
                std::min( desired_dist - size_t( 2), m_VdwBondLengthsByBondDist.GetSize() - size_t( 1))
              )[ atom_type_a_name + " XX"] += 1.0;
            }

            for( auto itr_ntc( not_the_closest.Begin()), itr_ntc_end( not_the_closest.End()); itr_ntc != itr_ntc_end; ++itr_ntc)
            {
              iterate::Generic< const chemistry::AtomConformationalInterface> itr_b( itr_ntc->First());

              std::string atom_type_b_name( itr_b->GetAtomType()->GetTwoLetterCode());

              const std::string bond_length_description
              (
                atom_type_a_name < atom_type_b_name
                ? atom_type_a_name + ' ' + atom_type_b_name
                : atom_type_b_name + ' ' + atom_type_a_name
              );

              double &st( m_ExpectedLengthsByBondDist[ bond_length_description]);
              if( st == 0.0 || !util::IsDefined( st))
              {
                st = chemistry::BondLengths::GetAverageCovalentRadius( *itr)
                     + chemistry::BondLengths::GetAverageCovalentRadius( *itr_b);
              }
              m_ExpectedVdw[ bond_length_description] =
                       itr->GetElementType()->GetProperty( chemistry::ElementTypeData::e_VDWaalsRadius)
                       + itr_b->GetElementType()->GetProperty( chemistry::ElementTypeData::e_VDWaalsRadius);
              if
              (
                ( itr->GetNumberCovalentlyBoundHydrogens() + itr_b->GetNumberCovalentlyBoundHydrogens())
                &&
                std::max
                (
                  itr->GetElementType()->GetProperty( chemistry::ElementTypeData::e_ElectroNegativity),
                  itr_b->GetElementType()->GetProperty( chemistry::ElementTypeData::e_ElectroNegativity)
                ) > 2.8
              )
              {
                continue;
              }
              if( itr_ntc->Second() < m_ExpectedVdw[ bond_length_description] + 0.1)
              {
                m_VdwBondLengths[ bond_length_description].PushBack( min_dist);
                m_VdwBondLengthsByBondDist
                (
                  std::min( desired_dist - size_t( 2), m_VdwBondLengthsByBondDist.GetSize() - size_t( 1))
                )[ bond_length_description].PushBack( min_dist);
              }
              else
              {
                m_TimesNotTheClosest
                (
                  std::min( desired_dist - size_t( 2), m_VdwBondLengthsByBondDist.GetSize() - size_t( 1))
                )[ bond_length_description] += 1.0 / double( not_the_closest.GetSize());
              }
            }
          }
        }
      }
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &GenerateAtomHybridizationDescriptors::GetReadMe() const
    {
      static std::string s_read_me =
        "GenerateAtomHybridizationDescriptors generates atom descriptors used for building models that determine"
        "the orbital hybridization state of atom types and, determining covalent and van der waals radii from a "
        "molecule structure library, such as the cambridge structural database";
      return s_read_me;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int GenerateAtomHybridizationDescriptors::Main() const
    {
      io::OFStream output;
      m_VdwBondLengthsByBondDist.Reset();
      m_VdwBondLengthsByBondDist.Resize( size_t( 3));
      m_TimesNotTheClosest.Reset();
      m_TimesNotTheClosest.Resize( size_t( 3));
      m_BondLengths.Reset();
      m_AtomTypesForSimpleBondType.Reset();
      m_BondLengths.Resize
      (
        s_NumberRingTypes,
        storage::Vector< storage::Map< std::string, math::RunningAverageSD< double> > >( s_NumberBondTypes)
      );
      m_AtomTypesForSimpleBondType.Resize
      (
        s_NumberRingTypes,
        storage::Vector< storage::Set< std::string> >( s_NumberBondTypes)
      );

      // iterate over the fragments from whatever feeds were given
      chemistry::FragmentEnsemble molecules_from_feed;
      for( chemistry::FragmentFeed feed; feed.NotAtEnd(); ++feed)
      {
        // generate all the descriptors

        // remove H first because they are improperly resolved by x-ray scattering experiments
        chemistry::FragmentComplete frag( *feed);
        frag.SaturateWithH();

        // Skip molecules with non-gasteiger atom types
        if( frag.HasNonGasteigerAtomTypes())
        {
          BCL_MessageStd
          (
            "Molecule #" + util::Format()( feed.GetPosition())
            + " has atoms that cannot be described with gasteiger atom types; skipping"
          );
          continue;
        }

        // Skip molecules with bad 3d coordinates
//        if( frag.HasBadGeometry())
//        {
//          BCL_MessageStd
//          (
//            "Molecule #" + util::Format()( feed.GetPosition()) + " has invalid 3d geometry; skipping"
//          );
//          continue;
//        }
        Process( frag);
        frag.SaturateWithH();
        molecules_from_feed.PushBack( frag);
      }

      // Collect stats
      if( m_WriteElementBondStatistics->GetFlag())
      {
        CountElementBonds( molecules_from_feed);
      }

      const std::string file_prefix
      (
        m_OutputFilenameBase->GetFirstParameter()->GetValue() + "."
      );

      if( m_WriteDescriptors->GetFlag())
      {
        const std::string file_suffix
        (
          "angle_descriptors.w" + std::string( m_SkipAromaticRingDescriptors->GetFlag() ? "o_aromatic.txt" : "_aromatic.txt")
        );

        // track all the txt filenames
        storage::Vector< std::string> all_filenames;

        for
        (
          storage::Map< std::string, storage::Vector< storage::VectorND< 2, linal::Vector< double> > > >::const_iterator
            itr( m_AngleDescriptorsByType.Begin()), itr_end( m_AngleDescriptorsByType.End());
          itr != itr_end;
          ++itr
        )
        {
          const std::string filename( file_prefix + itr->first + "." + file_suffix);
          all_filenames.PushBack( filename);
          io::File::MustOpenOFStream( output, filename);
          io::Serialize::Write( itr->second, output);
          io::File::CloseClearFStream( output);
        }

        for
        (
          storage::Vector< std::string>::iterator itr( all_filenames.Begin()), itr_end( all_filenames.End());
          itr != itr_end;
          ++itr
        )
        {
          // initialize the source dataset retriever
          util::Implementation< model::RetrieveDataSetBase> source( "File(filename=" + *itr + ")");

          // generate the output filename by changing the last 3 characters to bin
          itr->replace( itr->size() - 3, 3, "bin");

          // store the dataset
          model::RetrieveDatasetSubset::StoreMasterDataset( *itr, *source);
        }
      }

      if( m_WriteBondLengthsInitializer->GetFlag())
      {
        io::File::MustOpenOFStream
        (
          output,
          m_OutputFilenameBase->GetFirstParameter()->GetValue() + ".bond_lengths_initializer.txt"
        );
        WriteBondLengthsInitializer( output);
        io::File::CloseClearFStream( output);
      }
      // end
      return 0;
    }

    //! @brief generate counts of element-element bonds for statistics
    //! @param FRAG the fragment for which to calculate the counts
    void GenerateAtomHybridizationDescriptors::CountElementBonds( const chemistry::FragmentEnsemble &FRAG) const
    {
      // Initialize total bond counter for the current ensemble
      size_t n_total_bonds( 0);

      // ~N_total_eletypes^2 * 4
//      float normalization_factor( 1600.0);

      // ~N_total_eletypes^2 * 10 bond types (minus the undefined)
      float normalization_factor( 4000.0);

      // make a map and increment the value every time we see the key
      storage::Map< storage::Triplet< chemistry::ElementType, chemistry::ElementType, chemistry::ConfigurationalBondType>, size_t> counts;
      storage::Map< storage::Pair< chemistry::ElementType, chemistry::ConfigurationalBondType>, size_t> bg_counts;

      // add pseudo count for unknowns
      storage::Triplet< chemistry::ElementType, chemistry::ElementType, chemistry::ConfigurationalBondType> unknown
      (
        chemistry::GetElementTypes().e_Undefined,
        chemistry::GetElementTypes().e_Undefined,
        chemistry::GetConfigurationalBondTypes().e_Undefined
      );
      counts.Insert( std::make_pair( unknown, size_t( 0)));

      // collect statistics
      chemistry::ConfigurationalBondType canonical_bonds[ 11] =
      {
          chemistry::GetConfigurationalBondTypes().e_Undefined,
          chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBond,
          chemistry::GetConfigurationalBondTypes().e_ConjugatedSingleBond,
          chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBond,
          chemistry::GetConfigurationalBondTypes().e_ConjugatedTripleBond,
          chemistry::GetConfigurationalBondTypes().e_NonConjugatedSingleBondInRing,
          chemistry::GetConfigurationalBondTypes().e_ConjugatedSingleBondInRing,
          chemistry::GetConfigurationalBondTypes().e_ConjugatedDoubleBondInRing,
          chemistry::GetConfigurationalBondTypes().e_ConjugatedTripleBondInRing,
          chemistry::GetConfigurationalBondTypes().e_AmideSingleBond,
          chemistry::GetConfigurationalBondTypes().e_AromaticBond
      };

      for( auto ens_itr( FRAG.Begin()), ens_itr_end( FRAG.End()); ens_itr != ens_itr_end; ++ens_itr)
      {
        // increment total number of bonds across ensemble
        n_total_bonds += ens_itr->GetNumberBonds();

        // collect the number of each element-element-bond triplet instances and total element-element pair instances
        // iterate over atoms in molecule
        for( auto atom_itr( ens_itr->GetAtomsIterator()); atom_itr.NotAtEnd(); ++atom_itr)
        {
          // iterate over atoms bonded to the current atom
          for( auto bond_itr( atom_itr->GetBonds().Begin()), bond_itr_end( atom_itr->GetBonds().End()); bond_itr != bond_itr_end; ++bond_itr)
          {
            // make a key from current atom, bonded atom, and bond type
            storage::Triplet< chemistry::ElementType, chemistry::ElementType, chemistry::ConfigurationalBondType> key
            (
              atom_itr->GetElementType(),
              bond_itr->GetTargetAtom().GetElementType(),
              canonical_bonds[ bond_itr->GetBondType()->GetBondData( chemistry::ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness)]
            );
            // insert new key/value pairs, or increment existing
            auto insert_itr( counts.Insert( std::make_pair( key, size_t( 1))));
            if( !insert_itr.second)
            {
              insert_itr.first->second += 1;
            }
            // don't forget to grab background data
            storage::Pair< chemistry::ElementType, chemistry::ConfigurationalBondType> bg_key
            (
              atom_itr->GetElementType(),
              canonical_bonds[ bond_itr->GetBondType()->GetBondData( chemistry::ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness)]
            );
            //insert new key/value pairs for background data, or increment existing
            auto bg_insert_itr( bg_counts.Insert( std::make_pair( bg_key, size_t( 1))));
            if( !bg_insert_itr.second)
            {
              bg_insert_itr.first->second += 1;
            }
          }
        }
      }

      // iterate over our fancy new map to collect counts
      storage::Vector
      <
        storage::Pair
        <
          storage::Triplet< chemistry::ElementType, chemistry::ElementType, chemistry::ConfigurationalBondType>,
          size_t
        >
      > raw_counts( counts.GetSize());
      storage::Vector
      <
        storage::Pair
        <
          storage::Triplet< chemistry::ElementType, chemistry::ElementType, chemistry::ConfigurationalBondType>,
          float
        >
      > normalized_counts( counts.GetSize()), log_normalized_counts( counts.GetSize());
      storage::Map
      <
        storage::Triplet< chemistry::ElementType, chemistry::ElementType, chemistry::ConfigurationalBondType>,
        float
      > energies;
      size_t index( 0);
      for( auto counts_itr( counts.Begin()), counts_itr_end( counts.End()); counts_itr != counts_itr_end; ++counts_itr, ++index)
      {
        // accumulate raw counts
        raw_counts( index).First() = counts_itr->first;
        raw_counts( index).Second() = counts_itr->second;

        // accumulate normalized counts
        normalized_counts( index).First() = counts_itr->first;
        normalized_counts( index).Second() = float( counts_itr->second) /
            (
                float( bg_counts.Find( std::make_pair( counts_itr->first.First(), counts_itr->first.Third()))->second)
                * float( bg_counts.Find( std::make_pair( counts_itr->first.Second(), counts_itr->first.Third()))->second)
            ) / ( float( n_total_bonds) / normalization_factor);

        // accumulate log normalized counts
        log_normalized_counts( index).First() = counts_itr->first;
        log_normalized_counts( index).Second() = log10f
            (
              float( counts_itr->second) /
              (
                  float( bg_counts.Find( std::make_pair( counts_itr->first.First(), counts_itr->first.Third()))->second)
                  * float( bg_counts.Find( std::make_pair( counts_itr->first.Second(), counts_itr->first.Third()))->second)
              ) / ( float( n_total_bonds) / normalization_factor)
            );

        // accumulate final energies as a map
        energies[ counts_itr->first] = -log10f
            (
              float( counts_itr->second + 1) / ( float( n_total_bonds) / normalization_factor)
            );
      }
      // output counts
      if( m_OutputFilenameBase->GetFlag())
      {
        io::OFStream raw_counts_out;
        io::File::MustOpenOFStream( raw_counts_out, m_OutputFilenameBase->GetFirstParameter()->GetValue() + ".element_bond_counts.raw.txt");
        raw_counts_out << raw_counts << '\n';
        io::File::CloseClearFStream( raw_counts_out);

        io::OFStream counts_out;
        io::File::MustOpenOFStream( counts_out, m_OutputFilenameBase->GetFirstParameter()->GetValue() + ".element_bond_counts.normalized.txt");
        counts_out << normalized_counts << '\n';
        io::File::CloseClearFStream( counts_out);

        io::OFStream log_counts_out;
        io::File::MustOpenOFStream( log_counts_out, m_OutputFilenameBase->GetFirstParameter()->GetValue() + ".element_bond_counts.log_normalized.txt");
        log_counts_out << log_normalized_counts << '\n';
        io::File::CloseClearFStream( log_counts_out);

        io::OFStream energies_out;
        io::File::MustOpenOFStream( energies_out, m_OutputFilenameBase->GetFirstParameter()->GetValue() + ".element_bond_counts.energies.txt");
        energies_out << energies << '\n';
        io::File::CloseClearFStream( energies_out);
      }
    }

    //! @brief write the initializers for chemistry::BondLengths
    //! @param OUTPUT output stream to write to
    void GenerateAtomHybridizationDescriptors::WriteBondLengthsInitializer( std::ostream &OUTPUT) const
    {
      // maintain a map from atom type to the maximum covalent bond length observed for that atom type
      storage::Map< std::string, double> atom_type_to_max_value;
      std::ostringstream oss_obj, oss_model;
      bool isfirst( true);
      for( size_t i( 0), sz( m_VdwBondLengthsByBondDist.GetSize()); i < sz; ++i)
      {
        storage::Map< std::pair< size_t, size_t>, std::pair< double, size_t> > atom_types_to_vdw;
        storage::Map< std::string, size_t> atom_type_to_index;
        size_t curr_atom_type_index( 0);
        for( auto itr( m_VdwBondLengthsByBondDist( i).Begin()), itr_end( m_VdwBondLengthsByBondDist( i).End()); itr != itr_end; ++itr)
        {
          itr->second.Sort( std::less< double>());
          if( !itr->second.GetSize())
          {
            continue;
          }
          const size_t original_size( itr->second.GetSize());
          std::ostringstream oss;
          double expected_vdw( m_ExpectedVdw[ itr->first]);
          oss << "Vdw stuff: " << itr->first << " " << i + 2 << " " << original_size << " ";
          size_t sz_not_closest( m_TimesNotTheClosest( i)[ itr->first]);
          oss << itr->second(  itr->second.GetSize() / 200) << ' '
              << itr->second(  itr->second.GetSize() / 2000) << ' '
              << m_ExpectedLengthsByBondDist[ itr->first] << ' '
              << expected_vdw << ' ' << sz_not_closest << ' ';
          size_t n_below_vwd( 0);
          math::RunningAverageSD< double> avesd;
          for( ; n_below_vwd < original_size; ++n_below_vwd)
          {
            avesd += itr->second( n_below_vwd);
          }
          itr->second.Resize( n_below_vwd);
          if( !itr->second.GetSize())
          {
            continue;
          }
          double median( itr->second( itr->second.GetSize() / 2.0));
          oss << float( sz_not_closest) / float( original_size + sz_not_closest)
              << ' ' << avesd.GetAverage() << ' ' << avesd.GetStandardDeviation() << ' '
              << median << ' ';
          if( itr->second.GetSize() < size_t( 100))
          {
            //continue;
            for( auto itrb( itr->second.Begin()), itrb_end( itr->second.End()); itrb < itrb_end; itrb += 2)
            {
              oss << *itrb << ' ';
            }
          }
          else
          {
            for( size_t j( 0), jsz( 100); j < jsz; j += 2)
            {
              oss << itr->second( j) << ' ';
            }
          }
          if( original_size > 0.0)
          {
            auto split( util::SplitString( itr->first, " \t"));
            if( itr->first.find( "XX") == std::string::npos && itr->first.find( "HX") == std::string::npos)
            {
              if( !atom_type_to_index.Has( split( 0)))
              {
                atom_type_to_index[ split( 0)] = curr_atom_type_index++;
              }
              if( !atom_type_to_index.Has( split( 1)))
              {
                atom_type_to_index[ split( 1)] = curr_atom_type_index++;
              }
              atom_types_to_vdw[ std::make_pair( atom_type_to_index[ split( 0)], atom_type_to_index[ split( 1)])]
                = std::make_pair( itr->second( itr->second.GetSize() / 200), original_size);
              if( m_ExpectedLengthsByBondDist[ itr->first] + 0.2 < itr->second( itr->second.GetSize() / 200))
              {
                oss_model << m_ExpectedLengthsByBondDist[ itr->first] + 0.01 << " <= "
                          << split( 0) << i+2 << " + "
                          << split( 1) << i+2 << " <= "
                          << itr->second( itr->second.GetSize() / 200) << ";\n";
              }
              if( !isfirst)
              {
                oss_obj << " + ";
              }
              else
              {
                isfirst = false;
                oss_obj << "max: ";
              }
              oss_obj << split( 0) << i + 2 << " " << original_size << " + " << split( 1) << i + 2 << " " << original_size;
            }
          }
          BCL_MessageStd( oss.str());
        }
        size_t n_records( atom_types_to_vdw.GetSize());
        linal::Matrix< double> vdw_length( n_records, curr_atom_type_index, double( 0.0));
        linal::Vector< double> lengths( n_records);
        size_t row( 0);
        for( auto itr( atom_types_to_vdw.Begin()), itr_end( atom_types_to_vdw.End()); itr != itr_end; ++itr, ++row)
        {
          vdw_length( row, itr->first.first) += itr->second.second;
          vdw_length( row, itr->first.second) += itr->second.second;
          lengths( row) = itr->second.first * itr->second.second;
        }
        storage::Vector< std::string> atom_type_names( curr_atom_type_index);
        for( auto itr( atom_type_to_index.Begin()), itr_end( atom_type_to_index.End()); itr != itr_end; ++itr)
        {
          atom_type_names( itr->second) = itr->first;
        }
        const linal::Vector< double> orbital_sizes
        (
          math::LinearLeastSquares::SolutionAndChiSquared( vdw_length, lengths).First()
        );
        for( size_t k( 0), sz( orbital_sizes.GetSize()); k < sz; ++k)
        {
          BCL_MessageStd
          (
            "VDW for size " + util::Format()( i + 2) + " "
            + atom_type_names( k) + " = " + util::Format()( orbital_sizes( k))
          );
        }
      }
      oss_obj << ";\n";
      io::OFStream outm;
      io::File::MustOpenOFStream( outm, "vdwnewmodel.lp.model");
      outm << oss_obj.str();
      outm << oss_model.str();
      io::File::CloseClearFStream( outm);
      return;

//      BCL_Debug( m_BondLengths);
      for( size_t ring_size_minus_two( 0); ring_size_minus_two < s_NumberRingTypes; ++ring_size_minus_two)
      {
        for( size_t bond_type_id( 1); bond_type_id < s_NumberBondTypes; ++bond_type_id)
        {
          const chemistry::AtomTypeData::Properties
            bond_type_property
            (
              chemistry::AtomTypeData::Properties( chemistry::AtomTypeData::e_VdWaalsRadiusCSD + bond_type_id)
            );
          const std::string bond_property_name
          (
            "AtomTypeData::e_" + chemistry::AtomTypeData::GetPropertyName( bond_type_property)
          );
          const storage::Set< std::string> &atom_types
          (
            m_AtomTypesForSimpleBondType( ring_size_minus_two)( bond_type_id)
          );
          storage::Map< std::string, size_t> atom_type_to_index;
          size_t counter( 0);
          for
          (
            storage::Set< std::string>::const_iterator itr( atom_types.Begin()), itr_end( atom_types.End());
            itr != itr_end;
            ++itr, ++counter
          )
          {
            atom_type_to_index[ *itr] = counter;
          }
          const size_t number_atom_types( atom_types.GetSize());
          const storage::Map< std::string, math::RunningAverageSD< double> > &bond_lengths
          (
            m_BondLengths( ring_size_minus_two)( bond_type_id)
          );
          const size_t number_bond_lengths( bond_lengths.GetSize());

          if( !number_bond_lengths)
          {
            continue;
          }
          BCL_MessageVrb
          (
            "# Atom types, bond lengths, counts for ring size " + util::Format()( ring_size_minus_two)
            + ", BT= " + util::Format()( bond_type_id) + ": "
            + util::Format()( number_atom_types) + " "
            + util::Format()( number_bond_lengths)
          );
          if( number_atom_types > number_bond_lengths)
          {
            BCL_MessageStd( "Under-determined system, skipping");
            continue;
          }
          linal::Matrix< double> bond_length_incidences( number_bond_lengths, number_atom_types, double( 0.0));
          linal::Vector< double> bond_length_aves( number_bond_lengths, double( 0.0));
          counter = 0;
          size_t sum_counts( 0);
          for
          (
            storage::Map< std::string, math::RunningAverageSD< double> >::const_iterator
              itr( bond_lengths.Begin()), itr_end( bond_lengths.End());
            itr != itr_end;
            ++itr, ++counter
          )
          {
            const math::RunningAverageSD< double> &average( itr->second);
            double weight( average.GetWeight());
            sum_counts += weight;
//            double variance( average.GetVariance());
//            if( weight < 10.0)
//            {
//              // if < 5 samples are known, then do not use the variance to inflate the weight
//              variance = 1.0;
//            }
//            else if( variance < 0.01)
//            {
//              // avoid numerical issues by setting a maximum variance
//              variance = 0.01;
//            }
//            weight /= variance;
            bond_length_aves( counter) = average.GetAverage() * weight;
            storage::Vector< std::string> substrings( util::SplitString( itr->first, " "));
            BCL_MessageStd
            (
              itr->first + " " + util::Format()( atom_type_to_index[ substrings( 0)])
              + " " + util::Format()( atom_type_to_index[ substrings( 1)])
              + " " + util::Format()( weight)
            );
            bond_length_incidences( counter, atom_type_to_index[ substrings( 0)]) += weight;
            bond_length_incidences( counter, atom_type_to_index[ substrings( 1)]) += weight;
          }
          BCL_Debug( bond_length_incidences);

          // calculate linear least squares fit
          const linal::Vector< double> orbital_sizes
          (
            math::LinearLeastSquares::SolutionAndChiSquared( bond_length_incidences, bond_length_aves).First()
          );
          counter = 0;
          for
          (
            storage::Set< std::string>::const_iterator itr( atom_types.Begin()), itr_end( atom_types.End());
            itr != itr_end;
            ++itr, ++counter
          )
          {
            const double radius( orbital_sizes( counter));
            const std::string element_type( itr->substr( 0, itr->size() - 2));
            const std::string n_bonds_electrons( itr->substr( itr->size() - 2, 2));
            if( radius > 0.2 && radius < 2.4)
            {
              OUTPUT << "SetPropertyInfo( phenotype_to_atom_type[ \"" << element_type
                     << n_bonds_electrons << "\"], "
                     << bond_property_name << ", "
                     << util::Format().FFP( 2)( radius) << ");\n";
              double &max_radius( atom_type_to_max_value[ element_type]);
              max_radius = std::max( max_radius, radius);
            }
            else
            {
              BCL_MessageStd
              (
                "Discarding incorrect value for the covalent radius of " + element_type + " w/ "
                + n_bonds_electrons.substr( 0, 1) + " bonds, " + n_bonds_electrons.substr( 1) + " electrons in bonds" +
                ( ring_size_minus_two ? " in ring of size " + util::Format()( ring_size_minus_two + 2) : " in chain ")
                + " for bond type " + util::Format()( bond_type_id)
              );
            }
          }
        }
      }

      storage::Map< std::string, double> vdw_radii_computed
      (
        chemistry::BondLengths::ComputeVdwRadii( m_VdwBondLengths, atom_type_to_max_value)
      );

      for
      (
        auto itr( vdw_radii_computed.Begin()), itr_end( vdw_radii_computed.End());
        itr != itr_end;
        ++itr
      )
      {
        // output the proper initializer
        OUTPUT << "SetPropertyInfo( phenotype_to_atom_type[ \"" << itr->first << "\"], "
               << "AtomTypeData::e_VdWaalsRadiusCSD, "
               << util::Format().FFP( 2)( itr->second) << ");\n";
      }
    }

  } // namespace app
} // namespace bcl

