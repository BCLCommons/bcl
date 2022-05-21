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
#include "chemistry/bcl_chemistry_fragment_probability_score.h"
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "chemistry/bcl_chemistry_small_molecule_fragment_isomorphism.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_split_rings.h"
#include "descriptor/bcl_descriptor_cheminfo_properties.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_statistics.h"

using bcl::linal::Distance;

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SmallMoleculeFragmentIsomorphism::s_Instance
    (
      GetObjectInstances().AddInstance( new SmallMoleculeFragmentIsomorphism())
    );

    //! @brief get the pseudocount that is used for all possible rotamers
    //! @return the pseudocount that is used for all possible rotamers
    double &SmallMoleculeFragmentIsomorphism::GetPseudocount()
    {
      static double s_pseudocount( 4.0);
      return s_pseudocount;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SmallMoleculeFragmentIsomorphism::SmallMoleculeFragmentIsomorphism()
    :
      m_Fragment(),
      m_FragmentCSI(),
      m_ContainsIncompleteRings( false)
    {
    }

    namespace
    {
      //! @brief get the number of bound hydrogens for each atom
      storage::Vector< size_t> GetNumberHydrogensExceptAtRingsAndTerminalAtoms( const FragmentComplete &A)
      {
        storage::Vector< size_t> n_h( A.GetSize(), size_t( 0));
        for( auto itr( A.GetAtomsIterator()); itr.NotAtEnd(); ++itr)
        {
          if( itr->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1)))
          {
            // atom in ring, undefined
            n_h( itr.GetPosition()) = util::GetUndefinedSize_t();
            continue;
          }
          const size_t real_nh( itr->GetNumberCovalentlyBoundHydrogens());
          if( itr->GetBonds().GetSize() <= real_nh + 1)
          {
            // terminal atom with hydrogens
            n_h( itr.GetPosition()) = util::GetUndefinedSize_t();
            continue;
          }
          n_h( itr.GetPosition()) = real_nh;
        }
        return n_h;
      }
    }

    //! @brief Constructor from data members
    //! @param FRAGMENT fragment whose rotamer information this class stores
    //! @param FRAGMENT_CSI isomorphism object of fragment and molecule
    SmallMoleculeFragmentIsomorphism::SmallMoleculeFragmentIsomorphism
    (
      const FragmentComplete &MOLECULE,
      const FragmentComplete &FRAGMENT,
      const graph::SubgraphIsomorphism< size_t, size_t> &FRAGMENT_CSI,
      bool REMOVE_RING_ROTAMER_INFO_UNLESS_FRAGMENT_IS_JUST_A_RING
    ) :
      m_Molecule( MOLECULE.Clone()),
      m_Fragment( FRAGMENT),
      m_FragmentCSI( FRAGMENT_CSI.GetIsomorphisms()),
      m_ContainsIncompleteRings( m_FragmentCSI.GetSize(), int( 0)),
      m_ContainsIncompleteHydrogenation( m_FragmentCSI.GetSize(), int( 0))
    {
      // get stored properties of the fragment
      RetrieveStoredProperties( REMOVE_RING_ROTAMER_INFO_UNLESS_FRAGMENT_IS_JUST_A_RING);
      // test for incomplete rings
      if( m_ContainsRings)
      {
        graph::SubgraphIsomorphism< util::SiPtr< const AtomConformationalInterface>, size_t> subg_iso_atoms;
        auto graph( ConformationGraphConverter::CreateGraphWithAtoms( MOLECULE));
        auto subgraph( ConformationGraphConverter::CreateGraphWithAtoms( FRAGMENT));
        subg_iso_atoms.SetGraphExternalOwnership( graph);
        subg_iso_atoms.SetSubgraphExternalOwnership( subgraph);
        subg_iso_atoms.SetIsomorphisms( FRAGMENT_CSI.GetIsomorphisms());
        auto subgiso( subg_iso_atoms.GetSubgraphIsomorphisms());
        for( size_t iso( 0), niso( FRAGMENT_CSI.GetIsomorphisms().GetSize()); iso < niso; ++iso)
        {
          auto adjacent_edges( subgiso( iso).GetAdjacentEdgeIndices());
          auto interior_edges( subgiso( iso).GetEdgeIndices());
          linal::Vector< size_t> n_ring_bonds_in_iso( MOLECULE.GetSize(), size_t( 0));
          for
          (
            auto itr_edges( interior_edges.Begin()), itr_edges_end( interior_edges.End());
            itr_edges != itr_edges_end;
            ++itr_edges
          )
          {
            if( ConfigurationalBondType( graph.GetEdgeData( itr_edges->First(), itr_edges->Second()))->IsBondInRing())
            {
              n_ring_bonds_in_iso( itr_edges->First()) += 1;
              n_ring_bonds_in_iso( itr_edges->Second()) += 1;
            }
          }
          for
          (
            auto itr_adj_edges( adjacent_edges.Begin()), itr_adj_edges_end( adjacent_edges.End());
            itr_adj_edges != itr_adj_edges_end;
            ++itr_adj_edges
          )
          {
            if
            (
              ( n_ring_bonds_in_iso( itr_adj_edges->First()) || n_ring_bonds_in_iso( itr_adj_edges->Second()))
              && ConfigurationalBondType( graph.GetEdgeData( itr_adj_edges->First(), itr_adj_edges->Second()))->IsBondInRing()
            )
            {
              m_ContainsIncompleteRings( iso) = 1;
              break;
            }
          }
        }
      }

      FragmentSplitRings ringsplit( true, 2);
      auto ring_systems( ringsplit( FRAGMENT));
      m_NumberRingSystems = ring_systems.GetSize();
      m_NumberAromaticRingSystems = size_t( 0);
      for( auto itr( ring_systems.Begin()), itr_end( ring_systems.End()); itr != itr_end; ++itr)
      {
        if( itr->CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsAromatic, 1))
        {
          ++m_NumberAromaticRingSystems;
        }
      }
      FragmentSplitRings chains( false, 1);
      m_NumberChainSystems = chains( FRAGMENT).GetSize();

      if( FRAGMENT.GetSize() > size_t( 4))
      {
        auto nh_molecule( GetNumberHydrogensExceptAtRingsAndTerminalAtoms( MOLECULE));
        auto nh_fragment( GetNumberHydrogensExceptAtRingsAndTerminalAtoms( m_Fragment));
        for( size_t iso( 0), niso( m_FragmentCSI.GetSize()); iso < niso; ++iso)
        {
          const auto &iso_ref( m_FragmentCSI( iso));
          for( size_t index_frag( 0), sz_frag( nh_fragment.GetSize()); index_frag < sz_frag; ++index_frag)
          {
            if( nh_fragment( index_frag) < nh_molecule( iso_ref( index_frag)))
            {
              m_ContainsIncompleteHydrogenation( iso) = 1;
              break;
            }
          }
        }
      }
      m_AllIncompleteHydrogenation = size_t( m_ContainsIncompleteHydrogenation.Sum()) == m_ContainsIncompleteHydrogenation.GetSize();
    }

    //! @brief Clone function
    //! @return pointer to new SmallMoleculeFragmentIsomorphism
    SmallMoleculeFragmentIsomorphism *SmallMoleculeFragmentIsomorphism::Clone() const
    {
      return new SmallMoleculeFragmentIsomorphism( *this);
    }

    //! @brief get the free energy of a particular rotamer
    //! @param ROTAMER rotamer signature of interest
    //! @param CONSIDER_ISOMETRY_CHANGES true - consider changes in isometry
    double SmallMoleculeFragmentIsomorphism::GetRotamerFreeEnergy
    (
      const linal::Vector< int> &ROTAMER,
      const bool &CONSIDER_ISOMETRY_CHANGES
    ) const
    {
      auto find_index( m_RotamerIndexMap.Find( ROTAMER));
      if( find_index == m_RotamerIndexMap.End())
      {
        return CONSIDER_ISOMETRY_CHANGES
             ? ( m_FreeEnergyUnseenRotamer - m_ExpectedFreeEnergy.GetAverage())
             : ( m_FreeEnergyUnseenRotamerNoIsometry - m_ExpectedFreeEnergyNoIsometry.GetAverage());
      }
      double rotamer_energy
      (
        CONSIDER_ISOMETRY_CHANGES
        ? m_RotamerFreeEnergy( find_index->second) - m_ExpectedFreeEnergy.GetAverage()
        : m_RotamerFreeEnergyIsometryChangesForbidden( find_index->second) - m_ExpectedFreeEnergyNoIsometry.GetAverage()
      );
      return rotamer_energy;
    }

    //! @brief get the free energy of a particular rotamer relative to others with the same dihedral bin
    //! @param ROTAMER rotamer signature of interest
    //! @param CONSIDER_ISOMETRY_CHANGES true - consider changes in isometry
    double SmallMoleculeFragmentIsomorphism::GetRotamerDihedralFineEnergy
    (
      const linal::Vector< double> &ROTAMER_ANGLES,
      const linal::Vector< int> &ROTAMER
    ) const
    {
      auto find_index( m_RotamerIndexMap.Find( ROTAMER));
      if( find_index == m_RotamerIndexMap.End())
      {
        return 0.0;
      }
      return linal::Distance( ROTAMER_ANGLES, m_Rotamers( find_index->second)) / 15.0 / math::Sqrt( double( ROTAMER.GetSize())) - 1.0;
    }

    //! @brief get all stored properties or rotamer library data of the fragment of interest
    //! @return rotamer library data of the fragment
    void SmallMoleculeFragmentIsomorphism::RetrieveStoredProperties
    (
      const bool &REMOVE_RING_ROTAMER_INFO_UNLESS_FRAGMENT_IS_JUST_A_RING
    )
    {
      m_BinSize = double( 30.0);
      // get total number of rotatable bonds in the fragment of interest
      m_RotatableBondCount = descriptor::GetCheminfoProperties().calc_NRotBond->SumOverObject( m_Fragment)( 0);
      m_RingBondCount = m_Fragment.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, 1);
      m_ChainBondCount = m_Fragment.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, 0);
      if( m_RingBondCount < size_t( 3))
      {
        // handle dihedral fragments with partial rings
        m_ChainBondCount += m_RingBondCount;
        m_RingBondCount = 0;
      }
      // this check is here due to some fragment library corruption that occurred due to a bug in RingWithUnsaturatedSubstituents
      if( m_RingBondCount)
      {
        for( auto itr( m_Fragment.GetAtomsIterator()); itr.NotAtEnd(); ++itr)
        {
          if
          (
            itr->GetBonds().GetSize() == size_t( 1)
            && itr->GetBonds()( 0).GetBondType()->GetBondData( ConfigurationalBondTypeData::e_IsInRing)
          )
          {
            BCL_MessageStd( "Messed up fragment (fragmentation along ring bonds)");
            m_Fragment.WriteMDL( util::GetLogger());
            m_DihedralChainBondCount = 0;
            m_RingBondCount = m_ChainBondCount = 0;
            m_RotatableBondCount = 0;
            m_Fragment = FragmentComplete();
            m_FragmentCSI.Reset();
            return;
          }
        }
      }

      m_DihedralBonds = PriorityDihedralAngles::GetDihedralEdges( m_Fragment);
      m_DihedralChainBondCount = m_DihedralBonds.GetSize() - m_RingBondCount;

      // exocyclic bonds not involved in dihedrals cause some issues with generate_3D, so for now we just remove them fromst
      // the fragment at this step
      // Use REMOVE_RING_ROTAMER_INFO_UNLESS_FRAGMENT_IS_JUST_A_RING as a proxy for whether it's okay to modify the molecule
      // and properties and such
      bool remove_exocyclic_atoms
      (
        REMOVE_RING_ROTAMER_INFO_UNLESS_FRAGMENT_IS_JUST_A_RING
        && !m_DihedralChainBondCount && m_RingBondCount && m_ChainBondCount
      );
      const size_t original_frag_n_atoms( m_Fragment.GetSize());
      storage::Vector< size_t> atoms_to_keep, inverse_atoms_to_keep;
      if( remove_exocyclic_atoms)
      {
        atoms_to_keep.AllocateMemory( original_frag_n_atoms);
        inverse_atoms_to_keep.Resize( original_frag_n_atoms, util::GetUndefined< size_t>());
        for( size_t a( 0); a < original_frag_n_atoms; ++a)
        {
          if( m_Fragment.GetAtomVector()( a).CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, 1))
          {
            inverse_atoms_to_keep( a) = atoms_to_keep.GetSize();
            atoms_to_keep.PushBack( a);
          }
        }
        auto atm_vec_cpy( m_Fragment.GetAtomVector());
        atm_vec_cpy.Reorder( atoms_to_keep);
        // update the fragment
        m_Fragment = FragmentComplete( atm_vec_cpy, m_Fragment.GetName(), m_Fragment.GetStoredProperties().GetMDLProperties());
        m_DihedralBonds = PriorityDihedralAngles::GetDihedralEdges( m_Fragment);
        m_ChainBondCount = 0;
        // update the isomorphisms
        for( auto itr_iso( m_FragmentCSI.Begin()), itr_iso_end( m_FragmentCSI.End()); itr_iso != itr_iso_end; ++itr_iso)
        {
          itr_iso->Reorder( atoms_to_keep);
        }
      }

      m_DihedralAngleIndices = PriorityDihedralAngles()( m_Fragment).Second();

      const bool remove_ring_rotamer_info
      (
        REMOVE_RING_ROTAMER_INFO_UNLESS_FRAGMENT_IS_JUST_A_RING && m_DihedralChainBondCount && m_RingBondCount
      );
      storage::Vector< size_t> dihedrals_to_keep;
      if( remove_ring_rotamer_info)
      {
        dihedrals_to_keep.AllocateMemory( m_DihedralBonds.GetSize() - m_RingBondCount);
        for( size_t i( 0), n_dihedral( m_DihedralBonds.GetSize()); i < n_dihedral; ++i)
        {
          if( !m_DihedralBonds( i).GetEdgeData()->IsBondInRing())
          {
            dihedrals_to_keep.PushBack( i);
          }
        }
        m_DihedralBonds.Reorder( dihedrals_to_keep);
        m_DihedralAngleIndices.Reorder( dihedrals_to_keep);
        m_DihedralChainBondCount = dihedrals_to_keep.GetSize();
      }

      m_ContainsRings = !remove_ring_rotamer_info
                        && m_Fragment.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, 1) > size_t( 2);
      m_ContainsRingConformations = m_ContainsRings
                                    && !remove_ring_rotamer_info
                                    && !m_Fragment.GetStoredProperties().GetMDLProperty( "Rotamer1Coordinates").empty();

      // get number of rotamers of the fragment in question
      linal::Vector< float> rotamer_counts( m_Fragment.GetStoredProperties().GetMDLPropertyAsVector( "RotamerCounts"));
      if( rotamer_counts.IsEmpty())
      {
        return;
      }
      m_RotamerCounts = storage::Vector< double>( rotamer_counts.Begin(), rotamer_counts.End());
      m_TotalRotamers = m_RotamerCounts.GetSize();
      double implied_frag_counts( 0);
      double pseudocount( GetPseudocount());
      for( size_t i( 0); i < m_TotalRotamers; ++i)
      {
        m_RotamerCounts( i) = std::max( m_RotamerCounts( i), double( 1.0));
        implied_frag_counts += m_RotamerCounts( i);
        pseudocount = std::max( std::min( pseudocount, m_RotamerCounts( i)), 1.0);
      }

      m_SimulatedCount = 0.0;
      if( !m_Fragment.GetMDLProperty( "NovelRotamersSimulated").empty())
      {
        m_SimulatedCount = m_FractionSimulatedNotObserved = m_Fragment.GetMDLPropertyAsVector( "NovelRotamersSimulated")( 0);
        m_RotamerFractionSimulated = m_Fragment.GetMDLPropertyAsVector( "RotamerCountsSimulated");
        m_SimulatedCount += m_RotamerFractionSimulated.Sum();
        m_FractionSimulatedNotObserved /= m_SimulatedCount;
      }
      else
      {
        m_RotamerFractionSimulated = linal::Vector< float>( m_TotalRotamers, 0.0);
        m_FractionSimulatedNotObserved = 0.0;
      }

      // get total number of times a fragment is seen in a structure database
      m_TotalFragmentCounts = std::max
                              (
                                util::ConvertStringToNumericalValue< double>
                                (
                                  m_Fragment.GetStoredProperties().GetMDLProperty( "ScaffoldCount")
                                ),
                                implied_frag_counts
                              );

      linal::Vector< double> average( m_Fragment.GetStoredProperties().GetMDLPropertyAsVector( "Average"));
      linal::Vector< int> bins( m_Fragment.GetStoredProperties().GetMDLPropertyAsVector( "Bins"));
      size_t number_dihedrals( bins.GetSize() / m_TotalRotamers);

      linal::Vector< int>::const_iterator itr_bin( bins.Begin());
      linal::Vector< double>::const_iterator itr_average( average.Begin());

      // for each rotamer get rotamer counts, average angle, rotamer bin signature, coordinates
      for( size_t rotamer_count( 0); rotamer_count < m_TotalRotamers; ++rotamer_count)
      {
        linal::VectorConstReference< int> cur_bin( number_dihedrals, itr_bin);
        linal::VectorConstReference< double> cur_average( number_dihedrals, itr_average);
        m_RotamerBins.PushBack( cur_bin);
        m_Rotamers.PushBack( cur_average);
        std::advance( itr_bin, number_dihedrals);
        std::advance( itr_average, number_dihedrals);
      }

      if( remove_ring_rotamer_info)
      {
        size_t n_to_keep( dihedrals_to_keep.GetSize());
        for( size_t rotamer_count( 0); rotamer_count < m_TotalRotamers; ++rotamer_count)
        {
          auto &ref_bins( m_RotamerBins( rotamer_count));
          auto &ref_ave( m_Rotamers( rotamer_count));
          for( size_t i( 0); i < n_to_keep; ++i)
          {
            ref_bins( i) = ref_bins( dihedrals_to_keep( i));
            ref_ave( i) = ref_ave( dihedrals_to_keep( i));
          }
          ref_bins.Shrink( n_to_keep);
          ref_ave.Shrink( n_to_keep);
        }
        // merge rotamers that differed only by their ring bins
        for( size_t rotamer_count( m_TotalRotamers - 1); rotamer_count < m_TotalRotamers; --rotamer_count)
        {
          for( size_t rotamer_count_b( 0); rotamer_count_b < rotamer_count; ++rotamer_count_b)
          {
            if( m_RotamerBins( rotamer_count) == m_RotamerBins( rotamer_count_b))
            {
              // compute the new average rotamer
              math::RunningAverage< linal::Vector< double> > rotamer_aves;
              rotamer_aves.AddWeightedObservation( m_Rotamers( rotamer_count), m_RotamerCounts( rotamer_count));
              rotamer_aves.AddWeightedObservation( m_Rotamers( rotamer_count_b), m_RotamerCounts( rotamer_count_b));
              m_Rotamers( rotamer_count_b) = rotamer_aves.GetAverage();
              // update counts
              m_RotamerCounts( rotamer_count_b) += m_RotamerCounts( rotamer_count);
              m_RotamerBins.RemoveElements( rotamer_count, 1);
              m_RotamerCounts.RemoveElements( rotamer_count, 1);
              m_Rotamers.RemoveElements( rotamer_count, 1);
              std::copy
              (
                m_RotamerFractionSimulated.Begin() + rotamer_count + 1,
                m_RotamerFractionSimulated.End(),
                m_RotamerFractionSimulated.Begin() + rotamer_count
              );
              m_RotamerFractionSimulated.Shrink( m_RotamerFractionSimulated.GetSize() - size_t( 1));
              --m_TotalRotamers;
              break;
            }
          }
        }
      }

      if( m_SimulatedCount)
      {
        m_RotamerFractionSimulated += float( pseudocount);
        m_SimulatedCount += pseudocount * m_RotamerFractionSimulated.GetSize();
        m_RotamerFractionSimulated /= float( m_SimulatedCount);
      }
      // get rotamer library data if total number of rotamers is greater than 1
      if( !remove_ring_rotamer_info)
      {
        if( m_TotalRotamers > 0 && m_ContainsRingConformations)
        {
          // for each rotamer get rotamer counts, average angle, rotamer bin signature, coordinates
          for( size_t rotamer_count( 1); rotamer_count <= m_TotalRotamers; ++rotamer_count)
          {
            // rotamer coordinates for a particular rotamer
            std::string coordinate_string
            (
              m_Fragment.GetStoredProperties().GetMDLProperty
              (
                "Rotamer" + util::Format()( rotamer_count) + "Coordinates"
              )
            );

            // get coordinates of the cluster center of rotamers of the fragment
            const storage::Vector< std::string> result
            (
              util::SplitString( util::TrimString( coordinate_string), " \t\n\r,")
            );
            storage::Vector< linal::Vector3D> rotamer_coordinates;

            for( size_t atom_index( 0); atom_index < original_frag_n_atoms; ++atom_index)
            {
              size_t vector_index( 16 * atom_index);
              rotamer_coordinates.PushBack
              (
                linal::Vector3D
                (
                  util::ConvertStringToNumericalValue< double>( result( vector_index)),
                  util::ConvertStringToNumericalValue< double>( result( vector_index + 1)),
                  util::ConvertStringToNumericalValue< double>( result( vector_index + 2))
                )
              );
            }
            if( remove_exocyclic_atoms)
            {
              rotamer_coordinates.Reorder( atoms_to_keep);
            }

            // store the information in the member variables
            m_RotamerCoordinates.PushBack( rotamer_coordinates);
          }
        }
        else if( m_TotalRotamers == size_t( 1))
        {
          m_RotamerCoordinates.PushBack();
          for( size_t atom_index( 0); atom_index < m_Fragment.GetNumberAtoms(); ++atom_index)
          {
            m_RotamerCoordinates( 0).PushBack( m_Fragment.GetAtomVector()( atom_index).GetPosition());
          }
        }
        else if( m_TotalRotamers > size_t( 1) && !m_DihedralChainBondCount)
        {
          BCL_MessageCrt( "Missing ring conformations, error in fragment library!");
          m_Fragment.WriteMDL( util::GetLogger());
          m_DihedralChainBondCount = 0;
          m_RingBondCount = m_ChainBondCount = 0;
          m_RotatableBondCount = 0;
          m_Fragment = FragmentComplete();
          m_FragmentCSI.Reset();
          return;
        }
      }

      m_MaxRotamerCounts = math::Statistics::MaximumValue( m_RotamerCounts.Begin(), m_RotamerCounts.End());

      m_TheoreticalMaxRotamersNoIsometry = ( m_RotamerCoordinates.IsEmpty() ? 1 : m_TotalRotamers);
      m_TheoreticalMaxRotamersNoIsometry *= math::Pow( size_t( 360.0 / m_BinSize), m_RotatableBondCount);
      m_TheoreticalMaxRotamers = m_TheoreticalMaxRotamersNoIsometry;
      // amide-bonds are essentially like double bonds, but do not formally count as an isometry change for whatever reason
      double isometric_multiplier( 1);
      for( auto itr_atoms( m_Fragment.GetAtomsIterator()); itr_atoms.NotAtEnd(); ++itr_atoms)
      {
        if( itr_atoms->GetBonds().GetSize() <= size_t( 1))
        {
          continue;
        }
        for
        (
          auto itr_bonds( itr_atoms->GetBonds().Begin()), itr_bonds_end( itr_atoms->GetBonds().End());
          itr_bonds != itr_bonds_end;
          ++itr_bonds
        )
        {
          if( &itr_bonds->GetTargetAtom() < &*itr_atoms || itr_bonds->GetBondType()->IsBondInRing())
          {
            continue;
          }
          if( itr_bonds->GetTargetAtom().GetBonds().GetSize() <= size_t( 1))
          {
            continue;
          }
          if( itr_bonds->GetBondType()->GetConjugation() == ConstitutionalBondTypeData::e_Amide)
          {
            m_TheoreticalMaxRotamersNoIsometry *= 2.0;
            m_TheoreticalMaxRotamers *= 2.0;
          }
          if
          (
            itr_bonds->GetBondType()->GetBondData( ConfigurationalBondTypeData::e_IsIsometric)
          )
          {
            m_TheoreticalMaxRotamers *= 2.0;
            isometric_multiplier *= 2.0;
          }
        }
      }

      m_TheoreticalMaxRotamers = std::max( m_TheoreticalMaxRotamers, double( m_TotalRotamers));
      m_TheoreticalMaxRotamersNoIsometry = std::max( m_TheoreticalMaxRotamersNoIsometry, double( m_TotalRotamers));
      m_RotamerFreeEnergy.Resize( m_TotalRotamers);
      m_RotamerFreeEnergyIsometryChangesForbidden.Resize( m_TotalRotamers);

      m_ExpectedFreeEnergy.Reset();
      m_ExpectedFreeEnergyNoIsometry.Reset();

//      if( m_SimulatedCount)
//      {
//        BCL_Debug( m_SimulatedCount)
//        BCL_Debug( m_RotamerFractionSimulated.Sum())
//        BCL_Debug( m_RotamerCounts)
//        BCL_Debug( m_RotamerFractionSimulated)
//        BCL_Debug( n_uncovered);
//      }

      if( 1) /// m_RotamerFractionSimulated.Sum() < 0.8 || m_Fragment.HasChiralCenters() || m_Fragment.HasIsometry())
      {
        const double average_count_per_rotamer( implied_frag_counts / m_TheoreticalMaxRotamers + GetPseudocount());
        const double average_count_per_rotamer_no_isometry( implied_frag_counts / m_TheoreticalMaxRotamersNoIsometry + GetPseudocount());

        m_FreeEnergyUnseenRotamer = -log( GetPseudocount() / average_count_per_rotamer);
        m_FreeEnergyUnseenRotamerNoIsometry = -log( GetPseudocount() / average_count_per_rotamer_no_isometry);
        m_FreeEnergyBestRotamer = -log( ( GetPseudocount() + m_MaxRotamerCounts) / average_count_per_rotamer);
        m_FreeEnergyBestRotamerNoIsometry = -log( ( GetPseudocount() + m_MaxRotamerCounts) / average_count_per_rotamer_no_isometry);

        for( size_t rotamer_n( 0); rotamer_n < m_TotalRotamers; ++rotamer_n)
        {
          m_RotamerFreeEnergy( rotamer_n) =
            -log( ( m_RotamerCounts( rotamer_n) + GetPseudocount()) / average_count_per_rotamer);
          m_RotamerFreeEnergyIsometryChangesForbidden( rotamer_n) =
            -log( ( m_RotamerCounts( rotamer_n) + GetPseudocount()) / average_count_per_rotamer_no_isometry);
          m_ExpectedFreeEnergy.AddWeightedObservation( m_RotamerFreeEnergy( rotamer_n), m_RotamerCounts( rotamer_n));
          m_ExpectedFreeEnergyNoIsometry.AddWeightedObservation
          (
            m_RotamerFreeEnergyIsometryChangesForbidden( rotamer_n),
            m_RotamerCounts( rotamer_n)
          );
        }
      }
      else
      {
        // use of simulated counts *should* in principle be better, however, our current simulation counts are lacking
        // Probably need to redo them properly with exhaustive angle-based sampling. This code kept in here just to test
        // this experimental feature easily.
        const double average_sim_rat_per_rotamer_raw( m_RotamerFractionSimulated.Sum() / double( m_RotamerFractionSimulated.GetSize()));
        const double average_sim_count_per_rotamer_raw( m_SimulatedCount * average_sim_rat_per_rotamer_raw);
        const double average_sim_count_per_rotamer( m_SimulatedCount * m_RotamerFractionSimulated.Sum() / ( m_RotamerFractionSimulated.GetSize() + m_FractionSimulatedNotObserved / average_sim_count_per_rotamer_raw));
        const double average_sim_count_per_rotamer_no_isometry( average_sim_count_per_rotamer * isometric_multiplier);
        const double average_sim_rat_per_rotamer( average_sim_count_per_rotamer / m_SimulatedCount);
        const double average_sim_rat_per_rotamer_no_isometry( average_sim_rat_per_rotamer * isometric_multiplier);
//        BCL_Debug( average_sim_rat_per_rotamer_raw);
//        BCL_Debug( average_sim_count_per_rotamer_raw);
//        BCL_Debug( average_sim_count_per_rotamer);
//        BCL_Debug( average_sim_count_per_rotamer_no_isometry);
//        BCL_Debug( average_sim_rat_per_rotamer);
//        BCL_Debug( average_sim_rat_per_rotamer_no_isometry);

        m_FreeEnergyUnseenRotamer += -log( GetPseudocount() / average_sim_count_per_rotamer);
        m_FreeEnergyUnseenRotamerNoIsometry += -log( GetPseudocount() / average_sim_count_per_rotamer_no_isometry);
        double rotcntsum( math::Statistics::Sum( m_RotamerCounts.Begin(), m_RotamerCounts.End()));
        double rotcntave( GetPseudocount() + rotcntsum / double( m_RotamerCounts.GetSize()));
        for( size_t rotamer_n( 0); rotamer_n < m_TotalRotamers; ++rotamer_n)
        {
          m_RotamerFreeEnergy( rotamer_n) =
              -log( ( m_RotamerCounts( rotamer_n) + GetPseudocount()) / rotcntave / ( m_RotamerFractionSimulated( rotamer_n) / average_sim_rat_per_rotamer));
          m_RotamerFreeEnergyIsometryChangesForbidden( rotamer_n) =
              -log( ( m_RotamerCounts( rotamer_n) + GetPseudocount()) / rotcntave / ( m_RotamerFractionSimulated( rotamer_n) / average_sim_rat_per_rotamer_no_isometry));
          m_ExpectedFreeEnergy.AddWeightedObservation( m_RotamerFreeEnergy( rotamer_n), m_RotamerCounts( rotamer_n));
          m_ExpectedFreeEnergyNoIsometry.AddWeightedObservation
          (
            m_RotamerFreeEnergyIsometryChangesForbidden( rotamer_n),
            m_RotamerCounts( rotamer_n)
          );
        }
//        BCL_Debug( m_RotamerFreeEnergy);
//        BCL_Debug( m_RotamerCounts);
//        BCL_Debug( m_RotamerFractionSimulated);
//        BCL_Debug( m_FreeEnergyUnseenRotamer);
//        BCL_Debug( m_FreeEnergyUnseenRotamerNoIsometry);
      }
      m_ExpectedFreeEnergyNoIsometry.AddWeightedObservation( m_FreeEnergyUnseenRotamerNoIsometry, GetPseudocount());
      m_ExpectedFreeEnergy.AddWeightedObservation( m_FreeEnergyUnseenRotamer, GetPseudocount());

      if( m_ContainsRings && !m_ChainBondCount && !m_Molecule->HasBadGeometry())
      {
        FragmentProbabilityScore scorer;
        typename FragmentProbabilityScore::t_Map bondmap;

        // all isomorphisms will yield the same central bonds so take the first one
        auto itr_isomorphisms( m_FragmentCSI.Begin());

        storage::Map< storage::Set< size_t>, std::pair< bool, storage::Vector< size_t> > > have_seen_conformation;
        // store priority dihedral angle information of fragment in terms of fragment isomorphism with the molecule of interest
        for
        (
          auto itr_isomorphisms_end( m_FragmentCSI.End());
          itr_isomorphisms != itr_isomorphisms_end;
          ++itr_isomorphisms
        )
        {
          // get isomorphism between fragment and the molecule of interest
          const storage::Vector< size_t> &isomorphism( *itr_isomorphisms);

          // insert placeholder into the map from atom indices to whether conformation was already seen
          auto itr_ins
          (
            have_seen_conformation.Insert
            (
              std::make_pair( storage::Set< size_t>( isomorphism.Begin(), isomorphism.End()), std::make_pair( false, *itr_isomorphisms))
            )
          );
          // skip if conformation already seen
          if( itr_ins.first->second.first)
          {
            continue;
          }

          // create a vector to store  dihedral bond atom indices of molecule ( this is isomorphic to fragment dihedral)
          storage::Vector< storage::VectorND< 4, size_t> > molecule_dihedral;
          molecule_dihedral.AllocateMemory( m_DihedralAngleIndices.GetSize());
          // get atom indices from the frame of reference of the molecule of interest
          for
          (
            storage::Vector< storage::VectorND< 4, size_t> >::const_iterator
              itr( m_DihedralAngleIndices.Begin()), itr_end( m_DihedralAngleIndices.End());
            itr != itr_end;
            ++itr
          )
          {
            molecule_dihedral.PushBack
            (
              storage::VectorND< 4, size_t>
              (
                isomorphism( itr->First()),
                isomorphism( itr->Second()),
                isomorphism( itr->Third()),
                isomorphism( itr->Fourth())
              )
            );
          }

          auto bins_and_angles( scorer.GetMolecularBinFromFragment( *m_Molecule, molecule_dihedral, bondmap));
          if( m_RotamerBins.Find( bins_and_angles.first) < m_RotamerBins.GetSize())
          {
            itr_ins.first->second.first = true;
          }
        }
        double ave( math::Statistics::Mean( m_RotamerCounts.Begin(), m_RotamerCounts.End()));
        for
        (
          auto itr_map( have_seen_conformation.Begin()), itr_map_end( have_seen_conformation.End());
          itr_map != itr_map_end;
          ++itr_map
        )
        {
          // conformation was already seen
          if( itr_map->second.first)
          {
            continue;
          }
          // get isomorphism between fragment and the molecule of interest
          const storage::Vector< size_t> &isomorphism( itr_map->second.second);

          // create a vector to store  dihedral bond atom indices of molecule ( this is isomorphic to fragment dihedral)
          storage::Vector< storage::VectorND< 4, size_t> > molecule_dihedral;
          molecule_dihedral.AllocateMemory( m_DihedralAngleIndices.GetSize());
          // get atom indices from the frame of reference of the molecule of interest
          for
          (
            storage::Vector< storage::VectorND< 4, size_t> >::const_iterator
              itr( m_DihedralAngleIndices.Begin()), itr_end( m_DihedralAngleIndices.End());
            itr != itr_end;
            ++itr
          )
          {
            molecule_dihedral.PushBack
            (
              storage::VectorND< 4, size_t>
              (
                isomorphism( itr->First()),
                isomorphism( itr->Second()),
                isomorphism( itr->Third()),
                isomorphism( itr->Fourth())
              )
            );
          }
          typename FragmentProbabilityScore::t_Map bondmap;
          auto bins_and_angles( scorer.GetMolecularBinFromFragment( *m_Molecule, molecule_dihedral, bondmap));

          m_ContainsRingConformations = true;
          m_Rotamers.PushBack( bins_and_angles.second);
          m_RotamerBins.PushBack( bins_and_angles.first);
          m_RotamerFreeEnergy.PushBack( 0.01);
          m_RotamerFreeEnergyIsometryChangesForbidden.PushBack( 0.01);
          storage::Vector< linal::Vector3D> coords( m_Fragment.GetNumberAtoms());
          for( size_t atom_index( 0); atom_index < m_Fragment.GetNumberAtoms(); ++atom_index)
          {
            coords( atom_index) = m_Molecule->GetAtomVector()( isomorphism( atom_index)).GetPosition();
          }
          m_RotamerCoordinates.PushBack( coords);
          m_RotamerCounts.PushBack( ave * 0.99);
        }
      }

      for( size_t rotamer_n( 0); rotamer_n < m_TotalRotamers; ++rotamer_n)
      {
        m_RotamerIndexMap[ m_RotamerBins( rotamer_n)] = rotamer_n;
      }
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SmallMoleculeFragmentIsomorphism::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SmallMoleculeFragmentIsomorphism::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_Fragment, ISTREAM);
      io::Serialize::Read( m_FragmentCSI, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &SmallMoleculeFragmentIsomorphism::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      io::Serialize::Write( m_Fragment, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_FragmentCSI, OSTREAM, INDENT);

      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl

