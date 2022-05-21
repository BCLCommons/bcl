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
#include "chemistry/bcl_chemistry_rotamer_cluster_center.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "math/bcl_math_running_average_sd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief constructor
    //! @param MOLECULE the configuration of interest
    //! @param TOLERANCE the tolerance level to be used for dihedral binning
    RotamerClusterCenter::RotamerClusterCenter( const ConformationInterface &MOLECULE, double BIN_SIZE) :
      m_Molecule( new FragmentConfigurationShared( MOLECULE)),
      m_ConformationComparison( BIN_SIZE),
      m_BinSize( BIN_SIZE),
      m_NumberSimulatedButNotSeen( 0.0),
      m_IsPureRing
      (
        MOLECULE.GetNumberBonds()
        ==
        MOLECULE.CountNonValenceBondsWithProperty( ConfigurationalBondTypeData::e_IsInRing, size_t( 1))
      )
    {
    }

    //! @brief Clone function
    //! @return pointer to new RotamerClusterCenter
    RotamerClusterCenter *RotamerClusterCenter::Clone() const
    {
      return new RotamerClusterCenter( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RotamerClusterCenter::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the configuration whose conformation clusters this object calculates
    //! @return configuration whose conformation clusters this object calculates
    const FragmentConfigurationShared &RotamerClusterCenter::GetConfiguration() const
    {
      return *m_Molecule;
    }

    //! @brief returns reference to all conformation clusters
    //! @return reference to all conformation clusters
    const storage::List< RotamerEnsemble> &RotamerClusterCenter::GetClusters() const
    {
      return m_Clusters;
    }

    //! @brief returns reference to all conformation clusters
    //! @return reference to all conformation clusters
    math::RunningAverageSD< double> RotamerClusterCenter::GetStats() const
    {
      math::RunningAverageSD< double> stats;
      for
      (
        storage::List< RotamerEnsemble>::const_iterator itr( m_Clusters.Begin()), itr_end( m_Clusters.End());
        itr != itr_end;
        ++itr
      )
      {
        stats += itr->GetWeights();
      }
      return stats;
    }

    //! @brief returns reference to all conformation clusters
    //! @return reference to all conformation clusters
    double RotamerClusterCenter::GetMaxWeight() const
    {
      double max_weight( 0.0);
      for
      (
        storage::List< RotamerEnsemble>::const_iterator itr( m_Clusters.Begin()), itr_end( m_Clusters.End());
        itr != itr_end;
        ++itr
      )
      {
        max_weight = std::max( max_weight, itr->GetWeights());
      }
      return max_weight;
    }

    //! @brief get the count of rotamer seen most number of times
    //! @return count of rotamer seen most number of times
    double RotamerClusterCenter::GetSimulatedWeightUnseenRotamers() const
    {
      return m_NumberSimulatedButNotSeen;
    }

    //! @brief returns reference to all conformation clusters
    //! @return reference to all conformation clusters
    size_t RotamerClusterCenter::GetMaxInstance() const
    {
      size_t max_instances( 0);
      for
      (
        storage::List< RotamerEnsemble>::const_iterator itr( m_Clusters.Begin()), itr_end( m_Clusters.End());
        itr != itr_end;
        ++itr
      )
      {
        max_instances = std::max( max_instances, itr->GetNumberInstances());
      }
      return max_instances;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculates which bin signature the given conformation belongs to
    //! @param CONFORMATION conformation to be used for calculation
    //! @return void
    void RotamerClusterCenter::operator()
    (
      const FragmentComplete &CONFORMATION,
      const double WEIGHT,
      const bool &SIMULATED
    )
    {
      storage::Vector< storage::VectorND< 4, size_t> > dihedral_atom_indices_storage;
      storage::Vector< int> dihedral_keys_storage;

      // get the dihedral angles for priority dihedral bonds and priority atoms for the conformation which is passed in
      storage::Pair< storage::Vector< double>, storage::Vector< storage::VectorND< 4, size_t> > >
      dihedral_angles_indices( m_PriorityDihedral( CONFORMATION));

      // get reference to dihedral bins
      storage::Vector< double> &dihedral_angles( dihedral_angles_indices.First());

      // get reference to priority atoms
      dihedral_atom_indices_storage = dihedral_angles_indices.Second();
      dihedral_keys_storage
        = m_ConformationComparison.DetermineDihedralKeys( dihedral_angles, true, true);

      const storage::Vector< storage::VectorND< 4, size_t> > &dihedral_atom_indices
      (
        dihedral_atom_indices_storage
      );

      const storage::Vector< int> &dihedral_keys
      (
        dihedral_keys_storage
      );

      // create graph colored by dihedral bins for conformation
      graph::ConstGraph< size_t, size_t> graph
      (
        CreateGraphWithDihedralBinsOnEdges( CONFORMATION, dihedral_keys, dihedral_atom_indices)
      );

      // create isomorphism object and set the sub graph
      graph::SubgraphIsomorphism< size_t, size_t> isomorphism;
      isomorphism.SetSubgraphExternalOwnership( graph);

      // get atom vector of the conformation which is passed in
      AtomVector< AtomConformationalShared> atoms_reordered( CONFORMATION.GetAtomInfo(), CONFORMATION.GetBondInfo());

      // go through each cluster information object and get graphs associated with each to determine whether the given
      // conformation is new or whether it belongs to an existing cluster
      bool found_match( false);
      for
      (
        storage::List< RotamerEnsemble>::iterator itr( m_Clusters.Begin()), itr_end( m_Clusters.End());
        itr != itr_end;
        ++itr
      )
      {
        // set the graph from cluster information as the graph in isomorphism object
        isomorphism.SetGraphExternalOwnership( itr->GetGraph());

        // if isomorphism is found, conformation belongs to an existing cluster. So add to it
        if( isomorphism.FindIsomorphism())
        {
          BCL_Assert( !found_match, "Should have stopped after the first match!");
          // reorder the molecule according to this isomorphism and add it to the cluster
          found_match = true;

          storage::Vector< size_t> inverse_mapping( isomorphism.GetInverseIsomorphism());
          atoms_reordered.Reorder( inverse_mapping);
          itr->operator ()( FragmentConformationShared( m_Molecule, atoms_reordered), WEIGHT, SIMULATED);
          break;
        }
      }
      // if no isomorphism is found then create a new cluster information object
      if( !found_match)
      {
        if( !SIMULATED)
        {
          m_Clusters.PushBack
          (
            RotamerEnsemble( FragmentConformationShared( m_Molecule, atoms_reordered), graph, WEIGHT, m_BinSize, SIMULATED)
          );
        }
        else
        {
          m_NumberSimulatedButNotSeen += WEIGHT;
        }
      }
    }

    //! @brief calculates which bin signature the given conformation belongs to
    //! @param CONFORMATION conformation to be used for calculation
    //! @return void
    void RotamerClusterCenter::operator()
    (
      const graph::SubgraphIsomorphism< util::SiPtr< const AtomConformationalInterface>, size_t> &ATOM_ISO,
      const storage::Vector< storage::Vector< size_t> > &ISO,
      const std::string &NAME,
      const bool &SIMULATED
    )
    {
      const bool is_autoisomorphism( ATOM_ISO.GetGraph().GetSize() == ATOM_ISO.GetSubgraph().GetSize());
      storage::Vector< double> weights( GetIsomorphismWeights( ISO, ATOM_ISO.GetGraph().GetSize()));
      storage::Vector< storage::Vector< size_t> >::const_iterator itr_simple_iso( ISO.Begin());
      storage::Vector< double>::const_iterator itr_weight( weights.Begin());
      auto subgraph_isomorphisms( ATOM_ISO.GetSubgraphIsomorphisms());
      for
      (
        storage::Vector< graph::Subgraph< util::SiPtr< const AtomConformationalInterface>, size_t> >::const_iterator
          itr_subgraph( subgraph_isomorphisms.Begin()), itr_subgraph_end( subgraph_isomorphisms.End());
        itr_subgraph != itr_subgraph_end;
        ++itr_subgraph, ++itr_simple_iso, ++itr_weight
      )
      {
        //const storage::Vector< size_t> &cur_iso( *itr_simple_iso);

        // yep
        // convert the graph back into a small molecule
        FragmentComplete fragment
        (
          ConformationGraphConverter::CreateAtomsFromGraph( itr_subgraph->ToGraph(), !is_autoisomorphism),
          NAME
        );

        operator ()( fragment, *itr_weight, SIMULATED);
      }
    }

    //! @brief Create a graph for a conformation
    //! @param CONFORMATION a conformation of configuration of interest
    //! @param DIHEDRAL_BIN dihedral bin signature of the conformation being passed in
    //! @param ATOM_DIHEDRAL_INDICES the atoms that make the priority dihedral angles for the conformation being passed in
    //! @return The conformation converted into a graph with the coloring schemes dictated by dihedral bin signature
    graph::ConstGraph< size_t, size_t> RotamerClusterCenter::CreateGraphWithDihedralBinsOnEdges
    (
      const ConformationInterface &CONFORMATION,
      const storage::Vector< int> &DIHEDRAL_BINS,
      const storage::Vector< storage::VectorND< 4, size_t> > &ATOM_DIHEDRAL_INDICES
    ) const
    {
      // create a conformation graph for the conformation which is passed in
      graph::ConstGraph< size_t, size_t> graph
      (
        ConformationGraphConverter
        (
          ConformationGraphConverter::e_AtomTypeAndComplexRingChirality,
          ConfigurationalBondTypeData::e_BondOrderAmideOrAromaticWithRingness
        )
        ( CONFORMATION)
      );

      // compute the radix, which is the maximum # that any of the configurational bond types data can return
      const size_t radix( GetConstitutionalBondTypes().GetEnumCount());

      // get an iterator to dihedral bin vector for conformation of interest
      auto itr_bins( DIHEDRAL_BINS.Begin());

      for
      (
        storage::Vector< storage::VectorND< 4, size_t> >::const_iterator
          itr( ATOM_DIHEDRAL_INDICES.Begin()), itr_end( ATOM_DIHEDRAL_INDICES.End());
        itr != itr_end;
        ++itr, ++itr_bins
      )
      {
        // get atoms involved in the central bond
        const size_t &vertex_a( itr->Second());
        const size_t &vertex_b( itr->Third());

        // get edge data for central bond from graph
        const size_t initial_edge_data( graph.GetEdgeData( vertex_a, vertex_b));

        // skip bonds in rings unless we are dealing with a purely ring-based
        if( initial_edge_data < size_t( 5) || m_IsPureRing)
        {
          // update the edge with information from the dihedral bins
          graph.EditEdge( vertex_a, vertex_b, initial_edge_data + radix * ( *itr_bins));
        }
      }

      return graph;
    }

    //! @brief get weights associated with each isomorphism
    storage::Vector< double> RotamerClusterCenter::GetIsomorphismWeights
    (
      const storage::Vector< storage::Vector< size_t> > &ISO,
      const size_t &MOLECULE_SIZE
    )
    {
      // get isomorphism from atom graph
      storage::Vector< double> isomorphism_counts( MOLECULE_SIZE, double( 0.0));
      for
      (
        storage::Vector< storage::Vector< size_t> >::const_iterator
          itr_iso( ISO.Begin()), itr_iso_end( ISO.End());
        itr_iso != itr_iso_end;
        ++itr_iso
      )
      {
        const storage::Vector< size_t> &inner_vec( *itr_iso);
        for
        (
          storage::Vector< size_t>::const_iterator itr_vec( inner_vec.Begin()), itr_vec_end( inner_vec.End());
          itr_vec != itr_vec_end;
          ++itr_vec
        )
        {
          isomorphism_counts( *itr_vec) += 1.0;
        }
      }
      storage::Vector< double> atom_count_fraction( isomorphism_counts);
      for
      (
        storage::Vector< double>::iterator
          itr_atom_count( atom_count_fraction.Begin()), itr_atom_count_end( atom_count_fraction.End());
        itr_atom_count != itr_atom_count_end;
        ++itr_atom_count
      )
      {
        *itr_atom_count = 1.0 / double( std::max( 1.0, *itr_atom_count));
      }
      storage::Vector< storage::Vector< size_t> >::const_iterator itr_simple_iso( ISO.Begin());

      storage::Vector< double> weights;
      weights.AllocateMemory( ISO.GetSize());

      for
      (
        storage::Vector< storage::Vector< size_t> >::const_iterator
          itr_simple_iso( ISO.Begin()), itr_simple_iso_end( ISO.End());
        itr_simple_iso != itr_simple_iso_end;
        ++itr_simple_iso
      )
      {
        const storage::Vector< size_t> &cur_iso( *itr_simple_iso);

        // weight of each isomorphism = 1/fragment_size * Summation(1/atom isomorphism count)
        double iso_weight( 0.0);
        for
        (
          storage::Vector< size_t>::const_iterator itr_cur_iso( cur_iso.Begin()), itr_cur_iso_end( cur_iso.End());
          itr_cur_iso != itr_cur_iso_end;
          ++itr_cur_iso
        )
        {
          iso_weight += atom_count_fraction( *itr_cur_iso);
        }
        iso_weight /= double( ISO.FirstElement().GetSize());
        weights.PushBack( iso_weight);
      }
      return weights;
    }

    //! @brief calculates the centers of cluster from ensembles stored in m_ConformationEnsembles
    //! @return void
    void RotamerClusterCenter::CalculateClusterCenters( size_t MAX_COUNTS)
    {
      // for each cluster find the center
      for
      (
        storage::List< RotamerEnsemble>::iterator itr_ensemble( m_Clusters.Begin()), itr_ensemble_end( m_Clusters.End());
        itr_ensemble != itr_ensemble_end;
        ++itr_ensemble
      )
      {
        itr_ensemble->CalculateCenter( MAX_COUNTS);
      }
    }

    //! @brief get an average structure to get a good structure
    //! @return void
    FragmentConformationShared RotamerClusterCenter::GetAverageStucture( size_t MAX_COUNTS)
    {
      storage::List< RotamerEnsemble>::iterator max_number;
      double max_counts( 0);
      for
      (
        storage::List< RotamerEnsemble>::iterator itr( m_Clusters.Begin()), itr_end( m_Clusters.End());
          itr != itr_end;
        ++itr
      )
      {
        if( itr->GetWeights() > max_counts)
        {
          max_number = itr;
          max_counts = itr->GetWeights();
        }
      }
      max_number->CalculateCenter( MAX_COUNTS);
      return max_number->GetClusterCenter();
    }

    void RotamerClusterCenter::PruneRotamersByWeight( double CUT_OFF)
    {
      for
      (
        storage::List< RotamerEnsemble>::iterator itr( m_Clusters.Begin()), itr_end( m_Clusters.End());
          itr != itr_end;
      )
      {
        if( itr->GetWeights() < CUT_OFF)
        {
          storage::List< RotamerEnsemble>::iterator itr_copy( itr);
          ++itr;
          m_Clusters.Remove( itr_copy);
        }
        else
        {
          ++itr;
        }
      }
    }

    void RotamerClusterCenter::PruneRotamersByInstances( size_t CUT_OFF)
    {
      for
      (
        storage::List< RotamerEnsemble>::iterator itr( m_Clusters.Begin()), itr_end( m_Clusters.End());
          itr != itr_end;
      )
      {
        if( itr->GetNumberInstances() < CUT_OFF)
        {
          storage::List< RotamerEnsemble>::iterator itr_copy( itr);
          ++itr;
          m_Clusters.Remove( itr_copy);
        }
        else
        {
          ++itr;
        }
      }
    }

    void RotamerClusterCenter::PruneNonPlanarAromaticRings()
    {
      size_t max_sz_planar( 0), max_sz_nonplanar( 0), total_planar( 0), total_nonplanar( 0);
      // count total size of planar vs non-planar aromatic rings
      for
      (
        storage::List< RotamerEnsemble>::iterator itr( m_Clusters.Begin()), itr_end( m_Clusters.End());
          itr != itr_end;
          ++itr
      )
      {
        if( !( itr->GetClusterCenter().AreAromaticRingsPlaner())) // && itr->GetClusterCenter().AreAmideBondsPlaner()))
        {
          total_nonplanar += itr->GetNumberInstances();
          max_sz_nonplanar = std::max( max_sz_nonplanar, itr->GetNumberInstances());
        }
        else
        {
          total_planar += itr->GetNumberInstances();
          max_sz_planar = std::max( max_sz_planar, itr->GetNumberInstances());
        }
      }
      if( total_planar >= total_nonplanar || max_sz_planar >= max_sz_nonplanar)
      {
        for
        (
          storage::List< RotamerEnsemble>::iterator itr( m_Clusters.Begin()), itr_end( m_Clusters.End());
            itr != itr_end;
        )
        {
          // This commented out block used to be in here, but when the algorithm was benchmarked, a bug in AreAmideBondsPlaner
          // caused it to only return true for terminal nitrogen groups (e.g. non-hydrogenated, terminal nitrogen) two bonds
          // away from an oxygen double bond; and was not checking planarity at all. Therefore, this block needs to be tested
          // with the bug-fixed version of the code before using this feature in production setting.
          // ~1% of all amide bonds in chains in the CSD are substantially (>10 degrees) non-planar due to sterics,
          // so it may be best to retain the non-planar amide bonds. Coupled amide bonds ring substituents often have
          // deflections of up to 30-40 degrees. Even higher deflections are sometimes seen with sulfonamid, where
          // no particular range of degrees is forbidden
          if( !( itr->GetClusterCenter().AreAromaticRingsPlaner())) // && itr->GetClusterCenter().AreAmideBondsPlaner()))
          {
            storage::List< RotamerEnsemble>::iterator itr_copy( itr);
            ++itr;
            m_Clusters.Remove( itr_copy);
          }
          else
          {
            ++itr;
          }
        }
      }
      else
      {
        BCL_MessageStd
        (
          "Skipping removal of non-planar aromatic rings - for this ring type there were more non-planar instances than "
          "planar "
        );
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RotamerClusterCenter::Read( std::istream &ISTREAM)
    {
      BCL_Exit( "Cannot read " + GetClassIdentifier(), -1);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &RotamerClusterCenter::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      BCL_Exit( "Cannot write " + GetClassIdentifier(), -1);
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
