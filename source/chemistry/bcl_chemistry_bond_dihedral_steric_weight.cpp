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
#include "chemistry/bcl_chemistry_bond_dihedral_steric_weight.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_undirected_edge.h"
#include "linal/bcl_linal_vector_operations.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new BondDihedralStericWeight
    BondDihedralStericWeight *BondDihedralStericWeight::Clone() const
    {
      return new BondDihedralStericWeight( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &BondDihedralStericWeight::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Determines steric weight at each side of the bond separately
    //! @param CONFORMATION the conformations of interest
    //! @param N_SPHERES number of spheres to go out
    storage::Vector< storage::Triplet< size_t, size_t, linal::Vector< double> > >
    BondDihedralStericWeight::CalculateDualSidedWeight
    (
      const FragmentComplete &CONFORMATION,
      const size_t N_SPHERES
    )
    {
      auto basic_graph( ConformationGraphConverter()( CONFORMATION));

      // accumulate VDW radius onto each node
      storage::Vector< double> vdw_total_each_node( basic_graph.GetSize(), double( 0.0));
      for( auto itr( CONFORMATION.GetAtomsIterator()); itr.NotAtEnd(); ++itr)
      {
        vdw_total_each_node( itr.GetPosition())
          += itr->GetAtomType()->GetAtomTypeProperty( AtomTypeData::e_CovalentRadiusTypical);
      }

      // this is a O(N^3) but could be made O(N^2) (possibly even O(N log(N))) if it ever becomes a bottleneck by transiting from tree leaves inward
      storage::Vector< storage::Triplet< size_t, size_t, linal::Vector< double> > > steric_data;
      steric_data.AllocateMemory( CONFORMATION.GetNumberBonds() * 2);
      const size_t mol_sz( CONFORMATION.GetSize());
      for( size_t atom_id( 0); atom_id < mol_sz; ++atom_id)
      {
        const storage::Vector< size_t> &basic_graph_neigh( basic_graph.GetNeighborIndices( atom_id));
        for( auto itr_neigh( basic_graph_neigh.Begin()), itr_neigh_end( basic_graph_neigh.End()); itr_neigh != itr_neigh_end; ++itr_neigh)
        {
          linal::Vector< double> vdw_sum( N_SPHERES + 1, double( 0.0));
          size_t atom_id_b( *itr_neigh);
          util::ShPtr< storage::Vector< size_t> > sh_ptr_dist
          (
            graph::Connectivity::DirectedDistancesToOtherVertices( basic_graph, atom_id, atom_id_b)
          );
          const storage::Vector< size_t> &distances( *sh_ptr_dist);
          for( size_t atom_id_c( 0); atom_id_c < mol_sz; ++atom_id_c)
          {
            if( distances( atom_id_c) <= N_SPHERES)
            {
              vdw_sum( distances( atom_id_c)) += vdw_total_each_node( atom_id_c);
            }
          }
          steric_data.PushBack
          (
            storage::Triplet< size_t, size_t, linal::Vector< double> >
            (
              atom_id,
              atom_id_b,
              vdw_sum
            )
          );
        }
      }
      return steric_data;
    }

    //! @brief Determines steric weight at each side of the bond separately
    //! @param CONFORMATION the conformations of interest
    //! @param N_SPHERES number of spheres to go out
    storage::Vector< graph::UndirectedEdge< double> >
    BondDihedralStericWeight::CalculateBondWeights
    (
      const FragmentComplete &CONFORMATION,
      const size_t N_SPHERES,
      const double DEPLETION_FACTOR
    )
    {
      storage::Vector< storage::Triplet< size_t, size_t, linal::Vector< double> > > bond_steric_factors
      (
        CalculateDualSidedWeight( CONFORMATION, N_SPHERES)
      );
      linal::Vector< double> weighting( N_SPHERES + 1, 1.0);
      for( size_t i( 1); i <= N_SPHERES; ++i)
      {
        weighting( i) = weighting( i - 1) * DEPLETION_FACTOR;
      }
      storage::List< storage::Triplet< size_t, size_t, double> > weight_left_right;
      std::map< storage::Pair< size_t, size_t>, storage::List< storage::Triplet< size_t, size_t, double> >::iterator> pair_to_iterator;
      for
      (
        auto itr( bond_steric_factors.Begin()), itr_end( bond_steric_factors.End());
        itr != itr_end;
        ++itr
      )
      {
        size_t low( 0), high( 0);
        if( itr->First() < itr->Second())
        {
          low = itr->First();
          high = itr->Second();
        }
        else
        {
          low = itr->Second();
          high = itr->First();
        }
        storage::Pair< size_t, size_t> ordered_pair( low, high);
        auto itr_to_map( pair_to_iterator.find( ordered_pair));
        if( itr_to_map == pair_to_iterator.end())
        {
          weight_left_right.PushBack
          (
            storage::Triplet< size_t, size_t, double>
            (
              low, high, linal::ScalarProduct( itr->Third(), weighting)
            )
          );
          pair_to_iterator[ ordered_pair] = weight_left_right.Last();
        }
        else
        {
          itr_to_map->second->Third() = std::min( itr_to_map->second->Third(), linal::ScalarProduct( itr->Third(), weighting));
        }
      }
      storage::Vector< graph::UndirectedEdge< double> > edges;
      edges.AllocateMemory( weight_left_right.GetSize());
      for( auto itr( weight_left_right.Begin()), itr_end( weight_left_right.End()); itr != itr_end; ++itr)
      {
        edges.PushBack( graph::UndirectedEdge< double>( itr->First(), itr->Second(), itr->Third()));
      }
      return edges;
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &BondDihedralStericWeight::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &BondDihedralStericWeight::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
