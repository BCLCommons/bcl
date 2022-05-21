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
#include "assemble/bcl_assemble_voxel_grid_mutation.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "biol/bcl_biol_aa_base.h"
#include "graph/bcl_graph_const_graph.h"
#include "graph/bcl_graph_undirected_edge.h"
#include "math/bcl_math_running_average_sd.h"
#include "score/bcl_score_aa_pair_clash.h"
#include "score/bcl_score_aa_pair_hi_res_clash.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //! @brief extract the 3D coordinates of a given t_DataType input. TO BE IMPLEMENTED BY DERIVED CLASSES if one input can have multiple positions
    //! @param INPUT pointer to the t_DataType
    //! @return 3D Vector of the input's coordinates
    util::SiPtrVector< const linal::Vector3D> VoxelGridMutation::ExtractPositions
    (
      const util::SiPtr< const biol::Mutation> &INPUT
    ) const
    {
      util::SiPtrVector< const linal::Vector3D> positions;
      positions.AllocateMemory( INPUT->GetChainIDs().size());
      for( size_t i( 0), sz( INPUT->GetAAs().GetSize()); i < sz; ++i)
      {
        positions.PushBack( ExtractPositionFromChain( INPUT, i));
        if( !positions.LastElement().IsDefined())
        {
          positions.PopBack();
        }
      }
      return positions;
    }

    //! @brief Extract an AA-instance's coordinates (Ha2 for Glycine, Cb-Atom otherwise)
    //! @param AA pointer to the amino acid object
    //! @return Reference to the coordinates. null-ptr if coordinates are undefined
    util::SiPtr< const linal::Vector3D> VoxelGridMutation::ExtractPositionFromChain
    (
      const util::SiPtr< const biol::Mutation> &MUTATION,
      const size_t &CHAIN_INDEX
    ) const
    {
      if( m_PreferCenter)
      {
        auto coords( MUTATION->GetAAs()( CHAIN_INDEX)->GetAtomCoordinates());
        math::RunningAverage< linal::Vector3D> ave_coord;
        for( auto itr( coords.Begin()), itr_end( coords.End()); itr != itr_end; ++itr)
        {
          if( ( *itr)->IsDefined())
          {
            ave_coord += **itr;
          }
        }
        if( ave_coord.GetWeight())
        {
          m_Centers.PushBack( ave_coord.GetAverage());
          return util::SiPtr< const linal::Vector3D>( &m_Centers.LastElement());
        }
      }
      // get firstSideChainAtom
      const biol::Atom &atom( m_PreferCA ? MUTATION->GetAAs()( CHAIN_INDEX)->GetCA() : MUTATION->GetAAs()( CHAIN_INDEX)->GetFirstSidechainAtom());

      // Check if it's coordinates are defined
      if( atom.AllCoordinatesDefined())
      { // if so return them
        return util::SiPtr< const linal::Vector3D>( atom.GetCoordinates());
      }
      return util::SiPtr< const linal::Vector3D>();
    }
  } // namespace assemble
} // namespace bcl
