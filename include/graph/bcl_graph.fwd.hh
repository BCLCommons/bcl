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

#ifndef BCL_GRAPH_FWD_HH_
#define BCL_GRAPH_FWD_HH_

// include bcl_defines.h header
#include "bcl_defines.h"

// This file contains forward declarations for the graph namespace
// This file is mostly automatically generated
// Developers may add typedefs and template default parameters
// all other changes will be removed the next time this file is generated
namespace bcl
{
  namespace graph
  {
  /////////////////////
  // regular classes //
  /////////////////////

    class CSISubstructure;
    class CommonSubgraphIsomorphismBase;
    class Connectivity;
    class EdgeCoverRingPerception;
    class ExhaustiveRingPerception;
    class Path;
    class Ring;
    class SubgraphIsomorphismBase;

  //////////////////////
  // template classes //
  //////////////////////

    template< typename t_VertexData, typename t_EdgeData>
    class CommonSubgraphIsomorphism;

    template< typename t_VertexData, typename t_EdgeData>
    class ConstGraph;

    template< typename t_VertexDataType, typename t_EdgeDataType>
    class EdgeWithData;

    template< typename t_VertexDataType, typename t_EdgeDataType>
    class GraphWithData;

    template< typename t_VertexData, typename t_EdgeData>
    class Subgraph;

    template< typename t_VertexData, typename t_EdgeData>
    class SubgraphIsomorphism;

    template< typename t_VertexData, typename t_EdgeData>
    class TreeNode;

    template< typename t_EdgeData>
    class UndirectedEdge;

    template< typename t_VertexDataType, typename t_EdgeDataType>
    class VertexWithData;

  //////////////
  // typedefs //
  //////////////

  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_FWD_HH_
