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
#include "graph/bcl_graph_edge_cover_ring_perception.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_const_graph.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> EdgeCoverRingPerception::s_Instance
    (
      GetObjectInstances().AddInstance( new EdgeCoverRingPerception())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    EdgeCoverRingPerception::EdgeCoverRingPerception() :
      m_Rings()
    {
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &EdgeCoverRingPerception::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Rings, ISTREAM);

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &EdgeCoverRingPerception::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write
      io::Serialize::Write( m_Rings, OSTREAM, INDENT);

      return OSTREAM;
    }

  } // namespace graph
} // namespace bcl
