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

#ifndef BCL_GRAPH_UNDIRECTED_EDGE_H_
#define BCL_GRAPH_UNDIRECTED_EDGE_H_

// include the namespace header
#include "bcl_graph.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialize.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class UndirectedEdge
    //! @brief Manages generic undirected edge info
    //!
    //! @see @link example_graph_undirected_edge.cpp @endlink
    //! @author mendenjl
    //! @date Jun 06, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_EdgeData>
    class UndirectedEdge :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      size_t     m_IndexLow;   //!< Smaller atom index
      size_t     m_IndexHigh;  //!< Larger atom index
      t_EdgeData m_EdgeData;   //!< Data type of edge

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      UndirectedEdge() :
        m_IndexLow( util::GetUndefinedSize_t()),
        m_IndexHigh( util::GetUndefinedSize_t()),
        m_EdgeData()
      {
      }

      //! @brief constructor from members
      UndirectedEdge
      (
        const size_t &INDEX_A,
        const size_t &INDEX_B,
        const t_EdgeData &EDGE_TYPE
      ) :
        m_IndexLow( std::min( INDEX_A, INDEX_B)),
        m_IndexHigh( std::max( INDEX_A, INDEX_B)),
        m_EdgeData( EDGE_TYPE)
      {
      }

      //! @brief Clone function
      //! @return pointer to new UndirectedEdge
      UndirectedEdge *Clone() const
      {
        return new UndirectedEdge( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get the smaller of the two indices
      //! @return the smaller of the two indices
      const size_t &GetIndexLow() const
      {
        return m_IndexLow;
      }

      //! @brief get the larger of the two indices
      //! @return the larger of the two indices
      const size_t &GetIndexHigh() const
      {
        return m_IndexHigh;
      }

      //! @brief get the edge data
      //! @return the edge data
      const t_EdgeData &GetEdgeData() const
      {
        return m_EdgeData;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief compare this bond info to other bond info
      //! @param BOND_INFO the object to compare this bond info to
      //! @return true if the bond info is the same
      bool operator ==( const UndirectedEdge &EDGE) const
      {
        return m_IndexLow == EDGE.m_IndexLow && m_IndexHigh == EDGE.m_IndexHigh && m_EdgeData == EDGE.m_EdgeData;
      }

      //! @brief compare this bond info to other bond info
      //! @param BOND_INFO the object to compare this bond info to
      //! @return true if the bond info is the same
      bool operator <( const UndirectedEdge &EDGE) const
      {
        if( m_IndexLow < EDGE.m_IndexLow)
        {
          return true;
        }
        return m_IndexLow == EDGE.m_IndexLow && m_IndexHigh < EDGE.m_IndexHigh;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        io::Serialize::Read( m_IndexLow,  ISTREAM);
        io::Serialize::Read( m_IndexHigh, ISTREAM);
        io::Serialize::Read( m_EdgeData, ISTREAM);

        // return
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        io::Serialize::Write( m_IndexLow,  OSTREAM, INDENT) << '\t';
        io::Serialize::Write( m_IndexHigh, OSTREAM, 0)      << '\n';
        io::Serialize::Write( m_EdgeData,  OSTREAM, INDENT);

        return OSTREAM;
      }

    }; // class UndirectedEdge

  } // namespace graph
} // namespace bcl

#endif //BCL_GRAPH_UNDIRECTED_EDGE_H_

