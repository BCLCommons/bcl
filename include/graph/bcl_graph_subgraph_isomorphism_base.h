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

#ifndef BCL_GRAPH_SUBGRAPH_ISOMORPHISM_BASE_H_
#define BCL_GRAPH_SUBGRAPH_ISOMORPHISM_BASE_H_

// include the namespace header
#include "bcl_graph.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_graph_const_graph.h"

// external includes - sorted alphabetically
namespace bcl
{
  namespace graph
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SubgraphIsomorphismBase
    //! @brief a common base class for all SubgraphIsomorphism that provides data and access methods
    //!
    //! @see @link example_graph_subgraph_isomorphism_base.cpp @endlink
    //! @author mendenjl
    //! @date Jun 04, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SubgraphIsomorphismBase :
      public util::ObjectInterface
    {

    protected:

    //////////
    // data //
    //////////

      //! all isomorphisms between subgraph and supergraph, index is the corresponding index of graph a
      storage::Vector< storage::Vector< size_t> > m_Isomorphisms;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! clone the object
      SubgraphIsomorphismBase *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @return class name of the object behind a pointer or the current object
      const std::string &GetClassIdentifier() const;

      //! @brief Get the isomorphisms between the sub and super graphs
      //! @return the isomorphisms
      const storage::Vector< storage::Vector< size_t> > &GetIsomorphisms() const;

      //! @brief Get the first subgraph isomorphism
      //! @return the first subgraph isomorphism
      const storage::Vector< size_t> &GetIsomorphism() const;

      //! @brief Get the first subgraph isomorphism
      //! @return the first subgraph isomorphism
      const storage::Vector< size_t> GetInverseIsomorphism() const;

      //! @brief set the isomorphism
      //! @param ISOMORPHISMS the new isomorphisms to use
      virtual void SetIsomorphisms( const storage::Vector< storage::Vector< size_t> > &ISOMORPHISMS);

      //! @brief set the isomorphism
      //! @param ISOMORPHISM the new isomorphism to use
      void SetIsomorphism( const storage::Vector< size_t> &ISOMORPHISM)
      {
        this->SetIsomorphisms( storage::Vector< storage::Vector< size_t> >( 1, ISOMORPHISM));
      }

      //! @brief prune isomorphisms to those that are unique
      //! @param VERTICES_DESIRED_TO_DIFFER vertex indices of subgraph that should differ to consider an isomorphism unique
      //!        from one found already
      void PruneIsomorphismsToThoseUniqueInField( const storage::Vector< size_t> &VERTICES_THAT_SHOULD_DIFFER);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    };

  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_SUBGRAPH_ISOMORPHISM_BASE_H_

