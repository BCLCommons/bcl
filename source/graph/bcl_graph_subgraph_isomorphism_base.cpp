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
#include "graph/bcl_graph_subgraph_isomorphism_base.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace graph
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SubgraphIsomorphismBase::s_Instance
    (
      GetObjectInstances().AddInstance( new SubgraphIsomorphismBase())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! clone the object
    SubgraphIsomorphismBase *SubgraphIsomorphismBase::Clone() const
    {
      return new SubgraphIsomorphismBase( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @return class name of the object behind a pointer or the current object
    const std::string &SubgraphIsomorphismBase::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Get the isomorphisms between the sub and super graphs
    //! @return the isomorphisms
    const storage::Vector< storage::Vector< size_t> > &SubgraphIsomorphismBase::GetIsomorphisms() const
    {
      return m_Isomorphisms;
    }

    //! @brief Get the first subgraph isomorphism
    //! @return the first subgraph isomorphism
    const storage::Vector< size_t> &SubgraphIsomorphismBase::GetIsomorphism() const
    {
      static storage::Vector< size_t> empty;
      return m_Isomorphisms.IsEmpty() ? empty : m_Isomorphisms.FirstElement();
    }

    //! @brief Get the first subgraph isomorphism
    //! @return the first subgraph isomorphism
    const storage::Vector< size_t> SubgraphIsomorphismBase::GetInverseIsomorphism() const
    {
      static storage::Vector< size_t> empty;
      if( m_Isomorphisms.IsEmpty())
      {
        return empty;
      }

      const storage::Vector< size_t> &isomorphism( m_Isomorphisms.FirstElement());
      storage::Vector< size_t> inverse( isomorphism.GetSize());

      size_t count( 0);
      for
      (
        storage::Vector< size_t>::const_iterator itr( isomorphism.Begin()), itr_end( isomorphism.End());
        itr != itr_end;
        ++itr, ++count
      )
      {
        inverse( *itr) = count;
      }
      return inverse;
    }

    //! @brief set the isomorphisms
    //! @param ISOMORPHISMS the new isomorphisms to use
    void SubgraphIsomorphismBase::SetIsomorphisms
    (
      const storage::Vector< storage::Vector< size_t> > &ISOMORPHISMS
    )
    {
      m_Isomorphisms = ISOMORPHISMS;
    }

    //! @brief prune isomorphisms to those that are unique
    //! @param VERTICES_DESIRED_TO_DIFFER vertex indices of subgraph that should differ to consider an isomorphism unique
    //!        from one found already
    void SubgraphIsomorphismBase::PruneIsomorphismsToThoseUniqueInField( const storage::Vector< size_t> &VERTICES_THAT_SHOULD_DIFFER)
    {
      storage::Set< storage::Vector< size_t> > iso_uniq;
      auto itr_place( m_Isomorphisms.Begin()), itr_next( m_Isomorphisms.Begin());
      for( auto itr_end( m_Isomorphisms.End()); itr_next != itr_end; ++itr_next)
      {
        storage::Vector< size_t> iso_key( *itr_next);
        iso_key.Reorder( VERTICES_THAT_SHOULD_DIFFER);
        if( iso_uniq.Insert( iso_key).second)
        {
          if( itr_place != itr_next)
          {
            *itr_place = *itr_next;
          }
          ++itr_place;
        }
      }
      const size_t n_size( std::distance( m_Isomorphisms.Begin(), itr_place));
      m_Isomorphisms.Resize( n_size);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SubgraphIsomorphismBase::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Isomorphisms, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &SubgraphIsomorphismBase::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Isomorphisms, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  } // namespace graph
} // namespace bcl
