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
#include "graph/bcl_graph_common_subgraph_isomorphism_base.h"

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
    const util::SiPtr< const util::ObjectInterface> CommonSubgraphIsomorphismBase::s_Instance
    (
      GetObjectInstances().AddInstance( new CommonSubgraphIsomorphismBase())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a solution type (connected by default)
    //! @param SOLUTION_TYPE whether the isomorphism must be connected
    CommonSubgraphIsomorphismBase::CommonSubgraphIsomorphismBase( const SolutionType &SOLUTION_TYPE) :
      m_SolutionType( SOLUTION_TYPE)
    {
    }

    //! clone the object
    CommonSubgraphIsomorphismBase *CommonSubgraphIsomorphismBase::Clone() const
    {
      return new CommonSubgraphIsomorphismBase( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @return class name of the object behind a pointer or the current object
    const std::string &CommonSubgraphIsomorphismBase::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief Get the isomorphisms between graphs A and B
    //! @return the isomorphisms
    const storage::Vector< storage::Map< size_t, size_t> > &CommonSubgraphIsomorphismBase::GetIsomorphisms() const
    {
      return m_Isomorphisms;
    }

    //! @brief Get the largest common subgraph isomorphism between graphs A and B
    //! @return the largest common subgraph isomorphism
    const storage::Map< size_t, size_t> &CommonSubgraphIsomorphismBase::GetIsomorphism() const
    {
      static storage::Map< size_t, size_t> empty;
      return m_Isomorphisms.IsEmpty() ? empty : m_Isomorphisms.FirstElement();
    }

    //! @brief Get whether this isomorphism had to be connected
    //! @return bool, true indicates this isomorphism had to be connected
    CommonSubgraphIsomorphismBase::SolutionType CommonSubgraphIsomorphismBase::GetSolutionType() const
    {
      return m_SolutionType;
    }

    //! @brief set the isomorphisms
    //! @param ISOMORPHISMS the new isomorphisms to use
    void CommonSubgraphIsomorphismBase::SetIsomorphisms
    (
      const storage::Vector< storage::Map< size_t, size_t> > &ISOMORPHISMS
    )
    {
      m_Isomorphisms = ISOMORPHISMS;
    }

    //! @brief set the isomorphism, assuming there is only one
    //! @param ISOMORPHISM the new isomorphism to use
    void CommonSubgraphIsomorphismBase::SetIsomorphism( const storage::Map< size_t, size_t> &ISOMORPHISM)
    {
      m_Isomorphisms.Reset();
      m_Isomorphisms.PushBack( ISOMORPHISM);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CommonSubgraphIsomorphismBase::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Isomorphisms, ISTREAM);

      size_t enum_value;
      ISTREAM >> enum_value;
      m_SolutionType = SolutionType( enum_value);

      // return
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &CommonSubgraphIsomorphismBase::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Isomorphisms, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SolutionType, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief swaps the isomorphisms
    //! thus makes indices in graph B be the keys to the isomorphism which map to indices in graph a
    void CommonSubgraphIsomorphismBase::SwapIsomorphisms()
    {
      for
      (
        storage::Vector< storage::Map< size_t, size_t> >::iterator
          itr_isos( m_Isomorphisms.Begin()), itr_isos_end( m_Isomorphisms.End());
        itr_isos != itr_isos_end;
        ++itr_isos
      )
      {
        // swap the isomorphism; used if graph_a was originally the larger graph

        // create a  new isomorphism, which will hold the swapped values
        storage::Map< size_t, size_t> new_isomorphism;
        for
        (
          storage::Map< size_t, size_t>::const_iterator iso( itr_isos->Begin()), iso_end( itr_isos->End());
          iso != iso_end;
          ++iso
        )
        {
          new_isomorphism[ iso->second] = iso->first;
        }

        // swap the isomorphism maps
        itr_isos->Swap( new_isomorphism);
      }
    }

  } // namespace graph
} // namespace bcl
