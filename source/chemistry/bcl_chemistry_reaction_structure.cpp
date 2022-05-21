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
#include "chemistry/bcl_chemistry_reaction_structure.h"

// includes from bcl - sorted alphabetically
#include "graph/bcl_graph_subgraph_isomorphism.h"

// external includes - sorted alphabetically
//#include <stdio.h>

namespace bcl
{
  namespace chemistry
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! default constructor
    ReactionStructure::ReactionStructure() :
      m_ElementGraph()
    {
    }

    //! @brief constructor with molecular structure and aromaticity specification
    //! @param FRAGMENT to use
    //! @param AROMATICITY aromaticity specification of all atoms
    ReactionStructure::ReactionStructure
    (
      const FragmentComplete &FRAGMENT,
      const std::vector< Aromaticity> &ALLOWED_AROMATICITY
    ) :
      m_Fragment( FRAGMENT),
      m_ElementGraph(),
      m_AllowedAromaticity( ALLOWED_AROMATICITY)
    {
      size_t n_atoms( m_Fragment.GetNumberAtoms());

      // make the graph
      m_ElementGraph = s_Converter( FRAGMENT);

      if( m_AllowedAromaticity.empty())
      {
        // arbitrary aliphatic or aromatic character
        m_AllowedAromaticity = std::vector< Aromaticity>( n_atoms, Aromaticity( e_AliphaticOrAromatic));
      }
      else
      {
        BCL_Assert
        ( 
          m_AllowedAromaticity.size() == n_atoms, 
          "ReactionStructure: allowed aromaticity vector does not have a size that matches the provided molecule graph"
        );
      }
    }

    //! @brief Clone function
    //! @return pointer to new ReactionStructure
    ReactionStructure *ReactionStructure::Clone() const
    {
      return new ReactionStructure( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief helper function to determine if a query aromaticity is allowed according to another aromaticity
    //! @param QUERY the aromaticity to check, will likely be one of e_Aromatic or e_Aliphatic
    //! @param TEST the aromaticity to check against, may be any Aromaticity specification
    //! @return true if QUERY is an allowed aromaticity value accordin to TEST
    //! @details since Aromaticity is implemented as a bitfield, this means that all bits set in QUERY are also set in TEST
    bool ReactionStructure::AromaticityMatches( const Aromaticity &QUERY, const Aromaticity &TEST)
    {
      return QUERY & TEST;
    }

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ReactionStructure::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get the internally stored ElementType/BondOrderOrAromatic graph
    const graph::ConstGraph< size_t, size_t> &ReactionStructure::GetGraph() const
    {
      return m_ElementGraph;
    }

    //! @brief get the number of atoms in the structure
    const size_t &ReactionStructure::GetSize() const
    {
      return m_ElementGraph.GetSize();
    }

    //! @brief get internal molecular structure data
    const FragmentComplete &ReactionStructure::GetFragment() const
    {
      return m_Fragment;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief determine if this structure contains another
    //! @param OTHER the reaction structure to check for
    //! @return true if this structure contains OTHER (or is equal to)
    bool ReactionStructure::Contains( const ReactionStructure &OTHER) const
    {
      const graph::ConstGraph< size_t, size_t> &other_graph( OTHER.m_ElementGraph);
      if( other_graph.GetSize() < m_ElementGraph.GetSize())
      {
        return false;
      }

      return !GetMatchingSubstructures( other_graph, OTHER.m_AllowedAromaticity).IsEmpty();
    }

    //! @brief get subgraph isomorphisms of one graph with graph from this object, including matching aromaticity
    //! @param OTHER_GRAPH the other graph to search
    //! @param OTHER_AROMATICITY aromaticity specification of atoms in OTHER_GRAPH
    //! @return a list of isomorphisms for substructures of OTHER_GRAPH that are equivalent to this one
    storage::Vector< storage::Vector< size_t> > ReactionStructure::GetMatchingSubstructures
    ( 
      const graph::ConstGraph< size_t, size_t> &OTHER_GRAPH,
      const std::vector< Aromaticity> &OTHER_AROMATICITY
    ) const
    {
      storage::Vector< storage::Vector< size_t> > matched_isos;

      graph::SubgraphIsomorphism< size_t, size_t> iso_search;
      iso_search.SetGraphExternalOwnership( OTHER_GRAPH);
      iso_search.SetSubgraphExternalOwnership( m_ElementGraph);

      // find all isomorphisms that are not overlapping on an element level
      if( iso_search.FindDisparateIsomorphisms())
      {
        const storage::Vector< storage::Vector< size_t> > &isos( iso_search.GetIsomorphisms());

        // check aromaticity
        for( size_t iso_no( 0), end_iso( isos.GetSize()); iso_no < end_iso; ++iso_no)
        {
          const storage::Vector< size_t> &iso( isos( iso_no));
          bool keep_iso( true);
          for( size_t ano( 0), end_ano( iso.GetSize()); ano < end_ano; ++ano)
          {
            if( !AromaticityMatches( OTHER_AROMATICITY[ iso( ano)], m_AllowedAromaticity[ ano]))
            {
              keep_iso = false;
              break;
            }
          }
          if( keep_iso)
          {
            matched_isos.PushBack( iso);
          }
        }
      }
      return matched_isos;
    }

    //! @brief get a list of isomorphisms that match this structure
    //! @param CONFORMATION the conformation to search
    //! @return a list of atom-to-atom isomorphisms that match this structure
    storage::Vector< storage::Vector< size_t> > ReactionStructure::GetMatchingSubstructures
    ( 
      const ConformationInterface &CONFORMATION
    ) const
    {
      storage::Vector< storage::Vector< size_t> > matched_isos;

      // create the graph
      graph::ConstGraph< size_t, size_t> mol_graph( s_Converter( CONFORMATION));

      if( mol_graph.GetSize() < m_ElementGraph.GetSize())
      {
        return matched_isos; // empty
      }

      // Set up aromaticity
      std::vector< Aromaticity> mol_arom( mol_graph.GetSize(), Aromaticity( e_Aliphatic));
      storage::Vector< graph::UndirectedEdge< size_t> > mol_bonds
      ( 
        CONFORMATION.GetAdjacencyList( ConfigurationalBondTypeData::e_BondOrderOrAromatic)
      );
      const size_t aromatic_val( 4);
      for
      (
        storage::Vector< graph::UndirectedEdge< size_t> >::const_iterator itr_edge( mol_bonds.Begin()), 
          itr_edge_end( mol_bonds.End());
        itr_edge != itr_edge_end;
        ++itr_edge
      )
      {
        if( itr_edge->GetEdgeData() == aromatic_val)
        {
          mol_arom[ itr_edge->GetIndexLow()] = Aromaticity( e_Aromatic);
          mol_arom[ itr_edge->GetIndexHigh()] = Aromaticity( e_Aromatic);
        }
      }

      matched_isos = GetMatchingSubstructures( mol_graph, mol_arom);

      return matched_isos;
    }

    //! @brief determine if this structure matches parts of a molecule
    //! @param CONFORMATION the conformation to check
    //! @return true if this structure matches any parts of CONFORMATION
    bool ReactionStructure::ContainedIn( const ConformationInterface &CONFORMATION) const
    {
      return !GetMatchingSubstructures( CONFORMATION).IsEmpty();
    }

    //! @brief equality operator
    //! @details graphs are same size, isomorphic, and all properties (aromaticity) are equivalent
    bool ReactionStructure::operator ==( const ReactionStructure &OTHER) const
    {

      bool is_equal( true);
      if( m_ElementGraph.GetSize() == OTHER.m_ElementGraph.GetSize())
      {
        storage::Vector< storage::Vector< size_t> > iso1( GetMatchingSubstructures( OTHER.m_ElementGraph, OTHER.m_AllowedAromaticity));
        storage::Vector< storage::Vector< size_t> > iso2( OTHER.GetMatchingSubstructures( m_ElementGraph, m_AllowedAromaticity));
        if( iso1.IsEmpty() || iso2.IsEmpty())
        {
          is_equal = false;
        }

      }
      else
      {
        is_equal = false;
      }

      return is_equal;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ReactionStructure::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ReactionStructure::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

    // instantiate the converter with proper arguments
    const ConformationGraphConverter ReactionStructure::s_Converter
    ( 
      ConformationGraphConverter::e_ElementType,
      ConfigurationalBondTypeData::e_BondOrderOrAromatic
    );

  } // namespace chemistry
} // namespace bcl
