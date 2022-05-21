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
#include "chemistry/bcl_chemistry_fragment_graph_marker.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    //! @brief constructor
    //! @param BOND_TYPE the bond type to be used for graph isomorphism
    //! @param SPLIT_PARAMETER the type of splitting that needs to be done like getting only rings, chain, scaffolds etc
    FragmentGraphMarker::FragmentGraphMarker
    (
      const ConformationGraphConverter &BOND_TYPE,
      const util::Implementation< FragmentSplitInterface> &SPLIT_PARAMETER
    ) :
      m_Converter( BOND_TYPE),
      m_Splitter( SPLIT_PARAMETER)
    {
    }

    //! @brief Clone function
    //! @return pointer to new FragmentGraphMarker
    FragmentGraphMarker *FragmentGraphMarker::Clone() const
    {
      return new FragmentGraphMarker( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &FragmentGraphMarker::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the unique Constitution passed through this class
    //! @return unique Constitution passed through this class
    const util::ShPtrList< FragmentConstitutionShared> &FragmentGraphMarker::GetConstitution() const
    {
      return m_SpecialFragments.GetConstitutions();
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief Create a graph of a conformation
    //! @param CONFORMATION a conformation
    //! @return The conformation converted into a graph with the given atom/bond representations
    graph::ConstGraph< size_t, size_t> FragmentGraphMarker::operator()( const ConformationInterface &CONFORMATION) const
    {
      // create a graphs
      graph::ConstGraph< size_t, size_t> mol_graph( m_Converter( CONFORMATION));
      ConformationGraphConverter::t_AtomGraph mol_graph_atom( m_Converter.CreateGraphWithAtoms( CONFORMATION));

      // get rings, chains or scaffold depending on the splitting parameter
      storage::List< storage::Vector< size_t> > indices( m_Splitter->GetComponentVertices( CONFORMATION, mol_graph_atom));
      FragmentEnsemble components( m_Splitter->ConvertComponentsIntoEnsemble( CONFORMATION, indices, mol_graph_atom));

      // get iterator to components
      storage::List< storage::Vector< size_t> >::const_iterator itr_indices( indices.Begin());

      // edges need to be weighted with a value more than current highest value of edge
      const size_t radix( GetConstitutionalBondTypes().GetEnumCount());

      // for each component
      for
      (
        storage::List< FragmentComplete>::const_iterator itr( components.Begin()), itr_end( components.End());
        itr != itr_end;
        ++itr, ++itr_indices
      )
      {
        // if fragment has not been seen then continue
        ConstitutionSet::const_iterator itr_find( m_SpecialFragments.Find( FragmentConstitutionShared( *itr)));
        if( itr_find == m_SpecialFragments.End())
        {
          continue;
        }

        // get the fragment index
        size_t fragment_index( std::distance( m_SpecialFragments.Begin(), itr_find));

        // get subgraph of fragment in the molecule of interest
        graph::Subgraph< size_t, size_t> subgraph
        (
          util::OwnPtr< const graph::ConstGraph< size_t, size_t> >( &mol_graph, false),
          *itr_indices
        );

        // get edges of that fragment maps to in the molecule of interest and color the corresponding edges in molecule
        // graph with fragment_index
        storage::Vector< graph::UndirectedEdge< size_t> > subgraph_edges( subgraph.GetUndirectedEdges());
        for
        (
          storage::Vector< graph::UndirectedEdge< size_t> >::const_iterator
            itr_edge( subgraph_edges.Begin()), itr_edge_end( subgraph_edges.End());
          itr_edge != itr_edge_end;
          ++itr_edge
        )
        {
          const size_t initial_edge_data( mol_graph.GetEdgeData( itr_edge->GetIndexLow(), itr_edge->GetIndexHigh()));
          BCL_Assert( initial_edge_data < radix, "Out of bounds initial edge data: " + util::Format()( initial_edge_data));
          mol_graph.EditEdge( itr_edge->GetIndexLow(), itr_edge->GetIndexHigh(), initial_edge_data + radix * ( fragment_index + 1));
        }
      }
      return mol_graph;
    }

    //! @brief create a vector of graphs of conformations that are passed as an ensemble
    //! @param ENSEMBLE ensemble of conformations whose graphs are desired
    //! @return a vector of graphs of a conformations that are passed as an ensemble
    storage::Vector< graph::ConstGraph< size_t, size_t> > FragmentGraphMarker::operator()( const FragmentEnsemble &ENSEMBLE) const
    {
      storage::Vector< graph::ConstGraph< size_t, size_t> > graphs;
      size_t scaffold_number( 0);
      for
      (
        storage::List< FragmentComplete>::const_iterator itr( ENSEMBLE.Begin()), itr_end( ENSEMBLE.End());
        itr != itr_end;
        ++itr, ++scaffold_number
      )
      {
        graphs.PushBack( operator()( *itr));
      }
      return graphs;
    }

    //! @brief Create a graph of a conformation and store Constitution in a Constitutionset
    //! @param CONFORMATION a conformation
    //! @return The conformation converted into a graph with the given atom/bond representations
    graph::ConstGraph< size_t, size_t> FragmentGraphMarker::Insert( const ConformationInterface &CONFORMATION)
    {
      FragmentEnsemble components( ( *m_Splitter)( CONFORMATION));

      // insert conformation into Constitution set
      for
      (
        storage::List< FragmentComplete>::const_iterator itr( components.Begin()), itr_end( components.End());
        itr != itr_end;
        ++itr
      )
      {
        m_SpecialFragments.Insert( FragmentConstitutionShared( *itr));
      }
      return operator ()( CONFORMATION);
    }

    //! @brief create a vector of graphs of conformations that are passed as an ensemble and storing the ensemble in a configuratin set
    //! @param ENSEMBLE ensemble of conformations whose graphs are desired
    //! @return create a vector of graphs of conformations that are passed as an ensemble
    storage::Vector< graph::ConstGraph< size_t, size_t> > FragmentGraphMarker::Insert( const FragmentEnsemble &ENSEMBLE)
    {
      storage::Vector< graph::ConstGraph< size_t, size_t> > graphs;

      // insert each conformation into Constitution set
      for
      (
        storage::List< FragmentComplete>::const_iterator itr( ENSEMBLE.Begin()), itr_end( ENSEMBLE.End());
        itr != itr_end;
        ++itr
      )
      {
        graphs.PushBack( Insert( *itr));
      }
      return graphs;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FragmentGraphMarker::Read( std::istream &ISTREAM)
    {
      BCL_Exit( "Cannot read " + GetClassIdentifier(), -1);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &FragmentGraphMarker::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      BCL_Exit( "Cannot write " + GetClassIdentifier(), -1);
      return OSTREAM;
    }

  } // namespace chemistry
} // namespace bcl
