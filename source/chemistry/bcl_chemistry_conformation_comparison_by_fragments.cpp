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
#include "chemistry/bcl_chemistry_conformation_comparison_by_fragments.h"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_atom_constitutional_interface.h"
#include "chemistry/bcl_chemistry_constitution_set.h"
#include "graph/bcl_graph_subgraph.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace chemistry
  {

    // add the interface to the set of known implementations
    // single instance of the class that only counts the number of common fragments
    const util::SiPtr< const util::ObjectInterface> ConformationComparisonByFragments::s_RawInstance
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance
      (
        new ConformationComparisonByFragments( false)
      )
    );

    // single instance of the tanimoto class
    const util::SiPtr< const util::ObjectInterface> ConformationComparisonByFragments::s_TanimotoInstance
    (
      util::Enumerated< ConformationComparisonInterface>::AddInstance
      (
        new ConformationComparisonByFragments( true)
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor, takes whether to compute tanimoto vs raw values, and whether to enforce a connected solution
    //! @param TANIMOTO whether to compute tanimoto values (vs. raw size of substructure)
    //! @param BINARY if true, it only looks for a single match of a fragment
    ConformationComparisonByFragments::ConformationComparisonByFragments( const bool &TANIMOTO) :
      m_ComputeTanimoto( TANIMOTO)
    {
    }

    //! virtual copy constructor
    ConformationComparisonByFragments *ConformationComparisonByFragments::Clone() const
    {
      return new ConformationComparisonByFragments( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! the class name as const ref std::string
    const std::string &ConformationComparisonByFragments::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ConformationComparisonByFragments::GetAlias() const
    {
      static const std::string s_tanimoto( "CommonFragmentsTanimoto"),
                               s_raw( "CommonFragments");

      return m_ComputeTanimoto ? s_tanimoto : s_raw;
    }

    //! @brief gets the splits that this class will use
    //! @return a vector of split implementations
    const storage::Vector< util::Implementation< FragmentSplitInterface> > &ConformationComparisonByFragments::GetSplits() const
    {
      return m_Splits;
    }

    //! @brief adds a split to the class
    //! @param SPLIT the spliting method to add
    void ConformationComparisonByFragments::AddSplit( const FragmentSplitInterface &SPLIT)
    {
      m_Splits.PushBack( util::Implementation< FragmentSplitInterface>( SPLIT));
    }

  /////////////////
  //  operations //
  /////////////////

    //! @brief align two small molecule objects and find the RMSD
    //! @param MOLECULE_A - first molecule being aligned
    //! @param MOLECULE_A - second molecule being aligned
    //! @return the RMSD between first and second molecule
    double ConformationComparisonByFragments::operator()
    (
      const ConformationInterface &MOLECULE_A,
      const ConformationInterface &MOLECULE_B
    ) const
    {
      
      // iterators to map elements that contain graphs of molecule A and B
      storage::Map
      < 
        util::SiPtr< const ConformationInterface>, 
        storage::Vector< graph::ConstGraph< size_t, size_t> >
      >::const_iterator 
        mol_a_itr( m_FragmentGraphs.Find( util::SiPtr< const ConformationInterface>( MOLECULE_A))),
        mol_b_itr( m_FragmentGraphs.Find( util::SiPtr< const ConformationInterface>( MOLECULE_B)));

      util::OwnPtr< const storage::Vector< graph::ConstGraph< size_t, size_t> > > op_mol_a_graphs, op_mol_b_graphs;

      // Get a pointer to the vector of graphs for molecule A, or make one if it doesn't exist
      if( mol_a_itr == m_FragmentGraphs.End())
      {
        // pointer to a new vector of graphs
        op_mol_a_graphs = util::OwnPtr< const storage::Vector< graph::ConstGraph< size_t, size_t> > >
                          (
                            new storage::Vector< graph::ConstGraph< size_t, size_t> >( MakeFragmentGraphs( MOLECULE_A)),
                            true
                          ); 
      }
      else
      {
        // pointer to a pre-existing vector of graphs
        op_mol_a_graphs = util::OwnPtr< const storage::Vector< graph::ConstGraph< size_t, size_t> > >
                          (
                            &( mol_a_itr->second),
                            false
                          );
      }

      // Get a pointer to the vector of graphs for molecule B, or make one if it doesn't exist
      if( mol_b_itr == m_FragmentGraphs.End())
      {
        // new vector of graphs
        op_mol_b_graphs = util::OwnPtr< const storage::Vector< graph::ConstGraph< size_t, size_t> > >
                          (
                            new storage::Vector< graph::ConstGraph< size_t, size_t> >( MakeFragmentGraphs( MOLECULE_B)),
                            true
                          ); 
      }
      else
      {
        // pre-existing vector of graphs
        op_mol_b_graphs = util::OwnPtr< const storage::Vector< graph::ConstGraph< size_t, size_t> > >
                          (
                            &( mol_b_itr->second),
                            false
                          );
      }

      // boolean vectors of which fragments have already been matched by which in each molecule
      linal::Vector< size_t> matched_mol_a_graphs( op_mol_a_graphs->GetSize(), size_t( 0));
      linal::Vector< size_t> matched_mol_b_graphs( op_mol_b_graphs->GetSize(), size_t( 0));

      // check every fragment of molecule A for its occurrence in molecule B
      size_t n_b_graphs( op_mol_b_graphs->GetSize());
      for( size_t a_no( 0), end_a_no( op_mol_a_graphs->GetSize()); a_no < end_a_no; ++a_no)
      {
        // find the fragment that matches this one
        size_t matched_graph( GraphVectorFind( ( *op_mol_a_graphs)( a_no), *op_mol_b_graphs));
        if( matched_graph != n_b_graphs)
        {
          /* Skip this for now
          if( matched_mol_b_graphs( matched_graph) == 1)
          {
            BCL_MessageStd
            ( 
              "Fragment " + util::Format()( matched_graph) + 
              " in molecule B matched multiple fragments in molecule A, including " + 
              util::Format()( a_no) + ".  Fragment comparisons must be one-to-one!"
            );
            return util::GetUndefined< double>();
          }
          */
          matched_mol_b_graphs( matched_graph) = 1;
          matched_mol_a_graphs( a_no) = 1;
        }
      }

      size_t common( matched_mol_a_graphs.Sum());
      size_t a_unique( matched_mol_a_graphs.GetSize() - common);
      size_t b_unique( matched_mol_b_graphs.GetSize() - common);
      size_t unique( common + a_unique + b_unique);

      double tanimoto( unique != 0 ? double( common) / double( unique) : 0);

      return m_ComputeTanimoto ? tanimoto : common;
    }

    //! @brief prepare the class for comparing a conformation
    //! @param MOLECULE the molecule to prepare to compare
    void ConformationComparisonByFragments::Prepare( const ConformationInterface &MOLECULE) const
    {
      m_FragmentGraphs[ util::SiPtr< const ConformationInterface>( &MOLECULE)] = MakeFragmentGraphs( MOLECULE);
    }

    //! @brief makes graphs of fragments of a molecule by calling the split functions and checking uniqueness
    //! @param MOLECULE the molecule to fragment
    //! @return a vector of unique graphs that represent fragments
    storage::Vector< graph::ConstGraph< size_t, size_t> > ConformationComparisonByFragments::MakeFragmentGraphs
    (
      const ConformationInterface &MOLECULE
    ) const
    {
      storage::Vector< graph::ConstGraph< size_t, size_t> > graphs;
      FragmentEnsemble frags;

      // Generate fragments by calling the split operator()
      for( size_t split_no( 0), end_no( m_Splits.GetSize()); split_no < end_no; ++split_no)
      {
        frags.Append( ( *m_Splits( split_no))( MOLECULE));
      }

      graphs.AllocateMemory( frags.GetSize());

      // make a graph of each fragment and compare it to the graphs that have been generated
      for
      ( 
        FragmentEnsemble::const_iterator itr_mol( frags.Begin()), itr_mol_end( frags.End());
        itr_mol != itr_mol_end;
        ++itr_mol
      )
      {
        graph::ConstGraph< size_t, size_t> frag_graph( m_GraphConverter( *itr_mol));
        if( GraphVectorFind( frag_graph, graphs) == graphs.GetSize())
        {
          graphs.PushBack( frag_graph);
        }
      }
      return graphs;
    }

    //! @brief determine if a graph is contained in a vector of graphs, i.e. same size and an isomorphism exists
    //! @param GRAPH the graph to check for
    //! @param GRAPH_VECTOR the vector of graphs to check
    //! @return the index of the first graph that matches, or GRAPH_VECTOR.GetSize() if not found
    size_t ConformationComparisonByFragments::GraphVectorFind
    ( 
      const graph::ConstGraph< size_t, size_t> &GRAPH, 
      const storage::Vector< graph::ConstGraph< size_t, size_t> > &GRAPH_VECTOR
    ) const
    {
      graph::SubgraphIsomorphism< size_t, size_t> iso_search;

      size_t matched_index( GRAPH_VECTOR.GetSize());

      iso_search.SetGraphExternalOwnership( GRAPH);
      size_t graph_size( GRAPH.GetSize()), graph_n_edges( GRAPH.NumEdges());

      // check the query graph against all graphs in the vector for uniqueness
      for( size_t graph_no( 0), end_no( matched_index); graph_no < end_no; ++graph_no)
      {
        // if the graph has a different number of vertices than the one being compared then
        // there is no way it can be equivalent
        if
        ( 
          graph_size != GRAPH_VECTOR( graph_no).GetSize() 
          || graph_n_edges != GRAPH_VECTOR( graph_no).NumEdges()
        )
        {
          continue;
        }

        iso_search.SetSubgraphExternalOwnership( GRAPH_VECTOR( graph_no));

        // if an isomorphism is found and its size matches the whole molecule then it is a match
        if( iso_search.FindIsomorphism() && iso_search.GetIsomorphism().GetSize() == graph_size)
        {
          matched_index = graph_no;
          break;
        }
      }
      
      return matched_index;
    }

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write out errors to
    bool ConformationComparisonByFragments::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      m_GraphConverter = ConformationGraphConverter( m_AtomComparison, m_BondComparison);
      return true;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ConformationComparisonByFragments::GetSerializer() const
    {
      io::Serializer parameters;
      const std::string measure( m_ComputeTanimoto ? "the tanimoto coefficient of the (#shared / #total)" : "");
      std::string desc( "Computes " + measure + "the number of common unique fragments between molecules");

      parameters.SetClassDescription( desc);

      parameters.AddInitializer
      (
        "",
        "the fragmentation methods to use",
        io::Serialization::GetAgent( &m_Splits)
      );

      parameters.AddInitializer
      (
        "atom comparison",
        "atom data that is compared to determine whether atoms are equivalent",
        io::Serialization::GetAgent( &m_AtomComparison),
        "ElementType"
      );
      parameters.AddInitializer
      (
        "bond comparison",
        "bond data that is compared",
        io::Serialization::GetAgent( &m_BondComparison),
        "BondOrderInRingOrAromatic"
      );

      return parameters;
    }

  } // namespace chemistry
} // namespace bcl
