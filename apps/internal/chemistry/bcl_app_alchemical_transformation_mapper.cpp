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

// include headers from the bcl - sorted alphabetically
#include "app/bcl_app_groups.h"
#include "chemistry/bcl_chemistry_atom_complete.h"
#include "chemistry/bcl_chemistry_atom_vector.h"
#include "chemistry/bcl_chemistry_atoms_complete_standardizer.h"
#include "chemistry/bcl_chemistry_configurational_bond_type_data.h"
#include "chemistry/bcl_chemistry_conformation_comparison_psi_field.h"
#include "chemistry/bcl_chemistry_conformation_graph_converter.h"
#include "chemistry/bcl_chemistry_constitution_set.h"
#include "chemistry/bcl_chemistry_fragment_complete.h"
#include "chemistry/bcl_chemistry_fragment_constitution_shared.h"
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "chemistry/bcl_chemistry_fragment_feed.h"
#include "chemistry/bcl_chemistry_hydrogens_handler.h"
#include "chemistry/bcl_chemistry_rotamer_library_interface.h"
#include "chemistry/bcl_chemistry_search_fragment_library_from_tree.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "graph/bcl_graph_common_subgraph_isomorphism_base.h"
#include "graph/bcl_graph_connectivity.h"
#include "graph/bcl_graph_const_graph.h"
#include "graph/bcl_graph_edge_cover_ring_perception.h"
#include "graph/bcl_graph_path.h"
#include "graph/bcl_graph_subgraph.h"
#include "graph/bcl_graph_subgraph_isomorphism.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_vector.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_string_functions.h"
namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AlchemicalTransformationMapper
    //! @brief Application for identifying core scaffold from substituted molecules
    //!
    //! @author brownbp1
    //! @date Oct 23, 2019
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class AlchemicalTransformationMapper :
      public Interface
    {
    private:

    //////////
    // data //
    //////////

      // flag specifying maximum number of edges between vertices
      util::ShPtr< command::FlagInterface> m_MaxCycleSize;

      // flag specifying maximum number of edges between vertices
      util::ShPtr< command::FlagInterface> m_MinEdgesBetweenVertices;
      
      // flag specifying reference molecule to which component pieces of the graph will be connected
      util::ShPtr< command::FlagInterface> m_ConnectToReferenceMol;

      // flag indicating if the reverse reactions are also desired
      util::ShPtr< command::FlagInterface> m_PrintReverse;

      // flag indicating if independent atom mask numbering is desired
      util::ShPtr< command::FlagInterface> m_IndependentMaskIndices;

      // flag build alchemical transformation graph based on similarity matrix
      util::ShPtr< command::FlagInterface> m_BuildGraph;

      // flag for input SDF files
      util::ShPtr< command::FlagInterface> m_InputFilenamesFlag;

      // flag for molecule id prefix
      util::ShPtr< command::FlagInterface> m_MolPrefixFlag;

      // flag for substructure connectivity
      util::ShPtr< command::FlagInterface> m_ConnectedFlag;

      // flag for scaffold-based alignment of pairs
      util::ShPtr< command::FlagInterface> m_AlignPairs;

      // flag for output file
      util::ShPtr< command::FlagInterface> m_OutputPrefixFlag;

      // flag for input similarity matrix
      util::ShPtr< command::FlagInterface> m_SimilarityMatrix;

      // flag for input ddG csv file
      util::ShPtr< command::FlagInterface> m_DDGFile;

      //! Substructure connection type
      mutable graph::CommonSubgraphIsomorphismBase::SolutionType m_SolutionType;

      //! Graphs of all molecules
      mutable storage::Map< util::SiPtr< const chemistry::FragmentComplete>, graph::ConstGraph< size_t, size_t> > m_Graphs;

    //////////////////////////////////
    // Construction and destruction //
    //////////////////////////////////

    public:

      //! @brief default constructor
      AlchemicalTransformationMapper();

      AlchemicalTransformationMapper *Clone() const
      {
        return new AlchemicalTransformationMapper( *this);
      }

      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize command to be returned
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // insert all the flags and params
        sp_cmd->AddFlag( m_MaxCycleSize);
        sp_cmd->AddFlag( m_MinEdgesBetweenVertices);
        sp_cmd->AddFlag( m_ConnectToReferenceMol);
        sp_cmd->AddFlag( m_PrintReverse);
        sp_cmd->AddFlag( m_IndependentMaskIndices);
        sp_cmd->AddFlag( m_BuildGraph);
        sp_cmd->AddFlag( m_InputFilenamesFlag);
        sp_cmd->AddFlag( m_MolPrefixFlag);
        sp_cmd->AddFlag( m_ConnectedFlag);
        sp_cmd->AddFlag( m_AlignPairs);
        sp_cmd->AddFlag( m_OutputPrefixFlag);
        sp_cmd->AddFlag( m_SimilarityMatrix);
        sp_cmd->AddFlag( m_DDGFile);

        // default flags
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

      graph::ConstGraph< std::string, double> BuildTransformationGraph() const
      {
        // read in similarity matrix
        linal::Matrix< double> distance_matrix;
        io::IFStream file;
        io::File::MustOpenIFStream( file, m_SimilarityMatrix->GetFirstParameter()->GetValue());
        io::Serialize::Read( distance_matrix, file);
        io::File::CloseClearFStream( file);

        // initialize vertex data from matrix rows
        storage::Vector< std::string> vertex_mol_ids;
        std::string prefix( m_MolPrefixFlag->GetFirstParameter()->GetValue());
        for( size_t i( 0); i < distance_matrix.GetNumberRows(); ++i)
        {
          vertex_mol_ids.PushBack( util::Format()( prefix) + "_" + util::Format()( i));
        }

        // initialize graph and similarity cutoff
        graph::ConstGraph< std::string, double> alchemy_graph
        (
          vertex_mol_ids,
          linal::Matrix< double>( vertex_mol_ids.GetSize(), vertex_mol_ids.GetSize(), double( 0.0)),
          double( 0.0)
        );

        // get upper triangle edges for each vertex
        storage::Vector< storage::Vector< storage::Pair< size_t, double >>> per_vertex_possible_edges;
        for( size_t row( 0); row < distance_matrix.GetNumberRows(); ++row)
        {
          storage::Vector< storage::Pair< size_t, double >> columns;
          for( size_t col( 0); col < distance_matrix.GetNumberCols(); ++col)
          {
            if( row == col)
            {
              continue;
            }
            columns.PushBack( std::make_pair( col, distance_matrix[ row][ col]));
          }
          per_vertex_possible_edges.PushBack( columns);
        }

        // sort the possible edges for each vertex
        static struct ComparePairs
        {
          bool operator()( const storage::Pair< size_t, double> &FIRST, const storage::Pair< size_t, double> &SECOND) const
          {
            return FIRST.Second() > SECOND.Second();
          }
        } pair_comparer;
        for( size_t row_index( 0); row_index < per_vertex_possible_edges.GetSize(); ++row_index)
        {
          // sort
          per_vertex_possible_edges( row_index).Sort< ComparePairs>( pair_comparer);
          size_t to_add( 0);
          while( alchemy_graph.GetNeighborIndices( row_index).GetSize() < size_t( m_MinEdgesBetweenVertices->GetFirstParameter()->GetNumericalValue< size_t>()))
          {
            alchemy_graph.AddEdge
            (
              row_index,
              per_vertex_possible_edges( row_index)( to_add).First(),
              per_vertex_possible_edges( row_index)( to_add).Second()
            );
            ++to_add;
          }
        }

        // now we need to connect all the components of our graph
        storage::List< storage::Vector< size_t> > components( graph::Connectivity::GetComponents( alchemy_graph));

        // connect all components to nearest reference vertex
        // having multiple reference vertices is very helpful when multiple experimental values are known
        storage::Vector< storage::Pair< size_t, size_t>> ref_connections;
        if( m_ConnectToReferenceMol->GetFlag())
        {
          // get reference vertex
          //          size_t ref_node_index( m_ConnectToReferenceMol->GetFirstParameter()->GetNumericalValue< size_t>());
          //          storage::Triplet< size_t, size_t, double> best_edge_score(0,0,0);
          storage::Vector< size_t> ref_vertices( m_ConnectToReferenceMol->GetNumericalList< size_t>());

          // go over each component
          for
          (
              auto all_comp_itr( components.Begin()), all_comp_itr_end( components.End());
              all_comp_itr != all_comp_itr_end;
              ++all_comp_itr
          )
          {
            // go over each vertex in component
            // track best possible connection to reference node for each vertex in this component
            storage::Triplet< size_t, size_t, double> best_edge_score( 0, 0, 0);
            for
            (
                auto vertex_itr( all_comp_itr->Begin()), vertex_itr_end( all_comp_itr->End());
                vertex_itr != vertex_itr_end;
                ++vertex_itr
            )
            {
              // go over each possible reference node
              for
              (
                  auto itr_ref_node( ref_vertices.Begin()), itr_ref_node_end( ref_vertices.End());
                  itr_ref_node != itr_ref_node_end;
                  ++itr_ref_node
              )
              {
                // if this component contains the reference node then don't do anything
                if( *vertex_itr == *itr_ref_node)
                {
                  break;
                }

                if( distance_matrix( *itr_ref_node, *vertex_itr) > best_edge_score.Third())
                {
                  best_edge_score = storage::Triplet< size_t, size_t, double>( *itr_ref_node, *vertex_itr, distance_matrix( *itr_ref_node, *vertex_itr));
                }
              } // end reference nodes
            } // end vertices in component

            // add an edge to the reference node
            alchemy_graph.AddEdge
            (
              best_edge_score.First(),
              best_edge_score.Second(),
              best_edge_score.Third()
            );

            // keep track so we can add these edges to the necessary pairs list
            if( best_edge_score.First() != best_edge_score.Second())
            {
              ref_connections.PushBack( std::make_pair( best_edge_score.First(), best_edge_score.Second()));
            }
          } // end components
        }

        // connect all components by similarity
        while( components.GetSize() > 1)
        {
          storage::Vector< size_t> &first( components.FirstElement());
          storage::Vector< size_t> &second( *( ++components.Begin()));
          storage::Triplet< size_t, size_t, double> best_edge_score( 0, 0, 0);
          for( size_t i( 0), sz( first.GetSize()); i < sz; ++i)
          {
            for( size_t j( 0), jsz( second.GetSize()); j < jsz; ++j)
            {
              if( distance_matrix( first( i), second( j)) > best_edge_score.Third())
              {
                best_edge_score = storage::Triplet< size_t, size_t, double>( first( i), second( j), distance_matrix( first( i), second( j)));
              }
            }
          }
          second.Append( first);
          alchemy_graph.AddEdge
          (
            best_edge_score.First(),
            best_edge_score.Second(),
            best_edge_score.Third()
          );
          components.PopFront();
        }

        // make sure shortest path to full cycle is < N
        bool did_something( true);
        size_t max_cycle_size( m_MaxCycleSize->GetFirstParameter()->GetNumericalValue< size_t>());
        while( did_something)
        {
          did_something = false;
          auto rings( graph::EdgeCoverRingPerception( alchemy_graph).GetRings());
          for
          (
              auto rings_itr( rings.Begin()), rings_itr_end( rings.End()); rings_itr != rings_itr_end; ++rings_itr
          )
          {
            if( rings_itr->GetSize() > size_t( max_cycle_size))
            {
              const storage::Vector< size_t> &first( rings_itr->GetVertices());
              storage::Triplet< size_t, size_t, double> best_edge_score( 0, 0, 0);
              for( size_t i( 0), sz( first.GetSize()); i < sz; ++i)
              {
                for( size_t j( 0); j < sz; ++j)
                {
                  if( j == i)
                  {
                    continue;
                  }
                  if( distance_matrix( first( i), first( j)) > best_edge_score.Third() && !alchemy_graph.AreConnected( first( i), first( j)))
                  {
                    best_edge_score = storage::Triplet< size_t, size_t, double>( first( i), first( j), distance_matrix( first( i), first( j)));
                  }
                }
              }
              alchemy_graph.AddEdge( best_edge_score.First(), best_edge_score.Second(), best_edge_score.Third());

              did_something = true;
              break; // prevent edges from crossing
            }
          }

          // make sure every vertex is in a ring
          for
          (
              size_t graph_index( 0), graph_index_end( alchemy_graph.GetVertices().GetSize());
              graph_index < graph_index_end;
              ++graph_index
          )
          {
            if
            (
                graph::Connectivity::LengthOfSmallestCycleWithVertex( alchemy_graph, graph_index) >= graph_index_end
            )
            {
              for
              (
                  size_t vertices( 0), vertices_end( per_vertex_possible_edges( graph_index).GetSize());
                  vertices < vertices_end;
                  ++vertices
              )
              {
                if( !alchemy_graph.AreConnected( graph_index, per_vertex_possible_edges( graph_index)( vertices).First()))
                {
                  alchemy_graph.AddEdge
                  (
                    graph_index,
                    per_vertex_possible_edges( graph_index)( vertices).First(),
                    per_vertex_possible_edges( graph_index)( vertices).Second()
                  );
                  if( graph::Connectivity::LengthOfSmallestCycleWithVertex( alchemy_graph, graph_index) <= m_MaxCycleSize->GetFirstParameter()->GetNumericalValue< size_t>())
                  {
                    break;
                  }
                  else
                  {
                    alchemy_graph.RemoveEdge
                    (
                      graph_index,
                      per_vertex_possible_edges( graph_index)( vertices).First()
                    );
                  }
                }
              }
            }
          }

          // output molecule pairs in each ring
          io::OFStream output;
          io::File::MustOpenOFStream( output, m_OutputPrefixFlag->GetFirstParameter()->GetValue() + ".pairs.txt");
          if( !did_something)
          {
            util::Format formatter;
            formatter.W( 2).Fill( '0');
            for
            (
                auto ring_itr( rings.Begin()), ring_itr_end( rings.End()); ring_itr != ring_itr_end; ++ring_itr
            )
            {
              size_t ring_index( 0);
              for
              (
                  auto itr( ring_itr->Begin());
                  ring_index < ring_itr->GetSize();
                  ++ring_index
              )
              {
                if( ring_index < ring_itr->GetSize() - 1)
                {
                  // print forward reaction
                  output << formatter(*itr) << ',';
                  ++itr;
                  output << formatter(*itr) << '\n';

                  // print reverse reactions
                  if( m_PrintReverse->GetFlag())
                  {
                    output << formatter(*itr) << ',';
                    --itr;
                    output << formatter(*itr) << '\n';

                    // fix itr
                    ++itr;
                  }
                }
                else
                {
                  // print forward reaction
                  output << formatter( *itr) << ',';
                  output << formatter( *( ring_itr->Begin())) << '\n';

                  // print reverse reaction
                  if( m_PrintReverse->GetFlag())
                  {
                    output << formatter( *( ring_itr->Begin())) << ',';
                    output << formatter( *itr) << '\n';
                  }
                }
              }
            }

            // add reference connections
            if( ref_connections.GetSize())
            {
              for
              (
                  auto ref_con_itr( ref_connections.Begin()), ref_con_itr_end( ref_connections.End());
                  ref_con_itr != ref_con_itr_end;
                  ++ref_con_itr
              )
              {
                // print forward reaction
                output << formatter( ref_con_itr->First()) << ',' << formatter( ref_con_itr->Second()) << '\n';

                // print reverse reactions
                if( m_PrintReverse->GetFlag())
                {
                  output << formatter( ref_con_itr->Second()) << ',' << formatter( ref_con_itr->First()) << '\n';
                }
              }
            }

            // close output stream
            io::File::CloseClearFStream( output);
            if( util::GetMessenger().GetCurrentMessageLevel() == util::Message::e_Debug)
            {
              BCL_Debug( rings);
            }
          }
        }

        // output final graph
        io::OFStream output;
        io::File::MustOpenOFStream( output, m_OutputPrefixFlag->GetFirstParameter()->GetValue() + ".graph.txt");
        output << alchemy_graph.GetBasicConnectivity();
        io::File::CloseClearFStream( output);
        if( util::GetMessenger().GetCurrentMessageLevel() == util::Message::e_Debug)
        {
          BCL_Debug( alchemy_graph.GetBasicConnectivity());
        }

        // output unique pairs
        io::File::MustOpenOFStream( output, m_OutputPrefixFlag->GetFirstParameter()->GetValue() + ".unique_pairs.txt");
        for
        (
            size_t graph_index( 0), graph_index_end( alchemy_graph.GetVertices().GetSize());
            graph_index < graph_index_end;
            ++graph_index
        )
        {
          const storage::Vector< size_t> &nbrs( alchemy_graph.GetNeighborIndices( graph_index));
          for
          (
              auto nbr_itr( nbrs.Begin()), nbr_itr_end( nbrs.End());
              nbr_itr != nbr_itr_end;
              ++nbr_itr
          )
          {
            if( graph_index > *nbr_itr || !alchemy_graph.AreConnected( *nbr_itr, graph_index))
            {
              output << graph_index << "," << *nbr_itr << "\n";
            }
          }
        }

        io::File::CloseClearFStream( output);

        // return our final graph
        return alchemy_graph;
      }

      void AlignPairs
      (
        const storage::Vector< chemistry::FragmentComplete> &MOLECULES
      ) const
      {
        // read in molecule pair file
        io::IFStream read;
        const std::string pair_filename( m_OutputPrefixFlag->GetFirstParameter()->GetValue() + ".unique_pairs.txt");
        io::File::MustOpenIFStream( read, pair_filename);
        const storage::Vector< std::string> pair_list( util::StringListFromIStream( read));
        io::File::CloseClearFStream( read);

        // split string by commas
        storage::Vector< storage::Vector< std::string>> string_pairs;
        for( size_t ddg_index( 0), sz( pair_list.GetSize()); ddg_index < sz; ++ddg_index)
        {
          string_pairs.PushBack( util::SplitString( pair_list( ddg_index), ","));
        }
        BCL_Debug( string_pairs);
        BCL_Debug( string_pairs( 0));
        BCL_Debug( string_pairs( 0)( 0));

        // align each pair
        for( size_t pair_index( 0), pair_index_end( string_pairs.GetSize()); pair_index < pair_index_end; ++pair_index)
        {
          size_t index_a( util::ConvertStringToNumericalValue< size_t>(string_pairs( pair_index)( 0)));
          size_t index_b( util::ConvertStringToNumericalValue< size_t>(string_pairs( pair_index)( 1)));
          chemistry::FragmentComplete mol_a(MOLECULES( index_a));
          chemistry::FragmentComplete mol_b(MOLECULES( index_b));

          static chemistry::ConformationComparisonPsiField aligner;
          aligner.ArbitraryScaffoldAlignment
          (
            mol_a,
            mol_b,
            static_cast< chemistry::ConformationGraphConverter::AtomComparisonType>( 2),
            static_cast< chemistry::ConfigurationalBondTypeData::Data>( 11)
          );

          // output each pair into an SDF file
          io::OFStream out;
          io::File::MustOpenOFStream( out, string_pairs( pair_index)( 0) + "_" + string_pairs( pair_index)( 1) + ".sdf", std::ios::app);
          mol_a.WriteMDL( out);
          mol_b.WriteMDL( out);
          io::File::CloseClearFStream( out);
        }
      }

      void Analyze() const
      {
        // read in ddG file
        io::IFStream read;
        const std::string ddg_filename( m_DDGFile->GetFirstParameter()->GetValue());
        io::File::MustOpenIFStream( read, ddg_filename);
        const storage::Vector< std::string> ddg_list( util::StringListFromIStream( read));
        io::File::CloseClearFStream( read);

        // split string by commas
        storage::Vector< storage::Vector< std::string>> split_strings;
        for( size_t ddg_index( 0), sz( ddg_list.GetSize()); ddg_index < sz; ++ddg_index)
        {
          split_strings.PushBack( util::SplitString( ddg_list( ddg_index), ","));
        }

        // get ddG values as floats
        storage::Map< storage::Pair< std::string, std::string>, double> ddg_values;
        for( size_t ddg_index( 0), sz( split_strings.GetSize()); ddg_index < sz; ++ddg_index)
        {
          storage::Pair< std::string, std::string> key( std::make_pair( split_strings( ddg_index)( 0), split_strings( ddg_index)( 1)));
          ddg_values.Insert( std::make_pair( key, util::ConvertStringToNumericalValue< double>( split_strings( ddg_index)( 2))));
        }

        // initialize vertex names
        storage::Set< std::string> vertex_mol_ids;
        storage::Map< std::string, size_t> index_mapper;
        storage::Vector< std::string> vertices;
        size_t index( 0);
        std::string prefix( m_MolPrefixFlag->GetFirstParameter()->GetValue());
        for( auto itr( split_strings.Begin()), itr_end( split_strings.End()); itr != itr_end; ++itr)
        {
          // our string names may not match the index of the set, so we need to map them back
          if( vertex_mol_ids.Insert( util::Format()( *( itr->Begin()))).second)
          {
            index_mapper[ ( *itr)( 0)] = index;
            ++index;
            vertices.PushBack( ( *itr)( 0));
          }
          if( vertex_mol_ids.Insert( util::Format()( *( itr->Begin() + 1))).second)
          {
            index_mapper[ ( *itr)( 1)] = index;
            ++index;
            vertices.PushBack( ( *itr)( 1));
          }
          std::cout << "Adding vertices: " << util::Format()( *( itr->Begin())) << " and " << util::Format()( *( itr->Begin() + 1)) << "\n";
        }
        BCL_Debug( vertices);

        // initialize graph and similarity cutoff
        graph::ConstGraph< std::string, double> alchemy_graph
        (
          vertices,
          linal::Matrix< double>( vertices.GetSize(), vertices.GetSize(), double( 12345.6789)),
          double( 12345.6789),
          true,
          false
        );
        graph::ConstGraph< std::string, double> alchemy_graph_unweighted
        (
          vertices,
          linal::Matrix< double>( vertices.GetSize(), vertices.GetSize(), double( 0.0)),
          double( 0.0),
          false,
          true
        );
        // add ddG values as edges
        for
        (
            storage::Map< storage::Pair< std::string, std::string>, double>::iterator edge_itr( ddg_values.Begin()), edge_itr_end( ddg_values.End());
            edge_itr != edge_itr_end;
            ++edge_itr
        )
        {
          // this is the value from the output of the GPU-TI calculation
          alchemy_graph.AddEdge
          (
            index_mapper.Find( edge_itr->first.First())->second,
            index_mapper.Find( edge_itr->first.Second())->second,
            edge_itr->second
          );
          // this is the reverse of the output value to make the opposite direction edge; theoretically this should be true
          alchemy_graph.AddEdge
          (
            index_mapper.Find( edge_itr->first.Second())->second,
            index_mapper.Find( edge_itr->first.First())->second,
            -edge_itr->second
          );
          // this is the value from the output of the GPU-TI calculation
          alchemy_graph_unweighted.AddEdge
          (
            index_mapper.Find( edge_itr->first.First())->second,
            index_mapper.Find( edge_itr->first.Second())->second,
            1.0
          );
        }

        // We want to find the shortest path by connectivity, not by value
        // until we have a FindShortestPath function instead of FindMinimalPath,
        // we need to make a "shadow graph" with all of the same edges but set to a constant
        // value. Then, we need to set the graph to directed so that it maintains fidelity to ordering
        // the vertex connections on output
        alchemy_graph_unweighted.SetDirected();

        // find path to nearest reference vertex
        if( m_ConnectToReferenceMol->GetFlag())
        {
          // read in the file with the reference nodes
          const std::string ref_dg_filename( m_ConnectToReferenceMol->GetFirstParameter()->GetValue());
          io::File::MustOpenIFStream( read, ref_dg_filename);
          const storage::Vector< std::string> dg_list( util::StringListFromIStream( read));
          io::File::CloseClearFStream( read);

          // split string by commas
          storage::Vector< storage::Vector< std::string>> split_strings;
          for( size_t dg_index( 0), sz( dg_list.GetSize()); dg_index < sz; ++dg_index)
          {
            split_strings.PushBack( util::SplitString( dg_list( dg_index), ","));
          }

          // get reference dG values as floats
          storage::Map< std::string, double> dg_values;
          storage::Vector< std::string> ref_vertices;
          for( size_t dg_index( 0), sz( split_strings.GetSize()); dg_index < sz; ++dg_index)
          {
            std::string key( split_strings( dg_index)( 0));
            BCL_Debug( key);
            ref_vertices.PushBack( key);
            dg_values.Insert( std::make_pair( key, util::ConvertStringToNumericalValue< double>( split_strings( dg_index)( 1))));
          }

          // find vertices that user fallaciously included
          storage::Set< std::string> ref_set( ref_vertices.Begin(), ref_vertices.End());
          storage::Set< std::string> fake_ref( ref_set - vertex_mol_ids);
          BCL_MessageStd( "Provided reference vertices are not in graph and will be excluded: " + util::Format()( fake_ref));

          // remove fake vertices from ref_set
          ref_set = ref_set - fake_ref;
          for
          (
              auto fake_ref_itr( fake_ref.Begin()), fake_ref_itr_end( fake_ref.End());
              fake_ref_itr != fake_ref_itr_end;
              ++fake_ref_itr
          )
          {
            dg_values.Erase( *fake_ref_itr);
          }

          // go over each vertex
          for
          (
              storage::Vector< std::string>::const_iterator vertex_itr( alchemy_graph.GetVertices().Begin()), vertex_itr_end( alchemy_graph.GetVertices().End());
              vertex_itr != vertex_itr_end;
              ++vertex_itr
          )
          {
            // go over each reference vertex and get path length
            graph::Path shortest_path;
            size_t shortest_path_length( math::GetHighestBoundedValue< size_t>());
            for
            (
                storage::Set< std::string>::const_iterator ref_itr( ref_set.Begin()), ref_itr_end( ref_set.End());
                ref_itr != ref_itr_end;
                ++ref_itr
            )
            {
              // should add an assert statement to make sure reference is actually in graph
//              BCL_Debug( *vertex_itr);
//              BCL_Debug( index_mapper.Find(*vertex_itr)->second);
//              BCL_Debug( *ref_itr);
//              BCL_Debug( index_mapper.Find(*ref_itr)->second);
              graph::Path path
              (
                graph::Connectivity::FindMinimalPath
                (
                  alchemy_graph_unweighted,
                  index_mapper.Find( *vertex_itr)->second,
                  index_mapper.Find( *ref_itr)->second
                )
              );
//              BCL_Debug( path.GetSize());
//              BCL_Debug( shortest_path.GetSize());
//              BCL_MessageStd("\n");
              if( path.GetSize() < shortest_path_length && path.GetSize() > 1)
              {
                shortest_path = path;
//                BCL_Debug( shortest_path.GetVertices());
                shortest_path_length = shortest_path.GetSize();
              }
            } // end reference vertex itr

            // compute the relative ddG to the nearest reference node
            double ddg( 0.0);
            if( shortest_path.GetSize())
            {
              // go over the shortest path
              for
              (
                  auto path_itr_follow( shortest_path.GetVertices().Begin()), path_itr_lead( shortest_path.GetVertices().Begin() + 1), path_itr_end( shortest_path.GetVertices().End());
                  path_itr_lead != path_itr_end;
                  ++path_itr_follow, ++path_itr_lead
              )
              {
                // find edge value
                ddg += alchemy_graph.GetEdgeData( *path_itr_follow, *path_itr_lead);
              }

              // compute predicted dG values of non-reference vertices
              std::string current_target( alchemy_graph.GetVertexData( *( shortest_path.GetVertices().Begin())));
              std::string current_ref( alchemy_graph.GetVertexData( *( --shortest_path.GetVertices().End())));
              double dg_ref( dg_values.Find( current_ref)->second);
              double dg_exp( dg_ref - ddg);
              BCL_MessageStd( "Ref: " + util::Format()(current_ref) + " dG = " + util::Format()( dg_values.Find( current_ref)->second) + " kcal/mol and, "
                  "Exp: " + util::Format()( current_target) + " dG = " + util::Format()( dg_exp) + " kcal/mol");
              BCL_MessageStd( "ddG: " + util::Format()( ddg));
            }
          } // end graph vertex itr
        }
      }

      void WriteAtomMask
      (
        const graph::ConstGraph< std::string, double> &GRAPH,
        const storage::Vector< chemistry::FragmentComplete> &MOLECULES
      ) const
      {
        // formatter
        util::Format formatter;
        formatter.W( 2).Fill( '0');

        // get all closed cycles in the graph
        auto rings( graph::EdgeCoverRingPerception( GRAPH).GetRings());

        // setup output
        io::OFStream output;
        io::File::MustOpenOFStream( output, m_OutputPrefixFlag->GetFirstParameter()->GetValue() + ".atom_masks.txt");

        // go over closed cycle
        for
        (
            auto ring_itr( rings.Begin()), ring_itr_end( rings.End()); ring_itr != ring_itr_end; ++ring_itr
        )
        {
          // iterate over each molecule pair in closed cycle; only need to do forward direction
          size_t ring_index( 0);
          for
          (
              auto itr( ring_itr->Begin());
              ring_index < ring_itr->GetSize();
              ++ring_index
          )
          {
            // indicate start of new molecule pair atom mask section
            output << "START ";

            chemistry::FragmentComplete mol_a, mol_b;
            // get molecule at this index and print molecule IDs to atom mask file
            if( ring_index < ring_itr->GetSize() - 1)
            {
              mol_a = MOLECULES( *itr);
              output << formatter( *itr) << '_';
              ++itr;
              mol_b = MOLECULES( *itr);
              output << formatter( *itr) << '\n';
            }
            else
            {
              mol_a = MOLECULES( *itr);
              output << formatter( *itr) << '_';
              mol_b = MOLECULES( *(ring_itr->Begin()));
              output << formatter( *( ring_itr->Begin())) << '\n';
            }

            // prepare dehydrogenated copy
            chemistry::FragmentComplete mol_a_noh( mol_a);
            chemistry::FragmentComplete mol_b_noh( mol_b);
            mol_a_noh.RemoveH();
            mol_b_noh.RemoveH();

            // prepare graphs
            Prepare( mol_a);
            Prepare( mol_b);
            Prepare( mol_a_noh);
            Prepare( mol_b_noh);

            // get graphs for the molecules
            storage::Map< util::SiPtr< const chemistry::FragmentComplete>, graph::ConstGraph< size_t, size_t> >::iterator
            itr_graph_a( m_Graphs.Find( mol_a)), itr_graph_b( m_Graphs.Find( mol_b)), itr_graph_a_noh( m_Graphs.Find( mol_a_noh)), itr_graph_b_noh( m_Graphs.Find( mol_b_noh));

            // set up the isomorphisms with pointers to the graphs
            graph::CommonSubgraphIsomorphism< size_t, size_t> isomorphism( m_SolutionType);
            isomorphism.SetGraphA
            (
              util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &itr_graph_a_noh->second, false)
            );
            isomorphism.SetGraphB
            (
              util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &itr_graph_b_noh->second, false)
            );

            // get the isomorphism
            isomorphism.FindIsomorphism( isomorphism.EstimateUpperBounds());

            // get inverted subgraph of a
            auto vertices_a
            (
              m_ConnectedFlag->GetFlag()
              ? isomorphism.GetIsomorphism().GetKeysAsVector()
              : isomorphism.GetSubgraphIsomorphismsOfGraphA().FirstElement().GetIdsOfInteriorVertices()
            );
            graph::Subgraph< size_t, size_t> subgraph_a
            (
              util::OwnPtr< graph::ConstGraph< size_t, size_t> >( &itr_graph_a_noh->second, false),
              vertices_a
            );

            // create an isomorphism with hydrogen atoms included
            graph::SubgraphIsomorphism< size_t, size_t> proper_isomorphism_a, proper_isomorphism_b;

            // do A first
            proper_isomorphism_a.SetSubgraph( subgraph_a.ToGraph());
            proper_isomorphism_a.SetGraph( itr_graph_a->second);
            proper_isomorphism_a.FindIsomorphism();
            graph::Subgraph< size_t, size_t> subgraph_a_with_h
            (
              proper_isomorphism_a.GetSubgraphIsomorphisms().FirstElement()
            );

            // do B first
            proper_isomorphism_a.SetSubgraph( subgraph_a.ToGraph());
            proper_isomorphism_a.SetGraph( itr_graph_b->second);
            proper_isomorphism_a.FindIsomorphism();
            graph::Subgraph< size_t, size_t> subgraph_b_with_h
            (
              proper_isomorphism_a.GetSubgraphIsomorphisms().FirstElement()
            );
            auto edges_a_outtie( subgraph_a_with_h.GetOrderedAdjacentEdgeIndices());
            auto edges_b_outtie( subgraph_b_with_h.GetOrderedAdjacentEdgeIndices());

            storage::Vector< size_t> subgraph_a_with_h_indices( subgraph_a_with_h.GetVertexIndices());
            storage::Vector< size_t> subgraph_b_with_h_indices( subgraph_b_with_h.GetVertexIndices());
            for( size_t common_index( 0), ncom( edges_a_outtie.GetSize()); common_index < ncom; ++common_index)
            {
              // A
              storage::List< size_t> h_a_indices, heavy_a_indices;
              for
              (
                auto itr_edges_a_outtie( edges_a_outtie( common_index).Begin()), itr_edges_a_outtie_end( edges_a_outtie( common_index).End());
                itr_edges_a_outtie != itr_edges_a_outtie_end;
                ++itr_edges_a_outtie
              )
              {
                if( mol_a.GetAtomVector()( *itr_edges_a_outtie).GetElementType()->GetAtomicNumber() == size_t( 1))
                {
                  h_a_indices.PushBack( *itr_edges_a_outtie);
                }
                else
                {
                  heavy_a_indices.PushBack( *itr_edges_a_outtie);
                }
              }

              // B
              storage::List< size_t> h_b_indices, heavy_b_indices;
              for
              (
                auto itr_edges_b_outtie( edges_b_outtie( common_index).Begin()), itr_edges_b_outtie_end( edges_b_outtie( common_index).End());
                itr_edges_b_outtie != itr_edges_b_outtie_end;
                ++itr_edges_b_outtie
              )
              {
                if( mol_b.GetAtomVector()( *itr_edges_b_outtie).GetElementType()->GetAtomicNumber() == size_t( 1))
                {
                  h_b_indices.PushBack( *itr_edges_b_outtie);
                }
                else
                {
                  heavy_b_indices.PushBack( *itr_edges_b_outtie);
                }
              }

              while( !h_a_indices.IsEmpty() && !h_b_indices.IsEmpty())
              {
                subgraph_a_with_h_indices.PushBack( h_a_indices.FirstElement());
                subgraph_b_with_h_indices.PushBack( h_b_indices.FirstElement());
                h_a_indices.PopFront();
                h_b_indices.PopFront();
              }
//              // if they differ only in terms of hydrogen atoms or heavy atoms, skip
//              if( h_a_indices.IsEmpty() && h_b_indices.IsEmpty())
//              {
//                continue;
//              }
//              if( heavy_a_indices.IsEmpty() && heavy_b_indices.IsEmpty())
//              {
//                continue;
//              }
//              if( !h_a_indices.IsEmpty() && !heavy_b_indices.IsEmpty())
//              {
//                // hydrogens left over on A, which were replaced by heavy atoms on B. Need to include the hydrogen atom on the subgraph
//              }
//              else // if( h_b_indices.IsEmpty())
//              {
//
//              }
            }

            // redo isomorphism
            subgraph_a_with_h.SetVertexIndices( subgraph_a_with_h_indices);
            subgraph_b_with_h.SetVertexIndices( subgraph_b_with_h_indices);

            // get inverted subgraph of b
            graph::Subgraph< size_t, size_t> complement_subgraph_a( subgraph_a_with_h.GetComplement());
            graph::Subgraph< size_t, size_t> complement_subgraph_b( subgraph_b_with_h.GetComplement());

            // output atom masks in Amber format
            // make an array of the masks needed
            std::string atom_masks[ 4] =
            {
                "timask1='@",
                "scmask1='@",
                "timask2='@",
                "scmask2='@"
            };
            for( size_t mask_index( 0); mask_index < size_t( 4); ++mask_index)
            {
              output << util::Format()( atom_masks[ mask_index]);
              // output from molecule A
              if( mask_index < size_t( 2))
              {

                for
                (
                    auto itr_a( complement_subgraph_a.GetVertexIndices().Begin()), itr_a_end( complement_subgraph_a.GetVertexIndices().End());
                    itr_a != itr_a_end;
                    ++itr_a
                )
                {
                  // syntactical shit
                  if( itr_a == complement_subgraph_a.GetVertexIndices().Begin())
                  {
                    output << util::Format()( *itr_a + 1); // atom indices are 1-indexed in most other places
                  }
                  else
                  {
                    output << "," + util::Format()( *itr_a + 1);
                  }
                }
                output << "'," << '\n';
              }
              // output from molecule B
              else
              {
                size_t offset( mol_a.GetSize());
                if( m_IndependentMaskIndices->GetFlag())
                {
                  offset = size_t( 0);
                }
                for
                (
                    auto itr_b( complement_subgraph_b.GetVertexIndices().Begin()), itr_b_end( complement_subgraph_b.GetVertexIndices().End());
                    itr_b != itr_b_end;
                    ++itr_b
                )
                {
                  // syntactical shit
                  if( itr_b == complement_subgraph_b.GetVertexIndices().Begin())
                  {
                    output << util::Format()( *itr_b + 1 + offset);
                  }
                  else
                  {
                    output << "," + util::Format()( *itr_b + 1 + offset);
                  }
                }
                output << "'," << '\n';
              }
            }
            output << "END" << '\n';

            // print debug statements
            if( util::GetMessenger().GetCurrentMessageLevel() == util::Message::e_Debug)
            {
              BCL_Debug( complement_subgraph_a.GetVertexIndices());
              BCL_Debug( complement_subgraph_b.GetVertexIndices());
            }
          } // done with pair in closed cycle
        } // done with all cycles
      }

      int Main() const
      {
        if( m_ConnectedFlag->GetFlag())
        {
          m_SolutionType = graph::CommonSubgraphIsomorphismBase::SolutionType::e_Connected;
        }

        // build alchemical transformation graph based on similarity matrix
        graph::ConstGraph< std::string, double> similarity_graph;
        if( m_BuildGraph->GetFlag() && m_SimilarityMatrix->GetFlag())
        {
          similarity_graph = BuildTransformationGraph();
        }

        // create atom masks based on substructure similarity
        // use of atom masks is a more precise way to make the alchemical transformation
        if( m_BuildGraph->GetFlag() && m_InputFilenamesFlag->GetFlag() && similarity_graph.GetSize())
        {
          // read in the molecules
          storage::Vector< chemistry::FragmentComplete> input_mols;
          chemistry::FragmentFeed mol_feed( m_InputFilenamesFlag->GetStringList(), sdf::e_Maintain);
          for( ; mol_feed.NotAtEnd(); ++mol_feed)
          {
            const chemistry::FragmentComplete &mol( *mol_feed);
            input_mols.PushBack( mol);
          }

          // input molecules are assumed to be in the same order as the molecules in the similarity matrix
          // this allows an easy mapping between graph vertices and fragment ensemble indices
          // it also allows flexibility in the type of similarity matrix generated (at the risk of fucking everything with unmatched orders)
          // minor check to prevent obvious error
          BCL_Assert( input_mols.GetSize() == similarity_graph.GetSize(),
            "There are " + util::Format()( input_mols.GetSize()) + " input molecules, but " + util::Format()( similarity_graph.GetSize()) + " vertices in graph!");

          // align molecule pairs
          if( m_AlignPairs->GetFlag())
          {
            AlignPairs( input_mols);
          }

          // generate atom masks
          WriteAtomMask( similarity_graph, input_mols);

        }

        // write unique pairs for simpler system setup

        // analyze ddG to compute affinities relative to reference molecule
        // provide error estimate for every ddG via cycle closure hysteresis
        if( m_DDGFile->GetFlag())
        {
          Analyze();
        }

        return 0;
      } // Main()

    protected:

      //! @brief prepare the class for comparing a conformation
      //! @param MOLECULE the molecule to prepare to compare
      void Prepare( const chemistry::FragmentComplete &MOLECULE) const
      {
        // create a graph converter
        chemistry::ConformationGraphConverter graph_maker
        (
          chemistry::ConformationGraphConverter::AtomComparisonType::e_AtomType,
          chemistry::ConfigurationalBondTypeData::Data::e_BondOrderAmideOrAromaticWithRingness,
          false
        );

        // make the graph and append it
        m_Graphs[ MOLECULE] = graph_maker( MOLECULE);
      }

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    private:

      static const ApplicationType AlchemicalTransformationMapper_Instance;

    }; // AlchemicalTransformationMapper

      //! @brief standard constructor
    AlchemicalTransformationMapper::AlchemicalTransformationMapper() :
      m_MaxCycleSize
      (
        new command::FlagStatic
        (
          "max_cycle_size",
          "maximum number of edges between any set of three vertices in a cycle",
          command::Parameter
          (
            "maximum cycle size",
            "integer specifying max number of edges per cycle",
            command::ParameterCheckRanged< size_t>( 3, 1000),
            "4"
          )
        )
      ),
      m_MinEdgesBetweenVertices
      (
        new command::FlagStatic
        (
          "min_mol_edges",
          "minimum number of edges from each vertex",
          command::Parameter
          (
            "minimum molecule edges",
            "integer specifying min number of edges from each vertex",
            command::ParameterCheckRanged< size_t>( 1, 1000),
            "2"
          )
        )
      ),
      m_MolPrefixFlag
      (
        new command::FlagStatic
        (
          "prefix",
          "prefix for molecule IDs",
          command::Parameter
          (
            "prefix",
            "prefix for molecule IDs",
            ""
          )
        )
      ),
      m_PrintReverse
      (
        new command::FlagStatic
        (
          "include_reverse_reactions",
          "output file will also have reverse reactions",
          command::Parameter
          (
            "include reverse reactions",
            "useful if the equilibrated pocket of each of the two molecules substantially differs",
            ""
          )
        )
      ),
      m_IndependentMaskIndices
      (
        new command::FlagStatic
        (
          "independent_mask_indices",
          "output the atom masks with indices corresponding to just the individual molecules",
          command::Parameter
          (
            "default behavior adds to the total number of atoms in mol_a to the beginning of the indices for the mol_b atom mask;"
            "this is because it is smart to make the two ligands molecules 1 and 2 in the PDB file for a TI calculation, and the indices are"
            "continuous",
            "if independent numbering is desired (e.g. start from 1 for each mol) then enable this flag",
            ""
          )
        )
      ),
      m_BuildGraph
      (
        new command::FlagStatic
        (
          "build_graph",
          "build graph connecting molecules based on similarity matrix",
          command::Parameter
          (
            "build graph",
            "ensures that every molecule is enclosed in at least one cycle and limits size of any cycle to no greater than max_cycle_size edges",
            ""
          )
        )
      ),
      m_InputFilenamesFlag
      (
        new command::FlagDynamic
        (
          "input_filenames",
          "filename of an sdf file containing the molecules to use",
          command::Parameter
          (
            "file containing the molecules",
            "name of files containing structures",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_OutputPrefixFlag
      (
        new command::FlagStatic
        (
          "output", "flag selecting the output file name",
          command::Parameter
          (
            "filename", "filename for output sdf", "alchemical_transformation_mapper.default.txt"
          )
        )
      ),
      m_SimilarityMatrix
      (
        new command::FlagStatic
        (
          "similarity_matrix",
          "similarity matrix between molecules",
          command::Parameter
          (
            "similarity matrix",
            "similarity matrix between molecules generated with molecule:Compare",
            ""
          )
        )
      ),
      m_DDGFile
      (
        new command::FlagStatic
        (
          "ddg_file",
          "ddG values for each transformation in graph",
          command::Parameter
          (
            "ddg file",
            "comma-separated file indicating molecule_1,molecule_2,ddG for all pairs",
            ""
          )
        )
      ),
      m_ConnectToReferenceMol
      (
        new command::FlagStatic
        (
          "connect_to_reference_mol",
          "connect the separate components of the graph to a reference node",
          command::Parameter
          (
            "reference file",
            "file containing reference data such that each row is has two comma-separated columns; "
            "the first column is the mol ID, and the second column is the dG in kcal/mol",
            ""
          )
        )
      ),
      m_ConnectedFlag
      (
        new command::FlagStatic
        (
          "connected",
          "maximum common connected substructure comparison for atom mask generation",
          command::Parameter
          (
            "connected",
            "maximum common connected substructure comparison for atom mask generation",
            ""
          )
        )
      ),
      m_AlignPairs
      (
        new command::FlagStatic
        (
          "align_pairs",
          "perform align to scaffold of the second molecule in each pair to the first",
          command::Parameter
          (
            "align_pairs",
            "perform align to scaffold of the second molecule in each pair to the first",
            ""
          )
        )
      ),
      m_SolutionType( graph::CommonSubgraphIsomorphismBase::SolutionType::e_Unconnected)
    {
    }

    const ApplicationType AlchemicalTransformationMapper::AlchemicalTransformationMapper_Instance
    (
      GetAppGroups().AddAppToGroup( new AlchemicalTransformationMapper(), GetAppGroups().e_Molecule)
    );

  } // namespace app
} // namespace bcl
