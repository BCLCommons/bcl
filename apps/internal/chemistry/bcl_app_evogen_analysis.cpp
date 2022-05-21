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
#include "app/bcl_app_apps.h"
#include "app/bcl_app_groups.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "graph/bcl_graph_const_graph.h"
#include "graph/bcl_graph_undirected_edge.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_running_min_max.h"
#include "sdf/bcl_sdf_fragment_factory.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_table.hpp"
#include "storage/bcl_storage_table_header.h"
#include "util/bcl_util_message.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_si_ptr.h"
#include "util/bcl_util_string_numeric_conversion.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically
#include <cstdio>

namespace bcl
{
  namespace app
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EvoGenAnalysis
    //! @brief Application for analyzing output of the EvoGen application
    //!
    //! @author geanesar
    //! @date Dec 10 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class EvoGenAnalysis :
      public Interface
    {

    //////////
    // data //
    //////////

      //! the filename of any histories that should be read in
      util::ShPtr< command::FlagInterface> m_HistoriesFilenameFlag;

      //! the filename to write the DOT file to
      util::ShPtr< command::FlagInterface> m_DotOutputFilenameFlag;

      //! where images are stored
      util::ShPtr< command::FlagInterface> m_ImagesDirectoryFlag;

      //! a list of members to look at
      util::ShPtr< command::FlagInterface> m_MembersFlag;

      //! the prefix for member-specific dot output
      util::ShPtr< command::FlagInterface> m_MemberDotOutputPrefixFlag;

      //! the prefix for member-specific info
      util::ShPtr< command::FlagInterface> m_MemberInfoOutputPrefixFlag;

      //! statistics output information
      util::ShPtr< command::FlagInterface> m_RunStatisticsFlag;

    ////////////////////
    // helper classes //
    ////////////////////

      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //!
      //! @class GraphNode
      //! @brief A small class used to represent EvoGen history when in a graph format
      //!
      //! @author geanesar
      //! @date Dec 10 2014
      //!
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      class GraphNode :
        public util::ObjectInterface
      {

      private:

      //////////
      // data //
      //////////

          //! the name of the object
          std::string m_Name;

          //! the fitness of the object
          std::string m_Fitness;

          //! population and member numbers
          storage::Pair< size_t, size_t> m_PopAndMemberNumbers;

          //! the names of the parents of the object
          storage::Vector< std::string> m_ParentNames;

          storage::Vector< storage::Pair< size_t, size_t> > m_ParentValues;

          //! the operation used to get to this object
          std::string m_Operation;

          //! @brief sets up a GraphNode object from an ObjectDataLabel
          //! @param LABEL the label to use
          void SetUpObjectFromLabel( const util::ObjectDataLabel &LABEL)
          {
            if( LABEL.IsEmpty() || LABEL.GetValue().empty())
            {
              return;
            }

            util::ObjectDataLabel::const_iterator itr_end( LABEL.End());

            m_Name = LABEL.GetValue();
            m_PopAndMemberNumbers = GetPopAndMemberNumbersFromString( LABEL.GetValue());

            util::ObjectDataLabel::const_iterator itr_fitness( LABEL.FindName( "fitness"));
            if( itr_fitness != itr_end)
            {
              m_Fitness = itr_fitness->GetValue();
            }

            util::ObjectDataLabel::const_iterator itr_parents( LABEL.FindName( "parents"));
            if( itr_parents != itr_end)
            {
              const std::vector< util::ObjectDataLabel> &parent_labels( itr_parents->GetArguments());
              for
              (
                std::vector< util::ObjectDataLabel>::const_iterator itr_lab( parent_labels.begin()), itr_lab_end( parent_labels.end());
                itr_lab != itr_lab_end;
                ++itr_lab
              )
              {
                m_ParentValues.PushBack( GetPopAndMemberNumbersFromString( itr_lab->GetValue()));
              }
            }

            util::ObjectDataLabel::const_iterator itr_op( LABEL.FindName( "operation"));
            if( itr_op != itr_end)
            {
              m_Operation = itr_op->GetValue();
            }

          }

      public:

          static storage::Pair< size_t, size_t> GetPopAndMemberNumbersFromString( const std::string &STRING)
          {
            std::stringstream sp;
            std::stringstream sm;
            int state( 0);
            for( std::string::const_iterator itr_str( STRING.begin()), itr_str_end( STRING.end()); itr_str != itr_str_end; ++itr_str)
            {
              if( state == 0) // beginning of string, first letter must be P
              {
                BCL_Assert( *itr_str == 'P', "Object value \"" + STRING + "\" does not fit the PXMY format (character " + *itr_str + ")");
                state = 1;
                continue;
              }
              else if( state == 1) // population number
              {
                if( isdigit( *itr_str))
                {
                  sp << *itr_str;
                }
                else
                {
                  BCL_Assert( *itr_str == 'M', "Object value \"" + STRING + "\" does not fit the PXMY format (character " + *itr_str + ")");
                  state = 2;
                }
              }
              else if( state == 2)
              {
                if( isdigit( *itr_str))
                {
                  sm << *itr_str;
                }
                else
                {
                  break;
                }
              }
            }
            int pop_val;
            int member_val;
            BCL_Assert( util::TryConvertFromString( pop_val, sp.str(), util::GetLogger()), "Cannot convert the population value to a string");
            BCL_Assert( util::TryConvertFromString( member_val, sm.str(), util::GetLogger()), "Cannot convert the member value to a string");
            return storage::Pair< size_t, size_t>( pop_val, member_val);
          }

          //! @brief default constructor
          GraphNode()
          {
          }

          //! @brief constructor with label
          //! @param LABEL the label to convert into a GraphNode
          GraphNode( const util::ObjectDataLabel &LABEL)
          {
            SetUpObjectFromLabel( LABEL);
          }

          //! @brief clone operation
          //! @return a copy of this object
          GraphNode *Clone() const
          {
            return new GraphNode( *this);
          }

          //! @brief get the class name
          //! @return the class name
          const std::string &GetClassIdentifier() const
          {
            return GetStaticClassName( this);
          }

          //! @brief get the name of the node
          //! @return the name of the node
          const std::string &GetName() const
          {
            return m_Name;
          }

          //! @brief get the fitness of the node
          //! @return a string containing the fitness of the node
          const std::string &GetFitness() const
          {
            return m_Fitness;
          }

          //! @brief get the parents
          //! @return a vector of parent names
          const storage::Vector< std::string> &GetParentNames() const
          {
            return m_ParentNames;
          }

          const storage::Vector< storage::Pair< size_t, size_t> > &GetParentValues() const
          {
            return m_ParentValues;
          }

          const size_t &GetPopulationNumber() const
          {
            return m_PopAndMemberNumbers.First();
          }

          const size_t &GetMemberNumber() const
          {
            return m_PopAndMemberNumbers.Second();
          }

          //! @brief get the operation
          //! @return the name of the operation
          const std::string &GetOperation() const
          {
            return m_Operation;
          }

          //! @brief checks if this node has a parent
          //! @return true if the parent exists in this node
          bool HasParent( const std::string &PARENT_NAME) const
          {
            for
            (
              storage::Vector< std::string>::const_iterator itr_parent( m_ParentNames.Begin()),
                itr_parent_end( m_ParentNames.End());
              itr_parent != itr_parent_end;
              ++itr_parent
            )
            {
              if( *itr_parent == PARENT_NAME)
              {
                return true;
              }
            }
            return false;
          }

          //! @brief Gets the population of a member
          //! @return the population a member belongs to
          size_t GetPopulation() const
          {
            std::stringstream pop_stream;
            for( size_t i( 0); i < m_Name.length(); ++i)
            {
              // Skip starting characters
              for( ; !isdigit( m_Name[ i]); ++i);

              // Output the first numerical part (the population) to a string
              for( ; isdigit( m_Name[ i]); ++i)
              {
                pop_stream << m_Name[ i];
              }
              break;
            }
            return util::ConvertStringToNumericalValue< size_t>( pop_stream.str());
          }

          //! @brief gets a numerical version of the fitness
          //! @return numerical encoding of the fitness
          double GetNumericalFitness() const
          {
            double member_fitness;
            util::TryConvertFromString( member_fitness, m_Fitness, util::GetLogger());
            return member_fitness;
          }

      protected:

          //! @brief read from std::istream
          //! @param ISTREAM input stream
          //! @return istream which was read from
          std::istream &Read( std::istream &ISTREAM)
          {
            util::ObjectDataLabel new_label( ISTREAM);
            SetUpObjectFromLabel( new_label);
            return ISTREAM;
          }

          //! @brief write to std::ostream
          //! @param OSTREAM output stream to write to
          //! @param INDENT - number of indentations
          //! @return output stream which was written to
          std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
          {
            return OSTREAM;
          }

      };

    ////////////////
    // operations //
    ////////////////

      //! @brief gets the available colors to use in a dot diagram
      //! @return a vector of strings specifying colors
      const storage::Vector< std::string> &GetAvailableColors() const
      {
        static storage::Vector< std::string> s_colors;
        s_colors.PushBack( "black");
        s_colors.PushBack( "red");
        s_colors.PushBack( "blue");
        s_colors.PushBack( "green");
        s_colors.PushBack( "orange");
        s_colors.PushBack( "purple");
        s_colors.PushBack( "gray");
        return s_colors;
      }

      //! @brief locates a SiPtr to a graph node with a particular name
      //! @param NODE_VECTOR the node vector to search
      //! @param NAME the name to look for
      //! @return the index of the first instance of that node
      size_t FindIndexOfNodeWithName( const util::SiPtrVector< GraphNode> &NODE_VECTOR, const std::string &NAME) const
      {
        size_t i( 0);
        for( ; i < NODE_VECTOR.GetSize(); ++i)
        {
          if( NODE_VECTOR( i)->GetName() == NAME)
          {
            break;
          }
        }
        return i;
      }

      //! @brief convert a full graph to a dot chart
      //! @param GRAPH the graph to convert
      //! @param GRAPH_NAME the name of the graph
      //! @param OSTREAM where to write the graph to
      //! @return returns OSTREAM
      std::ostream &GraphToDOT
      (
        const graph::ConstGraph< util::SiPtr< const GraphNode>, int> &GRAPH,
        const std::string &GRAPH_NAME,
        std::ostream &OSTREAM,
        const std::string &IMAGE_DIR = std::string()
      ) const
      {

        //const storage::Vector< std::string> &colors( GetAvailableColors());

        // Get the vertices of the graph
        const storage::Vector< util::SiPtr< const GraphNode> > &vertices( GRAPH.GetVertices());

        // The deader part of the thing
        OSTREAM << "digraph " + GRAPH_NAME + " {" << std::endl;

        OSTREAM << "  node [shape=box];" << std::endl;

        // add vertices/connections
        for( size_t pos( 0); pos < vertices.GetSize() - 1; ++pos)
        {
          storage::Vector< size_t> neighbors( GRAPH.GetNeighborIndices( pos));
          for
          (
            storage::Vector< size_t>::const_iterator itr_neighbor( neighbors.Begin()), itr_neighbor_end( neighbors.End());
            itr_neighbor != itr_neighbor_end;
            ++itr_neighbor
          )
          {
            if( *itr_neighbor <= pos)
            {
              continue;
            }
            //int edge_data( GRAPH.GetEdgeData( pos, *itr_neighbor));
            //const std::string &color( edge_data < colors.GetSize() ? colors( edge_data) : "gray");
            OSTREAM << "  " << vertices( pos)->GetName() << " -> " << vertices( *itr_neighbor)->GetName();
            if( !vertices( *itr_neighbor)->GetOperation().empty())
            {
              OSTREAM << " [ label=\"" + vertices( *itr_neighbor)->GetOperation() + "\" ]";
            }
            OSTREAM << std::endl;
                    //<< " [color=" + color + "];" << std::endl;

          }
        }

        // Add vertices to the graph
        for( size_t pos( 0); pos < vertices.GetSize(); ++pos)
        {
          if( IMAGE_DIR.empty())
          {
            OSTREAM << "  " << vertices( pos)->GetName() << " [label=\"" << vertices( pos)->GetName() << "\\n"
                    << "fitness=" << vertices( pos)->GetFitness() << "\"];" << std::endl;
          }
          else
          {
            std::string img_name( IMAGE_DIR + "/" + vertices( pos)->GetName() + ".png");
            OSTREAM << "  " << vertices( pos)->GetName() << " [label=<<TABLE>";
            OSTREAM << "<TR><TD>" + vertices( pos)->GetName() + "</TD></TR>";
            if( !vertices( pos)->GetFitness().empty())
            {
              OSTREAM << "<TR><TD>fitness: " + vertices( pos)->GetFitness() + "</TD></TR>";
            }
            OSTREAM << "<TR><TD><IMG SRC=\"" + img_name + "\"/></TD></TR></TABLE>>]" << std::endl;
          }
        }

        // close the graph
        OSTREAM << "}" << std::endl;
        return OSTREAM;
      }

      void GetNodeLineage
      (
        const linal::Matrix< int> &EDGE_MATRIX,
        const size_t &CURRENT_INDEX,
        storage::Set< size_t> &CURRENT_LINEAGE
      ) const
      {
        //const size_t &rows( EDGE_MATRIX.GetNumberRows());
        const size_t &cols( EDGE_MATRIX.GetNumberCols());

        if( !CURRENT_LINEAGE.Insert( CURRENT_INDEX).second)
        {
          return;
        }

        storage::Set< size_t> parents;
        for( size_t c( 0); c < cols; ++c)
        {
          if( EDGE_MATRIX( CURRENT_INDEX, c) < 0)
          {
            parents.Insert( c);
          }
        }

        for
        (
          storage::Set< size_t>::const_iterator itr_parent( parents.Begin()), itr_parent_end( parents.End());
          itr_parent != itr_parent_end;
          ++itr_parent
        )
        {
          GetNodeLineage( EDGE_MATRIX, *itr_parent, CURRENT_LINEAGE);
        }
      }

      //! @brief recurses through a node's lineage and gets the parent indices associated with it
      //! @param GRAPH the graph to use
      //! @param GRAPH_INDEX the graph index to look at
      //! @param INDEX_VECTOR the vector to add the found indices to
      //! @return a vector containing the lineage of
      void RecursivelyGetLineageOfNode
      (
        const graph::ConstGraph< util::SiPtr< GraphNode>, size_t> &GRAPH,
        const size_t &GRAPH_INDEX,
        storage::Set< size_t> &INDEX_VECTOR
      ) const
      {

        const util::SiPtr< GraphNode> &node( GRAPH.GetVertexData( GRAPH_INDEX));
        const storage::Vector< size_t> &neighbors( GRAPH.GetNeighborIndices( GRAPH_INDEX));

        // Determine the parent indices
        storage::Vector< size_t> parent_indices;
        parent_indices.AllocateMemory( neighbors.GetSize());
        for( size_t i( 0); i < neighbors.GetSize(); ++i)
        {
          if( node->HasParent( GRAPH.GetVertexData( neighbors( i))->GetName()))
          {
            parent_indices.PushBack( neighbors( i));
            INDEX_VECTOR.Insert( neighbors( i));
          }
        }

        // Recurse through the parents
        for( size_t i( 0); i < parent_indices.GetSize(); ++i)
        {
          RecursivelyGetLineageOfNode( GRAPH, parent_indices( i), INDEX_VECTOR);
        }
      }

      //! @brief recursively gets indices of a subgraph up to a given label
      //! @param GRAPH the whole graph to look in
      //! @param NODE the node to trace the lineage of
      //! @return a vector of indices that make up the subgraph
      storage::Vector< size_t> GetLineageOfNode
      (
        const graph::ConstGraph< util::SiPtr< GraphNode>, size_t> &GRAPH,
        const size_t &NODE
      ) const
      {
        if( NODE > GRAPH.GetSize())
        {
          BCL_MessageStd( "NODE was outside of GRAPH");
          return storage::Vector< size_t>();
        }

        // set of indices
        storage::Set< size_t> include_indices;
        // Insert the desired node
        include_indices.Insert( NODE);
        RecursivelyGetLineageOfNode( GRAPH, NODE, include_indices);

        return storage::Vector< size_t>( include_indices.Begin(), include_indices.End());
      }

      void CalculateStatisticsOfTransitions
      (
        const graph::ConstGraph< util::SiPtr< GraphNode>, size_t> &GRAPH,
        const storage::Map< std::string, size_t> &STRING_TO_EDGE_MAP
      ) const
      {
        util::SiPtrVector< GraphNode> vertices( GRAPH.GetVertices().Begin(), GRAPH.GetVertices().End());

        // A map of edge type (key) to the total number of edges of that type (values)
        storage::Map< size_t, size_t> edge_total_count_map;

        // A map of edge type (key) to the
        storage::Map< size_t, size_t> edge_improved_count_map;

        for( size_t index( 0); index < vertices.GetSize(); ++index)
        {
          double member_fitness( vertices( index)->GetNumericalFitness());
          for
          (
              storage::Vector< std::string>::const_iterator itr_parent( vertices( index)->GetParentNames().Begin()),
                itr_parent_end( vertices( index)->GetParentNames().End());
              itr_parent != itr_parent_end;
              ++itr_parent
          )
          {
            size_t parent_index( FindIndexOfNodeWithName( vertices, *itr_parent));
            if( parent_index >= vertices.GetSize())
            {
              continue;
            }
            double parent_fitness( vertices( parent_index)->GetNumericalFitness());
            size_t edge_type( GRAPH.GetEdgeData( index, parent_index));
            ++edge_total_count_map[ edge_type];
            if( parent_fitness <= member_fitness)
            {
              ++edge_improved_count_map[ edge_type];
            }
            else
            {
              // Make sure the entry is there
              if( !edge_improved_count_map.Has( edge_type))
              {
                edge_improved_count_map[ edge_type] = 0;
              }
            }
          }
        }

        // Make the table for outputting data
        storage::Vector< std::string> header_strings;
//        header_strings.PushBack( "Operation");
        header_strings.PushBack( "# Improved");
        header_strings.PushBack( "# Total");
        header_strings.PushBack( "Success rate");
        util::ShPtr< storage::TableHeader> sp_table_header( new storage::TableHeader( header_strings));
        storage::Table< std::string> edge_table( sp_table_header);

        if( edge_improved_count_map.GetSize() != edge_total_count_map.GetSize())
        {
          BCL_MessageStd( "ERROR with making tables");
          return;
        }

        BCL_MessageStd( "maps have " + util::Format()( edge_improved_count_map.GetSize()));

        // Build the table up
        storage::Map< size_t, size_t>::const_iterator itr_improved( edge_improved_count_map.Begin());
        for
        (
          storage::Map< size_t, size_t>::const_iterator itr_total( edge_total_count_map.Begin()), itr_total_end( edge_total_count_map.End());
          itr_total != itr_total_end;
          ++itr_total, ++itr_improved
        )
        {

          storage::Vector< std::string> row_info( 3);

          // Find the string for the operation type
          std::string row_name;
          for
          (
            storage::Map< std::string, size_t>::const_iterator itr_type( STRING_TO_EDGE_MAP.Begin()), itr_type_end( STRING_TO_EDGE_MAP.End());
            itr_type != itr_type_end;
            ++itr_type
          )
          {
            if( itr_type->second == itr_total->first)
            {
              row_name = itr_type->first;
              break;
            }
          }
          if( !row_name.length())
          {
            break;
          }
          row_info( 0) = util::Format()( itr_improved->second);
          row_info( 1) = util::Format()( itr_total->second);
          if( itr_total->second != 0)
          {
            row_info( 2) = util::Format()( double( itr_improved->second) / double( itr_total->second));
          }
          else
          {
            row_info( 2) = "nan";
          }
          edge_table.InsertRow( row_name, row_info);
        }
        edge_table.WriteFormatted( std::cout, util::Format().W( 10).R(), "Operation");
      }

    public:

    //////////////////////////////////
    // Construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      EvoGenAnalysis();

      //! @brief clone function - copies this object
      //! @return a pointer to a copy of this object
      EvoGenAnalysis *Clone() const
      {
        return new EvoGenAnalysis( *this);
      }

      //! @brief get the name of the class
      //! @return the name of the class
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief initializes the command object for that executable
      //! @return a ShPtr to a command object that contains command line flags for this app
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize command to be returned
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // insert all the flags and params
        sp_cmd->AddFlag( m_HistoriesFilenameFlag);
        sp_cmd->AddFlag( m_DotOutputFilenameFlag);
        sp_cmd->AddFlag( m_ImagesDirectoryFlag);
        sp_cmd->AddFlag( m_MembersFlag);
        sp_cmd->AddFlag( m_MemberDotOutputPrefixFlag);
        sp_cmd->AddFlag( m_MemberInfoOutputPrefixFlag);
        sp_cmd->AddFlag( m_RunStatisticsFlag);

        // default flags
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

      //! @brief the main function
      //! @return 0 if the function ran well
      int Main() const
      {

        io::IFStream input;
        io::File::MustOpenIFStream( input, m_HistoriesFilenameFlag->GetFirstParameter()->GetValue());

        storage::Vector< util::ObjectDataLabel> histories;
        storage::Vector< GraphNode> graph_nodes;

        while( !input.eof())
        {
          util::ObjectDataLabel new_label( input);
          histories.PushBack( new_label);
          GraphNode new_node( new_label);
          graph_nodes.PushBack( new_node);
        }

        storage::Map< size_t, storage::Map< size_t, size_t> > pop_member_map;
        util::SiPtrVector< const GraphNode> graph_node_ptrs;
        graph_node_ptrs.AllocateMemory( graph_nodes.GetSize());

        for
        (
          storage::Vector< GraphNode>::const_iterator itr_node( graph_nodes.Begin()), itr_node_end( graph_nodes.End());
          itr_node != itr_node_end;
          ++itr_node
        )
        {
          const size_t &pop( itr_node->GetPopulationNumber());
          const size_t &mem( itr_node->GetMemberNumber());
          size_t member_index( graph_node_ptrs.GetSize());
          if( pop_member_map.Has( pop))
          {
            if( pop_member_map.GetValue( pop).Has( mem))
            {
              BCL_MessageVrb( "Warning: P" + util::Format()( pop) + "M" + util::Format()( mem) + " has duplicate entry");
              member_index = pop_member_map[ pop][ mem];
            }
          }

          if( member_index == graph_node_ptrs.GetSize())
          {
            graph_node_ptrs.PushBack( util::SiPtr< const GraphNode>( &( *itr_node)));
          }

          pop_member_map[ pop][ mem] = member_index;
        }

        // Determine the connectivity of nodes
        linal::Matrix< int> connectivity_matrix( graph_node_ptrs.GetSize(), graph_node_ptrs.GetSize(), int( 0));

        for
        (
          storage::Map< size_t, storage::Map< size_t, size_t> >::const_iterator itr_map_1( pop_member_map.Begin()), itr_map_1_end( pop_member_map.End());
          itr_map_1 != itr_map_1_end;
          ++itr_map_1
        )
        {
          for
          (
            storage::Map< size_t, size_t>::const_iterator itr_map_2( itr_map_1->second.Begin()), itr_map_2_end( itr_map_1->second.End());
            itr_map_2 != itr_map_2_end;
            ++itr_map_2
          )
          {
            const size_t &pop( itr_map_1->first);
            const size_t &mem( itr_map_2->first);

            const size_t siptr_index( pop_member_map[ pop][ mem]);
            util::SiPtr< const GraphNode> &node_ptr( graph_node_ptrs( siptr_index));

            const storage::Vector< storage::Pair< size_t, size_t> > &parent_locations( node_ptr->GetParentValues());

            for( size_t i( 0), end_i( parent_locations.GetSize()); i < end_i; ++i)
            {
              const size_t &parent_pop( parent_locations( i).First());
              const size_t &parent_mem( parent_locations( i).Second());
              const size_t &parent_siptr_index( pop_member_map[ parent_pop][ parent_mem]);

              connectivity_matrix( parent_siptr_index, siptr_index) = 1;
              connectivity_matrix( siptr_index, parent_siptr_index) = -1;
            }
          }
        }

        graph::ConstGraph< util::SiPtr< const GraphNode>, int> relation_graph
        (
          graph_node_ptrs,
          connectivity_matrix,
          size_t( 0), // disconnected edge type
          true, // directed
          false // don't make undirected
        );

        if( m_DotOutputFilenameFlag->GetFlag())
        {
          // output the full dot graph
          std::string img_dir;
          if( m_ImagesDirectoryFlag->GetFlag())
          {
            img_dir = m_ImagesDirectoryFlag->GetFirstParameter()->GetValue();
          }
          io::OFStream dot_output;
          io::File::MustOpenOFStream( dot_output, m_DotOutputFilenameFlag->GetFirstParameter()->GetValue());
          BCL_MessageStd( "Writing full relation graph to " + m_DotOutputFilenameFlag->GetFirstParameter()->GetValue());
          GraphToDOT( relation_graph, "full_graph", dot_output, img_dir);
          io::File::CloseClearFStream( dot_output);
        }

        storage::Vector< std::string> members( m_MembersFlag->GetStringList());
        if( !members.IsEmpty() && m_MemberDotOutputPrefixFlag->GetFlag())
        {
          std::string prefix( m_MemberDotOutputPrefixFlag->GetFirstParameter()->GetValue());

          const linal::Matrix< int> &edge_matrix( relation_graph.GetEdgeDataMatrix());

          for( size_t m( 0), end_m( members.GetSize()); m < end_m; ++m)
          {
            storage::Pair< size_t, size_t> member( GraphNode::GetPopAndMemberNumbersFromString( members( m)));
            const size_t &mem_index( pop_member_map[ member.First()][ member.Second()]);
            storage::Set< size_t> vertices;
            GetNodeLineage( edge_matrix, mem_index, vertices);
            storage::Vector< size_t> vertices_v( vertices.Begin(), vertices.End());
            graph::ConstGraph< util::SiPtr< const GraphNode>, int> subgraph( relation_graph.GetSubgraph( vertices_v));

            std::string member_name( "P" + util::Format()( member.First()) + "M" + util::Format()( member.Second()));
            std::string out_filename( prefix + "_" + member_name + ".dot");

            io::OFStream dot_output;
            io::File::MustOpenOFStream( dot_output, out_filename);

            BCL_MessageStd( "Writing lineage of " + member_name + " to " + out_filename);

            std::string img_dir;
            if( m_ImagesDirectoryFlag->GetFlag())
            {
              img_dir = m_ImagesDirectoryFlag->GetFirstParameter()->GetValue();
            }

            GraphToDOT( subgraph, member_name, dot_output, img_dir);
          }
        }

        return 0;
      } // Main()

    protected:

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

      static const ApplicationType EvoGenAnalysis_Instance;

    }; // EvoGenAnalysis

      //! @brief standard constructor
    EvoGenAnalysis::EvoGenAnalysis() :
      m_HistoriesFilenameFlag
      (
        new command::FlagStatic
        (
          "histories",
          "the name of the file that contains the histories from an EvoGen run",
          command::Parameter
          (
            "histories",
            "a filename containing histories data labels",
            command::ParameterCheckFileExistence()
          )
        )
      ),
      m_DotOutputFilenameFlag
      (
        new command::FlagDynamic
        (
          "output_full",
          "the filename that the DOT graph for the full EvoGen run will be written to",
          command::Parameter
          (
            "output_full",
            "the filename that DOT output should be written to"
          )
        )
      ),
      m_ImagesDirectoryFlag
      (
        new command::FlagDynamic
        (
          "images_directory",
          "the directory that contains images of the molecules of interest, named PXMY.png (where X,Y are numbers)",
          command::Parameter
          (
            "directory",
            "the directory path"
          )
        )
      ),
      m_MembersFlag
      (
        new command::FlagDynamic
        (
          "members",
          "identifiers for members that should be analyzed (specified as PXMY)",
          command::Parameter
          (
            "members",
            "member identifiers (usually PXMY X,Y are numbers) of members that should be investigated"
          )
        )
      ),
      m_MemberDotOutputPrefixFlag
      (
        new command::FlagDynamic
        (
          "output_dot_prefix",
          "the prefix that will be used for member-specific lineage DOT graphs",
          command::Parameter
          (
            "output_dot_prefix",
            "the prefix of the filename that will be generated for member-specific DOT output."
            "  Output will be to <prefix>_<member name>.dot"
          )
        )
      ),
      m_MemberInfoOutputPrefixFlag
      (
        new command::FlagDynamic
        (
          "output_info_prefix",
          "the prefix that will be used to output member-specific information e.g. statistics",
          command::Parameter
          (
            "output_info_prefix",
            "the prefix that filenames will have.  output will be to <prefix>_<member name>.txt"
          )
        )
      ),
      m_RunStatisticsFlag
      (
        new command::FlagDynamic
        (
          "run_statistics",
          "where statistics for the overall run should be written to",
          command::Parameter
          (
            "run_statistics",
            "where statistics for the overall run should be written to"
          )
        )
      )
    {
    }

    const ApplicationType EvoGenAnalysis::EvoGenAnalysis_Instance
    (
      GetAppGroups().AddAppToGroup( new EvoGenAnalysis(), GetAppGroups().e_ChemInfo)
    );

  } // namespace app
} // namespace bcl
