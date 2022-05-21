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

#ifndef BCL_CLUSTER_OUTPUT_ROWS_H_
#define BCL_CLUSTER_OUTPUT_ROWS_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_cluster_node.h"
#include "bcl_cluster_output_interface.h"
#include "io/bcl_io_file.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class OutputRows
    //! @brief is a type of Output interface which outputs the results of a Dendrogram in rows.
    //! @details Each cluster in the dendrogram is output to a row, so all the member for a given cluster are on the same
    //! line. Also, the distance between the two clusters is output.
    //!
    //! @see @link example_cluster_output_rows.cpp @endlink
    //! @author alexanns
    //! @date October 9, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class OutputRows :
      public OutputInterface< t_DataType, t_PrecisionType>
    {
    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Default constructor
      OutputRows()
      {
      }

      //! @brief copy constructor
      OutputRows< t_DataType, t_PrecisionType> *Clone() const
      {
        return new OutputRows< t_DataType, t_PrecisionType>( *this);
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

    ///////////////
    // operators //
    ///////////////

      //! @brief WriteOutput writes some data of a list of Nodes to some output stream
      //! @param FILENAME is the stream to which the output will be written
      //! @param NODE_LIST the list of nodes which is going to be output
      //! @return returns data construct holding the distances between all objects for use in a LinkageInterface
      void WriteOutput
      (
        const std::string &FILENAME, const util::SiPtrList< const Node< t_DataType, t_PrecisionType> > &NODE_LIST
      ) const
      {
        io::OFStream ofstream;
        io::File::MustOpenOFStream( ofstream, FILENAME);

        for
        (
          typename util::SiPtrList< const Node< t_DataType, t_PrecisionType> >::const_iterator
            node_list_itr( NODE_LIST.Begin()), node_list_itr_end( NODE_LIST.End());
          node_list_itr != node_list_itr_end;
          ++node_list_itr
        )
        {
          util::SiPtrList< const Node< t_DataType, t_PrecisionType> > expanded_node_list( ( *node_list_itr)->ExpandAllNodes());

          // iterate through the nodes of "NODE_LIST" so that the members of each can be output
          for
          (
            typename util::SiPtrList< const Node< t_DataType, t_PrecisionType> >::const_iterator
              node_itr( expanded_node_list.Begin()), node_itr_end( expanded_node_list.End());
            node_itr != node_itr_end;
            ++node_itr
          )
          {
            // Output the members of the node currently denoted by "node_itr"
            OutputNodeMembers( ofstream, **node_itr);
          }
        }
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

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief OutputNodeMembers outputs the individual members of a Node
      //! @param OFSTREAM the output stream to which the Node will be written
      //! @param NODE the node which will be output to OFSTREAM
      //! @param NODE_COUNTER the number of the current node
      //! @param OUTER_NODE_COUNTER
      //! @return returns io::OFStream
      virtual std::ostream &OutputNodeMembers
      (
        std::ostream &OFSTREAM, const Node< t_DataType, t_PrecisionType> &NODE
      ) const
      {
        // create SiPtrList "member_list" and initialize with the members of "NODE"
        const util::SiPtrList< t_DataType> member_list( NODE.GetMembers());

        bool leaf( NODE.IsLeaf());

        // iterate through "cluster_members" in order to output each one
        for
        (
          typename util::SiPtrList< t_DataType>::const_iterator
            member_itr( member_list.Begin()), member_itr_end( member_list.End());
          member_itr != member_itr_end;
          ++member_itr
        )
        {
          // output the node number and the girth of the current node denoted by "itr"
          OFSTREAM << "NODE ";
          // write out the identifier; it might be nan, so use io::Serialize
          io::Serialize::Write( NODE.GetIdentifier(), OFSTREAM, 0);
          OFSTREAM << " : Member : " << std::string( util::Format()( **member_itr)) << " : Size : "
          << NODE.GetMembers().GetSize() << " : Leaf : " << util::Format()( leaf) << " : Linkage : " <<
          util::Format()( NODE.GetGirth()) << '\n';
        }

        // end
        return OFSTREAM;
      }

    }; // OutputRows

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_OUTPUT_ROWS_H_ 
