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

#ifndef BCL_CLUSTER_NODE_H_
#define BCL_CLUSTER_NODE_H_

// include the namespace header
#include "bcl_cluster.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"
#include "storage/bcl_storage_set.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_binary_function_interface.h"
#include "util/bcl_util_sh_ptr_list.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Node
    //! @brief is used for holding a collection of objects which are related by one another within a
    //! threshold of some difference metric.
    //! @details For example, one node might have all protein models that are
    //! within 5angstrom RMSD to one another. A node knows about the nodes that comprise it. For example,
    //! the node that has all protein models that are within 5angstrom RMSD to one another might be made up
    //! of a node that has all protein models that are within 3angstrom RMSD and a node that has all protein
    //! models within 1anstrom RMSD to one another.
    //!
    //! @see @link example_cluster_node.cpp @endlink
    //! @author alexanns, woetzen
    //! @date May 29, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType, typename t_PrecisionType>
    class Node :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! member "m_Nodes" is a util::ShPtrList< Node> which holds all of the sub-Nodes that make up this current
      //! node
      util::ShPtrList< Node< t_DataType, t_PrecisionType> > m_Nodes;

      //! member "m_Members" is a util::SiPtrList< t_DataType> which holds all of the object identifiers that are
      //! members of this current Node
      util::SiPtrList< t_DataType> m_Members;

      //! member "m_Girth" is a t_PrecisionType which is internal girth threshold for this current Node
      t_PrecisionType m_Girth;

      //! numerical identifier of this node
      size_t m_Identifier;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Node() :
        m_Nodes(),
        m_Members(),
        m_Girth( util::GetUndefined< t_PrecisionType>()),
        m_Identifier( util::GetUndefinedSize_t())
      {
      }

      //! @brief constructor for a Node which is not made up of any other Nodes but has members
      Node( const util::SiPtrList< t_DataType> &MEMBERS, const size_t IDENTIFIER = util::GetUndefinedSize_t()) :
        m_Nodes(),
        m_Members( MEMBERS),
        m_Girth( util::GetUndefined< t_PrecisionType>()),
        m_Identifier( IDENTIFIER)
      {
      }

      //! @brief constructor for when Nodes are merged into a new cluster
      //! @param NODES util::ShPtrList< Node> of Nodes which will be merged
      //! @param GIRTH t_PrecisionType which is girth threshold the constructed Node will have
      Node( const util::ShPtrList< Node< t_DataType, t_PrecisionType> > &NODES, const t_PrecisionType &GIRTH, const size_t IDENTIFIER = util::GetUndefinedSize_t()) :
        m_Nodes( NODES),
        m_Members(),
        m_Girth( GIRTH),
        m_Identifier( IDENTIFIER)
      {
        // iterate through "NODES" and fill "m_Members" by appending the "m_Members" of the Nodes in "NODES"
        for
        (
          typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::const_iterator itr( NODES.Begin()), itr_end( NODES.End());
          itr != itr_end;
          ++itr
        )
        {
          // append the members of Node currently denoted by "itr" to "m_Members"
          m_Members.Append( ( *itr)->GetMembers());
        }
      }

      //! @brief constructor from all members
      //! @param NODES util::ShPtrList< Node> of Nodes which will be merged
      //! @param MEMBERS the members this node will be constructed from
      //! @param GIRTH t_PrecisionType which is girth threshold the constructed Node will have
      Node
      (
        const util::ShPtrList< Node< t_DataType, t_PrecisionType> > &NODES,
        const util::SiPtrList< t_DataType> &MEMBERS,
        const t_PrecisionType &GIRTH,
        const size_t IDENTIFIER = util::GetUndefinedSize_t()
      ) :
        m_Nodes( NODES),
        m_Members( MEMBERS),
        m_Girth( GIRTH),
        m_Identifier( IDENTIFIER)
      {
      }

      //! @brief Clone function
      //! @return pointer to new Node
      virtual Node< t_DataType, t_PrecisionType> *Clone() const
      {
        return new Node< t_DataType, t_PrecisionType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      virtual const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief GetNodes gives a const reference to "m_Nodes"
      //! @return a const reference to "m_Nodes"
      util::ShPtrList< Node< t_DataType, t_PrecisionType> > const &GetNodes() const
      {
        return m_Nodes;
      }

      //! @brief GetNodes gives a reference to "m_Nodes"
      //! @return a const reference to "m_Nodes"
      util::ShPtrList< Node< t_DataType, t_PrecisionType> > &GetNodes()
      {
        return m_Nodes;
      }

      //! @brief SetNodes changes "m_Nodes" to a different ShPtrList of Nodes
      //! @param NODES is the ShPtrList of Nodes which "m_Nodes" will be set to
      void SetNodes( const util::ShPtrList< Node< t_DataType, t_PrecisionType> > &NODES)
      {
        m_Nodes = NODES;
      }

      //! @brief GetMembers gives a const reference to "m_Members"
      //! @return returns a const reference to "m_Members"
      util::SiPtrList< t_DataType> const &GetMembers() const
      {
        return m_Members;
      }

      //! @brief GetMembers gives a reference to "m_Members"
      //! @return returns a const reference to "m_Members"
      util::SiPtrList< t_DataType> &GetMembers()
      {
        return m_Members;
      }

      //! @brief SetMembers changes "m_Members" to a new SiPtrList of strings
      //! @param MEMBERS is the SiPtrList of strings which "m_Members" will be changed to
      void SetMembers( const util::SiPtrList< t_DataType> &MEMBERS)
      {
        m_Members = MEMBERS;
      }

      //! @brief GetGirth gives a const reference to "m_Girth"
      //! @return returns a const reference to "m_Girth"
      const t_PrecisionType &GetGirth() const
      {
        return m_Girth;
      }

      //! @brief SetGirth changes "m_Girth" to a new t_PrecisionType
      //! @param GIRTH is the new t_PrecisionType which "m_Girth" will be changed to
      void SetGirth( const t_PrecisionType GIRTH)
      {
        m_Girth = GIRTH;
      }

      //! @brief gives "m_Identifier"
      //! @return returns "m_Identifier"
      size_t GetIdentifier() const
      {
        return m_Identifier;
      }

      //! @brief sets "m_Identifier" to a new value
      //! @param IDENTIFIER the value that "m_Identifier" will be set to
      void SetIdentifier( const size_t IDENTIFIER)
      {
        m_Identifier = IDENTIFIER;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief RemoveNodesBelowSize
      //! @param SIZE
      //! @return returns
      void RemoveNodesBelowSize( const size_t SIZE)
      {
        // iterate through "m_Nodes" to remove those nodes who are too small (too few members)
         for
         (
           typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator itr( m_Nodes.Begin()), itr_end( m_Nodes.End());
           itr != itr_end;

         )
         {
           if( ( *itr)->GetMembers().GetSize() < SIZE)
           {
             itr = m_Nodes.Remove( itr);
           }
           else
           {
             ( *itr)->RemoveNodesBelowSize( SIZE);
             ++itr;
           }
         }
         if( m_Nodes.GetSize() == size_t( 1) && GetMembers().GetSize() < SIZE)
         {
           m_Nodes.Reset();
         }
      }

      void RemoveSingularBranches()
      {
        for
        (
          typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
            itr( m_Nodes.Begin()), itr_end( m_Nodes.End());
          itr != itr_end;
          ++itr
        )
        {
          ( *itr)->RemoveSingularBranches();
        }

        if( m_Nodes.GetSize() == 1)
        {
          // this temporary variable is necessary otherwise a segfault will result because a ShPtr will go out of scope
          util::ShPtr< Node< t_DataType, t_PrecisionType> > tmp_node( m_Nodes.FirstElement());
          *this = *tmp_node;
        }
      }

      //! @brief RemoveNodesWithHighSimilarity removes inner nodes which have a high internal similarity
      //! @param GIRTH_CUTOFF the girth which is used as the cutoff point for removing nodes or not
      //! @param BINARY_PREDICATE the method that should be used to compare girths
      //!        (less = girth is a dissimilarity measure, ie RMSD of protein models)
      void RemoveNodesWithHighSimilarity
      (
        const t_PrecisionType GIRTH_CUTOFF,
        const util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> &BINARY_PREDICATE
      )
      {
        // remove all children nodes, if girth is undefined
        if( !util::IsDefined( m_Girth))
        {
          m_Nodes.Reset();
          BCL_MessageDbg
          (
            "my girth is undefined girth is " + util::Format()( m_Girth)
            + " m_nodes size after reset " + util::Format()( m_Nodes.GetSize())
          );

          // exit function
          return;
        }

        // remove all children nodes, if girth does not fullfill the cutoff
        if( BINARY_PREDICATE( m_Girth, GIRTH_CUTOFF))
        {
          m_Nodes.Reset();
          BCL_MessageDbg
          (
            "my girth does not fulfill cutoff girth " + util::Format()( m_Girth) + " cutoff " + util::Format()( GIRTH_CUTOFF)
            + " m_nodes size after reset " + util::Format()( m_Nodes.GetSize())
          );
          return;
        }

        // cutoff is fulfilled, recursively call RemoveNodesWithHighSimilarity on children nodes
        BCL_MessageDbg
        (
          "my girth does fulfill cutoff girth " + util::Format()( m_Girth) + " cutoff " + util::Format()( GIRTH_CUTOFF)
          + " calling Remove...on inner nodes I contain " + util::Format()( m_Nodes.GetSize()) + " nodes"
        );
        // iterate through "m_Nodes" to remove those nodes who are too similar
        for
        (
          typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator itr( m_Nodes.Begin()), itr_end( m_Nodes.End());
          itr != itr_end;
          ++itr
        )
        {
          ( *itr)->RemoveNodesWithHighSimilarity( GIRTH_CUTOFF, BINARY_PREDICATE);
        }
      }

      //! @brief RemoveNodesWithLowSimilarity removes leaf nodes which have a low internal similarity
      //! @param GIRTH_CUTOFF the girth which is used as the cutoff point for removing nodes or not
      //! @param BINARY_PREDICATE the method that should be used to compare girths
      //!        (less = girth is a dissimilarity measure, ie RMSD of protein models)
      //! @note will normally be called only after RemoveNodesBelowSize and/or RemoveNodesWithHighSimilarity, along
      //!       with remove
      void RemoveNodesWithLowSimilarity
      (
        const t_PrecisionType GIRTH_CUTOFF,
        const util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> &BINARY_PREDICATE
      )
      {
        // remove all children nodes, if girth is undefined
        if( !util::IsDefined( m_Girth))
        {
          m_Nodes.Reset();

          // exit function
          return;
        }
        else if( m_Nodes.GetSize() == 1) // if there is only one child node, replace this node with it
        {
          // this temporary variable is necessary otherwise a segfault will result because a ShPtr will go out of scope
          util::ShPtr< Node< t_DataType, t_PrecisionType> > tmp_node( m_Nodes.FirstElement());
          *this = *tmp_node;
          return;
        }

        // iterate through "m_Nodes" to remove those nodes who are too dissimilar
        for
        (
          typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator
            itr( m_Nodes.Begin()), itr_end( m_Nodes.End());
          itr != itr_end;
          // iteration in loop
        )
        {
          // recursive step
          ( *itr)->RemoveNodesWithLowSimilarity( GIRTH_CUTOFF, BINARY_PREDICATE);

          // if the node is a leaf and its members are too dissimilar, remove it
          if( ( *itr)->IsLeaf() && BINARY_PREDICATE( GIRTH_CUTOFF, ( *itr)->GetGirth()))
          {
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::iterator itr_old( itr++);
            m_Nodes.Remove( itr_old);
          }
          else
          {
            // non-leaf node or node with true clusters
            ++itr;
          }
        }

        if( m_Nodes.GetSize() == 1)
        {
          // this temporary variable is necessary otherwise a segfault will result because a ShPtr will go out of scope
          util::ShPtr< Node< t_DataType, t_PrecisionType> > tmp_node( m_Nodes.FirstElement());
          *this = *tmp_node;
        }
      }

      //! @brief CountNumberBaseNodes gives the number of inner nodes which have zero or 1 contained node
      //!        it is important for visualizing a node because it tell how many nodes need to be spaced along x-axis
      //! @param FOR_VISUALIZATION true if the counting is being done for visualization purposes - false by default
      //! @return size_t which gives the number of base nodes contained in this node
      size_t CountNumberBaseNodes( const bool FOR_VISUALIZATION = false) const
      {
        // if this node has no nodes it is a base node so return a count for this node
        if( m_Nodes.IsEmpty())
        {
          return 1;
        }
        else //< this node contains at least one other node so it is not a base node
        {
          size_t count( 0);

          // this will add an extra count which will the result in nodes being staggered (stair stepped) in pymol
          // visualization instead of the nodes being directly overlaid in x-coordinates.
          // Needs to be done here because there could be multiple nodes in this node which only contain one node so
          // if all these nodes are not recursively taken care of (i.e. adding an extra count to each) the bottom of
          // the dendrogram will clash.
          if( m_Nodes.GetSize() == 1 && FOR_VISUALIZATION)
          {
            ++count;
          }

          // iterate through "m_Nodes" to count how many inner nodes have no other nodes
          for
          (
            typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::const_iterator itr( m_Nodes.Begin()), itr_end( m_Nodes.End());
            itr != itr_end;
            ++itr
          )
          {
            count += ( *itr)->CountNumberBaseNodes( FOR_VISUALIZATION);
          }

          return count;
        }
      }

      //! @brief ExpandAllNodes returns a list of all the nodes that make up this node so that nodes buried deep
      //!        in the cluster hierarchy can be directly accessed.
      //! @return SiPtrList to the Nodes that make up this node
      util::SiPtrList< const Node< t_DataType, t_PrecisionType> > ExpandAllNodes() const
      {
        // created a util::SiPtrList< const Node< t_DataType, t_PrecisionType> > which will hold the Nodes
        util::SiPtrList< const Node< t_DataType, t_PrecisionType> > nodes;

        // create a copy of "this" Node
        util::SiPtr< const Node< t_DataType, t_PrecisionType> > this_node( this);

        // add "this_node" to "nodes"
        nodes.InsertElement( this_node);

        // iterate through "m_Nodes" to get the Nodes
         for
         (
           typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::const_iterator itr( m_Nodes.Begin()), itr_end( m_Nodes.End());
           itr != itr_end;
           ++itr
         )
         {
           // have the Node currently denoted by "itr" provide the Nodes that it has and Append these Nodes to "nodes"
           const util::SiPtrList< const Node< t_DataType, t_PrecisionType> > temp( ( *itr)->ExpandAllNodes());
           nodes.Append( temp);
         }

         // return "nodes" which has all of the Nodes whose "m_Girth" is less than "GIRTH"
         return nodes;
      }

      //! @brief get all members that are not contained in child nodes
      util::SiPtrList< const t_DataType> GetOrphanedMembers() const
      {
        // create a set with all members of this node
        storage::Set< util::SiPtr< t_DataType> > all_members( m_Members.Begin(), m_Members.End());

        // iterate through m_Nodes to remove all members held by other nodes from the list
        for
        (
          typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::const_iterator
            itr( m_Nodes.Begin()), itr_end( m_Nodes.End());
          itr != itr_end;
          ++itr
        )
        {
          // get the members of the node
          const util::SiPtrList< t_DataType> &node_members( ( *itr)->GetMembers());

          // iterate through the members of the node
          for
          (
            typename util::SiPtrList< t_DataType>::const_iterator
              itr_member( node_members.Begin()), itr_member_end( node_members.End());
            itr_member != itr_member_end;
            ++itr_member
          )
          {
            all_members.RemoveElement( all_members.Find( *itr_member));
          }
        }
        util::SiPtrList< const t_DataType> members;
        for
        (
          typename storage::Set< util::SiPtr< t_DataType> >::const_iterator
            itr( all_members.Begin()), itr_end( all_members.End());
          itr != itr_end;
          ++itr
        )
        {
          members.PushBack( util::SiPtr< const t_DataType>( **itr));
        }
        return members;
      }

      //! @brief determines the center of the cluster as the member with the smallest sum distance to all other members
      //! @param MEMBER_DISTANCE_FUNCTION
      //! @param BINARY_PREDICATE
      util::SiPtr< const t_DataType> GetCenter
      (
        const util::ShPtr< math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const t_DataType> >, t_PrecisionType> >
          &MEMBER_DISTANCE_FUNCTION,
        const util::ShPtr< util::BinaryFunctionInterface< t_PrecisionType, t_PrecisionType, bool> > &BINARY_PREDICATE
      ) const
      {
        // make sure the node has at least 1 member in it, otherwise doesn't make sense to call GetCenter
        BCL_Assert( !GetMembers().IsEmpty(), "trying to get the center of a node with no members");

        // t_PrecisionType "cluster_center_distance" to keep track of the smallest distance from all other members
        t_PrecisionType cluster_center_distance( util::GetUndefined< t_PrecisionType>());

        // create iterator "cluster_center_itr" to indicate which member is the cluster center
        typename util::SiPtrList< t_DataType>::const_iterator cluster_center_itr( m_Members.Begin());

        // iterate through "m_Members" and determine the cluster center
        for
        (
          typename util::SiPtrList< t_DataType>::const_iterator
            member_itr( m_Members.Begin()), member_itr_end( m_Members.End());
          member_itr != member_itr_end;
          ++member_itr
        )
        {
          // create t_PrecisionType "current_distance" to hold the distance the object denoted by "member_itr" is from
          // all other members of the node currently denoted by "node_itr" and initialize to zero
          t_PrecisionType current_distance( 0.0);

          // iterate through the "members" to add up "current_distance"
          for
          (
            typename util::SiPtrList< t_DataType>::const_iterator member_itr_b( m_Members.Begin());
            member_itr_b != member_itr_end;
            ++member_itr_b
          )
          {
            if( member_itr != member_itr_b)
            {
              BCL_Assert( MEMBER_DISTANCE_FUNCTION.IsDefined(), "Member distance function is not defined");

              // add the distance between the objects currently denoted by "member_itr"
              // and "member_itr_b" to "current_distance"
              current_distance += MEMBER_DISTANCE_FUNCTION->operator()
              (
                storage::VectorND< 2, util::SiPtr< const t_DataType> >( **member_itr, **member_itr_b)
              );
            }
          }

          BCL_MessageDbg( "current distance is " + util::Format()( current_distance));

          if( !util::IsDefined( cluster_center_distance))
          {
            // set "cluster_center_distance" to "current_distance"
            cluster_center_distance = current_distance;

            // set "cluster_center_itr" to "member_itr"
            cluster_center_itr = member_itr;
          }
          // true if the object currently denoted by "member_itr" is the cluster center so far
          else if( BINARY_PREDICATE->operator()( current_distance, cluster_center_distance))
          {
            // set "cluster_center_distance" to "current_distance"
            cluster_center_distance = current_distance;

            // set "cluster_center_itr" to "member_itr"
            cluster_center_itr = member_itr;
          }
          // true if "current_distance" is the same as a "cluster_center_distance"
          else if( current_distance == cluster_center_distance)
          {
            BCL_MessageDbg
            (
              util::Format()( **member_itr) + " is the same distance from all other objects in this node as "
              + util::Format()( **cluster_center_itr) + " with a distance of "
              + util::Format()( cluster_center_distance) + ". " + util::Format()( **cluster_center_itr)
              + " is being kept as the current cluster center."
            );
          }
        }

        // true if the node has more than one member, then the "cluster_center_distance"
        // should be defined
        if( GetMembers().GetSize() > 1)
        {
          // make sure that the "cluster_center_distance" is defined
          BCL_Assert
          (
            util::IsDefined( cluster_center_distance),
            "cluster center could not be determined"
          );
        }

        // return a SiPtr to the node denoted by "cluster_center_itr"
        return util::SiPtr< const t_DataType>( **cluster_center_itr);
      }

      //! @brief determines the largest and smallest defined girth that is observed in the node
      //! @return vectornd 2 with t_PrecisionTypes indicating the smallest and largest defined girth, respectively
      storage::VectorND< 2, t_PrecisionType> GetSmallestLargestDefinedGirth() const
      {
        storage::Set< t_PrecisionType> girths;
        girths.Insert( m_Girth);

        // iterate through "m_Nodes"
        for
        (
          typename util::ShPtrList< Node< t_DataType, t_PrecisionType> >::const_iterator itr( m_Nodes.Begin()), itr_end( m_Nodes.End());
          itr != itr_end;
          ++itr
        )
        {
          const t_PrecisionType current_girth( ( *itr)->GetGirth());
          if( util::IsDefined( current_girth))
          {
            girths.Insert( current_girth);
            const storage::VectorND< 2, t_PrecisionType> current_smallest_largest( ( *itr)->GetSmallestLargestDefinedGirth());
            if( util::IsDefined( current_smallest_largest.First()))
            {
              girths.Insert( current_smallest_largest.First());
            }
            if( util::IsDefined( current_smallest_largest.Second()))
            {
              girths.Insert( current_smallest_largest.Second());
            }
          }
        }

        storage::VectorND< 2, t_PrecisionType> smallest_largest_girth
        (
          util::GetUndefined< t_PrecisionType>(), util::GetUndefined< t_PrecisionType>()
        );

        // true if the set of girths is not empty
        if( !girths.IsEmpty())
        {
          // set the smallest and largest girths
          smallest_largest_girth.First() = *girths.Begin();
          smallest_largest_girth.Second() = *girths.ReverseBegin();
        }

        return smallest_largest_girth;
      }

      //! @brief determines if a node is a leaf i.e. at the very bottom of the dendrogram
      //! @return bool true if the node is a leaf - false otherwise
      bool IsLeaf() const
      {
        return m_Nodes.IsEmpty();
      }

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      virtual std::istream &Read( std::istream &ISTREAM)
      {
        // read in members
        io::Serialize::Read( m_Nodes, ISTREAM);
        // io::Serialize::Read( m_Members, ISTREAM);
        io::Serialize::Read( m_Girth, ISTREAM);
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write "m_Nodes"
        io::Serialize::Write( m_Nodes, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Members, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Girth, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // Node

    // instantiate s_Instance
    template< typename t_DataType, typename t_PrecisionType>
    const util::SiPtr< const util::ObjectInterface> Node< t_DataType, t_PrecisionType>::s_Instance
    (
      GetObjectInstances().AddInstance( new Node< t_DataType, t_PrecisionType>())
    );

    //! @brief operator == checks if two Nodes are the same
    //! @param NODE_A first Node
    //! @param NODE_B second Node
    //! @return true, if Nodes are identical
    template< typename t_DataType, typename t_PrecisionType>
    inline bool operator ==
    (
      const Node< t_DataType, t_PrecisionType> &NODE_A,
      const Node< t_DataType, t_PrecisionType> &NODE_B
    )
    {
      // create bool "girth_matches" to hold whether or not the girths of the two nodes matches
      bool girth_matches;

      // true if both girths are defined
      if( util::IsDefined( NODE_A.GetGirth()) && util::IsDefined( NODE_B.GetGirth()))
      {
        girth_matches = ( NODE_A.GetGirth() == NODE_B.GetGirth());
      }
      else if //< one girth is defined but the other is not
      (
        util::IsDefined( NODE_A.GetGirth()) != util::IsDefined( NODE_B.GetGirth())
      )
      {
        girth_matches = false;
      }
      else //< both girths are not defined
      {
        girth_matches = true;
      }

      // return true if all the members agree with one another; false otherwise
      return NODE_A.GetNodes() == NODE_B.GetNodes()     &&
             NODE_A.GetMembers() == NODE_B.GetMembers() &&
             girth_matches;
    }

  } // namespace cluster
} // namespace bcl

#endif // BCL_CLUSTER_NODE_H_
