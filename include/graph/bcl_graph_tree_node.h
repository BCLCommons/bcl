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

#ifndef BCL_GRAPH_TREE_NODE_H_
#define BCL_GRAPH_TREE_NODE_H_

// include the namespace header
#include "bcl_graph.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_graph_undirected_edge.h"
#include "linal/bcl_linal_matrix.h"
#include "math/bcl_math.h"
#include "storage/bcl_storage_list.h"
#include "storage/bcl_storage_set.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_own_ptr.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically
#include <numeric>

namespace bcl
{
  namespace graph
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class TreeNode
    //! @brief a class which represents a node of a tree: single parent, multiple children.  Acts as an auxillary class
    //!        for a higher tree class
    //!
    //! @see @link example_graph_tree_node.cpp @endlink
    //! @author geanesar
    //! @date Mar 22, 2016
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_VertexData, typename t_EdgeData>
    class TreeNode :
      public util::ObjectInterface
    {
    private:

      t_VertexData m_Data; //!< data stored in the node
      util::SiPtr< TreeNode< t_VertexData, t_EdgeData> > m_Parent; //!< parent node
      storage::Map< util::SiPtr< TreeNode< t_VertexData, t_EdgeData> >, t_EdgeData> m_Children; //!< children nodes

    //////////
    // data //
    //////////

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor, creates an empty graph
      TreeNode() :
        m_Data(),
        m_Parent(),
        m_Children()
      {
      }

      //! @brief constructor with data
      TreeNode( const t_VertexData &DATA) :
        m_Data( DATA),
        m_Parent(),
        m_Children()
      {
      }

      //! clone the object
      TreeNode *Clone() const
      {
        return new TreeNode( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @return class name of the object behind a pointer or the current object
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief return whether this tree node is a leaf, i.e. if it has no children
      bool IsLeaf() const
      {
        return m_Children.IsEmpty();
      }

      //! @brief return whether the tree node is a root, i.e. if it has no parents
      bool IsRoot() const
      {
        return !m_Parent.IsDefined();
      }

      //! @brief get whether the node is isolated, i.e. no parents or children
      bool IsIsolated() const
      {
        return IsRoot() && IsLeaf();
      }

      //! @brief get the contained data
      const t_VertexData &GetData() const
      {
        return m_Data;
      }

      //! @brief get the contained data
      t_VertexData &GetData()
      {
        return m_Data;
      }

      //! @brief set the data contained in the node
      void SetData( const t_VertexData &DATA)
      {
        m_Data = DATA;
      }

      //! @brief get the maximum depth of this node with regards to any of its parents
      size_t GetDepth() const
      {
        return IsRoot() ? 0 : m_Parent->GetDepth() + 1; // add 1 for the parent node
      }

      //! @brief get the distance from this node to the most distant leaf node under it
      size_t GetHeight() const
      {
        size_t height( 0);
        for
        (
          typename storage::Map< util::SiPtr< TreeNode< t_VertexData, t_EdgeData> >, t_EdgeData>::const_iterator itr_child( m_Children.Begin()), 
            itr_child_end( m_Children.End());
          itr_child != itr_child_end;
          ++itr_child
        )
        {
          if( itr_child->first.IsDefined())
          {
            size_t newheight( itr_child->first->GetHeight() + 1); // add one for moving to the child node
            height = std::max( newheight, height);
          }
        }
        return height;
      }

      //! @brief Get the parents of the node
      const util::SiPtr< TreeNode< t_VertexData, t_EdgeData> > &GetParent() const
      {
        return m_Parent;
      }

      //! @brief query if a parent is present
      bool HasParent( const util::SiPtr< TreeNode< t_VertexData, t_EdgeData> > &NODE) const
      {
        return m_Parent == NODE;
      }

      //! @brief add a parent to this node
      void SetParent( const util::SiPtr< TreeNode< t_VertexData, t_EdgeData> > &NODE)
      {
        if( NODE.IsDefined())
        {
          m_Parent = NODE;
        }
      }

      //! @brief remove a parent from this node
      void RemoveParent()
      {
        SetParent( util::SiPtr< TreeNode< t_VertexData, t_EdgeData> >());
      }

      //! @brief get the children of the node
      const storage::Map< util::SiPtr< TreeNode< t_VertexData, t_EdgeData> >, t_EdgeData> &GetInternalChildrenInfo() const
      {
        return m_Children;
      }

      //! @brief query if a child is present
      bool HasChild( const util::SiPtr< TreeNode< t_VertexData, t_EdgeData> > &NODE) const
      {
        return m_Children.Has( NODE);
      }

      //! @brief add a child node
      void AddChild( util::SiPtr< TreeNode< t_VertexData, t_EdgeData> > NODE, const t_EdgeData &EDGE_DATA)
      {
        if( NODE.IsDefined())
        {
          m_Children[ NODE] = EDGE_DATA;
          util::SiPtr< TreeNode< t_VertexData, t_EdgeData> > this_sp( this);
          NODE->SetParent( this_sp);
        }
      }

      //! @brief remove a child node
      void RemoveChild( util::SiPtr< TreeNode< t_VertexData, t_EdgeData> > NODE)
      {
        if( NODE.IsDefined())
        {
          m_Children.RemoveElement( NODE);
          util::SiPtr< TreeNode< t_VertexData, t_EdgeData> > this_sp( this);
          NODE->RemoveParent( this_sp);
        }
      }

      //! @brief get the size of the tree starting from this node
      size_t GetSize() const
      {
        size_t size( 1);
        for
        (
          typename storage::Map< util::SiPtr< TreeNode< t_VertexData, t_EdgeData> >, t_EdgeData>::const_iterator itr_child( m_Children.Begin()),
            itr_child_end( m_Children.End());
          itr_child != itr_child_end;
          ++itr_child
        )
        {
          if( itr_child->first.IsDefined())
          {
            size += itr_child->first->GetSize();
          }
        }
        return size;
      }

      //! @brief get a list of all leaf nodes in the tree
      util::SiPtrVector< TreeNode< t_VertexData, t_EdgeData> > GetLeaves() const
      {
        util::SiPtrVector< TreeNode< t_VertexData, t_EdgeData> > leaves;
        if( IsLeaf())
        {
          leaves.PushBack( util::SiPtr< TreeNode< t_VertexData, t_EdgeData> >( this));
        }
        else
        {
          for
          (
            typename storage::Map< util::SiPtr< TreeNode< t_VertexData, t_EdgeData> >, t_EdgeData>::const_iterator itr_child( m_Children.Begin()),
              itr_child_end( m_Children.End());
            itr_child != itr_child_end;
            ++itr_child
          )
          {
            if( itr_child->first.IsDefined())
            {
              leaves.Append( itr_child->first->GetLeaves());
            }
          }
        }
        return leaves;
      }

    ////////////////
    // operations //
    ////////////////

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

    //////////////////////
    // helper functions //
    //////////////////////

    public:

    }; // template class Graph

  } // namespace graph
} // namespace bcl

#endif // BCL_GRAPH_TREE_NODE_H_

