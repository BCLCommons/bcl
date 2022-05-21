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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "graph/bcl_graph_tree_node.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_graph_tree_node.cpp
  //!
  //! @author geanesar
  //! @date Mar 22, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleGraphTreeNode :
    public ExampleInterface
  {
  public:

    ExampleGraphTreeNode *Clone() const
    {
      return new ExampleGraphTreeNode( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      storage::Vector< graph::TreeNode< int, int> > nodes( 4);
      nodes( 0).SetData( 0);
      nodes( 1).SetData( 1);
      nodes( 2).SetData( 2);
      nodes( 3).SetData( 3);

      util::SiPtrVector< graph::TreeNode< int, int> > node_ptrs( 4);
      node_ptrs( 0) = &nodes( 0);
      node_ptrs( 1) = &nodes( 1);
      node_ptrs( 2) = &nodes( 2);
      node_ptrs( 3) = &nodes( 3);

      nodes( 0).AddChild( node_ptrs( 1), 1);
      nodes( 1).AddChild( node_ptrs( 2), 1);

      BCL_ExampleCheck( nodes( 0).IsRoot(), true);
      BCL_ExampleCheck( nodes( 1).IsRoot(), false);

      BCL_ExampleCheck( nodes( 2).IsLeaf(), true);
      BCL_ExampleCheck( nodes( 1).IsLeaf(), false);

      BCL_ExampleCheck( nodes( 3).IsIsolated(), true);
      BCL_ExampleCheck( nodes( 0).IsIsolated(), false);

      BCL_ExampleCheck( nodes( 0).GetDepth(), 0);
      BCL_ExampleCheck( nodes( 1).GetDepth(), 1);
      BCL_ExampleCheck( nodes( 2).GetDepth(), 2);

      BCL_ExampleCheck( nodes( 0).GetHeight(), 2);
      BCL_ExampleCheck( nodes( 1).GetHeight(), 1);
      BCL_ExampleCheck( nodes( 2).GetHeight(), 0);

      return 0;
    }
     // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleGraphTreeNode

  const ExampleClass::EnumType ExampleGraphTreeNode::s_Instance
  (
    GetExamples().AddEnum( ExampleGraphTreeNode())
  );

} // namespace bcl
