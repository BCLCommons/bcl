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
#include "graph/bcl_graph_common_subgraph_isomorphism_base.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_graph_common_subgraph_isomorphism_base.cpp
  //!
  //! @author mendenjl
  //! @date Nov 1, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleGraphCommonSubgraphIsomorphismBase :
    public ExampleInterface
  {
  public:

    //! @brief Clone
    ExampleGraphCommonSubgraphIsomorphismBase *Clone() const
    {
      return new ExampleGraphCommonSubgraphIsomorphismBase( *this);
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

    //! @brief Run tests the common subgraph isomorphism base class
    int Run() const
    {
      // create a default constructed csi base
      graph::CommonSubgraphIsomorphismBase default_csi_base;

      // create a csi base for an unconnected csi
      graph::CommonSubgraphIsomorphismBase unconnected_csi_base( graph::CommonSubgraphIsomorphismBase::e_Unconnected);

      // test clone
      BCL_ExampleIndirectCheck
      (
        util::ShPtr< graph::CommonSubgraphIsomorphismBase>( default_csi_base.Clone()).IsDefined(),
        true,
        "CommonSubgraphIsomorphismBase clone"
      );

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck( default_csi_base.GetSolutionType(), graph::CommonSubgraphIsomorphismBase::e_Connected);
      BCL_ExampleCheck( unconnected_csi_base.GetSolutionType(), graph::CommonSubgraphIsomorphismBase::e_Unconnected);

      // check that the initial isomorphism is empty
      BCL_ExampleCheck( default_csi_base.GetIsomorphisms().GetSize(), 0);

      // create a new isomorphism map and check the get/set functions
      storage::Map< size_t, size_t> new_isomorphism;

      // make vertex 0 map to vertex 5
      new_isomorphism[ 0] = 5;

      // set the isomorphism
      default_csi_base.SetIsomorphism( new_isomorphism);
      BCL_ExampleCheck( default_csi_base.GetIsomorphism().GetSize(), 1);
      BCL_ExampleCheck( default_csi_base.GetIsomorphism().GetValue( 0), 5);

    //////////////////////
    // input and output //
    //////////////////////

      BCL_ExampleIndirectCheck
      (
        TestBCLObjectIOForSymmetry( unconnected_csi_base, graph::CommonSubgraphIsomorphismBase()),
        true,
        "CommonSubgraphIsomorphismBase I/O"
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleGraphCommonSubgraphIsomorphismBase

  const ExampleClass::EnumType ExampleGraphCommonSubgraphIsomorphismBase::s_Instance
  (
    GetExamples().AddEnum( ExampleGraphCommonSubgraphIsomorphismBase())
  );
} // namespace bcl
