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
#include "math/bcl_math_mutate_decision_node.h"

// includes from bcl - sorted alphabetically
#include "mc/bcl_mc_mutate_loop_add.h"
#include "mc/bcl_mc_mutate_loop_remove.h"

// external includes - sorted alphabetically

#undef AddAtom

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_mutate_decision_node.cpp
  //!
  //! @author karakam, fischea
  //! @date Oct 23, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathMutateDecisionNode :
    public ExampleInterface
  {
  public:

    ExampleMathMutateDecisionNode *Clone() const
    { return new ExampleMathMutateDecisionNode( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief run routine
    //! @detail this is performing the execution of the example
    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create the decision node
      math::MutateDecisionNode< assemble::ProteinModel> node;

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( node.GetClassIdentifier(), ( GetStaticClassName< math::MutateDecisionNode< assemble::ProteinModel> >()));

    ///////////////
    // operators //
    ///////////////

      // create two mutates for testing
      const mc::MutateLoopAdd mutate_add;
      const mc::MutateLoopRemove mutate_remove;

      // add both mutates to the decision node with their corresponding probabilities
      node.AddMutate( mutate_add, 0.5);
      node.AddMutate( mutate_remove, 0.2);

    //////////////////////
    // input and output //
    //////////////////////

      // construct using serializer
      const std::string serialization_string( node.GetString());
      math::MutateDecisionNode< assemble::ProteinModel> node_serialized;
      node_serialized.AssertRead( serialization_string);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathMutateDecisionNode

  const ExampleClass::EnumType ExampleMathMutateDecisionNode::s_Instance
  (
    GetExamples().AddEnum( ExampleMathMutateDecisionNode())
  );

} // namespace bcl
