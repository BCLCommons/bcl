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
#include "fold/bcl_fold_mutate_tree.h"

// includes from bcl - sorted alphabetically
#include "fold/bcl_fold_default_mutates.h"

// external includes - sorted alphabetically

namespace bcl
{

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_fold_mutate_tree.cpp
  //! @brief this example tests the implementation of fold::MutateTree.
  //!
  //! @author fischea
  //! @date Aug 19, 2016
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  class ExampleFoldMutateTree :
    public ExampleInterface
  {

  //////////
  // data //
  //////////

  public:

    //! single instance of this class
    static const ExampleClass::EnumType s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief clone function
    //! @return pointer to a new ExampleFoldMutateTree
    ExampleFoldMutateTree *Clone() const
    {
      return new ExampleFoldMutateTree( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the name of this class
    //! @return the name of this class
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! @detail this is performing the execution of the example
    int Run() const
    {

    //////////////////////
    // data preparation //
    //////////////////////

      // initialize mutates
      fold::DefaultMutates &mutates( fold::DefaultMutates::GetInstance());
      mutates.InitializeMutates();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create an empty mutate tree
      fold::MutateTree tree;
      mutates.InitializeMutateTree( tree);

      // check if mutate probabilities can be set
      tree.SetMutateProbability( mutates.e_AddSSENextToSSE, 0.5);
      BCL_ExampleCheck( tree.GetMutateProbability( mutates.e_AddSSENextToSSE), 0.5);
      tree.SetMutateProbability( mutates.e_SSERotateLarge, 0.9);
      BCL_ExampleCheck( tree.GetMutateProbability( mutates.e_SSERotateLarge), 0.9);

      // store the string representation of the mutate tree
      const std::string parameters( tree.GetLabel().ToString());
      fold::MutateTree param_tree;
      param_tree.AssertRead( parameters);

      // check if both mutate trees are identical
      BCL_ExampleCheck( tree.GetMutateProbabilities().GetSize(), param_tree.GetMutateProbabilities().GetSize());
      BCL_ExampleCheck
      (
        tree.GetMutateTypeProbabilities().GetSize(), param_tree.GetMutateTypeProbabilities().GetSize()
      );

    /////////////////
    // data access //
    /////////////////

      // check class name
      BCL_ExampleCheck( tree.GetClassIdentifier(), ( GetStaticClassName< fold::MutateTree>()));

    ///////////////
    // operators //
    ///////////////

      return 0;
    }

  }; // class ExampleFoldMutateTree

  //! single instance of this class
  const ExampleClass::EnumType ExampleFoldMutateTree::s_Instance
  (
     GetExamples().AddEnum( ExampleFoldMutateTree())
  );

} // namespace bcl
