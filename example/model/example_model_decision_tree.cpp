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
#include "model/bcl_model_decision_tree.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_template_instantiations.h"
#include "model/bcl_model_feature_data_set.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_model_decision_tree.cpp
  //!
  //! @author lemmonwa, teixeipl, mendenjl
  //! @date Aug 26, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleModelDecisionTree :
    public ExampleInterface
  {
  public:

    ExampleModelDecisionTree *Clone() const
    {
      return new ExampleModelDecisionTree( *this);
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
      // create tree for tests
      util::ShPtr< model::DecisionTree> left_node( new model::DecisionTree( linal::Vector< float>( 1, 0.0)));
      util::ShPtr< model::DecisionTree> right_node( new model::DecisionTree( linal::Vector< float>( 1, 100.0)));
      util::ShPtr< model::DecisionTree> tree( new model::DecisionTree( linal::Vector< float>( 1, 50.0)));
      tree->ConfigureBranch( 1, 5.0, left_node, right_node);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      model::DecisionTree decision_tree_default;

      // constructor with all necessary parameters
      model::DecisionTree decision_tree_constructor( linal::Vector< float>( 1, 0.99));

      // check clone
      util::ShPtr< model::DecisionTree> decision_tree_clone( decision_tree_constructor.Clone());

      // check whether decision_node_constructor was initialized with appropriate first pair and correct Activity
      BCL_ExampleCheck( decision_tree_constructor.GetActivity(), float( 0.99));

      // check clone method construction by examining all the pairs within the clone
      util::ShPtr< model::DecisionTree> tree_clone( tree->Clone());

      // check whether decision_node_clone was initialized with appropriate first pair
      BCL_ExampleCheck( tree_clone->GetActivity(), linal::Vector< float>( 1, float( 50.0)));

      // test if the trees are separate entities
      BCL_ExampleIndirectCheck( tree_clone.IsDefined() && tree_clone != tree, true, "clone");
      BCL_ExampleCheck( tree->IsLeaf(), false);

    /////////////////
    // data access //
    /////////////////

      //GetProbabilityOfActivity
      BCL_ExampleCheck( tree_clone->GetActivity(), float( 50.0));

      //GetDecisionNodePairs
      BCL_ExampleCheck( tree_clone->IsLeaf(), false);

    ///////////////
    // operators //
    ///////////////

      //create test data
      linal::Matrix< float> descriptors( 1, 5);
      descriptors( 0, 0) = random::GetGlobalRandom().Random( float( 100));
      descriptors( 0, 1) = random::GetGlobalRandom().Random( float( 100));
      descriptors( 0, 2) = random::GetGlobalRandom().Random( float( 100));
      descriptors( 0, 3) = random::GetGlobalRandom().Random( float( 100));
      descriptors( 0, 4) = float( 5);
      model::FeatureDataSet< float> feature( descriptors);

      //operator()
      BCL_ExampleCheck( model::DecisionTree( linal::Vector< float>( 1, 0.5))( feature)( 0)( 0), linal::Vector< float>( 1, 0.5));

      // test operator()
      // check decision function on test input for left example
      linal::Matrix< float> three_descriptors( 1, 3);
      three_descriptors( 0, 0) = float( 0.1);
      three_descriptors( 0, 1) = float( 0.1);
      three_descriptors( 0, 2) = float( 0.9);
      BCL_ExampleCheck
      (
        tree->operator ()( model::FeatureDataSet< float>( three_descriptors))( 0)( 0),
        float( 0.0)
      );

      // test operator()
      // check decision function on test input for right example
      // in this example, tree is just looking to see whether the first descriptor is <= 10.0, so now try it with
      // an example where that descriptor is > 10.0 to see the other result
      three_descriptors( 0, 1) = 11.0;
      BCL_ExampleCheck
      (
        tree->operator ()( model::FeatureDataSet< float>( three_descriptors))( 0)( 0),
        float( 100.0)
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write bcl object
      WriteBCLObject( *tree);
      // create default object
      model::DecisionTree tree_read;
      // read bcl object
      ReadBCLObject( tree_read);
      // check get oversampling factor method
      BCL_ExampleCheck( tree_read.GetActivity(), float( 50.0));

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    } //end ExampleClass::DtreeNode

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleModelDecisionTree

  const ExampleClass::EnumType ExampleModelDecisionTree::s_Instance
  (
    GetExamples().AddEnum( ExampleModelDecisionTree())
  );

} // namespace bcl
