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
#include "cluster/bcl_cluster_node.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_wrapper_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    // explicit instantiation of Node< std::string, double>
    template class Node< std::string, double>;
  } // namespace cluster

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_cluster_node.cpp
  //!
  //! @author alexanns
  //! @date Sep 20, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleClusterNode :
    public ExampleInterface
  {
  public:

    ExampleClusterNode *Clone() const
    {
      return new ExampleClusterNode( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "Testing default constructor");

      // create Node "node_a"
      cluster::Node< std::string, double> node_a;

      // check that default constructor worked
      BCL_Example_Check
      (
        node_a.GetNodes().IsEmpty()
        && node_a.GetMembers().IsEmpty()
        && !util::IsDefined( node_a.GetGirth()),
          " m_Nodes should be zero but is " + util::Format()( node_a.GetNodes().GetSize())
        + " m_Members should have size zero but has size " + util::Format()( node_a.GetMembers().GetSize())
        + " m_Girth should be undefined but is " + util::Format()( node_a.GetGirth())
      );

      // test constructor taking a SiPtrList of members
      BCL_MessageStd( "Testing constructor taking a SiPtrList of members");

      // create std::string"second" and initialize with the string "second"
      std::string second( "second");

      // create SiPtrList "members_b" and initialize with "second"
      util::SiPtrList< std::string> members_b( 1, second);

      // call constructor taking a SiPtrList of members
      // create Node "node_b"
      cluster::Node< std::string, double> node_b( members_b);

      // check that constructor taking a SiPtrList of members worked
      BCL_Example_Check
      (
        node_b.GetNodes().IsEmpty()
        && node_b.GetMembers().GetSize() == 1
        && *node_b.GetMembers().FirstElement() == second
        && !util::IsDefined( node_b.GetGirth()),
          " m_Nodes should be zero but is " + util::Format()( node_b.GetNodes().GetSize())
        + " m_Members should have size one but has size " + util::Format()( node_b.GetMembers().GetSize())
        + " m_Member's one element should be " + util::Format()( second) + " but instead is "
        + util::Format()( *node_b.GetMembers().FirstElement())
        + " m_Distance should be undefined but is " + util::Format()( node_b.GetGirth())
      );

      // test constructor taking a ShPtrList of Nodes and a distance
      BCL_MessageStd( "Testing constructor taking a ShPtrList of Nodes and a distance");

      // create std::string"third" and initialize with the string "third"
      std::string third( "third");

      // create SiPtrList "members_c" and initialize with "first"
      util::SiPtrList< std::string> members_c( 1, third);

      // create Node "node_c"
      cluster::Node< std::string, double> node_c( members_c);

      // create ShPtrList of Nodes "nodes_b_c"
      util::ShPtrList< cluster::Node< std::string, double> > nodes_b_c;

      // add new copies of "node_b" and "node_c" to "nodes_b_c"
      nodes_b_c.PushBack( util::ShPtr< cluster::Node< std::string, double> >( new cluster::Node< std::string, double>( node_b)));
      nodes_b_c.PushBack( util::ShPtr< cluster::Node< std::string, double> >( new cluster::Node< std::string, double>( node_c)));

      // create double "dist_between_nodes_b_c" which is the distance between "node_b" and "node_c" and initialize to 1
      double dist_between_nodes_b_c( 1.0);

      // call constructor taking "nodes_b_c" and "dist_between_nodes_b_c"
      // create Node "node_d"
      cluster::Node< std::string, double> node_d( nodes_b_c, dist_between_nodes_b_c);

      // make sure constructor taking a ShPtrList of Nodes and a distance worked
      BCL_Example_Check
      (
        // check size of m_Nodes
           node_d.GetNodes().GetSize() == 2
           // check that member of first node in node_d is "second"
        && *( *node_d.GetNodes().Begin())->GetMembers().FirstElement() == second
           // check that member of second node in node_d is "third"
        && *( *++node_d.GetNodes().Begin())->GetMembers().FirstElement() == third
           // check that node_d has two members
        && node_d.GetMembers().GetSize() == 2
           // check that the first member of node_d is "second"
        && *node_d.GetMembers().FirstElement() == second
           // check that that second member of node_d is "third"
        && **( ++( node_d.GetMembers().Begin())) == third
           // check that the girth of node_d is 1.0
        && math::EqualWithinTolerance( node_d.GetGirth(), 1.0),
          " m_Nodes size should be 2 but is " + util::Format()( node_d.GetNodes().GetSize())
        + " \nthe member of the first node in node_d should be \"" + second + "\" but instead is \""
        + *( *node_d.GetNodes().Begin())->GetMembers().FirstElement()
        + "\" \nthe member of the second node in node_d should be \"" + third + "\" but instead is \""
        + *( *++node_d.GetNodes().Begin())->GetMembers().FirstElement()
        + "\" \nm_Members should have size 2 but has size " + util::Format()( node_b.GetMembers().GetSize())
        + " \nm_Member's first member should be \"" + second + "\" but instead is \"" + *node_d.GetMembers().FirstElement()
        + "\" m_Members second member should be \"third\" but instead is \"" + **( ++( node_d.GetMembers().Begin()))
        + "\" \nm_Girth should be 1.0 but is " + util::Format()( node_d.GetGirth())
      );

      // test copy constructor
      BCL_MessageStd( "Testing copy constructor");

      // call copy constructor
      cluster::Node< std::string, double> node_d_copy( node_d);

      // check the copy constructor worked correctly
      BCL_Example_Check
      (
        // check size of m_Nodes
           node_d_copy.GetNodes().GetSize() == 2
           // check that member of first node in node_d_copy is "second"
        && *( *node_d_copy.GetNodes().Begin())->GetMembers().FirstElement() == second
           // check that member of second node in node_d_copy is "third"
        && *( *++node_d_copy.GetNodes().Begin())->GetMembers().FirstElement() == third
           // check that node_d_copy has two members
        && node_d_copy.GetMembers().GetSize() == 2
           // check that the first member of node_d_copy is "second"
        && *node_d_copy.GetMembers().FirstElement() == second
           // check that that second member of node_d_copy is "third"
        && **( ++( node_d_copy.GetMembers().Begin())) == third
           // check that the girth of node_d_copy is 1.0
        && math::EqualWithinTolerance( node_d_copy.GetGirth(), 1.0),
          " m_Nodes size should be 2 but is " + util::Format()( node_d_copy.GetNodes().GetSize())
        + " \nthe member of the first node in node_d_copy should be \"" + second + "\" but instead is \""
        + *( *node_d_copy.GetNodes().Begin())->GetMembers().FirstElement()
        + "\" \nthe member of the second node in node_d_copy should be \"" + third + "\" but instead is \""
        + *( *++node_d_copy.GetNodes().Begin())->GetMembers().FirstElement()
        + "\" \nm_Members should have size 2 but has size " + util::Format()( node_d_copy.GetMembers().GetSize())
        + " \nm_Member's first member should be \"" + second + "\" but instead is \"" + *node_d_copy.GetMembers().FirstElement()
        + "\" m_Members second member should be \"third\" but instead is \"" + **( ++( node_d_copy.GetMembers().Begin()))
        + "\" \nm_Girth should be 1.0 but is " + util::Format()( node_d_copy.GetGirth())
      );

      // test clone constructor
      util::ShPtr< cluster::Node< std::string, double> > node_clone( node_d.Clone());
      BCL_Example_Check
      (
        *node_clone == node_d,
        "Clone function did not clone node_d : " + util::Format()( node_clone) + " properly to create node_clone " +
        util::Format()( node_d)
      );

    /////////////////
    // data access //
    /////////////////

      // test SetNodes and GetNodes
      BCL_MessageStd( "test SetNodes and GetNodes");

      // create std::string"fourth" and initialize with the string "fourth"
      std::string fourth( "fourth");

      // create SiPtrList "members_d" and initialize with "fourth"
      util::SiPtrList< std::string> members_e( 1, fourth);

      // call constructor taking a SiPtrList of members
      // create Node "node_e"
      cluster::Node< std::string, double> node_e( members_e);

      // create ShPtrList of Nodes "nodes_a_d"
      util::ShPtrList< cluster::Node< std::string, double> > nodes_d_e;

      // add new copies of "node_d" and "node_e" to "nodes_d_e"
      nodes_d_e.PushBack( util::ShPtr< cluster::Node< std::string, double> >( new cluster::Node< std::string, double>( node_d)));
      nodes_d_e.PushBack( util::ShPtr< cluster::Node< std::string, double> >( new cluster::Node< std::string, double>( node_e)));

      // call SetNodes function on node_a; now node_a will have node_d which consists of node_b and node_c, and node_a
      // will also have node_e
      node_a.SetNodes( nodes_d_e);

      // make sure SetNodes function worked by using the GetNodes function
      BCL_Example_Check
      (
        // make sure the size of m_Nodes for node_a is two
           node_a.GetNodes().GetSize() == 2
           // make sure that the first node of node_a ( which is node_d) has two nodes (node_b and node_c)
        && node_a.GetNodes().FirstElement()->GetNodes().GetSize() == 2
           // make sure that the second node of node_a ( which is node_e) has member "fourth"
        && *( node_a.GetNodes().LastElement()->GetMembers().FirstElement()) == fourth,
          " size of m_Nodes for node_a should be 2 but instead is " + util::Format()( node_a.GetNodes().GetSize())
        + "\n first node of node_a ( which is node_d) should have two nodes (node_b and node_c) but instead the size"
        + " of the first node of node_a is " + util::Format()( node_a.GetNodes().FirstElement()->GetNodes().GetSize())
        + "\n the second node of node_a (which is node_e) should have member \"fourth\" "
        + " but instead has \"" + *( node_a.GetNodes().LastElement()->GetMembers().FirstElement()) + "\"\n"
      );

      // test SetMembers and GetMembers
      BCL_MessageStd( "test SetMembers and GetMembers");

      // create SiPtrList "members_a"
      util::SiPtrList< std::string> members_a;

      // append the members of "node_d" to "members_a"
      members_a.Append( node_d.GetMembers());

      // append the members of "node_e" to "members_a"
      members_a.Append( node_e.GetMembers());

      // call SetMembers function on "node_a" with "members_a"
      node_a.SetMembers( members_a);

      // make sure that the members of "node_a" are correct by using GetMembers function
      BCL_Example_Check
      (
        // make sure that "node_a" has three members ( second, third, fourth)
           node_a.GetMembers().GetSize() == 3
           // make sure first member of "node_a" is "second"
        && *( node_a.GetMembers().FirstElement()) == second
           // make sure second member of "node_a" is "third"
        && **( ++node_a.GetMembers().Begin()) == third
           // make sure third member of "node_a" is fourth
        && *( node_a.GetMembers().LastElement()) == fourth,
          " node_a should have three members ( second, third, fourth) but instead has "
        + util::Format()( node_a.GetMembers().GetSize()) + " members \n"
        + " first member of node_a should be \"second\" but instead is "
        + *( node_a.GetMembers().FirstElement()) + "\n"
        + " second member of node_a should be \"third\" but instead is "
        + **( ++node_a.GetMembers().Begin()) + "\n"
        + " third member of node_a should be \"fourth\" but instead is "
        + *( node_a.GetMembers().LastElement()) + "\n"
      );

      // test SetGirth and GetGirth
      BCL_MessageStd( "test SetGirth and GetGirth");

      // call SetGirth function with node_a and set the girth of "node_a" to 2.0
      node_a.SetGirth( 2.0);

      // make sure the girth of "node_a" is 2.0 using the GetGirth function
      BCL_Example_Check
      (
        node_a.GetGirth() == 2.0,
        "node_a girth should be 2.0 but instead is " + util::Format()( node_a.GetGirth()) + "\n"
      );

      BCL_MessageStd( "node_a: \n" + util::Format()( node_a));

    ////////////////
    // operations //
    ////////////////

      // test GetCenter function
      // create SiPtrList "members_get_center" to hold the members for testing the "GetCenter" function
      util::SiPtrList< double> members_get_center;
      double five( 5.0);
      double eight( 8.0);
      double nine( 9.0);
      members_get_center.PushBack( util::SiPtr< double >( five));
      members_get_center.PushBack( util::SiPtr< double >( eight));
      members_get_center.PushBack( util::SiPtr< double >( nine));

      // instantiate the s_Instance needed for this example
      cluster::Node< double, double>::s_Instance.IsDefined();

      cluster::Node< double, double> node_get_center( members_get_center);
      util::SiPtr< const double> center
      (
        node_get_center.GetCenter
        (
          util::ShPtr
          <
            math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const double> >, double>
          >( new DoubleCompare()),
          util::ShPtr< util::BinaryFunctionInterface< double, double, bool> >( ( *math::Comparisons< double>::GetEnums().e_Less)->Clone())
        )
      );
      BCL_Example_Check( *center == eight, "center should be 8.0 but is " + util::Format()( *center));

      // test RemoveNodesBelowSize
      {
        // create SiPtrList "members_get_center" to hold the members for testing the "GetCenter" function
        util::SiPtrList< std::string> members_temp;
        std::string five( "five");
        std::string eight( "eight");
        std::string nine( "nine");
        util::SiPtr< std::string> sp_five( five);
        util::SiPtr< std::string> sp_eight( eight);
        util::SiPtr< std::string> sp_nine( nine);
        members_temp.PushBack( sp_five);
        members_temp.PushBack( sp_eight);
        members_temp.PushBack( sp_nine);
        cluster::Node< std::string, double> node_temp( members_temp);
        nodes_d_e.PushBack( util::ShPtr< cluster::Node< std::string, double> >( node_temp.Clone()));
        cluster::Node< std::string, double> node_remove_below_size( nodes_d_e, 9.4);
        BCL_MessageStd
        (
          "Before RemoveNodesBelowSize " +
          util::Format()( node_remove_below_size.GetNodes().GetSize())
        );
        node_remove_below_size.RemoveNodesBelowSize( 3);
        BCL_MessageStd
        (
          "After RemoveNodesBelowSize " +
          util::Format()( node_remove_below_size.GetNodes().GetSize())
        );
        BCL_ExampleIndirectCheck
        (
          node_remove_below_size.GetNodes().GetSize(),
          1,
          "node_remove_below_size.RemoveNodesBelowSize( 3)"
        );
      }

      // test GetSmallestLargestDefinedGirth
      {
        cluster::Node< std::string, double> node_get_min_max_girth( CreateNode( false));
        BCL_MessageStd
        (
          "GetSmallestLargestDefinedGirth : \n" +
          util::Format()( node_get_min_max_girth.GetSmallestLargestDefinedGirth())
        );
        BCL_Example_Check
        (
          node_get_min_max_girth.GetSmallestLargestDefinedGirth().First() == 3.0 &&
          node_get_min_max_girth.GetSmallestLargestDefinedGirth().Second() == 36.0,
          "min and max should be 3.0 and 36.0 but are " +
          util::Format()( node_get_min_max_girth.GetSmallestLargestDefinedGirth().First()) + " and " +
          util::Format()( node_get_min_max_girth.GetSmallestLargestDefinedGirth().Second())
        );
      }

      // test CountNumberBaseNodes
      {
        BCL_MessageStd( "CountNumberBaseNodes : \n");
        cluster::Node< std::string, double> node_count_base_nodes( CreateNode( false));

        // TODO: Figure out why this next, commented out, line causes the example to crash
        //BCL_MessageDbg( "test node for CountNumberBaseNodes node is " + util::Format()( node_count_base_nodes));
        BCL_MessageDbg( "number members is " + util::Format()( node_count_base_nodes.GetMembers().GetSize()));
        BCL_MessageDbg( "number inner nodes is " + util::Format()( node_count_base_nodes.GetNodes().GetSize()));

        BCL_Example_Check
        (
          node_count_base_nodes.CountNumberBaseNodes() == 8,
          "number of base nodes should be 8 but is " + util::Format()( node_count_base_nodes.CountNumberBaseNodes())
        );
      }

      // test RemoveNodesWithHighSimilarity with Greater and Less
      {
        // create SiPtrList "members_temp" to hold the members for testing the "RemoveNodesWithHighSimilarity" function
        cluster::Node< std::string, double> new_node_lt( CreateNode( false));
        cluster::Node< std::string, double> new_node_gt( CreateNode( true));
        new_node_lt.RemoveNodesWithHighSimilarity( 10.0, **math::Comparisons< double>::GetEnums().e_Less);
        new_node_gt.RemoveNodesWithHighSimilarity( -11.0, **math::Comparisons< double>::GetEnums().e_Greater);
        const size_t end_size_lt( new_node_lt.ExpandAllNodes().GetSize());
        const size_t correct_size_lt( 11);
        const size_t end_size_gt( new_node_gt.ExpandAllNodes().GetSize());
        const size_t correct_size_gt( 9);

        BCL_Example_Check
        (
          end_size_lt == correct_size_lt, "Some nodes were not removed according to girth cutoff and number of nodes left should be "
          + util::Format()( correct_size_lt) + " However number of nodes left is " + util::Format()( end_size_lt)
        )

        BCL_Example_Check
        (
          end_size_gt == correct_size_gt, "Some nodes were not removed according to negated girth cutoff and number of nodes left should be "
          + util::Format()( correct_size_gt) + " However number of nodes left is " + util::Format()( end_size_gt)
        )

//        cluster::Node< std::string, t_PrecisionType> node_temp_b( nodes_d_e, 4.0);
//        nodes_d_e.PushBack( util::ShPtr< cluster::Node< std::string, t_PrecisionType> >( node_temp_b.Clone()));
//        cluster::Node< std::string, t_PrecisionType> node_remove_with_high_similarity( nodes_d_e, 3.0);
//        BCL_Message
//        (
//          util::Message::e_Debug, "Before RemoveNodesWithHighSimilarity with Greater "
//          + util::Format()( node_remove_with_high_similarity)
//          + "\ntotal number of nodes (including top containing node) is "
//          + util::Format()( node_remove_with_high_similarity.ExpandAllNodes().GetSize())
//        );
//        node_remove_with_high_similarity.RemoveNodesWithHighSimilarity( 4.9, *math::Comparisons< double>::GetEnums().e_Greater);
//        BCL_Message
//        (
//          util::Message::e_Debug, "After RemoveNodesWithHighSimilarity with Greater " +
//          util::Format()( node_remove_with_high_similarity)
//          + "\ntotal number of nodes (including top containing node) is "
//          + util::Format()( node_remove_with_high_similarity.ExpandAllNodes().GetSize())
//        );
//        BCL_Example_Check
//        (
//          ExampleClass::ExampleResult::e_Trivial,
//          node_remove_with_high_similarity.ExpandAllNodes().GetSize() == 4,
//          "Should be 4 nodes total in node_remove_with_high_similarity (including itself) but there are " +
//          util::Format()( node_remove_with_high_similarity.ExpandAllNodes().GetSize())
//        );
      }

      BCL_MessageDbg( "node_a : \n" + util::Format()( node_a));

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      BCL_MessageStd( "Testing the write function");
      WriteBCLObject( node_a);

      // can't read simple pointers

    //////////////////////
    // helper functions //
    //////////////////////

      // end
      return 0;
    }

    const cluster::Node< std::string, double> CreateNode( bool NEGATE_GIRTHS) const
    {
      int negation( 1);
      if( NEGATE_GIRTHS)
      {
        negation = -1;
      }

      std::string one( "one");
      std::string two( "two");
      std::string three( "three");
      std::string four( "four");
      std::string five( "five");
      std::string six( "six");
      std::string seven( "seven");
      std::string eight( "eight");
      util::SiPtr< std::string> sp_one( one);
      util::SiPtr< std::string> sp_two( two);
      util::SiPtr< std::string> sp_three( three);
      util::SiPtr< std::string> sp_four( four);
      util::SiPtr< std::string> sp_five( five);
      util::SiPtr< std::string> sp_six( six);
      util::SiPtr< std::string> sp_seven( seven);
      util::SiPtr< std::string> sp_eight( eight);
      util::SiPtrList< std::string> members_one;
      util::SiPtrList< std::string> members_two;
      util::SiPtrList< std::string> members_three;
      util::SiPtrList< std::string> members_four;
      util::SiPtrList< std::string> members_five;
      util::SiPtrList< std::string> members_six;
      util::SiPtrList< std::string> members_seven;
      util::SiPtrList< std::string> members_eight;
      members_one.PushBack( sp_one);
      members_two.PushBack( sp_two);
      members_three.PushBack( sp_three);
      members_four.PushBack( sp_four);
      members_five.PushBack( sp_five);
      members_six.PushBack( sp_six);
      members_seven.PushBack( seven);
      members_eight.PushBack( sp_eight);
      cluster::Node< std::string, double> node_one( members_one);
      cluster::Node< std::string, double> node_two( members_two);
      cluster::Node< std::string, double> node_three( members_three);
      cluster::Node< std::string, double> node_four( members_four);
      cluster::Node< std::string, double> node_five( members_five);
      cluster::Node< std::string, double> node_six( members_six);
      cluster::Node< std::string, double> node_seven( members_seven);
      cluster::Node< std::string, double> node_eight( members_eight);

      util::ShPtrList< cluster::Node< std::string, double> > list_1_2;
      list_1_2.PushBack( util::ShPtr< cluster::Node< std::string, double> >( node_one.Clone()));
      list_1_2.PushBack( util::ShPtr< cluster::Node< std::string, double> >( node_two.Clone()));
      util::ShPtr< cluster::Node< std::string, double> > node_1_2
      ( new cluster::Node< std::string, double>( list_1_2, 3 * negation));

      util::ShPtrList< cluster::Node< std::string, double> > list_3_4;
      list_3_4.PushBack( util::ShPtr< cluster::Node< std::string, double> >( node_three.Clone()));
      list_3_4.PushBack( util::ShPtr< cluster::Node< std::string, double> >( node_four.Clone()));
      util::ShPtr< cluster::Node< std::string, double> > node_3_4
      ( new cluster::Node< std::string, double>( list_3_4, 7 * negation));

      util::ShPtrList< cluster::Node< std::string, double> > list_5_6;
      list_5_6.PushBack( util::ShPtr< cluster::Node< std::string, double> >( node_five.Clone()));
      list_5_6.PushBack( util::ShPtr< cluster::Node< std::string, double> >( node_six.Clone()));
      util::ShPtr< cluster::Node< std::string, double> > node_5_6
      ( new cluster::Node< std::string, double>( list_5_6, 11 * negation));

      util::ShPtrList< cluster::Node< std::string, double> > list_7_8;
      list_7_8.PushBack( util::ShPtr< cluster::Node< std::string, double> >( node_seven.Clone()));
      list_7_8.PushBack( util::ShPtr< cluster::Node< std::string, double> >( node_eight.Clone()));
      util::ShPtr< cluster::Node< std::string, double> > node_7_8
      ( new cluster::Node< std::string, double>( list_7_8, 15 * negation));
      util::ShPtrList< cluster::Node< std::string, double> > list_12_34;
      list_12_34.PushBack( node_1_2);
      list_12_34.PushBack( node_3_4);
      util::ShPtr< cluster::Node< std::string, double> > node_12_34
      ( new cluster::Node< std::string, double>( list_12_34, 10 * negation));

      util::ShPtrList< cluster::Node< std::string, double> > list_56_78;
      list_56_78.PushBack( node_5_6);
      list_56_78.PushBack( node_7_8);
      util::ShPtr< cluster::Node< std::string, double> > node_56_78
      ( new cluster::Node< std::string, double>( list_56_78, 25 * negation));
      util::ShPtrList< cluster::Node< std::string, double> > list_1234_5678;
      list_1234_5678.PushBack( node_12_34);
      list_1234_5678.PushBack( node_56_78);
      const cluster::Node< std::string, double> node_1234_5678( list_1234_5678, 36 * negation);

      return node_1234_5678;
    }

    //! DoubleCompare class is the function interface used for comparing the members of the Nodes in this example
    //! The operator() returns the absolute value of the difference between the two doubles
    class DoubleCompare :
      public math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const double> >, double>
    {
    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! default constructor
      DoubleCompare()
      {
      }

      // clone function
      DoubleCompare *Clone() const
      {
        return new DoubleCompare( *this);
      }

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! virtual operator taking an storage::VectorND< 2, double > and returning a double
      virtual double operator()
      (
        const storage::VectorND< 2, util::SiPtr< const double> > &DOUBLES
      ) const
      {
        return math::Absolute( *DOUBLES.First() - *DOUBLES.Second());
      }

      //! read from std::istream
      virtual std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! write to std::ostream
      virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }
    }; // helper class DoubleCompare

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleClusterNode

  const ExampleClass::EnumType ExampleClusterNode::s_Instance
  (
    GetExamples().AddEnum( ExampleClusterNode())
  );

} // namespace bcl
