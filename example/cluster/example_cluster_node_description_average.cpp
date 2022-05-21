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
#include "cluster/bcl_cluster_node_description_average.h"

// includes from bcl - sorted alphabetically
#include "cluster/bcl_cluster_node_description_from_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_cluster_node_description_average.cpp
  //!
  //! @author alexanns
  //! @date Sep 20, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleClusterNodeDescriptionAverage :
    public ExampleInterface
  {
  public:

    ExampleClusterNodeDescriptionAverage *Clone() const
    {
      return new ExampleClusterNodeDescriptionAverage( *this);
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
      // create const string "descriptions_filename" and initialize with the file that has the descriptions
      const std::string descriptions_filename( AddExampleInputPathToFilename( e_Cluster, "cluster_pdb_scores.txt"));

      // create method for getting individual member descriptors
      util::ShPtr< util::FunctionInterface< std::string, double> > descriptor_function
      (
        new cluster::NodeDescriptionFromFile( descriptions_filename)
      );

      // create a node to use as an example for calculating its descriptor value
      // create members which will be added to a SiPtrList which will be the members of the node
      std::string pdb_0000( "cluster_0000_final.pdb");
      std::string pdb_0001( "cluster_0001_final.pdb");
      std::string pdb_0002( "cluster_0002_final.pdb");
      std::string pdb_0003( "cluster_0003_final.pdb");
      std::string pdb_0004( "cluster_0004_final.pdb");
      // create SiPtrList "members" which will hold the members of the node
      util::SiPtrList< std::string> members( 1, pdb_0000);
      members.PushBack( pdb_0001);
      members.PushBack( pdb_0002);
      members.PushBack( pdb_0003);
      members.PushBack( pdb_0004);
      // create Node "node" and initialze with "members"
      cluster::Node< std::string, double> node( members);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // call default constructor
      cluster::NodeDescriptionAverage< std::string, double> def_constr;

      // test constructor taking parameter
      cluster::NodeDescriptionAverage< std::string, double> param_constr
      (
        descriptor_function
      );
      {
        const double expected_description( 1.0344708);
        const double given_description( param_constr( node));
        BCL_Example_Check
        (
          expected_description == given_description,
          "Constructor taking a parameter did not work.\n" + util::Format()( param_constr)
        );
      }

      // call copy constructor
      cluster::NodeDescriptionAverage< std::string, double> copy_constr( param_constr);
      {
        const double expected_description( 1.0344708);
        const double given_description( copy_constr( node));
        BCL_Example_Check
        (
          expected_description == given_description,
          "Copy constructor did not work.\n" + util::Format()( copy_constr)
        );
      }

      // virtual copy constructor
      util::ShPtr< cluster::NodeDescriptionAverage< std::string, double> > clone_constr
      (
        copy_constr.Clone()
      );
      {
        const double expected_description( 1.0344708);
        const double given_description( clone_constr->operator()( node));
        BCL_Example_Check
        (
          expected_description == given_description,
          "Clone did not work.\n" + util::Format()( *clone_constr)
        );
      }

    ///////////////////
    //// data access //
    ///////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name
      (
        "bcl::cluster::NodeDescriptionAverage<std::string,double>"
      );
      const std::string given_static_class_name
      (
        GetStaticClassName< cluster::NodeDescriptionAverage< std::string, double> >()
      );
      BCL_Example_Check
      (
        given_static_class_name == correct_static_class_name,
        "GetStaticClassName gives " + given_static_class_name + " but should give " + correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_Example_Check
      (
        correct_static_class_name == clone_constr->GetClassIdentifier(),
        "GetClassIdentifier gives " + clone_constr->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      // test operator()
      // test all member descriptor values defined
      {
        const double expected_description( 1.0344708);
        const double given_description( param_constr( node));
        BCL_Example_Check
        (
          expected_description == given_description,
          "Description value of " + util::Format()( given_description) + " was given instead of expected value\n" +
          util::Format()( expected_description)
        );
      }

      // test not all member descriptor values defined
      {
        // create a node to use as an example for calculating its descriptor value
        // create members which will be added to a SiPtrList which will be the members of the node
        std::string pdb_0000( "not_existant_member");
        std::string pdb_0001( "cluster_0001_final.pdb");
        std::string pdb_0002( "cluster_0002_final.pdb");
        std::string pdb_0003( "non_existent_member_b");
        std::string pdb_0004( "cluster_0004_final.pdb");
        // create SiPtrList "members" which will hold the members of the node
        util::SiPtrList< std::string> members( 1, pdb_0000);
        members.PushBack( pdb_0001);
        members.PushBack( pdb_0002);
        members.PushBack( pdb_0003);
        members.PushBack( pdb_0004);
        // create Node "node" and initialze with "members"
        cluster::Node< std::string, double> node( members);
        const double expected_description( 1.047071);
        const double given_description( param_constr( node));
        BCL_Example_Check
        (
          expected_description == given_description,
          "Description value of " + util::Format()( given_description) + " was given instead of expected value\n" +
          util::Format()( expected_description)
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      BCL_MessageStd( "Testing the write function");
      WriteBCLObject( param_constr);

      // read the file back
      BCL_MessageStd( "Testing the read function");
      cluster::NodeDescriptionAverage< std::string, double> read_output;
      ReadBCLObject( read_output);
      {
        const double expected_description( 1.0344708);
        const double given_description( read_output( node));
        BCL_Example_Check
        (
          expected_description == given_description,
          "Read did not work.\n" + util::Format()( read_output)
        );
      }

    //////////////////////
    // helper functions //
    //////////////////////

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleClusterNodeDescriptionAverage

  const ExampleClass::EnumType ExampleClusterNodeDescriptionAverage::s_Instance
  (
    GetExamples().AddEnum( ExampleClusterNodeDescriptionAverage())
  );

} // namespace bcl
