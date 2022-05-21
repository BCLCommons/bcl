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
#include "cluster/bcl_cluster_node_colorer.h"

// includes from bcl - sorted alphabetically
#include "cluster/bcl_cluster_node_description_average.h"
#include "cluster/bcl_cluster_node_description_from_file.h"
#include "linal/bcl_linal_vector_operations.h"
#include "util/bcl_util_color_gradient.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    // explicit instantiation of cluster::NodeDescriptionAverage< std::string>
    template class NodeDescriptionAverage< std::string, double>;
  } // namespace cluster

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_cluster_node_colorer.cpp
  //!
  //! @author alexanns
  //! @date Sep 17, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleClusterNodeColorer :
    public ExampleInterface
  {
  public:

    ExampleClusterNodeColorer *Clone() const
    {
      return new ExampleClusterNodeColorer( *this);
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
      // just have a constant color
      const storage::Vector< util::Color> constant_gradient_points( 2, util::GetColors().e_Red);

      // create coloring function "color_function" initialize and initialize with a OutputPymolColorGradient
      // that will give a constant color
      const util::ShPtr< util::FunctionInterface< double, linal::Vector3D> > color_function
      (
        new util::ColorGradient( math::Range< double>( 0.0, 1.0), constant_gradient_points)
      );

      // create const string "descriptions_filename" and initialize with the file that has member descriptions
      const std::string descriptions_filename( AddExampleInputPathToFilename( e_Cluster, "cluster_pdb_scores.txt"));

      // create method for getting individual member descriptors
      util::ShPtr< util::FunctionInterface< std::string, double> > member_descriptor_function
      (
        new cluster::NodeDescriptionFromFile( descriptions_filename)
      );

      // create descriptor function "node_descriptor_function"
      util::ShPtr< util::FunctionInterface< cluster::Node< std::string, double>, double> >
      node_descriptor_function
      (
        new cluster::NodeDescriptionAverage< std::string, double>( member_descriptor_function)
      );

      // create a node to use as an example for calculating its color value
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
      cluster::NodeColorer< std::string, double> def_constr;

      // call constructor taking parameters for each member variable
      cluster::NodeColorer< std::string, double> param_constr
      (
        node_descriptor_function, color_function
      );
      {
        const linal::Vector3D &expected_color( *constant_gradient_points.FirstElement());
        const linal::Vector3D given_color( param_constr( node));
        BCL_Example_Check
        (
          expected_color == given_color,
          "Constructor taking parameters did not work.\n" + util::Format()( param_constr)
        );
      }

      // call copy constructor
      cluster::NodeColorer< std::string, double> copy_constr( param_constr);
      {
        const linal::Vector3D &expected_color( *constant_gradient_points.FirstElement());
        const linal::Vector3D given_color( copy_constr( node));
        BCL_Example_Check
        (
          expected_color == given_color,
          "Copy constructor did not work.\n" + util::Format()( copy_constr)
        );
      }

      // virtual copy constructor
      util::ShPtr< cluster::NodeColorer< std::string, double> > clone_constr( copy_constr.Clone());
      {
        const linal::Vector3D &expected_color( *constant_gradient_points.FirstElement());
        const linal::Vector3D given_color( clone_constr->operator()( node));
        BCL_ExampleIndirectCheck( expected_color, given_color, "Clone constructor");
      }

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::cluster::NodeColorer<std::string,double>");
      const std::string given_static_class_name
      (
        GetStaticClassName< cluster::NodeColorer< std::string, double> >()
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
      {
        const linal::Vector3D expected_color( *constant_gradient_points.FirstElement());
        const linal::Vector3D given_color( copy_constr( node));
        BCL_Example_Check
        (
          expected_color == given_color,
          "Operator() did not work. Expected result :\n" + util::Format()( expected_color) + "Given result :\n" +
          util::Format()( given_color)
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      BCL_MessageStd( "Testing the write function");
      // initialize ofstream and check that it opened correctly
      WriteBCLObject( param_constr);

      // read the file back
      BCL_MessageStd( "Testing the read function");
      cluster::NodeColorer< std::string, double> read_output;
      ReadBCLObject( read_output);

      {
        const linal::Vector3D expected_color( *constant_gradient_points.FirstElement());
        const linal::Vector3D given_color( read_output( node));
        BCL_Example_Check
        (
          expected_color == given_color,
          "Reading in did not work.\n" + util::Format()( read_output)
        );
      }

    //////////////////////
    // helper functions //
    //////////////////////

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end NodeColorer

  const ExampleClass::EnumType ExampleClusterNodeColorer::s_Instance
  (
    GetExamples().AddEnum( ExampleClusterNodeColorer())
  );

} // namespace bcl
