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
#include "cluster/bcl_cluster_node_description_from_file.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_cluster_node_description_from_file.cpp
  //!
  //! @author alexanns, karakam
  //! @date Sep 20, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleClusterNodeDescriptionFromFile :
    public ExampleInterface
  {
  public:

    ExampleClusterNodeDescriptionFromFile *Clone() const
    {
      return new ExampleClusterNodeDescriptionFromFile( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // call default constructor
      cluster::NodeDescriptionFromFile def_constr;

      // call constructor taking parameters for each member variable
      cluster::NodeDescriptionFromFile param_constr( descriptions_filename);
      {
        const double expected_description( 0.549634);
        const double given_description( param_constr( "cluster_0005_final.pdb"));
        BCL_Example_Check
        (
          expected_description == given_description,
          "Constructor taking a filename parameters did not work.\n" + util::Format()( param_constr)
        );
      }

      // call copy constructor
      cluster::NodeDescriptionFromFile copy_constr( param_constr);
      {
        const double expected_description( -0.103498);
        const double given_description( copy_constr( "cluster_0020_final.pdb"));
        BCL_Example_Check
        (
          expected_description == given_description,
          "Copy constructor did not work.\n" + util::Format()( copy_constr)
        );
      }

      // virtual copy constructor
      util::ShPtr< cluster::NodeDescriptionFromFile> clone_constr( copy_constr.Clone());
      {
        const double expected_description( 0.429461);
        const double given_description
        (
          clone_constr->operator()( "cluster_0000_final.pdb")
        );
        BCL_Example_Check
        (
          expected_description == given_description,
          "Clone constructor did not work.\n" + util::Format()( *clone_constr)
        );
      }

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::cluster::NodeDescriptionFromFile");
      const std::string given_static_class_name
      (
        GetStaticClassName< cluster::NodeDescriptionFromFile>()
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
      // test with a member that does not exist or does not have a description
      {
        const double expected_description( util::GetUndefined< double>());
        const double given_description
        (
          clone_constr->operator()( "non_existant_member")
        );
        BCL_Example_Check
        (
          !util::IsDefined( given_description),
          "Description value of " + util::Format()( given_description) + " was given instead of expected value\n" +
          util::Format()( expected_description)
        );
      }
      // test with normal scenario
      {
        const double expected_description( 0.173507);
        const double given_description
        (
          clone_constr->operator()( "cluster_0019_final.pdb")
        );
        BCL_Example_Check
        (
          expected_description == given_description,
          "Description value of " + util::Format()( given_description) + " was given instead of expected value\n" +
          util::Format()( expected_description)
        );
      }
      // test what if there are too few columns on a line
      {
        const double expected_description( util::GetUndefined< double>());
        const double given_description
        (
          clone_constr->operator()( "too_few_columns")
        );
        BCL_Example_Check
        (
          !util::IsDefined( given_description),
          "Description value of " + util::Format()( given_description) + " was given instead of expected value\n" +
          util::Format()( expected_description)
        );
      }
      // test with too many columns on a line
      {
        const double expected_description( 0.983533);
        const double given_description
        (
          clone_constr->operator()( "cluster_0012_final.pdb")
        );
        BCL_Example_Check
        (
          expected_description == given_description,
          "Description value of " + util::Format()( given_description) + " was given instead of expected value\n" +
          util::Format()( expected_description)
        );
      }
      // test with same string object occuring more than once in file
      {
        const double expected_description( 0.320237);
        const double given_description
        (
          clone_constr->operator()( "cluster_0010_final.pdb")
        );
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
      // initialize ofstream and check that it opened correctly
      WriteBCLObject( param_constr);
      // read the file back
      BCL_MessageStd( "Testing the read function");
      cluster::NodeDescriptionFromFile read_output;
      ReadBCLObject( read_output);
      {
        const double expected_description( 1.60168);
        const double given_description( read_output( "cluster_0003_final.pdb"));
        BCL_Example_Check
        (
          expected_description == given_description,
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

  }; //end NodeDescriptionFromFile

  const ExampleClass::EnumType ExampleClusterNodeDescriptionFromFile::s_Instance
  (
    GetExamples().AddEnum( ExampleClusterNodeDescriptionFromFile())
  );

} // namespace bcl
