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
#include "cluster/bcl_cluster_dendrogram.h"

// includes from bcl - sorted alphabetically
#include "cluster/bcl_cluster_linkage_complete.h"
#include "util/bcl_util_wrapper_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_cluster_dendrogram.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleClusterDendrogram :
    public ExampleInterface
  {
  public:

    ExampleClusterDendrogram *Clone() const
    {
      return new ExampleClusterDendrogram( *this);
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
      // create ShPtrList of double "members" to hold the members to be clustered
      util::ShPtr< storage::List< double> > members
      (
        new storage::List< double>()
      );

      // create some cluster nodes and add to "members"
      members->PushBack( double( 1.0));
      members->PushBack( double( 11.0));
      members->PushBack( double( 6.0));
      members->PushBack( double( 8.0));
      members->PushBack( double( 2.0));
//      members->PushBack( double( 0.0));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "Testing default constructor");
      cluster::Dendrogram< double, double> def_constr;

      // test constructor taking a function interface and a list of members
      BCL_MessageStd( "test constructor taking a function interface");
      cluster::Dendrogram< double, double> param_constr
      (
        util::ShPtr< cluster::LinkageInterface< double, double> >
        (
          new cluster::LinkageComplete< double, double>
          (
            util::ShPtr
            <
              math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const double> >, double>
            >( new DoubleCompare()),
            util::ShPtr< util::BinaryFunctionInterface< double, double, bool> >( *math::Comparisons< double>::GetEnums().e_Less)
          )
        ),
        members
      );

      BCL_MessageStd( "param_const: finished");
      // print the dendrogram "param_constr"
      //BCL_MessageStd( "param_const: \n" + util::Format()( param_constr));

      // make sure "param_constr" has the correct number of nodes in it
      BCL_Example_Check
      (
        param_constr.GetNode().GetNodes().GetSize() == 2,
        "size of param_constr nodes should be 2 but is " + util::Format()( param_constr.GetNode().GetNodes().GetSize())
      );

      // make sure "param_constr" has the correct girth
      BCL_Example_Check
      (
        param_constr.GetNode().GetGirth() == 10,
        "girth of param_constr should be 10 but is " + util::Format()( param_constr.GetNode().GetGirth())
      );

      // make sure "param_constr" has the correct number of members
      BCL_Example_Check
      (
        param_constr.GetNode().GetMembers().GetSize() == 5,
        "the number of members in param_constr should be 5 but is "
        + util::Format()( param_constr.GetNode().GetMembers().GetSize())
      );

      // make sure the two nodes contained in the node of "param_constr" have the correct girth
      BCL_Example_Check
      (
        ( *param_constr.GetNode().GetNodes().Begin())->GetGirth() == 1.0,
        "girth of first node in the Nodes of param_constr should be 1 but is "
        + util::Format()( ( *param_constr.GetNode().GetNodes().Begin())->GetGirth())
      );

      // make sure the two nodes contained in the node of "param_constr" have the correct girth
      BCL_Example_Check
      (
        ( *++param_constr.GetNode().GetNodes().Begin())->GetGirth() == 5.0,
        "girth of first node in the Nodes of param_constr should be 5 but is "
        + util::Format()( ( *param_constr.GetNode().GetNodes().Begin())->GetGirth())
      );

      // test constructor taking a function interface and a list of members and taking a height cutoff
      BCL_MessageStd( "test constructor taking a height cutoff");
      cluster::Dendrogram< double, double> height_constr
      (
        util::ShPtr< cluster::LinkageInterface< double, double> >
        (
          new cluster::LinkageComplete< double, double>
          (
            util::ShPtr
            <
              math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const double> >, double>
            >( new DoubleCompare()),
            util::ShPtr< util::BinaryFunctionInterface< double, double, bool> >( *math::Comparisons< double>::GetEnums().e_Less)
          )
        ),
        members,
        2.0
      );

      // make sure "height_constr" has the correct number of nodes in it
      BCL_Example_Check
      (
        height_constr.GetNode().GetNodes().GetSize() == 3,
        "size of height_constr nodes should be 3 but is "
        + util::Format()( height_constr.GetNode().GetNodes().GetSize())
      );

      // make sure "height_constr" has the correct girth
      BCL_Example_Check
      (
        height_constr.GetNode().GetGirth() == 2.0,
        "girth of height_constr should be 2 but is " + util::Format()( height_constr.GetNode().GetGirth())
      );

      // make sure "height_constr" has the correct number of members
      BCL_Example_Check
      (
        height_constr.GetNode().GetMembers().GetSize() == 5,
        "the number of members in height_constr should be 5 but is "
        + util::Format()( height_constr.GetNode().GetMembers().GetSize())
      );

      // make sure the three nodes contained in the node of "height_constr" have the correct girth
      // first node (11) should have undefined girth
      BCL_Example_Check
      (
        !util::IsDefined( ( *height_constr.GetNode().GetNodes().Begin())->GetGirth()),
        "girth of first node in the Nodes of height_constr should be 1 but is "
        + util::Format()( ( *height_constr.GetNode().GetNodes().Begin())->GetGirth())
      );

      // second node (1 and 2) should have girth of 1
      BCL_Example_Check
      (
        ( *++height_constr.GetNode().GetNodes().Begin())->GetGirth() == 1,
        "girth of second node in the Nodes of height_constr should be 1 but is "
        + util::Format()( ( *++height_constr.GetNode().GetNodes().Begin())->GetGirth())
      );

      // third node (6 and 8) should have girth of 2
      BCL_Example_Check
      (
        ( *++++height_constr.GetNode().GetNodes().Begin())->GetGirth() == 2,
        "girth of third node in the Nodes of height_constr should be 2 but is "
        + util::Format()( ( *++++height_constr.GetNode().GetNodes().Begin())->GetGirth())
      );

      // test preclustering functionality
      {
        BCL_MessageStd( "test with preclustering");
        // create ShPtrList of double "members" to hold the members to be clustered
        util::ShPtr< storage::List< double> > members
        (
          new storage::List< double>()
        );

        // create some cluster nodes and add to "members"
        members->PushBack( double( 1.0));
        members->PushBack( double( 11.0));
        members->PushBack( double( 6.0));
        members->PushBack( double( 8.0));
        members->PushBack( double( 2.0));
        members->PushBack( double( 0.0));

        cluster::Dendrogram< double, double> precluster_constr
        (
          util::ShPtr< cluster::LinkageInterface< double, double> >
          (
            new cluster::LinkageComplete< double, double>
            (
              util::ShPtr
              <
                math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const double> >, double>
              >( new DoubleCompare()),
              util::ShPtr< util::BinaryFunctionInterface< double, double, bool> >( *math::Comparisons< double>::GetEnums().e_Less)
            )
          ),
          members,
          ( **math::Comparisons< double>::GetEnums().e_Less),
          1.1
        );

        BCL_Example_Check
        (
          ( *precluster_constr.GetNode().GetNodes().Begin())->GetNodes().GetSize() == 3,
          "Inner node should have 3 nodes in it but has\n"
          + util::Format()( ( *precluster_constr.GetNode().GetNodes().Begin())->GetNodes().GetSize()) + "\nthe node is\n"
          + util::Format()( *precluster_constr.GetNode().GetNodes().Begin())
        );

        BCL_Example_Check
        (
          ( *precluster_constr.GetNode().GetNodes().Begin())->GetGirth() == 2,
          "Inner node should have girth of 2 but has girth of\n"
          + util::Format()( ( *precluster_constr.GetNode().GetNodes().Begin())->GetGirth()) + "\nthe node is\n"
          + util::Format()( *precluster_constr.GetNode().GetNodes().Begin())
        );
      }

      // test preclustering functionality again
      {
        BCL_MessageStd( "test with preclustering again");
        // create ShPtrList of double "members" to hold the members to be clustered
        util::ShPtr< storage::List< double> > members
        (
          new storage::List< double>()
        );

        // create some cluster nodes and add to "members"
        members->PushBack( double( 2.0));
        members->PushBack( double( 5.0));
        members->PushBack( double( 3.0));
        members->PushBack( double( 4.0));
        members->PushBack( double( 10.0));
        members->PushBack( double( 12.0));

        cluster::Dendrogram< double, double> precluster_constr
        (
          util::ShPtr< cluster::LinkageInterface< double, double> >
          (
            new cluster::LinkageComplete< double, double>
            (
              util::ShPtr
              <
                math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const double> >, double>
              >( new DoubleCompare()),
              util::ShPtr< util::BinaryFunctionInterface< double, double, bool> >( *math::Comparisons< double>::GetEnums().e_Less)
            )
          ),
          members,
          ( **math::Comparisons< double>::GetEnums().e_Less),
          1.1
        );

        BCL_Example_Check
        (
          ( *precluster_constr.GetNode().GetNodes().Begin())->GetNodes().GetSize() == 4,
          "Inner node should have 3 nodes in it but has\n"
          + util::Format()( ( *precluster_constr.GetNode().GetNodes().Begin())->GetNodes().GetSize()) + "\nthe node is\n"
          + util::Format()( *precluster_constr.GetNode().GetNodes().Begin())
        );

        BCL_Example_Check
        (
          ( *precluster_constr.GetNode().GetNodes().Begin())->GetGirth() == 3,
          "Inner node should have girth of 2 but has girth of\n"
          + util::Format()( ( *precluster_constr.GetNode().GetNodes().Begin())->GetGirth()) + "\nthe node is\n"
          + util::Format()( *precluster_constr.GetNode().GetNodes().Begin())
        );
      }

      // test preclustering functionality yet again
      {
        BCL_MessageStd( "test with preclustering yet again");
        // create ShPtrList of double "members" to hold the members to be clustered
        util::ShPtr< storage::List< double> > members
        (
          new storage::List< double>()
        );

        // create some cluster nodes and add to "members"
        members->PushBack( double( 2.0));
        members->PushBack( double( 2.5));
        members->PushBack( double( 3.0));

        cluster::Dendrogram< double, double> precluster_constr
        (
          util::ShPtr< cluster::LinkageInterface< double, double> >
          (
            new cluster::LinkageComplete< double, double>
            (
              util::ShPtr
              <
                math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const double> >, double>
              >( new DoubleCompare()),
              util::ShPtr< util::BinaryFunctionInterface< double, double, bool> >( *math::Comparisons< double>::GetEnums().e_Less)
            )
          ),
          members,
          ( **math::Comparisons< double>::GetEnums().e_Less),
          1.1
        );

        BCL_Example_Check
        (
          precluster_constr.GetNode().GetNodes().GetSize() == 1,
          "Should only be one node in dendrogram node but there are\n"
          + util::Format()( precluster_constr.GetNode().GetNodes().GetSize())
        );

        BCL_Example_Check
        (
          ( *precluster_constr.GetNode().GetNodes().Begin())->GetNodes().GetSize() == 3,
          "Should only 3 nodes in the one node in the dendrogram but there are\n"
          + util::Format()( ( *precluster_constr.GetNode().GetNodes().Begin())->GetNodes().GetSize())
        );
      }

      // test copy constructor
      BCL_MessageStd( "test copy constructor");
      cluster::Dendrogram< double, double> copy_constr( param_constr);
      // make sure "copy_constr" was properly constructed
      BCL_Example_Check
      (
        copy_constr == param_constr,
        "copy constructor did not work properly : original is \n" + util::Format()( param_constr) +
        " \n and copy is \n" + util::Format()( copy_constr)
      );

      // test clone constructor
      BCL_MessageStd( "test clone constructor");
      util::ShPtr< cluster::Dendrogram< double, double> > clone( copy_constr.Clone());
      // make sure "clone" was created properly
      BCL_Example_Check
      (
        *clone == copy_constr,
        "clone constructor did not work properly : original is \n" + util::Format()( copy_constr) +
        " \n and clone is \n" + util::Format()( clone)
      );

    /////////////////
    // data access //
    /////////////////

      // test "GetClassIdentifier" function
      const std::string class_identifier( param_constr.GetClassIdentifier());
      const std::string correct_class_identifier( "bcl::cluster::Dendrogram<double,double>");
      BCL_Example_Check
      (
        class_identifier == correct_class_identifier,
        "GetClassIdentifier function returned wrong name : name should be \n" + correct_class_identifier +
        "\n but name given by GetClassIdentifier is :\n" + class_identifier
      );

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      BCL_MessageStd( "Testing the write function");
      WriteBCLObject( param_constr);

//    // read the file back
//    BCL_MessageStd( "Testing the read function");
//    cluster::Dendrogram< double > read_dendrogram;
//    ReadBCLObject( read_dendogram);

      // since the SiPtrs of the Node can't be read in, "read_dendrogram" will differ from "param_constr"

    //////////////////////
    // helper functions //
    //////////////////////

      // end
      return 0;
    }

  public:

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

  }; //end ExampleClusterDendrogram

  // instantiate s_Instance of DoubleCompare
  const util::SiPtr< const util::ObjectInterface> ExampleClusterDendrogram::DoubleCompare::s_Instance
  (
    GetObjectInstances().AddInstance( new DoubleCompare())
  );

  const ExampleClass::EnumType ExampleClusterDendrogram::s_Instance
  (
    GetExamples().AddEnum( ExampleClusterDendrogram())
  );

} // namespace bcl
