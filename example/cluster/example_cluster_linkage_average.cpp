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
#include "cluster/bcl_cluster_linkage_average.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_wrapper_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace cluster
  {
    // explicit instantiation of LinkageAverage< double >
    template class LinkageAverage< double, double>;
  } // namespace cluster

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_cluster_linkage_average.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleClusterLinkageAverage :
    public ExampleInterface
  {
  public:

    ExampleClusterLinkageAverage *Clone() const
    {
      return new ExampleClusterLinkageAverage( *this);
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

      // create some cluster nodes

      // create Node "node_a"
      double one( 1.0);
      cluster::Node< double, double> node_a
      (
        util::SiPtrList< double>( 1, one)
      );
      // create Node "node_b"
      double two( 2.0);
      cluster::Node< double, double> node_b
      (
        util::SiPtrList< double>( 1, two)
      );
      // create Node "node_c"
      double six( 6.0);
      cluster::Node< double, double> node_c
      (
        util::SiPtrList< double>( 1, six)
      );
      // create Node "node_d"
      double eight( 8.0);
      cluster::Node< double, double> node_d
      (
        util::SiPtrList< double>( 1, eight)
      );
      // create Node "node_e"
      double eleven( 11.0);
      cluster::Node< double, double> node_e
      (
        util::SiPtrList< double>( 1, eleven)
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "Testing default constructor");
      cluster::LinkageAverage< double, double> linkage_def;

      // test constructor taking a function interface
      BCL_MessageStd( "test constructor taking a function interface");
      cluster::LinkageAverage< double, double> linkage_func
      (
        util::ShPtr< math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const double> >, double> >( new DoubleCompare())
      );

      // test copy constructor
      BCL_MessageStd( "test copy constructor");
      cluster::LinkageAverage< double, double> linkage_copy( linkage_func);

      // test clone constructor
      BCL_MessageStd( "test clone constructor");
      util::ShPtr< cluster::LinkageAverage< double, double> > linkage_clone( linkage_copy.Clone());

      // make sure the function interface, copy, and clone constructors worked
      BCL_Example_Check
      (
        linkage_func( node_e, node_d) == 3.0
        && linkage_copy( node_a, node_b) == 1.0
        && linkage_clone->operator()( node_c, node_b) == 4.0,
          "using the operator() to check the construction worked correctly, linkage_func should return 3.0 and it returned"
        + util::Format()( linkage_func( node_e, node_d))
        + "\nlinkage_copy should return 1.0 and it returned " + util::Format()( linkage_copy( node_a, node_b))
        + "\nlinkage_clone should return 4.0 and it returned "
        + util::Format()( linkage_clone->operator()( node_c, node_b))
      );

    /////////////////
    // data access //
    /////////////////

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      BCL_MessageStd( "Testing operator()");

      // create ShPtrList of Nodes "nodes_b_c"
      util::ShPtrList< cluster::Node< double, double> > nodes_a_b_c;
      // add node_a to nodes_a_b_c
      nodes_a_b_c.PushBack( util::ShPtr< cluster::Node< double, double> >( new cluster::Node< double, double>( node_a)));
      // add node_b to nodes_a_b_c
      nodes_a_b_c.PushBack( util::ShPtr< cluster::Node< double, double> >( new cluster::Node< double, double>( node_b)));
      // add node_c to nodes_a_b_c
      nodes_a_b_c.PushBack( util::ShPtr< cluster::Node< double, double> >( new cluster::Node< double, double>( node_c)));

      // create ShPtrList of Nodes "nodes_d_e"
      util::ShPtrList< cluster::Node< double, double> > nodes_d_e;
      // add node_d to nodes_d_e
      nodes_d_e.PushBack( util::ShPtr< cluster::Node< double, double> >( new cluster::Node< double, double>( node_d)));
      // add node_e to nodes_d_e
      nodes_d_e.PushBack( util::ShPtr< cluster::Node< double, double> >( new cluster::Node< double, double>( node_e)));

      // create Node "node_a_b_c"
      cluster::Node< double, double> node_a_b_c
      (
        nodes_a_b_c, 5.0
      );

      // create Node "node_d_e"
      cluster::Node< double, double> node_d_e
      (
        nodes_d_e, 3.0
      );

      // test operator(): create double "distance" and initialize with the distance calculated by "linkage_copy"
      double distance( linkage_copy( node_a_b_c, node_d_e));

      // make sure that the correct distance was calculated
      BCL_Example_Check
      (
        math::EqualWithinTolerance( distance, 6.5),
        "distance between node_a_b_c and node_d_e should be 6.5 and is calculated as "
        + util::Format()( linkage_copy( node_a_b_c, node_d_e))
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      BCL_MessageStd( "Testing the write function");
      WriteBCLObject( linkage_func);

//    // read the file back
//    BCL_MessageStd( "Testing the read function");
//    cluster::LinkageComplete< double > linkage_read;
//    ReadBCLObject( linkage_read);

//    // check both collectors are same
//    BCL_MessageStd( "Checking the write and read worked");
//    BCL_Example_Check
//    (
//      ExampleClass::ExampleResult::e_Trivial,
//      math::EqualWithinTolerance( linkage_read( node_d_e, node_b), 9.0),
//        "read in object \"linkage_read\" should have calculated the distance between \"node_d_e\" and \"node_b\" as 9.0 and calcualted the distance as "
//      + util::Format()( linkage_read( node_d_e, node_b))
//    );

    //////////////////////
    // helper functions //
    //////////////////////

      // end
      return 0;
    }

  private:

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
      virtual double operator()( const storage::VectorND< 2, util::SiPtr< const double> > &DOUBLES) const
      {
        return math::Absolute( ( *DOUBLES.First()) - ( *DOUBLES.Second()));
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
    }; // helper class ExampleLinkageComplete

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleClusterLinkageAverage

  const ExampleClass::EnumType ExampleClusterLinkageAverage::s_Instance
  (
    GetExamples().AddEnum( ExampleClusterLinkageAverage())
  );
} // namespace bcl
