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
#include "cluster/bcl_cluster_output_centers.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_wrapper_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_cluster_output_centers.cpp
  //! @details example
  //!
  //! @author alexanns, karakam
  //! @date Sep 20, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleClusterOutputCenters :
    public ExampleInterface
  {
  public:

    ExampleClusterOutputCenters *Clone() const
    {
      return new ExampleClusterOutputCenters( *this);
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
      // create some nodes which can then be used to test the functions of a ClusterOutputCenters object

      util::SiPtrList< double> members_five;
      double five( 5.0);
      members_five.PushBack( util::SiPtr< double>( five));
      cluster::Node< double, double> node_five( members_five);

      util::SiPtrList< double> members_six;
      double six( 6.0);
      members_six.PushBack( util::SiPtr< double>( six));
      cluster::Node< double, double> node_six( members_six);

      util::SiPtrList< double> members_nine;
      double nine( 9.0);
      members_nine.PushBack( util::SiPtr< double>( nine));
      cluster::Node< double, double> node_nine( members_nine);

      util::SiPtrList< double> members_eight;
      double eight( 8.0);
      members_eight.PushBack( util::SiPtr< double>( eight));
      cluster::Node< double, double> node_eight( members_eight);

      util::ShPtrList< cluster::Node< double, double> > nodes_five_eight;
      nodes_five_eight.PushBack( util::ShPtr< cluster::Node< double, double> >( node_five.Clone()));
      nodes_five_eight.PushBack( util::ShPtr< cluster::Node< double, double> >( node_eight.Clone()));
      cluster::Node< double, double> node_five_eight( nodes_five_eight, 3.0);

      util::ShPtrList< cluster::Node< double, double> > nodes_five_six_eight
      (
        1, util::ShPtr< cluster::Node< double, double> >( node_five_eight.Clone())
      );
      nodes_five_six_eight.PushBack( util::ShPtr< cluster::Node< double, double> >( node_six.Clone()));
      cluster::Node< double, double> node_five_six_eight( nodes_five_six_eight, 3.1);

      util::ShPtrList< cluster::Node< double, double> > nodes_five_six_eight_nine
      (
        1, util::ShPtr< cluster::Node< double, double> >( node_five_six_eight.Clone())
      );
      nodes_five_six_eight_nine.PushBack
      (
        util::ShPtr< cluster::Node< double, double> >( node_nine.Clone())
      );

      // create Node "node_five_six_eight_nine" and initialize to contain all the nodes from above
      cluster::Node< double, double> node_five_six_eight_nine( nodes_five_six_eight_nine, 3.5);

      BCL_MessageDbg
      (
        "Node is " + util::Format()( node_five_six_eight_nine)
      );

      // create SiPtrList of nodes "nodes_to_output" and add "node_five_six_eight_nine" to it twice so that
      // outputting more than one node can be tested
      util::SiPtrList< const cluster::Node< double, double> > nodes_to_output;
      nodes_to_output.PushBack( util::SiPtr< const cluster::Node< double, double> >( node_five_six_eight_nine));
      nodes_to_output.PushBack( util::SiPtr< const cluster::Node< double, double> >( node_five_six_eight_nine));

      // create ShPtr to DoubleCompare to compare the members of a node
      util::ShPtr
      <
        math::FunctionInterfaceSerializable< storage::VectorND< 2, util::SiPtr< const double> >, double>
      > compare( new DoubleCompare());

      // create const std::string "output_filename"
      // initialize with the name of the file the cluster centers will be written to
      const std::string output_filename
      (
        AddExampleOutputPathToFilename( cluster::GetNamespaceIdentifier(), "bcl_cluster_output_centers_write_output.txt")
      );

      // create const std::string "output_filename_correct"
      // initialize with the name of the file that has the correct output
      const std::string output_filename_correct
      (
        AddExampleOutputPathToFilename( cluster::GetNamespaceIdentifier(), "bcl_cluster_output_centers_write_output_correct.txt")
      );

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // Default constructor
      cluster::OutputCenters< double, double> default_constr;

      // constructor taking all member variables as parameters
      cluster::OutputCenters< double, double> param_constr
      (
        compare,
        util::ShPtr< util::BinaryFunctionInterface< double, double, bool> >
        (
            ( *math::Comparisons< double>::GetEnums().e_Less)->Clone()
        )
      );
      // write the output
      param_constr.WriteOutput( output_filename, nodes_to_output);

      // check that the output is correct
      BCL_Example_Check
      (
        io::File::FilesMatch( output_filename_correct, output_filename),
        "files don't match"
      );

      // copy constructor
      cluster::OutputCenters< double, double> copy_constr( param_constr);
      // write the output
      copy_constr.WriteOutput( output_filename, nodes_to_output);
      // check that the output is correct
      BCL_Example_Check
      (
        io::File::FilesMatch( output_filename_correct, output_filename),
        "files don't match"
      );

      // clone constructor
      util::ShPtr< cluster::OutputCenters< double, double> > clone_constr( param_constr.Clone());
      // write the output
      clone_constr->WriteOutput( output_filename, nodes_to_output);
      // check that the output is correct
      BCL_Example_Check
      (
        io::File::FilesMatch( output_filename_correct, output_filename),
        "files don't match"
      );

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::cluster::OutputCenters<double,double>");
      const std::string classname( GetStaticClassName< cluster::OutputCenters< double, double> >());
      BCL_Example_Check
      (
        classname == correct_static_class_name,
        "GetStaticClassName gives " + classname + " but should give " +
        correct_static_class_name
      );

      // check GetClassIdentifier
      BCL_Example_Check
      (
        correct_static_class_name == clone_constr->GetClassIdentifier(),
        "GetClassIdentifier gives " + clone_constr->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      // test WriteOutput Function
      // write the output
      copy_constr.WriteOutput( output_filename, nodes_to_output);
      // check that the output is correct
      BCL_Example_Check
      (
        io::File::FilesMatch( output_filename_correct, output_filename),
        "files don't match"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      BCL_MessageStd( "testing Write");
      WriteBCLObject( copy_constr);
       // read the file back
       BCL_MessageStd( "testing Read");
       cluster::OutputCenters< double, double> read_output_center;
       ReadBCLObject( read_output_center);

       // make sure that the written outputter and the read-in outputter are the same
       // write the output
       read_output_center.WriteOutput( output_filename, nodes_to_output);
       // check that the output is correct
       BCL_Example_Check
       (
         io::File::FilesMatch( output_filename_correct, output_filename),
         "files don't match"
       );

    //////////////////////
    // helper functions //
    //////////////////////

      // OutputNodeMembers is checked implicitly with WriteOutput

      return 0;
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

  }; //end ExampleClusterOutputRows

  // instantiate s_Instance
  const util::SiPtr< const util::ObjectInterface> ExampleClusterOutputCenters::DoubleCompare::s_Instance
  (
    GetObjectInstances().AddInstance( new DoubleCompare())
  );

  const ExampleClass::EnumType ExampleClusterOutputCenters::s_Instance
  (
    GetExamples().AddEnum( ExampleClusterOutputCenters())
  );

} // namespace bcl
