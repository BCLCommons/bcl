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
#include "cluster/bcl_cluster_distances_stored.h"

// includes from bcl - sorted alphabetically
#include "cluster/bcl_cluster_input_table.h"
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_cluster_distances_stored.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleClusterDistancesStored :
    public ExampleInterface
  {
  public:

    ExampleClusterDistancesStored *Clone() const
    {
      return new ExampleClusterDistancesStored( *this);
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
      // instantiate the s_Instance that will be needed for this example
      cluster::DistancesStored< std::string, double>::s_Instance.IsDefined();

      // create string "file" and initialize with the path to the example table that will be used for input
      std::string file( AddExampleInputPathToFilename( e_Biology, "table_example.txt"));

      // create ifstream read
      io::IFStream read;

      io::File::CloseClearFStream( read);

      // open "file" for reading
      BCL_ExampleMustOpenInputFile( read, file);

      cluster::InputTable< double> input_table( true, false);

      // create data construct to hold all the distances for the example
      util::ShPtr< cluster::DistancesStored< std::string, double> > constr( input_table.HandleInput( read));

      storage::HashMap< size_t, storage::HashMap< size_t, double> > distances( constr->GetData());
      util::ShPtr< storage::List< std::string> > object_list( input_table.GetInputObjects());

      BCL_ExampleAssert( distances.GetSize(), 3);

      // make sure that "object_list" has the right size and is defined
      BCL_ExampleAssert( object_list.IsDefined(), true);
      BCL_ExampleCheck( object_list->GetSize(), 4);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      {
        // check the default constructor
        BCL_MessageStd( "Testing default constructor");

        // create DistancesStored "def_constr"
        util::ShPtr< cluster::DistancesStored< std::string, double> > constr( input_table.HandleInput( read));

        // check that "constr" was constructed properly
        BCL_ExampleIndirectCheck( constr->GetData().GetSize(), 0, "Constructor from empty table");
      }

      {
        // create DistancesStored "constr"
        cluster::DistancesStored< std::string, double> constr( distances);

        // check that "constr" constructed properly
        BCL_ExampleIndirectAssert( constr.GetData().GetSize(), 3, "Constructor from hashmap");
      }

      {
        // check copy constructor
        cluster::DistancesStored< std::string, double> dummy( distances);
        cluster::DistancesStored< std::string, double> constr( dummy);

        // check that "constr" constructed properly
        BCL_ExampleIndirectCheck( constr.GetData().GetSize(), 3, "Copy constructor");
      }

      {
        // check clone constructor
        cluster::DistancesStored< std::string, double> dummy( distances);
        util::ShPtr< cluster::DistancesStored< std::string, double> > constr( dummy.Clone());

        // check that "constr" was constructed properly
        BCL_ExampleIndirectCheck( constr->GetData().GetSize(), 3, "Clone");
      }

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck
      (
        constr->GetClassIdentifier(),
        ( GetStaticClassName< cluster::DistancesStored< std::string, double> >())
      );

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      // test the parentheses operator
      BCL_MessageStd( "test the parentheses operator");
      {
        const double correct_value( 1.2);
        storage::Pair< double, bool> result
        (
          CallParenthesesOperator( *object_list->Begin(), *--object_list->End(), *constr, correct_value)
        );
        BCL_Example_Check
        (
          result.Second(),
          "operator() gave " + util::Format()( result.First()) + " but should give " + util::Format()( correct_value)
        );
      }
      {
        const double correct_value( 6.0);
        storage::Pair< double, bool> result
        (
          CallParenthesesOperator( *object_list->Begin(), *----object_list->End(), *constr, correct_value)
        );
        BCL_Example_Check
        (
          result.Second(),
          "operator() gave " + util::Format()( result.First()) + " but should give " + util::Format()( correct_value)
        );
      }
      {
        const double correct_value( 3.5);
        storage::Pair< double, bool> result
        (
          CallParenthesesOperator( *++object_list->Begin(), *--object_list->End(), *constr, correct_value)
        );
        BCL_Example_Check
        (
          result.Second(),
          "operator() gave " + util::Format()( result.First()) + " but should give " + util::Format()( correct_value)
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

      // make sure that the written distances and the read-in distances are the same
      std::stringstream input_output;
      input_output << *constr;
      cluster::DistancesStored< std::string, double> read_distances;
      input_output >> read_distances;

      BCL_ExampleCheck( constr->GetData(), read_distances.GetData());

    //////////////////////
    // helper functions //
    //////////////////////

      // end
      return 0;
    }

    // @brief CallParenthesesOperator calls the parentheses operator and makes sure it works
    // @param STRING_A the first string which will be used to call the operator()
    // @param STRING_B the second string which will be used to call the operator()
    // @param DISTANCES_STORED the DistancesStored object whose operator() will be tested
    // @param EXPECTED_VALUE double which is the value expected to be returned by the operator()
    storage::Pair< double, bool> CallParenthesesOperator
    (
      const util::SiPtr< const std::string> &STRING_A,
      const util::SiPtr< const std::string> &STRING_B,
      const cluster::DistancesStored< std::string, double> &DISTANCES_STORED,
      const double EXPECTED_VALUE
    ) const;

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleClusterDistancesStored

  const ExampleClass::EnumType ExampleClusterDistancesStored::s_Instance
  (
    GetExamples().AddEnum( ExampleClusterDistancesStored())
  );

  // @brief CallParenthesesOperator calls the parentheses operator and makes sure it works
  // @param STRING_A the first string which will be used to call the operator()
  // @param STRING_B the second string which will be used to call the operator()
  // @param DISTANCES_STORED the DistancesStored object whose operator() will be tested
  // @param EXPECTED_VALUE double which is the value expected to be returned by the operator()
  storage::Pair< double, bool> ExampleClusterDistancesStored::CallParenthesesOperator
  (
    const util::SiPtr< const std::string> &STRING_A,
    const util::SiPtr< const std::string> &STRING_B,
    const cluster::DistancesStored< std::string, double> &DISTANCES_STORED,
    const double EXPECTED_VALUE
  ) const
  {
    BCL_MessageStd
    (
      "Testing parentheses operator with " + util::Format()( *STRING_A) + " and " + util::Format()( *STRING_B)
    );

    BCL_MessageDbg( "Distances stored object is\n" + util::Format()( DISTANCES_STORED));

    // create double "value" and initialize with value retrieved from DATA
    const double value
    (
      DISTANCES_STORED( storage::VectorND< 2, util::SiPtr< const std::string> >( STRING_A, STRING_B))
    );

    // make sure "value" matched "EXPECTED_VALUE"
    return storage::Pair< double, bool>( value, math::EqualWithinTolerance( EXPECTED_VALUE, value));
  }

} // namespace bcl
