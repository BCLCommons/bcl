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
#include "storage/bcl_storage_triplet.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_storage_triplet.cpp
  //!
  //! @author alexanns
  //! @date Nov 06, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleStorageTriplet :
    public ExampleInterface
  {
  public:

    ExampleStorageTriplet *Clone() const
    {
      return new ExampleStorageTriplet( *this);
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
      BCL_MessageStd
      (
        "this class is similar to the standard Library std::pair. Except that it is a triplet, it has the typical bcl Clone, Empty and read/write functionality."
      );
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      BCL_MessageStd( "Construct Triplet of int, double and string from int, double and string");
      storage::Triplet< int, double, std::string> residue;

      // test construction from three values and accessing the vaules by First, Second, and Third functions
      BCL_MessageStd
      (
        "construct from three values  " + util::Format()( 1) + ", " + util::Format()( 'A') + ", alpha"
      );
      storage::Triplet< int, double, std::string> start( 1, -1.5, "alpha");
      BCL_MessageStd( "test access to elements  " + util::Format()( 1) + ", -1.5, alpha");
      BCL_MessageStd( "start: first element:  " + util::Format()( start.First()));
      BCL_MessageStd( "start: second element: " + util::Format()( start.Second()));
      BCL_MessageStd( "start: third element: " + start.Third());
      BCL_Example_Check
      (
        start.First() == 1 && start.Second() == -1.5 && start.Third() == "alpha",
        "construction with three values and accessing with First,  Second, and Third functions did not work"
      );

      // test copy constructor
      BCL_MessageStd( "copy constructor: copy start  ");
      storage::Triplet< int, double, std::string> copy( start);
      BCL_MessageStd( "copy: first element:  " + util::Format()( copy.First()));
      BCL_MessageStd( "copy: second element: " + util::Format()( copy.Second()));
      BCL_MessageStd( "copy: third element: " + copy.Third());
      BCL_Example_Check
      (
        copy.First() == 1 && copy.Second() == -1.5 && copy.Third() == "alpha",
        "copy constructor did not work"
      );

    /////////////////
    // data access //
    /////////////////

      // test GetClassIdentifier function
      const std::string class_identifier( residue.GetClassIdentifier());
      const std::string class_identifier_expected( "bcl::storage::Triplet<int,double,std::string>");

      BCL_MessageStd( "GetClassIdentifier: " + class_identifier);
      BCL_Example_Check
      (
        class_identifier == class_identifier_expected,
        "GetClassIdentifier does not return the correct name: " +
        class_identifier_expected + " but " + class_identifier
      );

      // test data manipulation of First
      BCL_MessageStd( "test data manipulation of First: ");
      BCL_MessageStd( "start: first element: before " + util::Format()( start.First()));
      start.First() = 2;
      BCL_MessageStd( "start: first element: after " + util::Format()( start.First()));
      BCL_Example_Check
      (
        start.First() == 2, "data manipulation of First did not work correctly"
      );

      // test data manipulation of Second
      BCL_MessageStd( "test data manipulation of Second: ");
      BCL_MessageStd( "start: Second element: before " + util::Format()( start.Second()));
      start.Second() = -2.5;
      BCL_MessageStd( "start: Second element: after " + util::Format()( start.Second()));
      BCL_Example_Check
      (
        start.Second() == -2.5, "data manipulation of Second did not work correctly"
      );

      // test data manipulation of Third
      BCL_MessageStd( "test data manipulation of Third: ");
      BCL_MessageStd( "start: Third element: before " + util::Format()( start.Third()));
      start.Third() = "beta";
      BCL_MessageStd( "residue: Third element: after " + util::Format()( start.Third()));
      BCL_Example_Check
      (
        start.Third() == "beta", "data manipulation of Third did not work correctly"
      );

    ///////////////
    // operators //
    ///////////////

      // test operator == in falst scenario
      BCL_MessageStd
      (
        "test operator == with false scenario: " + util::Format()( start == copy)
      );
      BCL_Example_Check
      (
        !( copy == start), "operator == returns true when it should return false"
      );

      // test equal = operator
      BCL_MessageStd( "test operator = : start = copy");
      BCL_MessageStd( "start:  first element:  " + util::Format()( start.First()));
      BCL_MessageStd( "start: second element: " + util::Format()( start.Second()));
      BCL_MessageStd( "start: third element: " + start.Third());
      BCL_MessageStd( "copy: first element:  " + util::Format()( copy.First()));
      BCL_MessageStd( "copy: second element: " + util::Format()( copy.Second()));
      BCL_MessageStd( "copy: third element: " + copy.Third());
      start = copy;
      BCL_MessageStd( "after: start: first element:  " + util::Format()( start.First()));
      BCL_MessageStd( "after: start: second element: " + util::Format()( start.Second()));
      BCL_MessageStd( "start: third element: " + start.Third());
      BCL_Example_Check
      (
        start.First() == copy.First() && start.Second() == copy.Second() && start.Third() == copy.Third(),
        "operator = did not work properly"
      );

      // test operator == in true scenario
      BCL_MessageStd
      (
        "test operator == with true scenario: " + util::Format()( copy == start)
      );
      BCL_Example_Check
      (
        copy == start, "operator == returns false when it should return true"
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      // output
      BCL_MessageStd( "test output: this is the Triplet: " + util::Format()( start));

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleStorageTriplet

  const ExampleClass::EnumType ExampleStorageTriplet::s_Instance
  (
    GetExamples().AddEnum( ExampleStorageTriplet())
  );

} // namespace bcl
