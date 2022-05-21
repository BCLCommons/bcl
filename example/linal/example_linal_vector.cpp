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
#include "linal/bcl_linal_vector.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_vector.cpp
  //!
  //! @author woetzen
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalVector :
    public ExampleInterface
  {
  public:

    ExampleLinalVector *Clone() const
    {
      return new ExampleLinalVector( *this);
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
      // generate two Vector< double>
      double c[] = { double( 4), double( -3), double( 7)};

    //////////////////
    // construction //
    //////////////////

      // default constructed - empty vector
      linal::Vector< double> v0;

      // constructed from length - initialize all values to double (0.0);
      linal::Vector< double> v1( 3);

      // construct from length and pointer to data
      linal::Vector< double> v2( 3, c);

      // copy constructor
      linal::Vector< double> v3( v2);

      BCL_MessageStd( "this is an empty vector");
      BCL_MessageStd( util::Format()( v0));
      BCL_Example_Check
      (
        v0.GetSize() == 0,
        "default constructed vector is not empty, but has " + util::Format()( v0.GetSize()) + " elements!"
      );

      BCL_MessageStd( "this is a vector of three elements 0.0\n" + util::Format()( v1));
      BCL_Example_Check
      (
        v1.GetSize() == 3,
        "vector should have 3 elements, but has " + util::Format()( v1.GetSize()) + " elements!"
      );
      BCL_Example_Check
      (
        v1 == double( 0.0),
        "vector1 should only have elements 0.0"
      );

      BCL_MessageStd( "this is a vector initialized from an array\n" + util::Format()( v2));
      BCL_ExampleCheck( v2.GetSize(), 3);
      BCL_Example_Check
      (
        std::equal( v2.Begin(), v2.End(), c),
        "original data does not agree with data in vector"
      );

      BCL_MessageStd( "this is vector 3 copied from vector 2:\n" + util::Format()( v3));
      BCL_ExampleCheck( v2, v3);

    /////////////////
    // data access //
    /////////////////

      // names
      BCL_MessageStd( "this is the class identifier: " + v3.GetClassIdentifier());
      BCL_MessageStd( "this is the static class identifier: " + GetStaticClassName( v3));
      BCL_Example_Check
      (
        v3.GetClassIdentifier() == GetStaticClassName( v3), "identifier and class name should match"
      );

      // size
      BCL_MessageStd( "size of v3: " + util::Format()( v3.GetSize()));

      // begin and end
      BCL_MessageStd( "content of pointer on begin of v3: " + util::Format()( *v3.Begin()));
      BCL_Example_Check
      (
        v3.End() == v3.Begin() + v3.GetSize(), "end does not point to begin + size"
      );

    ////////////////
    // operations //
    ////////////////

      // norm and square norm = length and square of length
      const double norm( v3.Norm());
      const double square_norm( v3.SquareNorm());
      BCL_MessageStd
      (
        "norm and square norm of vector3: " + util::Format()( norm) + " " + util::Format()( square_norm)
      );

    ///////////////
    // operators //
    ///////////////

      // access operator()
      BCL_MessageStd( "this is the first element in v2(0): " + util::Format()( v2( 0)));
      BCL_Example_Check
      (
        v2( 0) == c[ 0], "first element in v2 does not match original array"
      );

      // assignment operator
      v1 = v2;
      BCL_MessageStd( "v1 after assigning it from v2: " + util::Format()( v1));
      BCL_Example_Check
      (
        v1 == v2, "assignment v1 = v2 was not successful"
      );

      // assign value
      v1 = double( 5.0);
      BCL_MessageStd( "v1 after assigning it from 5.0: " + util::Format()( v1));
      BCL_Example_Check
      (
        v1 == 5.0, "assignment v1 = 5.0 was not successful"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // filename
      BCL_MessageStd( "writing and reading from file");
      WriteBCLObject( v3);
      ReadBCLObject( v0);
      BCL_Example_Check
      (
        v0 == v3,
        "written and read vector are not identical" + util::Format()( v3) + " " + util::Format()( v0)
      );

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleLinalVector

  const ExampleClass::EnumType ExampleLinalVector::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalVector())
  );

} // namespace bcl
