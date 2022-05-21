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
#include "linal/bcl_linal_vector_reference.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_vector_reference.cpp
  //!
  //! @author mendenjl
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalVectorReference :
    public ExampleInterface
  {
  public:

    ExampleLinalVectorReference *Clone() const
    {
      return new ExampleLinalVectorReference( *this);
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
      // dimensions
      const size_t size( 3);

      // generate two Vector< double>
      double c[] = { double( 4), double( -3), double( 7)};

      // linal::Vector< double> to use as original matrix to reference
      linal::Vector< double> vector_orig( size, double( 3));
      linal::Vector< double> vector_c( size, c);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // checking GetSize() & default constructor
      BCL_ExampleCheck( linal::VectorReference< double>().GetSize(), 0);

      //! constructing from rows, cols, pointer to data
      linal::VectorReference< double> vector_b( 3, vector_orig.Begin());

      //! constructing from VectorInterface
      linal::VectorReference< double> vector_ref_c( vector_c);

    /////////////////
    // data access //
    /////////////////

      //! checking GetNumberCols()
      BCL_ExampleCheck( vector_b.GetSize(), 3);

      //! checking GetNumberOfElements()
      BCL_ExampleCheck( vector_ref_c.GetSize(), vector_c.GetSize());

      //! checking construction
      BCL_ExampleCheck( std::equal( vector_ref_c.Begin(), vector_ref_c.End(), vector_c.Begin()), true);

      //! checking Begin() and End()
      BCL_ExampleCheck( vector_c.Begin() + vector_c.GetSize(), vector_c.End());

    ///////////////
    // operators //
    ///////////////

      //! checking operator()( row, col)
      BCL_ExampleCheck( vector_ref_c( 1), vector_c( 1));

      //! checking operator = ( VectorInterface
      vector_ref_c = 1.0;
      BCL_ExampleCheck( vector_ref_c( 0), 1.0);

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleLinalVectorReference

  const ExampleClass::EnumType ExampleLinalVectorReference::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalVectorReference())
  );

} // namespace bcl
