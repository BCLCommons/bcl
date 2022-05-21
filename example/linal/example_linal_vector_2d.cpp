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
#include "linal/bcl_linal_vector_2d.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "linal/bcl_linal_vector_2d_operations.h"
#include "math/bcl_math_angle.h"
#include "math/bcl_math_statistics.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_vector_2d.cpp
  //! @details test the functionalities of linal::Vector2D
  //!
  //! @author mendenjl
  //! @date Dec 05, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalVector2D :
    public ExampleInterface
  {
  public:

    ExampleLinalVector2D *Clone() const
    {
      return new ExampleLinalVector2D( *this);
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
    //////////////////
    // constructors //
    //////////////////

      //default constructor
      linal::Vector2D default_vector2D;
      BCL_ExampleCheck( linal::Vector2D(), double( 0.0));

      //initialize Vector2D from one double - sets all elements to this value
      linal::Vector2D vector1( double( 3.5));
      BCL_ExampleCheck( vector1, double( 3.5));

      // initialize with array
      const double c[] = { double( 4), double( -3)};
      linal::Vector2D vector2( c);

      //randomize values
      math::Statistics::SetRand( vector1.Begin(), vector1.End(), 0.0, 1.0);

      BCL_MessageStd( "this is an random vector of two elements vector1");
      BCL_MessageStd( util::Format()( vector1));

      BCL_MessageStd( "this is vector of two elements vector2");
      BCL_MessageStd( util::Format()( vector2));

      // compute scalar product
      BCL_MessageStd( "this is the scalar product of the vectors, using two different functions");
      BCL_MessageStd
      (
        "linal::ScalarProduct( vector1, vector2): " + util::Format()( linal::ScalarProduct( vector1, vector2))
      );
      BCL_MessageStd( "vector1 * vector2: " + util::Format()( vector1 * vector2));

      // compute projection angle and output as formatted degree
      BCL_MessageStd
      (
        "compute projection angle between vector2 and vector1 and output as formatted degree"
      );
      BCL_MessageStd
      (
        util::Format().W( 7).FFP( 3)( math::Angle::Degree( linal::ProjAngle( vector2, vector1))) +
        std::string( 1, math::Angle::s_DegreeChar)
      );

      // initialize two doubles
      const double x( 1.0);
      const double y( 2.0);

      // construct an array of double
      double values[ 2] = { 2.0, 3.0};

      // construct a vector
      linal::Vector< double> linal_vector( 2, values);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      linal::Vector2D vector_def;

      // constructor from 3 values
      linal::Vector2D vector_a( x, y);
      const linal::Vector2D vector_const_a( x, y);

      // constructor from a data pointer
      linal::Vector2D vector_b( values);
      const linal::Vector2D vector_const_b( values);

      // constructor from vector interface
      linal::Vector2D vector_c( linal_vector);

      // test clone
      util::ShPtr< linal::VectorInterface< double> > sp_vector( vector_b.Clone());

    /////////////////
    // data access //
    /////////////////

      // check class descriptor
      BCL_ExampleCheck( vector_def.GetClassIdentifier(), GetStaticClassName< linal::Vector2D>());
      BCL_ExampleCheck( sp_vector->GetClassIdentifier(), GetStaticClassName< linal::Vector2D>());

      // test const X(), Y() and Z()
      BCL_ExampleCheck( vector_const_a.X(), x);
      BCL_ExampleCheck( vector_const_a.Y(), y);

      // test non-const X(), Y() and Z()
      const double diff( 1.5);
      vector_a.X() += diff;
      vector_a.Y() += diff;
      BCL_ExampleCheck( vector_a.X(), x + diff);
      BCL_ExampleCheck( vector_a.Y(), y + diff);

      // test const Begin() and End() functions
      double total( 0.0);
      for( const double *ptr( vector_const_b.Begin()), *ptr_end( vector_const_b.End()); ptr != ptr_end; ++ptr)
      {
        total += *ptr;
      }
      BCL_ExampleCheck( total, 5.0);

      // test non-const Begin() and End() functions
      total = 0.0;
      for( double *ptr( sp_vector->Begin()), *ptr_end( sp_vector->End()); ptr != ptr_end; ++ptr)
      {
        ( *ptr) += diff;
      }
      BCL_ExampleCheck( sp_vector->operator ()( 0), values[ 0] + diff);
      BCL_ExampleCheck( sp_vector->operator ()( 1), values[ 1] + diff);

      // test the values for the vector interface constructed vector2d
      BCL_ExampleCheck( vector_c( 0), linal_vector( 0));
      BCL_ExampleCheck( vector_c( 1), linal_vector( 1));

    ///////////////
    // operators //
    ///////////////

      // test const operator( POS)
      BCL_ExampleCheck( vector_const_b( 0), values[ 0]);
      BCL_ExampleCheck( vector_const_b( 1), values[ 1]);

      // test non-const X(), Y() and Z()
      vector_b( 0) += diff;
      vector_b( 1) += diff;
      BCL_ExampleCheck( vector_b( 0), values[ 0] + diff);
      BCL_ExampleCheck( vector_b( 1), values[ 1] + diff);

    //////////////////////
    // input and output //
    //////////////////////

      // write object
      WriteBCLObject( vector_a);
      linal::Vector2D vector_read;
      ReadBCLObject( vector_read);
      BCL_ExampleCheck( vector_read, vector_a);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleLinalVector2D

  const ExampleClass::EnumType ExampleLinalVector2D::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalVector2D())
  );

} // namespace bcl
