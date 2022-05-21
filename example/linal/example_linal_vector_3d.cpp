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
#include "linal/bcl_linal_vector_3d.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "linal/bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_angle.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "math/bcl_math_statistics.h"
#include "math/bcl_math_transformation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_vector_3d.cpp
  //! @details test the functionalities of linal::Vector3D
  //!
  //! @author karakam, woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalVector3D :
    public ExampleInterface
  {
  public:

    ExampleLinalVector3D *Clone() const
    {
      return new ExampleLinalVector3D( *this);
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
      linal::Vector3D default_vector3D;
      BCL_Example_Check
      (
        default_vector3D.X() == double( 0.0) &&
        default_vector3D.Y() == double( 0.0) &&
        default_vector3D.Z() == double( 0.0), "values should all be 0.0"
      );

      //initialize Vector3D from one double - sets all elements to this value
      linal::Vector3D vector1( double( 3.5));
      BCL_Example_Check
      (
        vector1.X() == double( 3.5) &&
        vector1.Y() == double( 3.5) &&
        vector1.Z() == double( 3.5), "values should all be 3.5"
      );

      // initialize with array
      const double c[] = { double( 4), double( -3), double( 7)};
      linal::Vector3D vector2( c);

      //randomize values
      math::Statistics::SetRand( vector1.Begin(), vector1.End(), 0.0, 1.0);

      BCL_MessageStd( "this is an random vector of three elements vector1");
      BCL_MessageStd( util::Format()( vector1));

      BCL_MessageStd( "this is vector of three elements vector2");
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

      // compute dihedral angle and output as formatted degree
      BCL_MessageStd( "compute dihedral angle and output as formatted degree");
      BCL_MessageStd
      (
        util::Format().W( 7).FFP( 3)( math::Angle::Degree( linal::Dihedral( vector1, vector2, vector2, vector1))) +
        std::string( 1, math::Angle::s_DegreeChar)
      );

      // generate two Vector3D as copies of vector1 and vector2
      linal::Vector3D v3( vector1), v4( vector2);

      // rotate v3 180 degree around X
      BCL_MessageStd( "Vector 3D");
      BCL_MessageStd( util::Format()( v3));
      BCL_MessageStd( "rotate 180 degree around X-axis");
      v3.Rotate( math::RotationMatrix3D( coord::GetAxes().e_X, math::Angle::Radian( 180.0)));
      BCL_MessageStd( util::Format()( v3));

      // initialize transformation matrix
      math::TransformationMatrix3D t1; //actually it is a matrix4D for a 3D space

      // add translation by -v4
      BCL_MessageStd( "transformation matrix gets a translation vector v4");

      BCL_MessageStd( util::Format()( t1));
      t1( -v4);
      BCL_MessageStd( util::Format()( t1));
      //add rotation
      BCL_MessageStd
      (
        "transformation matrix gets a translation vector rotation matrix around Z and 180 degree"
      );
      t1( coord::GetAxes().e_Z, math::Angle::Radian( 180.0));
      BCL_MessageStd( util::Format()( t1));

      // apply t1 on vector1
      BCL_MessageStd( "v3 is transformed");
      BCL_MessageStd( util::Format()( v3));
      v3.Transform( t1);
      BCL_MessageStd( util::Format()( v3));

      //generate a vector of a random translation in a sphere with the radius r=sphere_radius
      const double sphere_radius( 5);
      linal::Vector3D translation;
      for( size_t i( 0); i < 100; ++i)
      {
        translation.SetRandomTranslation( sphere_radius);
        BCL_Example_Check
        (
          translation.Norm() <= sphere_radius,
          "length of translation vector should be smaller than " + util::Format()( sphere_radius)
        );
      }

      // initialize three doubles
      const double x( 1.0);
      const double y( 2.0);
      const double z( 3.0);

      // construct an array of double
      double values[ 3] = { 2.0, 3.0, 4.0};

      // construct a vector
      linal::Vector< double> linal_vector( 3, values);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      linal::Vector3D vector_def;

      // constructor from 3 values
      linal::Vector3D vector_a( x, y, z);
      const linal::Vector3D vector_const_a( x, y, z);

      // constructor from a data pointer
      linal::Vector3D vector_b( values);
      const linal::Vector3D vector_const_b( values);

      // constructor from vector interface
      linal::Vector3D vector_c( linal_vector);

      // test clone
      util::ShPtr< linal::VectorInterface< double> > sp_vector( vector_b.Clone());

    /////////////////
    // data access //
    /////////////////

      // check class descriptor
      BCL_ExampleCheck( vector_def.GetClassIdentifier(), GetStaticClassName< linal::Vector3D>());
      BCL_ExampleCheck( sp_vector->GetClassIdentifier(), GetStaticClassName< linal::Vector3D>());

      // test const X(), Y() and Z()
      BCL_ExampleCheck( vector_const_a.X(), x);
      BCL_ExampleCheck( vector_const_a.Y(), y);
      BCL_ExampleCheck( vector_const_a.Z(), z);

      // test non-const X(), Y() and Z()
      const double diff( 1.5);
      vector_a.X() += diff;
      vector_a.Y() += diff;
      vector_a.Z() += diff;
      BCL_ExampleCheck( vector_a.X(), x + diff);
      BCL_ExampleCheck( vector_a.Y(), y + diff);
      BCL_ExampleCheck( vector_a.Z(), z + diff);

      // test const Begin() and End() functions
      const double expected_total( 9.0);
      double total( 0.0);
      for( const double *ptr( vector_const_b.Begin()), *ptr_end( vector_const_b.End()); ptr != ptr_end; ++ptr)
      {
        total += *ptr;
      }
      BCL_ExampleCheck( total, expected_total);

      // test non-const Begin() and End() functions
      total = 0.0;
      for( double *ptr( sp_vector->Begin()), *ptr_end( sp_vector->End()); ptr != ptr_end; ++ptr)
      {
        ( *ptr) += diff;
      }
      BCL_ExampleCheck( sp_vector->operator ()( 0), values[ 0] + diff);
      BCL_ExampleCheck( sp_vector->operator ()( 1), values[ 1] + diff);
      BCL_ExampleCheck( sp_vector->operator ()( 2), values[ 2] + diff);

      // test the values for the vector interface constructed vector3d
      BCL_ExampleCheck( vector_c( 0), linal_vector( 0));
      BCL_ExampleCheck( vector_c( 1), linal_vector( 1));
      BCL_ExampleCheck( vector_c( 2), linal_vector( 2));

    ///////////////
    // operators //
    ///////////////

      // test const operator( POS)
      BCL_ExampleCheck( vector_const_b( 0), values[ 0]);
      BCL_ExampleCheck( vector_const_b( 1), values[ 1]);
      BCL_ExampleCheck( vector_const_b( 2), values[ 2]);

      // test non-const X(), Y() and Z()
      vector_b( 0) += diff;
      vector_b( 1) += diff;
      vector_b( 2) += diff;
      BCL_ExampleCheck( vector_b( 0), values[ 0] + diff);
      BCL_ExampleCheck( vector_b( 1), values[ 1] + diff);
      BCL_ExampleCheck( vector_b( 2), values[ 2] + diff);

    //////////////////////
    // input and output //
    //////////////////////

      // write object
      WriteBCLObject( vector_a);
      linal::Vector3D vector_read;
      ReadBCLObject( vector_read);
      BCL_ExampleCheck( vector_read, vector_a);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleLinalVector3D

  const ExampleClass::EnumType ExampleLinalVector3D::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalVector3D())
  );

} // namespace bcl
