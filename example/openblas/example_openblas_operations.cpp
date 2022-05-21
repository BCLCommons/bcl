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

//c) Copyright BCL @ Vanderbilt University 2014
// (c) This file is part of the BCL software suite and is made available under license.
// (c) For commercial users:
// (c)
// (c)   The BCL copyright and license yields to non-BCL copyrights and licenses where indicated by code comments.
// include example header
#include "example.h"
// include the header of the class which this example is for
#include "openblas/bcl_openblas_operations.h"

// includes from bcl - sorted alphabetically
#include <math.h>

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_openblas_operations.cpp
  //!
  //! @author vuot2
  //! @date Mar 21, 2018
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOpenBlasOperations :
    public ExampleInterface
  {

  public:

    ExampleOpenBlasOperations *Clone() const
    {
      return new ExampleOpenBlasOperations( *this);
    }

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  /////////////////
  // data access //
  /////////////////

    int Run() const
    {
      // Test the double version
      openblas::Operations< double> blas_double;

      BCL_ExampleCheck( blas_double.GetClassIdentifier(), "bcl::openblas::Operations<double>");

      // create vectors and matrix for testing
      double d1[ 3] = { 1, 2, 3};
      double d2[ 3] = { 1, 5, 7};
      double d3[ 6] = { 1, 2, 3, 4, 5, 6};
      double d4[ 2] = { 32, 71};
      double d5[ 6] = { 7, 8, 9, 10, 11, 12};
      double d6[ 4] = { 58, 64, 139, 154};
      double d7[ 2] = { 37, 76};
      double d8[ 9] = { 1, 5, 7, 2, 10, 14, 3, 15, 21};

      linal::Vector< double> vd1( 3, d1);
      linal::Vector< double> vd2( 3, d2);
      linal::Matrix< double> md1( 2, 3, d3);
      linal::Vector< double> vd3( 2, d4);
      linal::Matrix< double> md2( 3, 2, d5);
      linal::Matrix< double> md3( 2, 2, d6);
      linal::Vector< double> vd4( 2, double(0) );
      linal::Vector< double> vd5( 2, double(5) );
      linal::Vector< double> vd6( 2, d7);
      linal::Vector< double> vd7( 2, double(0) );
      linal::Matrix< double> md4( 3, 3, d8);

      // Test the dot product
      BCL_ExampleCheckWithinAbsTolerance
      (
        blas_double.DotProduct( vd1, vd2), 32,
        float( 1.0e-6)
      );

      // Test the outer product
      BCL_ExampleCheck( blas_double.OuterProduct( vd1, vd2), md4);

      // Test the L2 norm of vector
      BCL_ExampleCheckWithinAbsTolerance
      (
        blas_double.Norm( vd2), std::sqrt( 75),
        float( 1.0e-6)
      );

      // Test the vector-matrix multiplication
      BCL_ExampleCheck( blas_double.Multiply( md1, vd2), vd3);

      // Test the matrix matrix multiplication
      BCL_ExampleCheck( blas_double.Multiply( md1, md2), md3);

      //VectorEqualsVectorTimesMatrix
      blas_double.VectorEqualsVectorTimesMatrix( vd4, vd2, md2);
      linal::VectorEqualsVectorTimesMatrix( vd7, vd2, md2);
      BCL_ExampleCheck( vd4, vd7);

      //VectorPlusEqualsMatrixTimesVector
      blas_double.VectorPlusEqualsMatrixTimesVector( vd5, md1, vd2);
      BCL_ExampleCheck( vd5, vd6);

      // Test the float version
      openblas::Operations< float> blas_float;

      BCL_ExampleCheck( blas_float.GetClassIdentifier(), "bcl::openblas::Operations<float>");

      // create vectors and matrix for testing
      float f1[ 3] = { 1, 2, 3};
      float f2[ 3] = { 1, 5, 7};
      float f3[ 6] = { 1, 2, 3, 4, 5, 6};
      float f4[ 2] = { 32, 71};
      float f5[ 6] = { 7, 8, 9, 10, 11, 12};
      float f6[ 4] = { 58, 64, 139, 154};
      float f7[ 2] = { 37, 76};
      float f8[ 9] = { 1, 5, 7, 2, 10, 14, 3, 15, 21};

      linal::Vector< float> vf1( 3, f1);
      linal::Vector< float> vf2( 3, f2);
      linal::Matrix< float> mf1( 2, 3, f3);
      linal::Vector< float> vf3( 2, f4);
      linal::Matrix< float> mf2( 3, 2, f5);
      linal::Matrix< float> mf3( 2, 2, f6);
      linal::Vector< float> vf4( 2, float(0) );
      linal::Vector< float> vf5( 2, float(5) );
      linal::Vector< float> vf6( 2, f7);
      linal::Vector< float> vf7( 2, float(0) );
      linal::Matrix< float> mf4( 3, 3, f8);

      // Test the dot product
      BCL_ExampleCheckWithinAbsTolerance
      (
        blas_float.DotProduct( vf1, vf2), 32,
        float( 1.0e-6)
      );

      // test the L2 norm of vector
      BCL_ExampleCheckWithinAbsTolerance
      (
        blas_float.Norm( vf2), std::sqrt( 75),
        float( 1.0e-6)
      );

      // Test the outer product
      BCL_ExampleCheck( blas_float.OuterProduct( vf1, vf2), mf4);
      
      // Test the vector-matrix multiplication
      BCL_ExampleCheck( blas_float.Multiply( mf1, vf2), vf3);

      // Test the matrix matrix multiplication
      BCL_ExampleCheck( blas_float.Multiply( mf1, mf2), mf3);

      //VectorEqualsVectorTimesMatrix
      blas_float.VectorEqualsVectorTimesMatrix( vf4, vf2, mf2);
      linal::VectorEqualsVectorTimesMatrix( vf7, vf2, mf2);
      BCL_ExampleCheck( vf4, vf7);

      //VectorPlusEqualsMatrixTimesVector
      blas_float.VectorPlusEqualsMatrixTimesVector( vf5, mf1, vf2);
      BCL_ExampleCheck( vf5, vf6);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleLinalOperations

  const ExampleClass::EnumType ExampleOpenBlasOperations::s_Instance
  (
    GetExamples().AddEnum( ExampleOpenBlasOperations())
  );

} // namespace bcl
