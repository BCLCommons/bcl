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
#include "linal/bcl_linal_matrix_operations.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_symmetric_eigensolver.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_matrix_operations.cpp
  //!
  //! @author woetzen, loweew
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalMatrixOperations :
    public ExampleInterface
  {
  public:

    ExampleLinalMatrixOperations *Clone() const
    {
      return new ExampleLinalMatrixOperations( *this);
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

    //! @brief create the matrix buy reading a matrix file with a dimension
    //! @param NAME : name of the file contain bcl table style matric in input/math
    //! @param SIZE : dimension of the square matrix
    linal::Matrix< double> MakeMatrix( const std::string &NAME, const size_t SIZE) const
    {
      linal::Matrix< double> matrix;
      const std::string filename
      (
        AddExampleInputPathToFilename( e_Math, NAME)
      );
      // get the score file
      io::IFStream input;
      io::File::TryOpenIFStream( input, filename);
      io::Serialize::Read( matrix, input);
      io::File::CloseClearFStream( input);
      return matrix;
    }

    //! @brief times the BCL's matrix matrix multiplication operation
    //! @param NAME : name of the file contain bcl table style matric in input/math
    //! @param SIZE : dimension of the square matrix
    void TimingMatrixOperation( const std::string &NAME, const size_t SIZE) const
    {
      // initiate all the matrix
      linal::Matrix< double> matrix = MakeMatrix( NAME, SIZE);

      // time the multiplication
      util::Stopwatch watch( false);
      watch.Reset();
      watch.Start();
      for( size_t i = 0; i < 1000; ++i) {
        matrix *matrix;
      }
      double time( watch.GetProcessDuration().GetSecondsFractional() / 1000.0);
      BCL_MessageStd
      (
        "matrix multiplication with dimension of " +
        util::Format()( SIZE) + "x" + util::Format()( SIZE) + " takes " +
        util::Format()( time) + " s."
      )
      watch.Stop();
    }

    //! @brief times the BCL's matrix matrix multiplication operation
    //! @param NAME : name of the file contain bcl table style matric in input/math
    //! @param SIZE : dimension of the square matrix
    void TimingMatrixVectorOperation( const std::string &NAME, const size_t SIZE) const
    {
      // initiate all the matrix and vector
      linal::Matrix< double> matrix( MakeMatrix( NAME, SIZE));
      linal::Vector< double> storage( SIZE, 0.0);
      linal::VectorConstReference< double> vec( SIZE, matrix[ 0]);
      // time the multiplication
      util::Stopwatch watch( false);
      watch.Reset();
      watch.Start();
      for( size_t i = 0; i < 1000; ++i) {
        linal::VectorPlusEqualsMatrixTimesVector( storage, matrix, vec);
      }
      double time( watch.GetProcessDuration().GetSecondsFractional() / 1000.0);
      BCL_MessageStd
      (
        "matrix vector multiplication with dimension of " +
        util::Format()( SIZE) + "x" + util::Format()( SIZE) + " takes " +
        util::Format()( time) + " s."
      )
      watch.Stop();
    }

    //! @brief times the BCL's matrix matrix multiplication operation
    //! @param NAME : name of the file contain bcl table style matric in input/math
    //! @param SIZE : dimension of the square matrix
    void TimingSVD( const std::string &NAME, const size_t SIZE) const
    {
      // initiate all the matrix
      linal::Matrix< double> matrix( MakeMatrix( NAME, SIZE));
      linal::Matrix< double> eig_v( SIZE, SIZE), eig_u( SIZE, SIZE);

      // time the SVD procedure
      util::Stopwatch watch( false);
      watch.Reset();
      watch.Start();
      for( size_t i = 0; i < 1000; ++i) {
        linal::SingularValueDecomposition( matrix, eig_v, eig_u);
      }
      double time( watch.GetProcessDuration().GetSecondsFractional() / 1000.0);
      BCL_MessageStd
      (
        "matrix multiplication with dimension of " +
        util::Format()( SIZE) + "x" + util::Format()( SIZE) + " takes " +
        util::Format()( time) + " s."
      )
      watch.Stop();
    }

    int Run() const
    {

      // creating data set
      float data1[ 9] = { 1.0, 2.0, 3.0, 2.0, 2.0, 4.0, 3.0, 4.0, 3.0};
      linal::Vector< float> vector( 3, 5.0);
      linal::Matrix< float> sqr_mat( 3, 3, 4.0);
      linal::Matrix< float> sym_mat( 3, 3, data1);
      linal::Matrix< float> eigenvec( 3, 3, 0.0f);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // testing replace diagonal
      linal::ReplaceDiagonal( sqr_mat, vector);
      BCL_ExampleCheck( sqr_mat( 2, 2), 5);

      // check find undefined rows
      const float nan_val( std::numeric_limits< float>::quiet_NaN());
      const float row_zero_and_two_nan[] = { 1.0, nan_val, 2.0, 3.0, nan_val, 4.0, 5.0, 6.0};
      BCL_ExampleCheck
      (
        linal::FindUndefinedRows( linal::Matrix< float>( size_t( 4), size_t( 2), row_zero_and_two_nan)),
        storage::Vector< size_t>::Create( 0, 2)
      );

      // helper variables
      static const double s_u[] = { 1.5, 2.5, 3.5};
      static const double s_v[] = { 4.0, 5.0, 6.0, 7.0};
      const linal::Vector< double> vec_u( 3, s_u);
      const linal::Vector< double> vec_v( 4, s_v);

      // outer product of two vectors
      linal::Matrix< double> outer_product_uv( vec_u.GetSize(), vec_v.GetSize(), 0.0);
      linal::AddOuterProductToMatrix( outer_product_uv, vec_u, vec_v);
      static const double sum_outer_product_uv_expected( 165.0);
      BCL_MessageDbg
      (
        "outer product of vector u and v:\n" + util::Format()( vec_u) + "\n*\n" + util::Format()( vec_v)
      );
      BCL_MessageStd
      (
        "this is the outer product u*v result: " + util::Format()( outer_product_uv)
      );
      BCL_ExampleCheck( outer_product_uv.Sum(), sum_outer_product_uv_expected);

      // try the matrix multiply with transpose functions, which alleviate the need to transpose the matrix (which may
      // barely fit in memory to begin with) for the purposes of multiplication
      linal::Matrix< double> random_5_by_4_matrix( 5, 4, double( 0.0));
      linal::Matrix< double> random_5_by_8_matrix( 5, 8, double( 0.0));
      random_5_by_4_matrix.AsVector().SetRand( 0.0, 10.0);
      random_5_by_8_matrix.AsVector().SetRand( 0.0, 10.0);

      linal::Matrix< double> explicit_mt_times_m( 4, 4), explicit_m_times_mt( 5, 5), explicit_mt_times_n( 4, 8);
      linal::MultiplySmallDenseMatrices( explicit_mt_times_m, random_5_by_4_matrix.Transposed(), random_5_by_4_matrix);
      linal::MultiplySmallDenseMatrices( explicit_m_times_mt, random_5_by_4_matrix, random_5_by_4_matrix.Transposed());
      linal::MultiplySmallDenseMatrices( explicit_mt_times_n, random_5_by_4_matrix.Transposed(), random_5_by_8_matrix);
      BCL_ExampleCheckWithinAbsTolerance
      (
        linal::MatrixTransposeTimesMatrix( random_5_by_4_matrix),
        explicit_mt_times_m,
        0.00001
      );
      BCL_ExampleCheckWithinAbsTolerance
      (
        linal::MatrixTimesItselfTransposed( random_5_by_4_matrix),
        explicit_m_times_mt,
        0.00001
      );
      BCL_ExampleCheckWithinAbsTolerance
      (
        linal::MatrixTransposeTimesMatrix( random_5_by_4_matrix, random_5_by_8_matrix),
        explicit_mt_times_n,
        0.00001
      );

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

      TimingSVD( "matrix.100", 100);
      TimingMatrixOperation( "matrix.100", 100);
      TimingMatrixVectorOperation( "matrix.100", 100);

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleLinalMatrixOperations

  const ExampleClass::EnumType ExampleLinalMatrixOperations::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalMatrixOperations())
  );

} // namespace bcl
