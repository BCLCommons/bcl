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
#include "linal/bcl_linal_matrix.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_matrix.cpp
  //!
  //! @author woetzen, putnamdk, fischea
  //! @date 10/13/2012
  //! @remarks status complete
  //! @remarks reviewed by putnamdk, fischea on 12/18/2012
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalMatrix :
    public ExampleInterface
  {
  public:

    ExampleLinalMatrix *Clone() const
    {
      return new ExampleLinalMatrix( *this);
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
      double c[] = { double( 4), double( -3), double( 7), double( 2)};

    //////////////////
    // construction //
    //////////////////

      // instantiate the s_Instance needed for this example
      linal::Matrix< double>::s_Instance.IsDefined();

      // default constructed - empty matrix
      linal::Matrix< double> matrix_empty;

      // constructed from dimensions - initialize all values to double (0.0);
      linal::Matrix< double> matrix_2x2( 2, 2);

      // construct from length and pointer to data
      linal::Matrix< double> matrix_array( 2, 2, c);

      // copy constructor
      linal::Matrix< double> matrix_copy( matrix_array);

      // clone
      util::ShPtr< linal::MatrixInterface< double> > ptr_matrix( matrix_array.Clone());

      // construct from const interface with paddings and fill value
      linal::Matrix< double> matrix_padding( matrix_2x2, 1, 2, 0.0);

      // check empty constructor
      BCL_MessageStd( "this is an empty matrix");
      BCL_MessageStd( util::Format()( matrix_empty));
      BCL_Example_Check
      (
        matrix_empty.GetNumberOfElements() == 0 && matrix_empty.GetNumberRows() == 0 && matrix_empty.GetNumberCols() == 0,
        "default constructed matrix is not empty, but has " + util::Format()( matrix_empty.GetNumberOfElements()) + " elements!"
      );

      // check 2x2 matrix with standard elements
      BCL_MessageStd( "this is a matrix of 2x2 elements 0.0");
      BCL_MessageStd( util::Format()( matrix_2x2));
      BCL_Example_Check
      (
        matrix_2x2.GetNumberOfElements() == 4 && matrix_2x2.GetNumberRows() == 2 && matrix_2x2.GetNumberCols() == 2,
        "matrix should have 2x2 elements, but has " + util::Format()( matrix_2x2) + " elements!"
      );
      BCL_Example_Check
      (
        matrix_2x2 == double( 0.0),
        "matrix1 should only have elements 0.0"
      );

      // check 2x2 matrix initialized from an array
      BCL_MessageStd( "this is a 2x2 matrix initialized from an array");
      BCL_MessageStd( util::Format()( matrix_array));
      BCL_Example_Check
      (
        matrix_2x2.GetNumberOfElements() == 4 && matrix_2x2.GetNumberRows() == 2 && matrix_2x2.GetNumberCols() == 2,
        "matrix2 should have 2x2 elements, but has " + util::Format()( matrix_2x2) + " elements!"
      );
      BCL_Example_Check
      (
        std::equal( matrix_array.Begin(), matrix_array.End(), c),
        "original data does not agree with data in matrix"
      );

      // check the copy constructor
      BCL_MessageStd( "this is matrix3 copied from matrix2:\n" + util::Format()( matrix_copy));
      BCL_Example_Check
      (
        matrix_array == matrix_copy, "matrix2 and matrix3 should be identical"
      );

      // check the clone
      BCL_MessageStd( "this is matrix4 cloned from matrix2:\n" + util::Format()( ptr_matrix));
      BCL_Example_Check
      (
        matrix_array == *ptr_matrix, "matrix2 and matrix4 should be identical"
      );

      // check constructor with paddings and fill value
      BCL_MessageStd
      (
        "this is matrix5 constructed from matrix2 with paddings:\n" + util::Format()( matrix_padding)
      );
      BCL_Example_Check
      (
        matrix_padding.GetNumberRows() == matrix_2x2.GetNumberRows() + 1 &&
        matrix_padding.GetNumberCols() == matrix_2x2.GetNumberCols() + 2,
        "matrix_padding should have " + util::Format()( matrix_2x2.GetNumberRows() + 1) + " rows and " +
        util::Format()( matrix_2x2.GetNumberCols() + 2) + " columns"
      );
      BCL_Example_Check
      (
        matrix_padding == ( double) 0.0,
        "matrix5 should only have elements 0.0"
      );

    /////////////////
    // data access //
    /////////////////

      // names
      BCL_ExampleCheck( ptr_matrix->GetClassIdentifier(), GetStaticClassName( matrix_array));

      // size
      BCL_ExampleIndirectCheck( matrix_array.GetNumberOfElements(), 4, "Constructor from array");
      BCL_ExampleCheck( matrix_array.GetNumberRows(), 2);
      BCL_ExampleCheck( matrix_array.GetNumberCols(), 2);

      BCL_MessageStd
      (
        util::Format()( matrix_array)
      )

      // pointers
      BCL_ExampleCheck( *matrix_array.Begin(), 4.0);

    ////////////////
    // operations //
    ////////////////

      // check method append
      storage::Vector< linal::Matrix< double> > matrices_to_append;
      matrix_array.Append( matrices_to_append);

      // check append with no matrices to append
      BCL_Example_Check
      (
        matrix_array.GetNumberRows() == 2, "number of rows should be 2"
      );

      // append an empty vector of matrices onto an empty matrix
      matrix_empty.Append( matrices_to_append);

      BCL_MessageStd( "starting with empty matrix and append: " + util::Format()( matrix_empty));
      BCL_ExampleCheck( linal::Matrix< double>().GetNumberRows(), 0);

      matrices_to_append.PushBack( matrix_array);
      matrices_to_append.PushBack( matrix_array);
      matrices_to_append.PushBack( matrix_array);

      matrix_array.Append( matrices_to_append);

      BCL_MessageStd( "appended matrices: \n" + util::Format()( matrix_array));

      BCL_Example_Check
      (
        matrix_array.GetNumberRows() == 8, "number of rows should be 8"
      );

      BCL_Example_Check
      (
        matrix_array.GetNumberCols() == 2, "number of cols should be 2"
      );

      // append an empty vector of matrices onto an empty matrix
      linal::Matrix< double> mtrx_empty;
      mtrx_empty.Append( matrices_to_append);
      BCL_MessageStd
      (
        "test: starting with empty matrix and append: " + util::Format()( mtrx_empty)
      );

      linal::Matrix< double> is_not_square( 3, 6);
      BCL_Example_Check
      (
        is_not_square.IsSquare() == false, "matrix is_not_square is square"
      )

      linal::Matrix< double> is_square( 5, 5);
      BCL_Example_Check
      (
        is_square.IsSquare() == true, "matrix is_square is not square"
      )

      linal::Matrix< double> diagonal( 3, 3, 0.0);
      diagonal( 0, 0) = 1.0;
      diagonal( 1, 1) = 2.0;
      diagonal( 2, 2) = 3.0;

      linal::Matrix< double> not_diagonal( 3, 3, 0.0);
      not_diagonal( 0, 0) = 1.0;
      not_diagonal( 1, 1) = 2.0;
      not_diagonal( 2, 2) = 3.0;
      not_diagonal( 0, 2) = 0.5;

      BCL_Example_Check
      (
        diagonal.IsDiagonal() == true, "matrix named diagonal is diagonal"
      );

      BCL_Example_Check
      (
        not_diagonal.IsDiagonal() == false,
        "matrix named not_diagonal is not diagonal "
      );

      linal::Matrix< double> symmetric( 3, 3, 0.0);
      symmetric( 0, 0) = 1.0;
      symmetric( 1, 1) = 2.0;
      symmetric( 2, 2) = 3.0;
      symmetric( 2, 0) = 4.0;
      symmetric( 0, 2) = 4.0;

      linal::Matrix< double> not_symmetric( 3, 3, 0.0);
      not_symmetric( 0, 0) = 1.0;
      not_symmetric( 0, 1) = 4.0;
      not_symmetric( 0, 2) = 0.0;

      not_symmetric( 1, 0) = 2.0;
      not_symmetric( 1, 1) = 1.0;
      not_symmetric( 1, 2) = 5.0;

      not_symmetric( 2, 0) = 0.0;
      not_symmetric( 2, 1) = 3.0;
      not_symmetric( 2, 2) = 1.0;

      BCL_Example_Check
      (
        symmetric.IsSymmetric() == true, "matrix named symmetric is not symmetric"
      );

      BCL_Example_Check
      (
        not_symmetric.IsSymmetric() == false,
        "matrix named not_symmetric is symmetric "
      );

      linal::Matrix< double> &tri_diagonal( not_symmetric);

      BCL_Example_Check
      (
        tri_diagonal.IsTriDiagonal() == true,
        "matrix named tri_diagonal is not tri_diagonal"
      );

      linal::Matrix< double> &not_tri_diagonal( not_symmetric);
      not_tri_diagonal( 0, 2) = -0.5;

      BCL_Example_Check
      (
        not_tri_diagonal.IsTriDiagonal() == false,
        "matrix named not_tri_diagonal is tri_diagonal"
      );

      // check if CreatePaddedMatrix works
      linal::Matrix< double> padded_matrix( symmetric.CreatePaddedMatrix( 2, 3, 1.0));

      for( size_t row( 0); row < symmetric.GetNumberRows(); ++row)
      {
        for( size_t column( 0); column < symmetric.GetNumberCols(); ++column)
        {
          BCL_Example_Check
          (
            symmetric( row, column) == padded_matrix( row, column),
            "elements don't match"
          );
        }
      }

      for( size_t row( symmetric.GetNumberRows()); row < padded_matrix.GetNumberRows(); ++row)
      {
        for( size_t column( symmetric.GetNumberCols()); column < padded_matrix.GetNumberCols(); ++column)
        {
          BCL_Example_Check
          (
            padded_matrix( row, column) == 1.0,
            "element should be 1.0"
          );
        }
      }

      BCL_Example_Check
      (
        padded_matrix.GetNumberRows() == symmetric.GetNumberRows() + 2 &&
        padded_matrix.GetNumberCols() == symmetric.GetNumberCols() + 3,
        "padded_matrix has wrong number od rows or columns"
      );

      // check if CreateSubMatrix works
      linal::Matrix< double> sub_matrix( padded_matrix.CreateSubMatrix( 3, 3, 0, 0));

      for( size_t row( 0); row < sub_matrix.GetNumberRows(); ++row)
      {
        for( size_t column( 0); column < sub_matrix.GetNumberCols(); ++column)
        {
          BCL_Example_Check
          (
            symmetric( row, column) == padded_matrix( row, column),
            "elements don't match"
          );
        }
      }

      BCL_Example_Check
      (
        sub_matrix.GetNumberRows() == 3 && sub_matrix.GetNumberCols() == 3,
        "sub_matrix has wrong number od rows or columns"
      );

      // check if SetSubMatrix works
      linal::Matrix< double> super_matrix( 5, 5, 0.0);
      super_matrix.SetSubMatrix( symmetric, 0, 0);

      for( size_t row( 0); row < symmetric.GetNumberRows(); ++row)
      {
        for( size_t column( 0); column < symmetric.GetNumberCols(); ++column)
        {
          BCL_Example_Check
          (
            symmetric( row, column) == super_matrix( row, column),
            "elements don't match"
          );
        }
      }

      for( size_t row( symmetric.GetNumberRows()); row < super_matrix.GetNumberRows(); ++row)
      {
        for( size_t column( symmetric.GetNumberCols()); column < super_matrix.GetNumberCols(); ++column)
        {
          BCL_Example_Check
          (
            super_matrix( row, column) == 0.0,
            "element should be 0.0"
          );
        }
      }

      // check if SetRow works
      linal::Matrix< double> test_matrix( 3, 4, 0.0);
      test_matrix( 0, 0) = 1.0;
      test_matrix( 0, 1) = 2.0;
      test_matrix( 0, 2) = 3.0;
      test_matrix( 0, 3) = 4.0;

      linal::Matrix< double> matrix_check( 3, 4, 0.0);

      linal::VectorND< double, 4> row_vector;
      row_vector( 0) = 1.0;
      row_vector( 1) = 2.0;
      row_vector( 2) = 3.0;
      row_vector( 3) = 4.0;

      matrix_check.ReplaceRow( 0, row_vector);

      for( size_t row( 0); row < matrix_check.GetNumberRows(); ++row)
      {
        for( size_t column( 0); column < matrix_check.GetNumberCols(); ++column)
        {
          BCL_Example_Check
          (
            matrix_check( row, column) == test_matrix( row, column),
            "elements don't match "
          );
        }
      }

      BCL_MessageStd( util::Format()( test_matrix));
      BCL_MessageStd( util::Format()( matrix_check));

      // Test ShrinkRows
      matrix_check.ShrinkRows( 2);

      BCL_Example_Check
      (
        matrix_check.GetNumberRows() == 2,
        "Shrink Rows did not properly execute"
      );

      // Test SetZero
      matrix_check.SetZero();

      for( size_t row( 0); row < matrix_check.GetNumberRows(); ++row)
      {
        for( size_t column( 0); column < matrix_check.GetNumberCols(); ++column)
        {
          BCL_Example_Check
          (
            matrix_check( row, column) == 0,
            "Set Zero did not properly execute"
          );
        }
      }

    /////////////////
    //// operators //
    /////////////////

      linal::Matrix< double> operator_matrix( 3, 4, 0.0);
      operator_matrix( 0, 0) = 1.0;
      operator_matrix( 0, 1) = 2.0;
      operator_matrix( 0, 2) = 3.0;
      operator_matrix( 0, 3) = 4.0;

      operator_matrix( 1, 0) = 5.0;
      operator_matrix( 1, 1) = 6.0;
      operator_matrix( 1, 2) = 7.0;
      operator_matrix( 1, 3) = 8.0;

      operator_matrix( 2, 0) =  9.0;
      operator_matrix( 2, 1) = 10.0;
      operator_matrix( 2, 2) = 11.0;
      operator_matrix( 2, 3) = 12.0;

      // access operator()
      BCL_Example_Check
      (
        operator_matrix( 1, 1) == 6.0,
        "Data Access operator did not properly execute"
      );

      // reference operator()
      const linal::Matrix< double> &const_operator_matrix( operator_matrix);

      BCL_Example_Check
      (
        const_operator_matrix( 1, 1) == 6.0,
        "Const Data Access operator did not properly execute"
      );

      // row operator()
      BCL_Example_Check
      (
        *operator_matrix[ 1] == 5.0,
        "Row operator did not point to the first member of the row " + util::Format()( *operator_matrix[ 1])
      );

      BCL_Example_Check
      (
        *const_operator_matrix[ 1] == 5.0,
        "const Row operator did not point to the first member of the row"
      );

      // Check assignment operator
      linal::Matrix< double> operator_matrix2 = operator_matrix;

      for( size_t row( 0); row < operator_matrix2.GetNumberRows(); ++row)
      {
        for( size_t column( 0); column < operator_matrix2.GetNumberCols(); ++column)
        {
          BCL_Example_Check
          (
            operator_matrix2( row, column) == operator_matrix( row, column),
            "elements don't match"
          );
        }
      }

      // Assignment from value operator()
      operator_matrix2 = 5.0;

      for( size_t row( 0); row < operator_matrix2.GetNumberRows(); ++row)
      {
        for( size_t column( 0); column < operator_matrix2.GetNumberCols(); ++column)
        {
          BCL_Example_Check
          (
            operator_matrix2( row, column) == 5.0,
            "Assigment from value operator does not work"
          );
        }
      }

    ////////////////////////
    //// input and output //
    ////////////////////////

      // filename
      BCL_MessageStd( "writing and reading from file");
      WriteBCLObject( matrix_copy);
      linal::Matrix< double> read_matrix;
      ReadBCLObject( read_matrix);
      BCL_Example_Check
      (
        read_matrix == matrix_copy,
        "written and read matrix are not identical" + util::Format()( matrix_copy) + " " + util::Format()( read_matrix)
      );

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleLinalMatrix

  const ExampleClass::EnumType ExampleLinalMatrix::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalMatrix())
  );

} // namespace bcl
