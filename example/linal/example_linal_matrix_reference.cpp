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
#include "linal/bcl_linal_matrix_reference.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_matrix_reference.cpp
  //! @details tests linal::MatrixReference< t_DataType>
  //!
  //! @author mendenjl
  //! @date Nov 07, 2012
  //! @remarks status complete
  //! @remarks reviewed by woetzen
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalMatrixReference :
    public ExampleInterface
  {
  public:

    ExampleLinalMatrixReference *Clone() const
    {
      return new ExampleLinalMatrixReference( *this);
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
      const size_t rows( 3);
      const size_t cols( 4);

      // linal::Matrix< double> to use as original matrix to reference
      linal::Matrix< double> matrix_orig( rows, cols, double( 3));
      linal::Matrix< double> squ_matrix( 3, 3, double( 2));

      // diagonal matrix
      double diag_elements[] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
      linal::Matrix< double> diag_matrix( 3, 3, diag_elements);

      // tridiagonal matrix
      double tridiag_elements[] = { 1.0, 4.0, 0.0, 0.0, 3.0, 4.0, 1.0, 0.0, 0.0, 2.0, 3.0, 4.0, 0.0, 0.0, 1.0, 3.0};
      linal::Matrix< double> tridiag_matrix( 4, 4, tridiag_elements);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! default constructor
      linal::MatrixReference< double> matrix_a;

      //! constructing from rows, cols, pointer to data
      linal::MatrixReference< double> matrix_b( rows, cols, matrix_orig.Begin());

      //! constructing from MatrixInterface
      linal::MatrixReference< double> matrix_c( matrix_orig);

      //! constructing from MatrixInterface
      linal::MatrixReference< double> matrix_ref_sq( squ_matrix);
      linal::MatrixReference< double> matrix_ref_diag( diag_matrix);
      linal::MatrixReference< double> matrix_ref_tridiag( tridiag_matrix);

    /////////////////
    // data access //
    /////////////////

      //! checking GetNumberRows()
      BCL_ExampleCheck( matrix_b.GetNumberRows(), matrix_orig.GetNumberRows());

      //! checking GetNumberCols()
      BCL_ExampleCheck( matrix_b.GetNumberCols(), matrix_orig.GetNumberCols());

      //! checking GetNumberOfElements()
      BCL_ExampleCheck( matrix_b.GetNumberOfElements(), matrix_orig.GetNumberCols() * matrix_orig.GetNumberRows());

      //! checking construction
      BCL_ExampleCheck( std::equal( matrix_b.Begin(), matrix_b.End(), matrix_orig.Begin()), true);

      //! checking Begin() and End()
      BCL_ExampleCheck( matrix_c.Begin() + matrix_c.GetNumberOfElements(), matrix_c.End());

    ///////////////
    // operators //
    ///////////////

      //! checking operator()( row, col)
      BCL_ExampleCheck( matrix_ref_diag( 1, 2), 0.0);

      //! checking operator = ( MatrixInterface
      matrix_a.Reference( matrix_orig);
      BCL_ExampleCheck( matrix_a( 1, 3), double( 3));

    ////////////////
    // operations //
    ////////////////

      //! IsSquare()
      BCL_ExampleCheck( matrix_b.IsSquare(), false);
      BCL_ExampleCheck( matrix_ref_sq.IsSquare(), true);

      //! IsDiagonal()
      BCL_ExampleCheck( matrix_b.IsDiagonal(), false);
      BCL_ExampleCheck( matrix_ref_diag.IsDiagonal(), true);

      //! IsTriDiagonal()
      BCL_ExampleCheck( matrix_b.IsTriDiagonal(), false);
      BCL_ExampleCheck( matrix_ref_tridiag.IsTriDiagonal(), true);

    //////////////////////
    // input and output //
    //////////////////////

      WriteBCLObject( matrix_b);

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleLinalMatrixReference

  const ExampleClass::EnumType ExampleLinalMatrixReference::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalMatrixReference())
  );

} // namespace bcl
