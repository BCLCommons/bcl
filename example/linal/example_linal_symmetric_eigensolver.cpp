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
#include "linal/bcl_linal_symmetric_eigensolver.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_linal_symmetric_eigensolver.cpp
  //!
  //! @author mendenjl
  //! @date Jun 19, 2014
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleLinalSymmetricEigenSolver :
    public ExampleInterface
  {
  public:

    ExampleLinalSymmetricEigenSolver *Clone() const
    {
      return new ExampleLinalSymmetricEigenSolver( *this);
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
      // create some random symmetric matrices
      linal::Matrix< double> random_5_by_4_matrix( 5, 4, double( 0.0));
      linal::Matrix< double> random_big_matrix( 10, 15, double( 0.0));
      random_5_by_4_matrix.AsVector().SetRand( 0.0, 10.0);
      random_big_matrix.AsVector().SetRand( 0.0, 10.0);

      linal::Matrix< double> four_by_four_symmetric( random_5_by_4_matrix.Transposed() * random_5_by_4_matrix);
      linal::Matrix< double> big_symmetric( random_big_matrix.Transposed() * random_big_matrix);

      linal::SymmetricEigenSolver< double> saes;

      BCL_ExampleCheck( saes.ComputeEigenvaluesOnly( four_by_four_symmetric), true);
      BCL_ExampleCheckWithinAbsTolerance
      (
        saes.GetSortedEigenvalues(),
        linal::MakeVector< double>( 480.411, 32.0611, 27.574, 5.44574),
        0.001
      );
      BCL_ExampleCheck( saes.ComputeEigenvaluesAndVectors( four_by_four_symmetric), true);
      double expected_eigenvectors[] =
      {
         0.45203 ,  -0.382426, -0.658758, -0.464173,
         0.361992,  0.687982 , 0.254356 , -0.575279,
         0.532664,  0.400492 , -0.335971, 0.665581 ,
         0.61717 ,  -0.469082, 0.623269 , 0.102947
      };
      linal::Matrix< double> expected_eigenvectors_mat( 4, 4, expected_eigenvectors);
      BCL_ExampleCheckWithinAbsTolerance
      (
        saes.GetSortedEigenvectors(),
        expected_eigenvectors_mat,
        0.001
      );
      BCL_ExampleIndirectCheckWithinAbsTolerance
      (
        saes.GetSortedEigenvalues(),
        linal::MakeVector< double>( 480.411, 32.0611, 27.574, 5.44574),
        0.001,
        "Computing eigenvectors does not change eigenvalues"
      );

      // create an object to hold the eigenvectors of the relatively big matrix; last five columns are zero because
      // this matrix will only have a rank of 10
      double expected_eigenvalues_big[] =
      {
        4396.26, 319.23, 244.443, 143.77, 127.803, 99.1342, 66.7642, 43.5377, 20.5639, 15.1134, 0, 0, 0, 0, 0
      };

      saes.ComputeEigenvaluesAndVectors( big_symmetric);
      BCL_ExampleCheckWithinAbsTolerance
      (
        saes.GetSortedEigenvalues(),
        linal::Vector< double>( 15, expected_eigenvalues_big),
        0.01
      );

      // example check copied over from linal::MatrixOperations for the old EigenSystemSquareSymmetricMatrixTridiagonal
      // testing tridiagonal
      double data1[ 9] = { 1.0, 2.0, 3.0, 2.0, 2.0, 4.0, 3.0, 4.0, 3.0};
      linal::Matrix< double> sym_mat( 3, 3, data1);
      saes.ComputeEigenvaluesOnly( sym_mat);
      BCL_ExampleCheckWithinTolerance( saes.GetSortedEigenvalues()( 2), -1.73042, 0.001);

    /////////////////
    // data access //
    /////////////////

    ///////////////
    // operators //
    ///////////////

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

  }; //end ExampleLinalSymmetricEigenSolver

  const ExampleClass::EnumType ExampleLinalSymmetricEigenSolver::s_Instance
  (
    GetExamples().AddEnum( ExampleLinalSymmetricEigenSolver())
  );

} // namespace bcl
