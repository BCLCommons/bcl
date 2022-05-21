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

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "math/bcl_math_linear_least_squares.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_matrix_inversion_moore_penrose.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> LinearLeastSquares::s_Instance
    (
      GetObjectInstances().AddInstance( new LinearLeastSquares())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LinearLeastSquares::LinearLeastSquares()
    {
    }

    //! @brief Clone is the virtual copy constructor
    LinearLeastSquares *LinearLeastSquares::Clone() const
    {
      return new LinearLeastSquares( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &LinearLeastSquares::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief SolutionAndChiSquared finds the solutions as well as the chi-squared
    //! chi-squared is the sum of residuals squared divided by the sum of of the y's squared
    storage::Pair< linal::Vector< double>, double> LinearLeastSquares::SolutionAndChiSquared
    (
      const linal::MatrixConstInterface< double> &X,
      const linal::VectorConstInterface< double> &Y
    )
    {
      // make a reference to X or X_Transposed, whichever has more rows
      // this way, the user does not have to remember whether X needs to be tall or wide, the function works either way
      if( X.GetNumberCols() > X.GetNumberRows())
      {
        return SolutionAndChiSquared( X.Transposed(), Y);
      }

      // Get the solutions
      linal::MatrixInversionMoorePenrose< double> inverter( X);
      const linal::Vector< double> solutions( inverter.Solve( Y));

      // calculate the sum of the solutions squared
      const double sum_solutions_squared( std::inner_product( Y.Begin(), Y.End(), Y.Begin(), 0.0));

      // calculate the sum of the residuals.
      // a residual is the experimental value minus the predicted value squared
      // the predicted value is the row of the matrix times the solutions
      double sum_residuals( 0.0);
      for( size_t index( 0), y_size( Y.GetSize()); index < y_size; index++)
      {
        // calculate the value given by the linear-least squares fit to the data
        const double linear_estimate( std::inner_product( solutions.Begin(), solutions.End(), X[ index], 0.0));

        // add the difference between the estimate and the actual value, squared, to the residual sum
        sum_residuals += std::pow( Y( index) - linear_estimate, 2.0);
      }

      return storage::Pair< linal::Vector< double>, double>
             (
               solutions,
               sum_residuals / sum_solutions_squared
             );
    }

    //! @brief SolutionAndChiSquared finds the solutions as well as the chi-squared
    //! chi-squared is the sum of residuals squared divided by the sum of of the y's squared
    storage::Pair< linal::Vector< float>, float> LinearLeastSquares::SolutionAndChiSquared
    (
      const linal::MatrixConstInterface< float> &X,
      const linal::VectorConstInterface< float> &Y
    )
    {
      // make a reference to X or X_Transposed, whichever has more rows
      // this way, the user does not have to remember whether X needs to be tall or wide, the function works either way
      if( X.GetNumberCols() > X.GetNumberRows())
      {
        return SolutionAndChiSquared( X.Transposed(), Y);
      }

      // Get the solutions
      linal::MatrixInversionMoorePenrose< float> inverter( X);
      const linal::Vector< float> solutions( inverter.Solve( Y));

      // calculate the sum of the solutions squared
      const float sum_solutions_squared( std::inner_product( Y.Begin(), Y.End(), Y.Begin(), 0.0));

      // calculate the sum of the residuals.
      // a residual is the experimental value minus the predicted value squared
      // the predicted value is the row of the matrix times the solutions
      float sum_residuals( 0.0);
      for( size_t index( 0), y_size( Y.GetSize()); index < y_size; index++)
      {
        // calculate the value given by the linear-least squares fit to the data
        const float linear_estimate( std::inner_product( solutions.Begin(), solutions.End(), X[ index], 0.0));

        // add the difference between the estimate and the actual value, squared, to the residual sum
        sum_residuals += std::pow( Y( index) - linear_estimate, float( 2.0));
      }

      return storage::Pair< linal::Vector< float>, float>
             (
               solutions,
               sum_residuals / sum_solutions_squared
             );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LinearLeastSquares::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &LinearLeastSquares::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace math
} // namespace bcl

