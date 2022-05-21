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

#ifndef BCL_MATH_LINEAR_LEAST_SQUARES_H_
#define BCL_MATH_LINEAR_LEAST_SQUARES_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_math_function_interface.h"
#include "linal/bcl_linal_matrix_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LinearLeastSquares
    //! @brief determining the parameters that best fit of a set of independent data to dependent data
    //!        where the parameters are determined using linear least squares
    //!
    //! @see @link example_math_linear_least_squares.cpp @endlink
    //! @author alexanns
    //! @date May 25, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API LinearLeastSquares :
      public util::ObjectInterface
    {
    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LinearLeastSquares();

      //! @brief Clone is the virtual copy constructor
      LinearLeastSquares *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief SolutionAndChiSquared finds the solutions as well as the chi-squared
      //! chi-squared is the sum of residuals squared divided by the sum of of the y's squared
      static storage::Pair< linal::Vector< double>, double> SolutionAndChiSquared
      (
        const linal::MatrixConstInterface< double> &X,
        const linal::VectorConstInterface< double> &Y
      );

      //! @brief SolutionAndChiSquared finds the solutions as well as the chi-squared
      //! chi-squared is the sum of residuals squared divided by the sum of of the y's squared
      static storage::Pair< linal::Vector< float>, float> SolutionAndChiSquared
      (
        const linal::MatrixConstInterface< float> &X,
        const linal::VectorConstInterface< float> &Y
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class LinearLeastSquares

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_LINEAR_LEAST_SQUARES_H_

