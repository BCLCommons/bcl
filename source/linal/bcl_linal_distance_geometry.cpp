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
#include "linal/bcl_linal_distance_geometry.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_principal_component_analysis.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  /////////////////
  // data access //
  /////////////////

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief virtual operator taking an ARGUMENT and returning a t_ResultType object
    //! @param SYMMETRIC_MATRIX containing pairwise distances for target
    //! @return function value of the given argument
    storage::Vector< VectorND< double, 3> > DistanceGeometry::operator()
    (
      const storage::SymmetricMatrix< double> &SYMMETRIC_MATRIX
    ) const
    {
      // Set result dimension
      const size_t result_dimension( 3);

      // Get input matrix size
      const size_t matrix_size( SYMMETRIC_MATRIX.GetSize());

      // Create matrix to be filled as the metric matrix
      Matrix< double> metric_matrix( matrix_size, matrix_size, 0.0);

      // TODO: Use math::Sqr()

      // Make a metric matrix using M = [m_ij] of dimension n, where m_ij = 1/2( p_0i^2 + p_0j^2 - p_ij^2)
      for( size_t i( 0); i < matrix_size; ++i)
      {
        for( size_t j( 0); j < matrix_size; ++j)
        {
          double p_0i = SYMMETRIC_MATRIX( i, 0);
          double p_0j = SYMMETRIC_MATRIX( 0, j);
          double p_ij = SYMMETRIC_MATRIX( i, j);
          metric_matrix( i, j) = ( pow( p_0i, 2) + pow( p_0j, 2) - pow( p_ij, 2)) / 2.0;
        }
      }

      // Create matrix to be filled along the diagonal with eigen values and empty eigen vectors matrix
      Matrix< double> eigen_values_matrix( metric_matrix);
      Matrix< double> eigen_vectors_matrix( matrix_size, matrix_size, 0.0);

      // Calculate eigen values
      SymmetricEigenSolver< double> eigensolver( metric_matrix, true);

      // Place diagonal containing values in a vector
      Vector< double> sorted_eigen_values( eigensolver.GetSortedEigenvalues());

      // Get sorted eigen_vectors matrix
      Matrix< double> sorted_eigen_vectors_matrix( eigensolver.GetSortedEigenvectors());

      // Create vector for checked and zeroed eigenvalues beyond the targeted dimension
      Matrix< double> clean_eigenvalues( 1, result_dimension, 0.0);

      // Take top three eigenvectors and check that they are all positive or zero, if negative set to zero
      for( size_t index( 0); index < result_dimension; ++index)
      {
        if( sorted_eigen_values( index) > 0)
        {
          // Take square root of each eigenvalue
          clean_eigenvalues( 0, index) = sqrt( sorted_eigen_values( index));
        }
      };

      // Multiply eigenvector by corresponding eigenvalue
      Matrix< double> result_matrix( matrix_size, result_dimension, 0.0);
      for( size_t col( 0); col < result_dimension; ++col)
      {
        for( size_t row( 0); row < matrix_size; ++row)
        {
          result_matrix( row, col) = clean_eigenvalues( 0, col) * sorted_eigen_vectors_matrix( row, col);
        }
      }

      // Create vector of size N to take all coordinates
      storage::Vector< VectorND< double, 3> > result_vector( matrix_size);

      // Iterate through n matrix placing first three positions into vector of VectorND of size result_dimension (3)
      VectorND< double, 3> temp_vect_nd;
      for( size_t i( 0); i < matrix_size; ++i)
      {
        for( size_t index( 0); index < result_dimension; ++index)
        {
          temp_vect_nd( index) = result_matrix( i, index);
        }
        result_vector( i) = temp_vect_nd;
      }

// TODO: Take first three using constructor from pointer

      // Return vector with coordinates
      return result_vector;

      // TODO: Add check later on for size mismatch ??
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DistanceGeometry::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DistanceGeometry::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // write members
      return OSTREAM;
    }

  } // namespace linal
} // namespace bcl
