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

#ifndef BCL_LINAL_MATRIX3X3_HPP_
#define BCL_LINAL_MATRIX3X3_HPP_

// include the header of this class
#include "bcl_linal_matrix3x3.h"

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix_inversion_gauss_jordan.h"
#include "bcl_linal_vector_3d.h"
#include "bcl_linal_vector_3d_operations.h"
#include "math/bcl_math_rotation_matrix_3d.h"
#include "math/bcl_math_statistics.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_logger_interface.h"

// external includes - sorted alphabetically
#include <math.h>

namespace bcl
{
  namespace linal
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    template< typename t_DataType>
    Matrix3x3< t_DataType>::Matrix3x3()
    {
      std::fill_n( &m_00a, s_NumberElements, t_DataType( 0));
    }

    //! @brief construct from filler
    //! @param FILL_VALUE assign every element to that value
    template< typename t_DataType>
    Matrix3x3< t_DataType>::Matrix3x3( const t_DataType &FILL_VALUE)
    {
      // set all values to FILL_VALUE
      std::fill_n( &m_00a, s_NumberElements, FILL_VALUE);
    }

    //! @brief construct from pointer to data
    //! @param DATA pointer to field of data
    template< typename t_DataType>
    Matrix3x3< t_DataType>::Matrix3x3( const t_DataType *DATA)
    {
      // copy data
      std::copy( DATA, DATA + s_NumberElements, &m_00a);
    }

    //! @brief copy constructor
    template< typename t_DataType>
    Matrix3x3< t_DataType>::Matrix3x3( const Matrix3x3< t_DataType> &MATRIX)
    {
      std::copy( &MATRIX.m_00a, &MATRIX.m_00a + s_NumberElements, &m_00a);
    }

    //! @brief copy constructor from MatrixInterface
    //! @param MATRIX_INTERFACE matrix interface to be copied from
    template< typename t_DataType>
    Matrix3x3< t_DataType>::Matrix3x3( const MatrixConstInterface< t_DataType> &MATRIX_INTERFACE)
    {
      // check if dimensions match
      BCL_Assert
      (
        MATRIX_INTERFACE.GetNumberRows() == s_Dimension || MATRIX_INTERFACE.GetNumberCols() == s_Dimension,
        "can only copy from 3x3 matrix!"
      );

      std::copy( MATRIX_INTERFACE.Begin(), MATRIX_INTERFACE.End(), &m_00a);
    }

    //! @brief Clone function
    //! @return pointer to new Matrix< t_DataType>
    template< typename t_DataType>
    Matrix3x3< t_DataType> *Matrix3x3< t_DataType>::Clone() const
    {
      return new Matrix3x3< t_DataType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    template< typename t_DataType>
    const std::string &Matrix3x3< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get number of rows
    //! @return number of rows
    template< typename t_DataType>
    size_t Matrix3x3< t_DataType>::GetNumberRows() const
    {
      return s_Dimension;
    }

    //! @brief get number of columns
    //! @return number of columns
    template< typename t_DataType>
    size_t Matrix3x3< t_DataType>::GetNumberCols() const
    {
      return s_Dimension;
    }

    //! @brief number of elements
    //! @return total number of elements in matrix
    template< typename t_DataType>
    size_t Matrix3x3< t_DataType>::GetNumberOfElements() const
    {
      return s_NumberElements;
    }

    //! @brief pointer to First Element
    //! @return const pointer to first element in range containing all elements of Matrix
    template< typename t_DataType>
    const t_DataType *Matrix3x3< t_DataType>::Begin() const
    {
      return &m_00a;
    }

    //! @brief pointer to First Element
    //! @return pointer to first element in range containing all elements of Matrix
    template< typename t_DataType>
    t_DataType *Matrix3x3< t_DataType>::Begin()
    {
      return &m_00a;
    }

    //! @brief pointer to end of range
    //! @return const pointer to address one after last element in Matrix
    template< typename t_DataType>
    const t_DataType *Matrix3x3< t_DataType>::End() const
    {
      return &m_00a + s_NumberElements;
    }

    //! @brief pointer to end of range
    //! @return pointer to address one after last element in Matrix
    template< typename t_DataType>
    t_DataType *Matrix3x3< t_DataType>::End()
    {
      return &m_00a + s_NumberElements;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief is matrix a square matrix
    //! @return true if number of cols and rows are identical
    template< typename t_DataType>
    bool Matrix3x3< t_DataType>::IsSquare() const
    {
      return true;
    }

    //! @brief is matrix a diagonal matrix
    //! @return true if all but the elements in the diagonal are 0
    template< typename t_DataType>
    bool Matrix3x3< t_DataType>::IsDiagonal() const
    {
      // return true if all elements but the diagonal is 0
      const t_DataType zero( 0);
      return math::EqualWithinMachineTolerance( m_01b, zero)
             && math::EqualWithinMachineTolerance( m_02c, zero)
             && math::EqualWithinMachineTolerance( m_10d, zero)
             && math::EqualWithinMachineTolerance( m_12f, zero)
             && math::EqualWithinMachineTolerance( m_20g, zero)
             && math::EqualWithinMachineTolerance( m_21h, zero);
    }

    //! @brief checks if matrix is defined
    //! @return bool - true if matrix is defined
    template< typename t_DataType>
    bool Matrix3x3< t_DataType>::IsDefined() const
    {
      return math::Statistics::IsDefined( &m_00a, &m_00a + s_NumberElements);
    }

    //! check whether matrix is symmetric
    template< typename t_DataType>
    bool Matrix3x3< t_DataType>::IsSymmetric() const
    {
      return math::EqualWithinMachineTolerance( m_01b, m_10d)
             && math::EqualWithinMachineTolerance( m_02c, m_20g)
             && math::EqualWithinMachineTolerance( m_12f, m_21h);
    }

    //! @brief is matrix a tridiagonal matrix
    //! @return if diagonal and adjacent diagonals are filled and the rest is 0
    template< typename t_DataType>
    bool Matrix3x3< t_DataType>::IsTriDiagonal() const
    {
      // return true if all elements but the inner three diagonals are 0
      return math::EqualWithinMachineTolerance( m_02c, t_DataType( 0))
             && math::EqualWithinMachineTolerance( m_20g, t_DataType( 0));
    }

    //! @brief swap elements between two rows
    //! @param ROW_A the first row
    //! @param ROW_B the second row
    template< typename t_DataType>
    void Matrix3x3< t_DataType>::SwapRows( const size_t ROW_A, const size_t ROW_B)
    {
      if( ROW_A >= s_Dimension || ROW_B >= s_Dimension)
      {
        return;
      }

      t_DataType tmp[ s_Dimension];
      t_DataType *row_a_beg( &m_00a + ROW_A * s_Dimension);
      t_DataType *row_a_end( row_a_beg + s_Dimension);
      t_DataType *row_b_beg( &m_00a + ROW_B * s_Dimension);
      t_DataType *row_b_end( row_b_beg + s_Dimension);

      std::copy( row_a_beg, row_a_end, tmp);
      std::copy( row_b_beg, row_b_end, row_a_beg);
      std::copy( tmp, tmp + s_Dimension, row_b_beg);
    }

    //! @brief replace row
    //! @param ROW the row to replace
    //! @param VECTOR the vector to replace the current row with
    template< typename t_DataType>
    void Matrix3x3< t_DataType>::ReplaceRow( const size_t ROW, const VectorConstInterface< t_DataType> &VECTOR)
    {
      if( ROW >= s_Dimension || VECTOR.GetSize() != s_Dimension)
      {
        return;
      }

      std::copy( VECTOR.Begin(), VECTOR.End(), &m_00a + ROW * s_Dimension);
    }

    //! @brief sort rows and given vectors (less than)
    //! @param VECTOR vector with 3 elements
    template< typename t_DataType>
    void Matrix3x3< t_DataType>::SortRowsAndVector( VectorInterface< t_DataType> &VECTOR)
    {
      // sort
      if( VECTOR( 0) < VECTOR( 1))
      {
        std::swap( VECTOR( 0), VECTOR( 1));
        SwapRows( 0, 1);
      }

      if( VECTOR( 1) < VECTOR( 2))
      {
        std::swap( VECTOR( 1), VECTOR( 2));
        SwapRows( 1, 2);
      }

      if( VECTOR( 0) < VECTOR( 1))
      {
        std::swap( VECTOR( 0), VECTOR( 1));
        SwapRows( 0, 1);
      }
    }

    //! @brief determinant of this matrix
    //! @return the determinant of the matrix
    template< typename t_DataType>
    t_DataType Matrix3x3< t_DataType>::Determinant() const
    {
      return
          m_00a * m_11e * m_22i
        + m_01b * m_12f * m_20g
        + m_02c * m_10d * m_21h
        - m_02c * m_11e * m_20g
        - m_01b * m_10d * m_22i
        - m_00a * m_12f * m_21h;
    }

    //! @brief trace of matrix
    //! @return the sum of the diagonal elements
    template< typename t_DataType>
    t_DataType Matrix3x3< t_DataType>::Trace() const
    {
      return m_00a + m_11e + m_22i;
    }

    //! @brief transpose the matrix
    template< typename t_DataType>
    void Matrix3x3< t_DataType>::Transpose()
    {
      std::swap( m_01b, m_10d);
      std::swap( m_02c, m_20g);
      std::swap( m_12f, m_21h);
    }

    //! @brief eigenvalues
    //! @return the eigenvalues
    template< typename t_DataType>
    Vector< t_DataType> Matrix3x3< t_DataType>::EigenValues() const
    {
      // simpler algorithm for symmetric matrix
      if( IsSymmetric())
      {
        return EigenValuesSymmetric();
      }

      Vector< t_DataType> eigenvalues( 3, t_DataType( 0));

      // matrix a b c \n d e f \n g h i
      // characteristic polynomial
      // a = a + e + i
      // b =
      //Use the equation from above to get your cubic equation and combine all constant terms possible to
      //give you a reduced equation we will use a, b, c and d to denote the coefficients of this equation.
      //Eqn = a*lambda^3 + b*lambda^2 + c*lambda + d = 0
      const t_DataType a( -1.0);
      const t_DataType b( m_00a + m_11e + m_22i); // diagonal
      const t_DataType c( m_01b * m_10d + m_02c * m_20g +
                          m_12f * m_21h - m_00a * m_11e -
                          m_00a * m_22i - m_11e * m_22i);
      const t_DataType d( Determinant());

      const t_DataType x( c / a - math::Sqr< t_DataType>( b) / ( 3 * math::Sqr< t_DataType>( a)));
      const t_DataType y( ( 2 * math::Pow< t_DataType>( b, 3) / math::Pow< t_DataType>( a, 3) - 9 * b * c / math::Sqr< t_DataType>( a) + 27 * d / a) / 27);
      const t_DataType z( math::Sqr< t_DataType>( y) / 4 + math::Pow< t_DataType>( x, 3) / 27);

      const t_DataType i( math::Sqrt< t_DataType>( math::Sqr< t_DataType>( y) / 4 - z));
      const t_DataType j( math::Pow< t_DataType>( i, 1.0 / 3.0));
      const t_DataType k( acos( -y / ( 2 * i)));
      const t_DataType m( cos( k / 3));
      const t_DataType n( math::Sqrt< t_DataType>( 3.0) * sin( k / 3));
      const t_DataType p( -b / ( 3 * a));

      eigenvalues( 0) = 2 * j * m + p;
      eigenvalues( 1) = -j * ( m + n) + p;
      eigenvalues( 2) = -j * ( m - n) + p;

      // end
      return eigenvalues;
    }

    //! @brief eigenvalues for a symmetric matrix
    //! @return the eigenvalues
    template< typename t_DataType>
    Vector< t_DataType> Matrix3x3< t_DataType>::EigenValuesSymmetric() const
    {
      Vector< t_DataType> eigenvalues( 3, t_DataType( 0));

      // Determine coefficients of characteristic poynomial. We write
      //       | a   d   f  |
      //  A =  | d*  b   e  |
      //       | f*  e*  c  |
      const t_DataType de( m_01b * m_12f);                  // d * e
      const t_DataType dd( math::Sqr< t_DataType>( m_01b)); // d^2
      const t_DataType ee( math::Sqr< t_DataType>( m_12f)); // e^2
      const t_DataType ff( math::Sqr< t_DataType>( m_02c)); // f^2
      const t_DataType m( Trace()); // sum of diagonal = trace
      const t_DataType c1( ( m_00a * m_11e + m_00a * m_22i + m_11e * m_22i) - dd - ee - ff); // a*b + a*c + b*c - d^2 - e^2 - f^2
      const t_DataType c0(   m_22i * dd    + m_00a * ee    + m_11e * ff     - m_00a * m_11e * m_22i - 2 * m_02c * de); // c*d^2 + a*e^2 + b*f^2 - a*b*c - 2*f*d*e)

      const t_DataType p( math::Sqr< t_DataType>( m) - 3 * c1);
      const t_DataType q( m * ( p - 1.5 * c1) - 13.5 * c0);
      const t_DataType sqrt_p( math::Sqrt< t_DataType>( math::Absolute( p)));

      t_DataType phi( 27.0 * ( 0.25 * math::Sqr< t_DataType>( c1) * ( p - c1) + c0 * ( q + 6.75 * c0)));
      phi = ( 1.0 / 3.0) * atan2( math::Sqrt< t_DataType>( math::Absolute( phi)), q);

      const t_DataType c( sqrt_p * cos( phi));
      const t_DataType s( ( 1.0 / math::Sqrt< t_DataType>( 3.0)) * sqrt_p * sin( phi));

      const double tmp( ( 1.0 / 3.0) * ( m - c));
      eigenvalues( 0) = tmp + c;
      eigenvalues( 1) = tmp - s;
      eigenvalues( 2) = tmp + s;

      // end
      return eigenvalues;
    }

    //! @brief transform matrix to tridiagonal form
    // Reduces a symmetric 3x3 matrix to tridiagonal form by applying
    // (unitary) Householder transformations:
    //            [ d[0]  e[0]       ]
    //    A = Q . [ e[0]  d[1]  e[1] ] . Q^T
    //            [       e[1]  d[2] ]
    //! @param EIGEN_VECTORS storage for Q
    //! @param EIGENVALUES storage for diagonal elements d
    //! @param OFF_DIAGONAL storage for off diagonal elements
    template< typename t_DataType>
    void Matrix3x3< t_DataType>::TridiagonalizeHouseholder( Matrix3x3< t_DataType> &EIGEN_VECTORS, Vector< t_DataType> &EIGEN_VALUES, t_DataType OFF_DIAGONAL[ s_Dimension])
    {
      t_DataType u[ s_Dimension], q[ s_Dimension];
      t_DataType omega, f;
      t_DataType K, h, g;

      // Initialize EIGEN_VECTORS to the identity matrix
      for( size_t i( 0); i < s_Dimension; ++i)
      {
        EIGEN_VECTORS( i, i) = 1.0;
        for( size_t j( 0); j < i; ++j)
        {
          EIGEN_VECTORS( i, j) = EIGEN_VECTORS( j, i) = 0.0;
        }
      }

      // Bring first row and column to the desired form
      h = math::Sqr< t_DataType>( m_01b) + math::Sqr< t_DataType>( m_02c);
      if( m_01b > 0.0)
      {
        g = -math::Sqrt< t_DataType>( h);
      }
      else
      {
        g = math::Sqrt< t_DataType>( h);
      }
      OFF_DIAGONAL[ 0] = g;
      f     = g * m_01b;
      u[ 1] = m_01b - g;
      u[ 2] = m_02c;

      omega = h - f;
      if( omega > 0.0)
      {
        omega = 1.0 / omega;
        K     = 0.0;
        for( size_t i( 1); i < s_Dimension; ++i)
        {
          f      = operator ()( 1, i) * u[ 1] + operator ()( i, 2) * u[2];
          q[ i]  = omega * f;                  // p
          K     += u[ i] * f;                  // u* A u
        }
        K *= 0.5 * math::Sqr< t_DataType>( omega);

        for( size_t i( 1); i < s_Dimension; ++i)
        {
          q[ i] = q[ i] - K * u[ i];
        }

        EIGEN_VALUES( 0) = m_00a;
        EIGEN_VALUES( 1) = m_11e - 2 * q[ 1] * u[ 1];
        EIGEN_VALUES( 2) = m_22i - 2 * q[ 2] * u[ 2];

        // Store inverse Householder transformation in Q
        for( size_t j( 1); j < s_Dimension; ++j)
        {
          f = omega * u[ j];
          for( size_t i( 1); i < s_Dimension; ++i)
          {
            EIGEN_VECTORS( i, j) -= f * u[ i];
          }
        }

        // Calculate updated A[1][2] and store it in OFF_DIAGONAL[1]
        OFF_DIAGONAL[ 1] = m_12f - q[ 1] * u[ 2] - u[ 1] * q[ 2];
      }
      else
      {
        for( size_t i( 0); i < s_Dimension; ++i)
        {
          EIGEN_VALUES( i) = operator ()( i, i);
        }
        OFF_DIAGONAL[ 1] = m_12f;
      }
    }

    //! @brief eigenvectors for symmetric matrix using QL over Householder
    //! @param EIGEN_VECTORS reference to the matrix that will hold the eigenvectors
    //! @param EIGEN_VALUES reference to the vector with eigenvalues
    //! @return true if successful
    template< typename t_DataType>
    bool Matrix3x3< t_DataType>::EigenVectorsSymmetricQL( Matrix3x3< t_DataType> &EIGEN_VECTORS, Vector< t_DataType> &EIGEN_VALUES) const
    {
      t_DataType off_diagonal[ 3];   // The third element is used only as temporary workspace
      double g, r, p, f, b, s, c, t; // Intermediate storage

      // Transform A to real tridiagonal form by the Householder method
      Matrix3x3< t_DataType> copy( *this);
      copy.TridiagonalizeHouseholder( EIGEN_VECTORS, EIGEN_VALUES, off_diagonal);

      // Calculate eigensystem of the remaining real symmetric tridiagonal matrix
      // with the QL method
      //
      // Loop over all off-diagonal elements
      for( int l( 0); l < int( s_Dimension) - 1; ++l)
      {
        size_t num_iter( 0); // number iterations
        while( true)
        {
          // Check for convergence and exit iteration loop if off-diagonal
          // element e(l) is zero
          int m;
          for( m = l; m <= int( s_Dimension) - 2; ++m)
          {
            g = math::Absolute( EIGEN_VALUES( m)) + math::Absolute( EIGEN_VALUES( m + 1));
            if( math::Absolute( off_diagonal[ m]) + g == g)
            {
              break;
            }
          }
          if( m == l)
          {
            break;
          }

          // should converge within 30 iterations
          if( num_iter++ >= 30)
          {
            return false;
          }

          // Calculate g = d_m - k
          g = ( EIGEN_VALUES( l + 1) - EIGEN_VALUES( l)) / ( off_diagonal[ l] + off_diagonal[ l]);
          r = math::Sqrt< t_DataType>( math::Sqr< t_DataType>( g) + 1.0);
          if( g > 0.0)
          {
            g = EIGEN_VALUES( m) - EIGEN_VALUES( l) + off_diagonal[ l] / ( g + r);
          }
          else
          {
            g = EIGEN_VALUES( m) - EIGEN_VALUES( l) + off_diagonal[ l] / ( g - r);
          }

          s = c = 1.0;
          p = 0.0;
          for( int i( m - 1); i >= l; --i)
          {
            f = s * off_diagonal[ i];
            b = c * off_diagonal[ i];
            if( math::Absolute( f) > math::Absolute( g))
            {
              c      = g / f;
              r      = math::Sqrt< t_DataType>( math::Sqr< t_DataType>( c) + 1.0);
              off_diagonal[ i + 1] = f * r;
              c     *= ( s = 1.0 / r);
            }
            else
            {
              s      = f / g;
              r      = math::Sqrt< t_DataType>( math::Sqr< t_DataType>( s) + 1.0);
              off_diagonal[ i + 1] = g * r;
              s     *= ( c = 1.0 / r);
            }

            g = EIGEN_VALUES( i + 1) - p;
            r = ( EIGEN_VALUES( i) - g) * s + 2 * c * b;
            p = s * r;
            EIGEN_VALUES( i + 1) = g + p;
            g = c * r - b;

            // Form eigenvectors
            for( size_t k( 0); k < s_Dimension; ++k)
            {
              t = EIGEN_VECTORS( k, i + 1);
              EIGEN_VECTORS( k, i + 1) = s * EIGEN_VECTORS( k, i) + c * t;
              EIGEN_VECTORS( k, i    ) = c * EIGEN_VECTORS( k, i) - s * t;
            }
          }
          EIGEN_VALUES( l) -= p;
          off_diagonal[ l] = g;
          off_diagonal[ m] = 0.0;
        }
      }

      return true;
    }

    //! @brief eigenvectors of symmetric matrix
    //! @param EIGEN_VECTORS reference to the matrix that will hold the eigenvectors
    //! @param EIGEN_VALUES reference to the vector that will hold the eigenvalues
    //! @return true is successful
    template< typename t_DataType>
    bool Matrix3x3< t_DataType>::EigenVectorsSymmetric( Matrix3x3< t_DataType> &EIGEN_VECTORS, Vector< t_DataType> &EIGEN_VALUES) const
    {
      if( !IsSymmetric())
      {
        return false;
      }

      Matrix3x3< t_DataType> copy( *this);
      EIGEN_VALUES = copy.EigenValuesSymmetric();

      const t_DataType n0( math::Sqr< t_DataType>( copy.m_00a) + math::Sqr< t_DataType>( copy.m_01b) + math::Sqr< t_DataType>( copy.m_02c));
      const t_DataType n1( math::Sqr< t_DataType>( copy.m_01b) + math::Sqr< t_DataType>( copy.m_11e) + math::Sqr< t_DataType>( copy.m_12f));

      t_DataType u, norm;
      t_DataType t = math::Absolute( EIGEN_VALUES( 0));
      if( ( u = math::Absolute( EIGEN_VALUES( 1))) > t)
      {
        t = u;
      }
      if( ( u = math::Absolute( EIGEN_VALUES( 2))) > t)
      {
        t = u;
      }
      if( t < 1.0)
      {
        u = t;
      }
      else
      {
        u = math::Sqr< t_DataType>( t);
      }
      t_DataType error( 256 * std::numeric_limits< t_DataType>::epsilon() * ( n0 + u) * ( n1 + u));

      EIGEN_VECTORS.m_01b = copy.m_01b * copy.m_12f - copy.m_02c * copy.m_11e;
      EIGEN_VECTORS.m_11e = copy.m_02c * copy.m_01b - copy.m_12f * copy.m_00a;
      EIGEN_VECTORS.m_21h = math::Sqr< t_DataType>( copy.m_01b);

      // Calculate first eigenvector by the formula
      //   v[0] = (A - EIGEN_VALUES( 0)).e1 x (A - EIGEN_VALUES( 0)).e2
      EIGEN_VECTORS.m_00a = EIGEN_VECTORS.m_01b + copy.m_02c * EIGEN_VALUES( 0);
      EIGEN_VECTORS.m_10d = EIGEN_VECTORS.m_11e + copy.m_12f * EIGEN_VALUES( 0);
      EIGEN_VECTORS.m_20g = ( copy.m_00a - EIGEN_VALUES( 0)) * ( copy.m_11e - EIGEN_VALUES( 0)) - EIGEN_VECTORS.m_21h;
      norm = math::Sqr< t_DataType>( EIGEN_VECTORS.m_00a) + math::Sqr< t_DataType>( EIGEN_VECTORS.m_10d) + math::Sqr< t_DataType>( EIGEN_VECTORS.m_20g);

      // If vectors are nearly linearly dependent, or if there might have
      // been large cancellations in the calculation of A[i][i] - EIGEN_VALUES( 0), fall
      // back to QL algorithm
      // Note that this simultaneously ensures that multiple eigenvalues do
      // not cause problems: If EIGEN_VALUES( 0) = EIGEN_VALUES( 1), then A - EIGEN_VALUES( 0) * I has rank 1,
      // i.e. all columns of A - EIGEN_VALUES( 0) * I are linearly dependent.
      if( norm <= error)
      {
        return EigenVectorsSymmetricQL( EIGEN_VECTORS, EIGEN_VALUES);
      }

      norm = math::Sqrt< t_DataType>( 1.0 / norm);
      for( size_t j( 0); j < s_Dimension; ++j)
      {
        EIGEN_VECTORS( j, 0) *= norm;
      }

      // calculate second eigenvector by the formula
      //   v[1] = (A - EIGEN_VALUES( 1)).e1 x (A - EIGEN_VALUES( 1)).e2
      EIGEN_VECTORS.m_01b = EIGEN_VECTORS.m_01b + copy.m_02c * EIGEN_VALUES( 1);
      EIGEN_VECTORS.m_11e = EIGEN_VECTORS.m_11e + copy.m_12f * EIGEN_VALUES( 1);
      EIGEN_VECTORS.m_21h = copy.m_00a * copy.m_11e - EIGEN_VECTORS.m_21h;
      norm = math::Sqr< t_DataType>( EIGEN_VECTORS.m_01b) + math::Sqr< t_DataType>( EIGEN_VECTORS.m_11e) + math::Sqr< t_DataType>( EIGEN_VECTORS.m_21h);

      error = n0 * n1;
      if( norm <= error)
      {
        return EigenVectorsSymmetricQL( EIGEN_VECTORS, EIGEN_VALUES);
      }

      norm = math::Sqrt< t_DataType>( 1.0 / norm);
      for( size_t j( 0); j < s_Dimension; ++j)
      {
        EIGEN_VECTORS( j, 1) *= norm;
      }

      // Calculate third eigenvector according to
      //   v[2] = v[0] x v[1]
      EIGEN_VECTORS.m_02c = EIGEN_VECTORS.m_10d * EIGEN_VECTORS.m_21h - EIGEN_VECTORS.m_20g * EIGEN_VECTORS.m_11e;
      EIGEN_VECTORS.m_12f = EIGEN_VECTORS.m_20g * EIGEN_VECTORS.m_01b - EIGEN_VECTORS.m_00a * EIGEN_VECTORS.m_21h;
      EIGEN_VECTORS.m_22i = EIGEN_VECTORS.m_00a * EIGEN_VECTORS.m_11e - EIGEN_VECTORS.m_10d * EIGEN_VECTORS.m_01b;

      // end
      return true;
    }

    //! @brief computes the euler angles between the two given frames.
    //! @param FRAME_1 the reference frame for the computation
    //! @param FRAME_2 the rotated frame for the computation
    //! @return the euler angles between the two given frames in the order precession, nutation, rotation
    template< typename t_DataType>
    Vector3D Matrix3x3< t_DataType>::ComputeEulerAngles
    (
      const Matrix3x3< t_DataType> &FRAME_1,
      const Matrix3x3< t_DataType> &FRAME_2
    )
    {
      // initialize projection matrix for coordinate system conversion
      Matrix3x3< t_DataType> proj_matrix( util::GetUndefined< t_DataType>());
      for( size_t x( 0); x <= 2; ++x)
      {
        for( size_t y( 0); y <= 2; ++y)
        {
          const Vector3D row_1( FRAME_1( x, 0), FRAME_1( x, 1), FRAME_1( x, 2));
          const Vector3D row_2( FRAME_2( y, 0), FRAME_2( y, 1), FRAME_2( y, 2));
          proj_matrix( x, y) = ScalarProduct( row_1, row_2);
        }
      }

      // compute the Euler angles in the order precession, nutation, rotation
      Vector3D euler_angles( util::GetUndefined< t_DataType>());
      const t_DataType z1xy( std::sqrt( std::pow( proj_matrix( 2, 0), 2) + std::pow( proj_matrix( 2, 1), 2)));
      if( z1xy > std::numeric_limits< t_DataType>::epsilon())
      {
        euler_angles( 0) = std::atan2
            (
              proj_matrix( 1, 0) * proj_matrix( 2, 1) - proj_matrix( 1, 1) * proj_matrix( 2, 0),
              proj_matrix( 0, 1) * proj_matrix( 2, 1) - proj_matrix( 0, 1) * proj_matrix( 2, 0)
            );
        euler_angles( 1) = std::atan2( z1xy, proj_matrix( 2, 2));
        euler_angles( 2) = -std::atan2( -proj_matrix( 2, 0), proj_matrix( 2, 1));
      }
      else
      {
        euler_angles( 0) = 0.0;
        euler_angles( 1) = proj_matrix( 2, 2) > 0.0 ? 0.0 : math::g_Pi;
        euler_angles( 2) = -std::atan2( proj_matrix( 0, 1), proj_matrix( 0, 0));
      }

      return euler_angles;
    }

    //! @brief computes the Euler angles according to x-y-z convention between the two given frames.
    //! @param FRAME_1 the reference frame for the computation
    //! @param FRAME_2 the rotated frame for the computation
    //! @return the euler angles between the two given frames according to x-y-z convention
    template< typename t_DataType>
    Vector3D Matrix3x3< t_DataType>::ComputeEulerAnglesXYZ
    (
      const Matrix3x3< double> &FRAME_1,
      const Matrix3x3< double> &FRAME_2
    )
    {
      MatrixInversionGaussJordan< double> mat_1( FRAME_1);
      const Matrix< double> inv_mat_1( mat_1.ComputeInverse());
      const Matrix3x3< double> rotation_matrix( FRAME_2 * inv_mat_1);
      const math::RotationMatrix3D tmp_mat( rotation_matrix);
      const Vector3D ea( tmp_mat.EulerAnglesXYZ());

      return ea;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief return reference to changeable element ( ROW, COL)
    //! @param ROW the row number, starting with 0
    //! @param COL the col number, starting with 0
    //! @return changeable reference to the element defined bey ROW and COL number
    template< typename t_DataType>
    t_DataType &Matrix3x3< t_DataType>::operator()( const size_t ROW, const size_t COL)
    {
      return ( &m_00a)[ ROW * s_Dimension + COL];
    }

    //! @brief return reference to const  element ( ROW, COL)
    //! @param ROW the row number, starting with 0
    //! @param COL the col number, starting with 0
    //! @return const element defined bey ROW and COL number
    template< typename t_DataType>
    const t_DataType &Matrix3x3< t_DataType>::operator()( const size_t ROW, const size_t COL) const
    {
      return ( &m_00a)[ ROW * s_Dimension + COL];
    }

    //! @brief access to a particular row
    //! @param ROW the row to access
    //! @return a pointer to the first member of that row
    template< typename t_DataType>
    t_DataType *Matrix3x3< t_DataType>::operator[]( const size_t ROW)
    {
      return &m_00a + ROW * s_Dimension;
    }

    //! @brief access to a particular row
    //! @param ROW the row to access
    //! @return a constant pointer to the first member of that row
    template< typename t_DataType>
    const t_DataType *Matrix3x3< t_DataType>::operator[]( const size_t ROW) const
    {
      return &m_00a + ROW * s_Dimension;
    }

    //! @brief assignment from Matrix3x3
    //! @param MATRIX the matrix used as source
    //! @return reference to this Matrix
    template< typename t_DataType>
    Matrix3x3< t_DataType> &Matrix3x3< t_DataType>::operator =( const Matrix3x3< t_DataType> &MATRIX)
    {
      // copy all elements
      if( &m_00a != &MATRIX.m_00a)
      {
        // copy elements
        std::copy( &MATRIX.m_00a, &MATRIX.m_00a + s_NumberElements, &m_00a);
      }

      // return reference to this Matrix
      return *this;
    }

    //! @brief assignment from MatrixInterface
    //! @param MATRIX_INTERFACE the matrix used as source
    //! @return reference to this Matrix
    template< typename t_DataType>
    Matrix3x3< t_DataType> &Matrix3x3< t_DataType>::operator =( const MatrixConstInterface< t_DataType> &MATRIX_INTERFACE)
    {
      // copy all elements
      if( &m_00a != MATRIX_INTERFACE.Begin())
      {
        // check if dimensions match
        BCL_Assert
        (
          MATRIX_INTERFACE.GetNumberRows() == s_Dimension || MATRIX_INTERFACE.GetNumberCols() == s_Dimension,
          "can only assign from 3x3 matrix!"
        );

        std::copy( MATRIX_INTERFACE.Begin(), MATRIX_INTERFACE.End(), &m_00a);
      }

      // return reference to this Matrix
      return *this;
    }

    //! @brief assignment from value
    //! @param VALUE all elements are set to that value
    //! @return reference to this assigned Matrix
    template< typename t_DataType>
    Matrix3x3< t_DataType> &Matrix3x3< t_DataType>::operator =( const t_DataType &VALUE)
    {
      // set all element to given VALUE
      std::fill( &m_00a, &m_00a + s_NumberElements, VALUE);

      // return reference to this Vector
      return *this;
    }

    //! @brief multiply with a second matrix3x3
    //! @param MATRIX other Matrix3x3
    //! @return reference to this matrix after multiplication
    template< typename t_DataType>
    Matrix3x3< t_DataType> &Matrix3x3< t_DataType>::operator *=( const Matrix3x3< t_DataType> &MATRIX)
    {
      Matrix3x3< t_DataType> copy( *this);

      m_00a = copy.m_00a * MATRIX.m_00a + copy.m_01b * MATRIX.m_10d + copy.m_02c * MATRIX.m_20g;
      m_01b = copy.m_00a * MATRIX.m_01b + copy.m_01b * MATRIX.m_11e + copy.m_02c * MATRIX.m_21h;
      m_02c = copy.m_00a * MATRIX.m_02c + copy.m_01b * MATRIX.m_12f + copy.m_02c * MATRIX.m_22i;

      m_10d = copy.m_10d * MATRIX.m_00a + copy.m_11e * MATRIX.m_10d + copy.m_12f * MATRIX.m_20g;
      m_11e = copy.m_10d * MATRIX.m_01b + copy.m_11e * MATRIX.m_11e + copy.m_12f * MATRIX.m_21h;
      m_12f = copy.m_10d * MATRIX.m_02c + copy.m_11e * MATRIX.m_12f + copy.m_12f * MATRIX.m_22i;

      m_20g = copy.m_20g * MATRIX.m_00a + copy.m_21h * MATRIX.m_10d + copy.m_22i * MATRIX.m_20g;
      m_21h = copy.m_20g * MATRIX.m_01b + copy.m_21h * MATRIX.m_11e + copy.m_22i * MATRIX.m_21h;
      m_22i = copy.m_20g * MATRIX.m_02c + copy.m_21h * MATRIX.m_12f + copy.m_22i * MATRIX.m_22i;

      return *this;
    }

    //! @brief multiply with a second matrix3x3
    //! @param MATRIX other Matrix3x3
    //! @return the result
    template< typename t_DataType>
    Matrix3x3< t_DataType> Matrix3x3< t_DataType>::operator *( const Matrix3x3< t_DataType> &MATRIX) const
    {
      Matrix3x3< t_DataType> copy( *this);

      copy.m_00a = m_00a * MATRIX.m_00a + m_01b * MATRIX.m_10d + m_02c * MATRIX.m_20g;
      copy.m_01b = m_00a * MATRIX.m_01b + m_01b * MATRIX.m_11e + m_02c * MATRIX.m_21h;
      copy.m_02c = m_00a * MATRIX.m_02c + m_01b * MATRIX.m_12f + m_02c * MATRIX.m_22i;

      copy.m_10d = m_10d * MATRIX.m_00a + m_11e * MATRIX.m_10d + m_12f * MATRIX.m_20g;
      copy.m_11e = m_10d * MATRIX.m_01b + m_11e * MATRIX.m_11e + m_12f * MATRIX.m_21h;
      copy.m_12f = m_10d * MATRIX.m_02c + m_11e * MATRIX.m_12f + m_12f * MATRIX.m_22i;

      copy.m_20g = m_20g * MATRIX.m_00a + m_21h * MATRIX.m_10d + m_22i * MATRIX.m_20g;
      copy.m_21h = m_20g * MATRIX.m_01b + m_21h * MATRIX.m_11e + m_22i * MATRIX.m_21h;
      copy.m_22i = m_20g * MATRIX.m_02c + m_21h * MATRIX.m_12f + m_22i * MATRIX.m_22i;

      return copy;
    }

    //! @brief subtract second matrix
    //! @param MATRIX the matrix to subtract
    //! @return reference to this matrix
    template< typename t_DataType>
    Matrix3x3< t_DataType> &Matrix3x3< t_DataType>::operator -=( const Matrix3x3< t_DataType> &MATRIX)
    {
      std::transform( &m_00a, &m_00a + s_NumberElements, &MATRIX.m_00a, &m_00a, std::minus< t_DataType>());
      return *this;
    }

    //! @brief subtract second matrix
    //! @param MATRIX the matrix to subtract
    //! @return result matrix
    template< typename t_DataType>
    Matrix3x3< t_DataType> Matrix3x3< t_DataType>::operator -( const Matrix3x3< t_DataType> &MATRIX) const
    {
      Matrix3x3< t_DataType> copy( *this);
      std::transform( &copy.m_00a, &copy.m_00a + s_NumberElements, &MATRIX.m_00a, &copy.m_00a, std::minus< t_DataType>());
      return copy;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_DataType>
    std::istream &Matrix3x3< t_DataType>::Read( std::istream &ISTREAM)
    {
      // read each element from ISTREAM
      for( t_DataType *ptr( Begin()), *ptr_end( End()); ptr != ptr_end; ++ptr)
      {
        BCL_Assert( io::Serialize::Read( *ptr, ISTREAM), "unable to read element");
      }

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    template< typename t_DataType>
    std::ostream &Matrix3x3< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      //write date
      for( size_t i = 0; i < s_Dimension;)
      {
        io::Serialize::InsertIndent( OSTREAM, INDENT);
        for( size_t j = 0; j < s_Dimension; ++j)
        {
          io::Serialize::Write( operator()( i, j), OSTREAM) << '\t';
        }
        ++i;
        // print line brake except after last line
        if( i != s_Dimension)
        {
          OSTREAM << '\n';
        }
      }

      //return
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! single instance of the Matrix class
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> Matrix3x3< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new Matrix3x3< t_DataType>())
    );

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_MATRIX3X3_HPP_
