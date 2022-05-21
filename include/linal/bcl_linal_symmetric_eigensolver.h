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

#ifndef BCL_LINAL_SYMMETRIC_EIGENSOLVER_H_
#define BCL_LINAL_SYMMETRIC_EIGENSOLVER_H_

// include the namespace header
#include "bcl_linal.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix.h"
#include "bcl_linal_matrix_operations.h"
#include "bcl_linal_matrix_reference.h"
#include "bcl_linal_operations.h"
#include "bcl_linal_operations_interface.h"
#include "bcl_linal_vector.h"
#include "bcl_linal_vector_const_reference.h"
#include "bcl_linal_vector_operations.h"
#include "math/bcl_math.h"
#include "math/bcl_math_rotation_matrix_2d.h"
#include "math/bcl_math_statistics.h"
#include "storage/bcl_storage_vector.h"
#include "type/bcl_type_enable_if.h"

// external includes - sorted alphabetically
#include <algorithm>
#include <numeric>

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//!
//! @file bcl_linal_symmetric_eigensolver.h
//!
//! Routines are adapted from Eigen,  lightweight C++ template library for linear algebra
//! @see @link http://eigen.tuxfamily.org/index.php @endlink
//! including the classes: SelfAdjointEigenSolver, HouseholderSequence, Tridiagonalization
//!
//! Eigen is subject to the following Copyright, which likewise applies to this file (only):
//! Copyright (C) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
//! Copyright (C) 2009 Gael Guennebaud <gael.guennebaud@inria.fr>
//!
//! This Source Code Form is subject to the terms of the Mozilla
//! Public License v. 2.0. If a copy of the MPL was not distributed
//! with this file, You can obtain one at http://mozilla.org/MPL/2.0/.
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace bcl
{
  namespace linal
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SymmetricEigenSolver
    //! @brief Computes eigenvalues and eigenvectors of symmetric, real matrices
    //!
    //! @tparam t_DataType the internal datatype of the matrix
    //!
    //! This class computes the eigenvalues and eigenvectors of a
    //! symmetric matrix. These are the scalars \f$ \lambda \f$ and vectors
    //! @f$ v \f$ such that \f$ Av = \lambda v \f$.  The eigenvalues of a
    //! selfadjoint matrix are always real. If \f$ D \f$ is a diagonal matrix with
    //! he eigenvalues on the diagonal, and \f$ V \f$ is a matrix with the
    //! eigenvectors as its columns, then \f$ A = V D V^{-1} \f$ (for selfadjoint
    //! matrices, the matrix \f$ V \f$ is always invertible). This is called the
    //! eigendecomposition.
    //!
    //! The algorithm exploits the fact that the matrix is symmetric, making it
    //! faster and more accurate than the general purpose eigenvalue algorithms
    //!
    //! Call the function Compute() to compute the eigenvalues and eigenvectors of
    //! a given matrix. Alternatively, you can use the
    //! SymmetricEigenSolver(const MatrixType&, bool) constructor which computes
    //! the eigenvalues and eigenvectors at construction time. Once the eigenvalue
    //! and eigenvectors are computed, they can be retrieved with the GetSortedEigenvalues()
    //! and GetSortedEigenvectors() functions.
    //!
    //! Differences between the BCL version of SelfAdjointEigenSolver's and Eigen's:
    //!   -- BCL version assumes row-major matrices, Eigen's version is optimized for column-major matrices
    //!   -- BCL version lacks support for complex types for simpler code and performance
    //!   -- Better use of object-oriented programming; less reliance on globally defined options
    //!   -- Somewhat slower because eigen has many optimizations using expression templates (not used in the BCL), that
    //!      require considerable care in implementation and maintenance (see Eigen topic aliasing); such as the need
    //!      to watch for matrices being on both the LHS and RHS of any expression.  In the BCL version; for loops are
    //!      occasionally used in place of Eigen's expressions to achieve similar performance
    //!   -- Eigen sorts eigenvalues/vectors in ascending order (smallest value first); the bcl sorts in
    //!      descending order to facilitate algorithms such as PCA where only the first N eigenvectors are desired,
    //!      since in this case the eigenmatrix may simply be truncated.
    //!   -- BCL uses STL's sorting algorithm O(nlog(n)) on the eigenvalues and vectors for somewhat better performance
    //!
    //! @see @link example_linal_symmetric_eigensolver.cpp @endlink
    //! @author mendenjl
    //! @date Jun 9, 2014
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    template< typename t_DataType>
    class SymmetricEigenSolver :
      public util::ObjectInterface
    {
    public:

    ///////////
    // enums //
    ///////////

      //! Status of the computation
      enum Status
      {
        e_Uninitialized,      //!< compute never called
        e_NoConvergence,      //!< Computation failed to converge
        e_Success,            //!< Successful computation
        s_NumberStatus        //!< Compute eigenvectors and sort them and the values
      };

    private:

      //! @brief Maximum number of iterations.
      //!
      //! The algorithm terminates if it does not converge within m_maxIterations * n iterations, where n
      //! denotes the size of the matrix. This value is currently set to 30 (copied from LAPACK).
      static const int s_MaxIterations = 30;

      Matrix< t_DataType> m_Eivec;          //!< Eigenvectors, if they were computed
      Matrix< t_DataType> m_Tmp;            //!< Temporary matrix used if eigenvalues are requested
      Vector< t_DataType> m_Eivalues;       //!< Eigenvalues
      Vector< t_DataType> m_Subdiag;        //!< Subdiagonal of the tridiagonal matrix
      Status              m_Info;           //!< Status of the computation
      bool                m_HaveVectors;    //!< True if eigenvectors were computed

    public:

      //  @brief Default constructor for fixed-size matrices.
      //
      //  The default constructor is useful in cases in which the user intends to
      // perform decompositions via compute().This constructor
      //  can only be used if \p _MatrixType is a fixed-size matrix; use
      //  SelfAdjointEigenSolver(size_t) for dynamic-size matrices.
      SymmetricEigenSolver() :
        m_Eivec(),
        m_Eivalues(),
        m_Subdiag(),
        m_Info( e_Uninitialized),
        m_HaveVectors( false)
      {
      }

      //! @brief Constructor, pre-allocates memory for dynamic-size matrices.
      //!
      //! @param [in]  SIZE  Positive integer, size of the matrix whose
      //! eigenvalues and eigenvectors will be computed.
      //!
      //! This constructor is useful for dynamic - size matrices, when the user
      //! intends to perform decompositions via compute(). The \p size
      //! parameter is only used as a hint. It is not an error to give a wrong
      //! @p size, but it may impair performance.
      //!
      //! @sa compute() for an example
      //!
      SymmetricEigenSolver( const size_t &SIZE) :
        m_Eivec( SIZE, SIZE),
        m_Tmp( m_Eivec),
        m_Eivalues( SIZE),
        m_Subdiag( SIZE > 1 ? SIZE - 1 : 1),
        m_Info( e_Uninitialized)
      {
      }

      //! @brief Constructor; computes eigendecomposition of given matrix.
      //!
      //! @param[in]  matrix  Selfadjoint matrix whose eigendecomposition is to
      //!    be computed. Only the lower triangular part of the matrix is referenced.
      //! @param[in]  EIGENVECTORS true to also compute eigenvectors
      //!
      //! This constructor calls compute(const MatrixType&, int) to compute the
      //! eigenvalues of the matrix \p matrix. The eigenvectors are computed if
      //! @p options equals #ComputeEigenvectors.
      //!
      //! Example: \include SelfAdjointEigenSolver_SelfAdjointEigenSolver_MatrixType.cpp
      //! Output: \verbinclude SelfAdjointEigenSolver_SelfAdjointEigenSolver_MatrixType.out
      //!
      //! @sa compute(const MatrixType&, int)
      //!
      SymmetricEigenSolver( const MatrixConstInterface< t_DataType> &MATRIX, const bool &EIGENVECTORS) :
        m_Eivec( MATRIX.GetNumberRows(), MATRIX.GetNumberCols()),
        m_Tmp( m_Eivec),
        m_Eivalues( MATRIX.GetNumberCols()),
        m_Subdiag( MATRIX.GetNumberRows() > 1 ? MATRIX.GetNumberRows() - 1 : 1),
        m_Info( e_Uninitialized),
        m_HaveVectors( true)
      {
        Compute( MATRIX, EIGENVECTORS);
      }

      //! @brief Clone function
      //! @return pointer to new MatrixInterface< t_DataType>
      SymmetricEigenSolver< t_DataType> *Clone() const
      {
        return new SymmetricEigenSolver< t_DataType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Returns the eigenvectors of given matrix.
      //!
      //! @returns  A const reference to the matrix whose columns are the eigenvectors.
      //!
      //! @pre The eigenvectors have been computed before.
      //!
      //! Row \f$ k \f$ of the returned matrix is an eigenvector corresponding
      //! to eigenvalue number \f$ k \f$ as returned by eigenvalues().  The
      //! eigenvectors are normalized to have (Euclidean) norm equal to one. If
      //! this object was used to solve the eigenproblem for the selfadjoint
      //! matrix \f$ A \f$, then the matrix returned by this function is the
      //! matrix \f$ V \f$ in the eigendecomposition \f$ A = V D V^{-1} \f$.
      //!
      //! @sa eigenvalues()
      //!
      const Matrix< t_DataType> &GetSortedEigenvectors() const
      {
        BCL_Assert( m_Info != e_Uninitialized, "Eigenvalues and vectors have not yet been computed!");
        BCL_Assert( m_Info != e_NoConvergence, "The computation failed to converge, no eigenvectors available!");
        BCL_Assert( m_HaveVectors, "The eigenvectors have not been computed together with the eigenvalues.");
        return m_Eivec;
      }

      //! @brief Returns the eigenvalues of given matrix, sorted (descending), so largest eigenvalues first
      //!
      //! @returns A const reference to the column vector containing the eigenvalues.
      //!
      //! @pre The eigenvalues have been computed before.
      //!
      //! The eigenvalues are repeated according to their algebraic multiplicity,
      //! so there are as many eigenvalues as rows in the matrix. The eigenvalues
      //! are sorted in desecending order.
      //!
      //! @sa GetSortedEigenvectors()
      //!
      const Vector< t_DataType> &GetSortedEigenvalues() const
      {
        BCL_Assert( m_Info != e_Uninitialized, "Eigenvalues have not yet been computed!");
        BCL_Assert( m_Info != e_NoConvergence, "The computation failed to converge, no eigenvalues available!");
        return m_Eivalues;
      }

      //! @brief compute the eigenvalues in O(4n^3/3) time
      //! @param MATRIX Selfadjoint matrix whose eigendecomposition is to be computed
      //!        Only the upper triangular part of the matrix is referenced.
      //! @return true on success (convergence)
      bool ComputeEigenvaluesOnly( const MatrixConstInterface< t_DataType> &MATRIX)
      {
        Compute( MATRIX, false);
        return m_Info == e_Success;
      }

      //! @brief compute the eigenvalues and vectors in O(9n^3) time
      //! @param MATRIX Selfadjoint matrix whose eigendecomposition is to be computed
      //!        Only the upper triangular part of the matrix is referenced.
      //! @return true on success (convergence)
      bool ComputeEigenvaluesAndVectors( const MatrixConstInterface< t_DataType> &MATRIX)
      {
        Compute( MATRIX, true);
        return m_Info == e_Success;
      }

      //! @brief Reports whether previous computation was successful.
      //!
      //! @returns \c Success if computation was succesful, \c Uninitialized if compute has not yet been called,
      //!          \c NoConvergence otherwise.
      Status Info() const
      {
        return m_Info;
      }

      //! @brief Determine whether this object tried to compute eigenvectors last time compute was called
      bool GetWereEigenvectorsRequested() const
      {
        return m_HaveVectors;
      }

    protected:

      //! @brief Computes eigendecomposition of given matrix.
      //!
      //! @param[in]  matrix  Selfadjoint matrix whose eigendecomposition is to
      //!    be computed. Only the lower triangular part of the matrix is referenced.
      //! @param[in]  WANT_EIGENVECTORS true to compute eigenvectors as well as values
      //! @returns    Reference to \c *this
      //!
      //! This function computes the eigenvalues of \p matrix.  The eigenvalues()
      //! function can be used to retrieve them.  If \p options equals #ComputeEigenvectors,
      //! then the eigenvectors are also computed and can be retrieved by
      //! calling eigenvectors().
      //!
      //! This implementation uses a symmetric QR algorithm. The matrix is first
      //! reduced to tridiagonal form using the Tridiagonalization class. The
      //! tridiagonal matrix is then brought to diagonal form with implicit
      //! symmetric QR steps with Wilkinson shift. Details can be found in
      //! Section 8.3 of Golub \& Van Loan, <i>%Matrix Computations</i>.
      //!
      //! The cost of the computation is about \f$ 9n^3 \f$ if the eigenvectors
      //! are required and \f$ 4n^3/3 \f$ if they are not required.
      //!
      //! This method reuses the memory in the SelfAdjointEigenSolver object that
      //! was allocated when the object was constructed, if the size of the
      //! matrix does not change.
      void Compute( const MatrixConstInterface< t_DataType> &MATRIX, const bool &WANT_EIGENVECTORS);

      //! @brief Performs a QR step on a tridiagonal symmetric matrix represented as a pair of two vectors
      //!        \a diag and \a subdiag.
      //!
      //! For compilation efficiency reasons, this procedure does not use eigen expression
      //! for its arguments.
      //!
      //! Implemented from Golub's "Matrix Computations", algorithm 8.3.2:
      //! "implicit symmetric QR step with Wilkinson shift"
      void TridiagonalQRStep( size_t START_INDEX, size_t END_INDEX);

      //! @brief Performs a full tridiagonalization in place
      //!
      //! @post * m_Eivec will have its tridiagonal decomposition is to be computed.  If computing eigenvectors, the
      //!         orthogonal matrix Q in the decomposition is stored in m_Eivec afterwards
      //!       * m_Eigenvalues is the diagonal of the tridiagonal matrix T in the decomposition.
      //!       * m_Subdiag  The subdiagonal of the tridiagonal matrix T in the decomposition.
      //!
      //! Computes the tridiagonal decomposition of the selfadjoint matrix \p m_Eivec in place
      //! such that \f$ m_Eivec = Q T Q^* \f$ where \f$ Q \f$ is unitary and \f$ T \f$ a real
      //! symmetric tridiagonal matrix.
      //!
      //! The tridiagonal matrix T is passed to the output parameters m_Eivalues and m_Subdiag. If
      //! extracting eigenvectors, then m_Eivec is set to the orthogonal matrix Q.  m_Eivec is otherwise consumed
      //! (rendered meaningless) by this operation
      //!
      //! Implemented from Golub's "Matrix Computations", algorithm 8.3.2:
      //! "implicit symmetric QR step with Wilkinson shift"
      void TridiagonalizationInPlace();

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        io::Serialize::Read( m_Eivec, ISTREAM);
        io::Serialize::Read( m_Eivalues, ISTREAM);
        io::Serialize::Read( m_Subdiag, ISTREAM);
        io::Serialize::Read( m_HaveVectors, ISTREAM);
        m_Tmp = m_Eivec;
        // end
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT indentation
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write NumberRows and -Cols
        io::Serialize::Write( m_Eivec, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Eivalues, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_Subdiag, OSTREAM, INDENT) << '\n';
        io::Serialize::Write( m_HaveVectors, OSTREAM, INDENT);

        //return
        return OSTREAM;
      }

    };

    namespace
    {

      //! @brief detect whether one value is much smaller than another
      //! @param A, B the values to compare
      //! @note specialization for float types
      inline bool IsMuchSmallerThan( const float &A, const float &B)
      {
        return math::Absolute( A) <= math::Absolute( B) * 1e-5f;
      }

      //! @brief detect whether one value is much smaller than another
      //! @param A, B the values to compare
      //! @note specialization for double types
      inline bool IsMuchSmallerThan( const double &A, const double &B)
      {
        return math::Absolute( A) <= math::Absolute( B) * 1e-12f;
      }

      //! @brief detect whether one value is much smaller than another
      //! @param A, B the values to compare
      //! @note specialization for integral types
      template< typename t_DataType>
      inline typename type::EnableIf< std::numeric_limits< t_DataType>::is_exact, bool>::Type
        IsMuchSmallerThan( const t_DataType &A, const t_DataType &B)
      {
        return !A && B;
      }

      //! @brief perform householder operations in-place on a vector
      //! @param VECTOR vector interface to make the householder coefficients out of
      //! @param TAU computed scaling factor of the Householder transform
      //! @param BETA result of householder * *this
      template< typename t_DataType>
      void MakeHouseholderInPlace( VectorInterface< t_DataType> &VECTOR, t_DataType &TAU, t_DataType &BETA)
      {
        VectorReference< t_DataType> essential_part( VECTOR.GetSize() - 1, VECTOR.Begin() + 1);
        const t_DataType tail_sq_norm( essential_part.SquareNorm());
        const t_DataType first( VECTOR( 0));
        if( tail_sq_norm == t_DataType( 0))
        {
          TAU = 0;
          BETA = first;
          essential_part = t_DataType( 0);
        }
        else
        {
          BETA = math::Sqrt( math::Sqr( first) + tail_sq_norm);
          if( first >= t_DataType( 0))
          {
            BETA = -BETA;
          }
          essential_part /= ( first - BETA);
          TAU = ( BETA - first) / BETA;
        }
      }
    }

    //! @brief Performs a full tridiagonalization in place
    //!
    //! @post * m_Eivec will have its tridiagonal decomposition is to be computed.  If computing eigenvectors, the
    //!         orthogonal matrix Q in the decomposition is stored in m_Eivec afterwards
    //!       * m_Eigenvalues is the diagonal of the tridiagonal matrix T in the decomposition.
    //!       * m_Subdiag  The subdiagonal of the tridiagonal matrix T in the decomposition.
    //!
    //! Computes the tridiagonal decomposition of the selfadjoint matrix \p m_Eivec in place
    //! such that \f$ m_Eivec = Q T Q^* \f$ where \f$ Q \f$ is unitary and \f$ T \f$ a real
    //! symmetric tridiagonal matrix.
    //!
    //! The tridiagonal matrix T is passed to the output parameters m_Eivalues and m_Subdiag. If
    //! extracting eigenvectors, then m_Eivec is set to the orthogonal matrix Q.  m_Eivec is otherwise consumed
    //! (rendered meaningless) by this operation
    //!
    //! Implemented from Golub's "Matrix Computations", algorithm 8.3.2:
    //! "implicit symmetric QR step with Wilkinson shift"
    template< typename t_DataType>
    void SymmetricEigenSolver< t_DataType>::TridiagonalizationInPlace()
    {
      const size_t n( m_Eivec.GetNumberRows());
      Vector< t_DataType> h_coeffs( n - 1);
      for( size_t i( 0); i < n - 1; ++i)
      {
        size_t remaining_size( n - i - 1);
        t_DataType beta;
        t_DataType h;
        VectorReference< t_DataType> mat_row_i_tail( m_Eivec.GetRow( i).Slice( i + 1));
        MakeHouseholderInPlace( mat_row_i_tail, h, beta);

        // Apply similarity transformation to remaining columns,
        // i.e., A = H A H' where H = I - h v v' and v = matA.col(i).tail(n-i-1)
        m_Eivec( i, i + 1) = 1;

        VectorReference< t_DataType> hcoeffs_tail( h_coeffs.Slice( i));
        for( size_t r( 0); r < remaining_size; ++r)
        {
          hcoeffs_tail( r) = ScalarProduct( m_Eivec.GetRow( r + i + 1).Slice( i + 1), mat_row_i_tail) * h;
        }
        hcoeffs_tail += ( h * t_DataType( -0.5) * ScalarProduct( hcoeffs_tail, mat_row_i_tail)) * mat_row_i_tail;
        for( size_t r( 0); r < remaining_size; ++r)
        {
          const t_DataType mat_row_i_col_r( mat_row_i_tail( r)), hcoeffs_tail_col_r( hcoeffs_tail( r));
          t_DataType *matrix_rpi_upper( m_Eivec[ r + i + 1] + i + 1);
          for( size_t s( 0); s < remaining_size; ++s, ++matrix_rpi_upper)
          {
            *matrix_rpi_upper -= mat_row_i_col_r * hcoeffs_tail( s) + mat_row_i_tail( s) * hcoeffs_tail_col_r;
          }
        }
        m_Subdiag( i) = m_Eivec( i, i + 1) = beta;
        h_coeffs( i) = h;
        m_Eivalues( i) = m_Eivec( i, i);
      }
      m_Eivalues( n - 1) = m_Eivec( n - 1, n - 1);
      if( m_HaveVectors)
      {
        // This entire scope comes from the following line in Eigen:
        // m_Eivec = HouseholderSequence< t_DataType>( m_Eivec, hCoeffs).setLength( n - 1).setShift( 1);
        Vector< t_DataType> workspace( n);
        m_Tmp.SetIdentity();
        // make the upper triangle of the matrix == the identity matrix
        // m_Eivec.diagonal().setOnes();
        // m_Eivec.template triangularView<StrictlyUpper>().setZero();
        // This code is translated to the equivalent BCL code
        for( size_t k( n - 2); k < n; --k)
        {
          const size_t corner_size( n - k - 1);
          // code for essential vector
          // Block<const VectorsType,Dynamic,1>(m_Eivec, k+2, k, n-k-2, 1);
          VectorConstReference< t_DataType> essential_vector_reference( m_Eivec.GetRow( k).Slice( k + 2));
          // Code from eigen:
          // m_Eivec.bottomRightCorner(cornerSize, cornerSize)
          //        .applyHouseholderOnTheLeft(essential_vector_reference, hCoeffs(k), workspace);
          const t_DataType tau( h_coeffs( k));
          if( corner_size == 1)
          {
            const t_DataType one_minus_tau( 1 - tau);
            m_Tmp( n - 1, n - 2) *= one_minus_tau;
          }
          else
          {
            MatrixReference< t_DataType> bottom_right_eig( corner_size, n, m_Tmp[ k + 1]);
            const size_t cornersize_next( corner_size - 1);
            // Block<Derived, EssentialPart::SizeAtCompileTime, Derived::ColsAtCompileTime> bottom(derived(), 1, 0, rows()-1, cols());
            VectorReference< t_DataType> workspace_ref( corner_size, workspace.Begin());
            workspace_ref = t_DataType( 0);

            // get a pointer to the beginning of the array
            const t_DataType *feature_ref( essential_vector_reference.Begin());

            // STORAGE(i) = ( sum of column i of MATRIX) * FEATURE( i)
            // start by setting STORAGE(i) = ( sum of column i of MATRIX)
            for( size_t i( 0); i < cornersize_next; ++i)
            {
              const t_DataType *itr_row( bottom_right_eig[ i + 1] + k + 1);
              const t_DataType vector_value( feature_ref[ i]);
              for
              (
                t_DataType *itr_result( workspace_ref.Begin()), *itr_result_end( workspace_ref.End());
                itr_result != itr_result_end;
                ++itr_result, ++itr_row
              )
              {
                *itr_result += vector_value * ( *itr_row);
              }
            }
            //VectorEqualsVectorTimesMatrix( workspace_ref, essential_vector_reference, bottom_right_bottom);
            workspace_ref += m_Tmp.GetRow( k + 1).Slice( k + 1);
            workspace_ref *= tau;
            VectorReference< t_DataType> eigenvector_bottom_right_upper_triangle_first_row
            (
              m_Tmp.GetRow( k + 1).Slice( k + 1)
            );
            eigenvector_bottom_right_upper_triangle_first_row -= workspace_ref;
            for( size_t x( 0); x < cornersize_next; ++x)
            {
              t_DataType *itr_tmp( m_Tmp[ x + k + 2] + k + 1);
              for( size_t y( 0); y < corner_size; ++y, ++itr_tmp)
              {
                *itr_tmp -= essential_vector_reference( x) * workspace( y);
              }
            }
          }
        }
        m_Eivec = m_Tmp;
      }
    }

    template< typename t_DataType>
    void SymmetricEigenSolver< t_DataType>::Compute
    (
      const MatrixConstInterface< t_DataType> &MATRIX,
      const bool &WANT_EIGENVECTORS
    )
    {
      BCL_Assert
      (
        MATRIX.IsSymmetric(),
        "SelfAdjointEigenSolver cannot compute eigendecomposition of non-symmetric matrix"
      );
      m_HaveVectors = WANT_EIGENVECTORS;
      const size_t n( MATRIX.GetNumberCols());
      if( m_Eivalues.GetSize() != n)
      {
        m_Eivalues = Vector< t_DataType>( n, t_DataType( 0));
      }
      else
      {
        m_Eivalues = t_DataType( 0);
      }
      if( WANT_EIGENVECTORS)
      {
        if( m_Eivec.GetNumberRows() != n)
        {
          m_Eivec = Matrix< t_DataType>( n, n, t_DataType( 0));
        }
        else
        {
          m_Eivec = t_DataType( 0);
        }
        if( m_Tmp.GetNumberRows() != n)
        {
          m_Tmp = Matrix< t_DataType>( n, n, t_DataType( 0));
        }
      }
      m_Info = e_Uninitialized;

      // trivial case of scalar matrix
      if( n == 1)
      {
        m_Eivalues( 0) = MATRIX( 0, 0);
        if( m_HaveVectors)
        {
          m_Eivec( 0, 0) = t_DataType( 1);
        }
        m_Info = e_Success;
        return;
      }

      // map the matrix coefficients to [-1:1] to avoid over- and underflow.
      m_Eivec = MATRIX;
      t_DataType scale = std::max( -MATRIX.AsVector().Min(), MATRIX.AsVector().Max());
      if( scale == t_DataType( 0))
      {
        scale = t_DataType( 1);
      }
      m_Eivec /= scale;
      m_Subdiag = Vector< t_DataType>( n - 1, t_DataType( 0));
      TridiagonalizationInPlace();

      for( size_t i( 0), last( n - 1); i < last; ++i)
      {
        if( IsMuchSmallerThan( m_Subdiag( i), t_DataType( math::Absolute(m_Eivalues(i))+math::Absolute(m_Eivalues(i+1)))))
        {
          m_Subdiag( i) = 0;
        }
      }

      bool finished( false);
      for( size_t start( 0), end( n - 1), iter( 0), max_iter( s_MaxIterations * n); iter <= max_iter; ++iter)
      {
        // find the largest unreduced block
        while( end > 0 && m_Subdiag( end - 1) == 0)
        {
          --end;
        }
        if( end == 0)
        {
          finished = true;
          break;
        }

        for( start = end - 1; start > 0 && m_Subdiag( start - 1) != 0; --start)
        {
        }

        TridiagonalQRStep( start, end);

        // zero out tiny components of m_Subdiag
        for( size_t i( start); i < end; ++i)
        {
          if( IsMuchSmallerThan( m_Subdiag( i), t_DataType( math::Absolute(m_Eivalues(i))+math::Absolute(m_Eivalues(i+1)))))
          {
            m_Subdiag( i) = 0;
          }
        }
      }

      m_Info = finished ? e_Success : e_NoConvergence;

      // Sort eigenvalues and corresponding vectors.
      if( m_Info == e_Success)
      {
        if( !m_HaveVectors)
        {
          // sort the eigenvalues directly
          std::sort( m_Eivalues.Begin(), m_Eivalues.End(), std::greater< t_DataType>());
        }
        else
        {
          // sort the by eigenvalue / index pair, then reorder the eigenvector matrix accordingly
          std::vector< std::pair< t_DataType, size_t> > diagonal;
          diagonal.reserve( n);
          for( size_t row( 0); row < n; ++row)
          {
            diagonal.push_back( std::make_pair( m_Eivalues( row), row));
          }
          std::stable_sort( diagonal.begin(), diagonal.end());
          storage::Vector< size_t> neworder;
          neworder.AllocateMemory( n);
          for( size_t row( 0), remaining_values( n - 1); row < n; ++row, --remaining_values)
          {
            neworder.PushBack( diagonal[ remaining_values].second);
            m_Eivalues( row) = diagonal[ remaining_values].first;
          }
          m_Eivec.Transpose();
          m_Eivec.ReorderRows( neworder);
          m_Eivec.Transpose();
        }
      }

      // scale back the eigen values
      m_Eivalues *= scale;
      return;
    }

    template< typename t_DataType>
    void SymmetricEigenSolver< t_DataType>::TridiagonalQRStep( size_t START_INDEX, size_t END_INDEX)
    {
      t_DataType *diag( m_Eivalues.Begin());
      t_DataType *subdiag( m_Subdiag.Begin());
      t_DataType td( ( diag[ END_INDEX - 1] - diag[ END_INDEX]) * t_DataType( 0.5));
      t_DataType e( subdiag[ END_INDEX - 1]);
      // Note that thanks to scaling, e^2 or td^2 cannot overflow, however they can still
      // underflow thus leading to inf/NaN values when using the following commented code:
      //   t_DataType e2 = math::Sqr(subdiag[end-1]);
      //   t_DataType mu = diag[end] - e2 / (td + (td>0 ? 1 : -1) * sqrt(td*td + e2));
      // This explain the following, somewhat more complicated, version:
      t_DataType mu( diag[ END_INDEX]);
      if(td==0)
      {
        mu -= math::Absolute(e);
      }
      else
      {
        t_DataType e2( math::Sqr( subdiag[ END_INDEX - 1]));
        t_DataType hypotneuse( math::Pythag( td, e));
        if( e2 == 0)
        {
          mu -= ( e / ( td + ( td > 0 ? 1 : -1))) * ( e / hypotneuse);
        }
        else
        {
          mu -= e2 / ( td + ( td > 0 ? hypotneuse : -hypotneuse));
        }
      }

      t_DataType x( diag[ START_INDEX] - mu);
      t_DataType z( subdiag[ START_INDEX]);
      math::RotationMatrix2D rot;
      for( size_t k( START_INDEX); k < END_INDEX; ++k)
      {
        rot.MakeGivens( x, z);
        const t_DataType c( rot.GetMatrix()( 0, 0)); // cosine component of transformation
        const t_DataType s( rot.GetMatrix()( 1, 0)); // sine component; (0,1) would be -sin component
        // do T = G' T G
        const t_DataType sdk( s * diag[ k] + c * subdiag[ k]);
        const t_DataType dkp1( s * subdiag[ k] + c * diag[ k + 1]);

        diag[ k] = c * ( c * diag[ k] - s * subdiag[ k]) - s * ( c * subdiag[ k] - s * diag[ k + 1]);
        diag[ k + 1] = s * sdk + c * dkp1;
        subdiag[ k] = c * sdk - s * dkp1;

        if( k > START_INDEX)
        {
          subdiag[ k - 1] = c * subdiag[ k - 1] - s * z;
        }

        x = subdiag[ k];

        if( k < END_INDEX - 1)
        {
          z = -s * subdiag[ k + 1];
          subdiag[ k + 1] = c * subdiag[ k + 1];
        }

        // apply the givens rotation to the unit matrix Q = Q * G
        if( m_HaveVectors)
        {
          // this operation could be more efficient if re-coded for cache performance on BCL's row-ordered matrices
          t_DataType *col_k( m_Eivec.Begin() + k), *col_kn( col_k + 1);
          for( size_t i( 0), n( m_Eivalues.GetSize()); i < n; ++i, col_k += n, col_kn += n)
          {
            t_DataType xi( *col_k), yi( *col_kn);
            *col_k = c * xi - s * yi;
            *col_kn = s * xi + c * yi;
          }
        }
      }
    }

    //! performs SVD of any matrix via MOORE-PENROSE-INVERSION of rectangular or non-
    //! symmetric square matrices by solving the eigen system of m * Transpose( m)
    //! this[rows,cols] = EigenfunctionsU[rows,cols] * Eigenvalues[cols,cols] * EigenfunctionsV[cols,cols]
    template< typename t_DataType>
    Matrix< t_DataType>
    SingularValueDecomposition
    (
      const MatrixConstInterface< t_DataType> &MATRIX,
      Matrix< t_DataType> &EIGEN_VECTORS_V,
      Matrix< t_DataType> &EIGEN_VECTORS_U
    )
    {
      // solve eigensystem of Transpose( *this) * ( *this)
      Matrix< t_DataType> mt_times_m( MatrixTransposeTimesMatrix( MATRIX));
      SymmetricEigenSolver< t_DataType> eigensolver( mt_times_m, !EIGEN_VECTORS_V.IsEmpty());
      if( eigensolver.GetWereEigenvectorsRequested())
      {
        EIGEN_VECTORS_V = eigensolver.GetSortedEigenvectors();
      }

      Vector< t_DataType> eigenvalues( eigensolver.GetSortedEigenvalues());

      // calculating eigenvalues
      for( size_t i( 0); i < eigenvalues.GetSize(); ++i)
      {
        eigenvalues( i) >= t_DataType( 0)
          ? eigenvalues( i) = math::Sqrt( eigenvalues( i))
          : eigenvalues( i) = math::Sqrt( -eigenvalues( i));
      }

      Matrix< t_DataType> &eigenvalue_diagonal_matrix( mt_times_m);
      eigenvalue_diagonal_matrix.SetIdentity();

      // calculate U
      if( !EIGEN_VECTORS_U.IsEmpty())
      {
        for( size_t i( 0), nr( eigenvalue_diagonal_matrix.GetNumberRows()); i < nr; ++i)
        {
          eigenvalue_diagonal_matrix( i, i) =
            eigenvalues( i) > std::numeric_limits< t_DataType>::min() ? 1.0 / eigenvalues( i) : 0.0;
        }
        EIGEN_VECTORS_U = MATRIX * EIGEN_VECTORS_V * eigenvalue_diagonal_matrix;
      }
      ReplaceDiagonal( eigenvalue_diagonal_matrix, eigenvalues);

      // invert V
      if( !EIGEN_VECTORS_V.IsEmpty())
      {
        EIGEN_VECTORS_V.Transpose();
      }

      // return
      return eigenvalue_diagonal_matrix;
    }

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_SYMMETRIC_EIGENSOLVER_H_
