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

#ifndef BCL_LINAL_PRINCIPAL_COMPONENT_ANALYSIS_H_
#define BCL_LINAL_PRINCIPAL_COMPONENT_ANALYSIS_H_

// include the namespace header
#include "bcl_linal.h"

// includes from bcl - sorted alphabetically
#include "bcl_linal_matrix.h"
#include "bcl_linal_matrix_operations.h"
#include "bcl_linal_operations.h"
#include "bcl_linal_symmetric_eigensolver.h"
#include "bcl_linal_vector.h"
#include "storage/bcl_storage_pair.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace linal
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PrincipalComponentAnalysis
    //! @brief performs principal component analysis on input matrix
    //!
    //! @tparam t_DataType can be int, float, double, or complex
    //!
    //! @see @link example_linal_principal_component_analysis.cpp @endlink
    //! @author loweew
    //! @date Jul 4, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class PrincipalComponentAnalysis :
      public util::ObjectInterface
    {

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new PrincipalComponentAnalysis< t_DataType>
      PrincipalComponentAnalysis< t_DataType> *Clone() const
      {
        return new PrincipalComponentAnalysis< t_DataType>( *this);
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

    ////////////////
    // operations //
    ////////////////

      //! @brief sorts eigenvalues and corresponding eigenvectors v
      //! @param EIGEN_VALUES the eigenvalues to be sorted
      //! @param EIGEN_VECTORS_V the corresponding eigenvectors to be sorted according to eigenvalues
      storage::VectorND< 2, Matrix< t_DataType> > SortEigenValuesAndReduce
      (
        const MatrixConstInterface< t_DataType> &EIGEN_VALUES,
        const MatrixConstInterface< t_DataType> &EIGEN_VECTORS_V,
        const t_DataType THRESHOLD
      )
      {
        BCL_Assert( EIGEN_VALUES.IsSquare(), "eigen values matrix is not square!");

        storage::Vector< storage::Pair< t_DataType, size_t> > diagonal;
        diagonal.AllocateMemory( EIGEN_VALUES.GetNumberRows());
        for( size_t row( 0), row_end( EIGEN_VALUES.GetNumberRows()); row < row_end; ++row)
        {
          diagonal.PushBack( storage::Pair< t_DataType, size_t>( EIGEN_VALUES( row, row), row));
        }

        diagonal.Sort( std::greater< storage::Pair< t_DataType, size_t> >());
        storage::Vector< storage::Pair< t_DataType, size_t> >
          reduced_eigen_values( ReduceEigenValues( diagonal, THRESHOLD));
        const size_t new_cols( reduced_eigen_values.GetSize());
        Matrix< t_DataType> reduced_eig_mat( new_cols, new_cols);
        Vector< t_DataType> sorted_eig_values( new_cols);
        for( size_t ctr( 0); ctr < new_cols; ++ctr)
        {
          sorted_eig_values( ctr) = reduced_eigen_values( ctr).First();
        }
        ReplaceDiagonal( reduced_eig_mat, sorted_eig_values);

        Matrix< t_DataType> trans_reduced_eigen_vectors( EIGEN_VECTORS_V.GetNumberRows(), new_cols);
        for( size_t itr( 0); itr < new_cols; ++itr)
        {
          const size_t position( reduced_eigen_values( itr).Second());
          std::copy( EIGEN_VECTORS_V[ position], EIGEN_VECTORS_V[ position] + new_cols, trans_reduced_eigen_vectors[ itr]);
        }
        return storage::VectorND< 2, Matrix< t_DataType> >( reduced_eig_mat, trans_reduced_eigen_vectors.Transposed());
      }

      //! @brief sorts eigenvalues and corresponding eigenvectors v
      //! @param EIGEN_VALUES the eigenvalues to be sorted
      //! @param EIGEN_VECTORS_V the corresponding eigenvectors sorted according to eigenvalues
      storage::VectorND< 2, Matrix< t_DataType> > SortEigenValues
      (
        const MatrixConstInterface< t_DataType> &EIGEN_VALUES,
        const MatrixConstInterface< t_DataType> &EIGEN_VECTORS_V
      )
      {
        BCL_Assert( EIGEN_VALUES.IsSquare(), "eigen values matrix is not square!");

        storage::Vector< storage::Pair< t_DataType, size_t> > diagonal( EIGEN_VALUES.GetNumberRows());
        const size_t rows( diagonal.GetSize());
        for( size_t row( 0); row < rows; ++row)
        {
          diagonal.PushBack( storage::Pair< t_DataType, size_t>( EIGEN_VALUES( row, row), row));
        }

        diagonal.Sort( std::greater< storage::Pair< t_DataType, size_t> >());
        Vector< t_DataType> sorted_eig_values( rows);
        for( size_t ctr( 0); ctr < rows; ++ctr)
        {
          sorted_eig_values( ctr) = diagonal( ctr).First();
        }
        Matrix< t_DataType> sorted_eig_mat( rows, rows);
        ReplaceDiagonal( sorted_eig_mat, sorted_eig_values);

        Matrix< t_DataType> trans_eigen_vectors( rows, rows);
        for( size_t itr( 0); itr < rows; ++itr)
        {
          const size_t position( diagonal( itr).Second());
          std::copy
          (
            EIGEN_VECTORS_V.Begin() + ( position * rows),
            EIGEN_VECTORS_V.Begin() + ( ( position + 1) * rows),
            trans_eigen_vectors.Begin() + itr * rows
          );
        }

        return storage::VectorND< 2, Matrix< t_DataType> >( sorted_eig_mat, trans_eigen_vectors.Transposed());
      }

      //! @brief determines number of eigenvalues to keep based on threshold
      //! @param EIGEN_VALUES storage vector of sorted eigenvalues with indices
      //! @param THRESHOLD the cutoff value for how many to keep in percent * 100.0
      //! @return reduced eigenvalues with corresponding indices
      storage::Vector< storage::Pair< t_DataType, size_t> > ReduceEigenValues
      (
        storage::Vector< storage::Pair< t_DataType, size_t> > &EIGEN_VALUES,
        const float THRESHOLD
      )
      {
        const size_t nr_rows( EIGEN_VALUES.GetSize());
        t_DataType variance( 0);
        for( size_t count( 0); count < nr_rows; ++count)
        {
          variance += EIGEN_VALUES( count).First();
        }

        const t_DataType desired_variance( THRESHOLD / 100.0 * variance);
        storage::Vector< storage::Pair< t_DataType, size_t> > reduced_eig;
        t_DataType sum( 0);
        for( size_t count( 0); count < nr_rows; ++count)
        {
          sum += EIGEN_VALUES( count).First();
          reduced_eig.PushBack( EIGEN_VALUES( count));
          if( sum >= desired_variance)
          {
            break;
          }
        }

        return reduced_eig;
      }

      //! @brief determines number of eigenvalues to keep based on threshold
      //! @param EIGEN_VALUES storage vector of sorted eigenvalues with indices
      //! @param N the number of desired eigen vectors
      //! @return reduced eigenvalues with corresponding indices
      storage::Vector< storage::Pair< t_DataType, size_t> > ReduceEigenValuesToN
      (
        storage::Vector< storage::Pair< t_DataType, size_t> > &EIGEN_VALUES,
        const size_t N
      )
      {
        return storage::Vector< storage::Pair< t_DataType, size_t> >( EIGEN_VALUES.Begin(), EIGEN_VALUES.Begin() + N);
      }

      //! @brief returns the reduced sorted eigenvectors and values
      static storage::Pair< Matrix< t_DataType>, Vector< t_DataType> >
      GetSortedEigenVectorsValues( const Matrix< t_DataType> &INPUT)
      {
        const size_t cols( INPUT.GetNumberCols());
        Matrix< t_DataType> eig_vec_v( cols, cols), eig_vec_u;
        Matrix< t_DataType> eig_values( SingularValueDecomposition( INPUT, eig_vec_v, eig_vec_u));

        BCL_Assert( eig_values.IsSquare(), "eigen values matrix is not square!");
        storage::Vector< storage::Pair< t_DataType, size_t> > diagonal;
        diagonal.AllocateMemory( cols);

        for( size_t row( 0); row < cols; ++row)
        {
          diagonal.PushBack( storage::Pair< t_DataType, size_t>( eig_values( row, row), row));
        }

        diagonal.Sort( std::greater< storage::Pair< t_DataType, size_t> >());

        Matrix< t_DataType> reduced_eig_mat( cols, cols);
        Vector< t_DataType> sorted_eig_values( cols);
        for( size_t ctr( 0); ctr < cols; ++ctr)
        {
          sorted_eig_values( ctr) = diagonal( ctr).First();
        }

        ReplaceDiagonal( reduced_eig_mat, sorted_eig_values);

        Matrix< t_DataType> trans_reduced_eigen_vectors( cols, cols);
        for( size_t itr( 0); itr < cols; ++itr)
        {
          const size_t position( diagonal( itr).Second());
          std::copy( eig_vec_v[ position], eig_vec_v[ position] + cols, trans_reduced_eigen_vectors[ itr]);
        }

        return
          storage::Pair< Matrix< t_DataType>, Vector< t_DataType> >( trans_reduced_eigen_vectors, sorted_eig_values);
      }

      //! @brief determines number of eigenvalues to keep based on threshold and max components to keep
      //! @param INPUT the input matrix
      //! @param STORAGE matrix in which to store the result
      //! @param EIGEN_VALUES storage vector of sorted eigenvalues with indices
      //! @param THRESHOLD the cutoff value for how many to keep in fraction of total (1.0 max)
      //! @param MAX_TO_KEEP maximum number of values to keep
      //! @return reduced eigenvalues with corresponding indices
      static void ReduceInputMatrix
      (
        const MatrixConstInterface< t_DataType> &INPUT,
        Matrix< t_DataType> &STORAGE,
        const Matrix< t_DataType> &EIGEN_VECTORS,
        const Vector< t_DataType> &EIGEN_VALUES,
        const t_DataType &THRESHOLD,
        const size_t &MAX_COMPONENTS
      )
      {
        const size_t nr_rows( std::min( MAX_COMPONENTS, EIGEN_VALUES.GetSize()));

        const t_DataType desired_variance( THRESHOLD * EIGEN_VALUES.Sum());

        size_t new_cols( 0);
        for( t_DataType sum( 0); new_cols < nr_rows && sum < desired_variance; ++new_cols)
        {
          sum += EIGEN_VALUES( new_cols);
        }

        Matrix< t_DataType> trans_reduced_eigen_vectors
        (
          Matrix< t_DataType>( new_cols, EIGEN_VECTORS.GetNumberRows(), EIGEN_VECTORS.Begin()).Transposed()
        );

        STORAGE = GetDefaultOperations< float>().Multiply( INPUT, trans_reduced_eigen_vectors);
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator returns the reduced sorted eigenvectors
      //! @param INPUT the input matrix
      //! @param THRESHOLD the threshold for how many eigenvectors to return
      //! @return the eigenvectors matrix
      Matrix< t_DataType> operator()( Matrix< t_DataType> &INPUT, const t_DataType &THRESHOLD)
      {
        Matrix< t_DataType> eig_vec_v( INPUT.GetNumberCols(), INPUT.GetNumberCols()), eig_vec_u;

        Matrix< t_DataType> eig_values( SingularValueDecomposition( INPUT, eig_vec_v, eig_vec_u));
        BCL_Assert( eig_values.IsDiagonal(), "eig values are not diagonal and should be!");

        storage::VectorND< 2, Matrix< t_DataType> > sorted_eig_vals_and_eig_vec_v
        (
          SortEigenValues( eig_values, eig_vec_v)
        );

        return sorted_eig_vals_and_eig_vec_v.Second();
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        // read members

        // return the stream
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write members

        // return the stream
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // template class PrincipalComponentAnalysis

    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> PrincipalComponentAnalysis< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( PrincipalComponentAnalysis< t_DataType>())
    );

  } // namespace linal
} // namespace bcl

#endif // BCL_LINAL_PRINCIPAL_COMPONENT_ANALYSIS_H_
