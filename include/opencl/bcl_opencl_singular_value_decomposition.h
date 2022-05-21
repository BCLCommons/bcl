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

#ifndef BCL_OPENCL_SINGULAR_VALUE_DECOMPOSITION_H_
#define BCL_OPENCL_SINGULAR_VALUE_DECOMPOSITION_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_kernel_source_file.h"
#include "bcl_opencl_kernel_source_interface.h"
#include "bcl_opencl_kernel_sources.h"
#include "bcl_opencl_matrix.h"
#include "bcl_opencl_matrix_multiply.h"
#include "bcl_opencl_matrix_transpose.h"
#include "bcl_opencl_vector.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_symmetric_eigensolver.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SingularValueDecomposition
    //! @brief performs singular value decomposition
    //! @details
    //! compute eigensystem of a square symmetric matrix via tridiagonalization
    //! this[cols,cols] = Eigenfunctions[cols,cols] * Eigenvalues[cols,cols] * Transpose( Eigenfunctions[cols,cols])
    //! numerical recipes in C *** 11.2 & 11.3 pages 469..481
    //!
    //! @tparam t_DataType float, int, double
    //!
    //! @see @link example_opencl_singular_value_decomposition.cpp @endlink
    //! @author loweew
    //! @date Aug 5, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class SingularValueDecomposition :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! host eigensystem
      linal::Matrix< t_DataType> m_EigenValues;
      linal::Matrix< t_DataType> m_EigenVectorsV;
      linal::Matrix< t_DataType> m_EigenVectorsU;

      //! command queue
      CommandQueue m_Queue;

      //! ocl matrix multiply object
      MatrixMultiply< t_DataType> m_MMult;

      //! ocl transpose object
      MatrixTranspose< t_DataType> m_MTranspose;

      //! ocl program
      cl::Program m_Program;

      //! blocksize
      static const cl_uint s_blocksize = 16;

      //! compiler options
      static const char *s_CLCompilerOptions;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SingularValueDecomposition()
      {
      }

      //! @brief constructor
      //! @param QUEUE the command queue
      SingularValueDecomposition( const CommandQueue &QUEUE) :
        m_Queue( QUEUE),
        m_MMult( QUEUE),
        m_MTranspose( QUEUE)
      {
        cl_int error_number( CL_SUCCESS);
        m_Program = KernelSources::Compile( GetKernelSources().e_Linal, util::CPPDataTypes::DataTypeFromTemplate< t_DataType>(), m_Queue, std::string( s_CLCompilerOptions), &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
      }

      //! @brief Clone function
      //! @return pointer to new SingularValueDecomposition< t_DataType>
      SingularValueDecomposition< t_DataType> *Clone() const
      {
        return new SingularValueDecomposition< t_DataType>( *this);
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

      //! @brief get eigen values
      //! @return eigen values
      const linal::Matrix< t_DataType> &GetEigenValues() const
      {
        return m_EigenValues;
      }

      //! @brief get eigen vectors v
      //! @return eigen vectors v
      const linal::Matrix< t_DataType> &GetEigenVectorsV() const
      {
        return m_EigenVectorsV;
      }

      //! @brief get eigen vectors u
      //! @return eigen vectors u
      const linal::Matrix< t_DataType> &GetEigenVectorsU() const
      {
        return m_EigenVectorsU;
      }

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief performs singular value decomposition
      //! @param MATRIX the matrix to operate on
      //! this[rows,cols] = EigenfunctionsU[rows,cols] * Eigenvalues[cols,cols] * EigenfunctionsV[cols,cols]
      void operator()
      (
        const linal::Matrix< t_DataType> &MATRIX
      )
      {
        // dims
        const size_t rows( MATRIX.GetNumberRows());
        const size_t cols( MATRIX.GetNumberCols());
        const size_t row_pad( ( s_blocksize - ( rows % s_blocksize)) % s_blocksize);
        const size_t col_pad( ( s_blocksize - ( cols % s_blocksize)) % s_blocksize);

        Matrix< t_DataType> device_eig_vals ( cols, cols, m_Queue, col_pad, col_pad);
        Matrix< t_DataType> device_eig_vec_u( rows, cols, m_Queue, row_pad, col_pad);

        // padded opencl matrix
        Matrix< t_DataType> input( MATRIX, m_Queue, row_pad, col_pad);

        // solve eigensystem of Transpose( mat) * ( mat)
        Matrix< t_DataType> trans_input( cols, rows, m_Queue, col_pad, row_pad);
        m_MTranspose( input, trans_input);
        m_MMult( trans_input, input, device_eig_vals);
        m_EigenValues = device_eig_vals.GetHostMatrix( col_pad, col_pad);

        // get eigen system
        m_EigenVectorsV = m_EigenValues;
        linal::SymmetricEigenSolver< t_DataType> eigensolver( m_EigenValues, true);
        m_EigenValues.SetIdentity();
        linal::ReplaceDiagonal( m_EigenValues, eigensolver.GetSortedEigenvalues());
        m_EigenVectorsV = eigensolver.GetSortedEigenvectors();

        // calculating eigenvalues
        for( size_t i( 0), end( m_EigenValues.GetNumberCols()); i < end; ++i)
        {
          m_EigenValues( i, i) >= t_DataType( 0)
            ? m_EigenValues( i, i) = math::Sqrt( m_EigenValues( i, i))
            : m_EigenValues( i, i) = math::Sqrt( -m_EigenValues( i, i));
        }

        // calculate U
        Matrix< t_DataType> device_eig_vec_v( m_EigenVectorsV, m_Queue, col_pad, col_pad);

        device_eig_vals = Matrix< t_DataType>( m_EigenValues, m_Queue, col_pad, col_pad);
        Matrix< t_DataType> tmp( rows, cols, m_Queue, row_pad, col_pad);

        m_MMult( input, device_eig_vec_v, tmp);
        m_MMult( tmp, device_eig_vals, device_eig_vec_u);
        m_EigenVectorsU = device_eig_vec_u.GetHostMatrix( row_pad, col_pad);

        // invert V
        Matrix< t_DataType> trans_eig_v( cols, cols, m_Queue, col_pad, col_pad);
        m_MTranspose( device_eig_vec_v, trans_eig_v);
        m_EigenVectorsV = trans_eig_v.GetHostMatrix( col_pad, col_pad);
        m_EigenValues = device_eig_vals.GetHostMatrix( col_pad, col_pad);
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

    }; // template class SingularValueDecomposition

    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> SingularValueDecomposition< t_DataType>::s_Instance;

    template< typename t_DataType>
    const char *SingularValueDecomposition< t_DataType>::s_CLCompilerOptions( "-cl-mad-enable -cl-fast-relaxed-math");

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_SINGULAR_VALUE_DECOMPOSITION_H_
