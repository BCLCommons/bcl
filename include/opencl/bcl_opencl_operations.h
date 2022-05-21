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

#ifndef BCL_OPENCL_OPERATIONS_H_
#define BCL_OPENCL_OPERATIONS_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_context.h"
#include "bcl_opencl_kernel_source_interface.h"
#include "linal/bcl_linal_operations_interface.h"
#include "linal/bcl_linal_vector_const_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Operations
    //!
    //! @author woetzen, loweew, vuot2
    //!
    //! @see @link example_opencl_operations.cpp @endlink
    //!
    //! @date Mar 14, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class Operations :
      public linal::OperationsInterface< t_DataType>
    {

    private:

    //////////
    // data //
    //////////

      CommandQueue   m_CommandQueue; //!< the command queue
      Context        m_Context;      //!< the context
      cl::Program    m_Program;      //!< the opencl program

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Operations();

      //! @brief Clone function
      //! @return pointer to new Operations
      Operations< t_DataType> *Clone() const;

      //! @brief is this class compatible with given command queue
      //! @param COMMAND_QUEUE the command queue this object would operate on
      //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
      bool IsCompatible( const CommandQueue &COMMAND_QUEUE) const;

      //! @brief initialize this class
      //! @brief COMMAND_QUEUE queue to use
      //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
      bool Initialize( const CommandQueue &COMMAND_QUEUE);

      //! @brief access to instance
      //! @return reference to instance of this class
      static const Operations< t_DataType> &GetInstance();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief access to the program
      //! @return reference to the program
      const cl::Program &GetProgram() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief dot product of two vectors
      //! @param VECTOR_A first vector
      //! @param VECTOR_B second vector
      //! @return dot product of VECTOR_A * VECTOR_B
      t_DataType DotProduct
      (
        const linal::VectorConstInterface< t_DataType> &VECTOR_A,
        const linal::VectorConstInterface< t_DataType> &VECTOR_B
      ) const;

      //! @brief norm of a vector
      //! @param VECTOR vector
      //! @return the euclidean norm of the VECTOR
      t_DataType Norm( const linal::VectorConstInterface< t_DataType> &VECTOR) const;

      //! @brief matrix-vector multiplication
      //! @param MATRIX matrix to be multiplied
      //! @param VECTOR vector to be multiplied
      //! @return resulting linal::Vector< t_DataType>
      linal::Vector< t_DataType> Multiply
      (
        const linal::MatrixConstInterface< t_DataType> &MATRIX,
        const linal::VectorConstInterface< t_DataType> &VECTOR
      ) const;

      //! @brief performs matrix-matrix multiplication using level 3 cublas
      //! @param MATRIX_A the first matrix in the multiplication
      //! @param MATRIX_B the second matrix in the multiplication
      //! @return the resulting linal::Matrix< t_DataType>
      linal::Matrix< t_DataType> Multiply
      (
        const linal::MatrixConstInterface< t_DataType> &MATRIX_A,
        const linal::MatrixConstInterface< t_DataType> &MATRIX_B
      ) const;

      //! @brief calculate pair-wise distances between lists of vectors
      //! @param LIST_VECTORS list of vectors
      //! @return triangular matrix of values of the distances
      linal::Matrix< t_DataType> DistanceMatrix
      (
        const util::SiPtrVector< const linal::VectorConstInterface< t_DataType> > &LIST_VECTORS
      ) const;

      //! @brief calculate pair-wise distances between lists of vectors
      //! @param LIST_VECTORS_A list of vectors
      //! @param LIST_VECTORS_B list of vectors
      //! @return matrix of values of the distances
      linal::Matrix< t_DataType> DistanceMatrix
      (
        const util::SiPtrVector< const linal::VectorConstInterface< t_DataType> > &LIST_VECTORS_A,
        const util::SiPtrVector< const linal::VectorConstInterface< t_DataType> > &LIST_VECTORS_B
      ) const;

      //! @brief reduction sum kernel
      //! @param VECTOR the vector to get the sum of
      //! @return the reduced result sum
      t_DataType Sum( const linal::VectorConstInterface< t_DataType> &VECTOR) const;

      //! @brief reduction sum kernel
      //! @param VECTOR the vector to get the sum of
      //! @return the reduced result sum
      t_DataType Sum( const Vector< t_DataType> &VECTOR) const;

      //! @brief reduction min kernel
      //! @param VECTOR the vector to get the min of
      //! @return the reduced resulting min
      t_DataType Min( const linal::VectorConstInterface< t_DataType> &VECTOR) const;

      //! @brief reduction max kernel
      //! @param VECTOR the vector to get the max of
      //! @return the reduced resulting max
      t_DataType Max( const linal::VectorConstInterface< t_DataType> &VECTOR) const;

      //! @brief gets the min and max for each column in a matrix
      //! @param MATRIX the matrix input
      //! @return a storage vector of math ranges with min and max for each column in the matrix
      storage::Vector< math::Range< t_DataType> > MinMax( const linal::MatrixConstInterface< t_DataType> &MATRIX) const;

      //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
      //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
      void VectorEqualsVectorTimesMatrix
      (
        linal::VectorInterface< t_DataType> &STORAGE,
        const linal::VectorConstInterface< t_DataType> &FEATURE,
        const linal::MatrixConstInterface< t_DataType> &MATRIX
      ) const;

      //! an optimized implementation of matrix * vector, with storage the storage (output) vector given externally
      //! STORAGE may optionally be the same as FEATURE.  Storage should be initialized externally
      void VectorPlusEqualsMatrixTimesVector
      (
        linal::VectorInterface< t_DataType> &STORAGE,
        const linal::MatrixConstInterface< t_DataType> &MATRIX,
        const linal::VectorConstInterface< t_DataType> &FEATURE
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief compile programs for given precision
      //! @param PRECISION float or double
      //! @return ERROR error that occured, CL_SUCCESS if no error
      cl_int CompilePrograms( const util::CPPDataTypes::Types &PRECISION);

    }; // template class Operations

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Operations< float>;
    BCL_EXPIMP_TEMPLATE template class BCL_API Operations< double>;

  } // namespace opencl
} // namespace bcl

#endif //BCL_OPENCL_OPERATIONS_H_
