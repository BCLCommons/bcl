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

#ifndef BCL_OPENCL_EUCLIDEAN_DISTANCE_H_
#define BCL_OPENCL_EUCLIDEAN_DISTANCE_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_buffer.h"
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_kernel_sources.h"
#include "bcl_opencl_matrix.h"
#include "linal/bcl_linal_matrix.h"
#include "model/bcl_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class EuclideanDistance
    //! @brief calculates pairwise euclidean distance between rows of two matrices
    //!
    //! @see @link example_opencl_euclidean_distance.cpp @endlink
    //! @author loweew
    //! @date Mar 23, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class EuclideanDistance :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! opencl queue
      CommandQueue              m_Queue;

      //! opencl program
      cl::Program               m_Program;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from command queue
      EuclideanDistance()
      {
      }

      //! @brief constructor from command queue
      //! @param QUEUE command queue
      EuclideanDistance( const CommandQueue &QUEUE) :
        m_Queue( QUEUE)
      {
        cl_int error_number( CL_SUCCESS);
        m_Program = KernelSources::Compile( GetKernelSources().e_EuclideanDistance, util::CPPDataTypes::DataTypeFromTemplate< t_DataType>(), m_Queue, std::string(), &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
      }

      //! @brief Clone function
      //! @return pointer to new EuclideanDistance
      EuclideanDistance< t_DataType> *Clone() const
      {
        return new EuclideanDistance< t_DataType>( *this);
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

    ///////////////
    // operators //
    ///////////////

      //! @brief pairwise euclidean distance calculation
      //! @param INPUT_A, INPUT_B opencl matrices to compare
      //! @return matrix containing euclidean distances
      Matrix< t_DataType> operator()
      (
        const Matrix< t_DataType> &INPUT_A,
        const Matrix< t_DataType> &INPUT_B
      ) const
      {
        const cl_uint     block_size( 16);
        // catch errors
        cl_int error_number = CL_SUCCESS;
        BCL_Assert
        (
          INPUT_A.GetNumberRows() % block_size == 0 &&
          INPUT_A.GetNumberCols() % block_size == 0 &&
          INPUT_B.GetNumberRows() % block_size == 0 &&
          INPUT_B.GetNumberCols() % block_size == 0,
          "input buffers are not padded correctly!"
        );

        const cl_uint rows_a( INPUT_A.GetNumberRows());
        const cl_uint cols_a( INPUT_A.GetNumberCols());
        const cl_uint rows_b( INPUT_B.GetNumberRows());
        const cl_uint cols_b( INPUT_B.GetNumberCols());

        // rmsd group and worksize
        const cl::NDRange block_dim( block_size, block_size);
        const cl::NDRange worksize( Tools::RoundUp( block_size, rows_b), Tools::RoundUp( block_size, rows_a));
        const cl::NDRange offset;

        // allocating buffer
        Matrix< t_DataType> output_buffer( rows_a - INPUT_A.GetRowPadding(), rows_b - INPUT_B.GetRowPadding(), m_Queue, INPUT_A.GetRowPadding(), INPUT_B.GetRowPadding());

        // construct kernel
        cl::Kernel kernel( m_Program, "EuclideanDistance", &error_number);

        // set args for rmsd
        error_number  = kernel.setArg( 0, output_buffer.GetData());
        error_number |= kernel.setArg( 1, INPUT_A.GetData());
        error_number |= kernel.setArg( 2, INPUT_B.GetData());
        error_number |= kernel.setArg( 3, sizeof( t_DataType) * block_size * block_size, 0);
        error_number |= kernel.setArg( 4, sizeof( t_DataType) * block_size * block_size, 0);
        error_number |= kernel.setArg( 5, rows_a);
        error_number |= kernel.setArg( 6, cols_a);
        error_number |= kernel.setArg( 7, rows_b);
        error_number |= kernel.setArg( 8, cols_b);
        BCL_Assert( error_number == CL_SUCCESS, "euclidean distance arg error: " + opencl::Tools::ErrorString( error_number));

        // launch kernel
        error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dim, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

        return output_buffer;
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
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; // class EuclideanDistance

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_EUCLIDEAN_DISTANCE_H_ 
