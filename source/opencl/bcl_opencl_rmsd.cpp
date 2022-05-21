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
#include "opencl/bcl_opencl_rmsd.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math.h"
#include "opencl/bcl_opencl_kernel_sources.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from command queue
    RMSD::RMSD()
    {
    }

    //! @brief constructor from command queue
    //! @param QUEUE command queue
    RMSD::RMSD( const CommandQueue &QUEUE) :
      m_Queue( QUEUE)
    {
      cl_int error_number( CL_SUCCESS);
      m_Program = KernelSources::Compile( GetKernelSources().e_RMSD, util::CPPDataTypes::e_Float, m_Queue, std::string(), &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
    }

    //! @brief Clone function
    //! @return pointer to new RMSD
    RMSD *RMSD::Clone() const
    {
      return new RMSD( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RMSD::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief rmsd calculation
    //! @param INPUT_A, INPUT_B cl buffers of matrices to compare
    //! @return rmsd of the INPUT_A vs INPUT_B
    float RMSD::operator()
    (
      const Matrix< float> &INPUT_A,
      const Matrix< float> &INPUT_B
    ) const
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint rows( INPUT_A.GetNumberRows());
      const cl_uint cols( INPUT_A.GetNumberCols());
      const cl_uint row_pad( INPUT_A.GetRowPadding());
      const cl_uint col_pad( INPUT_A.GetColPadding());
      const cl_uint rows_b( INPUT_B.GetNumberRows());
      const cl_uint cols_b( INPUT_B.GetNumberCols());

      BCL_Assert( rows == rows_b && cols == cols_b, "dimensions don't match!");

      // rmsd group and worksize
      const size_t      block_size( 256);
      const cl::NDRange block_dim ( block_size);
      const cl::NDRange worksize  ( Tools::RoundUp( block_size, rows * cols));
      const cl::NDRange offset;
      const size_t      num_groups((( rows * cols) % block_size == 0 ? 0 : 1 ) + (( rows * cols) / block_size));

      // construct kernel
      cl::Kernel kernel( m_Program, "RmsdKernel", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      Vector< float> vector( num_groups, m_Queue);

      // set args for rmsd
      error_number  = kernel.setArg( 0, INPUT_A.GetData());
      error_number |= kernel.setArg( 1, INPUT_B.GetData());
      error_number |= kernel.setArg( 2, vector.GetData());
      error_number |= kernel.setArg( 3, rows * cols);
      error_number |= kernel.setArg( 4, sizeof( float) * cl_uint( block_size), 0);
      BCL_Assert( error_number == CL_SUCCESS, "rmsd arg error: " + opencl::Tools::ErrorString( error_number));

      // launch kernel
      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dim, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // read back partially reduced sum of size num_groups
      linal::Vector< float> rmsd( vector.GetHostVector());

      BCL_MessageDbg( "printing rmsd partial sum: " + util::Format()( rmsd));

      return math::Sqrt( rmsd.Sum() / ( ( rows - row_pad) * ( cols - col_pad)));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &RMSD::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &RMSD::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace opencl
} // namespace bcl
