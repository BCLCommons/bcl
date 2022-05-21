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
#include "opencl/bcl_opencl_insertion_sort.h"

// includes from bcl - sorted alphabetically
#include "opencl/bcl_opencl_kernel_sources.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////
  // data //
  //////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from command queue
    InsertionSort::InsertionSort()
    {
    }

    //! @brief constructor from command queue
    //! @param QUEUE command queue
    InsertionSort::InsertionSort( const CommandQueue &QUEUE) :
      m_Queue( QUEUE)
    {
      cl_int error_number( CL_SUCCESS);
      m_Program = KernelSources::Compile( GetKernelSources().e_InsertionSort, util::CPPDataTypes::e_Float, m_Queue, std::string(), &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
    }

    //! @brief Clone function
    //! @return pointer to new InsertionSort
    InsertionSort *InsertionSort::Clone() const
    {
      return new InsertionSort( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &InsertionSort::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief sorts k values of columns in matrix into first k positions
    //!        also provides buffer of keys as "index buffer" through m_IndexMatrixOnDevice variable which can be retrieved
    //!        using Get function
    //! XXX NOTE: doing a partial sort over writes the values in the first k rows of the column instead of pushing them down the col
    //! @param DATA matrix buffer to sort
    //! @param NR_TO_SORT the k number of values you want sorted into the first k cols
    //! @return buffer with first k lowest values sorted in first k columns
    Matrix< float> InsertionSort::operator()
    (
      Matrix< float> &DATA, const size_t &NR_TO_SORT
    ) const
    {
      // error catching
      cl_int error_number = CL_SUCCESS;

      const cl_uint rows( DATA.GetNumberRows());
      const cl_uint cols( DATA.GetNumberCols());
      const cl_uint row_pad( DATA.GetRowPadding());
      const cl_uint col_pad( DATA.GetColPadding());

      // launch kernel
      cl::Kernel kernel( m_Program, "InsertionSort", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      m_IndexMatrixOnDevice = Matrix< int>( NR_TO_SORT, cols, m_Queue);

      error_number  = kernel.setArg(  0, DATA.GetData());
      error_number  = kernel.setArg(  1, cl_uint( cols));
      error_number |= kernel.setArg(  2, m_IndexMatrixOnDevice.GetData());
      error_number |= kernel.setArg(  3, cl_uint( cols));
      error_number |= kernel.setArg(  4, cl_uint( cols - col_pad));
      error_number |= kernel.setArg(  5, cl_uint( rows - row_pad));
      error_number |= kernel.setArg(  6, cl_uint( NR_TO_SORT));
      BCL_Assert( error_number == CL_SUCCESS, "setting kernel args error: " + Tools::ErrorString( error_number));

      const cl_uint block_size( 256);
      cl::NDRange offset;
      cl::NDRange block_dims( block_size);
      cl::NDRange worksize( Tools::RoundUp( block_size, cols));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
      BCL_Assert( error_number == CL_SUCCESS, "error launching kernel: " + Tools::ErrorString( error_number));

      m_Queue.finish();

      return DATA;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &InsertionSort::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &InsertionSort::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace opencl
} // namespace bcl
