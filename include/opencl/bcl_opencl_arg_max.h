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

#ifndef BCL_OPENCL_ARG_MAX_H_
#define BCL_OPENCL_ARG_MAX_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_kernel_sources.h"
#include "bcl_opencl_vector.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_limits.h"
#include "storage/bcl_storage_pair.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ArgMax
    //! @brief gets the max value in an array and returns it as well as its index. Optimized for the gpu
    //!
    //! @tparam t_DataType can be float, double, complex, int, etc...
    //!
    //! @see @link example_opencl_arg_max.cpp @endlink
    //! @author loweew
    //! @date Mar 25, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class ArgMax :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! command queue
      CommandQueue m_Queue;

      //! program
      cl::Program m_Program;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      ArgMax()
      {
      }

      //! @brief default constructor
      ArgMax( const CommandQueue &QUEUE) :
        m_Queue( QUEUE),
        m_Program()
      {
        cl_int error_number( CL_SUCCESS);
        m_Program = KernelSources::Compile( GetKernelSources().e_Linal, util::CPPDataTypes::DataTypeFromTemplate< t_DataType>(), m_Queue, std::string(), &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
      }

      //! @brief Clone function
      //! @return pointer to new ArgMax
      ArgMax< t_DataType> *Clone() const
      {
        return new ArgMax< t_DataType>( *this);
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

      //! @brief operator for performing arg max
      //! @param DATA data for which to find max and index of max
      //! @return the resulting max value and its index
      storage::Pair< t_DataType, int> operator()( Vector< t_DataType> &DATA) const
      {
        // error catching
        cl_int error_number = CL_SUCCESS;
        const cl_uint block_size( 256);

        const size_t nr_groups( ( DATA.GetSize() % block_size == 0 ? 0 : 1) + DATA.GetSize() / block_size);

        Vector< t_DataType> result_elements( nr_groups, m_Queue);
        Vector< int> result_indeces( nr_groups, m_Queue);
        linal::Vector< t_DataType> elements( nr_groups);
        linal::Vector< int> indeces( nr_groups);

        cl::Kernel kernel( m_Program, "ArgMax", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // set the args values
        error_number  = kernel.setArg( 0, DATA.GetData());
        error_number |= kernel.setArg( 1, result_indeces.GetData());
        error_number |= kernel.setArg( 2, result_elements.GetData());
        error_number |= kernel.setArg( 3, cl_uint( DATA.GetSize()));
        error_number |= kernel.setArg( 4, block_size * sizeof( t_DataType), NULL); //shared memory
        error_number |= kernel.setArg( 5, block_size * sizeof( int), NULL); //shared memory
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // Execute
        const cl::NDRange block_dims( block_size);
        const cl::NDRange offset;
        const cl::NDRange worksize( Tools::RoundUp( block_size, DATA.GetSize()));

        m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);

        elements = result_elements.GetHostVector();
        indeces = result_indeces.GetHostVector();

        // complete on cpu
        t_DataType final_element( math::GetLowestUnboundedValue< t_DataType>());
        size_t final_index( 0);
        for( size_t count( 0); count < nr_groups; ++count)
        {
          size_t greater_than( elements( count) > final_element ? 1 : 0);
          greater_than ? final_element = elements( count), final_index = indeces( count) : 0;
        }

        return storage::Pair< t_DataType, int>( final_element, final_index);
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

    }; // template class ArgMax

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_ARG_MAX_H_
