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

#ifndef BCL_OPENCL_COORDINATES_TRANSFORMER_H_
#define BCL_OPENCL_COORDINATES_TRANSFORMER_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_buffer.h"
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_kernel_source_file.h"
#include "bcl_opencl_kernel_source_interface.h"
#include "bcl_opencl_matrix.h"
#include "bcl_opencl_tools.h"
#include "math/bcl_math_transformation_matrix_3d.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CoordinateTransformer
    //! @brief given a buffer of coordinates and a transformation matrix, transforms the coordinates
    //! @details TODO: add an detailed description to this class
    //!
    //! @tparam t_DataType TODO: add a description for the template parameter
    //!
    //! @see @link example_opencl_coordinate_transformer.cpp @endlink
    //! @author woetzen
    //! @date Dec 5, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class CoordinateTransformer :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      CommandQueue m_CommandQueue; //!< command queue
      cl::Program  m_Program;      //!< OpenCL program that contains kernels for each precision

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CoordinateTransformer()
      {
      }

      //! @brief Clone function
      //! @return pointer to new CoordinateTransformer< t_DataType>
      CoordinateTransformer< t_DataType> *Clone() const
      {
        return new CoordinateTransformer< t_DataType>( *this);
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

      //! @brief initialize this class
      //! @brief COMMAND_QUEUE queue to use
      //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
      bool Initialize( const CommandQueue &COMMAND_QUEUE)
      {
        // update the command queue
        m_CommandQueue = COMMAND_QUEUE;

        // for precision type
        cl_int error_number = CompilePrograms( util::CPPDataTypes::DataTypeFromTemplate< t_DataType>());
        if( error_number != CL_SUCCESS)
        {
          BCL_MessageDbg( "error compiling programs:\n" + Tools::ErrorString( error_number));
          return false;
        }

        return true;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief transform coordinates in a given buffer
      //! @param COORDINATES buffer with matrix of 4 cols and multiple of 64 rows
      //! @param TRANSFORMATION_MATRIX the matrix to apply
      //! @return Buffer a new buffer with the transformed coordiantes
      Matrix< t_DataType> operator()( const Matrix< t_DataType> &COORDINATES, const math::TransformationMatrix3D &TRANSFORMATION_MATRIX, const size_t NR_COORDINATES) const
      {
        cl_int error_number( CL_SUCCESS);

        // create new atom buffer
        Matrix< t_DataType> new_atoms( COORDINATES.GetNumberRows(), COORDINATES.GetNumberCols(), m_CommandQueue);

        // enqueue kernel
        cl::Kernel kernel( m_Program, "CoordinateTransformation", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        const t_DataType *transformation_ptr( TRANSFORMATION_MATRIX.GetMatrix().Begin());
        error_number  = kernel.setArg(  0, COORDINATES.GetData());
        error_number |= kernel.setArg(  1, ( *transformation_ptr)); // a1
        ++transformation_ptr;
        error_number |= kernel.setArg(  2, ( *transformation_ptr)); // a2
        ++transformation_ptr;
        error_number |= kernel.setArg(  3, ( *transformation_ptr)); // a3
        ++transformation_ptr;
        ++transformation_ptr;
        error_number |= kernel.setArg(  4, ( *transformation_ptr)); // b1
        ++transformation_ptr;
        error_number |= kernel.setArg(  5, ( *transformation_ptr)); // b2
        ++transformation_ptr;
        error_number |= kernel.setArg(  6, ( *transformation_ptr)); // b3
        ++transformation_ptr;
        ++transformation_ptr;
        error_number |= kernel.setArg(  7, ( *transformation_ptr)); // c1
        ++transformation_ptr;
        error_number |= kernel.setArg(  8, ( *transformation_ptr)); // c2
        ++transformation_ptr;
        error_number |= kernel.setArg(  9, ( *transformation_ptr)); // c3
        ++transformation_ptr;
        ++transformation_ptr;
        error_number |= kernel.setArg( 10, ( *transformation_ptr)); // t1
        ++transformation_ptr;
        error_number |= kernel.setArg( 11, ( *transformation_ptr)); // t2
        ++transformation_ptr;
        error_number |= kernel.setArg( 12, ( *transformation_ptr)); // t3
        error_number |= kernel.setArg( 13, new_atoms.GetData());
        BCL_Assert( error_number == CL_SUCCESS, "setting gpu kernel args error: " + Tools::ErrorString( error_number));

        const cl::NDRange local_worksize( 64);
        const cl::NDRange offset;
        cl::NDRange global_worksize( NR_COORDINATES);

        // launching kernel
        error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, offset, global_worksize, local_worksize);
        BCL_Assert( error_number == CL_SUCCESS, "kernel launch error: " + Tools::ErrorString( error_number));

        // end
        return new_atoms;
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

      //! @brief compile programs for given precision
      //! @param PRECISION float or double
      //! @return ERROR error that occured, CL_SUCCESS if no error
      cl_int CompilePrograms( const util::CPPDataTypes::Types &PRECISION)
      {
        cl_int error_number( CL_SUCCESS);

        // compile the program
        cl::Program::Sources source;
        KernelSourceFile transform_source_file( "atom_transformation.cl");

        const Device device( m_CommandQueue.GetDevice( &error_number));
        BCL_Assert( error_number == CL_SUCCESS, "cannot get device for command queue");

        const std::string transform_source( transform_source_file.GetSource( PRECISION, device.Extensions()));
        if( transform_source.empty())
        {
          return CL_INVALID_KERNEL_DEFINITION;
        }
        source.push_back( std::make_pair( transform_source.c_str(), transform_source.length()));

        const Context context( m_CommandQueue.GetContext( &error_number));
        BCL_Assert( error_number == CL_SUCCESS, "cannot get context from command queue");

        // create the program
        cl::Program &current_program( m_Program);
        current_program = cl::Program( context, source, &error_number);
        if( error_number != CL_SUCCESS)
        {
          BCL_MessageCrt( "create program error: " + Tools::ErrorString( error_number));
          return error_number;
        }

        // build the program
        error_number = current_program.build( std::vector< cl::Device>( 1, device));
        if( error_number != CL_SUCCESS)
        {
          BCL_MessageCrt( "build program error: " + Tools::ErrorString( error_number));
          std::string build_info;
          error_number = current_program.getBuildInfo( device, CL_PROGRAM_BUILD_LOG, &build_info);
          if( error_number != CL_SUCCESS)
          {
            BCL_MessageCrt( "get build info error: " + Tools::ErrorString( error_number));
          }
          else
          {
            BCL_MessageCrt( "build log: " + build_info);
          }
          return error_number;
        }

        // end
        return error_number;
      }

    }; // template class CoordinateTransformer

    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> CoordinateTransformer< t_DataType>::s_Instance
    (
      GetObjectInstances().AddInstance( new CoordinateTransformer< t_DataType>())
    );

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_COORDINATES_TRANSFORMER_H_ 
