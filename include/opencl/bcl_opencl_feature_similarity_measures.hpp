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

#ifndef BCL_OPENCL_FEATURE_SIMILARITY_MEASURES_HPP_
#define BCL_OPENCL_FEATURE_SIMILARITY_MEASURES_HPP_

// include header of this class
#include "bcl_opencl_feature_similarity_measures.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl.h"
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_kernel_source_file.h"
#include "bcl_opencl_kernel_source_interface.h"
#include "bcl_opencl_matrix.h"
#include "bcl_opencl_tools.h"
#include "bcl_opencl_vector.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_matrix.h"
#include "model/bcl_model_feature_similarity_interface.h"
#include "util/bcl_util_binary_function_interface.h"
#include "util/bcl_util_object_data_label.h"

// external includes - sorted alphabetically
#include <CL/cl.hpp>

namespace bcl
{
  namespace opencl
  {
    // instantiate s_Instance
    template< typename t_DataType>
    const util::SiPtr< const util::ObjectInterface> FeatureSimilarityMeasures< t_DataType>::s_Instance
    (
      util::Enumerated< model::FeatureSimilarityMeasuresInterface< t_DataType> >::AddInstance( new FeatureSimilarityMeasures())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor
    template< typename t_DataType>
    FeatureSimilarityMeasures< t_DataType>::FeatureSimilarityMeasures()
    {
    }

    //! @brief constructor
    //! @param MEASURE the measurement
    //! @param QUEUE the command queue
    template< typename t_DataType>
    FeatureSimilarityMeasures< t_DataType>::FeatureSimilarityMeasures( std::string &MEASURE, const CommandQueue &QUEUE) :
      m_Measure( MEASURE),
      m_Queue( QUEUE)
    {
      std::string kernel_name( "feature_similarity_measures.cl");
      cl_int error_number = GetTools().CompilePrograms
          (
            util::CPPDataTypes::DataTypeFromTemplate< t_DataType>(),
            m_Queue,
            m_Program,
            kernel_name
          );
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
    }

    //! @brief Clone function
    //! @return pointer to new FeatureSimilarityMeasures< t_DataType>
    template< typename t_DataType>
    FeatureSimilarityMeasures< t_DataType> *FeatureSimilarityMeasures< t_DataType>::Clone() const
    {
      return new FeatureSimilarityMeasures< t_DataType>( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &FeatureSimilarityMeasures< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    template< typename t_DataType>
    const std::string &FeatureSimilarityMeasures< t_DataType>::GetAlias() const
    {
      static const std::string s_Name( "OpenCLSimilarityMeasures");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the similarity measure
    //! @param INPUT_A the input matrix for which to compare against the matrix from the constructor
    //! @return the Matrix of the resulting similarity measures
    template< typename t_DataType>
    linal::Matrix< t_DataType> FeatureSimilarityMeasures< t_DataType>::operator()
    (
      const linal::MatrixConstInterface< t_DataType> &INPUT_A
    ) const
    {
      // calculate padding
      const size_t row_pad_a( Tools::CalcPadding( s_blocksize, INPUT_A.GetNumberRows()));
      const size_t col_pad_a( Tools::CalcPadding( s_blocksize, INPUT_A.GetNumberCols()));

      // create device input matrices
      const Matrix< t_DataType> d_input_a( INPUT_A, m_Queue, row_pad_a, col_pad_a);

      BCL_MessageStd( "input_a rows: " + util::Format()( d_input_a.GetNumberRows()) + "\ninput_a cols: " + util::Format()( d_input_a.GetNumberCols()));

      // perform calculation
      Matrix< t_DataType> output;
      const size_t s_size_limit( 15000);
      linal::Matrix< t_DataType> result( INPUT_A.GetNumberRows(), INPUT_A.GetNumberRows());
      if( INPUT_A.GetNumberRows() > s_size_limit)
      {
        result = SetupAndLaunchKernelSlow( d_input_a, m_Measure);
      }
      else
      {
        output = SetupAndLaunchKernelFast( d_input_a, m_Measure);
        result = output.GetHostMatrix( row_pad_a, row_pad_a);
      }
      return result;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    template< typename t_DataType>
    std::istream &FeatureSimilarityMeasures< t_DataType>::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    template< typename t_DataType>
    std::ostream &FeatureSimilarityMeasures< t_DataType>::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer FeatureSimilarityMeasures< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "performs pairwise similarity coefficients between two inputs");

      parameters.AddInitializer
      (
        "measurement",
        "the measurement type i.e. Tanimoto, Cosine, Dice, Euclidean, Manhattan",
        io::Serialization::GetAgentWithCheck
        (
          &m_Measure,
          command::ParameterCheckAllowed
          (
            storage::Vector< std::string>::Create( "Tanimoto", "Dice", "Cosine", "Euclidean", "Manhattan", "RMSD")
          )
        ),
        "Tanimoto"
      );
      return parameters;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief sets up the kernel and launches it
    //! @param INPUT_A the input matrix
    //! @param INPUT_B the other input matrix
    //! @param KERNEL_NAME the name of the kernel as a string
    //! @return the error number for checking if launch is successful
    template< typename t_DataType>
    Matrix< t_DataType> FeatureSimilarityMeasures< t_DataType>::SetupAndLaunchKernelFast
    (
      const Matrix< t_DataType> &INPUT_A,
      const std::string &KERNEL_NAME
    ) const
    {
      // error catching
      cl_int error_number = CL_SUCCESS;
      const cl_uint rows_a( INPUT_A.GetNumberRows());
      const cl_uint row_pad_a( INPUT_A.GetRowPadding());
      const cl_uint cols_a( INPUT_A.GetNumberCols());

      Matrix< t_DataType> output
      (
        rows_a - row_pad_a,
        rows_a - row_pad_a,
        m_Queue,
        row_pad_a,
        row_pad_a
      );
      std::string rmsd( "RMSD");

      // construct kernels
      cl::Kernel kernel( m_Program, KERNEL_NAME.c_str(), &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      error_number =  kernel.setArg( 0, output.GetData());
      error_number |= kernel.setArg( 1, INPUT_A.GetData());
      error_number |= kernel.setArg( 2, s_blocksize * s_blocksize * sizeof( t_DataType), 0);
      error_number |= kernel.setArg( 3, s_blocksize * s_blocksize * sizeof( t_DataType), 0);
      error_number |= kernel.setArg( 4, rows_a);
      error_number |= kernel.setArg( 5, cols_a);
      if( KERNEL_NAME == rmsd)
      {
        const cl_uint nr_atoms( ( cols_a - INPUT_A.GetColPadding()) / 3);
        error_number |= kernel.setArg( 6, nr_atoms);
      }
      BCL_Assert( error_number == CL_SUCCESS, "kernel error: " + opencl::Tools::ErrorString( error_number));

      cl_uint range_x( Tools::RoundUp( s_blocksize, rows_a));
      cl_uint range_y( Tools::RoundUp( s_blocksize, rows_a));

      // thread launch dimensions
      const cl::NDRange block_dim( s_blocksize, s_blocksize);
      const cl::NDRange offset;
      const cl::NDRange worksize( range_x, range_y);

      // launch kernel
      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dim, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      return output;
    }

    //! @brief sets up the kernel and launches it
    //! @param INPUT_A the input matrix
    //! @param INPUT_B the other input matrix
    //! @param KERNEL_NAME the name of the kernel as a string
    //! @return the error number for checking if launch is successful
    template< typename t_DataType>
    linal::Matrix< t_DataType> FeatureSimilarityMeasures< t_DataType>::SetupAndLaunchKernelSlow
    (
      const Matrix< t_DataType> &INPUT_A,
      const std::string &KERNEL_NAME
    ) const
    {
      // error catching
      cl_int error_number = CL_SUCCESS;
      const cl_uint rows_a( INPUT_A.GetNumberRows());
      const cl_uint row_pad_a( INPUT_A.GetRowPadding());
      const cl_uint cols_a( INPUT_A.GetNumberCols());

      linal::Matrix< t_DataType> result( rows_a, rows_a);
      linal::Vector< t_DataType> tmp( size_t( rows_a), 0.0);
      Vector< t_DataType> output
      (
        rows_a - row_pad_a,
        m_Queue,
        row_pad_a
      );

      std::string orig( KERNEL_NAME);
      std::string ext( "Slow");
      orig.append( ext);

      // construct kernels
      cl::Kernel kernel( m_Program, orig.c_str(), &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      for( cl_uint row( 0); row < rows_a; ++row)
      {
        error_number =  kernel.setArg( 0, output.GetData());
        error_number |= kernel.setArg( 1, row);
        error_number |= kernel.setArg( 2, INPUT_A.GetData());
        error_number |= kernel.setArg( 3, s_blocksize * s_blocksize * sizeof( t_DataType), 0);
        error_number |= kernel.setArg( 4, rows_a);
        error_number |= kernel.setArg( 5, cols_a);
        BCL_Assert( error_number == CL_SUCCESS, "kernel error: " + opencl::Tools::ErrorString( error_number));

        cl_uint range_x( Tools::RoundUp( s_blocksize * s_blocksize, rows_a));

        // thread launch dimensions
        const cl::NDRange block_dim( s_blocksize * s_blocksize);
        const cl::NDRange offset;
        const cl::NDRange worksize( range_x);

        // launch kernel
        error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dim, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
        tmp = output.GetHostVector();
        std::copy( tmp.Begin(), tmp.End(), result.Begin() + ( row * rows_a));
      }

      return result.CreateSubMatrix( rows_a - row_pad_a, rows_a - row_pad_a);
    }

    //! @brief responsible for updating to a valid queue
    //! @param TOOLS opencl tools
    template< typename t_DataType>
    void FeatureSimilarityMeasures< t_DataType>::UpdateQueue( Tools &TOOLS)
    {
      if( !TOOLS.HasCommandQueues())
      {
        return;
      }
      std::string kernel_name( "feature_similarity_measures.cl");

      m_Queue = GetTools().GetFirstCommandQueue();
      cl_int error_number = GetTools().CompilePrograms
          (
            util::CPPDataTypes::DataTypeFromTemplate< t_DataType>(),
            m_Queue,
            m_Program,
            kernel_name
          );
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
    }

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_FEATURE_SIMILARITY_MEASURES_HPP_
