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

#ifndef BCL_OPENCL_FEATURE_SIMILARITY_MEASURES_H_
#define BCL_OPENCL_FEATURE_SIMILARITY_MEASURES_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_kernel_source_file.h"
#include "bcl_opencl_kernel_source_interface.h"
#include "bcl_opencl_matrix.h"
#include "bcl_opencl_tools.h"
#include "io/bcl_io_serialization.h"
#include "io/bcl_io_serialization_interface.h"
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
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FeatureSimilarityMeasures
    //! @brief provides several feature similarity measures that work on matrices of feature vectors
    //! @details
    //!
    //! @tparam t_DataType can be float, double, int
    //!
    //! @see @link example_opencl_feature_similarity_measures.cpp @endlink
    //! @author loweew
    //! @date Oct 03, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_DataType>
    class FeatureSimilarityMeasures :
      public model::FeatureSimilarityMeasuresInterface< t_DataType>
    {

    private:

    //////////
    // data //
    //////////

      std::string m_Measure; // enum

      //! command queue
      CommandQueue m_Queue;

      //! program
      cl::Program m_Program;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
      static const cl_uint s_blocksize = 16;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    public:

      //! @brief constructor
      FeatureSimilarityMeasures();

      //! @brief constructor
      //! @param MEASURE the measurement
      //! @param QUEUE the command queue
      FeatureSimilarityMeasures( std::string &MEASURE, const CommandQueue &QUEUE);

      //! @brief Clone function
      //! @return pointer to new FeatureSimilarityMeasures< t_DataType>
      FeatureSimilarityMeasures< t_DataType> *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate the similarity measure
      //! @param INPUT_A the input matrix for which to compare against the matrix from the constructor
      //! @return the Matrix of the resulting similarity measures
      linal::Matrix< t_DataType> operator()
      (
        const linal::MatrixConstInterface< t_DataType> &INPUT_A
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
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief sets up the kernel and launches it
      //! @param INPUT_A the input matrix
      //! @param KERNEL_NAME the name of the kernel as a string
      //! @return the error number for checking if launch is successful
      Matrix< t_DataType> SetupAndLaunchKernelFast
      (
        const Matrix< t_DataType> &INPUT_A,
        const std::string &KERNEL_NAME
      ) const;

      //! @brief sets up the kernel and launches it
      //! @param INPUT_A the input matrix
      //! @param KERNEL_NAME the name of the kernel as a string
      //! @return the error number for checking if launch is successful
      linal::Matrix< t_DataType> SetupAndLaunchKernelSlow
      (
        const Matrix< t_DataType> &INPUT_A,
        const std::string &KERNEL_NAME
      ) const;

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      //! @return result of any validation performed internally
      io::ValidationResult PreReadHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
      {
        UpdateQueue( GetTools());
        return io::ValidationResult( true);
      }

      //! @brief responsible for updating to a valid queue
      //! @param TOOLS opencl tools
      void UpdateQueue( Tools &TOOLS);

    }; // template class FeatureSimilarityMeasures

  ////////////////////////////
  // explicit instantiation //
  ////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API FeatureSimilarityMeasures< double>;
    BCL_EXPIMP_TEMPLATE template class BCL_API FeatureSimilarityMeasures< float>;

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_FEATURE_SIMILARITY_MEASURES_H_
