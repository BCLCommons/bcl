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

#ifndef BCL_OPENCL_SUPPORT_VECTOR_MACHINE_H_
#define BCL_OPENCL_SUPPORT_VECTOR_MACHINE_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_buffer.h"
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_context.h"
#include "bcl_opencl_kernel_source_interface.h"
#include "bcl_opencl_matrix.h"
#include "bcl_opencl_model_interface.h"
#include "bcl_opencl_vector.h"
#include "model/bcl_model_feature_data_set.h"
#include "model/bcl_model_interface.h"
#include "model/bcl_model_support_vector_kernel_base.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SupportVectorMachine
    //! @brief opencl implementation of support vector machine which can only be used for prediction and is optimized
    //!        for the gpu which uses float data type
    //!
    //! @see @link example_opencl_support_vector_machine.cpp @endlink
    //! @author loweew, woetzen
    //! @date Mar 12, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SupportVectorMachine :
      public ModelInterface
    {

    private:

    //////////
    // data //
    //////////

      //! bias
      float m_Bias;

      //! alphas
      linal::Vector< float> m_AlphasHost;
      Vector< float>       m_Alphas;

      //! support vectors
      model::FeatureDataSet< float> m_SupportVectorsHost;
      Matrix< float>       m_SupportVectors;
      size_t       m_NumberSupportVectors; //! number of support vectors in buffer including the padding

      //! Kernel that is used
      util::Implementation< model::SupportVectorKernelBase> m_Kernel;

      //! rescale in
      model::RescaleFeatureDataSet m_RescaleInput;

      //! rescale out
      model::RescaleFeatureDataSet m_RescaleOutput;

      //! queue
      CommandQueue m_Queue;

      //! program
      cl::Program m_Program;

      //! default size for a groups
      static const cl_uint s_GroupSize = 128;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      SupportVectorMachine();

      SupportVectorMachine( const CommandQueue &QUEUE);

      //! @brief default constructor
      //! @param BIAS the shift of the hyperplane from the origin
      //! @param ALPHAS vector of the alpha values
      //! @param SUPPORT_VECTORS the support vectors from the training
      //! @param KERNEL the kernel function to use
      //! @param RESCALE_IN, RESCALE_OUT the rescale functions
      SupportVectorMachine
      (
        const float BIAS,
        const linal::VectorConstInterface< float> &ALPHAS,
        const model::FeatureDataSetInterface< float> &SUPPORT_VECTORS,
        const util::Implementation< model::SupportVectorKernelBase> &KERNEL,
        const model::RescaleFeatureDataSet &RESCALE_IN,
        const model::RescaleFeatureDataSet &RESCALE_OUT,
        const CommandQueue &QUEUE
      );

      //! @brief Clone function
      //! @return pointer to new SupportVectorMachine
      SupportVectorMachine *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief method for getting m_Kernel
      //! @return the used SupportVectorKernelBase
      const util::Implementation< model::SupportVectorKernelBase> &GetKernel() const
      {
        return m_Kernel;
      }

      //! @brief set the alpha vector
      //! @param ALPHA the alphas
      void SetAlpha( const linal::VectorConstInterface< float> &ALPHA)
      {
        m_AlphasHost = ALPHA;
        m_Alphas = Vector< float>( ALPHA, m_Queue);
      }

      //! @brief set the bias
      //! @param BIAS the bias
      void SetBias( const float &BIAS)
      {
        m_Bias = BIAS;
      }

      //! @brief sets the support vectors
      //! @param SV the support vectors
      void SetSupportVectors( const model::FeatureDataSet< float> &SV)
      {
        m_SupportVectorsHost = SV;
        m_SupportVectors = Matrix< float>( SV.GetMatrix(), m_Queue);
      }

      //! @brief set number of sv
      //! @param NR_SV the number of support vectors
      void SetNumberSupportVectors( const size_t &NR_SV)
      {
        m_NumberSupportVectors = NR_SV;
      }

      //! @brief get number sv
      //! @return nr sv
      size_t GetNumberSupportVectors() const
      {
        return m_NumberSupportVectors;
      }

      //! @brief get the output feature size for this model
      //! @return the output feature size for this model
      size_t GetNumberOutputs() const
      {
        return size_t( 1);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Set the scaling of a feature set according to the model
      //! @param FEATURES feature set of interest
      //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
      //!       when operator() is called
      void Rescale( model::FeatureDataSet< float> &FEATURE) const;

      //! @brief predict result with model using a NOT rescaled feature vector
      //! @param FEATURE not rescaled feature vector
      //! @return predcited result vector using a model
      model::FeatureDataSet< float>
      PredictWithoutRescaling( const model::FeatureDataSetInterface< float> &FEATURE) const;

    ///////////////
    // operators //
    ///////////////

      //! @brief predict result with model using a rescaled feature vector
      //! @param FEATURE normalized or rescaled feature vector
      //! @return predcited result vector using a model
      model::FeatureDataSet< float>
      operator()( const model::FeatureDataSetInterface< float> &FEATURE) const;

      //! @brief predict based on opencl::Matrix data structure
      //! @param FEATURE the matrix of features
      //! @return the matrix of predicted outputs
      Matrix< float> operator()( const Matrix< float> &FEATURE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

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

    }; // class SupportVectorMachine

  } // namespace opencl
} // namespace bcl

#endif // BCL_OPENCL_SUPPORT_VECTOR_MACHINE_H_
