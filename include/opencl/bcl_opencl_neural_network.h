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

#ifndef BCL_OPENCL_NEURAL_NETWORK_H_
#define BCL_OPENCL_NEURAL_NETWORK_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_opencl_kernel_source_interface.h"
#include "bcl_opencl_matrix.h"
#include "bcl_opencl_matrix_multiply.h"
#include "bcl_opencl_matrix_transpose.h"
#include "bcl_opencl_model_interface.h"
#include "bcl_opencl_transfer_function_sigmoid.h"
#include "bcl_opencl_vector.h"
#include "bcl_opencl_vector_matrix_add.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_vector.h"
#include "math/bcl_math_range.h"
#include "model/bcl_model.h"
#include "model/bcl_model_feature_data_set.h"
#include "model/bcl_model_transfer_function_interface.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_object_data_label.h"

// external includes - sorted alphabetically
#include <iosfwd>

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class NeuralNetwork
    //! @brief is an artificial neural network in opencl optimized for the gpu
    //!
    //! @see @link example_opencl_neural_network.cpp @endlink
    //! @author loweew
    //! @date Mar 28, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API NeuralNetwork :
      public ModelInterface
    {
    public:

    /////////////
    // friends //
    /////////////

    //////////
    // data //
    //////////

      //! the default input range for neural network transfer functions
      static const math::Range< float> s_DefaultInputRange;

      //! the default input range for neural network transfer functions
      static const math::Range< float> s_DefaultOutputRange;

    private:

    //////////
    // data //
    //////////

      //! transfer function used in every neuron
      util::Implementation< model::TransferFunctionInterface> m_TransferFunction;

      //! rescale input to transfer function input range
      util::ShPtr< model::RescaleFeatureDataSet>       m_RescaleInput;

      //! rescale output from transfer function output range
      util::ShPtr< model::RescaleFeatureDataSet>       m_RescaleOutput;

      //! bias
      storage::Vector< linal::Vector< float> >         m_Bias;

      //! weight
      storage::Vector< linal::Matrix< float> >         m_Weight;

      //! architecture of the neural network, stored for convenience
      storage::Vector< size_t>                         m_Architecture;

      //! device bias buffers
      mutable storage::Vector< Vector< float> >        m_BiasBuffers;

      //! device weight buffers
      mutable storage::Vector< Matrix< float> >        m_WeightBuffers;

      //! tranposed device weight buffers
      mutable storage::Vector< Matrix< float> >        m_TransWeightBuffers;

      //! opencl queue
      CommandQueue                                     m_Queue;

      //! opencl program
      cl::Program                                      m_Program;

      //! gpu mmult
      MatrixMultiply< float>                           m_GpuMMult;

      //! gpu matrix transpose
      MatrixTranspose< float>                          m_GpuTranspose;

      //! gpu v-m add
      VectorMatrixAdd< float>                          m_GpuVMAdd;

      //! gpu sigmoid
      TransferFunctionSigmoid< float>                  m_GpuSigmoid;

      //! bool for whether or not the device data is set
      mutable bool                                     m_DeviceDataBool;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////
    // data //
    //////////

      //! block size
      static const cl_uint s_blocksize;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from queue
      NeuralNetwork();

      //! @brief constructor from queue
      NeuralNetwork( const CommandQueue &QUEUE);

      //! @brief construct from all necessary parameters
      //! @param RESCALE_INPUT function to rescale input to transfer function range
      //! @param RESCALE_OUTPUT function to rescale output from transfer function range
      //! @param BIAS the bias for each neuron of a hidden layer
      //! @param WEIGHT the weight for each neuron of a hidden layer
      //! @param TRANSFER_FUNCTION function to simulate the threshold in the neuron
      NeuralNetwork
      (
        const util::ShPtr< model::RescaleFeatureDataSet> &RESCALE_INPUT,
        const util::ShPtr< model::RescaleFeatureDataSet> &RESCALE_OUTPUT,
        const storage::Vector< linal::Vector< float> > &BIAS,
        const storage::Vector< linal::Matrix< float> > &WEIGHT,
        const util::Implementation< model::TransferFunctionInterface> &TRANSFER_FUNCTION,
        const CommandQueue &QUEUE
      );

      //! @brief construct from all necessary parameters for the gpu prediction
      //!        this assumes that the inputs are already padded appropriately and rescaled
      //! @param BIAS the bias for each neuron of a hidden layer
      //! @param WEIGHT the weight for each neuron of a hidden layer
      NeuralNetwork
      (
        const util::ShPtr< model::RescaleFeatureDataSet> &RESCALE_INPUT,
        const util::ShPtr< model::RescaleFeatureDataSet> &RESCALE_OUTPUT,
        const storage::Vector< Vector< float> > &BIAS,
        const storage::Vector< Matrix< float> > &WEIGHT,
        const CommandQueue &QUEUE
      );

      //! copy constructor
      NeuralNetwork *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief get number layers
      size_t GetNumberLayers() const
      {
        return m_Weight.GetSize();
      }

      //! @brief returns architecture of this neural network
      //! @return first entry is the number of input neurons
      //!         middle entries are the number of hidden neurons in each layer
      //!         last number in vector is the number of output neurons
      const storage::Vector< size_t> &GetArchitecture() const;

      //! @brief set the architecture architecture of this neural network
      //! @param ARCHITECTURE first entry is the number of input neurons
      //!         middle entries are the number of hidden neurons in each layer
      //!         last number in vector is the number of output neurons
      void SetArchitecture( const storage::Vector< size_t> &ARCHITECTURE);

      //! get number inputs
      size_t GetNumberInputs() const;

      //! get number output neurons
      size_t GetNumberOutputs() const;

      //! @brief get the bias
      const storage::Vector< linal::Vector< float> > &GetBias() const
      {
        return m_Bias;
      }

      //! @brief set the bias
      void SetBias( const storage::Vector< linal::Vector< float> > &BIAS)
      {
        m_Bias = BIAS;
      }

      //! @brief get the weight
      const storage::Vector< linal::Matrix< float> > &GetWeight() const
      {
        return m_Weight;
      }

      //! @brief set the weight
      void SetWeight( const storage::Vector< linal::Matrix< float> > &WEIGHT)
      {
        m_Weight = WEIGHT;
      }

      //! @brief access the transfer function
      const util::Implementation< model::TransferFunctionInterface> &GetTransferFunction() const
      {
        return m_TransferFunction;
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
      //! @return predicted result vector using a model
      model::FeatureDataSet< float> PredictWithoutRescaling
      (
        const model::FeatureDataSetInterface< float> &FEATURE
      ) const;

      //! @brief predict result with model using a rescaled feature vector
      //! @param FEATURE normalized or rescaled feature vector
      //! @return predicted result vector using a model
      model::FeatureDataSet< float> operator()( const model::FeatureDataSetInterface< float> &FEATURE) const;

      //! @brief predict with gpu without rescaling
      Matrix< float> operator()( const Matrix< float> &INPUT) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read NeuralNetwork from std::istream
      std::istream &Read( std::istream &ISTREAM);

      //! write NeuralNetwork into std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief determines the architecture already present in the weight matrix/vectors of the neural network
      void SetImplicitArchitecture();

      //! @brief sigmoid transfer function to all elements in a matrix
      //! @param INPUT the input matrix
      void Sigmoid
      (
        Matrix< float> &INPUT,
        Matrix< float> &OUTPUT
      );

      //! @brief adds single row to all rows of a matrix
      //! @param BIAS the bias to add
      //! @param MATRIX the data to add the bias to
      void AddBias
      (
        Vector< float> &BIAS,
        Matrix< float> &MATRIX
      );

      //! @brief sets up everything necessary for gpu computing
      void EnableGpu();

      //! @brief sets device data
      void SetDeviceData() const;

      //! @brief compile programs for given precision
      //! @param PRECISION float or double
      //! @return ERROR error that occured, CL_SUCCESS if no error
      cl_int CompilePrograms( const util::CPPDataTypes::Types &PRECISION);

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      //! @return result of any validation performed internally
      io::ValidationResult PreReadHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
      {
        UpdateQueue( GetTools());
        return io::ValidationResult( true);
      }

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
      {
        SetArchitecture( m_Architecture);
        return true;
      }

      //! @brief responsible for updating to a valid queue
      //! @param TOOLS opencl tools
      void UpdateQueue( Tools &TOOLS);

    }; // class NeuralNetwork

  } // namespace opencl
} // namespace bcl

#endif //BCL_OPENCL_NEURAL_NETWORK_H_
