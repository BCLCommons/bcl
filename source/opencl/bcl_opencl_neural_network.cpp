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
#include "opencl/bcl_opencl_neural_network.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "model/bcl_model_transfer_sigmoid.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////
  // data //
  //////////

    //! the default input range for neural network transfer functions
    const math::Range< float> NeuralNetwork::s_DefaultInputRange( -1, 1);

    //! the default input range for neural network transfer functions
    const math::Range< float> NeuralNetwork::s_DefaultOutputRange( 0.1, 0.9);

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> NeuralNetwork::s_Instance
    (
      GetObjectInstances().AddInstance( new NeuralNetwork())
    );

    const cl_uint NeuralNetwork::s_blocksize = 16;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from queue
    NeuralNetwork::NeuralNetwork() :
      m_TransferFunction(),
      m_RescaleInput(),
      m_RescaleOutput(),
      m_Bias(),
      m_Weight(),
      m_Queue(),
      m_DeviceDataBool( false)
    {
    }

    //! @brief constructor from queue
    NeuralNetwork::NeuralNetwork( const CommandQueue &QUEUE) :
      m_TransferFunction(),
      m_RescaleInput(),
      m_RescaleOutput(),
      m_Bias(),
      m_Weight(),
      m_Queue( QUEUE),
      m_DeviceDataBool( false)
    {
    }

    //! @brief construct from all necessary parameters
    //! @param RESCALE_INPUT
    //! @param RESCALE_OUTPUT
    //! @param BIAS
    //! @param WEIGHT
    //! @param TRANSFER_FUNCTION
    NeuralNetwork::NeuralNetwork
    (
      const util::ShPtr< model:: RescaleFeatureDataSet> &RESCALE_INPUT,
      const util::ShPtr< model:: RescaleFeatureDataSet> &RESCALE_OUTPUT,
      const storage::Vector< linal::Vector< float> > &BIAS,
      const storage::Vector< linal::Matrix< float> > &WEIGHT,
      const util::Implementation< model::TransferFunctionInterface> &TRANSFER_FUNCTION,
      const CommandQueue &QUEUE
    ) :
      m_TransferFunction( TRANSFER_FUNCTION),
      m_RescaleInput( RESCALE_INPUT),
      m_RescaleOutput( RESCALE_OUTPUT),
      m_Bias( BIAS),
      m_Weight( WEIGHT),
      m_Queue( QUEUE),
      m_DeviceDataBool( false)
    {
      SetImplicitArchitecture();
      EnableGpu();
    }

    //! @brief construct from all necessary parameters for the gpu prediction
    //!        this assumes that the inputs are already padded appropriately and rescaled
    //! @param BIAS the bias for each neuron of a hidden layer
    //! @param WEIGHT the weight for each neuron of a hidden layer
    NeuralNetwork::NeuralNetwork
    (
      const util::ShPtr< model:: RescaleFeatureDataSet> &RESCALE_INPUT,
      const util::ShPtr< model:: RescaleFeatureDataSet> &RESCALE_OUTPUT,
      const storage::Vector< Vector< float> > &BIAS,
      const storage::Vector< Matrix< float> > &WEIGHT,
      const CommandQueue &QUEUE
    ) :
      m_TransferFunction( util::Implementation< model::TransferFunctionInterface>( model::TransferSigmoid())),
      m_RescaleInput( RESCALE_INPUT),
      m_RescaleOutput( RESCALE_OUTPUT),
      m_BiasBuffers( BIAS),
      m_WeightBuffers( WEIGHT),
      m_Queue( QUEUE),
      m_GpuMMult( QUEUE),
      m_GpuTranspose( QUEUE),
      m_GpuVMAdd( QUEUE),
      m_GpuSigmoid( QUEUE),
      m_DeviceDataBool( true)
    {
      cl_int error_number = CompilePrograms( util::CPPDataTypes::DataTypeFromTemplate< float>());
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
      for( size_t count( 0), count_end( m_WeightBuffers.GetSize()); count < count_end; ++count)
      {
        m_TransWeightBuffers.PushBack
        (
          Matrix< float>
          (
            m_WeightBuffers( count).GetNumberCols() - m_WeightBuffers( count).GetColPadding(),
            m_WeightBuffers( count).GetNumberRows() - m_WeightBuffers( count).GetRowPadding(),
            m_Queue,
            m_WeightBuffers( count).GetColPadding(),
            m_WeightBuffers( count).GetRowPadding()
          )
        );
      }
    }

    //! copy constructor
    NeuralNetwork *NeuralNetwork::Clone() const
    {
      return new NeuralNetwork( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &NeuralNetwork::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &NeuralNetwork::GetAlias() const
    {
      static const std::string s_Name( "OpenCLNeuralNetwork");
      return s_Name;
    }

  //////////////
  // operator //
  //////////////

  /////////////////
  // data access //
  /////////////////

    //! @brief returns architecture of this neural network
    //! @return first entry is the number of input neurons
    //!         middle entries are the number of hidden neurons in each layer
    //!         last number in vector is the number of output neurons
    const storage::Vector< size_t> &NeuralNetwork::GetArchitecture() const
    {
      return m_Architecture;
    }

    //! @brief set the architecture architecture of this neural network
    //! @param ARCHITECTURE first entry is the number of input neurons
    //!         middle entries are the number of hidden neurons in each layer
    //!         last number in vector is the number of output neurons
    void NeuralNetwork::SetArchitecture( const storage::Vector< size_t> &ARCHITECTURE)
    {
      m_Architecture = ARCHITECTURE;
      m_Bias.Reset();
      m_Weight.Reset();

      for( size_t i( 1); i < ARCHITECTURE.GetSize(); ++i)
      {
        m_Bias.PushBack( linal::Vector< float>( ARCHITECTURE( i)));
        m_Weight.PushBack( linal::Matrix< float>( ARCHITECTURE( i), ARCHITECTURE( i - 1)));
      }
    }

    //! get number inputs
    size_t NeuralNetwork::GetNumberInputs() const
    {
      return m_Weight.FirstElement().GetNumberCols();
    }

    //! get number output neurons
    size_t NeuralNetwork::GetNumberOutputs() const
    {
      return m_Weight.LastElement().GetNumberRows();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief predict result with model using a NOT rescaled feature vector
    //! @param FEATURE not rescaled feature vector
    //! @return predicted result vector using a model
    model::FeatureDataSet< float> NeuralNetwork::PredictWithoutRescaling
    (
      const model::FeatureDataSetInterface< float> &FEATURE
    ) const
    {
      if( !m_DeviceDataBool)
      {
        SetDeviceData();
      }
      const cl_uint row_pad( ( s_blocksize - ( FEATURE.GetNumberFeatures() % s_blocksize)) % s_blocksize);
      const cl_uint col_pad( ( s_blocksize - ( FEATURE.GetFeatureSize() % s_blocksize)) % s_blocksize);
      const Matrix< float> feature( FEATURE.GetMatrix(), m_Queue, row_pad, col_pad);
      Matrix< float> result( operator()( feature));
      return
        model::FeatureDataSet< float>
        (
          result.GetHostMatrix( result.GetRowPadding(), result.GetColPadding()),
          *m_RescaleOutput
        );
    }

    //! @brief Set the scaling of a feature set according to the model
    //! @param FEATURES feature set of interest
    //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
    //!       when operator() is called
    void NeuralNetwork::Rescale( model::FeatureDataSet< float> &FEATURE) const
    {
      if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_RescaleInput)
      {
        FEATURE.DeScale();
        FEATURE.Rescale( *m_RescaleInput);
      }
    }

    //! @brief predict result with model using a rescaled feature vector
    //! @param FEATURE normalized or rescaled feature vector
    //! @return predicted result vector using a model
    model::FeatureDataSet< float> NeuralNetwork::operator()( const model::FeatureDataSetInterface< float> &FEATURE) const
    {
      // handle the case where rescaling is necessary
      if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_RescaleInput)
      {
        model::FeatureDataSet< float> feature( FEATURE);
        feature.Rescale( *m_RescaleInput);
        return PredictWithoutRescaling( feature).DeScale();
      }

      // data is already rescaled
      return PredictWithoutRescaling( FEATURE).DeScale();
    }

    //! @brief predict with gpu without rescaling
    //! @param INPUT the features to predict
    Matrix< float> NeuralNetwork::operator()( const Matrix< float> &INPUT) const
    {
      BCL_Assert
      (
        INPUT.GetNumberCols() == m_WeightBuffers( 0).GetNumberCols(),
        "Input features and weight buffers have non-matching dimensionality!"
      );
      m_GpuTranspose( m_WeightBuffers( 0), m_TransWeightBuffers( 0));

      Matrix< float> hidden( INPUT.GetNumberRows() - INPUT.GetRowPadding(), m_TransWeightBuffers( 0).GetNumberCols() - m_TransWeightBuffers( 0).GetColPadding(), m_Queue, INPUT.GetRowPadding(), m_TransWeightBuffers( 0).GetColPadding());

      m_GpuMMult( INPUT, m_TransWeightBuffers( 0), hidden);

      m_GpuVMAdd( m_BiasBuffers( 0), hidden);
      m_GpuSigmoid.F( hidden, hidden);

      Matrix< float> this_hidden = hidden;

      // remaining layers
      for( size_t k( 1), k_end( m_WeightBuffers.GetSize()); k != k_end; ++k)
      {
        Matrix< float> hidden_input( this_hidden.GetNumberRows() - this_hidden.GetRowPadding(), m_TransWeightBuffers( k).GetNumberCols() - m_TransWeightBuffers( k).GetColPadding(), m_Queue, this_hidden.GetRowPadding(), m_TransWeightBuffers( k).GetColPadding());
        m_GpuTranspose( m_WeightBuffers( k), m_TransWeightBuffers( k));

        m_GpuMMult( this_hidden, m_TransWeightBuffers( k), hidden_input);
        m_GpuVMAdd( m_BiasBuffers( k), hidden_input);
        m_GpuSigmoid.F( hidden_input, hidden_input);
        this_hidden = hidden_input;
      }

      return this_hidden;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read NeuralNetwork from std::istream
    std::istream &NeuralNetwork::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_TransferFunction, ISTREAM);
      io::Serialize::Read( m_RescaleInput, ISTREAM);
      io::Serialize::Read( m_RescaleOutput, ISTREAM);
      io::Serialize::Read( m_Bias, ISTREAM);
      io::Serialize::Read( m_Weight, ISTREAM);

      SetImplicitArchitecture();
      UpdateQueue( GetTools());
      EnableGpu();

      // end
      return ISTREAM;
    }

    //! write NeuralNetwork into std::ostream
    std::ostream &NeuralNetwork::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_TransferFunction, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RescaleInput, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RescaleOutput, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Bias, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Weight, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer NeuralNetwork::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "see http://en.wikipedia.org/wiki/Neural_network"
      );

      parameters.AddInitializer
      (
        "transfer function",
        "function that translates input from neurons in the prior layer into the output of each hidden layer",
        io::Serialization::GetAgent( &m_TransferFunction),
        "Sigmoid"
      );
      parameters.AddInitializer
      (
        "architecture",
        "# of neurons in each layer (input, hidden layer(s), output), e.g. (100, 4, 4, 1)",
        io::Serialization::GetAgentContainerWithCheck( &m_Architecture, io::Serialization::GetAgentWithMin( size_t( 1)))
      );

      return parameters;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief determines the architecture already present in the weight matrix/vectors of the neural network
    void NeuralNetwork::SetImplicitArchitecture()
    {
      // architecture
      m_Architecture.Reset();

      if( m_Weight.GetSize())
      {
        m_Architecture.AllocateMemory( GetNumberLayers() + 1);

        // number of input neurons
        m_Architecture.PushBack( GetNumberInputs());

        // iterate over all layers to get their size
        for
        (
          storage::Vector< linal::Matrix< float> >::const_iterator
            layer_itr( m_Weight.Begin()), layer_itr_end( m_Weight.End());
          layer_itr != layer_itr_end;
          ++layer_itr
        )
        {
          m_Architecture.PushBack( layer_itr->GetNumberRows());
        }
      }
    }

    //! @brief sigmoid transfer function to all elements in a matrix
    //! @param INPUT the input matrix
    void NeuralNetwork::Sigmoid
    (
      Matrix< float> &INPUT,
      Matrix< float> &OUTPUT
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint rows(    INPUT.GetNumberRows());
      const cl_uint cols(    INPUT.GetNumberCols());
      const cl_uint row_pad( INPUT.GetRowPadding());
      const cl_uint col_pad( INPUT.GetColPadding());

      // create kernel
      cl::Kernel kernel( m_Program, "SigmoidKernel", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( s_blocksize, s_blocksize); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( s_blocksize, cols), Tools::RoundUp( s_blocksize, rows));

      // set args for add bias
      error_number  = kernel.setArg( 0, OUTPUT.GetData());
      error_number |= kernel.setArg( 1, INPUT.GetData());
      error_number |= kernel.setArg( 2, rows);
      error_number |= kernel.setArg( 3, cols);
      error_number |= kernel.setArg( 4, col_pad);
      error_number |= kernel.setArg( 5, row_pad);
      BCL_Assert( error_number == CL_SUCCESS, "sigmoid error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    //! @brief adds single row to all rows of a matrix
    //! @param BIAS the bias to add
    //! @param MATRIX the data to add the bias to
    void NeuralNetwork::AddBias
    (
      Vector< float> &BIAS,
      Matrix< float> &MATRIX
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint rows( MATRIX.GetNumberRows());
      const cl_uint cols( MATRIX.GetNumberCols());
      const cl_uint row_pad( MATRIX.GetRowPadding());
      const cl_uint col_pad( MATRIX.GetColPadding());

      // create kernel
      cl::Kernel add_bias_kernel( m_Program, "AddBiasKernel", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( s_blocksize, s_blocksize); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange add_bias_worksize( Tools::RoundUp( s_blocksize, cols), Tools::RoundUp( s_blocksize, rows));

      // set args for add bias
      error_number  = add_bias_kernel.setArg( 0, MATRIX.GetData());
      error_number |= add_bias_kernel.setArg( 1, BIAS.GetData());
      error_number |= add_bias_kernel.setArg( 2, rows);
      error_number |= add_bias_kernel.setArg( 3, cols);
      error_number |= add_bias_kernel.setArg( 4, col_pad);
      error_number |= add_bias_kernel.setArg( 5, row_pad);
      BCL_Assert( error_number == CL_SUCCESS, "add bias error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( add_bias_kernel, offset, add_bias_worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      m_Queue.finish();
    }

    //! @brief sets up everything necessary for gpu computing
    void NeuralNetwork::EnableGpu()
    {
      cl_int error_number = CompilePrograms( util::CPPDataTypes::e_Float);
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));

      m_GpuMMult     = MatrixMultiply< float>( m_Queue);
      m_GpuTranspose = MatrixTranspose< float>( m_Queue);
      m_GpuVMAdd     = VectorMatrixAdd< float>( m_Queue);
      m_GpuSigmoid   = TransferFunctionSigmoid< float>( m_Queue);

      SetDeviceData();
    }

    //! @brief sets device data
    void NeuralNetwork::SetDeviceData() const
    {
      for( size_t count( 0), count_end( GetNumberLayers()); count < count_end; ++count)
      {
        size_t bias_pad( ( s_blocksize - ( m_Bias( count).GetSize() % s_blocksize)) % s_blocksize);
        size_t weight_row_pad( ( s_blocksize - ( m_Weight( count).GetNumberRows() % s_blocksize)) % s_blocksize);
        size_t weight_col_pad( ( s_blocksize - ( m_Weight( count).GetNumberCols() % s_blocksize)) % s_blocksize);

        m_BiasBuffers.PushBack( Vector< float>( m_Bias( count), m_Queue, bias_pad));
        m_WeightBuffers.PushBack( Matrix< float>( m_Weight( count), m_Queue, weight_row_pad, weight_col_pad));
        m_TransWeightBuffers.PushBack( Matrix< float>( m_Weight( count).GetNumberCols(), m_Weight( count).GetNumberRows(), m_Queue, weight_col_pad, weight_row_pad));
      }

      m_DeviceDataBool = true;
    }

    //! @brief compile programs for given precision
    //! @param PRECISION float or double
    //! @return ERROR error that occured, CL_SUCCESS if no error
    cl_int NeuralNetwork::CompilePrograms( const util::CPPDataTypes::Types &PRECISION)
    {
      cl_int error_number( CL_SUCCESS);

      // compile the program
      cl::Program::Sources source;
      const KernelSourceFile ann_kernels_source_file( "ann_kernels.cl");

      // construct opencl device
      const Device device( m_Queue.GetDevice( &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "cannot get device for command queue");

      // construct kernel source strings
      const std::string ann_kernels_source( ann_kernels_source_file.GetSource( PRECISION, device.Extensions()));

      // check if kernel strings are empty
      if
      (
          ann_kernels_source.empty()
      )
      {
        return CL_INVALID_KERNEL_DEFINITION;
      }

      // pushback strings to program sources vector
      source.push_back( std::make_pair( ann_kernels_source.c_str(), ann_kernels_source.length()));

      // create the program
      cl::Program &current_program( m_Program);
      current_program = cl::Program( m_Queue.GetContext(), source, &error_number);
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

    //! @brief responsible for updating to a valid queue
    //! @param TOOLS opencl tools
    void NeuralNetwork::UpdateQueue( Tools &TOOLS)
    {
      if( !TOOLS.HasCommandQueues())
      {
        return;
      }

      m_Queue = GetTools().GetFirstCommandQueue();

      cl_int error_number = CompilePrograms( util::CPPDataTypes::e_Float);
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));

      m_GpuMMult     = MatrixMultiply< float>( m_Queue);
      m_GpuTranspose = MatrixTranspose< float>( m_Queue);
      m_GpuVMAdd     = VectorMatrixAdd< float>( m_Queue);
      m_GpuSigmoid   = TransferFunctionSigmoid< float>( m_Queue);
    }

  } // namespace opencl
} // namespace bcl
