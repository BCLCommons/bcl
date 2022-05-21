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
#include "opencl/bcl_opencl_approximator_simple_propagation.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "math/bcl_math_statistics.h"
#include "model/bcl_model_transfer_sigmoid.h"
#include "opencl/bcl_opencl_neural_network.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    const char *ApproximatorSimplePropagation::s_CLCompilerOptions =
            "-cl-mad-enable -cl-fast-relaxed-math";

    const cl_uint ApproximatorSimplePropagation::s_Blocksize = 16;

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ApproximatorSimplePropagation::s_Instance
    (
      util::Enumerated< model::ApproximatorBase>::AddInstance( new ApproximatorSimplePropagation())
    );

    //! @brief default constructor
    ApproximatorSimplePropagation::ApproximatorSimplePropagation() :
      m_HiddenArchitecture(),
      m_Alpha(),
      m_Eta(),
      m_StepsPerCall()
    {
    }

    //! @brief constructor from training data, transfer, rescale, and objective functions
    //! @param TRAINING_DATA data to train the NeuralNetwork on
    //! @param ARCHITECTURE the architecture of the neural network, determines size of the matrices
    //! @param ALPHA learning momentum
    //! @param ETA learning rate
    //! @param STEPS_PER_CALL steps per call
    //! @param QUEUE command queue
    ApproximatorSimplePropagation::ApproximatorSimplePropagation
    (
      util::ShPtr< descriptor::Dataset> &TRAINING_DATA,
      const storage::Vector< size_t> &ARCHITECTURE,
      const float ALPHA,
      const float ETA,
      const size_t STEPS_PER_CALL,
      const CommandQueue &QUEUE
    ) :
      m_HiddenArchitecture( ARCHITECTURE),
      m_Alpha( ALPHA),
      m_Eta( ETA),
      m_StepsPerCall( STEPS_PER_CALL),
      m_Queue( QUEUE),
      m_GpuRMSD( m_Queue),
      m_GpuMMult( m_Queue),
      m_GpuTranspose( m_Queue),
      m_GpuMMAdd( m_Queue),
      m_GpuVMAdd( m_Queue),
      m_GpuSigmoid( m_Queue)
    {
      cl_int error_number( CL_SUCCESS);
      m_Program = KernelSources::Compile( GetKernelSources().e_ArtificialNeuralNetwork, util::CPPDataTypes::e_Float, m_Queue, std::string( s_CLCompilerOptions), &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
      if( !TRAINING_DATA.IsDefined())
      {
        return;
      }
      // set and rescale training data set
      SetTrainingData( TRAINING_DATA);
    }

    //! @brief copy constructor
    //! @return a new ApproximatorSimplePropagation copied from this instance
    ApproximatorSimplePropagation *ApproximatorSimplePropagation::Clone() const
    {
      return new ApproximatorSimplePropagation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ApproximatorSimplePropagation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ApproximatorSimplePropagation::GetAlias() const
    {
      static const std::string s_Name( "OpenCLSimplePropagation");
      return s_Name;
    }

    //! @brief set training data set for a specific iterate in approximater framework
    //! @param DATA training data set
    void ApproximatorSimplePropagation::SetTrainingData
    (
      util::ShPtr< descriptor::Dataset> &DATA
    )
    {
      // rescale function for in an output and denormalization
      m_TrainingData = DATA;
      DATA->GetFeatures().Rescale( NeuralNetwork::s_DefaultInputRange);
      DATA->GetResults().Rescale( model::TransferSigmoid().GetDynamicOutputRange());

      BCL_MessageDbg( "scaled training data: " + util::Format()( *m_TrainingData->GetFeaturesPtr()));

      // now that the training data has been defined, set up the internal data members for the desired architecture
      SetupArchitecture();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief gets model from gpu
    //! @return shptr to the new model interface
    util::ShPtr< model::Interface> ApproximatorSimplePropagation::GetCurrentModel() const
    {
      m_Bias.Reset();
      m_Weight.Reset();
      for( size_t ctr( 0), ctr_end( m_BiasBuffers.GetSize()); ctr < ctr_end; ++ctr)
      {
        m_Bias.PushBack( m_BiasBuffers( ctr).GetHostVector( m_BiasBuffers( ctr).GetPadding()));
        m_Weight.PushBack( m_WeightBuffers( ctr).GetHostMatrix( m_WeightBuffers( ctr).GetRowPadding(), m_WeightBuffers( ctr).GetColPadding()));
      }
      // make a new model out of the current data members
      return util::ShPtr< model::Interface>
      (
        new NeuralNetwork
        (
          GetRescaleFeatureDataSet(),
          GetRescaleResultDataSet(),
          m_Bias,
          m_Weight,
          util::Implementation< model::TransferFunctionInterface>( model::TransferSigmoid()),
          m_Queue
        )
      );
    }

    //! @brief construct a model from the current iterate
    //! @return shptr to the new model interface
    util::ShPtr< model::Interface> ApproximatorSimplePropagation::GetCurrentGPUModel() const
    {
      // make a new model out of the current data members
      return util::ShPtr< model::Interface>
      (
        new NeuralNetwork
        (
          GetRescaleFeatureDataSet(),
          GetRescaleResultDataSet(),
          m_BiasBuffers,
          m_WeightBuffers,
          m_Queue
        )
      );
    }

    //! @brief returns the current approximation
    //! @return current argument result pair
    const util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> >
      ApproximatorSimplePropagation::GetCurrentApproximation() const
    {
      util::ShPtr< model::Interface> model( GetCurrentModel());
      return
        util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> >
        (
          new storage::Pair< util::ShPtr< model::Interface>, float>
          (
            model,
            m_ObjectiveFunction->operator()( model)
          )
        );
    }

    //! @brief conducts the next approximation step and stores the approximation
    void ApproximatorSimplePropagation::Next()
    {
      if( m_MonitorFeaturesOnDevice.GetNumberOfElements() < 1)
      {
        const size_t mon_feat_pad_row( ( s_Blocksize - ( m_ObjectiveFunction->GetData()->GetSize() % s_Blocksize)) % s_Blocksize);
        const size_t mon_feat_pad_col( ( s_Blocksize - ( m_ObjectiveFunction->GetData()->GetFeatureSize() % s_Blocksize)) % s_Blocksize);
        m_MonitorFeaturesOnDevice = Matrix< float>( m_ObjectiveFunction->GetData()->GetFeaturesPtr()->GetMatrix(), m_Queue, mon_feat_pad_row, mon_feat_pad_col);

        const size_t mon_result_pad_row( ( s_Blocksize - ( m_ObjectiveFunction->GetData()->GetSize() % s_Blocksize)) % s_Blocksize);
        const size_t mon_result_pad_col( ( s_Blocksize - ( m_ObjectiveFunction->GetData()->GetResultSize() % s_Blocksize)) % s_Blocksize);
        m_MonitorTargetOnDevice   = Matrix< float>( GetRescaleResultDataSet()->operator()( *m_ObjectiveFunction->GetData()->GetResultsPtr()).GetMatrix(), m_Queue, mon_result_pad_row, mon_result_pad_col);
      }

      float rmsd_train( std::numeric_limits< float>::max());
      float rmsd_mon  ( std::numeric_limits< float>::max());
      for( size_t count( 0), count_end( m_StepsPerCall); count < count_end; ++count)
      {
        Train();
        UpdateWeights();
      }

      util::ShPtr< ModelInterface> current_model( GetCurrentGPUModel());
      Matrix< float> train_result = current_model->operator()( m_TrainingFeaturesOnDevice);
      rmsd_train = m_GpuRMSD( train_result, m_TrainingTargetOnDevice);

      Matrix< float> mon_result( current_model->operator()( m_MonitorFeaturesOnDevice));
      rmsd_mon = m_GpuRMSD( mon_result, m_MonitorTargetOnDevice);

      BCL_MessageStd( "Training rmsd: " + util::Format()( rmsd_train) + "\tmonitor rmsd: " + util::Format()( rmsd_mon))

      util::ShPtr< model::Interface> current( GetCurrentModel());

      this->GetTracker().Track
      (
        util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> >
        (
          new storage::Pair< util::ShPtr< model::Interface>, float>( current, m_ObjectiveFunction->operator()( current))
        )
      );
    }

    //! read ApproximatorSimplePropagation from std::istream
    std::istream &ApproximatorSimplePropagation::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_TrainingData, ISTREAM);
      io::Serialize::Read( m_HiddenArchitecture, ISTREAM);
      io::Serialize::Read( m_Bias, ISTREAM);
      io::Serialize::Read( m_Weight, ISTREAM);

      // return
      return ISTREAM;
    }

    //! write ApproximatorSimplePropagation into std::ostream
    std::ostream &ApproximatorSimplePropagation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_TrainingData, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_HiddenArchitecture, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Bias, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Weight, OSTREAM, INDENT) << '\n';

      // end
      return OSTREAM;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ApproximatorSimplePropagation::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "trains a neural network using backpropagation (see http://en.wikipedia.org/wiki/Backpropagation)"
      );

      parameters.AddInitializer
      (
        "hidden architecture",
        "# of neurons in each hidden layer, e.g. (100) declares that there is 1 hidden layer w/ 100 neurons",
        io::Serialization::GetAgentContainerWithCheck
        (
          &m_HiddenArchitecture,
          io::Serialization::GetAgentWithMin( size_t( 1))
        )
      );
      parameters.AddInitializer
      (
        "eta",
        "learning rate",
        io::Serialization::GetAgentWithRange( &m_Eta, 0.0, 1.0),
        "10.0"
      );
      parameters.AddInitializer
      (
        "alpha",
        "momentum; how long a change in weights persists",
        io::Serialization::GetAgentWithRange( &m_Alpha, 0.0, 1.0),
        "0.5"
      );
      parameters.AddInitializer
      (
        "steps per call",
        "after how many steps to check rmsd",
        io::Serialization::GetAgent( &m_StepsPerCall),
        "500"
      );
      parameters.AddInitializer
      (
        "objective function",
        "function that evaluates the model after each batch step",
        io::Serialization::GetAgent( &m_ObjectiveFunction->GetImplementation()),
        "RMSD"
      );
      parameters.AddInitializer
      (
        "initial network",
        "provides an initial network model from which to start training from",
        io::Serialization::GetAgentInputFilename( &m_NetworkFilename),
        ""
      );
      return parameters;
    }

    //! run forward through ANN and compute hidden terms m_Hidden (test)
    void ApproximatorSimplePropagation::CalcHiddenTerms()
    {
      // cpu equivalent
      // // input layer
      // m_HiddenInput( 0) = m_Weight( 0) * FEATURE + m_Bias( 0);
      // m_Hidden( 0) = m_TransferFunction->F( m_HiddenInput( 0));
      //
      // // remaining layers
      // for( size_t k( 1); k != m_Hidden.GetSize(); ++k)
      // {
      //   m_HiddenInput( k) = m_Weight( k) * m_Hidden( k - 1) + m_Bias( k);
      //   m_Hidden( k) = m_TransferFunction->F( m_HiddenInput( k));
      // }

      m_GpuTranspose( m_WeightBuffers( 0), m_TransWeightBuffers( 0));

      m_GpuMMult( m_TrainingFeaturesOnDevice, m_TransWeightBuffers( 0), m_HiddenInputBuffers( 0));

      m_GpuVMAdd( m_BiasBuffers( 0), m_HiddenInputBuffers( 0));

      m_GpuSigmoid.F( m_HiddenInputBuffers( 0), m_HiddenBuffers( 0));

      // remaining layers
      for( size_t k( 1), k_end( m_HiddenBuffers.GetSize()); k != k_end; ++k)
      {
        m_GpuTranspose( m_WeightBuffers( k), m_TransWeightBuffers( k));

        m_GpuMMult( m_HiddenBuffers( k - 1), m_TransWeightBuffers( k), m_HiddenInputBuffers( k));

        m_GpuVMAdd( m_BiasBuffers( k), m_HiddenInputBuffers( k));

        m_GpuSigmoid.F( m_HiddenInputBuffers( k), m_HiddenBuffers( k));
      }
    }

    //! run backward through ANN and compute err terms m_Errors (train)
    void ApproximatorSimplePropagation::CalcErrorTerms()
    {
      // cpu eqivalent

      // // output layer
      // m_Errors.LastElement() = RESULT;
      // m_Errors.LastElement() -= m_Hidden.LastElement();
      // const float *hid( m_Hidden.LastElement().Begin());
      // const float *hid_input( m_HiddenInput.LastElement().Begin());
      // for
      // (
      //   float *err( m_Errors.LastElement().Begin()), *err_end( m_Errors.LastElement().End());
      //   err != err_end;
      //   ++err, ++hid, ++hid_input
      // )
      // {
      //   ( *err) *= m_TransferFunction->dF( *hid_input, *hid);
      // }
      //
      // // all other layers
      // for( size_t i( m_Hidden.GetSize() - 1); i > 0; --i)
      // {
      //   const float *hid( m_Hidden( i - 1).Begin());
      //   const float *hid_input( m_HiddenInput( i - 1).Begin());
      //
      //   m_Errors( i - 1) = m_Errors( i) * m_Weight( i);
      //
      //   for
      //   (
      //     float *err( m_Errors( i - 1).Begin()), *err_end( m_Errors( i - 1).End());
      //     err != err_end;
      //     ++err, ++hid, ++hid_input
      //   )
      //   {
      //     ( *err) *= m_TransferFunction->dF( *hid_input, *hid);
      //   }
      // }

      OutputLayerErrorTerms();

      for( size_t i( m_HiddenBuffers.GetSize() - 1); i > 0; --i)
      {
        m_GpuMMult( m_ErrorsBuffers( i), m_WeightBuffers( i), m_ErrorsBuffers( i - 1));
        OtherLayerErrorTerms( i);
      }

    }

    //! compute changes m_SlopeBias/Weight to be applied on ANN (train)
    void ApproximatorSimplePropagation::CalcChangeTerms()
    {
      // cpu equivalent

      // // first layer
      // m_SlopesBias( 0) += m_Errors( 0);
      //
      // math::OuterProductOperation( m_SlopesWeight( 0), math::PlusEquals< float>(), m_Errors( 0), FEATURE);
      //
      // // all other layers
      // for( size_t i( 1); i < m_Hidden.GetSize(); ++i)
      // {
      //   m_SlopesBias( i) += m_Errors( i);
      //   math::OuterProductOperation( m_SlopesWeight( i), math::PlusEquals< float>(), m_Errors( i), m_Hidden( i - 1));
      // }

      ReduceErrorsToVector( m_ErrorsBuffers( 0), m_SlopesBiasBuffers( 0));
      m_GpuTranspose( m_ErrorsBuffers( 0), m_TransErrorsBuffers( 0));
      m_GpuMMult( m_TransErrorsBuffers( 0), m_TrainingFeaturesOnDevice, m_SlopesWeightBuffers( 0));

      for( size_t i( 1); i < m_HiddenBuffers.GetSize(); ++i)
      {
        ReduceErrorsToVector( m_ErrorsBuffers( i), m_SlopesBiasBuffers( i));
        m_GpuTranspose( m_ErrorsBuffers( i), m_TransErrorsBuffers( i));

        m_GpuMMult( m_TransErrorsBuffers( i), m_HiddenBuffers( i - 1), m_SlopesWeightBuffers( i));
      }
    }

    //! train ANN with a feature
    void ApproximatorSimplePropagation::Train()
    {
      // training step
      CalcHiddenTerms();
      CalcErrorTerms();
      CalcChangeTerms();

    }

    void ApproximatorSimplePropagation::UpdateWeights()
    {
      for( size_t i( 0); i < m_HiddenBuffers.GetSize(); ++i)
      {
        ApplyBiasChanges( i);
        ApplyWeightChanges( i);
        AddBiasChangeTerms( i);
        m_GpuMMAdd( m_WeightBuffers( i), m_ChangeWeightBuffers( i));
      }
    }

    //! @brief calculates the output layer error terms
    void ApproximatorSimplePropagation::OutputLayerErrorTerms()
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;
      BCL_MessageDbg( "OutputLayerErrorTerms()");
      // dimensions
      const cl_uint rows(    m_ErrorsBuffers.LastElement().GetNumberRows());
      const cl_uint cols(    m_ErrorsBuffers.LastElement().GetNumberCols());
      const cl_uint row_pad( m_ErrorsBuffers.LastElement().GetRowPadding());
      const cl_uint col_pad( m_ErrorsBuffers.LastElement().GetColPadding());

      // create kernel
      cl::Kernel kernel( m_Program, "ErrorKKernel", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( s_Blocksize, s_Blocksize); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( s_Blocksize, cols), Tools::RoundUp( s_Blocksize, rows));

      BCL_MessageDbg( "m_TrainingTargetOnDevice: " + util::Format()( m_TrainingTargetOnDevice.GetHostMatrix()));
      BCL_MessageDbg( "m_ErrorsBuffers: " + util::Format()( m_ErrorsBuffers.LastElement().GetHostMatrix()));
      BCL_MessageDbg( "m_HiddenBuffers: " + util::Format()( m_HiddenBuffers.LastElement().GetHostMatrix()));

      error_number  = kernel.setArg( 0, m_ErrorsBuffers.LastElement().GetData());
      error_number |= kernel.setArg( 1, m_TrainingTargetOnDevice.GetData());
      error_number |= kernel.setArg( 2, m_HiddenBuffers.LastElement().GetData());
      error_number |= kernel.setArg( 3, rows);
      error_number |= kernel.setArg( 4, cols);
      error_number |= kernel.setArg( 5, col_pad);
      error_number |= kernel.setArg( 6, row_pad);
      BCL_Assert( error_number == CL_SUCCESS, "error_k error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      BCL_MessageDbg( "output layer error terms: " + util::Format()( m_ErrorsBuffers.LastElement().GetHostMatrix()));
    }

    //! @brief calculates the hidden layer error terms
    //! @param LAYER the layer number
    void ApproximatorSimplePropagation::OtherLayerErrorTerms( const size_t &LAYER)
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint rows(    m_ErrorsBuffers( LAYER - 1).GetNumberRows());
      const cl_uint cols(    m_ErrorsBuffers( LAYER - 1).GetNumberCols());
      const cl_uint row_pad( m_ErrorsBuffers( LAYER - 1).GetRowPadding());
      const cl_uint col_pad( m_ErrorsBuffers( LAYER - 1).GetColPadding());

      // create kernel
      cl::Kernel kernel( m_Program, "ErrorJKernel", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( s_Blocksize, s_Blocksize); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( s_Blocksize, cols), Tools::RoundUp( s_Blocksize, rows));

      // output, i, y
      error_number  = kernel.setArg( 0, m_ErrorsBuffers( LAYER - 1).GetData());
      error_number |= kernel.setArg( 1, m_ErrorsBuffers( LAYER - 1).GetData());
      error_number |= kernel.setArg( 2, m_HiddenBuffers( LAYER - 1).GetData());
      error_number |= kernel.setArg( 3, rows);
      error_number |= kernel.setArg( 4, cols);
      error_number |= kernel.setArg( 5, col_pad);
      error_number |= kernel.setArg( 6, row_pad);
      BCL_Assert( error_number == CL_SUCCESS, "error_k error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

    }

    //! @brief adds bias change terms to previous bias
    //! @param LAYER the layer number
    void ApproximatorSimplePropagation::AddBiasChangeTerms( const size_t &LAYER)
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint rows(    1);
      const cl_uint cols(    m_BiasBuffers( LAYER).GetSize());
      const cl_uint row_pad( 0);
      const cl_uint col_pad( m_BiasBuffers( LAYER).GetPadding());

      // create kernel
      cl::Kernel kernel( m_Program, "AddKernel", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( s_Blocksize, s_Blocksize); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( s_Blocksize, cols), Tools::RoundUp( s_Blocksize, rows));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, m_BiasBuffers( LAYER).GetData());
      error_number |= kernel.setArg( 1, m_ChangeBiasBuffers( LAYER).GetData());
      error_number |= kernel.setArg( 2, rows);
      error_number |= kernel.setArg( 3, cols);
      error_number |= kernel.setArg( 4, col_pad);
      error_number |= kernel.setArg( 5, row_pad);
      BCL_Assert( error_number == CL_SUCCESS, "update bias error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

    }

    //! @brief helper function to do a column-wise reduction on the bias errors
    //! @param ERRORS the errors buffer
    //! @param REDUCED_ERRORS the reduced error vector
    void ApproximatorSimplePropagation::ReduceErrorsToVector
    (
      Matrix< float> &ERRORS, Vector< float> &REDUCED_ERRORS
    )
    {
      // error catching
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint rows( ERRORS.GetNumberRows());
      const cl_uint cols( ERRORS.GetNumberCols());

      // blocksize
      const cl_uint block_size( 256);

      // construct kernel
      cl::Kernel kernel( m_Program, "SumBiasColumnsKernel", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // num_groups equals size of partially reduced array returned from gpu partial sum
      const size_t num_groups( ( rows - 1) / block_size + 1);

      // construct partial sum result matrix
      linal::Matrix< float> temp_result( cols, num_groups);

      // thread launch dimensions
      const cl::NDRange block_dim( block_size, 1);
      const cl::NDRange offset;
      const cl::NDRange kernel_worksize( Tools::RoundUp( block_size, rows), cols);

      // allocate result array on device
      Matrix< float> temp_array( temp_result, m_Queue);

      // set args for kernel
      error_number  = kernel.setArg( 0, ERRORS.GetData());
      error_number |= kernel.setArg( 1, temp_array.GetData());
      error_number |= kernel.setArg( 2, rows);
      error_number |= kernel.setArg( 3, cols);
      error_number |= kernel.setArg( 4, block_size * sizeof( float), 0);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // launch kernel
      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, kernel_worksize, block_dim, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // read back partially reduced result array
      temp_result = temp_array.GetHostMatrix();

      // final reduction result vector
      linal::Vector< float> result_vector( size_t( cols), 0.0);

      // finish reduction of vectors/cols
      for( cl_uint i( 0); i < cols; ++i)
      {
        result_vector( i) = std::accumulate( temp_result.Begin() + i * num_groups, temp_result.Begin() + ( i + 1) * num_groups, float( 0));
      }

      // write back fully reduced vector to device
      REDUCED_ERRORS = Vector< float>( result_vector.CreateSubVector( REDUCED_ERRORS.GetSize() - REDUCED_ERRORS.GetPadding()), m_Queue, REDUCED_ERRORS.GetPadding());
    }

    //! @brief applies change terms
    void ApproximatorSimplePropagation::ApplyWeightChanges( const size_t &LAYER)
    {
      // error catching
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint rows(    m_SlopesWeightBuffers( LAYER).GetNumberRows());
      const cl_uint cols(    m_SlopesWeightBuffers( LAYER).GetNumberCols());
      const cl_uint row_pad( m_SlopesWeightBuffers( LAYER).GetRowPadding());
      const cl_uint col_pad( m_SlopesWeightBuffers( LAYER).GetColPadding());

      cl::Kernel kernel( m_Program, "DeltaKernel", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( s_Blocksize, s_Blocksize); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( s_Blocksize, cols), Tools::RoundUp( s_Blocksize, rows));

      // set args for kernel
      error_number =  kernel.setArg( 0, m_ChangeWeightBuffers( LAYER).GetData());
      error_number |= kernel.setArg( 1, m_SlopesWeightBuffers( LAYER).GetData());
      error_number |= kernel.setArg( 2, m_Eta);
      error_number |= kernel.setArg( 3, m_Alpha);
      error_number |= kernel.setArg( 4, rows);
      error_number |= kernel.setArg( 5, cols);
      error_number |= kernel.setArg( 6, col_pad);
      error_number |= kernel.setArg( 7, row_pad);
      BCL_Assert( error_number == CL_SUCCESS, "delta kernel error: " + opencl::Tools::ErrorString( error_number));

      // launch kernel
      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

    }

    //! @brief applies change terms
    void ApproximatorSimplePropagation::ApplyBiasChanges( const size_t &LAYER)
    {
      // error catching
      cl_int error_number = CL_SUCCESS;
      BCL_MessageDbg( "m_SlopesWeightBuffers( LAYER): a" + util::Format()( m_SlopesWeightBuffers( LAYER).GetHostMatrix()));

      // dimensions
      const cl_uint rows(    1);
      const cl_uint cols(    m_SlopesBiasBuffers( LAYER).GetSize());
      const cl_uint row_pad( 0);
      const cl_uint col_pad( m_SlopesBiasBuffers( LAYER).GetPadding());

      cl::Kernel kernel( m_Program, "DeltaKernel", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( s_Blocksize, s_Blocksize); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( s_Blocksize, cols), Tools::RoundUp( s_Blocksize, rows));

      // set args for kernel
      error_number =  kernel.setArg( 0, m_ChangeBiasBuffers( LAYER).GetData());
      error_number |= kernel.setArg( 1, m_SlopesBiasBuffers( LAYER).GetData());
      error_number |= kernel.setArg( 2, m_Eta);
      error_number |= kernel.setArg( 3, m_Alpha);
      error_number |= kernel.setArg( 4, rows);
      error_number |= kernel.setArg( 5, cols);
      error_number |= kernel.setArg( 6, col_pad);
      error_number |= kernel.setArg( 7, row_pad);
      BCL_Assert( error_number == CL_SUCCESS, "delta kernel error: " + opencl::Tools::ErrorString( error_number));

      // launch kernel
      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
    }

    //! @brief sets up the neural network with a particular architecture, after training data was set
    void ApproximatorSimplePropagation::SetupArchitecture()
    {
      BCL_Assert
      (
        m_TrainingData.IsDefined(),
        "SetupArchitecture requires valid training data!"
      );

      // create the complete architecture, including input and output neurons
      // start by making a vector with just the input neurons
      util::ShPtr< NeuralNetwork> network;
      storage::Vector< size_t> architecture;
      bool read_net( false);
      if( m_NetworkFilename.size() > 3)
      {
        read_net = true;
        io::IFStream in;
        io::File::TryOpenIFStream( in, m_NetworkFilename);
        double rmsd( util::GetUndefined< double>());
        std::string paren;
        io::Serialize::Read( rmsd, in);
        BCL_MessageStd( "initial network had a score of " + util::Format()( rmsd));
        io::Serialize::Read( paren, in);
        io::Serialize::Read( network, in);
        io::File::CloseClearFStream( in);
        architecture = network->GetArchitecture();
      }
      else
      {
        architecture.PushBack( m_TrainingData->GetFeatureSize());
        // append the hidden neurons
        if( !m_HiddenArchitecture.IsEmpty())
        {
          architecture.Append( m_HiddenArchitecture);
        }

        // add the output neurons
        architecture.PushBack( m_TrainingData->GetResultSize());
      }

      m_Bias.Reset();
      m_Weight.Reset();
      m_Bias.AllocateMemory( architecture.GetSize() - 1);
      m_Weight.AllocateMemory( architecture.GetSize() - 1);

      if( m_TrainingData->GetFeaturesPtr()->GetNumberFeatures() <= 0)
      {
        return;
      }

      // ensure that there are no empty layers
      BCL_Assert
      (
        math::Statistics::MinimumValue( architecture.Begin(), architecture.End()) > 0,
        "Each layer must have at least one neuron!"
      );

      // initialize data
      for( size_t i( 1); i < architecture.GetSize(); ++i)
      {
        m_Bias.PushBack( linal::Vector< float>( architecture( i), 0.0));
        m_Weight.PushBack( linal::Matrix< float>( architecture( i), architecture( i - 1), 0.0));
      }

      if( read_net)
      {
        m_Weight = network->GetWeight();
        m_Bias   = network->GetBias();
        BCL_MessageStd( "using weights and bias from the provided initial network file")
      }

      // randomize the initial values of bias and weight
      for
      (
        storage::Vector< linal::Vector< float> >::iterator itr( m_Bias.Begin()), itr_end( m_Bias.End());
        itr != itr_end;
        ++itr
      )
      {
        if( !read_net)
        {
          math::Statistics::SetRand( itr->Begin(), itr->End(), -0.1, 0.1);
        }
        const size_t pad( ( s_Blocksize - ( itr->GetSize() % s_Blocksize)) % s_Blocksize);
        m_BiasBuffers.PushBack( Vector< float>( *itr, m_Queue, pad));
      }

      for
      (
        storage::Vector< linal::Matrix< float> >::iterator itr( m_Weight.Begin()), itr_end( m_Weight.End());
        itr != itr_end;
        ++itr
      )
      {
        if( !read_net)
        {
          math::Statistics::SetRand( itr->Begin(), itr->End(), -0.1, 0.1);
        }
        const size_t pad_row( ( s_Blocksize - ( itr->GetNumberRows() % s_Blocksize)) % s_Blocksize);
        const size_t pad_col( ( s_Blocksize - ( itr->GetNumberCols() % s_Blocksize)) % s_Blocksize);
        m_WeightBuffers.PushBack( Matrix< float>( *itr, m_Queue, pad_row, pad_col));
        m_TransWeightBuffers.PushBack( Matrix< float>( itr->GetNumberCols(), itr->GetNumberRows(), m_Queue, pad_col, pad_row));
      }

      BCL_MessageDbg( "starting weights: " + util::Format()( m_Weight) + "\nstarting biases: " + util::Format()( m_Bias));
      const size_t train_feat_pad_row( ( s_Blocksize - ( m_TrainingData->GetSize() % s_Blocksize)) % s_Blocksize);
      const size_t train_feat_pad_col( ( s_Blocksize - ( m_TrainingData->GetFeatureSize() % s_Blocksize)) % s_Blocksize);
      m_TrainingFeaturesOnDevice = Matrix< float>( m_TrainingData->GetFeaturesPtr()->GetMatrix(), m_Queue, train_feat_pad_row, train_feat_pad_col);

      const size_t train_result_pad_row( ( s_Blocksize - ( m_TrainingData->GetSize() % s_Blocksize)) % s_Blocksize);
      const size_t train_result_pad_col( ( s_Blocksize - ( m_TrainingData->GetResultSize() % s_Blocksize)) % s_Blocksize);
      m_TrainingTargetOnDevice   = Matrix< float>( m_TrainingData->GetResultsPtr()->GetMatrix(), m_Queue, train_result_pad_row, train_result_pad_col);

      // initialize data
      for( size_t i( 1); i < architecture.GetSize(); ++i)
      {
        m_HiddenBuffers.PushBack(       Matrix< float>( m_TrainingFeaturesOnDevice.GetNumberRows() - m_TrainingFeaturesOnDevice.GetRowPadding(), m_BiasBuffers( i - 1).GetSize() - m_BiasBuffers( i - 1).GetPadding(), m_Queue, m_TrainingFeaturesOnDevice.GetRowPadding(), m_BiasBuffers( i - 1).GetPadding(), 0));
        m_HiddenInputBuffers.PushBack(  Matrix< float>( m_TrainingFeaturesOnDevice.GetNumberRows() - m_TrainingFeaturesOnDevice.GetRowPadding(), m_BiasBuffers( i - 1).GetSize() - m_BiasBuffers( i - 1).GetPadding(), m_Queue, m_TrainingFeaturesOnDevice.GetRowPadding(), m_BiasBuffers( i - 1).GetPadding(), 0));
        m_ErrorsBuffers.PushBack(       Matrix< float>( m_TrainingFeaturesOnDevice.GetNumberRows() - m_TrainingFeaturesOnDevice.GetRowPadding(), m_BiasBuffers( i - 1).GetSize() - m_BiasBuffers( i - 1).GetPadding(), m_Queue, m_TrainingFeaturesOnDevice.GetRowPadding(), m_BiasBuffers( i - 1).GetPadding(), 0));
        m_TransErrorsBuffers.PushBack(  Matrix< float>( m_BiasBuffers( i - 1).GetSize() - m_BiasBuffers( i - 1).GetPadding(), m_TrainingFeaturesOnDevice.GetNumberRows() - m_TrainingFeaturesOnDevice.GetRowPadding(), m_Queue, m_BiasBuffers( i - 1).GetPadding(), m_TrainingFeaturesOnDevice.GetRowPadding(), 0));
        m_SlopesBiasBuffers.PushBack(   Vector< float>( m_BiasBuffers( i - 1).GetSize() - m_BiasBuffers( i - 1).GetPadding(), m_Queue, m_BiasBuffers( i - 1).GetPadding(), 0));
        m_ChangeBiasBuffers.PushBack(   Vector< float>( m_BiasBuffers( i - 1).GetSize() - m_BiasBuffers( i - 1).GetPadding(), m_Queue, m_BiasBuffers( i - 1).GetPadding(), 0));
        m_SlopesWeightBuffers.PushBack( Matrix< float>( m_WeightBuffers( i - 1).GetNumberRows() - m_WeightBuffers( i - 1).GetRowPadding(), m_WeightBuffers( i - 1).GetNumberCols() - m_WeightBuffers( i - 1).GetColPadding(), m_Queue, m_WeightBuffers( i - 1).GetRowPadding(), m_WeightBuffers( i - 1).GetColPadding(), 0));
        m_ChangeWeightBuffers.PushBack( Matrix< float>( m_WeightBuffers( i - 1).GetNumberRows() - m_WeightBuffers( i - 1).GetRowPadding(), m_WeightBuffers( i - 1).GetNumberCols() - m_WeightBuffers( i - 1).GetColPadding(), m_Queue, m_WeightBuffers( i - 1).GetRowPadding(), m_WeightBuffers( i - 1).GetColPadding(), 0));
      }
    }

    //! @brief responsible for updating to a valid queue
    //! @param TOOLS opencl tools
    void ApproximatorSimplePropagation::UpdateQueue( Tools &TOOLS)
    {
      if( !TOOLS.HasCommandQueues())
      {
        return;
      }

      m_Queue = GetTools().GetFirstCommandQueue();

      cl_int error_number( CL_SUCCESS);
      m_Program = KernelSources::Compile( GetKernelSources().e_ArtificialNeuralNetwork, util::CPPDataTypes::e_Float, m_Queue, std::string( s_CLCompilerOptions), &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));

      m_GpuRMSD      = RMSD( m_Queue);
      m_GpuMMult     = MatrixMultiply< float>( m_Queue);
      m_GpuTranspose = MatrixTranspose< float>( m_Queue);
      m_GpuMMAdd     = MatrixAdd< float>( m_Queue);
      m_GpuVMAdd     = VectorMatrixAdd< float>( m_Queue);
      m_GpuSigmoid   = TransferFunctionSigmoid< float>( m_Queue);
    }

  } // namespace opencl
} // namespace bcl
