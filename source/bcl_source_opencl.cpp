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
#include "opencl/bcl_opencl_approximator_resilient_propagation.h"

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
    const char *ApproximatorResilientPropagation::s_CLCompilerOptions =
            "-cl-mad-enable -cl-fast-relaxed-math";

    const cl_uint ApproximatorResilientPropagation::s_Blocksize = 16;

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ApproximatorResilientPropagation::s_Instance
    (
      util::Enumerated< model::ApproximatorBase>::AddInstance( new ApproximatorResilientPropagation())
    );

    //! @brief default constructor
    ApproximatorResilientPropagation::ApproximatorResilientPropagation() :
      m_HiddenArchitecture(),
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
    ApproximatorResilientPropagation::ApproximatorResilientPropagation
    (
      util::ShPtr< descriptor::Dataset> &TRAINING_DATA,
      const storage::Vector< size_t> &ARCHITECTURE,
      const size_t STEPS_PER_CALL,
      const CommandQueue &QUEUE
    ) :
      m_HiddenArchitecture( ARCHITECTURE),
      m_StepsPerCall( STEPS_PER_CALL),
      m_Queue( QUEUE),
      m_GpuRMSD     ( m_Queue),
      m_GpuMMult    ( m_Queue),
      m_GpuTranspose( m_Queue),
      m_GpuMMAdd    ( m_Queue),
      m_GpuVMAdd    ( m_Queue),
      m_GpuSigmoid  ( m_Queue)
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
    //! @return a new ApproximatorResilientPropagation copied from this instance
    ApproximatorResilientPropagation *ApproximatorResilientPropagation::Clone() const
    {
      return new ApproximatorResilientPropagation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ApproximatorResilientPropagation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ApproximatorResilientPropagation::GetAlias() const
    {
      static const std::string s_Name( "OpenCLResilientPropagation");
      return s_Name;
    }

    //! @brief set training data set for a specific iterate in approximater framework
    //! @param DATA training data set
    void ApproximatorResilientPropagation::SetTrainingData
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
    util::ShPtr< model::Interface> ApproximatorResilientPropagation::GetCurrentModel() const
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
    util::ShPtr< model::Interface> ApproximatorResilientPropagation::GetCurrentGPUModel() const
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
      ApproximatorResilientPropagation::GetCurrentApproximation() const
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
    void ApproximatorResilientPropagation::Next()
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

      Matrix< float> train_result( current_model->operator()( m_TrainingFeaturesOnDevice));
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

    //! read ApproximatorResilientPropagation from std::istream
    std::istream &ApproximatorResilientPropagation::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_TrainingData, ISTREAM);
      io::Serialize::Read( m_HiddenArchitecture, ISTREAM);
      io::Serialize::Read( m_Bias, ISTREAM);
      io::Serialize::Read( m_Weight, ISTREAM);

      // return
      return ISTREAM;
    }

    //! write ApproximatorResilientPropagation into std::ostream
    std::ostream &ApproximatorResilientPropagation::Write( std::ostream &OSTREAM, const size_t INDENT) const
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
    io::Serializer ApproximatorResilientPropagation::GetSerializer() const
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
      parameters.AddOptionalInitializer
      (
        "initial network",
        "provides an initial network model from which to start training from",
        io::Serialization::GetAgent( &m_NetworkFilename)
      );
      return parameters;
    }

    //! run forward through ANN and compute hidden terms m_Hidden (test)
    void ApproximatorResilientPropagation::CalcHiddenTerms()
    {
      // cpu equivalent
      //  // input layer
      //  m_HiddenInput( 0) = m_Weight( 0) * FEATURE + m_Bias( 0);
      //  m_Hidden( 0) = m_TransferFunction->F( m_HiddenInput( 0));
      //
      //  // remaining layers
      //  for( size_t k( 1); k != m_Hidden.GetSize(); ++k)
      //  {
      //    m_HiddenInput( k) = m_Weight( k) * m_Hidden( k - 1) + m_Bias( k);
      //    m_Hidden( k) = m_TransferFunction->F( m_HiddenInput( k));
      //  }

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

      BCL_MessageDbg( "hidden layer: " + util::Format()( m_HiddenBuffers( 0).GetHostMatrix()));

    }

    //! run backward through ANN and compute err terms m_Errors (train)
    void ApproximatorResilientPropagation::CalcErrorTerms()
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
    void ApproximatorResilientPropagation::CalcChangeTerms()
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
    void ApproximatorResilientPropagation::Train()
    {
      // training step
      CalcHiddenTerms();
      CalcErrorTerms();
      CalcChangeTerms();

    }

    void ApproximatorResilientPropagation::UpdateWeights()
    {
      for( size_t i( 0); i < m_HiddenBuffers.GetSize(); ++i)
      {
        PerformWeightsUpdate( i);
        PerformBiasUpdate( i);
      }
    }

    //! @brief performs the resilient update algorithm
    //! @param LAYER the layer number
    void ApproximatorResilientPropagation::PerformWeightsUpdate( const size_t &LAYER)
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;
      // dimensions
      const cl_uint rows(    m_WeightBuffers( LAYER).GetNumberRows());
      const cl_uint cols(    m_WeightBuffers( LAYER).GetNumberCols());
      const cl_uint row_pad( m_WeightBuffers( LAYER).GetRowPadding());
      const cl_uint col_pad( m_WeightBuffers( LAYER).GetColPadding());

      // create kernel
      cl::Kernel kernel( m_Program, "ResilientUpdate", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( s_Blocksize, s_Blocksize); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( s_Blocksize, cols), Tools::RoundUp( s_Blocksize, rows));

      error_number  = kernel.setArg( 0, m_WeightBuffers( LAYER).GetData());
      error_number |= kernel.setArg( 1, m_ChangeWeightBuffers( LAYER).GetData());
      error_number |= kernel.setArg( 2, m_SlopesWeightBuffers( LAYER).GetData());
      error_number |= kernel.setArg( 3, m_PrevSlopesWeightBuffers( LAYER).GetData());
      error_number |= kernel.setArg( 4, rows);
      error_number |= kernel.setArg( 5, cols);
      error_number |= kernel.setArg( 6, col_pad);
      error_number |= kernel.setArg( 7, row_pad);
      error_number |= kernel.setArg( 8, std::numeric_limits< float>::min());
      BCL_Assert( error_number == CL_SUCCESS, "resilient setarg error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

    }

    //! @brief performs the resilient update algorithm
    //! @param LAYER the layer number
    void ApproximatorResilientPropagation::PerformBiasUpdate( const size_t &LAYER)
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;
      // dimensions
      const cl_uint rows( 1);
      const cl_uint cols(    m_BiasBuffers( LAYER).GetSize());
      const cl_uint row_pad( 0);
      const cl_uint col_pad( m_BiasBuffers( LAYER).GetPadding());

      // create kernel
      cl::Kernel kernel( m_Program, "ResilientUpdate", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( s_Blocksize, s_Blocksize); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( s_Blocksize, cols), Tools::RoundUp( s_Blocksize, rows));

      error_number  = kernel.setArg( 0, m_BiasBuffers( LAYER).GetData());
      error_number |= kernel.setArg( 1, m_ChangeBiasBuffers( LAYER).GetData());
      error_number |= kernel.setArg( 2, m_SlopesBiasBuffers( LAYER).GetData());
      error_number |= kernel.setArg( 3, m_PrevSlopesBiasBuffers( LAYER).GetData());
      error_number |= kernel.setArg( 4, rows);
      error_number |= kernel.setArg( 5, cols);
      error_number |= kernel.setArg( 6, col_pad);
      error_number |= kernel.setArg( 7, row_pad);
      error_number |= kernel.setArg( 8, std::numeric_limits< float>::min());
      BCL_Assert( error_number == CL_SUCCESS, "resilient setarg error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
    }

    //! @brief calculates the output layer error terms
    void ApproximatorResilientPropagation::OutputLayerErrorTerms()
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
    }

    //! @brief calculates the hidden layer error terms
    //! @param LAYER the layer number
    void ApproximatorResilientPropagation::OtherLayerErrorTerms( const size_t &LAYER)
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

    //! @brief helper function to do a column-wise reduction on the bias errors
    //! @param ERRORS the errors buffer
    //! @param REDUCED_ERRORS the reduced error vector
    void ApproximatorResilientPropagation::ReduceErrorsToVector
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
//        linal::Matrix< float> errors_tmp( ROWS, COLS);
//        m_Queue.enqueueReadBuffer( ERRORS, CL_TRUE, 0, sizeof( float) * COLS * ROWS, errors_tmp.Begin());
//        BCL_MessageDbg( "bias to be reduced: " + util::Format()( errors_tmp));

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

      return;
    }

    //! @brief sets up the neural network with a particular architecture, after training data was set
    void ApproximatorResilientPropagation::SetupArchitecture()
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
        io::Serialize::Read( rmsd, in);
        BCL_MessageStd( "initial network had a score of " + util::Format()( rmsd));
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
        m_PrevSlopesBiasBuffers.PushBack( Vector< float>( linal::Vector< float>( itr->GetSize(), 0.0125), m_Queue, pad));
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
        m_PrevSlopesWeightBuffers.PushBack( Matrix< float>( linal::Matrix< float>( itr->GetNumberRows(), itr->GetNumberCols(), 0.0125), m_Queue, pad_row, pad_col));
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
    void ApproximatorResilientPropagation::UpdateQueue( Tools &TOOLS)
    {
      if( !TOOLS.HasCommandQueues())
      {
        return;
      }

      m_Queue = TOOLS.GetFirstCommandQueue();

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
#include "opencl/bcl_opencl_approximator_sequential_minimial_optimization.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_const_interface.hpp"
#include "model/bcl_model_support_vector_kernel_rbf.h"
#include "model/bcl_model_support_vector_machine.h"
#include "opencl/bcl_opencl_kernel_sources.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    const cl_uint ApproximatorSequentialMinimialOptimization::s_Lower_Bound( 0);    //!<  0 indicates lower boundary
    const cl_uint ApproximatorSequentialMinimialOptimization::s_Upper_Bound( 1);    //!<  1 indicates upper boundary
    const cl_uint ApproximatorSequentialMinimialOptimization::s_Free( 2);           //!<  2 indicates no boundary
    const float ApproximatorSequentialMinimialOptimization::m_EPS_A( 1e-12);
    const float ApproximatorSequentialMinimialOptimization::m_P( 0.1);
    const char *ApproximatorSequentialMinimialOptimization::s_CLCompilerOptions =
            "-cl-mad-enable -cl-fast-relaxed-math";

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ApproximatorSequentialMinimialOptimization::s_Instance
    (
      util::Enumerated< model::ApproximatorBase>::AddInstance( new ApproximatorSequentialMinimialOptimization())
    );

    //! @brief default constructor
    ApproximatorSequentialMinimialOptimization::ApproximatorSequentialMinimialOptimization() :
      m_CostParameterC( 0.0),
      m_Status(),
      m_Alpha(),
      m_Gradient(),
      m_GradientBar(),
      m_Bias(),
      m_Signs(),
      m_Labels(),
      m_ActiveSize(),
      m_ProbLength(),
      m_OptimizationGapThreshold( 0.1),
      m_OptimizationGap( 1.0),
      m_Model( new SupportVectorMachine()),
      m_NumberIterations( 0),
      m_NumberCurrentSupportVectors( 0)
    {
    }

    //! @brief Iterate for Sequential Minimal Optimization Learning Algorithm
    //! @param COST_PARAMETER_C penalty parameter c for svm regression training
    //! @param MODEL initial support vector model
    //! @param TRAINING_DATA training data set of choice
    //! @param NUMBER_ITERATIONS
    ApproximatorSequentialMinimialOptimization::ApproximatorSequentialMinimialOptimization
    (
      const float COST_PARAMETER_C,
      const util::ShPtr< SupportVectorMachine> &MODEL,
      util::ShPtr< descriptor::Dataset> &TRAINING_DATA,
      const size_t NUMBER_ITERATIONS,
      const CommandQueue &QUEUE
    ) :
      m_CostParameterC( COST_PARAMETER_C),
      m_Status(),
      m_Alpha(),
      m_Gradient(),
      m_GradientBar(),
      m_Bias(),
      m_Signs(),
      m_Labels(),
      m_ActiveSize( 0),
      m_ProbLength( 2 * TRAINING_DATA->GetSize()),
      m_Queue( QUEUE),
      m_OptimizationGapThreshold( 0.1),
      m_OptimizationGap( 1.0),
      m_Model( MODEL),
      m_NumberIterations( NUMBER_ITERATIONS),
      m_NumberCurrentSupportVectors( 0)
    {
      cl_int error_number( CL_SUCCESS);
      m_Program = KernelSources::Compile( GetKernelSources().e_SequentialMinimalOptimization, util::CPPDataTypes::e_Float, m_Queue, std::string( s_CLCompilerOptions), &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));

      SetTrainingData( TRAINING_DATA);
    }

    //! @brief copy constructor
    //! @return a new ApproximatorSequentialMinimialOptimization copied from this instance
    ApproximatorSequentialMinimialOptimization *ApproximatorSequentialMinimialOptimization::Clone() const
    {
      return new ApproximatorSequentialMinimialOptimization( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ApproximatorSequentialMinimialOptimization::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ApproximatorSequentialMinimialOptimization::GetAlias() const
    {
      static const std::string s_Name( "OpenclApproximatorSequentialMinimialOptimization");
      return s_Name;
    }

    //! @brief set training data set for a specific iterate in approximater framework
    //! @param DATA training data set
    void ApproximatorSequentialMinimialOptimization::SetTrainingData( util::ShPtr< descriptor::Dataset> &DATA)
    {
      m_TrainingData = DATA;
      DATA->GetFeatures().Rescale( model::SupportVectorMachine::s_DefaultInputRange);
      DATA->GetResults().Rescale( model::SupportVectorMachine::s_DefaultInputRange);

      // set rescale functions and kernel svm model used in iterative training process
      if( m_Model->GetNumberSupportVectors() == 0)
      {
        m_Model = util::ShPtr< SupportVectorMachine>
        (
          new SupportVectorMachine
          (
            0.0,
            linal::Vector< float>( 1),
            model::FeatureDataSet< float>( linal::Matrix< float>( 1, m_TrainingData->GetFeatureSize())),
            util::Implementation< model::SupportVectorKernelBase>( model::SupportVectorKernelRBF( m_Gamma)),
            *GetRescaleFeatureDataSet(),
            *GetRescaleResultDataSet(),
            m_Queue
          )
        );
      }

      m_TrainingFeaturesOnDevice = Matrix< float>( m_TrainingData->GetFeaturesPtr()->GetMatrix(), m_Queue);
      m_TrainingResultsOnDevice = Matrix< float>( m_TrainingData->GetResultsPtr()->GetMatrix(), m_Queue);

      // initialize model
      InitializeMemberVectorsForTraining();
    }

    //! @brief construct a model from the current iterate
    //! @return shptr to the new model interface
    util::ShPtr< model::Interface> ApproximatorSequentialMinimialOptimization::GetCurrentModel() const
    {
      return m_Model;
    }

    //! @brief returns the current approximation
    //! @return current argument result pair
    const util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> >
      ApproximatorSequentialMinimialOptimization::GetCurrentApproximation() const
    {
      util::ShPtr< model::Interface> model( GetCurrentModel());
      return
        util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> >
        (
          new storage::Pair< util::ShPtr< model::Interface>, float>
          (
            model.HardCopy(),
            m_OptimizationGap
          )
        );
    }

    //! @brief conducts the next approximation step and stores the approximation
    void ApproximatorSequentialMinimialOptimization::Next()
    {
      // optimization gap for indicating the error difference delta epsilon to a given epsilon
      // iterate for a number of internal iterations
      if( m_OptimizationGap > m_OptimizationGapThreshold && m_NumberCurrentSupportVectors < float( 0.9) * m_TrainingData->GetSize())
      {
        for( size_t counter( 0); counter < m_NumberIterations; ++counter)
        {
          m_OptimizationGap = IterationStep();
        }
      }

      // postprocess model and determine final support vectors
      FinalizeSupportVectorModel();

      m_NumberCurrentSupportVectors = m_Model->GetNumberSupportVectors();
      // create final pair with model and objective function evaluation
      util::ShPtr< storage::Pair< util::ShPtr< model::Interface>, float> > current_model
      (
        new storage::Pair< util::ShPtr< model::Interface>, float>( m_Model.HardCopy(), m_OptimizationGap)
      );
      this->GetTracker().Track( current_model);

      BCL_MessageStd
      (
        " #SV: " + util::Format()( m_Model->GetNumberSupportVectors())
        + " gap: " + util::Format()( m_OptimizationGap)
      );
    }

    //! @brief iterates one cycle and returns ShPtr to pair of resultant argument and corresponding score
    float ApproximatorSequentialMinimialOptimization::IterationStep()
    {
      size_t first_vector_index( 0);
      size_t second_vector_index( 0);

      // determine index i and j for the two data vectors for which the quadratic problem has to be solved
      const float optimization_gap
      (
        DetermineFeatureVectorCombination( first_vector_index, second_vector_index)
      );

      // Solving Quadratic Problem sub problems for two given data vectors
      SolveQuadraticProblemSubProblem( first_vector_index, second_vector_index);

      // return optimization difference to error threshold epsilon
      return optimization_gap;
    }

    //! @brief evaluates whether the approximation can continue
    //! @return true, if the approximation can continue - otherwise false
    bool ApproximatorSequentialMinimialOptimization::CanContinue() const
    {
      return m_OptimizationGap > m_OptimizationGapThreshold
             && m_NumberCurrentSupportVectors < float( 0.9) * m_TrainingData->GetSize();
    }

    //! @brief finalize support vector model and determine final SVs, Alphas and Bias
    void ApproximatorSequentialMinimialOptimization::FinalizeSupportVectorModel()
    {
      // reconstruct the whole gradient
      ReconstructGradient();

      // counter for iteration
      size_t vector_counter( 0);

      // final vector with all alphas of support vectors
      storage::Vector< float> alpha_final;

      storage::List< size_t> sv_indices;

      Vector< float> device_final_alphas( m_AlphaOnDevice.GetSize() / 2, m_Queue);
      Vector< cl_uint> device_sv_indices( m_AlphaOnDevice.GetSize() / 2, m_Queue);

      oclAssembleFinalAlphaVector( device_final_alphas, device_sv_indices);

      linal::Vector< float> tmp_alphas( device_final_alphas.GetHostVector());
      linal::Vector< cl_uint> tmp_sv_indices( device_sv_indices.GetHostVector());

      for( size_t count( 0), count_end( m_AlphaOnDevice.GetSize() / 2); count < count_end; ++count)
      {
        if( tmp_alphas( count) != 0)
        {
          alpha_final.PushBack( tmp_alphas( count));
          sv_indices.PushBack( tmp_sv_indices( count));
        }
      }

      // matrix containing support vectors
      linal::Matrix< float> sv_matrix
      (
        sv_indices.GetSize(),                             // rows
        m_TrainingData->GetFeatureSize(), // cols
        float( 0)                                  // default value
      );

      vector_counter = 0;

      // fill support vector matrix with vectors by support vector index
      for
      (
        storage::List< size_t>::const_iterator itr_sv( sv_indices.Begin()), itr_sv_end( sv_indices.End());
        itr_sv != itr_sv_end;
        ++itr_sv, ++vector_counter
      )
      {
        const model::FeatureReference< float> &support_vector( m_TrainingData->GetFeaturesPtr()->operator ()( *itr_sv));

        std::copy( support_vector.Begin(), support_vector.End(), sv_matrix[ vector_counter]);
      }

      // assemble SV Model
      m_Model->SetAlpha( linal::Vector< float>( alpha_final.Begin(), alpha_final.End()));
      m_Model->SetBias( CalculateBias());
      m_Model->SetSupportVectors( model::FeatureDataSet< float>( sv_matrix));
      m_Model->SetNumberSupportVectors( sv_matrix.GetNumberRows());
    }

    //! read NeuralNetwork from std::istream
    std::istream &ApproximatorSequentialMinimialOptimization::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_NumberIterations, ISTREAM);
      io::Serialize::Read( m_CostParameterC, ISTREAM);
      io::Serialize::Read( m_Status, ISTREAM);
      io::Serialize::Read( m_Alpha, ISTREAM);
      io::Serialize::Read( m_Gradient, ISTREAM);
      io::Serialize::Read( m_GradientBar, ISTREAM);
      io::Serialize::Read( m_Bias, ISTREAM);
      io::Serialize::Read( m_Signs, ISTREAM);
      io::Serialize::Read( m_Labels, ISTREAM);
      io::Serialize::Read( m_ActiveSize, ISTREAM);
      io::Serialize::Read( m_ProbLength, ISTREAM);
      io::Serialize::Read( m_Model, ISTREAM);
      io::Serialize::Read( m_TrainingData, ISTREAM);
      // return
      return ISTREAM;
    }

    //! write NeuralNetwork into std::ostream
    std::ostream &ApproximatorSequentialMinimialOptimization::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_NumberIterations, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_CostParameterC, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Status, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Alpha, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Gradient, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_GradientBar, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Bias, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Signs, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Labels, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ActiveSize, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ProbLength, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Model, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_TrainingData, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief method checks whether a LaGrange Multiplier reached a certain boundary or not
    //! @param ALPHA a LaGrange Multiplier
    //! @return const int that indicates whether ALPHA reached a certain boundary or not
    int ApproximatorSequentialMinimialOptimization::AlphaToStatus( const float &ALPHA) const
    {
      if( ALPHA >= m_CostParameterC - m_EPS_A)
      {
        return s_Upper_Bound;
      }
      else if( ALPHA <= m_EPS_A)
      {
        return s_Lower_Bound;
      }
      else
      {
        return s_Free;
      }
    }

    //! @brief
    void ApproximatorSequentialMinimialOptimization::UpdateAlphaStatus( const int &FEATURE_VECTOR_I)
    {
      if( m_Alpha( FEATURE_VECTOR_I) >= m_CostParameterC)
      {
        m_Status( FEATURE_VECTOR_I) = s_Upper_Bound;
      }
      else if( m_Alpha( FEATURE_VECTOR_I) <= 0)
      {
        m_Status( FEATURE_VECTOR_I) = s_Lower_Bound;
      }
      else
      {
        m_Status( FEATURE_VECTOR_I) = s_Free;
      }
    }

    //! @brief method checks whether a feature_vector i reached the upper boundary
    //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
    //! @return const bool that indicates whether FEATURE_VECTOR_I reached upper boundary
    inline bool ApproximatorSequentialMinimialOptimization::IsUpperBound( const size_t &FEATURE_VECTOR_I) const
    {
      return m_Status( FEATURE_VECTOR_I) == s_Upper_Bound;
    }

    //! @brief checks whether a feature_vector i reached the lower boundary
    //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
    //! @return const bool that indicates whether FEATURE_VECTOR_I reached lower boundary
    inline bool ApproximatorSequentialMinimialOptimization::IsLowerBound( const size_t &FEATURE_VECTOR_I) const
    {
      return m_Status( FEATURE_VECTOR_I) == s_Lower_Bound;
    }

    //! @brief checks whether a feature_vector i is not bound and between the boundaries
    //! @param FEATURE_VECTOR_I a position of a feature vector in m_Status
    //! @return const bool that indicates whether FEATURE_VECTOR_I reached no boundary
    inline bool ApproximatorSequentialMinimialOptimization::IsFree( const size_t &FEATURE_VECTOR_I) const
    {
      return m_Status( FEATURE_VECTOR_I) == s_Free;
    }

    //! @brief Initialize member vectors and variables to prepare SVR training
    //! @param TRAINING_DATA feature vector data set with labels
    //! @param SVR_MODEL model to be trained
    void ApproximatorSequentialMinimialOptimization::InitializeMemberVectorsForTraining()
    {
      // Initialization of variables and vectors
      const size_t number_feature_vectors( m_TrainingData->GetSize());

      // has to be set in case of this class is read from file and constructor was not explicitly applied
      m_ActiveSize = 0;
      m_ProbLength = 2 * number_feature_vectors;

      // initialize vector for lagrange multipliers
      m_Alpha = linal::Vector< float>( m_ProbLength, 0.0);

      // initialize vector for labeling of
      m_Labels = linal::Vector< float>( m_ProbLength, float( 1.0));
      m_Signs = linal::Vector< cl_int>( m_ProbLength, 1);

      for( size_t count( number_feature_vectors), count_end( m_ProbLength); count < count_end; ++count)
      {
        m_Signs( count) = -1;
        m_Labels( count) = -1.0;
      }

      // status indication of KKT determination for every vector
      m_Status = linal::Vector< cl_uint>( m_ProbLength);
      m_GradientBar = linal::Vector< float>( m_ProbLength, 0.0);

      // initialization of bias for support vector model
      m_Bias = linal::Vector< float>( m_ProbLength, m_P);

      // initializing
      for( size_t progress( 0); progress < number_feature_vectors; ++progress)
      {
        const float result( m_TrainingData->GetResultsPtr()->operator()( progress)( 0));
        m_Bias( progress) -= result;
        m_Bias( number_feature_vectors + progress) += result;
      }

      m_Gradient = m_Bias;

      m_Q_i                 = Vector< float>( m_ProbLength, m_Queue);
      m_Q_d                 = Vector< float>( m_ProbLength, m_Queue, 0, 1.0);
      m_Q_j                 = Vector< float>( m_ProbLength, m_Queue);
      m_GradientOnDevice    = Vector< float>( m_Gradient, m_Queue);
      m_BiasOnDevice        = Vector< float>( m_Bias, m_Queue);
      m_AlphaOnDevice       = Vector< float>( m_Alpha, m_Queue);
      m_LabelsOnDevice      = Vector< float>( m_Labels, m_Queue);
      m_SignsOnDevice       = Vector< cl_int>( m_Signs, m_Queue);
      m_GradientBarOnDevice = Vector< float>( m_GradientBar, m_Queue);
      m_StatusOnDevice      = Vector< cl_uint>( m_Status, m_Queue);

      // update alpha_status
      oclUpdateAllAlphaStatus();

      m_ActiveSize = m_ProbLength;

      Vector< float> kernel_vector( m_ProbLength / 2, m_Queue);

      // initializing Gradients
      for( size_t progress( 0); progress < m_ActiveSize; ++progress)
      {
        if( !IsLowerBound( progress))
        {
          oclGetInputIKernelMatrixResultingVector( progress, m_Gamma, m_Q_i, kernel_vector);

          oclInitializeGradient( m_Q_i);

          if( IsUpperBound( progress))
          {
            oclInitializeGradientBar( m_Q_i);
          }
        }
      }

      BCL_MessageDbg( " Initialization.. complete");
    }

    //! @brief determine the bias value for Support Vector Regression Model
    //! @return calculated bias value
    float ApproximatorSequentialMinimialOptimization::CalculateBias()
    {
      float bias;

      const cl_uint block_size( 128);
      const cl_uint num_groups( ( m_Alpha.GetSize() % block_size == 0 ? 0 : 1) + ( m_Alpha.GetSize() / block_size));
      Vector< float> upper_output( num_groups, m_Queue);
      Vector< float> lower_output( num_groups, m_Queue);
      Vector< cl_uint> nr_free_output( num_groups, m_Queue);
      Vector< float> sum_free_output( num_groups, m_Queue);
      float final_upper, final_lower, final_sum_free;
      cl_uint final_nr_free;
      oclCalculateBias( upper_output, lower_output, nr_free_output, sum_free_output, final_upper, final_lower, final_nr_free, final_sum_free);

      // calculate bias
      if( final_nr_free > 0)
      {
        bias = final_sum_free / final_nr_free;
      }
      else
      {
        bias = ( final_upper + final_lower) / 2;
      }

      // return calculated bias value
      return bias;
    }

    //! @brief heuristic approach to find a feature_vector i and j
    //! @brief depending on m_Gradients for examination of SMO Classification
    //! @param TRAINING_DATA data set of feature vectors inclusive labels
    //! @param SVR_MODEL
    //! @param FIRST_VECTOR_INDEX
    //! @param SECOND_VECTOR_INDEX
    //! @return pair of two indices for feature_vector i and j
    float ApproximatorSequentialMinimialOptimization::DetermineFeatureVectorCombination
    (
      size_t &FIRST_VECTOR_INDEX,
      size_t &SECOND_VECTOR_INDEX
    )
    {
      // return i,j such that
      // i: maximizes -y_i * grad(f)_i, i in I_up(\alpha)
      // j: minimizes the decrease of obj value
      //    (if quadratic coefficient <= 0, replace it with tau)
      //    -y_j*grad(f)_j < -y_i*grad(f)_i, j in I_low(\alpha)

      // input for kernel as std in opencl doesn't work
      float maximum_gradient( -std::numeric_limits< float>::infinity());
      cl_int maximum_gradient_index( -1);
      cl_int gmin_ind( -1);
      float gmax2( -std::numeric_limits< float>::infinity());

      oclFindMaxGradient( maximum_gradient, maximum_gradient_index);

      const cl_int index_i( maximum_gradient_index);

      Vector< float> kernel_vector( m_ProbLength / 2, m_Queue);

      if( index_i != -1) // null Q_i not accessed: Gmax=-INF if index_i=-1
      {
        oclGetInputIKernelMatrixResultingVector( index_i, m_Gamma, m_Q_i, kernel_vector);
      }

      const cl_uint block_size( 128);
      const cl_uint num_groups( ( m_Labels.GetSize() % block_size == 0 ? 0 : 1) + ( m_Labels.GetSize() / block_size));
      Vector< float> tmp_gmax( num_groups, m_Queue);
      Vector< float> tmp_obj_diff_min( num_groups, m_Queue);
      Vector< cl_uint> tmp_gmin_ind( num_groups, m_Queue);

      oclGetGminGmax2( m_Q_i, m_Q_d, index_i, maximum_gradient, tmp_gmax, tmp_obj_diff_min, tmp_gmin_ind, gmax2, gmin_ind);

      FIRST_VECTOR_INDEX  = maximum_gradient_index;
      SECOND_VECTOR_INDEX = gmin_ind;

      return maximum_gradient + gmax2; // optimization gap
    }

    //! @brief solve sub problem of quadratic problem according two feature vectors
    //! @param SVR_MODEL support vector machine model of interest
    //! @param FEATURE_VECTOR_I index of feature vector in MATRIX
    //! @param FEATURE_VECTOR_J index of feature vector in MATRIX
    //! @return a pair of vector of float with computed values of kernel function
    void ApproximatorSequentialMinimialOptimization::SolveQuadraticProblemSubProblem
    (
      const int &FEATURE_VECTOR_I,
      const int &FEATURE_VECTOR_J
    )
    {
      Vector< float> kernel_vector_i( m_Alpha.GetSize() / 2, m_Queue);

      oclGetInputIKernelMatrixResultingVector( FEATURE_VECTOR_I, m_Gamma, m_Q_i, kernel_vector_i);

      Vector< float> kernel_vector_j( m_Alpha.GetSize() / 2, m_Queue);
      oclGetInputIKernelMatrixResultingVector( FEATURE_VECTOR_J, m_Gamma, m_Q_j, kernel_vector_j);

      float quad_coef( m_Q_i.GetHostVector()( FEATURE_VECTOR_I) + m_Q_j.GetHostVector()( FEATURE_VECTOR_J) + 2 * m_Q_i.GetHostVector()( FEATURE_VECTOR_J));

      if( quad_coef <= 0)
      {
        quad_coef = m_EPS_A;
      }

      m_Gradient = m_GradientOnDevice.GetHostVector();
      m_Alpha    = m_AlphaOnDevice.GetHostVector();
      m_Labels   = m_LabelsOnDevice.GetHostVector();

      const float old_m_alpha_i( m_Alpha( FEATURE_VECTOR_I));
      const float old_m_alpha_j( m_Alpha( FEATURE_VECTOR_J));

      // if one of the two classLabels is negative
      // computing LaGrange Multipliers m_Alpha
      if( m_Labels( FEATURE_VECTOR_I) * m_Labels( FEATURE_VECTOR_J) < 0)
      {
        const float delta( ( -1 * m_Gradient( FEATURE_VECTOR_I) - m_Gradient( FEATURE_VECTOR_J)) / quad_coef);

        const float diff( m_Alpha( FEATURE_VECTOR_I) - m_Alpha( FEATURE_VECTOR_J));

        m_Alpha( FEATURE_VECTOR_I) += delta;
        m_Alpha( FEATURE_VECTOR_J) += delta;

        if( diff > 0)
        {
          if( m_Alpha( FEATURE_VECTOR_J) < 0)
          {
            m_Alpha( FEATURE_VECTOR_J) = 0;
            m_Alpha( FEATURE_VECTOR_I) = diff;
          }
          if( m_Alpha( FEATURE_VECTOR_I) > m_CostParameterC)
          {
            m_Alpha( FEATURE_VECTOR_I) = m_CostParameterC;
            m_Alpha( FEATURE_VECTOR_J) = m_CostParameterC - diff;
          }
        }
        else
        {
          if( m_Alpha( FEATURE_VECTOR_I) < 0)
          {
            m_Alpha( FEATURE_VECTOR_I) = 0;
            m_Alpha( FEATURE_VECTOR_J) = -diff;
          }
          if( m_Alpha( FEATURE_VECTOR_J) > m_CostParameterC)
          {
            m_Alpha( FEATURE_VECTOR_J) = m_CostParameterC;
            m_Alpha( FEATURE_VECTOR_I) = m_CostParameterC + diff;
          }
        }
      }
      else // if both classLabels are positive
      {
        const float delta( ( m_Gradient( FEATURE_VECTOR_I) - m_Gradient( FEATURE_VECTOR_J)) / quad_coef);
        const float sum( m_Alpha( FEATURE_VECTOR_I) + m_Alpha( FEATURE_VECTOR_J));

        m_Alpha( FEATURE_VECTOR_I) -= delta;
        m_Alpha( FEATURE_VECTOR_J) += delta;

        if( sum > m_CostParameterC)
        {
          if( m_Alpha( FEATURE_VECTOR_I) > m_CostParameterC)
          {
            m_Alpha( FEATURE_VECTOR_I) = m_CostParameterC;
            m_Alpha( FEATURE_VECTOR_J) = sum - m_CostParameterC;
          }
          if( m_Alpha( FEATURE_VECTOR_J) > m_CostParameterC)
          {
            m_Alpha( FEATURE_VECTOR_J) = m_CostParameterC;
            m_Alpha( FEATURE_VECTOR_I) = sum - m_CostParameterC;
          }
        }
        else
        {
          if( m_Alpha( FEATURE_VECTOR_J) < 0)
          {
            m_Alpha( FEATURE_VECTOR_J) = 0;
            m_Alpha( FEATURE_VECTOR_I) = sum;
          }
          if( m_Alpha( FEATURE_VECTOR_I) < 0)
          {
            m_Alpha( FEATURE_VECTOR_I) = 0;
            m_Alpha( FEATURE_VECTOR_J) = sum;
          }
        }
      }

      // update Gradient
      const float delta_alpha_i( m_Alpha( FEATURE_VECTOR_I) - old_m_alpha_i);
      const float delta_alpha_j( m_Alpha( FEATURE_VECTOR_J) - old_m_alpha_j);

      oclUpdateGradient( delta_alpha_i, delta_alpha_j, m_Q_i, m_Q_j);

      // update alpha_status and m_GradientBar
      m_Status = m_StatusOnDevice.GetHostVector();
      const bool feature_i_is_upper_bound( IsUpperBound( FEATURE_VECTOR_I));
      const bool feature_j_is_upper_bound( IsUpperBound( FEATURE_VECTOR_J));

      UpdateAlphaStatus( FEATURE_VECTOR_I);
      UpdateAlphaStatus( FEATURE_VECTOR_J);

      m_StatusOnDevice = Vector< cl_uint>( m_Status, m_Queue);
      m_AlphaOnDevice = Vector< float>( m_Alpha, m_Queue);

      short cost_multiplier( 0);

      if( feature_i_is_upper_bound != IsUpperBound( FEATURE_VECTOR_I))
      {
        cost_multiplier += ( feature_i_is_upper_bound ? -1 : 1);
      }
      if( feature_j_is_upper_bound != IsUpperBound( FEATURE_VECTOR_J))
      {
        cost_multiplier += ( feature_j_is_upper_bound ? -1 : 1);
      }

      if( cost_multiplier != 0)
      {
        const double effective_cost( cost_multiplier * m_CostParameterC);

        oclUpdateGradientBar( effective_cost, m_Q_i);
      }
    }

    void ApproximatorSequentialMinimialOptimization::ReconstructGradient()
    {
      // reconstruct inactive elements of m_Gradient from m_GradientBar and free variables
      if( m_ActiveSize == m_ProbLength)
      {
        return;
      }

      oclUpdateGradientWithBias();

      // vector of kernel values of vector i and all other vectors in training set
      Vector< float> kernel_vector( m_Alpha.GetSize() / 2, m_Queue);

      // progress counter accessing elements at a certain position in vector
      size_t progress( 0);

      // iterate over all LaGrange multipliers
      m_Alpha = m_AlphaOnDevice.GetHostVector();

      for
      (
        linal::Vector< float>::iterator
          iter_begin_alpha( m_Alpha.Begin()),
          iter_end_alpha( m_Alpha.End());
        iter_begin_alpha != iter_end_alpha;
        ++iter_begin_alpha, ++progress
      )
      {
        // if correspondent status values indicates that it is not between upper or lower bound
        if( IsFree( progress))
        {
          oclGetInputIKernelMatrixResultingVector( progress, m_Gamma, m_Q_i, kernel_vector);

          oclAddAlphaKernelToGradient( m_Q_i, *iter_begin_alpha);
        }
      }
    }

    void ApproximatorSequentialMinimialOptimization::oclFindMaxGradient( float &MAX_GRAD, cl_int &MAX_IND)
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( m_GradientOnDevice.GetSize());
      const cl_uint block_size( 128);
      const size_t num_groups( ( elements % block_size == 0 ? 0 : 1) + ( elements / block_size));

      Vector< float> tmp_max( num_groups, m_Queue);
      Vector< cl_uint> tmp_indexes( num_groups, m_Queue);

      // create kernel
      cl::Kernel kernel( m_Program, "FindMaxGradient", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, tmp_max.GetData());
      error_number |= kernel.setArg( 1, tmp_indexes.GetData());
      error_number |= kernel.setArg( 2, m_GradientOnDevice.GetData());
      error_number |= kernel.setArg( 3, m_LabelsOnDevice.GetData());
      error_number |= kernel.setArg( 4, m_StatusOnDevice.GetData());
      error_number |= kernel.setArg( 5, s_Upper_Bound);
      error_number |= kernel.setArg( 6, s_Lower_Bound);
      error_number |= kernel.setArg( 7, block_size * sizeof( float), 0);
      error_number |= kernel.setArg( 8, block_size * sizeof( cl_uint), 0);
      error_number |= kernel.setArg( 9, elements);
      BCL_Assert( error_number == CL_SUCCESS, "FindMaxGradient error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      linal::Vector< float> max_vector( tmp_max.GetHostVector());
      linal::Vector< cl_uint> ind_vector( tmp_indexes.GetHostVector());

      // complete on cpu
      float max_element( -std::numeric_limits< float>::infinity());
      size_t final_index( 0);
      for( size_t count( 0); count < num_groups; ++count)
      {
        size_t greater_than( max_vector( count) >= max_element ? 1 : 0);
        greater_than ? max_element = max_vector( count), final_index = ind_vector( count) : 0;
      }
      MAX_GRAD = max_element;
      MAX_IND  = final_index;
      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclGetInputIKernelMatrixResultingVector
    (
      const cl_uint &VECTOR_ID,
      const float &GAMMA,
      Vector< float> &OUTPUT,
      Vector< float> &KERNEL_VECTOR
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint rows( m_TrainingFeaturesOnDevice.GetNumberRows());
      const cl_uint cols( m_TrainingFeaturesOnDevice.GetNumberCols());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel_a( m_Program, "ComputeInputIKernelMatrix", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      cl::Kernel kernel_b( m_Program, "GetInputIKernelMatrixResultingVector", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      cl_uint index_vector;

      // adjust index for input vector i
      if( VECTOR_ID >= rows)
      {
        index_vector = VECTOR_ID - ( m_ProbLength / 2);
      }
      else
      {
        index_vector = VECTOR_ID;
      }

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize_a( Tools::RoundUp( block_size, rows));
      const cl::NDRange worksize_b( Tools::RoundUp( block_size, m_ProbLength));

      // add change in bias to previous bias
      error_number  = kernel_a.setArg( 0, index_vector);
      error_number |= kernel_a.setArg( 1, m_TrainingFeaturesOnDevice.GetData());
      error_number |= kernel_a.setArg( 2, cols);
      error_number |= kernel_a.setArg( 3, rows);
      error_number |= kernel_a.setArg( 4, KERNEL_VECTOR.GetData());
      error_number |= kernel_a.setArg( 5, block_size * sizeof( float), 0);
      error_number |= kernel_a.setArg( 6, GAMMA);
      BCL_Assert( error_number == CL_SUCCESS, "ComputeInputIKernelMatrix error: " + opencl::Tools::ErrorString( error_number));

      // add change in bias to previous bias
      error_number  = kernel_b.setArg( 0, VECTOR_ID);
      error_number |= kernel_b.setArg( 1, KERNEL_VECTOR.GetData());
      error_number |= kernel_b.setArg( 2, m_SignsOnDevice.GetData());
      error_number |= kernel_b.setArg( 3, cl_uint( m_ProbLength));
      error_number |= kernel_b.setArg( 4, OUTPUT.GetData());
      BCL_Assert( error_number == CL_SUCCESS, "GetInputIKernelMatrixResultingVector error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel_a, offset, worksize_a, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel_b, offset, worksize_b, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclGetGminGmax2
    (
      Vector< float> &Q_I,
      Vector< float> &Q_D,
      const cl_int &INDEX_I,
      float &MAX_GRADIENT,
      Vector< float> &TMP_GMAX2_OUTPUT,
      Vector< float> &TMP_OBJ_DIFF_MIN_OUTPUT,
      Vector< cl_uint> &TMP_GMIN_INDEX_OUTPUT,
      float &FINAL_GMAX2,
      cl_int       &FINAL_GMIN_IND
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( m_GradientOnDevice.GetSize());
      const cl_uint block_size( 128);
      const size_t num_groups( ( elements % block_size == 0 ? 0 : 1) + ( elements / block_size));

      // create kernel
      cl::Kernel kernel( m_Program, "GetGminGmax2", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, Q_I.GetData());
      error_number |= kernel.setArg( 1, Q_D.GetData());
      error_number |= kernel.setArg( 2, m_GradientOnDevice.GetData());
      error_number |= kernel.setArg( 3, m_LabelsOnDevice.GetData());
      error_number |= kernel.setArg( 4, m_StatusOnDevice.GetData());
      error_number |= kernel.setArg( 5, s_Upper_Bound);
      error_number |= kernel.setArg( 6, s_Lower_Bound);
      error_number |= kernel.setArg( 7, cl_uint( INDEX_I));
      error_number |= kernel.setArg( 8, MAX_GRADIENT);
      error_number |= kernel.setArg( 9, m_EPS_A);
      error_number |= kernel.setArg( 10, elements);
      error_number |= kernel.setArg( 11, TMP_GMAX2_OUTPUT.GetData());
      error_number |= kernel.setArg( 12, TMP_OBJ_DIFF_MIN_OUTPUT.GetData());
      error_number |= kernel.setArg( 13, TMP_GMIN_INDEX_OUTPUT.GetData());
      error_number |= kernel.setArg( 14, block_size * sizeof( float), 0);
      error_number |= kernel.setArg( 15, block_size * sizeof( float), 0);
      error_number |= kernel.setArg( 16, block_size * sizeof( cl_uint), 0);
      BCL_Assert( error_number == CL_SUCCESS, "GetGminGmax2 error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      linal::Vector< float> gmax2_vector( TMP_GMAX2_OUTPUT.GetHostVector());
      linal::Vector< float> obj_diff_min_vector( TMP_OBJ_DIFF_MIN_OUTPUT.GetHostVector());
      linal::Vector< cl_uint> gmin_ind_vector( TMP_GMIN_INDEX_OUTPUT.GetHostVector());

      // complete on cpu
      float gmax_element( -std::numeric_limits< float>::infinity());
      float obj_diff_min( std::numeric_limits< float>::infinity());
      int gmin_ind( 0);
      for( size_t count( 0); count < num_groups; ++count)
      {
        if( gmax2_vector( count) > gmax_element)
        {
          gmax_element = gmax2_vector( count);
        }

        if( obj_diff_min_vector( count) < obj_diff_min)
        {
          obj_diff_min = obj_diff_min_vector( count);
          gmin_ind = gmin_ind_vector( count);
        }
      }

      FINAL_GMAX2 = gmax_element;
      FINAL_GMIN_IND = gmin_ind;

      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclCalculateBias
    (
      Vector< float> &UPPER_OUTPUT,
      Vector< float> &LOWER_OUTPUT,
      Vector< cl_uint>      &NR_FREE_OUTPUT,
      Vector< float> &SUM_FREE_OUTPUT,
      float &FINAL_UPPER_OUTPUT,
      float &FINAL_LOWER_OUTPUT,
      cl_uint      &FINAL_NR_FREE_OUTPUT,
      float &FINAL_SUM_FREE_OUTPUT
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( m_StatusOnDevice.GetSize());
      const cl_uint block_size( 128);
      const size_t num_groups( ( elements % block_size == 0 ? 0 : 1) + ( elements / block_size));

      // create kernel
      cl::Kernel kernel( m_Program, "CalculateBias", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, m_StatusOnDevice.GetData());
      error_number |= kernel.setArg( 1, m_GradientOnDevice.GetData());
      error_number |= kernel.setArg( 2, m_LabelsOnDevice.GetData());
      error_number |= kernel.setArg( 3, s_Lower_Bound);
      error_number |= kernel.setArg( 4, s_Upper_Bound);
      error_number |= kernel.setArg( 5, elements);
      error_number |= kernel.setArg( 6, block_size * sizeof( float), 0);
      error_number |= kernel.setArg( 7, block_size * sizeof( float), 0);
      error_number |= kernel.setArg( 8, block_size * sizeof( cl_uint), 0);
      error_number |= kernel.setArg( 9, block_size * sizeof( float), 0);
      error_number |= kernel.setArg( 10, UPPER_OUTPUT.GetData());
      error_number |= kernel.setArg( 11, LOWER_OUTPUT.GetData());
      error_number |= kernel.setArg( 12, NR_FREE_OUTPUT.GetData());
      error_number |= kernel.setArg( 13, SUM_FREE_OUTPUT.GetData());
      BCL_Assert( error_number == CL_SUCCESS, "CalculateBias error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      linal::Vector< cl_uint> nr_free_vector( NR_FREE_OUTPUT.GetHostVector());
      linal::Vector< float> sum_free_vector( SUM_FREE_OUTPUT.GetHostVector());
      linal::Vector< float> upper_vector( UPPER_OUTPUT.GetHostVector());
      linal::Vector< float> lower_vector( LOWER_OUTPUT.GetHostVector());

      FINAL_NR_FREE_OUTPUT = nr_free_vector.Sum();
      FINAL_SUM_FREE_OUTPUT = sum_free_vector.Sum();

      // complete on cpu
      float final_upper( std::numeric_limits< float>::infinity());
      float final_lower( -std::numeric_limits< float>::infinity());
      for( size_t count( 0); count < num_groups; ++count)
      {
        if( upper_vector( count) < final_upper)
        {
          final_upper = upper_vector( count);
        }

        if( lower_vector( count) > final_lower)
        {
          final_lower = lower_vector( count);
        }
      }

      FINAL_LOWER_OUTPUT = final_lower;
      FINAL_UPPER_OUTPUT = final_upper;

      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclUpdateAllAlphaStatus
    (
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( m_AlphaOnDevice.GetSize());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel( m_Program, "UpdateAllAlphaStatus", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number |= kernel.setArg( 0, m_CostParameterC);
      error_number |= kernel.setArg( 1, s_Upper_Bound);
      error_number |= kernel.setArg( 2, s_Lower_Bound);
      error_number |= kernel.setArg( 3, s_Free);
      error_number |= kernel.setArg( 4, m_AlphaOnDevice.GetData());
      error_number |= kernel.setArg( 5, m_StatusOnDevice.GetData());
      error_number |= kernel.setArg( 6, elements);
      BCL_Assert( error_number == CL_SUCCESS, "UpdateAllAlphaStatus error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclUpdateGradient
    (
      const float &DELTA_ALPHA_I,
      const float &DELTA_ALPHA_J,
      const Vector< float> &Q_I,
      const Vector< float> &Q_J
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( m_GradientOnDevice.GetSize());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel( m_Program, "UpdateGradient", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, DELTA_ALPHA_I);
      error_number |= kernel.setArg( 1, DELTA_ALPHA_J);
      error_number |= kernel.setArg( 2, m_GradientOnDevice.GetData());
      error_number |= kernel.setArg( 3, Q_I.GetData());
      error_number |= kernel.setArg( 4, Q_J.GetData());
      error_number |= kernel.setArg( 5, elements);
      BCL_Assert( error_number == CL_SUCCESS, "UpdateGradient error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();

    }

    void ApproximatorSequentialMinimialOptimization::oclUpdateGradientBar
    (
      const float &EFFECTIVE_COST,
      const Vector< float> &Q_I
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( m_GradientBarOnDevice.GetSize());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel( m_Program, "UpdateGradientBar", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, EFFECTIVE_COST);
      error_number |= kernel.setArg( 1, m_GradientBarOnDevice.GetData());
      error_number |= kernel.setArg( 2, Q_I.GetData());
      error_number |= kernel.setArg( 3, elements);
      BCL_Assert( error_number == CL_SUCCESS, "UpdateGradientBar error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclUpdateGradientWithBias()
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( m_GradientBarOnDevice.GetSize());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel( m_Program, "UpdateGradientWithBias", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, m_GradientBarOnDevice.GetData());
      error_number |= kernel.setArg( 1, m_BiasOnDevice.GetData());
      error_number |= kernel.setArg( 2, elements);
      error_number |= kernel.setArg( 3, m_GradientOnDevice.GetData());
      BCL_Assert( error_number == CL_SUCCESS, "UpdateGradientBar error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclAddAlphaKernelToGradient
    (
      const Vector< float> &Q_I,
      const float &ALPHA
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( Q_I.GetSize());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel( m_Program, "AddAlphaKernelToGradient", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, cl_uint( m_ActiveSize));
      error_number |= kernel.setArg( 1, m_GradientOnDevice.GetData());
      error_number |= kernel.setArg( 2, Q_I.GetData());
      error_number |= kernel.setArg( 3, ALPHA);
      error_number |= kernel.setArg( 4, elements);
      BCL_Assert( error_number == CL_SUCCESS, "AddAlphaKernelToGradient error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclAssembleFinalAlphaVector
    (
      Vector< float> &FINAL_ALPHAS,
      Vector< cl_uint> &SV_INDECES
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( FINAL_ALPHAS.GetSize());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel( m_Program, "AssembleFinalAlphaVector", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, cl_uint( m_ProbLength / 2));
      error_number |= kernel.setArg( 1, m_AlphaOnDevice.GetData());
      error_number |= kernel.setArg( 2, FINAL_ALPHAS.GetData());
      error_number |= kernel.setArg( 3, SV_INDECES.GetData());
      BCL_Assert( error_number == CL_SUCCESS, "AssembleFinalAlphaVector error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclInitializeGradient
    (
      Vector< float> &Q_I
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( Q_I.GetSize());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel( m_Program, "InitializeGradient", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, m_AlphaOnDevice.GetData());
      error_number |= kernel.setArg( 1, Q_I.GetData());
      error_number |= kernel.setArg( 2, m_GradientOnDevice.GetData());
      error_number |= kernel.setArg( 3, elements);
      BCL_Assert( error_number == CL_SUCCESS, "InitializeGradient error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    void ApproximatorSequentialMinimialOptimization::oclInitializeGradientBar
    (
      Vector< float> &Q_I
    )
    {
      // catch errors
      cl_int error_number = CL_SUCCESS;

      // dimensions
      const cl_uint elements( Q_I.GetSize());

      const cl_uint block_size( 128);

      // create kernel
      cl::Kernel kernel( m_Program, "InitializeGradientBar", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

      // set thread block dimensions
      const cl::NDRange block_dimensions( block_size); //< all thread blocks have same dimensions
      const cl::NDRange offset;
      const cl::NDRange worksize( Tools::RoundUp( block_size, elements));

      // add change in bias to previous bias
      error_number  = kernel.setArg( 0, m_CostParameterC);
      error_number |= kernel.setArg( 1, Q_I.GetData());
      error_number |= kernel.setArg( 2, m_GradientBarOnDevice.GetData());
      error_number |= kernel.setArg( 3, elements);
      BCL_Assert( error_number == CL_SUCCESS, "InitializeGradientBar error: " + opencl::Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
      m_Queue.finish();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ApproximatorSequentialMinimialOptimization::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "trains a support vector machine using sequential-minimal-optimization "
        "(see http://en.wikipedia.org/wiki/Sequential_Minimal_Optimization)"
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
        "gamma",
        "currently hardcoded to RBF, so this is the gamma value",
        io::Serialization::GetAgent( &m_Gamma),
        "0.5"
      );

      parameters.AddInitializer
      (
        "iterations",
        "# of iterations used internally to improve the optimization gap",
        io::Serialization::GetAgent( &m_NumberIterations),
        "0"
      );
      parameters.AddInitializer
      (
        "cost",
        "controls the trade off between allowing training errors and forcing rigid margins; high values may lead to better training values, but risks overtraining",
        io::Serialization::GetAgent( &m_CostParameterC),
        "0.0"
      );
      return parameters;
    }

    //! @brief responsible for updating to a valid queue
    //! @param TOOLS opencl tools
    void ApproximatorSequentialMinimialOptimization::UpdateQueue( Tools &TOOLS)
    {
      if( !TOOLS.HasCommandQueues())
      {
        return;
      }

      m_Queue = GetTools().GetFirstCommandQueue();

      cl_int error_number( CL_SUCCESS);
      m_Program = KernelSources::Compile( GetKernelSources().e_SequentialMinimalOptimization, util::CPPDataTypes::e_Float, m_Queue, std::string( s_CLCompilerOptions), &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
    }

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_buffer.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Buffer::Buffer()
    {
    }

    //! @brief construct from cl::Buffer
    //! @param BUFFER the cl::Buffer
    Buffer::Buffer( const cl::Buffer &BUFFER) :
      cl::Buffer( BUFFER)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Buffer
    Buffer *Buffer::Clone() const
    {
      return new Buffer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Buffer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief size of memory in bytes
    //! @param ERROR_PTR error will be written to this location
    //! @return size_t size of buffer
    size_t Buffer::GetSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_MEM_SIZE>( ERROR_PTR);
    }

    //! @brief get the context for this buffer
    //! @param ERROR_PTR error will be written to this location
    //! @return Context the context
    Context Buffer::GetContext( cl_int *ERROR_PTR) const
    {
      return Context( getInfo< CL_MEM_CONTEXT>( ERROR_PTR));
    }

    //! @brief the associated host pointer
    //! @param ERROR_PTR error will be written to this location
    //! @return void* the host pointer if any is associated, NULL otherwise
    void *Buffer::GetHostPtr( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_MEM_HOST_PTR>( ERROR_PTR);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Buffer::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Buffer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_command_queue.h"

// includes from bcl - sorted alphabetically
#include "opencl/bcl_opencl_context.h"
#include "opencl/bcl_opencl_device.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CommandQueue::CommandQueue()
    {
    }

    //! @brief construct from context and device
    //! @param CONTEXT the context
    //! @param DEVICE the device for this command queue
    //! @param ERROR_PTR error will be written to this location
    CommandQueue::CommandQueue( const Context &CONTEXT, const Device &DEVICE, cl_int *ERROR_PTR) :
      cl::CommandQueue( CONTEXT, DEVICE, 0, ERROR_PTR)
    {
    }

    //! @brief construct from cl::CommandQueue
    //! @param COMMAND_QUEUE the cl::CommandQueue
    CommandQueue::CommandQueue( const cl::CommandQueue &COMMAND_QUEUE) :
      cl::CommandQueue( COMMAND_QUEUE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new CommandQueue
    CommandQueue *CommandQueue::Clone() const
    {
      return new CommandQueue( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CommandQueue::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get the context for this command queue
    //! @param ERROR_PTR error will be written to this location
    //! @return Context the context associated with this queue
    Context CommandQueue::GetContext( cl_int *ERROR_PTR) const
    {
      return Context( cl::CommandQueue::getInfo< CL_QUEUE_CONTEXT>( ERROR_PTR));
    }

    //! @brief get the devices for this command queue
    //! @param ERROR_PTR error will be written to this location
    //! @return Device the device associated with this queue
    Device CommandQueue::GetDevice( cl_int *ERROR_PTR) const
    {
      return Device( cl::CommandQueue::getInfo< CL_QUEUE_DEVICE>( ERROR_PTR));
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief comparison less than
    //! @param RHS right hand side queue
    //! @return true is this queue id is smaller than the rhs queue
    bool CommandQueue::operator <( const CommandQueue &RHS) const
    {
      return cl::CommandQueue::operator ()() < RHS();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CommandQueue::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &CommandQueue::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_context.h"

// includes from bcl - sorted alphabetically
#include "opencl/bcl_opencl_tools.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Context::Context()
    {
    }

    //! @brief construct from cl::Context
    //! @param CONTEXT the cl::Context to construct from
    Context::Context( const cl::Context &CONTEXT) :
      cl::Context( CONTEXT)
    {
    }

    //! @brief Context from devices
    //! @param ERROR_PTR pointer to error storage
    Context::Context( const storage::Vector< Device> &DEVICES, cl_int *ERROR_PTR)
    {
      // convert to cl::Devices
      std::vector< cl::Device> devices( DEVICES.Begin(), DEVICES.End());
      BCL_MessageVrb( "Creating context");
      cl_int error( CL_SUCCESS);
      cl::Context new_context( devices, NULL, NULL, NULL, &error);
      if( error == CL_SUCCESS)
      {
        BCL_MessageVrb( "Created new_context");
        cl::Context::operator =( new_context);
      }
      else
      {
        BCL_MessageCrt( "Error creating context: " + Tools::ErrorString( error));
        Tools::AssignError( ERROR_PTR, error);
      }
    }

    //! @brief Clone function
    //! @return pointer to new Context
    Context *Context::Clone() const
    {
      return new Context( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Context::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Context::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Context::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_class_descriptor.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_density_fit_protein_minimizer_powell.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_fit_protein_minimizer_powell.h"
#include "density/bcl_density_fit_protein_minimizers.h"
#include "density/bcl_density_map.h"
#include "opti/bcl_opti_approximator_powell.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_number_iterations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  ///////////
  // types //
  ///////////

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param RESOLUTION the resolution to simulate for
    DensityFitProteinMinimzerPowell::PositionCorrelation::PositionCorrelation()
    {
    }

    //! @brief Clone function
    //! @return pointer to new PositionCorrelation
    DensityFitProteinMinimzerPowell::PositionCorrelation *DensityFitProteinMinimzerPowell::PositionCorrelation::Clone() const
    {
      return new PositionCorrelation( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DensityFitProteinMinimzerPowell::PositionCorrelation::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the protein model
    //! @param PROTEIN_MODEL
    void DensityFitProteinMinimzerPowell::PositionCorrelation::SetProtein( const assemble::ProteinModel &PROTEIN_MODEL)
    {
      m_Atoms            = m_Simulator.AtomsToDevice( PROTEIN_MODEL.GetAtoms());
      m_NrAtoms          = m_Atoms.GetNumberRows();
    }

    //! @brief set the map
    //! @param DENSITY_MAP
    void DensityFitProteinMinimzerPowell::PositionCorrelation::SetDensityMap( const util::SiPtr< const density::Map> &DENSITY_MAP)
    {
      m_SpMap            = DENSITY_MAP;
      m_Dimensions       = m_Simulator.RoundUpDimensions( m_SpMap->GetDimensions());
      m_DensityMapBuffer = m_Correlation.MapToDevice( *m_SpMap, m_Dimensions);
    }

    //! @brief set the resolution
    //! @param RESOLUTION
    void DensityFitProteinMinimzerPowell::PositionCorrelation::SetResolution( const double RESOLUTION)
    {
      m_Simulator.SetResolution( RESOLUTION);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief is this class compatible with given command queue
    //! @param COMMAND_QUEUE the command queue this object would operate on
    //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
    bool DensityFitProteinMinimzerPowell::PositionCorrelation::IsCompatible( const CommandQueue &COMMAND_QUEUE) const
    {
      return m_Simulator.IsCompatible( COMMAND_QUEUE) &&
             m_Correlation.IsCompatible( COMMAND_QUEUE);
    }

    //! @brief initialize this class
    //! @brief COMMAND_QUEUE queue to use
    //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
    bool DensityFitProteinMinimzerPowell::PositionCorrelation::Initialize( const CommandQueue &COMMAND_QUEUE)
    {
      m_CommandQueue = COMMAND_QUEUE;
      m_MinMax = DataSetMinMax< double>( m_CommandQueue);
      // initialize opencl classes
      return m_Simulator.Initialize( m_CommandQueue) &&
             m_Correlation.Initialize( m_CommandQueue) &&
             m_Transformer.Initialize( m_CommandQueue);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief calculate the correlation
    //! @param VECTOR containing 6 elements - three rotations, three translations
    double DensityFitProteinMinimzerPowell::PositionCorrelation::operator()( const linal::Vector< double> &VECTOR) const
    {
      const math::TransformationMatrix3D transformation( VECTOR);

      // transform the coordinates
      const Matrix< double> atoms( m_Transformer( m_Atoms, transformation, m_NrAtoms));

      linal::Vector3D mincoord( m_MinMax.Min( atoms).Begin());
      linal::Vector3D maxcoord( m_MinMax.Max( atoms).Begin());

      // add margin
      mincoord -= 2 * m_SpMap->GetCellWidth();
      maxcoord += 2 * m_SpMap->GetCellWidth();

      // determine index
      const storage::VectorND< 3, int> index
      (
        int( std::floor( mincoord.X() / m_SpMap->GetCellWidth().X())),
        int( std::floor( mincoord.Y() / m_SpMap->GetCellWidth().Y())),
        int( std::floor( mincoord.Z() / m_SpMap->GetCellWidth().Z()))
      );

      // dimensions of grid
      const storage::VectorND< 3, size_t> exact_dimensions
      (
        size_t( std::ceil( maxcoord.X() / m_SpMap->GetCellWidth().X())) - index.First() ,
        size_t( std::ceil( maxcoord.Y() / m_SpMap->GetCellWidth().Y())) - index.Second(),
        size_t( std::ceil( maxcoord.Z() / m_SpMap->GetCellWidth().Z())) - index.Third()
      );
      // dimensions of grid
      const storage::VectorND< 3, size_t> dimensions( m_Simulator.RoundUpDimensions( exact_dimensions));

      // relative index of argument density relative to this
      linal::Vector< int> this_start( 3, int( 0));
      linal::Vector< int> arg_start(  3, int( 0));
      linal::Vector< int> this_end(   3, int( 0));
      linal::Vector< int> arg_end(    3, int( 0));
      // find common start and end
      linal::Vector< int> common_start( 3, int( 0));
      linal::Vector< int> common_end( 3, int( 0));
      linal::Vector< int> extent( 3, int( 0));

      // iterate over dimensions
      for( size_t i( 0); i < 3; ++i)
      {
        this_start( i) = m_SpMap->GetIndex()( i);
        arg_start(  i) =               index( i);
        const int this_end      = this_start( i) + m_Dimensions( i);
        const int arg_end       = arg_start(  i) + dimensions( i);
        const int common_start  = std::max( this_start( i), arg_start( i));
        const int common_end    = std::min( this_end      , arg_end);
        extent( i) = common_end - common_start - 1;

        // no common sub density
        if( extent( i) <= 0)
        {
          return 0.0;
        }

        this_start( i) = common_start - m_SpMap->GetIndex()( i);
        arg_start(  i) = common_start -               index( i);
      }

      // buffer for the density map
      const Vector< double> density( m_Simulator.Simulate( atoms, m_NrAtoms, index, dimensions));

      // calculate cross correlation
      const double ccc
      (
        -m_Correlation.CrossCorrelationCoefficient
        (
          m_DensityMapBuffer,
          density,
          this_start,
          m_Dimensions,
          arg_start,
          dimensions,
          extent,
          0.0
        )
      );

      // end
      return ccc;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DensityFitProteinMinimzerPowell::PositionCorrelation::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DensityFitProteinMinimzerPowell::PositionCorrelation::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> DensityFitProteinMinimzerPowell::s_Instance
    (
      GetObjectInstances().AddInstance( new DensityFitProteinMinimzerPowell())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DensityFitProteinMinimzerPowell::DensityFitProteinMinimzerPowell() :
      m_CorrelationFunction( new PositionCorrelation())
    {
    }

    //! @brief Clone function
    //! @return pointer to new DensityFitProteinMinimzerPowell
    DensityFitProteinMinimzerPowell *DensityFitProteinMinimzerPowell::Clone() const
    {
      return new DensityFitProteinMinimzerPowell( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DensityFitProteinMinimzerPowell::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the resolution of the density map
    //! @param RESOLUTION density map and simulation resolution
    void DensityFitProteinMinimzerPowell::SetResolution( const double RESOLUTION)
    {
      m_CorrelationFunction->SetResolution( RESOLUTION);
    }

    //! @brief set max translation and rotation
    //! @param MAX_TRANSLATION max translation in any direction for a single iteration
    //! @param MAX_ROTATION max rotation in radians in any direction for a single iteration
    void DensityFitProteinMinimzerPowell::SetMaxTranslationAndRotation( const double MAX_TRANSLATION, const double MAX_ROTATION)
    {
      m_MaxTranslation = MAX_TRANSLATION;
      m_MaxRotation    = MAX_ROTATION;
    }

    //! @brief set the max number of iterations for minimization
    //! @param MAX_NUMBER_ITERATIONS maximum number of iterations for minimization
    void DensityFitProteinMinimzerPowell::SetMaxIterations( const size_t MAX_NUMBER_ITERATIONS)
    {
      m_MaxIterations = MAX_NUMBER_ITERATIONS;
    }

    //! @brief set protein agreement measure to be used
    //! @param AGREEMENT protein agreement enumerator
    void DensityFitProteinMinimzerPowell::SetProteinAgreement( const density::ProteinAgreement &AGREEMENT)
    {
      if( AGREEMENT != density::GetProteinAgreements().e_CCC)
      {
        BCL_MessageStd( "currently, only CCC as ProteinAgreement is supported via Opencl");
      }
    }

    //! @brief simulator to use
    //! @param DENSITY_SIMULATOR simulator enumerator
    void DensityFitProteinMinimzerPowell::SetSimulator( const density::Simulator &DENSITY_SIMULATOR)
    {
      if( DENSITY_SIMULATOR != density::GetSimulators().e_GaussianSphere)
      {
        BCL_MessageStd( "currently, only GaussianSphere as density simulator is supported via Opencl");
      }
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief is this class compatible with given command queue
    //! @param COMMAND_QUEUE the command queue this object would operate on
    //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
    bool DensityFitProteinMinimzerPowell::IsCompatible( const CommandQueue &COMMAND_QUEUE) const
    {
      return m_CorrelationFunction->IsCompatible( COMMAND_QUEUE);
    }

    //! @brief initialize this class
    //! @brief COMMAND_QUEUE queue to use
    //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
    bool DensityFitProteinMinimzerPowell::Initialize( const CommandQueue &COMMAND_QUEUE)
    {
      return m_CorrelationFunction->Initialize( COMMAND_QUEUE);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator minimizing the position of a protein model within a given density map
    //! @param PROTEIN_MODEL start position of given protein model
    //! @param DENSITY_MAP the density map to fit the PROTEIN_MODEL into
    //! @return the fitted protein model
    assemble::ProteinModel DensityFitProteinMinimzerPowell::operator()( const assemble::ProteinModel &PROTEIN_MODEL, const density::Map &DENSITY_MAP) const
    {
      // set members of correlation function
      m_CorrelationFunction->SetProtein( PROTEIN_MODEL);
      m_CorrelationFunction->SetDensityMap( util::ToSiPtr( DENSITY_MAP));

      storage::Vector< linal::Vector< double> > search_directions;
      linal::Vector< double> direction( 6, 0.0);
      const linal::Vector< double> start( direction);
      direction( 0) = m_MaxRotation;
      search_directions.PushBack( direction);
      direction( 0) = 0.0;
      direction( 1) = m_MaxRotation;
      search_directions.PushBack( direction);
      direction( 1) = 0.0;
      direction( 2) = m_MaxRotation * 0.5;
      search_directions.PushBack( direction);
      direction( 2) = 0.0;
      direction( 3) = m_MaxTranslation;
      search_directions.PushBack( direction);
      direction( 3) = 0.0;
      direction( 4) = m_MaxTranslation;
      search_directions.PushBack( direction);
      direction( 4) = 0.0;
      direction( 5) = m_MaxTranslation;
      search_directions.PushBack( direction);

      // create termination criteria for the approximation
      opti::CriterionCombine< linal::Vector< double>, double> criterion_combine;
      criterion_combine.InsertCriteria
      (
        opti::CriterionNumberIterations< linal::Vector< double>, double>( m_MaxIterations / 10)
      );
      criterion_combine.InsertCriteria
      (
        opti::CriterionConvergenceResult< linal::Vector< double>, double>( 1, 0.001)
      );

      // create powell approximator from its members
      opti::ApproximatorPowell< linal::Vector< double>, double> approximator
      (
        *m_CorrelationFunction, criterion_combine, search_directions, start
      );

      // do the actual approximation
      approximator.Approximate();

      // create the final model
      const util::ShPtr< assemble::ProteinModel> best_model
      (
        density::FitProteinMinimizerPowell::PositionCorrelation::TransformedHardCopy
        (
          PROTEIN_MODEL, approximator.GetTracker().GetBest()->First()
        )
      );

      return *best_model;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DensityFitProteinMinimzerPowell::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DensityFitProteinMinimzerPowell::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DensityFitProteinMinimzerPowellEnumHandler
    //! @brief handler class for adding the operations enum handler
    //! @author woetzen, loweew
    //! @date Mar 14, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DensityFitProteinMinimzerPowellEnumHandler :
      public signal::Slots
    {

    public:

    //////////
    // data //
    //////////

      //! the enum in the density::FitProteinMinimizers
      density::FitProteinMinimizer e_Minimzer;

      //! the only instance of this class
      static const DensityFitProteinMinimzerPowellEnumHandler s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    private:

      //! @brief default constructor
      DensityFitProteinMinimzerPowellEnumHandler();

    public:

    ////////////////
    // operations //
    ////////////////

      //! @brief update the enum with the command queue from the Tools
      //! @param TOOLS the tolls to get the commandqueue from
      void UpdateEnum( Tools &TOOLS);

    }; // template class DensityFitProteinMinimzerPowellEnumHandler

    //! instance of DensityFitProteinMinimzerPowellEnumHandler
    const DensityFitProteinMinimzerPowellEnumHandler DensityFitProteinMinimzerPowellEnumHandler::s_Instance = DensityFitProteinMinimzerPowellEnumHandler();

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    DensityFitProteinMinimzerPowellEnumHandler::DensityFitProteinMinimzerPowellEnumHandler() :
      e_Minimzer( density::GetFitProteinMinimizers().AddEnum( "PowellOpencl", util::ShPtr< DensityFitProteinMinimzerPowell>()))
    {
      GetTools().GetQueueUpdateSignal().Connect( this, &DensityFitProteinMinimzerPowellEnumHandler::UpdateEnum);
    }

    //! @brief update the enum with the command queue from the Tools
    //! @param TOOLS the tolls to get the commandqueue from
    void DensityFitProteinMinimzerPowellEnumHandler::UpdateEnum( Tools &TOOLS)
    {
      util::ShPtr< DensityFitProteinMinimzerPowell> sp_minimizer( new DensityFitProteinMinimzerPowell());
      if( !TOOLS.HasCommandQueues())
      {
        *e_Minimzer = util::ShPtr< DensityFitProteinMinimzerPowell>();
        return;
      }

      // try to initialize
      if( sp_minimizer->Initialize( TOOLS.GetFirstCommandQueue()))
      {
        // just update the existing one with the new one
        *e_Minimzer = sp_minimizer;
      }
      else
      {
        BCL_MessageVrb( "unable to initialize enum: OpenCL");
      }
    }

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_density_simulate_gaussian_sphere.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "density/bcl_density_simulators.h"
#include "opencl/bcl_opencl_dataset_min_max.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

    //! @class DensitySimulateEnumHandler
    //! @brief handler class for adding the density simulate enum handler
    class BCL_API DensitySimulateEnumHandler :
      public signal::Slots
    {

    private:

    //////////
    // data //
    //////////

      //! the enum in the density::Simulators
      density::Simulator e_DensitySimulateGaussianSphere;

      //! the only instance of this class
      static const DensitySimulateEnumHandler s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DensitySimulateEnumHandler() :
        e_DensitySimulateGaussianSphere( density::GetSimulators().AddEnum( "OpenclGaussianSphere", util::ShPtr< DensitySimulateGaussianSphere>()))
      {
        // register enum with opencl queue update signal
        GetTools().GetQueueUpdateSignal().Connect( this, &DensitySimulateEnumHandler::UpdateEnum);
      }

      //! @brief update the enum with the command queue from the Tools
      //! @param TOOLS the tolls to get the commandqueue from
      void UpdateEnum( Tools &TOOLS)
      {
        util::ShPtr< DensitySimulateGaussianSphere> sp_simulator( new DensitySimulateGaussianSphere());
        if( !TOOLS.HasCommandQueues())
        {
          *e_DensitySimulateGaussianSphere = util::ShPtr< DensitySimulateGaussianSphere>();
          return;
        }

        // try to initialize
        if( sp_simulator->Initialize( TOOLS.GetFirstCommandQueue()))
        {
          // just update the existing one with the new one
          *e_DensitySimulateGaussianSphere = sp_simulator;
        }
        else
        {
          BCL_MessageVrb( "unable to initialize enum: OpenclGaussianSphere");
        }
      }

    }; // class DensitySimulateEnumHandler

    //! instance of DensitySimulateEnumHandler
    const DensitySimulateEnumHandler DensitySimulateEnumHandler::s_Instance = DensitySimulateEnumHandler();

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from grid spacing and resolution
    //! @param GRID_SPACING the spacing for the density grid
    //! @param RESOLUTION the resolution to simulate for
    DensitySimulateGaussianSphere::DensitySimulateGaussianSphere() :
      m_GridSpacing( density::Simulators::GetDefaultGridSpacing()),
      m_Resolution( density::Simulators::GetDefaultResolution()),
      m_Margin( 2),
      m_CommandQueue()
    {
    }

    //! @brief Clone function
    //! @return pointer to new DensitySimulateGaussianSphere
    DensitySimulateGaussianSphere *DensitySimulateGaussianSphere::Clone() const
    {
      return new DensitySimulateGaussianSphere( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DensitySimulateGaussianSphere::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set the resolution
    //! @param RESOLUTION the resolution for the density map to be generated
    void DensitySimulateGaussianSphere::SetResolution( const double RESOLUTION)
    {
      m_Resolution = RESOLUTION;
    }

    //! @brief set the resolution
    double DensitySimulateGaussianSphere::GetResolution() const
    {
      return m_Resolution;
    }

    //! @brief set the grid spacing
    //! @param GRID_SPACING the width of a grid element in x, y and z
    void DensitySimulateGaussianSphere::SetGridSpacing( const linal::Vector3D &GRID_SPACING)
    {
      m_GridSpacing = GRID_SPACING;
    }

    //! @brief set the margin
    //! @param MARGIN number of additional cells next to last atom occupied cells
    void DensitySimulateGaussianSphere::SetMargin( const size_t MARGIN)
    {
      m_Margin = MARGIN;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief is this class compatible with given command queue
    //! @param COMMAND_QUEUE the command queue this object would operate on
    //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
    bool DensitySimulateGaussianSphere::IsCompatible( const CommandQueue &COMMAND_QUEUE) const
    {
      cl_int error_number( CL_SUCCESS);
      const Device device( COMMAND_QUEUE.GetDevice( &error_number));

      // can get device
      if( error_number != CL_SUCCESS)
      {
        return false;
      }

      const storage::Set< Extension> extensions( device.Extensions( &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "unable to get extensions from device");

      return KernelSourceInterface::PrecisionCompatibleWithExtensions
             (
               util::CPPDataTypes::DataTypeFromTemplate< double>(),
               extensions
             );
    }

    //! @brief initialize this class
    //! @brief COMMAND_QUEUE queue to use
    //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
    bool DensitySimulateGaussianSphere::Initialize( const CommandQueue &COMMAND_QUEUE)
    {
      // check if this is a compatible command queue
      if( !IsCompatible( COMMAND_QUEUE))
      {
        BCL_MessageDbg( "command queue is not compatible");
        return false;
      }

      // update the command queue
      m_CommandQueue = COMMAND_QUEUE;

      // for precision type
      const cl_int error_number = CompilePrograms( util::CPPDataTypes::DataTypeFromTemplate< double>());
      if( error_number != CL_SUCCESS)
      {
        BCL_MessageDbg( "error compiling programs:\n" + Tools::ErrorString( error_number));
        return false;
      }

      return true;
    }

    //! @brief simulate a density map for a buffer of atoms into given grid dimensions
    //! @param ATOMS a Buffer of atoms with their weight
    //! @param NR_ATOMS number of atoms in Buffer with size 4 * NR_ATOMS
    //! @param INDEX the index of the grid
    //! @param DIMENSION the dimension of the grid
    Vector< double> DensitySimulateGaussianSphere::Simulate
    (
      const Matrix< double> &ATOMS,
      const size_t NR_ATOMS,
      const storage::VectorND< 3, int> &INDEX,
      const storage::VectorND< 3, size_t> &DIMENSIONS
    ) const
    {
      cl_int error_number( CL_SUCCESS);

      BCL_Assert
      (
        DIMENSIONS.First() % 4 == 0 && DIMENSIONS.Second() % 4 == 0 && DIMENSIONS.Third() % 4 == 0,
        "dimensions need to be a multiple of 4"
      );

      BCL_Assert( NR_ATOMS % 64 == 0, "NR_ATOMS need to be a multiple of 64");

      // constants describing gaussian blob shape
      const double blob_k( math::Sqr( math::g_Pi / ( 2.4 + 0.8 * m_Resolution)));

      // constant for square distance cutoff
      const double cutoff_square( math::Sqr( 3 * ( 1.0 / math::Sqrt( 2)) * ( ( 2.4 + 0.8 * m_Resolution) / math::g_Pi)));

      // number of elements in grid
      const size_t grid_size( DIMENSIONS.First() * DIMENSIONS.Second() * DIMENSIONS.Third());

      // store the real space index for easier access
      const linal::Vector3D realspaceindex
        (
          INDEX.First()  * m_GridSpacing.X(),
          INDEX.Second() * m_GridSpacing.Y(),
          INDEX.Third()  * m_GridSpacing.Z()
        );

      // flat array for the density map
      Vector< double> device_grid( grid_size, m_CommandQueue);

      cl::Kernel kernel( m_Program, "SimulateDensityGaussianSphere", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number  = kernel.setArg(  0, ATOMS.GetData());
      error_number |= kernel.setArg(  1, cl_uint( NR_ATOMS));
      error_number |= kernel.setArg(  2, realspaceindex.X());
      error_number |= kernel.setArg(  3, realspaceindex.Y());
      error_number |= kernel.setArg(  4, realspaceindex.Z());
      error_number |= kernel.setArg(  5, m_GridSpacing.X());
      error_number |= kernel.setArg(  6, m_GridSpacing.Y());
      error_number |= kernel.setArg(  7, m_GridSpacing.Z());
      error_number |= kernel.setArg(  8, blob_k);
      error_number |= kernel.setArg(  9, cutoff_square);
      error_number |= kernel.setArg( 10, device_grid.GetData());
      error_number |= kernel.setArg( 11, 64 * 4 * sizeof( double), 0);
      BCL_Assert( error_number == CL_SUCCESS, "setting gpu kernel args error: " + Tools::ErrorString( error_number));

      const cl::NDRange local_worksize( 4, 4, 4);
      const cl::NDRange offset;
      const cl::NDRange global_worksize( DIMENSIONS.First(), DIMENSIONS.Second(), DIMENSIONS.Third());

      // launching kernel
      error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, offset, global_worksize, local_worksize);
      BCL_Assert( error_number == CL_SUCCESS, "kernel launch error: " + Tools::ErrorString( error_number));

      // end
      return device_grid;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief generate simulated density from given list of atoms
    //! @param ATOMS siptrvector of atoms
    //! @return a simulated density map
    density::Map DensitySimulateGaussianSphere::operator()( const util::SiPtrVector< const biol::Atom> &ATOMS) const
    {
      cl_int error_number( CL_SUCCESS);

      // copy atoms to device
      Matrix< double> device_matrix( AtomsToDevice( ATOMS));

      DataSetMinMax< double> min_max( m_CommandQueue);
      const int number_rows( Tools::RoundUp( 64, ATOMS.GetSize()));

      linal::Vector3D mincoord( min_max.Min( device_matrix).Begin());
      linal::Vector3D maxcoord( min_max.Max( device_matrix).Begin());

      // add margin
      mincoord -= m_Margin * m_GridSpacing;
      maxcoord += m_Margin * m_GridSpacing;

      // determine index
      const linal::VectorND< int, 3> index
      (
        int( std::floor( mincoord.X() / m_GridSpacing.X())),
        int( std::floor( mincoord.Y() / m_GridSpacing.Y())),
        int( std::floor( mincoord.Z() / m_GridSpacing.Z()))
      );

      // dimensions of grid
      const storage::VectorND< 3, size_t> exact_dimensions
      (
        size_t( std::ceil( maxcoord.X() / m_GridSpacing.X())) - index( 0) ,
        size_t( std::ceil( maxcoord.Y() / m_GridSpacing.Y())) - index( 1),
        size_t( std::ceil( maxcoord.Z() / m_GridSpacing.Z())) - index( 2)
      );
      // dimensions of grid
      const storage::VectorND< 3, size_t> dimensions( RoundUpDimensions( exact_dimensions));

      const storage::VectorND< 3, int> index_nd
        (
         int( std::floor( mincoord.X() / m_GridSpacing.X())),
         int( std::floor( mincoord.Y() / m_GridSpacing.Y())),
         int( std::floor( mincoord.Z() / m_GridSpacing.Z()))
         );

      // flat array (vector) for the density map since we don't currently have an opencl::Tensor class
      Vector< double> device_grid( Simulate( device_matrix, number_rows, index_nd, dimensions));

      // allocate mask size and set values to 0
      math::Tensor< double> grid
      (
        dimensions.Third(), dimensions.Second(), dimensions.First(), double( 0.0)
      );

      // read result
      error_number = m_CommandQueue.enqueueReadBuffer( device_grid.GetData(), CL_TRUE, 0, sizeof( double) * grid.GetSize(), grid.Begin());
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // intervals
      linal::VectorND< int, 3> intervals;
      for( size_t i( 0); i < 3; ++i)
      {
        if( index( i) <= -int( index( i)))
        {
          intervals( i) = math::Absolute( index( i));
        }
        else if( index( i) >= 0)
        {
          intervals( i) = index( i) + index( i) - 1;
        }
        else
        {
          intervals( i) = index( i) - 1;
        }
      }

      // length
      const linal::Vector3D length
      (
        m_GridSpacing.X() * double( intervals( 0)),
        m_GridSpacing.Y() * double( intervals( 1)),
        m_GridSpacing.Z() * double( intervals( 2))
      );

      // end
      return density::Map
             (
               grid,
               index,
               intervals,
               length,
               m_GridSpacing,
               density::Map::GetDefaultAngle(),
               density::Map::GetDefaultAxis(),
               linal::Vector3D( 0.0)
             );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DensitySimulateGaussianSphere::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DensitySimulateGaussianSphere::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief atoms to device buffer
    //! @param ATOMS siptrvector of atoms
    //! @return Buffer a buffer on the device that is associated with this command queue
    Matrix< double> DensitySimulateGaussianSphere::AtomsToDevice( const util::SiPtrVector< const biol::Atom> &ATOMS) const
    {
      // constants describing gaussian blob shape
      const double blob_k( math::Sqr( math::g_Pi / ( 2.4 + 0.8 * m_Resolution)));
      const double blob_c( math::Pow( blob_k / math::g_Pi, 1.5));

      const linal::Matrix< double> atoms_matrix( AtomsToMatrix( ATOMS, blob_c));

      Matrix< double> device_matrix( atoms_matrix, m_CommandQueue);

      // end
      return device_matrix;
    }

    //! @brief round up dimensions so that they are compatible with work group size
    //! @param DIMENSIONS the original dimensions
    //! @return storage::VectorND< 3, size_t> the new diemsions rounded up
    storage::VectorND< 3, size_t> DensitySimulateGaussianSphere::RoundUpDimensions
    (
      const storage::VectorND< 3, size_t> &DIMENSIONS
    ) const
    {
      return storage::VectorND< 3, size_t>
      (
        Tools::RoundUp( 4, DIMENSIONS.First()) ,
        Tools::RoundUp( 4, DIMENSIONS.Second()),
        Tools::RoundUp( 4, DIMENSIONS.Third())
      );
    }

    //! @brief convert atoms to padded matrix
    //! @param ATOMS siptrvector of atoms
    //! @return linal::Matrix< double> a matrix with 4 cols (3 coordinates with 1 weight) and number of atoms + x rows
    //!         ( which have the coordinates of the last row with 0 weight) so that they are a multiple of 64
    linal::Matrix< double> DensitySimulateGaussianSphere::AtomsToMatrix( const util::SiPtrVector< const biol::Atom> &ATOMS, const double BLOB_C)
    {
      linal::Matrix< double> atom_matrix( Tools::RoundUp( 64, ATOMS.GetSize()), 4, 0.0);

      double *row_ptr( atom_matrix.Begin());
      for
      (
        util::SiPtrVector< const biol::Atom>::const_iterator atom_itr( ATOMS.Begin()), atom_itr_end( ATOMS.End());
        atom_itr != atom_itr_end;
        ++atom_itr, row_ptr += 4
      )
      {
        const biol::Atom &current_atom( **atom_itr);
        std::copy( current_atom.GetCoordinates().Begin(), current_atom.GetCoordinates().End(), row_ptr);
        row_ptr[ 3] = BLOB_C * current_atom.GetType()->GetElementType()->GetProperty( chemistry::ElementTypeData::e_Mass);
      }

      // fill all remaining rows with the coordinates of the last atoms
      const linal::Vector3D &last_atom_coordinates( ATOMS.LastElement()->GetCoordinates());

      for( double *row_ptr_end( atom_matrix.End()); row_ptr != row_ptr_end; row_ptr += 4)
      {
        std::copy( last_atom_coordinates.Begin(), last_atom_coordinates.End(), row_ptr);
      }

      return atom_matrix;
    }

    //! @brief compile programs for given precision
    //! @param PRECISION float or double
    //! @return ERROR error that occured, CL_SUCCESS if no error
    cl_int DensitySimulateGaussianSphere::CompilePrograms( const util::CPPDataTypes::Types &PRECISION)
    {
      cl_int error_number( CL_SUCCESS);

      // compile the program
      cl::Program::Sources source;
      KernelSourceFile simulate_source_file( "simulate_density_gaussian_sphere.cl");

      const Device device( m_CommandQueue.GetDevice( &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "cannot get device for command queue");

      const std::string simulate_source( simulate_source_file.GetSource( PRECISION, device.Extensions()));
      if( simulate_source.empty())
      {
        return CL_INVALID_KERNEL_DEFINITION;
      }
      source.push_back( std::make_pair( simulate_source.c_str(), simulate_source.length()));

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

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_device.h"

// includes from bcl - sorted alphabetically
#include "opencl/bcl_opencl_tools.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

    //! @brief DataType as string
    //! @param DATA_TYPE the data type
    //! @return the DataType as string
    const std::string &Device::GetDataTypeString( const DataType &DATA_TYPE)
    {
      static const std::string s_data_type_strings[] =
      {
        "CHAR", "SHORT", "INT", "LONG", "FLOAT", "DOUBLE", "HALF",
        GetStaticClassName< DataType>()
      };

      return s_data_type_strings[ DATA_TYPE];
    }

    //! @brief convert DataType to cl device info preferred vector width
    //! @param DATA_TYPE the data type
    //! @return cl_device_info preferred vector width for data type
    cl_device_info Device::DataTypeToDeviceInfoPreferredVectorWidth( const DataType &DATA_TYPE)
    {
      switch( DATA_TYPE)
      {
        case e_Char:   return CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR;
        case e_Short:  return CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT;
        case e_Int:    return CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT;
        case e_Long:   return CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG;
        case e_Float:  return CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT;
        case e_Double: return CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE;
        case e_Half:   return CL_DEVICE_PREFERRED_VECTOR_WIDTH_HALF;
        default:       return 0;
      }
    }

    //! @brief convert DataType to cl device info native vector width
    //! @param DATA_TYPE the data type
    //! @return cl_device_info native vector width for data type
    cl_device_info Device::DataTypeToDeviceInfoNativeVectorWidth( const DataType &DATA_TYPE)
    {
      switch( DATA_TYPE)
      {
        case e_Char:   return CL_DEVICE_NATIVE_VECTOR_WIDTH_CHAR;
        case e_Short:  return CL_DEVICE_NATIVE_VECTOR_WIDTH_SHORT;
        case e_Int:    return CL_DEVICE_NATIVE_VECTOR_WIDTH_INT;
        case e_Long:   return CL_DEVICE_NATIVE_VECTOR_WIDTH_LONG;
        case e_Float:  return CL_DEVICE_NATIVE_VECTOR_WIDTH_FLOAT;
        case e_Double: return CL_DEVICE_NATIVE_VECTOR_WIDTH_DOUBLE;
        case e_Half:   return CL_DEVICE_NATIVE_VECTOR_WIDTH_HALF;
        default:       return 0;
      }
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Device::Device() :
      cl::Device()
    {
    }

    //! @brief construct from cl::Device
    //! @param DEVICE the cl::Device
    Device::Device( const cl::Device &DEVICE) :
      cl::Device( DEVICE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new Device
    Device *Device::Clone() const
    {
      return new Device( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Device::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief access the name defined by CL_DEVICE_NAME
    //! @param ERROR_PTR error will be written to this location
    //! @return name of device
    std::string Device::Name( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_NAME>( ERROR_PTR);
    }

    //! @brief get the vendor defined by CL_DEVICE_VENDOR
    //! @param ERROR_PTR error will be written to this location
    //! @return the vendor of the device
    std::string Device::Vendor( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_VENDOR>( ERROR_PTR);
    }

    //! @brief get the version defined by CL_DEVICE_VERSION
    //! @param ERROR_PTR error will be written to this location
    //! @return the version of the device
    std::string Device::Version( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_VERSION>( ERROR_PTR);
    }

    //! @brief get the driver version defined by CL_DRIVER_VERSION
    //! @param ERROR_PTR error will be written to this location
    //! @return the vendor of the device
    std::string Device::DriverVersion( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DRIVER_VERSION>( ERROR_PTR);
    }

    //! @brief returns device type defined by CL_DEVICE_TYPE
    //! @param ERROR_PTR error will be written to this location
    //! @return device type
    cl_device_type Device::DeviceType( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_TYPE>( ERROR_PTR);
    }

    //! @brief vendor id as defined by CL_DEVICE_VENDOR_ID
    //! @param ERROR_PTR error will be written to this location
    //! @return vendor id
    cl_uint Device::VendorID( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_VENDOR_ID>( ERROR_PTR);
    }

    //! @brief the platform this device is associated with
    //! @param ERROR_PTR error will be written to this location
    //! @return the platform
    Platform Device::GetPlatform( cl_int *ERROR_PTR) const
    {
      return Platform( getInfo< CL_DEVICE_PLATFORM>( ERROR_PTR));
    }

    //! @brief return number of compute units on device
    //! @param ERROR_PTR error will be written to this location
    //! @return number of compute units
    cl_uint Device::MaxComputeUnits( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_VENDOR_ID>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS
    //! @param ERROR_PTR error will be written to this location
    //! @return max work item dimension
    cl_uint Device::MaxWorkItemDimension( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_WORK_ITEM_SIZES
    //! @param ERROR_PTR error will be written to this location
    //! @return the max work item size for all dimensions
    std::vector< size_t> Device::MaxWorkItemSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_WORK_ITEM_SIZES>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_WORK_GROUP_SIZE
    //! @param ERROR_PTR error will be written to this location
    //! @return the max group size - which means the max sum of work item sizes
    size_t Device::MaxWorkGroupSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_WORK_GROUP_SIZE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_CLOCK_FREQUENCY in MHz
    //! @param ERROR_PTR error will be written to this location
    //! @return max clock for the device
    cl_uint Device::MaxClockFrequency( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_CLOCK_FREQUENCY>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_ADDRESS_BITS
    //! @param ERROR_PTR error will be written to this location
    //! @return a bitfiled for the address bits
    cl_bitfield Device::AddressBits( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_ADDRESS_BITS>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_MEM_ALLOC_SIZE
    //! @param ERROR_PTR error will be written to this location
    //! @return maximal size for allocated memory
    cl_ulong Device::MaxMemAllocSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_MEM_ALLOC_SIZE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_GLOBAL_MEM_SIZE
    //! @param ERROR_PTR error will be written to this location
    //! @return total memory size for global memory
    cl_ulong Device::GlobalMemSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_GLOBAL_MEM_SIZE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_ERROR_CORRECTION_SUPPORT
    //! @param ERROR_PTR error will be written to this location
    //! @return is error correcte memory supported (ECC)
    cl_bool Device::ErrorCorrection( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_ERROR_CORRECTION_SUPPORT>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_LOCAL_MEM_TYPE
    //! @param ERROR_PTR error will be written to this location
    //! @return type of local memory
    cl_device_local_mem_type Device::LocalMemType( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_LOCAL_MEM_TYPE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_LOCAL_MEM_SIZE
    //! @param ERROR_PTR error will be written to this location
    //! @return size of local memory
    cl_ulong Device::LocalMemSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_LOCAL_MEM_SIZE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE
    //! @param ERROR_PTR error will be written to this location
    //! @return size of constant memory buffer
    cl_ulong Device::MaxConstantBufferSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MEM_BASE_ADDR_ALIGN
    //! @param ERROR_PTR error will be written to this location
    //! @return address alignment
    cl_uint Device::MemBaseAddrAlign( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MEM_BASE_ADDR_ALIGN>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE
    //! @param ERROR_PTR error will be written to this location
    //! @return
    cl_uint Device::MinDataTypeAlignSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MIN_DATA_TYPE_ALIGN_SIZE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_QUEUE_PROPERTIES & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE
    //! @param ERROR_PTR error will be written to this location
    //! @return does device support out of order execution with multiple command queues
    cl_bool Device::QueueOutOfOrderExecution( cl_int *ERROR_PTR) const
    {
      const cl_command_queue_properties device_command_queue_properties( getInfo< CL_DEVICE_QUEUE_PROPERTIES>( ERROR_PTR));
      return device_command_queue_properties & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE;
    }

    //! @brief defined by CL_DEVICE_QUEUE_PROPERTIES & CL_QUEUE_PROFILING_ENABLE
    //! @param ERROR_PTR error will be written to this location
    //! @return is queue profiling available
    cl_bool Device::QueueProfiling( cl_int *ERROR_PTR) const
    {
      const cl_command_queue_properties device_command_queue_properties( getInfo< CL_DEVICE_QUEUE_PROPERTIES>( ERROR_PTR));
      return device_command_queue_properties & CL_QUEUE_PROFILING_ENABLE;
    }

    //! @brief defined by CL_DEVICE_PROFILING_TIMER_RESOLUTION
    //! @param ERROR_PTR error will be written to this location
    //! @return timer resolution for queue profiling
    size_t Device::ProfilingTimerResolution( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_PROFILING_TIMER_RESOLUTION>();
    }

    //! @brief defined by CL_DEVICE_ENDIAN_LITTLE
    //! @param ERROR_PTR error will be written to this location
    //! @return is device little endian (false = big endian)
    cl_bool Device::EndianLittle( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_ENDIAN_LITTLE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_AVAILABLE
    //! @param ERROR_PTR error will be written to this location
    //! @return is device available
    cl_bool Device::Available( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_AVAILABLE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_COMPILER_AVAILABLE
    //! @param ERROR_PTR error will be written to this location
    //! @return is compiler for device available
    cl_bool Device::CompilerAvailable( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_AVAILABLE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_EXECUTION_CAPABILITIES & CL_EXEC_KERNEL
    //! @param ERROR_PTR error will be written to this location
    //! @return can kernel be executed
    cl_bool Device::ExecKernel( cl_int *ERROR_PTR) const
    {
      const cl_device_exec_capabilities device_exec_capabilities( getInfo< CL_DEVICE_EXECUTION_CAPABILITIES>( ERROR_PTR));
      return device_exec_capabilities & CL_EXEC_KERNEL;
    }

    //! @brief defined by CL_DEVICE_EXECUTION_CAPABILITIES & CL_EXEC_NATIVE_KERNEL
    //! @param ERROR_PTR error will be written to this location
    //! @return does device execute native kernels
    cl_bool Device::ExecNativeKernel( cl_int *ERROR_PTR) const
    {
      const cl_device_exec_capabilities device_exec_capabilities( getInfo< CL_DEVICE_EXECUTION_CAPABILITIES>( ERROR_PTR));
      return device_exec_capabilities & CL_EXEC_NATIVE_KERNEL;
    }

    //! @brief defined by CL_DEVICE_IMAGE_SUPPORT
    //! @param ERROR_PTR error will be written to this location
    //! @return does device support iamge data
    cl_bool Device::ImageSupport( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_IMAGE_SUPPORT>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_READ_IMAGE_ARGS
    //! @param ERROR_PTR error will be written to this location
    //! @return max number of image arguments
    cl_uint Device::MaxReadImageArgs( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_READ_IMAGE_ARGS>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_WRITE_IMAGE_ARGS
    //! @param ERROR_PTR error will be written to this location
    //! @return maximal number of image write arguments
    cl_uint Device::MaxWriteImageArgs( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_WRITE_IMAGE_ARGS>( ERROR_PTR);
    }

    //! @brief defined by CL_FP_DENORM
    //! @param ERROR_PTR error will be written to this location
    //! @param DATA_TYPE the denorm support for this data type (Float, Double, Half)
    //! @return floating point denorms supported
    cl_bool Device::FPDenorms( const DataType DATA_TYPE, cl_int *ERROR_PTR) const
    {
      const cl_int config_type( FpConfigFromDataType( DATA_TYPE));
      if( config_type == 0)
      {
        return CL_FALSE;
      }

      cl_device_fp_config device_fp_config;
      Tools::AssignError( ERROR_PTR, getInfo( config_type, &device_fp_config));

      return device_fp_config & CL_FP_DENORM;
    }

    //! @brief defined by CL_FP_INF_NAN
    //! @param ERROR_PTR error will be written to this location
    //! @param DATA_TYPE the inf nan support for this data type (Float, Double, Half)
    //! @return floating point inf and nan supported
    cl_bool Device::FPInfNan( const DataType DATA_TYPE, cl_int *ERROR_PTR) const
    {
      const cl_int config_type( FpConfigFromDataType( DATA_TYPE));
      if( config_type == 0)
      {
        return CL_FALSE;
      }

      cl_device_fp_config device_fp_config;
      Tools::AssignError( ERROR_PTR, getInfo( config_type, &device_fp_config));

      return device_fp_config & CL_FP_INF_NAN;
    }

    //! @brief defined by CL_FP_ROUND_TO_NEAREST
    //! @param ERROR_PTR error will be written to this location
    //! @param DATA_TYPE the round to nearest support for this data type (Float, Double, Half)
    //! @return floating point round to nearest supported
    cl_bool Device::FPRoundToNearest( const DataType DATA_TYPE, cl_int *ERROR_PTR) const
    {
      const cl_int config_type( FpConfigFromDataType( DATA_TYPE));
      if( config_type == 0)
      {
        return CL_FALSE;
      }

      cl_device_fp_config device_fp_config;
      Tools::AssignError( ERROR_PTR, getInfo( config_type, &device_fp_config));

      return device_fp_config & CL_FP_ROUND_TO_NEAREST;
    }

    //! @brief defined by CL_FP_ROUND_TO_ZERO
    //! @param ERROR_PTR error will be written to this location
    //! @param DATA_TYPE the round to zero support for this data type (Float, Double, Half)
    //! @return floating point round to sero supported
    cl_bool Device::FPRoundToZero( const DataType DATA_TYPE, cl_int *ERROR_PTR) const
    {
      const cl_int config_type( FpConfigFromDataType( DATA_TYPE));
      if( config_type == 0)
      {
        return CL_FALSE;
      }

      cl_device_fp_config device_fp_config;
      Tools::AssignError( ERROR_PTR, getInfo( config_type, &device_fp_config));

      return device_fp_config & CL_FP_ROUND_TO_ZERO;
    }

    //! @brief defined by CL_FP_ROUND_TO_INF
    //! @param ERROR_PTR error will be written to this location
    //! @param DATA_TYPE the round to inf support for this data type (Float, Double, Half)
    //! @return floating point round to +inf or -inf supported
    cl_bool Device::FPRoundToInf( const DataType DATA_TYPE, cl_int *ERROR_PTR) const
    {
      const cl_int config_type( FpConfigFromDataType( DATA_TYPE));
      if( config_type == 0)
      {
        return CL_FALSE;
      }

      cl_device_fp_config device_fp_config;
      Tools::AssignError( ERROR_PTR, getInfo( config_type, &device_fp_config));

      return device_fp_config & CL_FP_ROUND_TO_INF;
    }

    //! @brief defined by CL_FP_FMA
    //! @param ERROR_PTR error will be written to this location
    //! @param DATA_TYPE the fma support for this data type (Float, Double, Half)
    //! @return floating mointed fused multiply add supported
    cl_bool Device::FPFusedMultiplyAdd( const DataType DATA_TYPE, cl_int *ERROR_PTR) const
    {
      const cl_int config_type( FpConfigFromDataType( DATA_TYPE));
      if( config_type == 0)
      {
        return CL_FALSE;
      }

      cl_device_fp_config device_fp_config;
      Tools::AssignError( ERROR_PTR, getInfo( config_type, &device_fp_config));

      return device_fp_config & CL_FP_FMA;
    }

    //! @brief defined by CL_DEVICE_IMAGE2D_MAX_WIDTH
    //! @param ERROR_PTR error will be written to this location
    //! @return max width for 2d image
    size_t Device::Image2DMaxWidth( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_IMAGE2D_MAX_WIDTH>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_IMAGE2D_MAX_HEIGHT
    //! @param ERROR_PTR error will be written to this location
    //! @return max height for 2d image
    size_t Device::Image2DMaxHeight( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_IMAGE2D_MAX_HEIGHT>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_IMAGE3D_MAX_WIDTH
    //! @param ERROR_PTR error will be written to this location
    //! @return max width for 3d image
    size_t Device::Image3DMaxWidth( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_IMAGE3D_MAX_WIDTH>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_IMAGE3D_MAX_HEIGHT
    //! @param ERROR_PTR error will be written to this location
    //! @return max height for 3d image
    size_t Device::Image3DMaxHeight( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_IMAGE3D_MAX_HEIGHT>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_IMAGE3D_MAX_DEPTH
    //! @param ERROR_PTR error will be written to this location
    //! @return max depth for 3d image
    size_t Device::Image3DMaxDepth( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_IMAGE3D_MAX_DEPTH>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_SAMPLERS
    //! @param ERROR_PTR error will be written to this location
    //! @return max samlers for image
    cl_uint Device::MaxSamplers( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_SAMPLERS>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_MAX_PARAMETER_SIZE
    //! @param ERROR_PTR error will be written to this location
    //! @return max kernel parameter size
    size_t Device::MaxParameterSize( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_MAX_PARAMETER_SIZE>( ERROR_PTR);
    }

    //! @brief defined by CL_DEVICE_EXTENSIONS
    //! @param ERROR_PTR error will be written to this location
    //! @return all extensions the device supports
    storage::Set< Extension> Device::Extensions( cl_int *ERROR_PTR) const
    {
      return GetExtensions().ExtensionsFromString( getInfo< CL_DEVICE_EXTENSIONS>( ERROR_PTR));
    }

    //! @brief defined if extensions contain "cl_nv_device_attribute_query"
    //! @param ERROR_PTR error will be written to this location
    //! @return is nvidia device
    cl_bool Device::NVDevice( cl_int *ERROR_PTR) const
    {
      const storage::Set< Extension> extensions( Extensions( ERROR_PTR));
      return extensions.Find( GetExtensions().e_nv_device_attribute_query) != extensions.End();
    }

    //! @brief defined by CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV
    //! @param ERROR_PTR error will be written to this location
    //! @return major compute capability for nvidia device
    cl_uint Device::ComputeCapabilityMajorNV( cl_int *ERROR_PTR) const
    {
      cl_uint cc_major( 999);
      Tools::AssignError( ERROR_PTR, getInfo( CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV, &cc_major));
      return cc_major;
    }

    //! @brief defined by CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV
    //! @param ERROR_PTR error will be written to this location
    //! @return minor compute capability for nvidia device
    cl_uint Device::ComputeCapabilityMinorNV( cl_int *ERROR_PTR) const
    {
      cl_uint cc_minor( 999);
      Tools::AssignError( ERROR_PTR, getInfo( CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV, &cc_minor));
      return cc_minor;
    }

    //! @brief defined by CL_DEVICE_REGISTERS_PER_BLOCK_NV
    //! @param ERROR_PTR error will be written to this location
    //! @return registers per block for nvidia device
    cl_uint Device::RegistersPerBlockNV( cl_int *ERROR_PTR) const
    {
      cl_uint registers( 0);
      Tools::AssignError( ERROR_PTR, getInfo( CL_DEVICE_REGISTERS_PER_BLOCK_NV, &registers));
      return registers;
    }

    //! @brief defined by CL_DEVICE_WARP_SIZE_NV
    //! @param ERROR_PTR error will be written to this location
    //! @return warp size for nvidia device
    cl_uint Device::WarpSizeNV( cl_int *ERROR_PTR) const
    {
      cl_uint warp_size( 0);
      Tools::AssignError( ERROR_PTR, getInfo( CL_DEVICE_WARP_SIZE_NV, &warp_size));
      return warp_size;
    }

    //! @brief defined by CL_DEVICE_GPU_OVERLAP_NV
    //! @param ERROR_PTR error will be written to this location
    //! @return overlapping gpu for nvidia device
    cl_bool Device::GPUOverlapNV( cl_int *ERROR_PTR) const
    {
      cl_bool gpu_overlap( false);
      Tools::AssignError( ERROR_PTR, getInfo( CL_DEVICE_GPU_OVERLAP_NV, &gpu_overlap));
      return gpu_overlap;
    }

    //! @brief defined by CL_DEVICE_KERNEL_EXEC_TIMEOUT_NV
    //! @param ERROR_PTR error will be written to this location
    //! @return kernel execution timeout for nvidia device
    cl_bool Device::KernelExecTimeoutNV( cl_int *ERROR_PTR) const
    {
      cl_bool exec_timeout( false);
      Tools::AssignError( ERROR_PTR, getInfo( CL_DEVICE_KERNEL_EXEC_TIMEOUT_NV, &exec_timeout));
      return exec_timeout;
    }

    //! @brief defined by CL_DEVICE_INTEGRATED_MEMORY_NV
    //! @param ERROR_PTR error will be written to this location
    //! @return integrated memory for nvidia device
    cl_bool Device::IntegratedMemoryNV( cl_int *ERROR_PTR) const
    {
      cl_bool integrated_memory( false);
      Tools::AssignError( ERROR_PTR, getInfo( CL_DEVICE_INTEGRATED_MEMORY_NV, &integrated_memory));
      return integrated_memory;
    }

    //! @brief preferred and native vector width - native and DataType half only available from opencl 1.1
    //! @param DATA_TYPE for which datatype
    //! @param ERROR_PTR error will be written to this location
    //! @return pair of preferred and native vector width, undefined values for opencl < 1.1
    std::pair< cl_uint, cl_uint> Device::PreferredNativeVectorWidth( const DataType DATA_TYPE, cl_int *ERROR_PTR) const
    {
      std::pair< cl_uint, cl_uint> preferred_native( util::GetUndefined< cl_uint>(), util::GetUndefined< cl_uint>());

      Tools::AssignError( ERROR_PTR, getInfo( DataTypeToDeviceInfoPreferredVectorWidth( DATA_TYPE), &preferred_native.first));
      const cl_bool opencl11( OpenclCVersion( ERROR_PTR).find( "1.1") != std::string::npos);

      // TODO make this if statement test whether half type is available
      if( opencl11)
      {
        Tools::AssignError( ERROR_PTR, getInfo( DataTypeToDeviceInfoNativeVectorWidth( DATA_TYPE), &preferred_native.second));
      }

      return preferred_native;
    }

  ////////////////
  // opencl 1.1 //
  ////////////////

    //! defined by CL_DEVICE_OPENCL_C_VERSION
    //! @param ERROR_PTR error will be written to this location
    //! @return OPENCLCVERSION - empty for opencl < 1.1
    std::string Device::OpenclCVersion( cl_int *ERROR_PTR) const
    {
      std::string opencl_version;
      Tools::AssignError( ERROR_PTR, getInfo( CL_DEVICE_OPENCL_C_VERSION, &opencl_version));

      return opencl_version;
    }

    //! defined by CL_DEVICE_HOST_UNIFIED_MEMORY
    //! @param ERROR_PTR error will be written to this location
    //! @return
    cl_bool Device::HostUnifiedMemory( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_DEVICE_HOST_UNIFIED_MEMORY>( ERROR_PTR);
    }

    //! @brief get a description string for the device
    //! @param ERROR_PTR error will be written to this location
    //! @return string with a description for the device
    std::string Device::GetDescription( cl_int *ERROR_PTR) const
    {
      std::string description;

      description +=   "NAME                       " + Name( ERROR_PTR);
      description += "\nVENDOR                     " + Vendor( ERROR_PTR);
      description += "\nDRIVER_VERSION             " + DriverVersion( ERROR_PTR);

      description += "\nTYPE                       ";
      const cl_device_type type( DeviceType( ERROR_PTR));
      if( type & CL_DEVICE_TYPE_CPU)         description += "TYPE_CPU ";
      if( type & CL_DEVICE_TYPE_GPU)         description += "TYPE_GPU ";
      if( type & CL_DEVICE_TYPE_ACCELERATOR) description += "TYPE_ACCELERATOR ";
      if( type & CL_DEVICE_TYPE_DEFAULT)     description += "TYPE_DEFAULT ";

      description += "\nVERSION                    " + Version( ERROR_PTR);
      description += "\nVENDOR_ID                  " + util::Format()( VendorID( ERROR_PTR));
      description += "\nMAX_COMPUTE_UNITS          " + util::Format()( MaxComputeUnits( ERROR_PTR));
      description += "\nMAX_WORK_ITEM_DIMENSIONS   " + util::Format()( MaxWorkItemDimension( ERROR_PTR));

      description += "\nMAX_WORK_ITEM_SIZES        ";
      const std::vector< size_t> max_workitem_size( MaxWorkItemSize( ERROR_PTR));
      for( std::vector< size_t>::const_iterator itr( max_workitem_size.begin()), itr_end( max_workitem_size.end()); itr != itr_end; ++itr)
      {
        description += util::Format()( *itr) + " ";
      }

      description += "\nMAX_WORK_GROUP_SIZE        " + util::Format()( MaxWorkGroupSize( ERROR_PTR));
      description += "\nMAX_CLOCK_FREQUENCY        " + util::Format()( MaxClockFrequency( ERROR_PTR)) + " MHz";
      description += "\nADDRESS_BITS               " + util::Format()( AddressBits( ERROR_PTR));
      description += "\nMAX_MEM_ALLOC_SIZE         " + util::Format()( MaxMemAllocSize( ERROR_PTR) / ( 1024 * 1024)) + " MByte";
      description += "\nGLOBAL_MEM_SIZE            " + util::Format()( GlobalMemSize( ERROR_PTR) / ( 1024 * 1024)) + " MByte";
      description += "\nERROR_CORRECTION_SUPPORT   " + std::string(    ErrorCorrection( ERROR_PTR)   ? "yes" : "no");

      description += "\nLOCAL_MEM_TYPE             ";
      switch( LocalMemType( ERROR_PTR))
      {
        case CL_LOCAL:
          description += "local";
          break;
        case CL_GLOBAL:
          description += "global";
          break;
        default:
          description += "unknown";
          break;
      }

      description += "\nLOCAL_MEM_SIZE             " + util::Format()( LocalMemSize( ERROR_PTR) / 1024) + " kByte";
      description += "\nMAX_CONSTANT_BUFFER_SIZE   " + util::Format()( MaxConstantBufferSize( ERROR_PTR) / 1024) + " kByte";
      description += "\nMEM_BASE_ADDR_ALIGN        " + util::Format()( MemBaseAddrAlign( ERROR_PTR));
      description += "\nMIN_DATA_TYPE_ALIGN_SIZE   " + util::Format()( MinDataTypeAlignSize( ERROR_PTR));
      description += "\nQUEUE_OUT_OF_ORDER_EXEC    " + std::string(    QueueOutOfOrderExecution( ERROR_PTR) ? "yes" : "no");
      description += "\nQUEUE_PROFILING            " + std::string(    QueueProfiling( ERROR_PTR) ? "yes" : "no");
      description += "\nPROFILING_TIMER_RESOLUTION " + util::Format()( ProfilingTimerResolution( ERROR_PTR)) + " ns";

      description += "\nENDIAN_LITTLE              " + std::string( EndianLittle( ERROR_PTR)      ? "yes" : "no");
      description += "\nAVAILABLE                  " + std::string( Available( ERROR_PTR)         ? "yes" : "no");
      description += "\nCOMPILER_AVAILABLE         " + std::string( CompilerAvailable( ERROR_PTR) ? "yes" : "no");

      description += "\nEXEC_KERNEL                " + std::string( ExecKernel( ERROR_PTR)       ? "yes" : "no");
      description += "\nEXEC_NATIVE_KERNEL         " + std::string( ExecNativeKernel( ERROR_PTR) ? "yes" : "no");

      description += "\nIMAGE_SUPPORT              " + std::string(    QueueProfiling( ERROR_PTR) ? "yes" : "no");
      description += "\nMAX_READ_IMAGE_ARGS        " + util::Format()( MaxReadImageArgs( ERROR_PTR));
      description += "\nMAX_WRITE_IMAGE_ARGS       " + util::Format()( MaxWriteImageArgs( ERROR_PTR));

      for( size_t i( e_Float); i <= e_Half; ++i)
      {
        description += '\n' + GetDataTypeString( DataType( i));
        description += "\n  FP_DENORM                " + std::string( FPDenorms(          DataType( i), ERROR_PTR) ? "yes" : "no");
        description += "\n  FP_INF_NAN               " + std::string( FPInfNan(           DataType( i), ERROR_PTR) ? "yes" : "no");
        description += "\n  FP_ROUND_TO_NEAREST      " + std::string( FPRoundToNearest(   DataType( i), ERROR_PTR) ? "yes" : "no");
        description += "\n  FP_ROUND_TO_ZERO         " + std::string( FPRoundToZero(      DataType( i), ERROR_PTR) ? "yes" : "no");
        description += "\n  FP_ROUND_TO_INF          " + std::string( FPRoundToInf(       DataType( i), ERROR_PTR) ? "yes" : "no");
        description += "\n  FP_FMA                   " + std::string( FPFusedMultiplyAdd( DataType( i), ERROR_PTR) ? "yes" : "no");
      }

      description += "\nIMAGE2D_MAX_WIDTH          " + util::Format()( Image2DMaxWidth( ERROR_PTR) );
      description += "\nIMAGE2D_MAX_HEIGHT         " + util::Format()( Image2DMaxHeight( ERROR_PTR));
      description += "\nIMAGE3D_MAX_WIDTH          " + util::Format()( Image3DMaxWidth( ERROR_PTR) );
      description += "\nIMAGE3D_MAX_HEIGHT         " + util::Format()( Image3DMaxHeight( ERROR_PTR));
      description += "\nIMAGE3D_MAX_DEPTH          " + util::Format()( Image3DMaxDepth( ERROR_PTR) );

      description += "\nMAX_SAMPLERS               " + util::Format()( MaxSamplers( ERROR_PTR));

      description += "\nMAX_PARAMETER_SIZE         " + util::Format()( MaxParameterSize( ERROR_PTR));

      description += "\nEXTENSIONS                 " + Extensions::ExtensionsToString( Extensions( ERROR_PTR));

      // for nvidia device
      if( NVDevice( ERROR_PTR))
      {
        description += "\nREGISTERS_PER_BLOCK_NV     " + util::Format()( RegistersPerBlockNV( ERROR_PTR));
        description += "\nWARP_SIZE_NV               " + util::Format()( WarpSizeNV( ERROR_PTR));
        description += "\nGPU_OVERLAP_NV             " + std::string(    GPUOverlapNV( ERROR_PTR)        ? "yes" : "no");
        description += "\nKERNEL_EXEC_TIMEOUT_NV     " + std::string(    KernelExecTimeoutNV( ERROR_PTR) ? "yes" : "no");
        description += "\nINTEGRATED_MEMORY_NV       " + std::string(    IntegratedMemoryNV( ERROR_PTR)  ? "yes" : "no");
      }

      description += "\nPREFERRED_VECTOR_WIDTH_CHAR   " + util::Format()( PreferredNativeVectorWidth( e_Char  , ERROR_PTR).first);
      description += "\nPREFERRED_VECTOR_WIDTH_SHORT  " + util::Format()( PreferredNativeVectorWidth( e_Short , ERROR_PTR).first);
      description += "\nPREFERRED_VECTOR_WIDTH_INT    " + util::Format()( PreferredNativeVectorWidth( e_Int   , ERROR_PTR).first);
      description += "\nPREFERRED_VECTOR_WIDTH_LONG   " + util::Format()( PreferredNativeVectorWidth( e_Long  , ERROR_PTR).first);
      description += "\nPREFERRED_VECTOR_WIDTH_FLOAT  " + util::Format()( PreferredNativeVectorWidth( e_Float , ERROR_PTR).first);
      description += "\nPREFERRED_VECTOR_WIDTH_DOUBLE " + util::Format()( PreferredNativeVectorWidth( e_Double, ERROR_PTR).first);

      if( OpenclCVersion( ERROR_PTR).find( "1.1") == std::string::npos)
      {
        return description;
      }

      // opencl 1.1
      description += "\nPREFFERED_VECTOR_WIDTH_HALF   " + util::Format()( PreferredNativeVectorWidth( e_Half  , ERROR_PTR).first);
      description += "\nNATIVE_VECTOR_WIDTH_CHAR      " + util::Format()( PreferredNativeVectorWidth( e_Char  , ERROR_PTR).second);
      description += "\nNATIVE_VECTOR_WIDTH_SHORT     " + util::Format()( PreferredNativeVectorWidth( e_Short , ERROR_PTR).second);
      description += "\nNATIVE_VECTOR_WIDTH_INT       " + util::Format()( PreferredNativeVectorWidth( e_Int   , ERROR_PTR).second);
      description += "\nNATIVE_VECTOR_WIDTH_LONG      " + util::Format()( PreferredNativeVectorWidth( e_Long  , ERROR_PTR).second);
      description += "\nNATIVE_VECTOR_WIDTH_FLOAT     " + util::Format()( PreferredNativeVectorWidth( e_Float , ERROR_PTR).second);
      description += "\nNATIVE_VECTOR_WIDTH_DOUBLE    " + util::Format()( PreferredNativeVectorWidth( e_Double, ERROR_PTR).second);
      description += "\nNATIVE_VECTOR_WIDTH_HALF      " + util::Format()( PreferredNativeVectorWidth( e_Half  , ERROR_PTR).second);
      description += "\nOPENCL_C_VERSION              " + OpenclCVersion( ERROR_PTR);
      description += "\nHOST_UNIFIED_MEMORY           " + std::string(    HostUnifiedMemory( ERROR_PTR) ? "yes" : "no");

      // end
      return description;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Device::Read( std::istream &ISTREAM)
    {
      // write

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return output stream which was written to
    std::ostream &Device::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief the config flag to get the info for fp config
    //! @param DATA_TYPE float, Double or Half
    //! @return  CL_DEVICE_SINGLE_FP_CONFIG, CL_DEVICE_DOUBLE_FP_CONFIG, CL_DEVICE_HALF_FP_CONFIG, 0 for non available config
    cl_int Device::FpConfigFromDataType( const DataType &DATA_TYPE)
    {
      switch( DATA_TYPE)
      {
        case e_Float : return CL_DEVICE_SINGLE_FP_CONFIG;
        case e_Double: return CL_DEVICE_DOUBLE_FP_CONFIG;
        case e_Half  : return CL_DEVICE_HALF_FP_CONFIG;
        default      : return 0;
      }

      // end
      return 0;
    }

    //! @brief convert device type from string
    //! @param TYPE device type
    //! @return string for that device type
    const std::string &Device::TypeToString( const cl_device_type TYPE)
    {
      switch( TYPE)
      {
        case CL_DEVICE_TYPE_CPU:
        {
          static const std::string s_type_string_cpu( "TYPE_CPU");
          return s_type_string_cpu;
        }
        case CL_DEVICE_TYPE_GPU:
        {
          static const std::string s_type_string_gpu( "TYPE_GPU");
          return s_type_string_gpu;
        }
        case CL_DEVICE_TYPE_ACCELERATOR:
        {
          static const std::string s_type_string_accelerator( "TYPE_ACCELERATOR");
          return s_type_string_accelerator;
        }
        case CL_DEVICE_TYPE_ALL:
        {
          static const std::string s_type_string_all( "TYPE_ALL");
          return s_type_string_all;
        }
        case CL_DEVICE_TYPE_DEFAULT:
        default:
        {
          static const std::string s_type_string_default( "TYPE_DEFAULT");
          return s_type_string_default;
        }
      }
    }

    //! @brief convert string to device type
    //! @param TYPE_STRING string for device type
    //! @return string for that device type
    cl_device_type Device::TypeFromString( const std::string &TYPE_STRING)
    {
      if( TYPE_STRING == "TYPE_CPU")
      {
        return CL_DEVICE_TYPE_CPU;
      }
      else if( TYPE_STRING == "TYPE_GPU")
      {
        return CL_DEVICE_TYPE_GPU;
      }
      else if( TYPE_STRING == "TYPE_ACCELERATOR")
      {
        return CL_DEVICE_TYPE_ACCELERATOR;
      }
      else if( TYPE_STRING == "TYPE_ALL")
      {
        return CL_DEVICE_TYPE_ALL;
      }

      // end
      return CL_DEVICE_TYPE_DEFAULT;
    }

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_extension_data.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_instances.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////
  // data //
  //////////

    //! @brief prefix to all extensions "cl"
    //! @return prefix string for opencl extensions
    const std::string &ExtensionData::GetPrefix()
    {
      static const std::string s_prefix( "cl");
      return s_prefix;
    }

    //! @brief separator for prefix, vendor and name
    //! @return separator '_'
    char ExtensionData::GetSeparator()
    {
      return '_';
    }

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ExtensionData::s_Instance;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ExtensionData::ExtensionData()
    {
    }

    //! @brief construct from extension string
    //! @param EXTENSION_STRING string for extension with prefix, vendor and name, e.g. "cl_amd_fp64"
    ExtensionData::ExtensionData( const std::string &EXTENSION_STRING)
    {
      BCL_Assert( SetMembersFromExtensionString( EXTENSION_STRING), "incorrect extension string: " + EXTENSION_STRING);
    }

    //! @brief Clone function
    //! @return pointer to new Extension
    ExtensionData *ExtensionData::Clone() const
    {
      return new ExtensionData( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ExtensionData::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the vendor, eg. "amd"
    //! @return name of vendor
    const std::string &ExtensionData::GetVendor() const
    {
      return m_Vendor;
    }

    //! @brief return the name of extension, eg. "nv"
    //! @return name of extension
    const std::string &ExtensionData::GetName() const
    {
      return m_Name;
    }

    //! @brief get the extension string made up of prefix, vendor and name
    //! @return full extension string, e.g. "cl_amd_fp64"
    std::string ExtensionData::GetString() const
    {
      return GetPrefix() + GetSeparator() + m_Vendor + GetSeparator() + m_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ExtensionData::Read( std::istream &ISTREAM)
    {
      // read extension string
      std::string extension_string;
      io::Serialize::Read( extension_string, ISTREAM);
      BCL_Assert( SetMembersFromExtensionString( extension_string), "incorrect extension string: " + extension_string);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ExtensionData::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write extension string
      io::Serialize::Write( GetString(), OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief set vendor name and extension name form extension string
    //! @param EXTENSION_STRING string for extension with prefix, vendor and name, e.g. "cl_amd_fp64"
    //! @return true, if successful
    bool ExtensionData::SetMembersFromExtensionString( const std::string &EXTENSION_STRING)
    {
      // check prefix
      static const std::string prefix( GetPrefix() + GetSeparator());
      if( EXTENSION_STRING.compare( 0, prefix.length(), prefix) != 0)
      {
        return false;
      }

      // position of second separator
      const std::string::size_type pos_2nd_separator( EXTENSION_STRING.find( GetSeparator(), prefix.length()));
      if( pos_2nd_separator == std::string::npos)
      {
        return false;
      }

      // extract vendor name between prefix and second separator
      const std::string::const_iterator itr_beg( EXTENSION_STRING.begin());
      m_Vendor = std::string( itr_beg + prefix.length(), itr_beg + pos_2nd_separator);

      // extract name
      m_Name = std::string( itr_beg + pos_2nd_separator + 1, EXTENSION_STRING.end());

      // end
      return true;
    }

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_extensions.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_enumerate.hpp"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Extensions::Extensions() :
      e_khr_3d_image_writes              ( AddExtension( "cl_khr_3d_image_writes"              )),
      e_khr_byte_addressable_store       ( AddExtension( "cl_khr_byte_addressable_store"       )),
      e_khr_d3d9_sharing                 ( AddExtension( "cl_khr_d3d9_sharing"                 )),
      e_khr_d3d10_sharing                ( AddExtension( "cl_khr_d3d10_sharing"                )),
      e_khr_d3d11_sharing                ( AddExtension( "cl_khr_d3d11_sharing"                )),
      e_khr_icd                          ( AddExtension( "cl_khr_icd"                          )),
      e_khr_fp16                         ( AddExtension( "cl_khr_fp16"                         )),
      e_khr_fp64                         ( AddExtension( "cl_khr_fp64"                         )),
      e_khr_gl_event                     ( AddExtension( "cl_khr_gl_event"                     )),
      e_khr_gl_sharing                   ( AddExtension( "cl_khr_gl_sharing"                   )),
      e_khr_global_int32_base_atomics    ( AddExtension( "cl_khr_global_int32_base_atomics"    )),
      e_khr_global_int32_extended_atomics( AddExtension( "cl_khr_global_int32_extended_atomics")),
      e_khr_int64_base_atomics           ( AddExtension( "cl_khr_int64_base_atomics"           )),
      e_khr_int64_extended_atomics       ( AddExtension( "cl_khr_int64_extended_atomics"       )),
      e_khr_local_int32_base_atomics     ( AddExtension( "cl_khr_local_int32_base_atomics"     )),
      e_khr_local_int32_extended_atomics ( AddExtension( "cl_khr_local_int32_extended_atomics" )),
      e_khr_select_fprounding_mode       ( AddExtension( "cl_khr_select_fprounding_mode"       )),
      e_ext_device_fission               ( AddExtension( "cl_ext_device_fission"               )),
      e_amd_device_attribute_query       ( AddExtension( "cl_amd_device_attribute_query"       )),
      e_amd_event_callback               ( AddExtension( "cl_amd_event_callback"               )),
      e_amd_fp64                         ( AddExtension( "cl_amd_fp64"                         )),
      e_amd_media_ops                    ( AddExtension( "cl_amd_media_ops"                    )),
      e_amd_popcnt                       ( AddExtension( "cl_amd_popcnt"                       )),
      e_amd_printf                       ( AddExtension( "cl_amd_printf"                       )),
      e_amd_vec3                         ( AddExtension( "cl_amd_vec3"                         )),
      e_apple_gl_sharing                 ( AddExtension( "cl_apple_gl_sharing"                 )),
      e_nv_compiler_options              ( AddExtension( "cl_nv_compiler_options"              )),
      e_nv_d3d9_sharing                  ( AddExtension( "cl_nv_d3d9_sharing"                  )),
      e_nv_d3d10_sharing                 ( AddExtension( "cl_nv_d3d10_sharing"                 )),
      e_nv_d3d11_sharing                 ( AddExtension( "cl_nv_d3d11_sharing"                 )),
      e_nv_device_attribute_query        ( AddExtension( "cl_nv_device_attribute_query"        )),
      e_nv_pragma_unroll                 ( AddExtension( "cl_nv_pragma_unroll"                 )),
      e_intel_printf                     ( AddExtension( "cl_intel_printf"                     )),
      e_intel_immediate_execution        ( AddExtension( "cl_intel_immediate_execution"        )),
      e_intel_exec_by_local_thread       ( AddExtension( "cl_intel_exec_by_local_thread"       ))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Extensions::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add extension enum
    //! @param EXTENSION_STRING string for extension, e.g. "cl_amd_fp64"
    //! @return the enum that was added, e_Undefined if there was an error
    Extension Extensions::AddExtension( const std::string &EXTENSION_STRING)
    {
      return AddEnum( EXTENSION_STRING, ExtensionData( EXTENSION_STRING));
    }

    //! @brief get set of extensions from a string, that was queried CL_PLATFORM_EXTENSIONS or CL_DEVICE_EXTENSIONS
    //! @param EXTENSIONS_STRING string containing multiple extensions, whitespace separated
    //! @return set of extensions
    storage::Set< Extension> Extensions::ExtensionsFromString( const std::string &EXTENSIONS_STRING) const
    {
      // all extension strings
      const storage::Vector< std::string> extension_strings( util::SplitString( EXTENSIONS_STRING));
      const storage::Vector< Extension> extension_vector( extension_strings.Begin(), extension_strings.End());
      // set of extensions
      storage::Set< Extension> extensions;
      extensions.InsertElements( extension_vector.Begin(), extension_vector.End());

      // there will be an undefined entry if one of the extension strings is unknown
      if( extensions.Find( e_Undefined) != extensions.End())
      {
        BCL_MessageDbg( "found at least one unknown extension: " + EXTENSIONS_STRING);
        extensions.Erase( e_Undefined);
      }

      // end
      return extensions;
    }

    //! @brief convert set of extensions into one string for output
    //! @param EXTENSIONS set of extensions
    //! @return string like it would result form query CL_PLATFORM_EXTENSIONS
    std::string Extensions::ExtensionsToString( const storage::Set< Extension> &EXTENSIONS)
    {
      std::string descriptor;

      // iterate over all extensions
      for
      (
        storage::Set< Extension>::const_iterator itr( EXTENSIONS.Begin()), itr_end( EXTENSIONS.End());
        itr != itr_end;
        ++itr
      )
      {
        descriptor += itr->GetName() + " ";
      }

      // end
      return descriptor;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief construct on access function for all Extensions
    //! @return reference to only instances of Extensions
    const Extensions &GetExtensions()
    {
      return Extensions::GetEnums();
    }

  } // namespace opencl

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< opencl::ExtensionData, opencl::Extensions>;

  } // namespace util
} // namespace bcl
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
#include "opencl/bcl_opencl_feature_similarity_measures.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
  ////////////////////////////
  // explicit instantiation //
  ////////////////////////////

    template class BCL_API FeatureSimilarityMeasures< double>;
    template class BCL_API FeatureSimilarityMeasures< float>;
  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_kappa_nearest_neighbor.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////
  // data //
  //////////

    //! the default input range
    const math::Range< float> KappaNearestNeighbor::s_DefaultInputRange( 0, 1);

    //! the default output range
    const math::Range< float> KappaNearestNeighbor::s_DefaultOutputRange( 0, 1);

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> KappaNearestNeighbor::s_Instance
    (
      GetObjectInstances().AddInstance( new KappaNearestNeighbor())
    );

    const char *KappaNearestNeighbor::s_CLCompilerOptions =
            "-cl-mad-enable -cl-fast-relaxed-math";

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    KappaNearestNeighbor::KappaNearestNeighbor() :
      m_ReferenceData(),
      m_Kappa( 1),
      m_Queue(),
      m_GpuDist(),
      m_GpuSort()
    {
    }

    //! @brief default constructor
    KappaNearestNeighbor::KappaNearestNeighbor( const CommandQueue &QUEUE) :
      m_ReferenceData(),
      m_Kappa( 1),
      m_Queue( QUEUE),
      m_GpuDist( QUEUE),
      m_GpuSort( QUEUE)
    {
    }

    //! @brief constructor from training data, query data, and rescale functions
    //! @param REFERENCE_DATA data to train the NeuralNetwork on
    //! @param KAPPA kappa value for number of nearest neighbors to consider
    //! @param RESCALE_INPUT function for rescaling the input
    //! @param RESCALE_OUTPUT function for rescaling the output
    KappaNearestNeighbor::KappaNearestNeighbor
    (
      util::ShPtr< descriptor::Dataset> &REFERENCE_DATA,
      const size_t KAPPA,
      const CommandQueue &QUEUE
    ) :
      m_ReferenceData( REFERENCE_DATA),
      m_Kappa( KAPPA),
      m_Queue( QUEUE),
      m_GpuDist( QUEUE),
      m_GpuSort( QUEUE)
    {
      cl_int error_number = CompilePrograms( util::CPPDataTypes::DataTypeFromTemplate< float>());
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
      SetReferenceDataOnDevice();
    }

    //! @brief copy constructor
    KappaNearestNeighbor *KappaNearestNeighbor::Clone() const
    {
      return new KappaNearestNeighbor( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &KappaNearestNeighbor::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &KappaNearestNeighbor::GetAlias() const
    {
      static const std::string s_Name( "OpenCLKappaNearestNeighbor");
      return s_Name;
    }

  //////////////
  // operator //
  //////////////

    //! @brief Set the scaling of a feature set according to the model
    //! @param FEATURES feature set of interest
    //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
    //!       when operator() is called
    void KappaNearestNeighbor::Rescale( model::FeatureDataSet< float> &FEATURE) const
    {
      if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_ReferenceData->GetFeaturesPtr()->GetScaling())
      {
        FEATURE.DeScale();
        FEATURE.Rescale( *m_ReferenceData->GetFeaturesPtr()->GetScaling());
      }
    }

    //! @brief predict result with model using a NOT rescaled feature vector
    //! @param FEATURE not rescaled feature vector
    //! @return predicted result vector using a model
    model::FeatureDataSet< float> KappaNearestNeighbor::PredictWithoutRescaling( const model::FeatureDataSetInterface< float> &FEATURE) const
    {
      const size_t block( 16);
      const size_t row_pad( ( block - ( FEATURE.GetNumberFeatures() % block)) % block);
      const size_t col_pad( ( block - ( FEATURE.GetFeatureSize() % block)) % block);
      Matrix< float> feature( FEATURE.GetMatrix(), m_Queue, row_pad, col_pad);
      return model::FeatureDataSet< float>( operator()( feature).GetHostMatrix(), *m_ReferenceData->GetResultsPtr()->GetScaling());
    }

    //! @brief predict result with model using a rescaled feature vector
    //! @param FEATURE normalized or rescaled feature vector
    //! @return predicted result vector using a model
    model::FeatureDataSet< float> KappaNearestNeighbor::operator()( const model::FeatureDataSetInterface< float> &FEATURE) const
    {
      // handle the case where rescaling is necessary
      if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != *m_ReferenceData->GetFeaturesPtr()->GetScaling())
      {
        model::FeatureDataSet< float> feature( FEATURE);
        feature.Rescale( *m_ReferenceData->GetFeaturesPtr()->GetScaling());
        return PredictWithoutRescaling( feature).DeScale();
      }

      // data is already rescaled
      return PredictWithoutRescaling( FEATURE).DeScale();
    }

    //! @brief predict result with model using a rescaled feature vector
    //! @param FEATURE normalized or rescaled Matrix
    //! @return predicted result vector using a model
    Matrix< float> KappaNearestNeighbor::operator()( const Matrix< float> &FEATURE) const
    {
      Matrix< float> distances( m_GpuDist( m_ReferenceFeatures, FEATURE));
      Matrix< float> sorted_distances( m_GpuSort( distances, m_Kappa));
      Matrix< int> sorted_result_indeces( m_GpuSort.GetIndexMatrix());

      Matrix< float> predicted_results( CalculateWeightedResults( sorted_distances, sorted_result_indeces));

      return predicted_results;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! read KappaNearestNeighbor from std::istream
    std::istream &KappaNearestNeighbor::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_ReferenceData, ISTREAM);
      io::Serialize::Read( m_Kappa,         ISTREAM);

      // end
      return ISTREAM;
    }

    //! write KappaNearestNeighbor into std::ostream
    std::ostream &KappaNearestNeighbor::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_ReferenceData, OSTREAM) << '\n';
      io::Serialize::Write( m_Kappa,         OSTREAM);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief sets reference data on the device
    void KappaNearestNeighbor::SetReferenceDataOnDevice()
    {
      const size_t block( 16);
      const size_t row_pad( ( block - ( m_ReferenceData->GetFeaturesPtr()->GetNumberFeatures() % block)) % block);
      const size_t col_pad( ( block - ( m_ReferenceData->GetFeaturesPtr()->GetFeatureSize() % block)) % block);
      m_ReferenceFeatures = Matrix< float>( m_ReferenceData->GetFeaturesPtr()->GetMatrix(), m_Queue, row_pad, col_pad);
      m_ReferenceResults = Matrix< float>( m_ReferenceData->GetResultsPtr()->GetMatrix(), m_Queue);
    }

    //! @brief calculates the weighted results based on distance
    //! @param SORTED_DISTANCES the kappa closest distances
    //! @param SORTED_INDECES the index corresponding to that distance for the result
    //! @return the output matrix of results
    Matrix< float> KappaNearestNeighbor::CalculateWeightedResults
    (
      const Matrix< float> &SORTED_DISTANCES,
      const Matrix< int> &SORTED_INDECES
    ) const
    {
      cl_int error_number = CL_SUCCESS;

      const cl_uint rows( SORTED_DISTANCES.GetNumberRows());
      const cl_uint cols( SORTED_DISTANCES.GetNumberCols());
      const cl_uint row_pad( SORTED_DISTANCES.GetRowPadding());
      const cl_uint col_pad( SORTED_DISTANCES.GetColPadding());

      // launch kernel
      cl::Kernel kernel( m_Program, "CalculateWeightedResults", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      Matrix< float> results( SORTED_DISTANCES.GetNumberCols(), size_t( 1), m_Queue);

      error_number  = kernel.setArg(  0, SORTED_DISTANCES.GetData());
      error_number |= kernel.setArg(  1, SORTED_INDECES.GetData());
      error_number |= kernel.setArg(  2, m_ReferenceResults.GetData());
      error_number |= kernel.setArg(  3, results.GetData());
      error_number |= kernel.setArg(  4, cl_uint( m_Kappa));
      error_number |= kernel.setArg(  5, cl_uint( rows));
      error_number |= kernel.setArg(  6, cl_uint( cols));
      error_number |= kernel.setArg(  7, cl_uint( col_pad));
      error_number |= kernel.setArg(  8, cl_uint( row_pad));
      BCL_Assert( error_number == CL_SUCCESS, "setting kernel args error: " + Tools::ErrorString( error_number));

      const cl_uint block_size( 128);
      cl::NDRange offset;
      cl::NDRange block_dims( block_size);
      cl::NDRange worksize( Tools::RoundUp( block_size, cols));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
      BCL_Assert( error_number == CL_SUCCESS, "error launching kernel: " + Tools::ErrorString( error_number));

      m_Queue.finish();

      return results;
    }

    //! @brief compile programs for given precision
    //! @param PRECISION float or double
    //! @return ERROR error that occured, CL_SUCCESS if no error
    cl_int KappaNearestNeighbor::CompilePrograms( const util::CPPDataTypes::Types &PRECISION)
    {
      cl_int error_number( CL_SUCCESS);

      // compile the program
      cl::Program::Sources source;

      // construct opencl device
      const Device device( m_Queue.GetDevice( &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "cannot get device for command queue");

      BCL_MessageDbg( "device description: " + util::Format()( device.GetDescription()));

      // construct kernel source strings
      const std::string knn_kernels_source( ( *GetKernelSources().e_KappaNearestNeighbor)->GetSource( PRECISION, device.Extensions()));

      // check if kernel strings are empty
      if( knn_kernels_source.empty())
      {
        return CL_INVALID_KERNEL_DEFINITION;
      }

      // pushback strings to program sources vector
      source.push_back( std::make_pair( knn_kernels_source.c_str(), knn_kernels_source.length()));

      // create the program
      cl::Program &current_program( m_Program);
      current_program = cl::Program( m_Queue.GetContext(), source, &error_number);
      if( error_number != CL_SUCCESS)
      {
        BCL_MessageCrt( "create program error: " + Tools::ErrorString( error_number));
        return error_number;
      }

      // build the program
      error_number = current_program.build( std::vector< cl::Device>( 1, device), s_CLCompilerOptions);
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

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer KappaNearestNeighbor::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "see http://en.wikipedia.org/wiki/K-nearest_neighbor_algorithm"
      );

      parameters.AddInitializer
      (
        "kappa",
        "number of nearest neighbors to use for computing average",
        io::Serialization::GetAgentWithMin( &m_Kappa, size_t( 1)),
        "3"
      );

      return parameters;
    }

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_kernel_source_alternative.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from two alternatives
    //! @param ALTERNATIVE_SOURCE_1
    //! @param ALTERNATIVE_SOURCE_2
    KernelSourceAlternative::KernelSourceAlternative
    (
      const KernelSourceInterface &ALTERNATIVE_SOURCE_1,
      const KernelSourceInterface &ALTERNATIVE_SOURCE_2
    ) :
      m_Alternatives()
    {
      m_Alternatives.PushBack( util::ShPtr< KernelSourceInterface>( ALTERNATIVE_SOURCE_1.Clone()));
      m_Alternatives.PushBack( util::ShPtr< KernelSourceInterface>( ALTERNATIVE_SOURCE_2.Clone()));
    }

    //! @brief Clone function
    //! @return pointer to new KernelSourceAlternative
    KernelSourceAlternative *KernelSourceAlternative::Clone() const
    {
      return new KernelSourceAlternative( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &KernelSourceAlternative::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief identifier for kernel source, so that the progam build can be cached based on it
    //! @return identifier like the filename or the function name
    const std::string &KernelSourceAlternative::GetIdentifier() const
    {
      static const std::string s_identifier;
      return m_Alternatives.IsEmpty() ? s_identifier : m_Alternatives.FirstElement()->GetIdentifier();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief get the source for compilation
    //! @param PRECISION precision of the kernel
    //! @param EXTENSIONS extensions of the devices, this kernel is compiled for
    //! @return source with precision set correctly
    std::string KernelSourceAlternative::GetSource( const util::CPPDataTypes::Types PRECISION, const storage::Set< Extension> &EXTENSIONS) const
    {
      // iterate over all alternatives
      for
      (
        util::ShPtrVector< KernelSourceInterface>::const_iterator
          itr( m_Alternatives.Begin()), itr_end( m_Alternatives.End());
        itr != itr_end;
        ++itr
      )
      {
        const std::string source( ( *itr)->GetSource( PRECISION, EXTENSIONS));
        if( !source.empty())
        {
          return source;
        }
      }

      // no alternative gave source
      return std::string();
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &KernelSourceAlternative::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Alternatives, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &KernelSourceAlternative::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Alternatives, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_kernel_source_file.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"
#include "opencl/bcl_opencl_extensions.h"
#include "opencl/bcl_opencl_kernel_sources.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from filename
    //! @param FILE_NAME filename of the kernel
    KernelSourceFile::KernelSourceFile( const std::string &FILE_NAME) :
      m_FileName( FILE_NAME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new KernelSourceFile
    KernelSourceFile *KernelSourceFile::Clone() const
    {
      return new KernelSourceFile( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &KernelSourceFile::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief identifier for kernel source, so that the progam build can be cached based on it
    //! @return identifier like the filename or the function name
    const std::string &KernelSourceFile::GetIdentifier() const
    {
      return m_FileName;
    }

    //! @brief get the source for compilation
    //! @param PRECISION precision of the kernel
    //! @param EXTENSIONS extensions of the devices, this kernel is compiled for
    //! @return source with precision set correctly
    std::string KernelSourceFile::GetSource( const util::CPPDataTypes::Types PRECISION, const storage::Set< Extension> &EXTENSIONS) const
    {
      // source
      std::string source;

      // open file and acquire content
      io::IFStream read;
      if( !io::File::TryOpenIFStream( read, KernelSources::AddKernelPath( m_FileName)))
      {
        BCL_MessageDbg( "cannot find kernel source file: " + KernelSources::AddKernelPath( m_FileName));
        return source;
      }

      source += "#define PRECISION " + util::CPPDataTypes::GetCPPString( PRECISION) + '\n';
      if( PRECISION == util::CPPDataTypes::e_Double)
      {
        if( EXTENSIONS.Find( GetExtensions().e_amd_fp64) != EXTENSIONS.End())
        {
          source += "#pragma OPENCL EXTENSION " + GetExtensions().e_amd_fp64.GetName() + " : enable\n";
        }
        else if( EXTENSIONS.Find( GetExtensions().e_khr_fp64) != EXTENSIONS.End())
        {
          source += "#pragma OPENCL EXTENSION " + GetExtensions().e_khr_fp64.GetName() + " : enable\n";
        }
        else
        {
          BCL_MessageCrt( "double precision is not supported by given extensions");
          return std::string();
        }
      }

      // filecontent
      std::stringstream file_content;
      file_content << read.rdbuf();
      io::File::CloseClearFStream( read);

      // add file content to source
      source += file_content.str() + "\n";

      // end
      return source;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &KernelSourceFile::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_FileName, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &KernelSourceFile::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_FileName, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_kernel_source_interface.h"

// includes from bcl - sorted alphabetically
#include "opencl/bcl_opencl_extensions.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

    //! @brief is the given Precision compatible with the given extension
    //! @param PRECISION the precision in question
    //! @param EXTENSIONS extensions of the devices
    //! @return true if the precision is supported
    bool KernelSourceInterface::PrecisionCompatibleWithExtensions( const util::CPPDataTypes::Types PRECISION, const storage::Set< Extension> &EXTENSIONS)
    {
      switch( PRECISION)
      {
        case util::CPPDataTypes::e_Float:
          return true;
        case util::CPPDataTypes::e_Double:
          return EXTENSIONS.Find( GetExtensions().e_khr_fp64) != EXTENSIONS.End() || EXTENSIONS.Find( GetExtensions().e_amd_fp64) != EXTENSIONS.End();
        default: break;
      }

      return false;
    }

    //! @brief get compiler flags necessary for precision
    //! @param PRECISION the precision of the kernel to be compiled
    //! @return string containing the compiler options that should be used
    std::string KernelSourceInterface::GetPrecisionCompilerOptions( const util::CPPDataTypes::Types PRECISION)
    {
      switch( PRECISION)
      {
        case util::CPPDataTypes::e_Float:
          // prevent warning
          // double-precision constant is represented as
          // single-precision constant because double is not enabled
          // PRECISION c_sub = 0.0;
          //                   ^
          return std::string( " -cl-single-precision-constant");
        default:
          return std::string();
      }

      return std::string();
    }

    //! @brief get additional compiler flags
    //! @return string containing the additional compiler options that are desired by the user
    std::string KernelSourceInterface::GetAdditionalCompilerOptions()
    {
      // kernels have to be compiled with -g to add debug information, so that it can be debugged
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        return std::string( " -g");
      }

      return std::string();
    }

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_kernel_sources.h"

// includes from bcl - sorted alphabetically
#include "bcl_version.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_file_in_search_path.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "opencl/bcl_opencl_kernel_source_file.h"
#include "opencl/bcl_opencl_tools.h"
#include "util/bcl_util_enumerate.hpp"

// external includes - sorted alphabetically

// path for opencl kernels
#if defined (__MINGW32__)
  #define BCL_KERNEL_PATH "opencl_kernels/"
#elif defined (__GNUC__)
  #define BCL_KERNEL_PATH "/dors/meilerlab/apps/bcl/opencl_kernels/rev_4941/"
#elif defined (_MSC_VER)
  #define BCL_KERNEL_PATH "../../../opencl_kernels/"
#endif

namespace bcl
{
  namespace opencl
  {

    const std::string &GetDefaultKernelsSourcePath()
    {
      static const std::string s_default_kernels_path
      (
        GetVersion().IsLicense() ?
          std::string( GetVersion().GetInstallPrefix() + "/opencl_kernels/")
          :
          std::string( BCL_KERNEL_PATH)
      );

      return s_default_kernels_path;
    }

    //! flag to change opencl kernels path
    util::ShPtr< command::FlagInterface> &KernelSources::GetKernelsSourcePathFlag()
    {
      static util::ShPtr< command::FlagInterface> s_kernels_path_flag
      (
        new command::FlagStatic
        (
          "opencl_kernels_path",
          "Path from which to retrieve opencl kernel (.cl) files",
          command::Parameter
          (
            "path",
            "relative or absolute path to a directory containing opencl kernels",
            command::ParameterCheckFileInSearchPath
            (
              "opencl_kernels",
              GetDefaultKernelsSourcePath(),
              io::Directory::e_Dir
            ),
            ""
          )
        )
      );

      return s_kernels_path_flag;
    }

    //! flag for path to to change opencl kernel binaries
    util::ShPtr< command::FlagInterface> &KernelSources::GetKernelsBinaryPathFlag()
    {
      static util::ShPtr< command::FlagInterface> s_kernels_path_flag
      (
        new command::FlagStatic
        (
          "opencl_kernels_binary_path",
          "if path is defined, all compiled binaries are stored and retrieved in the given directory after compilation",
          command::Parameter
          (
            "path",
            "relative or absolute path to read/write compiled binaries from/to",
            ""
          ),
          // This signal must be made after all opencl flags have been set
          // TODO: add SetSignal to flag interface so that this can be setup from where the flags are added, rather
          // than forcing it to be here
          &Tools::UpdateCurrentPlatformDevicesQueuesFromCommandLineFlag
        )
      );

      return s_kernels_path_flag;
    }

    KernelSources::ProgramMap &KernelSources::GetBuildPrograms()
    {
      static ProgramMap s_build_programs;
      return s_build_programs;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    KernelSources::KernelSources() :
      e_EuclideanDistance            ( AddEnum( "EuclideanDistance" , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "euclidean_distance.cl"  )))),
      e_SvmKernelFunctions           ( AddEnum( "SvmKernelFunctions", util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "svm_kernel_functions.cl")))),
      e_InsertionSort                ( AddEnum( "InsertionSort"     , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "insertion_sort.cl"      )))),
      e_RMSD                         ( AddEnum( "RMSD"              , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "rmsd.cl"                )))),
      e_SequentialMinimalOptimization( AddEnum( "SMO"               , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "svr_kernels.cl"         )))),
      e_KappaNearestNeighbor         ( AddEnum( "KNN"               , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "knn_kernels.cl"         )))),
      e_ArtificialNeuralNetwork      ( AddEnum( "ANN"               , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "ann_kernels.cl"         )))),
      e_Linal                        ( AddEnum( "Linal"             , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "linal.cl"               )))),
      e_Quality                      ( AddEnum( "Quality"           , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "quality.cl"             )))),
      e_Saxs                         ( AddEnum( "SAXS"              , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "saxs_kernels.cl"        )))),
      e_Buffer                       ( AddEnum( "Buffer"            , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "buffer.cl"              )))),
      e_ClusterCenters               ( AddEnum( "ClusterCenters"    , util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "cluster_centers.cl"     ))))
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &KernelSources::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add the Opencl kernel path to filename
    //! @param FILE_NAME filename of kernel
    //! @return kernel_path/FILE_NAME
    std::string KernelSources::AddKernelPath( const std::string &FILE_NAME)
    {
      const std::string resolved_filename
      (
        util::GetRuntimeEnvironment().ResolveFileName( GetKernelsSourcePathFlag()->GetFirstParameter()->GetValue() + FILE_NAME)
      );

      BCL_Assert
      (
        !resolved_filename.empty(),
        "unable to resolve filename: " + GetKernelsSourcePathFlag()->GetFirstParameter()->GetValue() + FILE_NAME
      );

      return resolved_filename;
    }

    //! @brief compile a kernel and return a binary
    //! @param KERNEL the kernel to compile
    //! @param PRECISION the precision to use
    //! @param QUEUE the command queue to program will run in
    //! @param OPTIONS compiler options
    //! @param ERROR_PTR location to store the ERROR at
    //! @return the cl::Program, either the stored already compiled version, or the newly compiled version
    cl::Program
    KernelSources::Compile
    (
      const KernelSource &KERNEL,
      const util::CPPDataTypes::Types &PRECISION,
      const CommandQueue &QUEUE,
      const std::string &OPTIONS,
      cl_int *ERROR_PTR
    )
    {
      cl_int error_number( CL_SUCCESS);

      // program
      cl::Program program;

      // device the program will run on
      const Device device( QUEUE.GetDevice( &error_number));
      if( error_number != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error_number);
        return program;
      }

      // platform
      const Platform platform( device.GetPlatform( &error_number));
      if( error_number != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error_number);
        return program;
      }

      // context
      const Context context( QUEUE.GetContext( &error_number));
      if( error_number != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error_number);
        return program;
      }

      // source identifier
      std::string identifier( ( *KERNEL)->GetIdentifier());
      // identifier refined with platform and device name
      if( !identifier.empty())
      {
        identifier = platform.Name() + device.Name() + OPTIONS + identifier;
        identifier += '_' + util::CPPDataTypes::GetCPPDatatypeName( PRECISION) + ".ptx";
        identifier = util::RemoveSpacesFromString( identifier);
      }

      // check if for given identifier, a program was already build
      ProgramMap::const_iterator itr_ident( GetBuildPrograms().find( identifier));
      if( itr_ident != GetBuildPrograms().end())
      {
        std::map< CommandQueue, cl::Program>::const_iterator queue_itr( itr_ident->second.find( QUEUE));
        if( queue_itr != itr_ident->second.end())
        {
          return queue_itr->second;
        }
      }

      io::DirectoryEntry bin_file;
      if( GetKernelsBinaryPathFlag()->GetFirstParameter()->GetWasSetInCommandLine() && !identifier.empty())
      {
        // file for that kernel, platform and device and options
        bin_file = io::DirectoryEntry( GetKernelsBinaryPathFlag()->GetFirstParameter()->GetValue() + identifier);

        // if such a file exists, compilation might not be necessary
        if( bin_file.DoesExist())
        {
          // read the file and make a program out of it
          cl::Program::Binaries binaries;
          io::IFStream read;
          io::File::MustOpenIFStream( read, bin_file.GetFullName(), std::ios::binary);
          const std::string file_content( ( std::istreambuf_iterator< char>( read)), std::istreambuf_iterator< char>());
          io::File::CloseClearFStream( read);
          binaries.push_back( std::make_pair( file_content.c_str(), file_content.size()));
          std::vector< cl_int> binary_status( 1, CL_SUCCESS); // number of devices - error for each device
          program = cl::Program( context, std::vector< cl::Device>( 1, device), binaries, &binary_status, &error_number);

          // was program creation successful
          if( error_number == CL_SUCCESS)
          {
            for( size_t i( 0); i < binary_status.size(); ++i)
            {
              if( binary_status[ i] != CL_SUCCESS)
              {
                error_number = binary_status[ i];
                BCL_MessageCrt( "unable to load binary for device: " + util::Format()( i) + " with error: " + Tools::ErrorString( binary_status[ i]));
                break;
              }
            }

            // no problem loading binaries
            if( error_number == CL_SUCCESS)
            {
              // compile
              error_number = program.build( std::vector< cl::Device>( 1, device), OPTIONS.c_str());
              if( error_number == CL_SUCCESS)
              {
                BCL_MessageVrb
                (
                  "use opencl binary file instead of sources: " + bin_file.GetFullName()
                );
                GetBuildPrograms()[ identifier][ QUEUE] = program;
                return program;
              }
            }
            // some error occurred
            error_number = CL_SUCCESS;
            program = cl::Program();
          }

          // file could not be used
          BCL_MessageCrt
          (
            "unable to use opencl binary file, recompiling for: " + bin_file.GetFullName()
          );
          identifier.clear();
        }
      }

      // actual source of the kernel
      const std::string ccc_source( ( *KERNEL)->GetSource( PRECISION, device.Extensions( NULL)));
      if( ccc_source.empty())
      {
        Tools::AssignError( ERROR_PTR, CL_INVALID_KERNEL_DEFINITION);
        return program;
      }

      // create the program
      cl::Program::Sources source;
      source.push_back( std::make_pair( ccc_source.c_str()      , ccc_source.length()));
      program = cl::Program( context, source, &error_number);
      if( error_number != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error_number);
        return program;
      }

      // compile
      error_number = program.build( std::vector< cl::Device>( 1, device), OPTIONS.c_str());

      // error in compilation
      if( error_number != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error_number); // assign the build error
        BCL_MessageCrt( "build program error: " + Tools::ErrorString( error_number));
        BCL_MessageCrt( "build program source:\n" + ccc_source);
        std::string build_info;
        error_number = program.getBuildInfo( device, CL_PROGRAM_BUILD_LOG, &build_info);
        if( error_number != CL_SUCCESS)
        {
          BCL_MessageCrt( "get build info error: " + Tools::ErrorString( error_number));
        }
        else
        {
          BCL_MessageCrt( "build log: " + build_info);
        }
        program = cl::Program();
        return program;
      }

      // successful compilation - store the binary if identifier is not empty
      if( GetKernelsBinaryPathFlag()->GetFirstParameter()->GetWasSetInCommandLine() && !identifier.empty() && !bin_file.DoesExist())
      {
        Tools::LogPtx( program, bin_file.GetFullName());
      }

      // cache the program build
      if( !identifier.empty())
      {
        GetBuildPrograms()[ identifier][ QUEUE] = program;
      }

      // end
      return program;
    }

    //! @brief construct on access function for all KernelSources
    //! @return reference to only instances of KernelSources
    KernelSources &GetKernelSources()
    {
      return KernelSources::GetEnums();
    }

  } // namespace opencl

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Enumerate< ShPtr< opencl::KernelSourceInterface>, opencl::KernelSources>;

  } // namespace util
} // namespace bcl
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
#include "opencl/bcl_opencl_kernel_source_string.h"

// includes from bcl - sorted alphabetically
#include "opencl/bcl_opencl_extensions.h"
#include "storage/bcl_storage_set.h"
#include "util/bcl_util_message.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    KernelSourceString::KernelSourceString() :
      m_Source()
    {
    }

    //! @brief construct from string
    //! @param SOURCE_STRING kernel source string
    KernelSourceString::KernelSourceString( const std::string &SOURCE_STRING) :
      m_Source( SOURCE_STRING)
    {
    }

    //! @brief Clone function
    //! @return pointer to new KernelSourceString
    KernelSourceString *KernelSourceString::Clone() const
    {
      return new KernelSourceString( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &KernelSourceString::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief identifier for kernel source, so that the progam build can be cached based on it
    //! @return identifier like the filename or the function name
    const std::string &KernelSourceString::GetIdentifier() const
    {
      static const std::string s_identifier;
      return s_identifier;
    }

    //! @brief get the source for compilation
    //! @param PRECISION precision of the kernel
    //! @param EXTENSIONS extensions of the devices, this kernel is compiled for
    //! @return source with precision set correctly
    std::string KernelSourceString::GetSource( const util::CPPDataTypes::Types PRECISION, const storage::Set< Extension> &EXTENSIONS) const
    {
      std::string source;
      source += "#define PRECISION " + util::CPPDataTypes::GetCPPString( PRECISION) + '\n';
      if( PRECISION == util::CPPDataTypes::e_Double)
      {
        if( EXTENSIONS.Find( GetExtensions().e_amd_fp64) != EXTENSIONS.End())
        {
          source += "#pragma OPENCL EXTENSION " + GetExtensions().e_amd_fp64.GetName() + " : enable\n";
        }
        else if( EXTENSIONS.Find( GetExtensions().e_khr_fp64) != EXTENSIONS.End())
        {
          source += "#pragma OPENCL EXTENSION " + GetExtensions().e_khr_fp64.GetName() + " : enable\n";
        }
        else
        {
          BCL_MessageCrt( "double precision is not supported by given extensions");
          return std::string();
        }
      }
      source += m_Source + "\n";
      return source;
    }

  ////////////////
  // operations //
  ////////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &KernelSourceString::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Source, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return output stream which was written to
    std::ostream &KernelSourceString::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Source, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_matrix3x3.hpp"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Matrix3x3< float>;
    template class BCL_API Matrix3x3< double>;

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_operations.hpp"

// includes from bcl - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    template class BCL_API Operations< float>;
    template class BCL_API OperationsEnumHandler< float>;
    template class BCL_API Operations< double>;
    template class BCL_API OperationsEnumHandler< double>;

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_platform.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command_state.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "opencl/bcl_opencl_kernel_sources.h"
#include "opencl/bcl_opencl_tools.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
  //////////
  // data //
  //////////

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PlatformFlagHelper
    //! @brief Adds Opencl flags to default flags, but only if -opencl Disable was not given on the command line
    //!        This is necessary to avoid querying the opencl platform if opencl was disabled, because even querying the
    //!        platform causes the cuda runtime to initialize, which on cluster jobs results in the processes' virtual
    //!        memory appearing to be equal to the total (phys + virtual + GPU) memory available on the given node due
    //!        to cuda's universal virtual addressing system. This is bad because it can cause pbs cluster jobs to crash
    //!        because they don't (or can't) request all virtual memory on the given node
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Oct 09, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    struct BCL_API PlatformFlagHelper :
      public signal::Slots
    {
      static const PlatformFlagHelper s_Instance;

      PlatformFlagHelper()
      {
        command::CommandState::GetGlobalCommandState().GetParseArgumentsSignal().Connect
        (
          this,
          &PlatformFlagHelper::AddOpenclFlagsToDefaultFlags
        );
      }

      //! @brief function to add opencl flags in a well defined order to the app default flags enum, but only if opencl
      //!        was not explicitly disabled
      void AddOpenclFlagsToDefaultFlags( const command::CommandState &STATE)
      {
        const storage::Vector< std::string> &opencl_args( STATE.GetArguments( "opencl"));
        if
        (
          ( opencl_args.GetSize() == size_t( 1) && opencl_args( 0) == "Disable")
          || ( opencl_args.IsEmpty() && !STATE.GetState().Has( "opencl"))
        )
        {
          Platform::GetIsOpenclDisabled() = true;
        }
        if( Platform::GetPlatformFlag().IsDefined())
        {
          // add opencl flags
          command::GetAppDefaultFlags().AddDefaultFlag( Platform::GetPlatformFlag(), command::e_Opencl);
          // add these flags only if opencl functionality is not prohibited due to incomplete opencl installation
          if( Platform::GetPlatformFlag()->GetParameterList().GetSize() > size_t( 1))
          {
            command::GetAppDefaultFlags().AddDefaultFlag( KernelSources::GetKernelsSourcePathFlag(), command::e_Opencl);
            command::GetAppDefaultFlags().AddDefaultFlag( KernelSources::GetKernelsBinaryPathFlag(), command::e_Opencl);

            // Ensure that the last flag has the signal to update the compiled opencl kernels and enums
            BCL_Assert
            (
              command::GetAppDefaultFlags().GetDefaultFlagsOfType( command::e_Opencl).LastElement()->GetSignal()
              == &Tools::UpdateCurrentPlatformDevicesQueuesFromCommandLineFlag,
              "The last opencl flag must have the update signal on it"
            );
          }
        }
      }
    };

    //! @brief add the opencl flags to the default flags
    const PlatformFlagHelper PlatformFlagHelper::s_Instance;

    //! @brief access to the commandline flag
    //! @return flag to select platform and processor type
    const util::ShPtr< command::FlagInterface> &Platform::GetPlatformFlag()
    {
      static util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "opencl",
          "choice of opencl platform and device type"
        )
      );
      static bool s_attempted_initialization( false); // keep track of whether we have attempted initialization

      if( !s_attempted_initialization)
      {
        // cast to FlagStatic for pushback
        util::ShPtr< command::FlagStatic> flag( s_flag);
        s_attempted_initialization = true;
        cl_int error( CL_SUCCESS);
        storage::Vector< std::string> names
        (
          GetIsOpenclDisabled()
          ? storage::Vector< std::string>()
          : QueryPlatformNamesStandardized( &error)
        );
        // warn if platform names could not be retrieved
        if( error != CL_SUCCESS)
        {
          BCL_MessageCrt( "unable to get opencl platform names: " + Tools::ErrorString( error));
        }
        names.PushBack( "Disable");
        if( GetIsOpenclDisabled() || error != CL_SUCCESS)
        {
          flag->PushBack
          (
            util::ShPtr< command::Parameter>
            (
              new command::Parameter
              (
                "platform",
                "opencl platform; Disabled because this machine lacks libOpenCL or does not have the appropriate "
                "/etc/OpenCL/vendors .icd files, or Disable was already given",
                command::ParameterCheckAllowed( names),
                names.LastElement() // Disable
              )
            )
          );
          return s_flag;
        }
        cl_device_type type;
        Platform optimal;
        optimal = GetOptimalPlatform( type, &error);
        if( error != CL_SUCCESS)
        {
          BCL_MessageCrt( "unable to get optimal platform: " + Tools::ErrorString( error));
          flag->PushBack
          (
            util::ShPtr< command::Parameter>
            (
              new command::Parameter
              (
                "platform",
                "opencl platform; Disabled because this machine lacks libOpenCL or does not have the appropriate "
                "/etc/OpenCL/vendors .icd files",
                command::ParameterCheckAllowed( names),
                names.LastElement() // Disable
              )
            )
          );
        }
        else
        {
          BCL_MessageVrb( "optimal platform: " + optimal.Name() + " " + Device::TypeToString( type));
        }

        // only insert parameters if there is at least one platform with either CPU or GPU
        if( error == CL_SUCCESS)
        {
          flag->PushBack
          (
            util::ShPtr< command::Parameter>
            (
              new command::Parameter
              (
                "platform",
                "opencl platform; select Disable on older machines with little or buggy opencl support",
                command::ParameterCheckAllowed( names),
                StandardizeName( optimal.Name())
              )
            )
          );
          flag->PushBack( GetDeviceTypeParam( &type));
          flag->PushBack( GetDeviceIDsParam());
        }
      }

      // end
      return s_flag;
    }

    //! @brief access to the command line param for processor type
    //! @param DEVICE_TYPE_PTR the device type to be default, if NULL, use the GPU or last setting
    //! @return param to select processor type
    util::ShPtr< command::ParameterInterface> &Platform::GetDeviceTypeParam( const cl_device_type *DEVICE_TYPE_PTR)
    {
      static util::ShPtr< command::ParameterInterface> s_param;

      // first time
      if( !s_param.IsDefined())
      {
        s_param = util::ShPtr< command::ParameterInterface>
        (
          new command::Parameter
          (
            "device_type",
            "choice of device type",
            command::ParameterCheckAllowed
            (
              storage::Vector< std::string>::Create
              (
                Device::TypeToString( CL_DEVICE_TYPE_GPU),
                Device::TypeToString( CL_DEVICE_TYPE_CPU),
                Device::TypeToString( CL_DEVICE_TYPE_ACCELERATOR),
                Device::TypeToString( CL_DEVICE_TYPE_ALL),
                Device::TypeToString( CL_DEVICE_TYPE_DEFAULT)
              )
            ),
            Device::TypeToString( CL_DEVICE_TYPE_GPU)
          )
        );
      }

      if( DEVICE_TYPE_PTR != NULL)
      {
        s_param->SetDefaultParameter( Device::TypeToString( *DEVICE_TYPE_PTR));
      }

      // end
      return s_param;
    }

    //! @brief access to the command line param for device vendor id's
    //! @return param to select devices by vendor id's
    const util::ShPtr< command::ParameterInterface> &Platform::GetDeviceIDsParam()
    {
      static util::ShPtr< command::ParameterInterface> s_param
      (
        new command::Parameter
        (
          "device_ids",
          "comma separated list of device ids to be used for selected platform, e.g \"0,1\"",
          ""
        )
      );

      // end
      return s_param;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Platform::Platform()
    {
    }

    //! @brief construct from name
    //! @param NAME platform name e.g. "ATI", "NVIDIA", "ATI Stream", "NVIDIA CUDA"
    Platform::Platform( const std::string &NAME, cl_int *ERROR_PTR)
    {
      const std::string standardized_name( StandardizeName( NAME));
      cl_int error( CL_SUCCESS);

      // get platforms
      const storage::Vector< Platform> &platforms( QueryPlatforms( &error));

      if( error != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error);
        return;
      }

      // iterate over platforms and find the one with the correct name
      for( storage::Vector< Platform>::const_iterator itr( platforms.Begin()), itr_end( platforms.End()); itr != itr_end; ++itr)
      {
        if
        (
          StandardizeName( itr->Name()) == standardized_name
        )
        {
          *this = *itr;
          return;
        }
      }

      Tools::AssignError( ERROR_PTR, CL_INVALID_VALUE);
    }

    //! @brief construct from an opencl platform
    Platform::Platform( const cl::Platform &PLATFORM, cl_int *ERROR_PTR) :
      cl::Platform( PLATFORM)
    {}

    //! @brief Clone function
    //! @return pointer to new Platform
    Platform *Platform::Clone() const
    {
      return new Platform( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Platform::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief access the name defined by CL_PLATFORM_NAME
    //! @param ERROR_PTR error will be written to this location
    //! @return the name of that platform
    std::string Platform::Name( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_PLATFORM_NAME>( ERROR_PTR);
    }

    //! @brief get the name defined by CL_PLATFORM_NAME as standardized name
    //! @param ERROR_PTR error will be written to this location
    //! @return bcl standardized platform name (without spaces)
    std::string Platform::StandardizedName( cl_int *ERROR_PTR) const
    {
      return StandardizeName( Name( ERROR_PTR));
    }

    //! @brief get the version defined by CL_PLATFORM_VERSION
    //! @param ERROR_PTR error will be written to this location
    //! @return the version of the platform
    std::string Platform::Version( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_PLATFORM_VERSION>( ERROR_PTR);
    }

    //! @brief get the vendor defined by CL_PLATFORM_VENDOR
    //! @param ERROR_PTR error will be written to this location
    //! @return the vendor of the platform
    std::string Platform::Vendor( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_PLATFORM_VENDOR>( ERROR_PTR);
    }

    //! @brief get the profile defined by CL_PLATFORM_PROFILE
    //! @param ERROR_PTR error will be written to this location
    //! @return the profile of the platform
    std::string Platform::Profile( cl_int *ERROR_PTR) const
    {
      return getInfo< CL_PLATFORM_PROFILE>( ERROR_PTR);
    }

    //! @brief get the extensions defined by CL_PLATFORM_EXTENSIONS
    //! @param ERROR_PTR error will be written to this location
    //! @return set of extensions
    storage::Set< Extension> Platform::Extensions( cl_int *ERROR_PTR) const
    {
      return GetExtensions().ExtensionsFromString( getInfo< CL_PLATFORM_EXTENSIONS>( ERROR_PTR));
    }

    //! @brief does it support ICD extension and what is the suffix
    //! @param ERROR_PTR error will be written to this location
    //! @return pair of bool, wheter is supports the icd extensions, and a string representing the icd suffix
    std::pair< bool, std::string> Platform::ICDExtension( cl_int *ERROR_PTR) const
    {
      std::pair< bool, std::string> icd_extension;
      storage::Set< Extension> extensions( Extensions( ERROR_PTR));

      icd_extension.first = extensions.Find( GetExtensions().e_khr_icd) != extensions.End();
      if( icd_extension.first)
      {
        Tools::AssignError( ERROR_PTR, getInfo( CL_PLATFORM_ICD_SUFFIX_KHR, &icd_extension.second));
      }

      // end
      return icd_extension;
    }

    //! @brief get a description string for the device
    //! @return string with a description for the device
    std::string Platform::Description( cl_int *ERROR_PTR) const
    {
      std::string description;

      description +=   "PROFILE        " + Profile( ERROR_PTR);
      description += "\nVERSION        " + Version( ERROR_PTR);
      description += "\nNAME           " + Name( ERROR_PTR);
      description += "\nVENDOR         " + Vendor( ERROR_PTR)    ;
      description += "\nEXTENSIONS     " + Extensions::ExtensionsToString( Extensions( ERROR_PTR));

      const std::pair< bool, std::string> icd_extension( ICDExtension( ERROR_PTR));
      if( icd_extension.first)
      {
        description += "\nICD_SUFFIX_KHR " + icd_extension.second;
      }

      // end
      return description;
    }

    //! @brief get devices for that platform
    //! @param DEVICE_TYPE the device type to be queried
    //! @return list of devices that match the queried types
    storage::Vector< Device> Platform::Devices( cl_device_type DEVICE_TYPE, cl_int *ERROR_PTR) const
    {
      cl_int error( CL_SUCCESS);
      storage::Vector< Device> final_devices;
      std::vector< cl::Device> devices;

      // get devices of queried type
      error = getDevices( DEVICE_TYPE, &devices);

      // check error
      if( error == CL_DEVICE_NOT_FOUND)
      {
        BCL_MessageVrb
        (
          "cannot find device for platform and device type: " +
          Name() + ' ' + Device::TypeToString( DEVICE_TYPE)
        );

        return final_devices;
      }
      if( error != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error);
      }

      // iterate over devices and insert
      for( std::vector< cl::Device>::const_iterator itr( devices.begin()), itr_end( devices.end()); itr != itr_end; ++itr)
      {
        final_devices.PushBack( Device( *itr));
      }

      // end
      return final_devices;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Platform::Read( std::istream &ISTREAM)
    {
      // read the name
      std::string tmp_name;
      io::Serialize::Read( tmp_name, ISTREAM);
      if( !tmp_name.empty()) // plaform name was given - check if there is one with that name
      {
        cl_int error( CL_SUCCESS);
        Platform tmp_platform( tmp_name, &error);
        BCL_Assert( error == CL_SUCCESS, "unable to construct platform with name: " + tmp_name);
        *this = tmp_platform;
      }
      else // empty platform
      {
        *this = Platform();
      }

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT - number of indentations
    //! @return output stream which was written to
    std::ostream &Platform::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write the platform name
      io::Serialize::Write( Name(), OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief query all platforms available
    //! @param ERROR_PTR error will be written to this location
    //! @return vector of all platforms
    const storage::Vector< Platform> &Platform::QueryPlatforms( cl_int *ERROR_PTR)
    {
      static cl_int s_error( CL_SUCCESS);

      // final list
      static storage::Vector< Platform> s_final_platform_list;

      // if we already have a result, return it
      if( !s_final_platform_list.IsEmpty() || s_error != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, s_error);
        return s_final_platform_list;
      }

      // get all platforms
      std::vector< cl::Platform> platform_list;
      s_error = cl::Platform::get( &platform_list);
      if( s_error != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, s_error);
        return s_final_platform_list;
      }

      // final list
      // iterate over platforms and insert
      for( std::vector< cl::Platform>::const_iterator itr( platform_list.begin()), itr_end( platform_list.end()); itr != itr_end; ++itr)
      {
        s_final_platform_list.PushBack( Platform( *itr, ERROR_PTR));
      }

      // end
      return s_final_platform_list;
    }

    //! @brief query all platform names
    //! @param ERROR_PTR error will be written to this location
    //! @return vector of all platform names
    const storage::Vector< std::string> &Platform::QueryPlatformNamesStandardized( cl_int *ERROR_PTR)
    {
      cl_int error( CL_SUCCESS);

      // get platforms
      const storage::Vector< Platform> &platforms( QueryPlatforms( &error));

      // name
      static storage::Vector< std::string> s_names;

      if( error != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error);
        return s_names;
      }
      else if( !s_names.IsEmpty() || platforms.IsEmpty())
      {
        // already initialized names
        return s_names;
      }

      // iterate over platforms
      for( storage::Vector< Platform>::const_iterator itr( platforms.Begin()), itr_end( platforms.End()); itr != itr_end; ++itr)
      {
        s_names.PushBack( StandardizeName( itr->Name()));
      }

      // end
      return s_names;
    }

    //! @brief get default device type as set in command line
    //! @brief return device type
    cl_device_type Platform::CommandLineDeviceType()
    {
      return Device::TypeFromString( GetDeviceTypeParam()->GetValue());
    }

    //! @brief get platform and devices from the command line
    //! @param PLATFORM platform, will be set to the command line platform
    //! @param DEVICES vector of devices to initialize
    //! @param ERROR_PTR error will be written to this location
    void Platform::InitializeFromCommandLine
    (
      Platform &PLATFORM,
      storage::Vector< Device> &DEVICES,
      cl_int *ERROR_PTR
    )
    {
      if( !GetPlatformFlag().IsDefined())
      {
        Tools::AssignError( ERROR_PTR, CL_PLATFORM_NOT_FOUND_KHR);
        PLATFORM = Platform();
        DEVICES.Reset();
        return;
      }

      const std::string standardized_name( StandardizeName( GetPlatformFlag()->GetFirstParameter()->GetValue()));
      cl_int error( CL_SUCCESS);

      // get platforms
      const storage::Vector< Platform> &platforms( QueryPlatforms( &error));

      if( error != CL_SUCCESS)
      {
        PLATFORM = Platform();
        DEVICES.Reset();
        Tools::AssignError( ERROR_PTR, error);
        return;
      }

      // get the command line device type
      cl_device_type device_type( CommandLineDeviceType());

      // keep track of whether a platform with the given name was found
      bool found_platform_with_name( false);

      // if a platform with the given name was found, assign error accordingly
      cl_int devices_err( CL_SUCCESS);

      // iterate over platforms and find the one with the correct name
      for
      (
        storage::Vector< Platform>::const_iterator itr( platforms.Begin()), itr_end( platforms.End());
        itr != itr_end;
        ++itr
      )
      {
        if( StandardizeName( itr->Name()) != standardized_name)
        {
          continue;
        }
        found_platform_with_name = true;
        cl_int local_devices_err( CL_SUCCESS);
        DEVICES = itr->Devices( device_type, &local_devices_err);

        // if no devices are available for that platform, issue a warning message that the default type has been updated
        if( DEVICES.IsEmpty() && !GetDeviceTypeParam()->GetWasSetInCommandLine())
        {
          DEVICES = itr->Devices( CL_DEVICE_TYPE_ALL, &local_devices_err);
          if( local_devices_err != CL_SUCCESS)
          {
            BCL_MessageStd
            (
              "Warning: Found platform with no devices (bad OpenCL driver setup?)"
            );
          }
          else if( !DEVICES.IsEmpty())
          {
            BCL_MessageStd
            (
              "Note: no devices with the default type = " + GetDeviceTypeParam()->GetDefaultValue()
              + " was found; instead taking the first device of the given platform"
            );
            GetDeviceTypeParam()->SetDefaultParameter( Device::TypeToString( CL_DEVICE_TYPE_ALL));
          }
        }
        if( local_devices_err != CL_SUCCESS)
        {
          devices_err = local_devices_err;
          continue;
        }
        if( !DEVICES.IsEmpty())
        {
          PLATFORM = *itr;
          break;
        }
      }

      // platform not found
      if( !found_platform_with_name)
      {
        Tools::AssignError( ERROR_PTR, CL_INVALID_VALUE);
        PLATFORM = Platform();
        DEVICES.Reset();
        return;
      }

      // no devices found for the platform
      if( devices_err != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, devices_err);
        PLATFORM = Platform();
        DEVICES.Reset();
        return;
      }

      // filter devices, if device_vendor_ids parameter was given
      if( GetDeviceIDsParam()->GetWasSetInCommandLine())
      {
        // vendor ids from commandline
        const storage::Vector< size_t> ids_vector( util::SplitStringToNumerical< size_t>( GetDeviceIDsParam()->GetValue(), ","));

        // create set to condense to unique ids
        const storage::Set< size_t> ids( ids_vector.Begin(), ids_vector.End());

        storage::Vector< Device> new_devices;

        // iterate over ids
        for( storage::Set< size_t>::const_iterator id_itr( ids.Begin()), id_itr_end( ids.End()); id_itr != id_itr_end; ++id_itr)
        {
          if( *id_itr >= DEVICES.GetSize())
          {
            BCL_MessageCrt
            (
              "there is no device with id: " + util::Format()( *id_itr) + " for platform: " +
              PLATFORM.Name() + " => ignoring this id"
            );
          }
          else
          {
            new_devices.PushBack( DEVICES( *id_itr));
          }
        }

        // end
        DEVICES = new_devices;
      }
    }

    //! @brief first platform with at least one device of given type
    //! @param DEVICE_TYPE the device type to be queried
    //! @param ERROR_PTR error will be written to this location, CL_INVALID_PLATFORM if no platform was found
    //! @return platform with gpu
    Platform Platform::FirstPlatformWithDeviceType( const cl_device_type DEVICE_TYPE, cl_int *ERROR_PTR)
    {
      cl_int error( CL_SUCCESS);

      // get all platforms
      const storage::Vector< Platform> &platforms( QueryPlatforms( &error));

      if( error != CL_SUCCESS)
      {
        Tools::AssignError( ERROR_PTR, error);
        return Platform();
      }

      const Platform *optimal_platform( NULL);
      // iterate over platforms and find first with given device type
      // iterate over platforms
      for( storage::Vector< Platform>::const_iterator itr( platforms.Begin()), itr_end( platforms.End()); itr != itr_end; ++itr)
      {
        storage::Vector< Device> devices( itr->Devices( DEVICE_TYPE));
        if( !devices.IsEmpty())
        {
          Tools::AssignError( ERROR_PTR, CL_SUCCESS);
          optimal_platform = &( *itr);

          // prefer Intel over AMD, since AMD writes all binaries to /tmp/OCL*.so, which causes the tmp folder to overflow
          if( optimal_platform->Name() == "Intel(R) OpenCL")
          {
            return *optimal_platform;
          }
        }
      }

      // no optimal platform found
      if( optimal_platform != NULL)
      {
        Tools::AssignError( ERROR_PTR, CL_SUCCESS);
        return *optimal_platform;
      }

      Tools::AssignError( ERROR_PTR, CL_INVALID_PLATFORM);
      // end - no platform found
      return Platform();
    }

    //! @brief identify optimal platform
    //! searches first for plaform with gpu, if there is non, the first one with cpu
    //! @param DEVICE_TYPE location where the device type of the optimal platform is stored
    //! @param ERROR_PTR error will be written to this location, CL_INVALID_PLATFORM if no platform was found
    //! @return optimal Platform
    Platform Platform::GetOptimalPlatform( cl_device_type &DEVICE_TYPE, cl_int *ERROR_PTR)
    {
      cl_int error( CL_SUCCESS);

      Platform optimal;

      // try gpu first
      optimal = FirstPlatformWithDeviceType( CL_DEVICE_TYPE_GPU, &error);
      if( error == CL_SUCCESS)
      {
        DEVICE_TYPE = CL_DEVICE_TYPE_GPU;
        return optimal;
      }

      // try cpu
      error = CL_SUCCESS;
      optimal = FirstPlatformWithDeviceType( CL_DEVICE_TYPE_CPU, &error);
      if( error == CL_SUCCESS)
      {
        DEVICE_TYPE = CL_DEVICE_TYPE_CPU;
        return optimal;
      }

      // no such platform
      Tools::AssignError( ERROR_PTR, CL_INVALID_PLATFORM);
      return optimal;
    }

    //! @brief standardize platform names, by replacing spaces
    //! @param PLATFORM_NAME name of platform that contains spaces like "ATI Stream"
    //! @return NAME with '_' instead of ' ' like "ATI_Stream"
    std::string Platform::StandardizeName( const std::string &PLATFORM_NAME)
    {
      std::string new_name( PLATFORM_NAME);
      // iterate over string
      for( std::string::iterator itr( new_name.begin()), itr_end( new_name.end()); itr != itr_end; ++itr)
      {
        if( *itr == ' ')
        {
          *itr = '_';
        }
      }

      // end
      return new_name;
    }

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_protein_agreement_ccc.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_side_chain_factory.h"
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "density/bcl_density_protein_agreements.h"
#include "opencl/bcl_opencl_kernel_source_file.h"
#include "opencl/bcl_opencl_kernel_sources.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

    //! @class DensityProteinAgreementCCCEnumHandler
    //! @brief handler class for adding the density simulate enum handler
    class BCL_API DensityProteinAgreementCCCEnumHandler :
      public signal::Slots
    {

    private:

    //////////
    // data //
    //////////

      //! the enum in the density::ProteinAgreement
      density::ProteinAgreement e_DensityProteinAgreementCCCGaussianSphere;

      //! the only instance of this class
      static const DensityProteinAgreementCCCEnumHandler s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DensityProteinAgreementCCCEnumHandler() :
        e_DensityProteinAgreementCCCGaussianSphere( density::GetProteinAgreements().AddEnum( "ProteinAgreementCCCOpencl", util::ShPtr< ProteinAgreementCCC>()))
      {
        GetTools().GetQueueUpdateSignal().Connect( this, &DensityProteinAgreementCCCEnumHandler::UpdateEnum);
      }

      //! @brief update the enum with the command queue from the Tools
      //! @param TOOLS the tolls to get the commandqueue from
      void UpdateEnum( Tools &TOOLS)
      {
        util::ShPtr< ProteinAgreementCCC> sp_agreement( new ProteinAgreementCCC());
        if( !TOOLS.HasCommandQueues())
        {
          *e_DensityProteinAgreementCCCGaussianSphere = util::ShPtr< ProteinAgreementCCC>();
          return;
        }

        if( sp_agreement->Initialize( TOOLS.GetFirstCommandQueue()))
        {
          // just update the existing one with the new one
          *e_DensityProteinAgreementCCCGaussianSphere = sp_agreement;
        }
        else
        {
          BCL_MessageVrb( "unable to initialize enum: ProteinAgreementCCCOpencl");
        }
      }

    }; // class DensityProteinAgreementCCCEnumHandler

    //! instance of DensityProteinAgreementCCCEnumHandler
    const DensityProteinAgreementCCCEnumHandler DensityProteinAgreementCCCEnumHandler::s_Instance = DensityProteinAgreementCCCEnumHandler();

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param ADD_SIDECHAIN_ATOMS protein model will get side chains before the density simulation
    ProteinAgreementCCC::ProteinAgreementCCC( const bool ADD_SIDECHAIN_ATOMS) :
      m_SimulatedContourLevelCutoff( 0),
      m_UseSideChains( ADD_SIDECHAIN_ATOMS)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new ProteinAgreementCCC copied from this one
    ProteinAgreementCCC *ProteinAgreementCCC::Clone() const
    {
      return new ProteinAgreementCCC( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinAgreementCCC::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to the density simulator
    //! @return the density simulator used
    const util::ShPtr< density::SimulateInterface> &ProteinAgreementCCC::GetSimulator() const
    {
      return m_Simulate;
    }

    //! @brief set the simulator used for density agreement, e.g. for simulating density from the protein model
    //! @param SP_SIMULATOR ShPtr to SimulatInterface
    void ProteinAgreementCCC::SetSimulator( const util::ShPtr< density::SimulateInterface> &SP_SIMULATOR)
    {
      m_Simulate = SP_SIMULATOR;
    }

    //! @brief access to the density used for agreement calculation
    //! @return SiPtr to the density
    const util::SiPtr< const density::Map> &ProteinAgreementCCC::GetDensity() const
    {
      return m_Map;
    }

    //! @brief set the density used for agreement calculation
    //! @param SP_DENSITY SiPtr to the density map
    void ProteinAgreementCCC::SetDensityMap( const util::SiPtr< const density::Map> &SP_DENSITY)
    {
      m_Map = SP_DENSITY;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief is this class compatible with given command queue
    //! @param COMMAND_QUEUE the command queue this object would operate on
    //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
    bool ProteinAgreementCCC::IsCompatible( const CommandQueue &COMMAND_QUEUE) const
    {
      cl_int error_number( CL_SUCCESS);
      const Device device( COMMAND_QUEUE.GetDevice( &error_number));

      // can get device
      if( error_number != CL_SUCCESS)
      {
        return false;
      }

      const storage::Set< Extension> extensions( device.Extensions( &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "unable to get extensions from device");

      return KernelSourceInterface::PrecisionCompatibleWithExtensions
             (
               util::CPPDataTypes::DataTypeFromTemplate< double>(),
               extensions
             );
    }

    //! @brief initialize this class
    //! @brief COMMAND_QUEUE queue to use
    //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
    bool ProteinAgreementCCC::Initialize( const CommandQueue &COMMAND_QUEUE)
    {
      // check if this is a compatible command queue
      if( !IsCompatible( COMMAND_QUEUE))
      {
        return false;
      }

      // update the command queue
      m_CommandQueue = COMMAND_QUEUE;

      // for precision type
      cl_int error_number( CL_SUCCESS);
      const std::string options;
      m_Program = KernelSources::Compile
          (
            GetCCCKernel(),
            util::CPPDataTypes::DataTypeFromTemplate< double>(),
            m_CommandQueue,
            options,
            &error_number
          );
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

    //! @brief operator() taking the a ProteinModel and returning the cross correlation coefficient
    //! @param PROTEIN_MODEL
    //! @return correlation between the member density map and a simulated density map for PROTEIN_MODEL
    double ProteinAgreementCCC::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // use protein as it is
      if( !m_UseSideChains)
      {
        return -CrossCorrelationCoefficient( *m_Map, m_Simulate->operator ()( PROTEIN_MODEL.GetAtoms()), m_SimulatedContourLevelCutoff);
      }

      // generate new protein model with side chains from PROTEINMODEL
      util::ShPtr< assemble::ProteinModel> protein_model_with_side_chains
      (
        biol::AASideChainFactory( false, true).ProteinModelWithSideChains( PROTEIN_MODEL)
      );

      const density::Map simulated_map( m_Simulate->operator ()( protein_model_with_side_chains->GetAtoms()));

      // return correlation between protein model with side chains and member density map
      return -CrossCorrelationCoefficient( *m_Map, simulated_map, m_SimulatedContourLevelCutoff);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &ProteinAgreementCCC::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Map               , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Simulate          , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SimulatedContourLevelCutoff, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UseSideChains     , OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &ProteinAgreementCCC::Read( std::istream &ISTREAM)
    {
      // write member
      io::Serialize::Read( m_Map               , ISTREAM);
      io::Serialize::Read( m_Simulate          , ISTREAM);
      io::Serialize::Read( m_SimulatedContourLevelCutoff, ISTREAM);
      io::Serialize::Read( m_UseSideChains     , ISTREAM);

      // end
      return ISTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    const KernelSource &ProteinAgreementCCC::GetCCCKernel()
    {
      static const KernelSource e_ccc_kernel( GetKernelSources().AddEnum( "DensityCorrelation", util::ShPtr< KernelSourceInterface>( new KernelSourceFile( "density_correlation.cl"))));
      return e_ccc_kernel;
    }

    //! @brief copy the density map to the device
    //! @param DENSITY_MAP the map to be copied
    //! @param NEW_DIMENSIONS new dimensions for map - used for padding
    //! @return Buffer the buffer containing the map
    Vector< double> ProteinAgreementCCC::MapToDevice
    (
      const density::Map &DENSITY_MAP,
      const storage::VectorND< 3, size_t> &NEW_DIMENSIONS
    ) const
    {
      const storage::VectorND< 3, size_t> map_dimension( DENSITY_MAP.GetDimensions());

      // new dimensions are the same
      if
      (
        map_dimension.First() == NEW_DIMENSIONS.First() &&
        map_dimension.Second() == NEW_DIMENSIONS.Second() &&
        map_dimension.Third() == NEW_DIMENSIONS.Third()
      )
      {
        return TensorToDevice( DENSITY_MAP.GetData());
      }

      // create padded tensor
      const math::Tensor< double> padded_tensor
      (
        DENSITY_MAP.GetData().CreatePaddedTensor
        (
          NEW_DIMENSIONS.Third()  - map_dimension.Third(),
          NEW_DIMENSIONS.Second() - map_dimension.Second(),
          NEW_DIMENSIONS.First()  - map_dimension.First()
        )
      );

      // copy padded tensor
      return TensorToDevice( padded_tensor);
    }

    //! @brief copy the tensor to the device
    //! @param TENSOR the tensor to be copied
    //! @return Buffer the buffer containing the tensor
    Vector< double> ProteinAgreementCCC::TensorToDevice( const math::Tensor< double> &TENSOR) const
    {
      Vector< double> device_tensor( linal::Vector< double>( TENSOR.GetSize(), TENSOR.Begin()), m_CommandQueue);

      return device_tensor;
    }

    //! @brief calculate the cross correlation between experimental and simulated density map
    //! this ccc measure only considers voxels, if the intensity in the simulated map is above the given contour level
    //! this prevents that adjacent protein densities in experimental maps have a negative contribution to the measure
    //! @param EXPERIMENTAL_DENSITY_MAP map from experiment
    //! @param SIMULATED_DENSITY_MAP map simulated from protein structure
    //! @param CONTOUR_LEVEL_SIMULATED the minimal intensity in the simulated map voxel, to be considered for the CCC
    //! @return cross correlation coefficient
    double ProteinAgreementCCC::CrossCorrelationCoefficient
    (
      const density::Map &EXPERIMENTAL_DENSITY_MAP,
      const density::Map &SIMULATED_DENSITY_MAP,
      const double CONTOUR_LEVEL_SIMULATED
    ) const
    {
      // create common sub tensor
      const storage::VectorND< 2, math::Tensor< double> > exp_sim_sub_tensor( EXPERIMENTAL_DENSITY_MAP.CommonSubTensor( SIMULATED_DENSITY_MAP));

      // calculate cross correlation
      return CrossCorrelationCoefficient( exp_sim_sub_tensor.First(), exp_sim_sub_tensor.Second(), CONTOUR_LEVEL_SIMULATED);
    }

    //! @brief calculate the cross correlation between experimental and simulated density map
    //! this ccc measure only considers voxels, if the intensity in the simulated map is above the given contour level
    //! this prevents that adjacent protein densities in experimental maps have a negative contribution to the measure
    //! @param EXPERIMENTAL_SUBDENSITY_MAP map from experiment
    //! @param SIMULATED_SUBDENSITY_MAP map simulated from protein structure
    //! @param CONTOUR_LEVEL_SIMULATED the minimal intensity in the simulated map voxel, to be considered for the CCC
    //! @return cross correlation coefficient
    double ProteinAgreementCCC::CrossCorrelationCoefficient
    (
      const math::Tensor< double> &EXPERIMENTAL_SUBDENSITY_MAP,
      const math::Tensor< double> &SIMULATED_SUBDENSITY_MAP,
      const double CONTOUR_LEVEL_SIMULATED
    ) const
    {
      Vector< double> device_exp_input( TensorToDevice( EXPERIMENTAL_SUBDENSITY_MAP));
      Vector< double> device_sim_input( TensorToDevice( SIMULATED_SUBDENSITY_MAP));

      return CrossCorrelationCoefficient( device_exp_input, device_sim_input, EXPERIMENTAL_SUBDENSITY_MAP.GetSize(), CONTOUR_LEVEL_SIMULATED);
    }

    //! @brief calculate the cross correlation between experimental and simulated density map
    //! @see CrossCorrelationCoefficient
    //! @param EXPERIMENTAL_BUFFER map from experiment
    //! @param SIMULATED_BUFFER map simulated from protein structure
    //! @param GRID_SIZE number of elements in buffer
    //! @param CONTOUR_LEVEL_SIMULATED the minimal intensity in the simulated map voxel, to be considered for the CCC
    //! @return cross correlation coefficient
    double ProteinAgreementCCC::CrossCorrelationCoefficient
    (
      const Vector< double> &EXPERIMENTAL_BUFFER,
      const Vector< double> &SIMULATED_BUFFER,
      const size_t GRID_SIZE,
      const double CONTOUR_LEVEL_SIMULATED
    ) const
    {
      BCL_Assert( EXPERIMENTAL_BUFFER.GetSize() == SIMULATED_BUFFER.GetSize(), "map sizes do not match!");

      cl_int error_number = CL_SUCCESS;

      const Context context( m_CommandQueue.GetContext( &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "cannot get context from command queue");

      const Device device( m_CommandQueue.GetDevice( &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "cannot get device for command queue");

      const size_t block_size( 64);
      const size_t num_groups
      (
        device.DeviceType() == CL_DEVICE_TYPE_GPU ?
            ( ( GRID_SIZE % block_size == 0 ? 0 : 1) + GRID_SIZE / block_size) :
            ( ( GRID_SIZE % device.MaxComputeUnits() == 0 ? 0 : 1) + GRID_SIZE / device.MaxComputeUnits())
      );

      Vector< double> device_exp_sum_output    ( num_groups, m_CommandQueue);
      Vector< double> device_sim_sum_output    ( num_groups, m_CommandQueue);
      Vector< double> device_exp_sim_sum_output( num_groups, m_CommandQueue);
      Vector< double> device_exp2_sum_output   ( num_groups, m_CommandQueue);
      Vector< double> device_sim2_sum_output   ( num_groups, m_CommandQueue);
      Vector< int> device_count_output         ( num_groups, m_CommandQueue);

      cl::NDRange local_worksize;
      const cl::NDRange offset;
      cl::NDRange global_worksize;
      cl::Kernel kernel;

      // Create the kernel
      if( device.DeviceType() == CL_DEVICE_TYPE_GPU)
      {
        kernel = cl::Kernel( m_Program, "DensityCorrelation", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
        local_worksize  = cl::NDRange( block_size); // all thread blocks have same dimensions
        global_worksize = cl::NDRange( Tools::RoundUp( block_size, GRID_SIZE));
        error_number  = kernel.setArg(  0, EXPERIMENTAL_BUFFER.GetData());
        error_number |= kernel.setArg(  1, SIMULATED_BUFFER.GetData());
        error_number |= kernel.setArg(  2, CONTOUR_LEVEL_SIMULATED);
        error_number |= kernel.setArg(  3, cl_uint( GRID_SIZE));
        error_number |= kernel.setArg(  4, device_exp_sum_output.GetData());
        error_number |= kernel.setArg(  5, device_sim_sum_output.GetData());
        error_number |= kernel.setArg(  6, device_exp_sim_sum_output.GetData());
        error_number |= kernel.setArg(  7, device_exp2_sum_output.GetData());
        error_number |= kernel.setArg(  8, device_sim2_sum_output.GetData());
        error_number |= kernel.setArg(  9, device_count_output.GetData());
        error_number |= kernel.setArg( 10, block_size * sizeof( double), 0);
        error_number |= kernel.setArg( 11, block_size * sizeof( double), 0);
        error_number |= kernel.setArg( 12, block_size * sizeof( double), 0);
        error_number |= kernel.setArg( 13, block_size * sizeof( double), 0);
        error_number |= kernel.setArg( 14, block_size * sizeof( double), 0);
        error_number |= kernel.setArg( 15, block_size * sizeof( int), 0);
        BCL_Assert( error_number == CL_SUCCESS, "setting gpu kernel args error: " + Tools::ErrorString( error_number));
      }
      else
      {
        kernel = cl::Kernel( m_Program, "DensityCorrelationCPU", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
        BCL_MessageStd( "using special cpu optimized opencl kernel");
        local_worksize  = cl::NDRange( 1); // all thread blocks have same dimensions
        global_worksize = cl::NDRange( num_groups);
        error_number  = kernel.setArg(  0, EXPERIMENTAL_BUFFER.GetData());
        error_number |= kernel.setArg(  1, SIMULATED_BUFFER.GetData());
        error_number |= kernel.setArg(  2, CONTOUR_LEVEL_SIMULATED);
        error_number |= kernel.setArg(  3, cl_uint( GRID_SIZE));
        error_number |= kernel.setArg(  4, cl_uint( num_groups));
        error_number |= kernel.setArg(  5, device_exp_sum_output.GetData());
        error_number |= kernel.setArg(  6, device_sim_sum_output.GetData());
        error_number |= kernel.setArg(  7, device_exp_sim_sum_output.GetData());
        error_number |= kernel.setArg(  8, device_exp2_sum_output.GetData());
        error_number |= kernel.setArg(  9, device_sim2_sum_output.GetData());
        error_number |= kernel.setArg( 10, device_count_output.GetData());
        BCL_Assert( error_number == CL_SUCCESS, "setting cpu kernel args error: " + Tools::ErrorString( error_number));
      }

      // launching kernel
      error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, offset, global_worksize, local_worksize);
      BCL_Assert( error_number == CL_SUCCESS, "kernel launch error: " + Tools::ErrorString( error_number));

      linal::Vector< double> exp_sum_output    ( device_exp_sum_output.GetHostVector());
      linal::Vector< double> sim_sum_output    ( device_sim_sum_output.GetHostVector());
      linal::Vector< double> exp_sim_sum_output( device_exp_sim_sum_output.GetHostVector());
      linal::Vector< double> exp2_sum_output   ( device_exp2_sum_output.GetHostVector());
      linal::Vector< double> sim2_sum_output   ( device_sim2_sum_output.GetHostVector());
      linal::Vector< int>    count_output      ( device_count_output.GetHostVector());

      const double count_voxel( count_output.Sum());
      const double sum_exp_sim( exp_sim_sum_output.Sum());
      const double sum_exp( exp_sum_output.Sum());
      const double sum_sim( sim_sum_output.Sum());
      const double sum_exp2( exp2_sum_output.Sum());
      const double sum_sim2( sim2_sum_output.Sum());

      double correlation( count_voxel * sum_exp_sim - sum_exp * sum_sim);
      correlation /= math::Sqrt( count_voxel * sum_exp2 - math::Sqr( sum_exp)) * math::Sqrt( count_voxel * sum_sim2 - math::Sqr( sum_sim));

      return correlation;
    }

    //! @brief calculate the cross correlation between experimental and simulated density map
    //! @see CrossCorrelationCoefficient
    //! @param EXPERIMENTAL_BUFFER map from experiment
    //! @param SIMULATED_BUFFER map simulated from protein structure
    //! @param GRID_SIZE number of elements in buffer
    //! @param CONTOUR_LEVEL_SIMULATED the minimal intensity in the simulated map voxel, to be considered for the CCC
    //! @return cross correlation coefficient
    double ProteinAgreementCCC::CrossCorrelationCoefficient
    (
      const Vector< double> &EXPERIMENTAL_BUFFER,
      const Vector< double> &SIMULATED_BUFFER,
      const linal::Vector< int> &EXP_START,
      const storage::VectorND< 3, size_t> &EXP_DIMENSIONS,
      const linal::Vector< int> &SIM_START,
      const storage::VectorND< 3, size_t> &SIM_DIMENSIONS,
      const linal::Vector< int> &EXTENSION,
      const double CONTOUR_LEVEL_SIMULATED
    ) const
    {
      cl_int error_number = CL_SUCCESS;

      const size_t block_size_x( 4);
      const size_t block_size_y( 4);
      const size_t block_size_z( 4);
      const size_t block_size( block_size_x * block_size_y * block_size_z);
      const size_t num_groups_x( ( ( EXTENSION( 0) % block_size_x == 0 ? 0 : 1) + EXTENSION( 0) / block_size_x));
      const size_t num_groups_y( ( ( EXTENSION( 1) % block_size_y == 0 ? 0 : 1) + EXTENSION( 1) / block_size_y));
      const size_t num_groups_z( ( ( EXTENSION( 2) % block_size_z == 0 ? 0 : 1) + EXTENSION( 2) / block_size_z));
      const size_t num_groups( num_groups_x * num_groups_y * num_groups_z);

      Vector< double> device_exp_sum_output    ( num_groups, m_CommandQueue);
      Vector< double> device_sim_sum_output    ( num_groups, m_CommandQueue);
      Vector< double> device_exp_sim_sum_output( num_groups, m_CommandQueue);
      Vector< double> device_exp2_sum_output   ( num_groups, m_CommandQueue);
      Vector< double> device_sim2_sum_output   ( num_groups, m_CommandQueue);
      Vector< int> device_count_output         ( num_groups, m_CommandQueue);

      cl::NDRange local_worksize;
      const cl::NDRange offset;
      cl::NDRange global_worksize;
      cl::Kernel kernel;

      // Create the kernel
      kernel = cl::Kernel( m_Program, "DensityCorrelationOverlap", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      local_worksize = cl::NDRange( block_size_x, block_size_y, block_size_y); // all thread blocks have same dimensions
      global_worksize = cl::NDRange( Tools::RoundUp( block_size_x, EXTENSION( 0)), Tools::RoundUp( block_size_y, EXTENSION( 1)), Tools::RoundUp( block_size_z, EXTENSION( 2)));
      error_number  = kernel.setArg(  0, EXPERIMENTAL_BUFFER.GetData());
      error_number |= kernel.setArg(  1, SIMULATED_BUFFER.GetData());
      error_number |= kernel.setArg(  2, CONTOUR_LEVEL_SIMULATED);
      error_number |= kernel.setArg(  3, cl_int( EXP_START( 0)));
      error_number |= kernel.setArg(  4, cl_int( EXP_START( 1)));
      error_number |= kernel.setArg(  5, cl_int( EXP_START( 2)));
      error_number |= kernel.setArg(  6, cl_uint( EXP_DIMENSIONS( 0)));
      error_number |= kernel.setArg(  7, cl_uint( EXP_DIMENSIONS( 1)));
      error_number |= kernel.setArg(  8, cl_uint( EXP_DIMENSIONS( 2)));
      error_number |= kernel.setArg(  9, cl_int( SIM_START( 0)));
      error_number |= kernel.setArg( 10, cl_int( SIM_START( 1)));
      error_number |= kernel.setArg( 11, cl_int( SIM_START( 2)));
      error_number |= kernel.setArg( 12, cl_uint( SIM_DIMENSIONS( 0)));
      error_number |= kernel.setArg( 13, cl_uint( SIM_DIMENSIONS( 1)));
      error_number |= kernel.setArg( 14, cl_uint( SIM_DIMENSIONS( 2)));
      error_number |= kernel.setArg( 15, cl_uint( EXTENSION( 0)));
      error_number |= kernel.setArg( 16, cl_uint( EXTENSION( 1)));
      error_number |= kernel.setArg( 17, cl_uint( EXTENSION( 2)));
      error_number |= kernel.setArg( 18, device_exp_sum_output.GetData());
      error_number |= kernel.setArg( 19, device_sim_sum_output.GetData());
      error_number |= kernel.setArg( 20, device_exp_sim_sum_output.GetData());
      error_number |= kernel.setArg( 21, device_exp2_sum_output.GetData());
      error_number |= kernel.setArg( 22, device_sim2_sum_output.GetData());
      error_number |= kernel.setArg( 23, device_count_output.GetData());
      error_number |= kernel.setArg( 24, block_size * sizeof( double), 0);
      error_number |= kernel.setArg( 25, block_size * sizeof( double), 0);
      error_number |= kernel.setArg( 26, block_size * sizeof( double), 0);
      error_number |= kernel.setArg( 27, block_size * sizeof( double), 0);
      error_number |= kernel.setArg( 28, block_size * sizeof( double), 0);
      error_number |= kernel.setArg( 29, block_size * sizeof( int), 0);
      BCL_Assert( error_number == CL_SUCCESS, "setting gpu kernel args error: " + Tools::ErrorString( error_number));

      // launching kernel
      error_number = m_CommandQueue.enqueueNDRangeKernel( kernel, offset, global_worksize, local_worksize);
      BCL_Assert( error_number == CL_SUCCESS, "kernel launch error: " + Tools::ErrorString( error_number));

      linal::Vector< double> exp_sum_output    ( device_exp_sum_output.GetHostVector());
      linal::Vector< double> sim_sum_output    ( device_sim_sum_output.GetHostVector());
      linal::Vector< double> exp_sim_sum_output( device_exp_sim_sum_output.GetHostVector());
      linal::Vector< double> exp2_sum_output   ( device_exp2_sum_output.GetHostVector());
      linal::Vector< double> sim2_sum_output   ( device_sim2_sum_output.GetHostVector());
      linal::Vector< int> count_output         ( device_count_output.GetHostVector());

      const double count_voxel( count_output.Sum());
      const double sum_exp_sim( exp_sim_sum_output.Sum());
      const double sum_exp( exp_sum_output.Sum());
      const double sum_sim( sim_sum_output.Sum());
      const double sum_exp2( exp2_sum_output.Sum());
      const double sum_sim2( sim2_sum_output.Sum());

      double correlation( count_voxel * sum_exp_sim - sum_exp * sum_sim);
      correlation /= math::Sqrt( count_voxel * sum_exp2 - math::Sqr( sum_exp)) * math::Sqrt( count_voxel * sum_sim2 - math::Sqr( sum_sim));

      return correlation;
    }

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_quality_gdt.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_range.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "opencl/bcl_opencl_matrix.h"
#include "opencl/bcl_opencl_operations.h"
#include "opencl/bcl_opencl_vector.h"
#include "quality/bcl_quality_average.h"
#include "quality/bcl_quality_measures.h"
#include "quality/bcl_quality_superimpose_measures.h"
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

    //! @brief handler class for adding the quality superimpose enum handler
    class BCL_API GDTEnumHandler :
      public signal::Slots
    {

    private:

    //////////
    // data //
    //////////

      quality::Measure e_GDT_HA; //!< enum for high accuracy measure
      quality::Measure e_GDT_TS; //!< enum for ts measure

      //! the enum in the quality::SuperimposeMeasure
      storage::Map< double, quality::SuperimposeMeasure> m_GDTSuperImposeMeasures;
      //! the enum in the quality::Measure
      storage::Map< double, quality::Measure> m_GDTMeasures;

      //! the only instance of this class
      static const GDTEnumHandler s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      GDTEnumHandler() :
        e_GDT_HA( quality::GetMeasures().AddEnum( "OpenclGDT_HA", util::ShPtr< quality::MeasureInterface>())),
        e_GDT_TS( quality::GetMeasures().AddEnum( "OpenclGDT_TS", util::ShPtr< quality::MeasureInterface>()))
      {
        for( storage::Set< double>::const_iterator itr( quality::Measures::GetDistanceCutoffsTS().Begin()), itr_end( quality::Measures::GetDistanceCutoffsTS().End()); itr != itr_end; ++itr)
        {
          const std::string cutoff( util::Format().FFP( 0)( *itr));
          m_GDTSuperImposeMeasures[ *itr] = quality::GetSuperimposeMeasures().AddEnum( "OpenclGDT_" + cutoff + "A", util::ShPtr< QualityGDT>());
          m_GDTMeasures[ *itr] = quality::GetMeasures().AddEnum( "OpenclGDT_" + cutoff + "A", util::ShPtr< QualityGDT>());
        }
        // register enum with opencl queue update signal
        GetTools().GetQueueUpdateSignal().Connect( this, &GDTEnumHandler::UpdateEnum);
      }

      //! @brief update the enum with the command queue from the Tools
      //! @param TOOLS the tolls to get the commandqueue from
      void UpdateEnum( Tools &TOOLS)
      {
        if( !TOOLS.HasCommandQueues())
        {
          // iterate through all enums
          for( storage::Map< double, quality::SuperimposeMeasure>::iterator itr( m_GDTSuperImposeMeasures.Begin()), itr_end( m_GDTSuperImposeMeasures.End()); itr != itr_end; ++itr)
          {
            *itr->second = util::ShPtr< QualityGDT>();
          }
          for( storage::Map< double, quality::Measure>::iterator itr( m_GDTMeasures.Begin()), itr_end( m_GDTMeasures.End()); itr != itr_end; ++itr)
          {
            *itr->second = util::ShPtr< QualityGDT>();
          }

          *e_GDT_HA = util::ShPtr< QualityGDT>();
          *e_GDT_TS = util::ShPtr< QualityGDT>();

          return;
        }

        // iterate through all enums
        for( storage::Map< double, quality::SuperimposeMeasure>::iterator itr( m_GDTSuperImposeMeasures.Begin()), itr_end( m_GDTSuperImposeMeasures.End()); itr != itr_end; ++itr)
        {
          util::ShPtr< QualityGDT> sp_quality( new QualityGDT( itr->first));
          if( sp_quality->Initialize( TOOLS.GetFirstCommandQueue()))
          {
            *itr->second = sp_quality;
            *m_GDTMeasures.Find( itr->first)->second = sp_quality;
          }
          else
          {
            BCL_MessageVrb( "unable to initialize enum: " + itr->second.GetName());
          }
        }

        // HA
        {
          util::ShPtrVector< quality::SuperimposeInterface> gdts;
          for( storage::Set< double>::const_iterator itr( quality::Measures::GetDistanceCutoffsHA().Begin()), itr_end( quality::Measures::GetDistanceCutoffsHA().End()); itr != itr_end; ++itr)
          {
            const storage::Map< double, quality::SuperimposeMeasure>::const_iterator map_itr( m_GDTSuperImposeMeasures.Find( *itr));
            if( map_itr != m_GDTSuperImposeMeasures.End())
            {
              if( map_itr->second.IsDefined())
              {
                gdts.PushBack( *map_itr->second);
              }
              else
              {
                break;
              }
            }
            else
            {
              util::ShPtr< QualityGDT> sp_quality( new QualityGDT( *itr));
              if( sp_quality->Initialize( TOOLS.GetFirstCommandQueue()))
              {
                gdts.PushBack( sp_quality);
              }
              else
              {
                BCL_MessageVrb( "unable to initialize gdt for cutoff: " + util::Format()( *itr));
              }
            }
          }
          if( gdts.GetSize() == quality::Measures::GetDistanceCutoffsHA().GetSize())
          {
            *e_GDT_HA = util::ShPtr< quality::MeasureInterface>( new quality::Average( gdts));
          }
          else
          {
            BCL_MessageVrb( "unable to initialize enum: " + e_GDT_HA.GetName());
          }
        }

        // TS
        {
          util::ShPtrVector< quality::SuperimposeInterface> gdts;
          for( storage::Map< double, quality::SuperimposeMeasure>::const_iterator itr( m_GDTSuperImposeMeasures.Begin()), itr_end( m_GDTSuperImposeMeasures.End()); itr != itr_end; ++itr)
          {
            if( itr->second.IsDefined())
            {
              gdts.PushBack( *itr->second);
            }
            else
            {
              break;
            }
          }
          if( gdts.GetSize() == m_GDTMeasures.GetSize())
          {
            *e_GDT_TS = util::ShPtr< quality::MeasureInterface>( new quality::Average( gdts));
          }
          else
          {
            BCL_MessageVrb( "unable to initialize enum: " + e_GDT_TS.GetName());
          }
        }
      }

    }; // class GDTEnumHandler

    //! instance of DensitySimulateEnumHandler
    const GDTEnumHandler GDTEnumHandler::s_Instance = GDTEnumHandler();

  //////////
  // data //
  //////////

    const size_t QualityGDT::s_NumberRowsTransformationMatrix = 4;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from a single distance cutoff and seed length
    //! @param QUEUE command queue
    //! @param DISTANCE_CUTOFF distance cutoff to be used
    //! @param SEED_LENGTH length of seed
    QualityGDT::QualityGDT
    (
      const double &DISTANCE_CUTOFF,
      const size_t SEED_LENGTH
    ) :
      m_Queue(),
      m_QualityLCS( DISTANCE_CUTOFF, SEED_LENGTH)
    {
    }

    //! @brief construct from a single distance cutoff and seed length
    //! @param QUEUE command queue
    //! @param DISTANCE_CUTOFF distance cutoff to be used
    //! @param SEED_LENGTH length of seed
    QualityGDT::QualityGDT
    (
      const CommandQueue &QUEUE,
      const double &DISTANCE_CUTOFF,
      const size_t SEED_LENGTH
    ) :
      m_Queue( QUEUE),
      m_QualityLCS( m_Queue, DISTANCE_CUTOFF, SEED_LENGTH)
    {
    }

    //! @brief Clone function
    //! @return pointer to new QualityGDT
    QualityGDT *QualityGDT::Clone() const
    {
      return new QualityGDT( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &QualityGDT::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &QualityGDT::GetComparisonFunction() const
    {
      return **math::Comparisons< double>::GetEnums().e_Greater;
    }

    //! @brief get seed length
    //! @return seed length
    size_t QualityGDT::GetSeedLength() const
    {
      return m_QualityLCS.GetSeedLength();
    }

    //! @brief get distance cutoff
    //! @return distance cutoff
    double QualityGDT::GetDistanceCutoff() const
    {
      return m_QualityLCS.GetCutoff();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief is this class compatible with given command queue
    //! @param COMMAND_QUEUE the command queue this object would operate on
    //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
    bool QualityGDT::IsCompatible( const CommandQueue &COMMAND_QUEUE) const
    {
      return m_QualityLCS.IsCompatible( COMMAND_QUEUE);
    }

    //! @brief initialize this class
    //! @brief COMMAND_QUEUE queue to use
    //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
    bool QualityGDT::Initialize( const CommandQueue &COMMAND_QUEUE)
    {
      if( m_QualityLCS.Initialize( COMMAND_QUEUE))
      {
        m_Queue = COMMAND_QUEUE;
        return true;
      }
      m_Queue = CommandQueue();
      return false;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculates GDT between COORDINATES and REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return GDT between COORDINATES and REFERENCE_COORDINATES
    double QualityGDT::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // calculate GDT
      return CalculateGDTAndSuperimposition( COORDINATES, REFERENCE_COORDINATES).First();
    }

    //! @brief calculates the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    math::TransformationMatrix3D QualityGDT::CalculateSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // calculate GDT and return the superimposition
      return CalculateGDTAndSuperimposition( COORDINATES, REFERENCE_COORDINATES).Second();
    }

    //! @brief calculates GDT and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return pair of GDT and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    storage::Pair< double, math::TransformationMatrix3D> QualityGDT::CalculateGDTAndSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      // if the coordinates are empty
      if( COORDINATES.IsEmpty() || REFERENCE_COORDINATES.IsEmpty())
      {
        return storage::Pair< double, math::TransformationMatrix3D>
        (
          util::GetUndefinedDouble(), math::TransformationMatrix3D()
        );
      }

      // create opencl matrices
      const Matrix< double> coordinates( m_QualityLCS.GetQualityRMSD().MatrixFromCoordinates( COORDINATES, m_QualityLCS.GetQualityRMSD().s_BlockSize));
      const Matrix< double> ref_coordinates( m_QualityLCS.GetQualityRMSD().MatrixFromCoordinates( REFERENCE_COORDINATES, m_QualityLCS.GetQualityRMSD().s_BlockSize));

      const storage::Pair< double, Matrix< double> > gdt_transformation
      (
        CalculateGDTAndSuperimposition( coordinates, ref_coordinates)
      );

      // return the pair of GDT value and the corresponding transformation matrix
      return storage::Pair< double, math::TransformationMatrix3D>
             (
               gdt_transformation.First(),
               math::TransformationMatrix3D( gdt_transformation.Second().GetHostMatrix( 0, s_BlockSize - 4))
             );
    }

    //! @brief calculates GDT and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES matrix of coordinates of interest
    //! @param REFERENCE_COORDINATES matrix of reference coordinates
    //! @return pair of GDT and the transformation that superimposes COORDINATES onto REFERENCE_COORDINATES
    storage::Pair< double, Matrix< double> > QualityGDT::CalculateGDTAndSuperimposition
    (
      const Matrix< double> &COORDINATES,
      const Matrix< double> &REFERENCE_COORDINATES
    ) const
    {
      const size_t num_coords( COORDINATES.GetNumberRows());
      // if the coordinates are empty
      if( num_coords < 3)
      {
        BCL_MessageCrt( "less than 3 coordinates: " + util::Format()( num_coords));
        return storage::Pair< double, Matrix< double> >
        (
          util::GetUndefinedDouble(), Matrix< double>( s_NumberRowsTransformationMatrix, m_QualityLCS.GetQualityRMSD().s_BlockSize, m_Queue, 0, m_QualityLCS.GetQualityRMSD().s_BlockSize - 4)
        );
      }

      // seed fragments from lcs
      storage::List< math::Range< size_t> > lcs_ranges( m_QualityLCS.CalculateRanges( COORDINATES, REFERENCE_COORDINATES));
      if( !lcs_ranges.IsEmpty() && lcs_ranges.FirstElement().GetWidth() < 3)
      {
        lcs_ranges.Reset();
      }

      const size_t block_size( m_QualityLCS.GetQualityRMSD().s_BlockSize);
      const size_t number_seed_fragments( num_coords - m_QualityLCS.GetSeedLength() + 1 + lcs_ranges.GetSize());
      const size_t number_seed_fragments_rnd( Tools::RoundUp( block_size, number_seed_fragments));

      Matrix< int> best_selections( num_coords, number_seed_fragments_rnd, m_Queue);
      Matrix< int> start_selections( GetSeedSelections( num_coords, number_seed_fragments_rnd));
      // insert lcs seeds
      {
        size_t fragment_nr( num_coords - m_QualityLCS.GetSeedLength() + 1);
        for
        (
          storage::List< math::Range< size_t> >::const_iterator itr( lcs_ranges.Begin()), itr_end( lcs_ranges.End());
          itr != itr_end;
          ++itr, ++fragment_nr
        )
        {
          start_selections.Fill( 1, itr->GetMin(), itr->GetWidth() + 1, fragment_nr, 1);
        }
      }

      // centers
      // transformation matrices
      Matrix< double> start_transformation_matrices( s_NumberRowsTransformationMatrix * number_seed_fragments, block_size, m_Queue);
      Matrix< double> best_transformation_matrices( s_NumberRowsTransformationMatrix * number_seed_fragments, block_size, m_Queue);

      // calculate transformations for all the seeds
      CalculateTransformations
      (
        COORDINATES,
        REFERENCE_COORDINATES,
        start_selections,
        start_transformation_matrices,
        number_seed_fragments
      );

      // update selection
      UpdateSelections
      (
        COORDINATES,
        REFERENCE_COORDINATES,
        best_selections,
        start_transformation_matrices,
        number_seed_fragments
      );

      // initial best transformations
      CalculateTransformations
      (
        COORDINATES,
        REFERENCE_COORDINATES,
        best_selections,
        best_transformation_matrices,
        number_seed_fragments
      );

      Vector< int> counts( number_seed_fragments, m_Queue);
      Vector< int> differences( number_seed_fragments, m_Queue);
      size_t nr_iterations( 0);

      size_t best_frag_num( util::GetUndefined< size_t>());
      int best_frag_length( 0);

      while( true)
      {
        ++nr_iterations;

        // update selection
        UpdateSelections
        (
          COORDINATES,
          REFERENCE_COORDINATES,
          start_selections, // will be overwritten with selections for current best transformations
          best_transformation_matrices,
          number_seed_fragments
        );

        // calculate transformation matrices
        // initial best transformations
        CalculateTransformations
        (
          COORDINATES,
          REFERENCE_COORDINATES,
          start_selections,
          start_transformation_matrices, // will be updated with transformation matrices for the updated selections
          number_seed_fragments
        );

        // calculate the differences in the selections
        {
          const size_t block_size( m_QualityLCS.GetQualityRMSD().s_BlockSize);

          // collect new selections
          const cl::NDRange kernel_group_dims( block_size, block_size);
          const cl::NDRange kernel_offset;
          const cl::NDRange kernel_work_size( Tools::RoundUp( block_size, num_coords), number_seed_fragments_rnd);

          // create kernel
          cl_int error_number( CL_SUCCESS);
          cl::Kernel kernel( m_QualityLCS.GetQualityRMSD().GetProgram(), "UpdateSelections", &error_number);

          // set the args values
          error_number  = kernel.setArg(  0, start_transformation_matrices.GetData());
          error_number |= kernel.setArg(  1, start_selections.GetData());
          error_number |= kernel.setArg(  2, best_transformation_matrices.GetData());
          error_number |= kernel.setArg(  3, best_selections.GetData());
          error_number |= kernel.setArg(  4, cl_uint( start_transformation_matrices.GetNumberCols()));
          error_number |= kernel.setArg(  5, cl_uint( number_seed_fragments));
          error_number |= kernel.setArg(  6, cl_uint( num_coords));
          error_number |= kernel.setArg(  7, cl_uint( s_NumberRowsTransformationMatrix));
          error_number |= kernel.setArg(  8, counts.GetData());
          error_number |= kernel.setArg(  9, differences.GetData());
          error_number |= kernel.setArg( 10, block_size * block_size * sizeof( cl_int), 0); //shared memory
          error_number |= kernel.setArg( 11, block_size * block_size * sizeof( cl_int), 0); //shared memory
          BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

          // enqueue the kernel
          error_number = m_Queue.enqueueNDRangeKernel( kernel, kernel_offset, kernel_work_size, kernel_group_dims);
          BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
        }

        const linal::Vector< int> difference_host( differences.GetHostVector());
        const linal::Vector< int> counts_host( counts.GetHostVector());

        size_t num_improvements( 0);
        for( size_t frag_num( 0); frag_num < number_seed_fragments; ++frag_num)
        {
          if( best_frag_length < counts_host( frag_num))
          {
            best_frag_length = counts_host( frag_num);
            best_frag_num = frag_num;
          }
          if( difference_host( frag_num) > 0)
          {
            ++num_improvements;
          }
        }

//        BCL_MessageStd( "difference:\n" + util::Format()( difference_host));
        BCL_MessageDbg( "counts:\n" + util::Format()( counts_host));
        BCL_MessageDbg( "iteration: " + util::Format()( nr_iterations) + "\tbest_frag: " + util::Format()( best_frag_num) + "\tlength: " + util::Format()( best_frag_length));
        if( num_improvements == 0)
        {
          break;
        }
      }

      // return the pair of GDT value and the corresponding transformation matrix
      if( util::IsDefined( best_frag_num))
      {
        return storage::Pair< double, Matrix< double> >( double( best_frag_length) / num_coords * OptimalValue(), best_transformation_matrices.SubMatrix( best_frag_num * s_NumberRowsTransformationMatrix, s_NumberRowsTransformationMatrix));
      }
      else
      {
        return storage::Pair< double, Matrix< double> >( 0.0, Matrix< double>( s_NumberRowsTransformationMatrix, block_size, m_Queue));
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &QualityGDT::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_QualityLCS, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &QualityGDT::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_QualityLCS, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief returns all possible seed ranges as selections in cols (not necessarily valid)
    //! @param NUMBER_OF_FRAGMENTS number of coordinates and rows in the selection matrix
    //! @param NUMBER_OF_FRAGMENTS_RND rounded number of fragments - number of cols in selections matrix
    //! @return Matrix of selections for coordinate vectors
    Matrix< int> QualityGDT::GetSeedSelections
    (
      const size_t NUMBER_OF_COORDINATES,
      const size_t NUMBER_OF_FRAGMENTS_RND
    ) const
    {
      Matrix< int> template_selection( NUMBER_OF_COORDINATES, NUMBER_OF_FRAGMENTS_RND, m_Queue);

      // kernel
      // Create the kernel
      cl_int error_number( CL_SUCCESS);
      cl::Kernel kernel( GetTools().GetBufferProgram( util::CPPDataTypes::DataTypeFromTemplate< int>(), m_Queue), "FillSeeds", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      const size_t block_size( 16);

      const cl::NDRange block_dimensions( block_size, block_size);
      const cl::NDRange offset;
      const cl::NDRange kernel_worksize( Tools::RoundUp( block_size, NUMBER_OF_COORDINATES), Tools::RoundUp( block_size, m_QualityLCS.GetSeedLength()));

      error_number  = kernel.setArg( 0, template_selection.GetData());
      error_number |= kernel.setArg( 1, cl_uint( m_QualityLCS.GetSeedLength()));
      error_number |= kernel.setArg( 2, cl_uint( NUMBER_OF_FRAGMENTS_RND));
      error_number |= kernel.setArg( 3, cl_uint( NUMBER_OF_COORDINATES));
      error_number |= kernel.setArg( 4, cl_uint( 1));
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, kernel_worksize, block_dimensions, NULL, NULL);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      return template_selection;
    }

    //! @param COORDINATES matrix of coordinates of interest
    //! @param REFERENCE_COORDINATES matrix of reference coordinates
    //! @param SELECTIONS matrix of elements indicating which rows to use for each fragment
    //! @param TRANSFORMATIONS transformation matrices for all selections
    //! @param NUMBER_OF_FRAGMENTS number of fragments
    void QualityGDT::CalculateTransformations
    (
      const Matrix< double> &COORDINATES,
      const Matrix< double> &REFERENCE_COORDINATES,
      const Matrix< int>    &SELECTIONS,
            Matrix< double> &TRANSFORMATIONS,
      const size_t NUMBER_OF_FRAGMENTS
    ) const
    {
      const size_t block_size( m_QualityLCS.GetQualityRMSD().s_BlockSize);
      const size_t number_seed_fragments_rnd( SELECTIONS.GetNumberCols());
      const size_t num_coords( COORDINATES.GetNumberRows());
      // centers
      Matrix< double> centers_coords( NUMBER_OF_FRAGMENTS, block_size, m_Queue);
      Matrix< double> centers_ref_coords( NUMBER_OF_FRAGMENTS, block_size, m_Queue);

      // calculate transformations for all the seeds
      // centers
      {
        const cl::NDRange kernel_group_dims( block_size, block_size);
        const cl::NDRange kernel_offset;
        const cl::NDRange kernel_work_size( block_size, number_seed_fragments_rnd);

        // create kernel
        cl_int error_number( CL_SUCCESS);
        cl::Kernel kernel_coords( m_QualityLCS.GetQualityRMSD().GetProgram(), "CoordinatesCenterSelections", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
        cl::Kernel kernel_ref_coords( m_QualityLCS.GetQualityRMSD().GetProgram(), "CoordinatesCenterSelections", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // coordinates to consider
        error_number = kernel_coords.setArg( 0, COORDINATES.GetData());
        error_number |= kernel_ref_coords.setArg( 0, REFERENCE_COORDINATES.GetData());
        error_number |= kernel_coords.setArg( 1, SELECTIONS.GetData());
        error_number |= kernel_ref_coords.setArg( 1, SELECTIONS.GetData());

        // num cols
        error_number |= kernel_coords.setArg( 2, cl_uint( COORDINATES.GetNumberCols()));
        error_number |= kernel_ref_coords.setArg( 2, cl_uint( COORDINATES.GetNumberCols()));

        // num coordinates
        error_number |= kernel_coords.setArg( 3, cl_uint( num_coords));
        error_number |= kernel_ref_coords.setArg( 3, cl_uint( num_coords));

        // centers output
        error_number |= kernel_coords.setArg( 4, centers_coords.GetData());
        error_number |= kernel_ref_coords.setArg( 4, centers_ref_coords.GetData());
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // enqueue the kernel
        error_number = m_Queue.enqueueNDRangeKernel( kernel_coords, kernel_offset, kernel_work_size, kernel_group_dims);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
        error_number = m_Queue.enqueueNDRangeKernel( kernel_ref_coords, kernel_offset, kernel_work_size, kernel_group_dims);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      }

      // covariance matrices
      {
        const cl::NDRange kernel_group_dims( block_size, block_size, 1);
        const cl::NDRange kernel_offset;
        const cl::NDRange kernel_work_size( block_size, block_size, NUMBER_OF_FRAGMENTS);

        // create kernel
        cl_int error_number( CL_SUCCESS);
        cl::Kernel kernel( m_QualityLCS.GetQualityRMSD().GetProgram(), "BuildCovarianceMatrixSelections", &error_number);

        // set the args values
        error_number  = kernel.setArg(  0, COORDINATES.GetData());
        error_number |= kernel.setArg(  1, REFERENCE_COORDINATES.GetData());
        error_number |= kernel.setArg(  2, SELECTIONS.GetData());
        error_number |= kernel.setArg(  3, cl_uint( COORDINATES.GetNumberCols()));
        error_number |= kernel.setArg(  4, cl_uint( COORDINATES.GetNumberRows()));
        error_number |= kernel.setArg(  5, cl_uint( SELECTIONS.GetNumberCols()));
        error_number |= kernel.setArg(  6, centers_coords.GetData());
        error_number |= kernel.setArg(  7, centers_ref_coords.GetData());
        error_number |= kernel.setArg(  8, TRANSFORMATIONS.GetData());
        error_number |= kernel.setArg(  9, cl_uint( s_NumberRowsTransformationMatrix));
        error_number |= kernel.setArg( 10, block_size * block_size * sizeof( double), 0); //shared memory
        error_number |= kernel.setArg( 11, block_size * block_size * sizeof( double), 0); //shared memory
        error_number |= kernel.setArg( 12, block_size * sizeof( double), 0); //shared memory
        error_number |= kernel.setArg( 13, block_size * sizeof( double), 0); //shared memory
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // enqueue the kernel
        error_number = m_Queue.enqueueNDRangeKernel( kernel, kernel_offset, kernel_work_size, kernel_group_dims);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      }

      // transformations from covariance matrices
      {
        const cl::NDRange kernel_group_dims( block_size * block_size);
        const cl::NDRange kernel_offset;
        const cl::NDRange kernel_work_size( Tools::RoundUp( block_size * block_size, NUMBER_OF_FRAGMENTS));

        // create kernel
        cl_int error_number( CL_SUCCESS);
        cl::Kernel kernel( Operations< double>::GetInstance().GetProgram(), "TransformationsFromCovarianceMatrices", &error_number);

        // set the args values
        error_number  = kernel.setArg(  0, TRANSFORMATIONS.GetData());
        error_number |= kernel.setArg(  1, centers_coords.GetData());
        error_number |= kernel.setArg(  2, centers_ref_coords.GetData());
        error_number |= kernel.setArg(  3, cl_uint( COORDINATES.GetNumberCols()));
        error_number |= kernel.setArg(  4, cl_uint( s_NumberRowsTransformationMatrix));
        error_number |= kernel.setArg(  5, cl_uint( NUMBER_OF_FRAGMENTS));
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // enqueue the kernel
        error_number = m_Queue.enqueueNDRangeKernel( kernel, kernel_offset, kernel_work_size, kernel_group_dims);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      }
    }

    //! @brief update the selection with the given transformation
    //! @param COORDINATES matrix of coordinates of interest
    //! @param REFERENCE_COORDINATES matrix of reference coordinates
    //! @param SELECTIONS matrix of elements indicating which rows to use for each fragment
    //! @param TRANSFORMATIONS transformation matrices for all selections
    //! @param NUMBER_OF_FRAGMENTS number of fragments
    void QualityGDT::UpdateSelections
    (
      const Matrix< double> &COORDINATES,
      const Matrix< double> &REFERENCE_COORDINATES,
            Matrix< int>    &SELECTIONS,
      const Matrix< double> &TRANSFORMATIONS,
      const size_t           NUMBER_OF_FRAGMENTS
    ) const
    {
      const size_t block_size( m_QualityLCS.GetQualityRMSD().s_BlockSize);
      const size_t num_coords( COORDINATES.GetNumberRows());
      // collect new selections
      const cl::NDRange kernel_group_dims( block_size * s_NumberRowsTransformationMatrix, 1);
      const cl::NDRange kernel_offset;
      const cl::NDRange kernel_work_size( Tools::RoundUp( block_size * s_NumberRowsTransformationMatrix, num_coords), NUMBER_OF_FRAGMENTS);

      // create kernel
      cl_int error_number( CL_SUCCESS);
      cl::Kernel kernel( m_QualityLCS.GetQualityRMSD().GetProgram(), "CoordinateSelectionsBelowCutoff", &error_number);

      // set the args values
      error_number  = kernel.setArg( 0, COORDINATES.GetData());
      error_number |= kernel.setArg( 1, REFERENCE_COORDINATES.GetData());
      error_number |= kernel.setArg( 2, SELECTIONS.GetData());
      error_number |= kernel.setArg( 3, cl_uint( COORDINATES.GetNumberCols()));
      error_number |= kernel.setArg( 4, cl_uint( num_coords));
      error_number |= kernel.setArg( 5, cl_uint( SELECTIONS.GetNumberCols()));
      error_number |= kernel.setArg( 6, TRANSFORMATIONS.GetData());
      error_number |= kernel.setArg( 7, cl_uint( s_NumberRowsTransformationMatrix));
      error_number |= kernel.setArg( 8, math::Sqr( m_QualityLCS.GetCutoff()));
      error_number |= kernel.setArg( 9, block_size * s_NumberRowsTransformationMatrix * sizeof( double), 0); //shared memory
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // enqueue the kernel
      error_number = m_Queue.enqueueNDRangeKernel( kernel, kernel_offset, kernel_work_size, kernel_group_dims);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
    }

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_quality_lcs.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_range.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "opencl/bcl_opencl_matrix.h"
#include "opencl/bcl_opencl_vector.h"
#include "quality/bcl_quality_measures.h"
#include "quality/bcl_quality_superimpose_measures.h"
#include "storage/bcl_storage_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

    //! @brief handler class for adding the quality superimpose enum handler
    class BCL_API LCSEnumHandler :
      public signal::Slots
    {

    private:

    //////////
    // data //
    //////////

      //! the enum in the quality::SuperimposeMeasure
      quality::SuperimposeMeasure e_LCSSuperImposeMeasure;
      //! the enum in the quality::Measure
      quality::Measure e_LCSMeasure;

      //! the only instance of this class
      static const LCSEnumHandler s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LCSEnumHandler() :
        e_LCSSuperImposeMeasure( quality::GetSuperimposeMeasures().AddEnum( "OpenclLCS", util::ShPtr< QualityLCS>())),
        e_LCSMeasure( quality::GetMeasures().AddEnum( "OpenclLCS", util::ShPtr< QualityLCS>()))
      {
        // register enum with opencl queue update signal
        GetTools().GetQueueUpdateSignal().Connect( this, &LCSEnumHandler::UpdateEnum);
      }

      //! @brief update the enum with the command queue from the Tools
      //! @param TOOLS the tolls to get the commandqueue from
      void UpdateEnum( Tools &TOOLS)
      {
        util::ShPtr< QualityLCS> sp_quality( new QualityLCS());
        if( !TOOLS.HasCommandQueues())
        {
          *e_LCSSuperImposeMeasure = util::ShPtr< QualityLCS>();
          *e_LCSMeasure = util::ShPtr< QualityLCS>();
          return;
        }

        // try to initialize
        if( sp_quality->Initialize( TOOLS.GetFirstCommandQueue()))
        {
          // just update the existing one with the new one
          *e_LCSSuperImposeMeasure = sp_quality;
          *e_LCSMeasure = sp_quality;
        }
        else
        {
          BCL_MessageVrb( "unable to initialize enum: OpenclLCS");
        }
      }

    }; // class LCSEnumHandler

    //! instance of DensitySimulateEnumHandler
    const LCSEnumHandler LCSEnumHandler::s_Instance = LCSEnumHandler();

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from a RMSD cutoff and a seed length
    //! @param RMSD_CUTOFF distance cutoff
    //! @param SEED_LENGTH length of seeds
    QualityLCS::QualityLCS
    (
      const double RMSD_CUTOFF,
      const size_t SEED_LENGTH
    ) :
      m_Queue(),
      m_QualityRmsd(),
      m_RmsdCutoff( RMSD_CUTOFF),
      m_SeedLength( SEED_LENGTH)
    {
    }

    //! @brief construct from a RMSD cutoff and a seed length and queue
    //! @param RMSD_CUTOFF distance cutoff
    //! @param SEED_LENGTH length of seeds
    //! @param QUEUE command queue to use
    QualityLCS::QualityLCS
    (
      const CommandQueue &QUEUE,
      const double RMSD_CUTOFF,
      const size_t SEED_LENGTH
    ) :
      m_Queue( QUEUE),
      m_QualityRmsd( m_Queue),
      m_RmsdCutoff( RMSD_CUTOFF),
      m_SeedLength( SEED_LENGTH)
    {
    }

    //! @brief Clone function
    //! @return pointer to new LCS
    QualityLCS *QualityLCS::Clone() const
    {
      return new QualityLCS( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &QualityLCS::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the optimal value for that quality measurement
    //! @return the best value by which two sets of coordinates can agree
    double QualityLCS::OptimalValue() const
    {
      return util::GetUndefined< double>();
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &QualityLCS::GetComparisonFunction() const
    {
      return **math::Comparisons< double>::GetEnums().e_Greater;
    }

    //! @brief return rmsd cutoff
    //! @return rmsd cutoff
    double QualityLCS::GetCutoff() const
    {
      return m_RmsdCutoff;
    }

    //! @brief get seed length
    //! @return seed length
    size_t QualityLCS::GetSeedLength() const
    {
      return m_SeedLength;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief is this class compatible with given command queue
    //! @param COMMAND_QUEUE the command queue this object would operate on
    //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
    bool QualityLCS::IsCompatible( const CommandQueue &COMMAND_QUEUE) const
    {
      return m_QualityRmsd.IsCompatible( COMMAND_QUEUE);
    }

    //! @brief initialize this class
    //! @brief COMMAND_QUEUE queue to use
    //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
    bool QualityLCS::Initialize( const CommandQueue &COMMAND_QUEUE)
    {
      if( m_QualityRmsd.Initialize( COMMAND_QUEUE))
      {
        m_Queue = COMMAND_QUEUE;
        return true;
      }
      m_Queue = CommandQueue();
      return false;
    }

    //! @brief calculates LCS between given coordinates
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return LCS between COORDINATES and REFERENCE_COORDINATES
    double QualityLCS::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      return CalculateMeasure
             (
               m_QualityRmsd.MatrixFromCoordinates( COORDINATES, m_QualityRmsd.s_BlockSize),
               m_QualityRmsd.MatrixFromCoordinates( REFERENCE_COORDINATES, m_QualityRmsd.s_BlockSize)
             );
    }

    //! @brief calculates root mean square deviation between given coordinates
    //! @param COORDINATES matrix of coordinates of interest
    //! @param REFERENCE_COORDINATES matrix of reference coordinates
    //! @return root mean square deviation between given coordinates
    double QualityLCS::CalculateMeasure
    (
      const Matrix< double> &COORDINATES,
      const Matrix< double> &REFERENCE_COORDINATES
    ) const
    {
      // calculate the LCS and store the range for the first one
      const storage::List< math::Range< size_t> > longest_ranges
      (
        CalculateRanges( COORDINATES, REFERENCE_COORDINATES)
      );
      if( longest_ranges.IsEmpty())
      {
        return 0.0;
      }

      // take first of the longest ranges, since all are equally long, they might just be at different ranges
      const math::Range< size_t> lcs( longest_ranges.FirstElement());

      // return the length of the segment
      return lcs.GetWidth() + 1;
    }

    //! @brief calculates the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
    math::TransformationMatrix3D QualityLCS::CalculateSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      return CalculateSuperimposition
             (
               m_QualityRmsd.MatrixFromCoordinates( COORDINATES, m_QualityRmsd.s_BlockSize),
               m_QualityRmsd.MatrixFromCoordinates( REFERENCE_COORDINATES, m_QualityRmsd.s_BlockSize)
             );
    }

    //! @brief calculates the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES matrix of coordinates of interest
    //! @param REFERENCE_COORDINATES matrix of reference coordinates
    //! @return the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
    math::TransformationMatrix3D QualityLCS::CalculateSuperimposition
    (
      const Matrix< double> &COORDINATES,
      const Matrix< double> &REFERENCE_COORDINATES
    ) const
    {
      // calculate the LCS and store the range for the first one
      const storage::List< math::Range< size_t> > longest_ranges( CalculateRanges( COORDINATES, REFERENCE_COORDINATES));
      if( longest_ranges.IsEmpty())
      {
        return math::TransformationMatrix3D( util::UndefinedObject());
      }

      // take first of the longest ranges, since all are equally long, they might just be at different ranges
      const math::Range< size_t> &lcs( longest_ranges.FirstElement());

      // calculate and return the transformation
      return m_QualityRmsd.CalculateSuperimposition
             (
               COORDINATES.SubMatrix( lcs.GetMin(), lcs.GetWidth() + 1),
               REFERENCE_COORDINATES.SubMatrix( lcs.GetMin(), lcs.GetWidth() + 1)
             );
    }

    //! @brief returns the ranges of longest continuous segments that can be superimposed below cutoff for given coordinates
    //! @param COORDINATES matrix of coordinates of interest
    //! @param REFERENCE_COORDINATES matrix of reference coordinates
    //! @return the range of longest continuous segment that can be superimposed below cutoff for given coordinates
    storage::List< math::Range< size_t> > QualityLCS::CalculateRanges
    (
      const Matrix< double> &COORDINATES,
      const Matrix< double> &REFERENCE_COORDINATES
    ) const
    {
      // initialize variable to store the longest fragments
      storage::List< math::Range< size_t> > longest_fragments;

      // consider all fragment lengths
      for
      (
        size_t fragment_length( COORDINATES.GetNumberRows());
        fragment_length >= 3 && longest_fragments.IsEmpty();
        --fragment_length
      )
      {
        const linal::Vector< double> rmsds
        (
          RMSDOfEachFragment( COORDINATES, REFERENCE_COORDINATES, fragment_length).GetHostVector()
        );

        // iterate through vector
        for( size_t index( 0); index < COORDINATES.GetNumberRows() + 1 - fragment_length; ++index)
        {
          if( rmsds( index) < m_RmsdCutoff)
          {
            longest_fragments.PushBack( math::Range< size_t>( index, index + fragment_length - 1));
          }
        }
      }

      // end
      return longest_fragments;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &QualityLCS::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_RmsdCutoff, ISTREAM);
      io::Serialize::Read( m_SeedLength, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &QualityLCS::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_RmsdCutoff, OSTREAM, INDENT) << '\t';
      io::Serialize::Write( m_SeedLength, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief check if a given range superimposes below the cutoff
    //! @param COORDINATES matrix of coordinates of interest
    //! @param REFERENCE_COORDINATES matrix of reference coordinates
    //! @return true, if coordinates within the given range are superimposable below the cutoff
    bool QualityLCS::IsGoodRange
    (
      const math::Range< size_t> &RANGE,
      const Matrix< double> &COORDINATES,
      const Matrix< double> &REFERENCE_COORDINATES
    ) const
    {
      // if the current RMSD is smaller
      return m_QualityRmsd.CalculateMeasure
             (
               COORDINATES.SubMatrix( RANGE.GetMin(), RANGE.GetWidth() + 1),
               REFERENCE_COORDINATES.SubMatrix( RANGE.GetMin(), RANGE.GetWidth() + 1)
             ) < m_RmsdCutoff;
    }

    //! @brief find a larger range by extending the given one that has a RMSD below the cutoff
    //! @param RANGE range of coordinates that are used to be as seed and to be extended
    //! @param COORDINATES matrix of coordinates of interest
    //! @param REFERENCE_COORDINATES matrix of reference coordinates
    //! @return extended range
    math::Range< size_t> QualityLCS::ExtendRange
    (
      const math::Range< size_t> &RANGE,
      const Matrix< double> &COORDINATES,
      const Matrix< double> &REFERENCE_COORDINATES
    ) const
    {
      // make a copy of the range
      math::Range< size_t> longest_range( RANGE);
      math::Range< size_t> current_range( RANGE);

      // while the range can still be extended
      while( current_range.GetMax() < COORDINATES.GetNumberRows() - 1)
      {
        // increment the current range length
        current_range.SetMax( current_range.GetMax() + 1);

        // if the current RMSD is smaller
        if
        (
          IsGoodRange
          (
            current_range,
            COORDINATES,
            REFERENCE_COORDINATES
          )
        )
        {
          // update the longest range
          longest_range = current_range;
        }
      }

      // return the longest range found so far
      return longest_range;
    }

    //! @brief returns the indices to the coordinates of longest continuous segments that can be superimposed below cutoff for given coordinates
    //! @param COORDINATES matrix of coordinates of interest
    //! @param REFERENCE_COORDINATES matrix of reference coordinates
    //! @return the indices of longest continuous segments that can be superimposed below cutoff for given coordinates
    storage::List< storage::List< size_t> > QualityLCS::CalculateIndices
    (
      const Matrix< double> &COORDINATES,
      const Matrix< double> &REFERENCE_COORDINATES
    ) const
    {
      return quality::LCS::ConvertRangesToLists( CalculateRanges( COORDINATES, REFERENCE_COORDINATES));
    }

    //! @brief calculate RMSDs for all fragments of given fragment length
    //! @param COORDINATES matrix of coordinates of interest
    //! @param REFERENCE_COORDINATES matrix of reference coordinates
    //! @param FRAGMENT_LENGTH number of coordinates in each fragment
    Vector< double> QualityLCS::RMSDOfEachFragment
    (
      const Matrix< double> &COORDINATES,
      const Matrix< double> &REFERENCE_COORDINATES,
      const size_t FRAGMENT_LENGTH
    ) const
    {
      const size_t number_of_fragments( COORDINATES.GetNumberRows() - FRAGMENT_LENGTH + 1);
      const size_t number_of_fragments_rnd( Tools::RoundUp( m_QualityRmsd.s_BlockSize, number_of_fragments));

      // centers
      Matrix< double> centers_coords( number_of_fragments, m_QualityRmsd.s_BlockSize, m_Queue);
      Matrix< double> centers_ref_coords( number_of_fragments, m_QualityRmsd.s_BlockSize, m_Queue);

      // calculate the centers
      {
        const cl::NDRange kernel_group_dims( m_QualityRmsd.s_BlockSize, m_QualityRmsd.s_BlockSize);
        const cl::NDRange kernel_offset;
        const cl::NDRange kernel_work_size( m_QualityRmsd.s_BlockSize, number_of_fragments_rnd);

        // create kernel
        cl_int error_number( CL_SUCCESS);
        cl::Kernel kernel_coords( m_QualityRmsd.GetProgram(), "FragmentCenters", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
        cl::Kernel kernel_ref_coords( m_QualityRmsd.GetProgram(), "FragmentCenters", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // coordinates to consider
        error_number = kernel_coords.setArg( 0, COORDINATES.GetData());
        error_number |= kernel_ref_coords.setArg( 0, REFERENCE_COORDINATES.GetData());

        // fragment length
        error_number |= kernel_coords.setArg( 1, cl_uint( FRAGMENT_LENGTH));
        error_number |= kernel_ref_coords.setArg( 1, cl_uint( FRAGMENT_LENGTH));
        error_number |= kernel_coords.setArg( 2, cl_uint( number_of_fragments));
        error_number |= kernel_ref_coords.setArg( 2, cl_uint( number_of_fragments));

        // output
        error_number |= kernel_coords.setArg( 3, centers_coords.GetData());
        error_number |= kernel_ref_coords.setArg( 3, centers_ref_coords.GetData());
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // enqueue the kernel
        error_number = m_Queue.enqueueNDRangeKernel( kernel_coords, kernel_offset, kernel_work_size, kernel_group_dims);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
        error_number = m_Queue.enqueueNDRangeKernel( kernel_ref_coords, kernel_offset, kernel_work_size, kernel_group_dims);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      }

      // build covariance matrix for each fragment
      // square norm of centers
      Vector< double> square_norm_centered_coords( number_of_fragments, m_Queue);
      Vector< double> square_norm_centered_ref_coords( number_of_fragments, m_Queue);

      // covariance matrix
      const size_t num_rows_cov( 3);
      Matrix< double> covariance_matrices( 3 * number_of_fragments, m_QualityRmsd.s_BlockSize, m_Queue);
      {
        const cl::NDRange kernel_group_dims( m_QualityRmsd.s_BlockSize, m_QualityRmsd.s_BlockSize, 1);
        const cl::NDRange kernel_offset;
        const cl::NDRange kernel_work_size( m_QualityRmsd.s_BlockSize, m_QualityRmsd.s_BlockSize, number_of_fragments);

        // create kernel
        cl_int error_number( CL_SUCCESS);
        cl::Kernel kernel( m_QualityRmsd.GetProgram(), "BuildCovarianceMatrixFragments", &error_number);

        // set the args values
        error_number  = kernel.setArg(  0, COORDINATES.GetData());
        error_number |= kernel.setArg(  1, REFERENCE_COORDINATES.GetData());
        error_number |= kernel.setArg(  2, centers_coords.GetData());
        error_number |= kernel.setArg(  3, centers_ref_coords.GetData());
        error_number |= kernel.setArg(  4, cl_uint( COORDINATES.GetNumberCols()));
        error_number |= kernel.setArg(  5, cl_uint( FRAGMENT_LENGTH));
        error_number |= kernel.setArg(  6, covariance_matrices.GetData());
        error_number |= kernel.setArg(  7, cl_uint( num_rows_cov));
        error_number |= kernel.setArg(  8, square_norm_centered_coords.GetData());
        error_number |= kernel.setArg(  9, square_norm_centered_ref_coords.GetData());
        error_number |= kernel.setArg( 10, m_QualityRmsd.s_BlockSize * m_QualityRmsd.s_BlockSize * sizeof( double), 0); //shared memory
        error_number |= kernel.setArg( 11, m_QualityRmsd.s_BlockSize * m_QualityRmsd.s_BlockSize * sizeof( double), 0); //shared memory
        error_number |= kernel.setArg( 12, m_QualityRmsd.s_BlockSize * sizeof( double), 0); //shared memory
        error_number |= kernel.setArg( 13, m_QualityRmsd.s_BlockSize * sizeof( double), 0); //shared memory
        error_number |= kernel.setArg( 14, m_QualityRmsd.s_BlockSize * m_QualityRmsd.s_BlockSize * sizeof( double), 0); //shared memory
        error_number |= kernel.setArg( 15, m_QualityRmsd.s_BlockSize * m_QualityRmsd.s_BlockSize * sizeof( double), 0); //shared memory
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // enqueue the kernel
        error_number = m_Queue.enqueueNDRangeKernel( kernel, kernel_offset, kernel_work_size, kernel_group_dims);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      }

      // calculate the rmsds from the covariance matrices
      Vector< double> rmsds( number_of_fragments, m_Queue);
      {
        const cl::NDRange kernel_group_dims( m_QualityRmsd.s_BlockSize);
        const cl::NDRange kernel_offset;
        const cl::NDRange kernel_work_size( number_of_fragments_rnd);

        // create kernel
        cl_int error_number( CL_SUCCESS);
        cl::Kernel kernel( m_QualityRmsd.GetProgram(), "RMSDFromCovarianceMatrix", &error_number);

        // set arguments
        error_number  = kernel.setArg( 0, covariance_matrices.GetData());
        error_number |= kernel.setArg( 1, cl_uint( FRAGMENT_LENGTH));
        error_number |= kernel.setArg( 2, cl_uint( number_of_fragments));
        error_number |= kernel.setArg( 3, cl_uint( num_rows_cov));
        error_number |= kernel.setArg( 4, cl_uint( m_QualityRmsd.s_BlockSize));
        error_number |= kernel.setArg( 5, square_norm_centered_coords.GetData());
        error_number |= kernel.setArg( 6, square_norm_centered_ref_coords.GetData());
        error_number |= kernel.setArg( 7, rmsds.GetData());
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // enqueue the kernel
        error_number = m_Queue.enqueueNDRangeKernel( kernel, kernel_offset, kernel_work_size, kernel_group_dims);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      }

      return rmsds;
    }

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_quality_rmsd.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_comparisons.h"
#include "math/bcl_math_transformation_matrix_3d.h"
#include "opencl/bcl_opencl_kernel_sources.h"
#include "opencl/bcl_opencl_matrix.h"
#include "opencl/bcl_opencl_matrix3x3.h"
#include "opencl/bcl_opencl_operations.h"
#include "opencl/bcl_opencl_vector.h"
#include "quality/bcl_quality_measures.h"
#include "quality/bcl_quality_rmsd.h"
#include "quality/bcl_quality_superimpose_measures.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

    //! @brief handler class for adding the quality superimpose enum handler
    class BCL_API RMSDEnumHandler :
      public signal::Slots
    {

    private:

    //////////
    // data //
    //////////

      //! the enum in the quality::SuperimposeMeasure
      quality::SuperimposeMeasure e_RMSDSuperImposeMeasure;
      //! the enum in the quality::Measure
      quality::Measure e_RMSDMeasure;

      //! the only instance of this class
      static const RMSDEnumHandler s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      RMSDEnumHandler() :
        e_RMSDSuperImposeMeasure( quality::GetSuperimposeMeasures().AddEnum( "OpenclRMSD", util::ShPtr< QualityRMSD>())),
        e_RMSDMeasure( quality::GetMeasures().AddEnum( "OpenclRMSD", util::ShPtr< QualityRMSD>()))
      {
        // register enum with opencl queue update signal
        GetTools().GetQueueUpdateSignal().Connect( this, &RMSDEnumHandler::UpdateEnum);
      }

      //! @brief update the enum with the command queue from the Tools
      //! @param TOOLS the tolls to get the commandqueue from
      void UpdateEnum( Tools &TOOLS)
      {
        util::ShPtr< QualityRMSD> sp_quality( new QualityRMSD());
        if( !TOOLS.HasCommandQueues())
        {
          *e_RMSDSuperImposeMeasure = util::ShPtr< QualityRMSD>();
          *e_RMSDMeasure = util::ShPtr< QualityRMSD>();
          return;
        }

        // try to initialize
        if( sp_quality->Initialize( TOOLS.GetFirstCommandQueue()))
        {
          // just update the existing one with the new one
          *e_RMSDSuperImposeMeasure = sp_quality;
          *e_RMSDMeasure = sp_quality;
        }
        else
        {
          BCL_MessageVrb( "unable to initialize enum: OpenclRMSD");
        }
      }

    }; // class RMSDEnumHandler

    //! instance of DensitySimulateEnumHandler
    const RMSDEnumHandler RMSDEnumHandler::s_Instance = RMSDEnumHandler();

  //////////
  // data //
  //////////

    const size_t QualityRMSD::s_BlockSize = 16;

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a boolean to whether superimpose coordinates (set to true by default)
    //! @param SUPERIMPOSE_COORDINATES boolean to whether superimpose coordinates before calculating RMSD
    QualityRMSD::QualityRMSD
    (
      const bool SUPERIMPOSE_COORDINATES
    ) :
      m_Queue(),
      m_SuperimposeCoordinates( SUPERIMPOSE_COORDINATES)
    {
    }

    //! @brief constructor from a boolean to whether superimpose coordinates (set to true by default)
    //! @param QUEUE command queue
    //! @param SUPERIMPOSE_COORDINATES boolean to whether superimpose coordinates before calculating RMSD
    QualityRMSD::QualityRMSD( const CommandQueue &QUEUE, const bool SUPERIMPOSE_COORDINATES) :
      m_Queue( QUEUE),
      m_SuperimposeCoordinates( SUPERIMPOSE_COORDINATES)
    {
      BCL_Assert( Initialize( m_Queue), "cannot initialize from given command queue");
    }

    //! @brief virtual copy constructor
    //! @return pointer to new QualityRMSD
    QualityRMSD *QualityRMSD::Clone() const
    {
      return new QualityRMSD( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get name of the current class
    //! @return name of the class
    const std::string &QualityRMSD::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the comparison function for better quality
    //! @return binary function to compare two quality measure values
    const util::BinaryFunctionInterface< double, double, bool> &QualityRMSD::GetComparisonFunction() const
    {
      return **math::Comparisons< double>::GetEnums().e_Less;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief is this class compatible with given command queue
    //! @param COMMAND_QUEUE the command queue this object would operate on
    //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
    bool QualityRMSD::IsCompatible( const CommandQueue &COMMAND_QUEUE) const
    {
      cl_int error_number( CL_SUCCESS);
      const Device device( COMMAND_QUEUE.GetDevice( &error_number));

      // can get device
      if( error_number != CL_SUCCESS)
      {
        return false;
      }

      const storage::Set< Extension> extensions( device.Extensions( &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "unable to get extensions from device");

      return KernelSourceInterface::PrecisionCompatibleWithExtensions
             (
               util::CPPDataTypes::DataTypeFromTemplate< double>(),
               extensions
             );
    }

    //! @brief initialize this class
    //! @brief COMMAND_QUEUE queue to use
    //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
    bool QualityRMSD::Initialize( const CommandQueue &COMMAND_QUEUE)
    {
      // check if this is a compatible command queue
      if( !IsCompatible( COMMAND_QUEUE))
      {
        BCL_MessageDbg( "command queue is not compatible");
        return false;
      }

      // update the command queue
      m_Queue = COMMAND_QUEUE;

      // for precision type
      cl_int error_number( CL_SUCCESS);
      m_Program = KernelSources::Compile( GetKernelSources().e_Quality, util::CPPDataTypes::e_Double, m_Queue, std::string(), &error_number);

      if( error_number != CL_SUCCESS)
      {
        BCL_MessageDbg( "error compiling programs:\n" + Tools::ErrorString( error_number));
        return false;
      }

      return true;
    }

    //! @brief calculates root mean square deviation between given coordinates
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return root mean square deviation between given coordinates
    double QualityRMSD::CalculateMeasure
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      return CalculateMeasure
             (
               MatrixFromCoordinates( COORDINATES, s_BlockSize),
               MatrixFromCoordinates( REFERENCE_COORDINATES, s_BlockSize)
             );
    }

    //! @brief calculates root mean square deviation between given coordinates
    //! @param COORDINATES matrix of coordinates of interest
    //! @param REFERENCE_COORDINATES matrix of reference coordinates
    //! @return root mean square deviation between given coordinates
    double QualityRMSD::CalculateMeasure
    (
      const Matrix< double> &COORDINATES,
      const Matrix< double> &REFERENCE_COORDINATES
    ) const
    {
      // if superimpose coordinates is set
      if( m_SuperimposeCoordinates)
      {
        return SuperimposedRMSD( COORDINATES, REFERENCE_COORDINATES);
      }
      // don't superimpose coordinates
      else
      {
        return RealSpaceRMSD( COORDINATES, REFERENCE_COORDINATES);
      }
    }

    //! @brief calculates the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES coordinates of interest
    //! @param REFERENCE_COORDINATES reference coordinates
    //! @return the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
    math::TransformationMatrix3D QualityRMSD::CalculateSuperimposition
    (
      const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
      const util::SiPtrVector< const linal::Vector3D> &REFERENCE_COORDINATES
    ) const
    {
      return CalculateSuperimposition
             (
               MatrixFromCoordinates( COORDINATES, s_BlockSize),
               MatrixFromCoordinates( REFERENCE_COORDINATES, s_BlockSize)
             );
    }

    //! @brief calculates the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
    //! @param COORDINATES matrix of coordinates of interest
    //! @param REFERENCE_COORDINATES matrix of reference coordinates
    //! @return the transformation matrix that superimposes COORDINATES onto REFERENCE_COORDINATES
    math::TransformationMatrix3D QualityRMSD::CalculateSuperimposition
    (
      const Matrix< double> &COORDINATES,
      const Matrix< double> &REFERENCE_COORDINATES
    ) const
    {
      return math::TransformationMatrix3D
            (
              SuperimposeCoordinates
              (
                COORDINATES,
                REFERENCE_COORDINATES
              ).GetHostMatrix( s_BlockSize - 4, s_BlockSize - 4)
            );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &QualityRMSD::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SuperimposeCoordinates, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read from std::ostream
    //! @param OSTREAM input stream
    //! @param INDENT indentation
    //! @return ostream which was read from
    std::ostream &QualityRMSD::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SuperimposeCoordinates, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief determine the transformation matrix to optimally (lowest RMSD) superimpose two sets of coordinates
    //! @param COORDINATES matrix of coordinates of interest
    //! @param REFERENCE_COORDINATES matrix of reference coordinates
    //! @return Transformation matrix that superimposes B onto A
    Matrix< double> QualityRMSD::SuperimposeCoordinates
    (
      const Matrix< double> &REFERENCE_COORDINATES,
      const Matrix< double> &COORDINATES
    ) const
    {
      const Vector< double> center_coordinates( Center( COORDINATES));
      const Vector< double> center_reference_coordinates( Center( REFERENCE_COORDINATES));

      Vector< double> square_norm_centered_a;
      Vector< double> square_norm_centered_b;

      // Calculate the covariance matrix
      const Matrix3x3< double> moment_device
      (
        BuildCovarianceMatrix
        (
          COORDINATES,
          REFERENCE_COORDINATES,
          center_coordinates,
          center_reference_coordinates,
          square_norm_centered_a,
          square_norm_centered_b
        )
      );

      // return transformation matrix calculated
      const math::TransformationMatrix3D transformation
      (
        quality::RMSD::CovarianceToTransformationMatrix
        (
          moment_device.GetHostMatrix(),
          linal::Vector3D( center_coordinates.GetHostVector().Begin()),
          linal::Vector3D( center_reference_coordinates.GetHostVector().Begin())
        )
      );

      return Matrix< double>( linal::Matrix< double>( 4, 4, transformation.GetMatrix().Begin()), m_Queue, s_BlockSize - 4, s_BlockSize - 4);

//      // determine transformation
//      return CovarianceToTransformationMatrix( moment_device, center_coordinates, center_reference_coordinates);
    }

    //! @brief calculate the real space rmsd of two sets of coordinates
    //! uses the coordinates as they are passed
    //! @param COORDINATES matrix of coordinates of interest
    //! @param REFERENCE_COORDINATES matrix of reference coordinates
    //! @return the rmsd between the passed coordinates
    double QualityRMSD::RealSpaceRMSD
    (
      const Matrix< double> &COORDINATES,
      const Matrix< double> &REFERENCE_COORDINATES
    ) const
    {
      // error catching
      cl_int error_number = CL_SUCCESS;

      // setup kernel
      cl::Kernel kernel( m_Program, "RealSpaceRMSD", &error_number);

      // result
      Vector< double> rmsd( 1, m_Queue);

      // arguments
      error_number  = kernel.setArg( 0, COORDINATES.GetData());
      error_number |= kernel.setArg( 1, REFERENCE_COORDINATES.GetData());
      error_number |= kernel.setArg( 2, cl_uint( COORDINATES.GetNumberOfElements()));
      error_number |= kernel.setArg( 3, rmsd.GetData());
      error_number |= kernel.setArg( 4, s_BlockSize * s_BlockSize * sizeof( double), 0); // shared memory
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // launch kernel
      const cl::NDRange block_dims( s_BlockSize * s_BlockSize);
      const cl::NDRange offset;
      const cl::NDRange worksize( s_BlockSize * s_BlockSize);

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      return math::Sqrt( rmsd.GetHostVector()( 0) / COORDINATES.GetNumberRows());
    }

    //! @brief calculate the rmsd of two sets of coordinates if they are optimally
    //! they are never superimposed - the rmsd is derived from the eigenvalues of the SVD of the covariance matrix
    //! @param COORDINATES matrix of coordinates of interest
    //! @param REFERENCE_COORDINATES matrix of reference coordinates
    //! @return the rmsd of the coordinates
    double QualityRMSD::SuperimposedRMSD
    (
      const Matrix< double> &COORDINATES,
      const Matrix< double> &REFERENCE_COORDINATES
    ) const
    {
      const Vector< double> center_a( Center( COORDINATES));
      const Vector< double> center_b( Center( REFERENCE_COORDINATES));

      Vector< double> square_norm_centered_a;
      Vector< double> square_norm_centered_b;

      // Calculate the covariance matrix
      Matrix3x3< double> moment_device
      (
        BuildCovarianceMatrix
        (
          COORDINATES,
          REFERENCE_COORDINATES,
          center_a,
          center_b,
          square_norm_centered_a,
          square_norm_centered_b
        )
      );

//      if( !moment.IsDefined())
//      {
//        BCL_MessageCrt( "covariance matrix is undefined");
//
//        return util::GetUndefinedDouble();
//      }

      static const double s_chi_threshold( 1e-10);

      // determine sign of last element
      const Vector< double> determinant( moment_device.Determinant());

      moment_device.MultiplyWithTransposed();
      // sort diagonal
      linal::Vector< double> eigenvalues( moment_device.EigenValues( true).GetHostVector( s_BlockSize - 3));
      std::sort( eigenvalues.Begin(), eigenvalues.End());

      const int chi( determinant( 0) < s_chi_threshold ? -1 : 1);
      eigenvalues( 0) *= chi;

      // calculate the square deviation
      double square_deviation( 0.0);
      square_deviation += square_norm_centered_a( 0);
      square_deviation += square_norm_centered_b( 0);
      square_deviation -= 2 * eigenvalues.Sum();

      // root mean and return
      return math::Sqrt( std::max( square_deviation, double( 0)) / double( COORDINATES.GetNumberRows()));
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief create a padded matrix for a vector of coordinates
    //! @param COORDINATES the coordinates
    //! @param BLOCK_SIZE block size for kernels
    //! @return Matrix that are padded to desired blocksize and contain the coordinates as rows
    Matrix< double> QualityRMSD::MatrixFromCoordinates( const util::SiPtrVector< const linal::Vector3D> &COORDINATES, const size_t BLOCK_SIZE) const
    {
      const size_t num_vectors( COORDINATES.GetSize());
      const size_t vector_size( 3);

      linal::Matrix< double> input_matrix( num_vectors, vector_size);
      double *ptr( input_matrix.Begin());

      for( size_t ctr( 0); ctr < num_vectors; ++ctr)
      {
        const double *vec_itr( COORDINATES( ctr)->Begin());
        for( size_t vec( 0); vec < vector_size; ++vec, ++ptr, ++vec_itr)
        {
          ( *ptr) = ( *vec_itr);
        }
      }

      const size_t pad_col( ( BLOCK_SIZE - ( input_matrix.GetNumberCols() % BLOCK_SIZE)) % BLOCK_SIZE);
      return Matrix< double>( input_matrix, m_Queue, 0, pad_col);
    }

    //! @brief compute covariance matrix of two sets of coordinates COORDINATES_A on COORDINATES_B
    //! both coordinate sets are translated to the center of mass
    //! @param COORDINATES matrix of coordinates of interest
    //! @param REFERENCE_COORDINATES matrix of reference coordinates
    //! @param CENTER_A the center of COORDINATES_A
    //! @param CENTER_B the center of COORDINATES_B
    //! @param SQUARE_NORM_CENTERED_COORDINATES_A optional pointer to which the square norm of the centered coordinates a will be depsosited
    //! @param SQUARE_NORM_CENTERED_COORDINATES_B optional pointer to which the square norm of the centered coordinates b will be depsosited
    //! @return COORDINATES_A * COORDINATES_B
    Matrix3x3< double> QualityRMSD::BuildCovarianceMatrix
    (
      const Matrix< double> &COORDINATES,
      const Matrix< double> &REFERENCE_COORDINATES,
      const Vector< double> &CENTER_A,
      const Vector< double> &CENTER_B,
      Vector< double> &SQUARE_NORM_CENTERED_COORDINATES_A,
      Vector< double> &SQUARE_NORM_CENTERED_COORDINATES_B
    ) const
    {
      // setup kernel
      // error catching
      cl_int error_number = CL_SUCCESS;

      BCL_Assert
      (
        COORDINATES.GetNumberCols() % s_BlockSize == 0 &&
        REFERENCE_COORDINATES.GetNumberCols() % s_BlockSize == 0,
        "input buffers are not padded correctly!"
      );

      // output
      Matrix3x3< double> covariance_matrix( m_Queue, double( 0));
      SQUARE_NORM_CENTERED_COORDINATES_A = Vector< double>( 1, m_Queue);
      SQUARE_NORM_CENTERED_COORDINATES_B = Vector< double>( 1, m_Queue);

      cl::Kernel kernel( m_Program, "BuildCovarianceMatrix", &error_number);

      // set the args values
      error_number  = kernel.setArg(  0, COORDINATES.GetData());
      error_number |= kernel.setArg(  1, REFERENCE_COORDINATES.GetData());
      error_number |= kernel.setArg(  2, cl_uint( COORDINATES.GetNumberRows()));
      error_number |= kernel.setArg(  3, cl_uint( COORDINATES.GetNumberCols()));
      error_number |= kernel.setArg(  4, CENTER_A.GetData());
      error_number |= kernel.setArg(  5, CENTER_B.GetData());
      error_number |= kernel.setArg(  6, covariance_matrix.GetData());
      error_number |= kernel.setArg(  7, SQUARE_NORM_CENTERED_COORDINATES_A.GetData());
      error_number |= kernel.setArg(  8, SQUARE_NORM_CENTERED_COORDINATES_B.GetData());
      error_number |= kernel.setArg(  9, s_BlockSize * s_BlockSize * sizeof( double), 0); //shared memory
      error_number |= kernel.setArg( 10, s_BlockSize * s_BlockSize * sizeof( double), 0); //shared memory
      error_number |= kernel.setArg( 11, s_BlockSize * sizeof( double), 0); //shared memory
      error_number |= kernel.setArg( 12, s_BlockSize * sizeof( double), 0); //shared memory
      error_number |= kernel.setArg( 13, s_BlockSize * s_BlockSize * sizeof( double), 0); //shared memory
      error_number |= kernel.setArg( 14, s_BlockSize * s_BlockSize * sizeof( double), 0); //shared memory
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Execute
      const cl::NDRange block_dims( s_BlockSize, s_BlockSize);
      const cl::NDRange offset;
      const cl::NDRange worksize( s_BlockSize, s_BlockSize);

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // end
      return covariance_matrix;
    }

    //! @brief calculate the center of a given matrix of coordinates
    //! @param COORDINATES matrix of coordinates in rows
    //! @return the center as vector
    Vector< double> QualityRMSD::Center
    (
      const Matrix< double> &COORDINATES
    ) const
    {
      Vector< double> center( COORDINATES.GetNumberCols(), m_Queue, 3 % s_BlockSize);

      cl_int error_number = CL_SUCCESS;

      // setup kernel
      cl::Kernel kernel( m_Program, "CoordinatesCenter", &error_number);

      // set the args values
      error_number  = kernel.setArg(  0, COORDINATES.GetData());
      error_number |= kernel.setArg(  1, cl_uint( COORDINATES.GetNumberRows()));
      error_number |= kernel.setArg(  2, center.GetData());
      error_number |= kernel.setArg(  3, s_BlockSize * s_BlockSize * sizeof( double), 0); // shared memory
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // Execute
      const cl::NDRange block_dims( s_BlockSize, s_BlockSize);
      const cl::NDRange offset;
      const cl::NDRange worksize( s_BlockSize, s_BlockSize);

      error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      // end
      return center;
    }

    //! @brief Transformation matrix from Covariance matrix
    //! @param MOMENT covariance matrix
    //! @param CENTER_COORDINATES of coordinates
    //! @param CENTER_REFERENCE_COORDINATES center of reference coordinates
    //! @return transformation matrix
    Matrix< double> QualityRMSD::CovarianceToTransformationMatrix
    (
      const Matrix3x3< double> &MOMENT,
      const Vector< double> &CENTER_COORDINATES,
      const Vector< double> &CENTER_REFERENCE_COORDINATES
    ) const
    {
      // diagonalization
      Matrix3x3< double> rotate_device( MOMENT.HardCopy());
      rotate_device.MultiplyWithTransposed();

      // solve Eigensystem
      Vector< double> eigenvalues( 3, m_Queue, s_BlockSize - 3, 0.0);
      Matrix3x3< double> eigenvectors( m_Queue, 0.0);
      rotate_device.EigenVectorsSymmetric( eigenvectors, eigenvalues);
      eigenvectors.Transpose();
      eigenvectors.SortRowsAndVector( eigenvalues);
      eigenvectors.Orthogonalize( 2);
//      linal::Vector< double> eigenvalues( eigenvalues.GetHostVector( s_BlockSize - 3));
//
//      // check second eigenvalue
//      if( eigenvalues( 1) <= 0.0 || eigenvalues( 0) <= 0.0)
//      {
//        return Matrix< double>( 4, 4, m_Queue);
//      } //error

      //build rotation matrix
      Matrix3x3< double> rotate( eigenvectors.HardCopy());
      rotate *= MOMENT;
      rotate.NormalizeRows( eigenvalues);
      rotate.Orthogonalize( 2);
      rotate.Transpose();
      rotate *= eigenvectors;

      // shift and rotate molecule
      Matrix< double> transformation( 4, 4, m_Queue, s_BlockSize - 4, s_BlockSize - 4);
      transformation.SetDiagonal( double( 1));

      // apply translation
      {
        cl_int error_number( CL_SUCCESS);
        cl::Kernel kernel( Operations< double>::GetInstance().GetProgram(), "TranslationOnTransformation", &error_number);

        // set the args values
        error_number  = kernel.setArg(  0, transformation.GetData());
        error_number |= kernel.setArg(  1, CENTER_REFERENCE_COORDINATES.GetData());
        error_number |= kernel.setArg(  2, double( -1));
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // Execute
        const cl::NDRange block_dims( s_BlockSize);
        const cl::NDRange offset;
        const cl::NDRange worksize( s_BlockSize);

        error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      }
      {
        cl_int error_number( CL_SUCCESS);
        cl::Kernel kernel( Operations< double>::GetInstance().GetProgram(), "RotationOnTransformation", &error_number);

        // set the args values
        error_number  = kernel.setArg(  0, transformation.GetData());
        error_number |= kernel.setArg(  1, rotate.GetData());
        error_number |= kernel.setArg(  2, s_BlockSize * s_BlockSize * sizeof( double), 0);
        error_number |= kernel.setArg(  3, s_BlockSize * s_BlockSize * sizeof( double), 0);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // Execute
        const cl::NDRange block_dims( s_BlockSize, s_BlockSize);
        const cl::NDRange offset;
        const cl::NDRange worksize( s_BlockSize, s_BlockSize);

        error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      }
      {
        cl_int error_number( CL_SUCCESS);
        cl::Kernel kernel( Operations< double>::GetInstance().GetProgram(), "TranslationOnTransformation", &error_number);

        // set the args values
        error_number  = kernel.setArg(  0, transformation.GetData());
        error_number |= kernel.setArg(  1, CENTER_COORDINATES.GetData());
        error_number |= kernel.setArg(  2, double( 1));
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

        // Execute
        const cl::NDRange block_dims( s_BlockSize);
        const cl::NDRange offset;
        const cl::NDRange worksize( s_BlockSize);

        error_number = m_Queue.enqueueNDRangeKernel( kernel, offset, worksize, block_dims);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      }

      // return transformation matrix calculated
      return transformation;
    }

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_saxs_debye.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math.h"
#include "math/bcl_math_sum_function.h"
#include "opencl/bcl_opencl_vector.h"
#include "restraint/bcl_restraint_sas_data_parameters.h"
#include "restraint/bcl_restraint_sas_debye.h"
#include "restraint/bcl_restraint_sas_experimental_and_calculated_data.h"
#include "restraint/bcl_restraint_sas_scattering_data.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_stopwatch.h"
#include "util/bcl_util_undefined.h"

// external includes - sorted alphabetically
#include <CL/cl_platform.h>

namespace bcl
{
  namespace opencl
  {

    class BCL_API FormFactorParameters
    {

    private:

    //////////
    // data //
    //////////

      //! float (c1 parameter)
      float m_ExcludedVolumeParameter;

      //! float (c2 parameter)
      float m_HydrationShellParameter;

      //! float ( sasa)
      float m_SolventAccessableSurfaceArea;

      //! float DisplacedSolventVolume
      float m_DisplacedSolventVolume;

      //! float VanderWaals Radius
      float m_Radius;

      //! float Bound Hydrogen
      float m_BoundHydrogen;

      //! float Crommer Mann Coefficient A1
      float m_A1;

      //! float Crommer Mann Coefficient A2
      float m_A2;

      //! float Crommer Mann Coefficient A3
      float m_A3;

      //! float Crommer Mann Coefficient A4
      float m_A4;

      //! float Crommer Mann Coefficient B1
      float m_B1;

      //! float Crommer Mann Coefficient B2
      float m_B2;

      //! float Crommer Mann Coefficient B3
      float m_B3;

      //! float Crommer Mann Coefficient B4
      float m_B4;

      //! float Crommer Mann Coefficient C
      float m_C;

      //! float X Coordinate
      float m_X;

      //! float Y Coordinate
      float m_Y;

      //! float Z Coordinate
      float m_Z;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @breif Default Constructor
      FormFactorParameters();

      //! @brief Clone function
      //! @return pointer to new FormFactorParameters
      FormFactorParameters *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief set ExcludedVolumeParameter (c1) member variable
      void SetExcludedVolumeParameter( const float &C1)
      {
        m_ExcludedVolumeParameter = C1;
      }

      //! @brief set HydrationShellParameter (c2) member variable
      void SetHydrationShellParameter( const float &C2)
      {
        m_HydrationShellParameter = C2;
      }

      //! @brief set SolventAccessableSurfaceArea ( sasa) member variable
      void SetSolventAccessableSurfaceArea( const double &SASA)
      {
        m_SolventAccessableSurfaceArea = ( float)SASA;
      }

      //! @brief set DisplacedSolventVolume member variable
      void SetDisplacedSolventVolume( const double &DSV)
      {
        m_DisplacedSolventVolume = ( float)DSV;
      }

      //! @brief set VanderWaals Radius member variable
      void SetRadius( const double &RADIUS)
      {
        m_Radius = ( float)RADIUS;
      }

      //! @brief set Bound Hydrogen member variable
      void SetBoundHydrogen( const size_t &BOUND_HYDROGEN)
      {
        m_BoundHydrogen = ( float)BOUND_HYDROGEN;
      }

      //! @brief set Crommer Mann Coefficient A1 member variable
      void SetA1( const double &A1)
      {
        m_A1 = ( float)A1;
      }

      //! @brief set Crommer Mann Coefficient A2 member variable
      void SetA2( const double &A2)
      {
        m_A2 = ( float)A2;
      }

      //! @brief set Crommer Mann Coefficient A3 member variable
      void SetA3( const double &A3)
      {
        m_A3 = ( float)A3;
      }

      //! @brief set Crommer Mann Coefficient A4 member variable
      void SetA4( const double &A4)
      {
        m_A4 = ( float)A4;
      }

      //! @brief set Crommer Mann Coefficient B1 member variable
      void SetB1( const double &B1)
      {
        m_B1 = ( float)B1;
      }

      //! @brief set Crommer Mann Coefficient B2 member variable
      void SetB2( const double &B2)
      {
        m_B2 = ( float)B2;
      }

      //! @brief set Crommer Mann Coefficient B3 member variable
      void SetB3( const double &B3)
      {
        m_B3 = ( float)B3;
      }

      //! @brief set Crommer Mann Coefficient B4 member variable
      void SetB4( const double &B4)
      {
        m_B4 = ( float)B4;
      }

      //! @brief set Crommer Mann Coefficient C member variable
      void SetC( const double &C)
      {
        m_C = ( float)C;
      }

      //! @brief set X Coordinate
      void SetX( const double &X)
      {
        m_X = ( float)X;
      }

       //! @brief set Y Coordinate
      void SetY( const double &Y)
      {
        m_Y = ( float)Y;
      }

      //! @brief set X Coordinate
      void SetZ( const double &Z)
      {
        m_Z = ( float)Z;
      }

      void SetParameters
      (
        const storage::Vector< std::string> &ATOMGROUP,
        const storage::Vector< linal::Vector3D> &COORDINATES,
        const storage::Vector< double> &SASA_VALUE,
        size_t LOCATION,
        float EXCLUDED_VOLUME,
        float HYDRATION_SHELL
      );

      void ShowValues();

      cl_float16 GetParametersAsFloat16()
      {
        cl_float16 host;

        host.s[0] = m_ExcludedVolumeParameter;
        host.s[1] = m_HydrationShellParameter;
        host.s[2] = m_SolventAccessableSurfaceArea;
        host.s[3] = m_DisplacedSolventVolume;
        host.s[4] = m_Radius;
        host.s[5] = m_BoundHydrogen;
        host.s[6] = m_A1;
        host.s[7] = m_A2;
        host.s[8] = m_A3;
        host.s[9] = m_A4;
        host.s[10] = m_B1;
        host.s[11] = m_B2;
        host.s[12] = m_B3;
        host.s[13] = m_B4;
        host.s[14] = m_C;
        host.s[15] = float( 0.0);

        return host;

      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class FormFactorParameters

    FormFactorParameters::FormFactorParameters() :
      m_ExcludedVolumeParameter( util::GetUndefined< float>()),
      m_HydrationShellParameter( util::GetUndefined< float>()),
      m_SolventAccessableSurfaceArea( util::GetUndefined< float>()),
      m_DisplacedSolventVolume( util::GetUndefined< float>()),
      m_Radius( util::GetUndefined< float>()),
      m_BoundHydrogen( util::GetUndefined< float>()),
      m_A1( util::GetUndefined< float>()),
      m_A2( util::GetUndefined< float>()),
      m_A3( util::GetUndefined< float>()),
      m_A4( util::GetUndefined< float>()),
      m_B1( util::GetUndefined< float>()),
      m_B2( util::GetUndefined< float>()),
      m_B3( util::GetUndefined< float>()),
      m_B4( util::GetUndefined< float>()),
      m_C( util::GetUndefined< float>()),
      m_X( util::GetUndefined< float>()),
      m_Y( util::GetUndefined< float>()),
      m_Z( util::GetUndefined< float>())
    {
    }

    void FormFactorParameters::SetParameters
    (
      const storage::Vector< std::string> &ATOMGROUP,
      const storage::Vector< linal::Vector3D> &COORDINATES,
      const storage::Vector< double> &SASA_VALUE,
      size_t LOCATION,
      float EXCLUDED_VOLUME,
      float HYDRATION_SHELL
    )
    {
      SetExcludedVolumeParameter( EXCLUDED_VOLUME);
      SetHydrationShellParameter( HYDRATION_SHELL);
      SetSolventAccessableSurfaceArea( SASA_VALUE( LOCATION));
      SetX( COORDINATES( LOCATION).X());
      SetY( COORDINATES( LOCATION).Y());
      SetZ( COORDINATES( LOCATION).Z());

      if( ATOMGROUP( LOCATION) == "H")
      {
        SetDisplacedSolventVolume( 5.15);
        SetRadius( 1.07);
        SetBoundHydrogen( 0);
        SetA1( 0.493002);
        SetA2( 0.322912);
        SetA3( 0.140191);
        SetA4( 0.040810);
        SetB1( 10.510900);
        SetB2( 26.125700);
        SetB3(  3.142360);
        SetB4( 57.799700);
        SetC( 0.003038);
      }
      else if
      (
          ATOMGROUP( LOCATION) == "C" ||
          ATOMGROUP( LOCATION) == "CH" ||
          ATOMGROUP( LOCATION) == "CH2" ||
          ATOMGROUP( LOCATION) == "CH3"
      )
      {
        SetA1( 2.310000);
        SetA2( 1.020000);
        SetA3( 1.588600);
        SetA4( 0.865000);
        SetB1( 20.843900);
        SetB2( 10.207500);
        SetB3( 0.568700);
        SetB4( 51.651200);
        SetC( 0.215600);

        if( ATOMGROUP( LOCATION) == "C")
        {
          SetDisplacedSolventVolume( 16.44);
          SetRadius( 1.58);
          SetBoundHydrogen( 0);
        }
        else if( ATOMGROUP( LOCATION) == "CH")
        {
          SetDisplacedSolventVolume( 21.59);
          SetRadius( 1.73);
          SetBoundHydrogen( 1);
        }
        else if( ATOMGROUP( LOCATION) == "CH2")
        {
          SetDisplacedSolventVolume( 26.74);
          SetRadius( 1.85);
          SetBoundHydrogen( 2);
        }
        else
        {
          SetDisplacedSolventVolume( 31.89);
          SetRadius( 1.97);
          SetBoundHydrogen( 3);
        }
      }
      else if
      (
          ATOMGROUP( LOCATION) == "N" ||
          ATOMGROUP( LOCATION) == "NH" ||
          ATOMGROUP( LOCATION) == "NH2" ||
          ATOMGROUP( LOCATION) == "NH3"
      )
      {
        SetA1( 12.212600);
        SetA2( 3.132200);
        SetA3( 2.012500);
        SetA4( 1.166300);
        SetB1( 0.005700);
        SetB2( 9.893300);
        SetB3( 28.997500);
        SetB4( 0.582600);
        SetC( -11.529000);

        if( ATOMGROUP( LOCATION) == "N")
        {
          SetDisplacedSolventVolume( 2.49);
          SetRadius( 0.84);
          SetBoundHydrogen( 0);
        }
        else if( ATOMGROUP( LOCATION) == "NH")
        {
          SetDisplacedSolventVolume( 7.64);
          SetRadius( 1.22);
          SetBoundHydrogen( 1);
        }
        else if( ATOMGROUP( LOCATION) == "NH2")
        {
          SetDisplacedSolventVolume( 12.79);
          SetRadius( 1.45);
          SetBoundHydrogen( 2);
        }
        else
        {
          SetDisplacedSolventVolume( 17.94);
          SetRadius( 1.62);
          SetBoundHydrogen( 3);
        }
      }
      else if
      (
          ATOMGROUP( LOCATION) == "O" ||
          ATOMGROUP( LOCATION) == "OH"
      )
      {
        SetA1( 3.048500);
        SetA2( 2.286800);
        SetA3( 1.546300);
        SetA4(  0.867000);
        SetB1( 13.277100);
        SetB2(  5.701100);
        SetB3( 0.323900);
        SetB4( 32.908900);
        SetC(   0.250800);

        if( ATOMGROUP( LOCATION) == "O")
        {
          SetDisplacedSolventVolume( 9.13);
          SetRadius( 1.30);
          SetBoundHydrogen( 0);
        }
        else
        {
          SetDisplacedSolventVolume( 14.28);
          SetRadius( 1.5);
          SetBoundHydrogen( 1);
        }
      }
      else if
      (
         ATOMGROUP( LOCATION) == "S" ||
         ATOMGROUP( LOCATION) == "SH"
      )
      {
        SetA1( 6.905300);
        SetA2( 5.203400);
        SetA3( 1.437900);
        SetA4( 1.586300);
        SetB1( 1.467900);
        SetB2( 22.215100);
        SetB3( 0.253600);
        SetB4( 56.172000);
        SetC( 0.866900);

        if( ATOMGROUP( LOCATION) == "S")
        {
          SetDisplacedSolventVolume( 19.86);
          SetRadius( 1.68);
          SetBoundHydrogen( 0);
        }
        else
        {
          SetDisplacedSolventVolume( 25.10);
          SetRadius( 1.81);
          SetBoundHydrogen( 1);
        }
      }
      else if
      (
         ATOMGROUP( LOCATION) == "SE"
      )
      {
        SetA1( 17.000600);
        SetA2( 5.819600);
        SetA3( 3.973100);
        SetA4( 4.354300);
        SetB1( 2.409800);
        SetB2( 0.272600);
        SetB3( 15.237200);
        SetB4( 43.816300);
        SetC( 2.840900);
        SetDisplacedSolventVolume( 28.73);
        SetRadius( 1.90);
        SetBoundHydrogen( 0);
      }
      else
      {
        std::string grouptype( ATOMGROUP( LOCATION));
        BCL_Message( util::Message::e_Standard, " Something is wrong, the Atomtype is: " + util::Format()( grouptype));
      }
    }

    void FormFactorParameters::ShowValues()
    {
      BCL_Message( util::Message::e_Standard, " C1 scaling parameter : " + util::Format()(  m_ExcludedVolumeParameter ));
      BCL_Message( util::Message::e_Standard, " C2 scaling parameter : " + util::Format()(  m_HydrationShellParameter ));
      BCL_Message( util::Message::e_Standard, " Sasa parameter : " + util::Format()(  m_SolventAccessableSurfaceArea ));
      BCL_Message( util::Message::e_Standard, " Solvent parameter : " + util::Format()(  m_DisplacedSolventVolume ));
      BCL_Message( util::Message::e_Standard, " Radius : " + util::Format()(  m_Radius ));
      BCL_Message( util::Message::e_Standard, " Bound Hydrogen :http://www.sltrib.com " + util::Format()(  m_BoundHydrogen ));
      BCL_Message( util::Message::e_Standard, " A1 : " + util::Format()(  m_A1 ));
      BCL_Message( util::Message::e_Standard, " A2 : " + util::Format()(  m_A2 ));
      BCL_Message( util::Message::e_Standard, " A3 : " + util::Format()(  m_A3 ));
      BCL_Message( util::Message::e_Standard, " A4 : " + util::Format()(  m_A4 ));
      BCL_Message( util::Message::e_Standard, " B1 : " + util::Format()(  m_B1 ));
      BCL_Message( util::Message::e_Standard, " B2 : " + util::Format()(  m_B2 ));
      BCL_Message( util::Message::e_Standard, " B3 : " + util::Format()(  m_B3 ));
      BCL_Message( util::Message::e_Standard, " B4 : " + util::Format()(  m_B4 ));
      BCL_Message( util::Message::e_Standard, " C : " + util::Format()(  m_C ));
      BCL_Message( util::Message::e_Standard, " X : " + util::Format()(  m_X ));
      BCL_Message( util::Message::e_Standard, " Y : " + util::Format()(  m_Y ));
      BCL_Message( util::Message::e_Standard, " Z : " + util::Format()(  m_Z ));
    }

  //////////
  // data //
  //////////

    //! single instance of that class
    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SaxsDebye::s_Instance
    (
      util::Enumerated< restraint::SasDebyeInterface>::AddInstance( new SaxsDebye())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SaxsDebye::SaxsDebye() :
      m_ShouldApproximateLoops( false),
      m_DetermineAnalyticNormFactor( false),
      m_ExcludedVolumeParameter( 1.0),
      m_HydrationShellParameter( 0.0),
      m_ShouldApproximateSideChains( true),
      m_ReducedExpData( util::ShPtr< storage::Vector< restraint::SasScatteringPoint> >()),
      m_Queue(),
      m_Program()
    {
    }

    //! @brief Constructor that takes a bool
    //! @param LOOPS bool value to represent loops that are not present in the protein model
    SaxsDebye::SaxsDebye
    (
      const CommandQueue &QUEUE,
      const bool LOOPS,
      const bool USE_REGULA_FALSI_APPROXIMATION,
      const float EXCLUDED_VOLUME_PARAMETER,
      const float HYDRATION_SHELL_PARAMETER,
      const bool SIDE_CHAIN_APPROXIMATION,
      const util::ShPtr< storage::Vector< restraint::SasScatteringPoint> > REDUCED_EXPERIMENTAL_DATA
    ) :
      m_ShouldApproximateLoops( LOOPS),
      m_DetermineAnalyticNormFactor( USE_REGULA_FALSI_APPROXIMATION),
      m_ExcludedVolumeParameter( EXCLUDED_VOLUME_PARAMETER),
      m_HydrationShellParameter( HYDRATION_SHELL_PARAMETER),
      m_ShouldApproximateSideChains( SIDE_CHAIN_APPROXIMATION),
      m_ReducedExpData( REDUCED_EXPERIMENTAL_DATA),
      m_Queue( QUEUE)
    {
      cl_int error_number( CL_SUCCESS);
      m_Program = KernelSources::Compile( GetKernelSources().e_Saxs, util::CPPDataTypes::e_Float, m_Queue, std::string(), &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
    }

    //! @brief Clone function
    //! @return pointer to new SaxsDebye
    SaxsDebye *SaxsDebye::Clone() const
    {
      return new SaxsDebye( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SaxsDebye::GetAlias() const
    {
      static const std::string s_Name( "OpenCLSaxsDebye");
      return s_Name;
    }

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SaxsDebye::GetClassIdentifier() const
    {
      // Get BCL Standardized Class name
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes a proteinModel as input and calculates an intensity using the debye formula
    //! @param PROTEIN_MODEL
    //! @return the intensity for a given q value
    restraint::SasExperimentalAndCalculatedData SaxsDebye::operator()
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      util::ShPtr< restraint::SasScatteringData> sp_experimental_data;

      if( !this->GetReducedExperimentalData().IsDefined())
      {

        // get the experimental SAXS data
        sp_experimental_data = this->GetExperimentalData();

        // Verify you have experimental Saxs Data
        if( !sp_experimental_data.IsDefined())
        {
          // warn user and return empty data
          BCL_Message( util::Message::e_Critical, "No experimental SAXS data found, returning empty data");
          return restraint::SasExperimentalAndCalculatedData();
        }

      }
      else
      {

        restraint::SasScatteringData experimental_data;

        for
        (
          storage::Vector< restraint::SasScatteringPoint>::const_iterator
           exp_data_itr( this->GetReducedExperimentalData()->Begin()),
           exp_data_itr_end( this->GetReducedExperimentalData()->End());
          exp_data_itr != exp_data_itr_end;
          ++exp_data_itr
        )
        {
          experimental_data.PushBackScattering( *exp_data_itr);
        }

        // use the Clone to ShPtr to create a shared pointer to experimental data
        util::ShPtr< restraint::SasScatteringData> sp_reduced_data( util::CloneToShPtr( experimental_data));

        // set sp_experimental_data to the reduced data set
        sp_experimental_data = sp_reduced_data;
      }

      restraint::SasDebye saxs_debye_object
      (
        m_ShouldApproximateLoops,
        m_DetermineAnalyticNormFactor,
        m_ExcludedVolumeParameter,
        m_HydrationShellParameter,
        m_ShouldApproximateSideChains
      );

      saxs_debye_object.GetAtomsAndFormFactors( PROTEIN_MODEL);

       // create object to hold calculated data
      restraint::SasScatteringData calculated_data;

      calculated_data.AllocateScatteringMemory( sp_experimental_data->GetScatteringData().GetSize());

      const size_t number_of_atoms( saxs_debye_object.GetCoordinates().GetSize());

      // host_coords is pointer to a vector of 4 floats
      cl_float4 *host_coords;
      cl_float16 *host_params;

      // calloc allocates a block of memory for an array of num elements, each of them size bytes long
      // and initializes all its bits to zero

      // point host_coords to a block of memory the size of number of atoms by the size of cl_float4
      host_coords = ( cl_float4 *)calloc( number_of_atoms, sizeof( cl_float4));
      host_params = ( cl_float16 *)calloc( number_of_atoms, sizeof( cl_float16));

      FormFactorParameters data_set;

      // get the coordinates of each atom
      for( size_t row( 0), row_end( number_of_atoms); row < row_end; ++row)
      {
        host_coords[ row].s[0] = float( saxs_debye_object.GetCoordinates()( row)( 0));
        host_coords[ row].s[1] = float( saxs_debye_object.GetCoordinates()( row)( 1));
        host_coords[ row].s[2] = float( saxs_debye_object.GetCoordinates()( row)( 2));
        host_coords[ row].s[3] = float( 0);

        data_set.SetParameters
        (
          saxs_debye_object.GetAtomGroups(),
          saxs_debye_object.GetCoordinates(),
          saxs_debye_object.GetSASAPoint(),
          row,
          m_ExcludedVolumeParameter,
          m_HydrationShellParameter
        );

        host_params[ row] = data_set.GetParametersAsFloat16();
      }

      // initialize error number with success
      cl_int error_number = CL_SUCCESS;

      // bcl buffer - wrapper around opencl bindings
      // allocating memory on GPU, just an allocation
      Buffer device_coords( cl::Buffer( m_Queue.GetContext(), CL_FALSE, sizeof( cl_float4) * number_of_atoms, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "buffer allocation/copy failed with: " + opencl::Tools::ErrorString( error_number));

      // transfers host allocation of memory to the GPU
      error_number = m_Queue.enqueueWriteBuffer( device_coords, CL_FALSE, 0, sizeof( cl_float4) * number_of_atoms, host_coords);
      BCL_Assert( error_number == CL_SUCCESS, "buffer allocation/copy failed with: " + opencl::Tools::ErrorString( error_number));

      // allocating for parameters
      Buffer device_params( cl::Buffer( m_Queue.GetContext(), CL_FALSE, sizeof( cl_float16) * number_of_atoms, NULL, &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "buffer allocation/copy failed with: " + opencl::Tools::ErrorString( error_number));

      // transfers host allocation of memory to the GPU
      error_number = m_Queue.enqueueWriteBuffer( device_params, CL_FALSE, 0, sizeof( cl_float16) * number_of_atoms, host_params);
      BCL_Assert( error_number == CL_SUCCESS, "buffer allocation/copy failed with: " + opencl::Tools::ErrorString( error_number));

      // get the number of q values
      const size_t number_q_values( sp_experimental_data->GetScatteringData().GetSize());

      // set up a vector of floats and doubles for the number of q values
      linal::Vector< float> q_values( number_q_values);
      linal::Vector< double> q_doubles( number_q_values);

      // reset row to 0
      size_t row( 0);

      // precalculate water factors and H for each q value
      // get the structure factor function for hydrogen

      // get the crommer mann constants for hydrogen
      static const math::SumFunction< restraint::SasDataParameters, double>
      h_form_factor( biol::GetAtomTypes().H->GetElementType()->GetStructureFactor(), 1.0, 0.0);

      // get the crommer mann constants for oxygen
      static const math::SumFunction< restraint::SasDataParameters, double>
      o_form_factor( biol::GetAtomTypes().O->GetElementType()->GetStructureFactor(), 1.0, 0.0);

      // initialize structure factors for water
      math::SumFunction< restraint::SasDataParameters, double> water_factors;
      math::SumFunction< restraint::SasDataParameters, double> h_factors;

      water_factors += double( 2.0) * h_form_factor;
      water_factors += o_form_factor;

      h_factors += h_form_factor;

      storage::Vector< float> water_factors_vector( number_q_values);

      storage::Vector< float> h_factors_vector( number_q_values);
      // iterate over experimental data to get q-values
      for
      (
        storage::Vector< restraint::SasScatteringPoint>::const_iterator
          data_itr( sp_experimental_data->GetScatteringData().Begin()),
          data_itr_end( sp_experimental_data->GetScatteringData().End());
        data_itr != data_itr_end;
        ++data_itr, ++row
      )
      {
        // variable to hold q-value both double and float forms
        q_values( row) = float( data_itr->GetQvalue());
        q_doubles( row) = data_itr->GetQvalue();

        restraint::SasDataParameters q_value( q_doubles( row));

        water_factors_vector( row) = float( water_factors( q_value));
        h_factors_vector( row) = float( h_factors( q_value));
      }

      //BCL_MessageDbg( " Water factor vector: " + util::Format()( water_factors_vector));
      //BCL_MessageDbg( " hydrogen factor vector: " + util::Format()( h_factors_vector));

      Vector< float> res_ff( number_of_atoms, m_Queue);

      // initialize
      Vector< float> inner_sum_matrix( number_of_atoms, m_Queue);

      // allocate memory to be transfered back
      linal::Vector< float> calculated_intensities( number_q_values);

      cl::Kernel inner_sum_kernel, inner_sum_kernel_zero;

      cl::Kernel res_ff_kernel( m_Program, "CalculateResFF", &error_number);

      BCL_Assert( error_number == CL_SUCCESS, "CalculateResFF arg error: " + opencl::Tools::ErrorString( error_number));
      if( GetTools().GetFirstCommandQueue().GetDevice( NULL).DeviceType( NULL) == CL_DEVICE_TYPE_GPU)
      {

        BCL_MessageDbg( " GPU Platform");
        // relies of implicit synchronization because of the warps
        inner_sum_kernel = cl::Kernel( m_Program, "InnerSumsOTFGPU", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "inner sum kernel arg error: " + opencl::Tools::ErrorString( error_number));

        inner_sum_kernel_zero = cl::Kernel( m_Program, "InnerSumsOTFZeroGPU", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "inner sum kernel arg error: " + opencl::Tools::ErrorString( error_number));
      }
      else
      {
        BCL_MessageDbg( " CPU Platform");
        // does not rely on implicit synchronization
        inner_sum_kernel = cl::Kernel( m_Program, "InnerSumsOTFCPU", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "inner sum kernel arg error: " + opencl::Tools::ErrorString( error_number));

        inner_sum_kernel_zero = cl::Kernel( m_Program, "InnerSumsOTFZeroCPU", &error_number);
        BCL_Assert( error_number == CL_SUCCESS, "inner sum kernel arg error: " + opencl::Tools::ErrorString( error_number));
      }

      // sum of intensities
      cl::Kernel intensity_sums_kernel( m_Program, "ReductionSumOfIntensity", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "inner sum kernel arg error: " + opencl::Tools::ErrorString( error_number));

      // rmsd group and worksize
      const cl_uint block_size = 128;
      const cl::NDRange block_dim( block_size);
      const cl::NDRange nr_atoms_worksize( Tools::RoundUp( block_size, number_of_atoms));
      const cl::NDRange offset;

      const size_t number_elements_partial_reduction( ( number_of_atoms / block_size) + 1);

      size_t q_index( 0);

      Vector< float> partial_reduction_output( number_elements_partial_reduction, m_Queue);

      if( q_values( q_index) == 0)
      {
        //BCL_MessageDbg( " q_values is zero ");

        error_number  = res_ff_kernel.setArg( 0, device_params);
        error_number  = res_ff_kernel.setArg( 1, q_values( q_index));
        error_number  = res_ff_kernel.setArg( 2, res_ff.GetData());
        error_number  = res_ff_kernel.setArg( 3, cl_uint( number_of_atoms));
        error_number  = res_ff_kernel.setArg( 4, h_factors_vector( q_index));
        error_number  = res_ff_kernel.setArg( 5, water_factors_vector( q_index));
        BCL_Assert( error_number == CL_SUCCESS, "ffs arg error: " + opencl::Tools::ErrorString( error_number));

        error_number = m_Queue.enqueueNDRangeKernel( res_ff_kernel, offset, nr_atoms_worksize, block_dim, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

        //BCL_MessageDbg( " gpu ffs: " + util::Format()( linal::Matrix< float>( res_ff.GetSize(), 1, res_ff.GetHostVector().Begin())));

        error_number  = inner_sum_kernel_zero.setArg( 0, device_coords);
        error_number |= inner_sum_kernel_zero.setArg( 1, res_ff.GetData());
        error_number |= inner_sum_kernel_zero.setArg( 2, q_values( q_index));
        error_number |= inner_sum_kernel_zero.setArg( 3, inner_sum_matrix.GetData());
        error_number |= inner_sum_kernel_zero.setArg( 4, cl_uint( number_of_atoms));
        BCL_Assert( error_number == CL_SUCCESS, "inner kernel zero arg error: " + opencl::Tools::ErrorString( error_number));

        // set
        error_number  = intensity_sums_kernel.setArg( 0, inner_sum_matrix.GetData());
        error_number |= intensity_sums_kernel.setArg( 1, cl_uint( number_of_atoms));
        error_number |= intensity_sums_kernel.setArg( 2, partial_reduction_output.GetData());
        error_number |= intensity_sums_kernel.setArg( 3, block_size * sizeof( float), 0);
        BCL_Assert( error_number == CL_SUCCESS, "inner sum kernel arg error: " + opencl::Tools::ErrorString( error_number));

        // launch kernel
        error_number = m_Queue.enqueueNDRangeKernel( inner_sum_kernel_zero, offset, nr_atoms_worksize, block_dim, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

        error_number = m_Queue.enqueueNDRangeKernel( intensity_sums_kernel, offset, nr_atoms_worksize, block_dim, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

        calculated_data.PushBackScattering( restraint::SasScatteringPoint( q_doubles( q_index), double( partial_reduction_output.GetHostVector().Sum()), 0.0));

        ++q_index;
      }

      for( ; q_index < number_q_values; ++q_index)
      {
        //BCL_MessageDbg( " q_values is not zero ");
        error_number  = res_ff_kernel.setArg( 0, device_params);
        error_number  = res_ff_kernel.setArg( 1, q_values( q_index));
        error_number  = res_ff_kernel.setArg( 2, res_ff.GetData());
        error_number  = res_ff_kernel.setArg( 3, cl_uint( number_of_atoms));
        error_number  = res_ff_kernel.setArg( 4, h_factors_vector( q_index));
        error_number  = res_ff_kernel.setArg( 5, water_factors_vector( q_index));
        BCL_Assert( error_number == CL_SUCCESS, "ffs arg error: " + opencl::Tools::ErrorString( error_number));

        error_number = m_Queue.enqueueNDRangeKernel( res_ff_kernel, offset, nr_atoms_worksize, block_dim, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

        //BCL_MessageDbg( " gpu ffs: " + util::Format()( linal::Matrix< float>( res_ff.GetSize(), 1, res_ff.GetHostVector().Begin())));

        // set
        error_number  = inner_sum_kernel.setArg( 0, device_coords);
        error_number |= inner_sum_kernel.setArg( 1, res_ff.GetData());
        error_number |= inner_sum_kernel.setArg( 2, q_values( q_index));
        error_number |= inner_sum_kernel.setArg( 3, inner_sum_matrix.GetData());
        error_number |= inner_sum_kernel.setArg( 4, cl_uint( number_of_atoms));
        BCL_Assert( error_number == CL_SUCCESS, "inner kernel arg error: " + opencl::Tools::ErrorString( error_number));

        // set
        error_number  = intensity_sums_kernel.setArg( 0, inner_sum_matrix.GetData());
        error_number |= intensity_sums_kernel.setArg( 1, cl_uint( number_of_atoms));
        error_number |= intensity_sums_kernel.setArg( 2, partial_reduction_output.GetData());
        error_number |= intensity_sums_kernel.setArg( 3, block_size * sizeof( float), 0);
        BCL_Assert( error_number == CL_SUCCESS, "inner sum kernel arg error: " + opencl::Tools::ErrorString( error_number));

        // launch kernel
        error_number = m_Queue.enqueueNDRangeKernel( inner_sum_kernel, offset, nr_atoms_worksize, block_dim, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));

        error_number = m_Queue.enqueueNDRangeKernel( intensity_sums_kernel, offset, nr_atoms_worksize, block_dim, NULL, NULL);
        BCL_Assert( error_number == CL_SUCCESS, "error: " + opencl::Tools::ErrorString( error_number));
        calculated_data.PushBackScattering( restraint::SasScatteringPoint( q_doubles( q_index), double( partial_reduction_output.GetHostVector().Sum()), 0.0));
      }

      restraint::SasExperimentalAndCalculatedData saxs_data( *sp_experimental_data, calculated_data);

      // Free the memory used by calloc
      delete [] host_coords;
      delete [] host_params;

      return saxs_data;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief Set the members of this property from the given LABEL
    //! @param LABEL the label to parse
    //! @param ERR_STREAM stream to write errors out to
    bool SaxsDebye::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERR_STREAM)
    {
      // only return true if the opencl kernels are initialized; this might be false if -opencl Disable was passed
      return GetTools().HasCommandQueues();
    }

    //! @brief responsible for updating to a valid queue
    //! @param TOOLS opencl tools
    void SaxsDebye::UpdateQueue( Tools &TOOLS)
    {
      if( !TOOLS.HasCommandQueues())
      {
        return;
      }

      m_Queue = TOOLS.GetFirstCommandQueue();

      cl_int error_number( CL_SUCCESS);
      m_Program = KernelSources::Compile( GetKernelSources().e_Saxs, util::CPPDataTypes::e_Float, m_Queue, std::string(), &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SaxsDebye::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "performs saxs debye calculation using opencl for massively parallel architectures"
      );

      parameters.AddInitializer
      (
        "consider loops",
        "should loops be considered",
        io::Serialization::GetAgent( &m_ShouldApproximateLoops),
        "0"
      );
      parameters.AddInitializer
      (
        "analytic",
        "whether to determine the norm factor with regula falsi (1) or pythagorean approximation (0)",
        io::Serialization::GetAgent( &m_DetermineAnalyticNormFactor),
        "0"
      );
      parameters.AddInitializer
      (
        "excluded volume",
        "c1 adjustment parameter",
        io::Serialization::GetAgent( &m_ExcludedVolumeParameter),
        "1.0"
      );
      parameters.AddInitializer
      (
        "hydration shell",
        "c2 adjustment parameter",
        io::Serialization::GetAgent( &m_HydrationShellParameter),
        "0.0"
      );
      parameters.AddOptionalInitializer
      (
        "approximate_sidechains",
        "sum up form factor contribution on cb position",
        io::Serialization::GetAgent( &m_ShouldApproximateSideChains)
      );

      return parameters;
    }
  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_support_vector_machine.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "opencl/bcl_opencl_kernel_sources.h"
#include "util/bcl_util_stopwatch.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SupportVectorMachine::s_Instance
    (
      GetObjectInstances().AddInstance( new SupportVectorMachine())
    );

    SupportVectorMachine::SupportVectorMachine() :
      m_Bias(),
      m_NumberSupportVectors(),
      m_Kernel(),
      m_RescaleInput(),
      m_RescaleOutput(),
      m_Queue()
    {
    }

    SupportVectorMachine::SupportVectorMachine( const CommandQueue &QUEUE) :
      m_Bias(),
      m_NumberSupportVectors(),
      m_Kernel(),
      m_RescaleInput(),
      m_RescaleOutput(),
      m_Queue( QUEUE)
    {
    }

    //! @brief default constructor
    //! @param BIAS the shift of the hyperplane from the origin
    //! @param ALPHAS vector of the alpha values
    //! @param SUPPORT_VECTORS the support vectors from the training
    //! @param KERNEL the kernel function to use
    //! @param RESCALE_IN, RESCALE_OUT the rescale functions
    SupportVectorMachine::SupportVectorMachine
    (
      const float BIAS,
      const linal::VectorConstInterface< float> &ALPHAS,
      const model::FeatureDataSetInterface< float> &SUPPORT_VECTORS,
      const util::Implementation< model::SupportVectorKernelBase> &KERNEL,
      const model::RescaleFeatureDataSet &RESCALE_IN,
      const model::RescaleFeatureDataSet &RESCALE_OUT,
      const CommandQueue &QUEUE
    ) :
      m_Bias( BIAS),
      m_NumberSupportVectors( SUPPORT_VECTORS.GetNumberFeatures()),
      m_Kernel( KERNEL),
      m_RescaleInput( RESCALE_IN),
      m_RescaleOutput( RESCALE_OUT),
      m_Queue( QUEUE)
    {
      BCL_Assert( ALPHAS.GetSize() == SUPPORT_VECTORS.GetNumberFeatures(), "Alphas vector and number support vectors should match!");
      const size_t row_padding( 0); //( s_GroupSize - ( ALPHAS.GetSize() % s_GroupSize)) % s_GroupSize);
      linal::Vector< float> padded_alphas( ALPHAS, row_padding);

      m_Alphas = Vector< float>( padded_alphas, m_Queue);
      linal::Matrix< float> padded_svs( SUPPORT_VECTORS.GetMatrix(), row_padding, 0);
      m_SupportVectors = Matrix< float>( padded_svs, m_Queue);
      // for precision type
      const cl_int error_number( CompilePrograms( util::CPPDataTypes::DataTypeFromTemplate< float>()));
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));
    }

    //! @brief Clone function
    //! @return pointer to new SupportVectorMachine
    SupportVectorMachine *SupportVectorMachine::Clone() const
    {
      return new SupportVectorMachine( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SupportVectorMachine::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SupportVectorMachine::GetAlias() const
    {
      static const std::string s_Name( "OpenCLSupportVectorMachine");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief predict result with model using a NOT rescaled feature vector
    //! @param FEATURE not rescaled feature vector
    //! @return predcited result vector using a model
    model::FeatureDataSet< float> SupportVectorMachine::PredictWithoutRescaling( const model::FeatureDataSetInterface< float> &FEATURE) const
    {
      // number of elements in one feature
      const size_t feature_length( FEATURE.GetFeatureSize());
      const size_t number_features( FEATURE.GetNumberFeatures());
      const size_t row_padding( 0); //( s_GroupSize - ( m_NumberSupportVectors % s_GroupSize)) % s_GroupSize);
      const size_t nr_groups( ( m_NumberSupportVectors % s_GroupSize == 0 ? 0 : 1) + ( m_NumberSupportVectors / s_GroupSize)); // ( m_NumberSupportVectors + row_padding) / s_GroupSize + 1);

      cl_int error_number( CL_SUCCESS);

      linal::Matrix< float> result_matrix( number_features, 1, 0.0);

      // buffer for result
      Vector< float> result_buffer( m_NumberSupportVectors + row_padding, m_Queue);
      Vector< float> reduction_result_buffer( nr_groups, m_Queue);
      linal::Vector< float> result( nr_groups, 0.0);

      // launch kernel
      cl::Kernel svm_kernel( m_Program, "ClassifyFeatureSupportVectors", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));
      cl::Kernel reduction_kernel( m_Program, "ReductionSum", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number  = svm_kernel.setArg(  1, m_SupportVectors.GetData());
      error_number |= svm_kernel.setArg(  2, m_Alphas.GetData());
      error_number |= svm_kernel.setArg(  3, cl_uint( feature_length));
      error_number |= svm_kernel.setArg(  4, cl_uint( m_NumberSupportVectors));
      error_number |= svm_kernel.setArg(  5, result_buffer.GetData());
      BCL_Assert( error_number == CL_SUCCESS, "setting gpu svm_kernel args error: " + Tools::ErrorString( error_number));

      error_number |= reduction_kernel.setArg(  1, cl_uint( m_NumberSupportVectors + row_padding));
      error_number |= reduction_kernel.setArg(  2, reduction_result_buffer.GetData());
      error_number |= reduction_kernel.setArg(  3, cl_uint( s_GroupSize * sizeof( float)), 0);
      BCL_Assert( error_number == CL_SUCCESS, "setting gpu reduction_kernel args error: " + Tools::ErrorString( error_number));

      const cl::NDRange local_worksize( s_GroupSize);
      const cl::NDRange offset;
      const cl::NDRange global_worksize( Tools::RoundUp( s_GroupSize, m_NumberSupportVectors + row_padding));
      util::Stopwatch timer;
      timer.Reset();
      timer.Start();
      // iterate over given features
      for( size_t i( 0); i < number_features; ++i)
      {
        // create buffer for feature vector
        const model::FeatureReference< float> &current_feature( FEATURE( i));
        Vector< float> device_feature( current_feature, m_Queue);

        error_number = svm_kernel.setArg( 0, device_feature.GetData());
        BCL_Assert( error_number == CL_SUCCESS, "setting gpu svm_kernel args error: " + Tools::ErrorString( error_number));

        error_number = m_Queue.enqueueNDRangeKernel( svm_kernel, offset, global_worksize, local_worksize);
        BCL_Assert( error_number == CL_SUCCESS, "svm_kernel enqueue error: " + Tools::ErrorString( error_number));

        error_number = reduction_kernel.setArg( 0, result_buffer.GetData());
        BCL_Assert( error_number == CL_SUCCESS, "setting gpu reduction_kernel args error: " + Tools::ErrorString( error_number));

        error_number = m_Queue.enqueueNDRangeKernel( reduction_kernel, offset, global_worksize, local_worksize);
        BCL_Assert( error_number == CL_SUCCESS, "reduction_kernel enqueue error: " + Tools::ErrorString( error_number));

        // read result
        result = reduction_result_buffer.GetHostVector();

        // reduce result
        result_matrix( i, 0) = result.Sum() - m_Bias;
      }
      timer.Stop();

      // return feature data set
      return model::FeatureDataSet< float>( result_matrix, m_RescaleOutput);
    }

    //! @brief Set the scaling of a feature set according to the model
    //! @param FEATURES feature set of interest
    //! @note this allows for external classes that own a dataset to ensure that a new dataset is never created
    //!       when operator() is called
    void SupportVectorMachine::Rescale( model::FeatureDataSet< float> &FEATURE) const
    {
      if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != m_RescaleInput)
      {
        FEATURE.DeScale();
        FEATURE.Rescale( m_RescaleInput);
      }
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief predict result with model using a rescaled feature vector
    //! @param FEATURE normalized or rescaled feature vector
    //! @return predcited result vector using a model
    model::FeatureDataSet< float> SupportVectorMachine::operator()( const model::FeatureDataSetInterface< float> &FEATURE) const
    {
      // handle the case where rescaling is necessary
      if( !FEATURE.IsRescaled() || *FEATURE.GetScaling() != m_RescaleInput)
      {
        model::FeatureDataSet< float> feature( FEATURE);
        feature.Rescale( m_RescaleInput);
        return PredictWithoutRescaling( feature).DeScale();
      }

      // data is already rescaled
      return PredictWithoutRescaling( FEATURE).DeScale();
    }

    //! @brief predict based on opencl::Matrix data structure
    //! @param FEATURE the matrix of features
    //! @return the matrix of predicted outputs
    Matrix< float> SupportVectorMachine::operator()( const Matrix< float> &FEATURE) const
    {
      // number of elements in one feature
      const size_t feature_length( FEATURE.GetNumberCols());
      const size_t number_features( FEATURE.GetNumberRows());
      const size_t row_padding( 0); //( s_GroupSize - ( m_NumberSupportVectors % s_GroupSize)) % s_GroupSize);
      const size_t nr_groups( ( m_NumberSupportVectors % s_GroupSize == 0 ? 0 : 1) + ( m_NumberSupportVectors / s_GroupSize)); // ( m_NumberSupportVectors + row_padding) / s_GroupSize + 1);

      cl_int error_number( CL_SUCCESS);

      linal::Matrix< float> result_matrix( number_features, 1, 0.0);

      // buffer for result
      Vector< float> result_buffer( m_NumberSupportVectors + row_padding, m_Queue);
      Vector< float> reduction_result_buffer( nr_groups, m_Queue);
      linal::Vector< float> result( nr_groups, 0.0);

      // launch kernel
      cl::Kernel svm_kernel( m_Program, "ClassifyFeatureSupportVectors", &error_number);
      cl::Kernel reduction_kernel( m_Program, "ReductionSum", &error_number);
      BCL_Assert( error_number == CL_SUCCESS, "error: " + Tools::ErrorString( error_number));

      error_number  = svm_kernel.setArg(  1, m_SupportVectors.GetData());
      error_number |= svm_kernel.setArg(  2, m_Alphas.GetData());
      error_number |= svm_kernel.setArg(  3, cl_uint( feature_length));
      error_number |= svm_kernel.setArg(  4, cl_uint( m_NumberSupportVectors));
      error_number |= svm_kernel.setArg(  5, result_buffer.GetData());
      BCL_Assert( error_number == CL_SUCCESS, "setting gpu svm_kernel args error: " + Tools::ErrorString( error_number));

      error_number |= reduction_kernel.setArg(  1, cl_uint( m_NumberSupportVectors + row_padding));
      error_number |= reduction_kernel.setArg(  2, reduction_result_buffer.GetData());
      error_number |= reduction_kernel.setArg(  3, s_GroupSize * sizeof( float), 0);
      BCL_Assert( error_number == CL_SUCCESS, "setting gpu reduction_kernel args error: " + Tools::ErrorString( error_number));

      const cl::NDRange local_worksize( s_GroupSize);
      const cl::NDRange offset;
      const cl::NDRange global_worksize( Tools::RoundUp( s_GroupSize, m_NumberSupportVectors + row_padding));
      util::Stopwatch timer;
      timer.Reset();
      timer.Start();

      // iterate over given features
      Buffer device_feature = Buffer::AllocateBufferOfSize< float>( feature_length, m_Queue);
      for( size_t i( 0); i < number_features; ++i)
      {
        // create buffer for feature vector
        m_Queue.enqueueCopyBuffer( FEATURE.GetData(), device_feature, i * feature_length, 0, sizeof( float) * feature_length);

        error_number = svm_kernel.setArg( 0, device_feature);
        BCL_Assert( error_number == CL_SUCCESS, "setting gpu svm_kernel args error: " + Tools::ErrorString( error_number));

        error_number = m_Queue.enqueueNDRangeKernel( svm_kernel, offset, global_worksize, local_worksize);
        BCL_Assert( error_number == CL_SUCCESS, "svm_kernel enqueue error: " + Tools::ErrorString( error_number));

        error_number = reduction_kernel.setArg( 0, result_buffer.GetData());
        BCL_Assert( error_number == CL_SUCCESS, "setting gpu reduction_kernel args error: " + Tools::ErrorString( error_number));

        error_number = m_Queue.enqueueNDRangeKernel( reduction_kernel, offset, global_worksize, local_worksize);
        BCL_Assert( error_number == CL_SUCCESS, "reduction_kernel enqueue error: " + Tools::ErrorString( error_number));

        // read result
        result = reduction_result_buffer.GetHostVector();

        // reduce result
        result_matrix( i, 0) = result.Sum() - m_Bias;
      }
      timer.Stop();

      return Matrix< float>( result_matrix, m_Queue);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SupportVectorMachine::Read( std::istream &ISTREAM)
    {

      io::Serialize::Read( m_AlphasHost          , ISTREAM);
      io::Serialize::Read( m_Bias                , ISTREAM);
      io::Serialize::Read( m_Kernel              , ISTREAM);
      io::Serialize::Read( m_NumberSupportVectors, ISTREAM);
      io::Serialize::Read( m_RescaleInput        , ISTREAM);
      io::Serialize::Read( m_RescaleOutput       , ISTREAM);
      io::Serialize::Read( m_SupportVectorsHost  , ISTREAM);

      m_Queue = GetTools().GetFirstCommandQueue();

      m_Alphas = Vector< float>( m_AlphasHost, m_Queue);
      m_SupportVectors = Matrix< float>( m_SupportVectorsHost.GetMatrix(), m_Queue);

      // for precision type
      cl_int error_number = CompilePrograms( util::CPPDataTypes::DataTypeFromTemplate< float>());
      BCL_Assert( error_number == CL_SUCCESS, "error compiling programs:\n" + Tools::ErrorString( error_number));

      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SupportVectorMachine::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_AlphasHost          , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Bias                , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Kernel              , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_NumberSupportVectors, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RescaleInput        , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_RescaleOutput       , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SupportVectorsHost  , OSTREAM, INDENT);

      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SupportVectorMachine::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "see http://en.wikipedia.org/wiki/Support_vector_machine"
      );

      parameters.AddInitializer
      (
        "kernel",
        "kernel used to map pairs of features onto a hyperplane",
        io::Serialization::GetAgent( &m_Kernel),
        "RBF(gamma=0.5)"
      );

      return parameters;
    }

    //! @brief compile programs for given precision
    //! @param PRECISION float or double
    //! @return ERROR error that occured, CL_SUCCESS if no error
    cl_int SupportVectorMachine::CompilePrograms( const util::CPPDataTypes::Types &PRECISION)
    {
      cl_int error_number( CL_SUCCESS);

      // compile the program
      cl::Program::Sources source;

      const Device device( m_Queue.GetDevice( &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "cannot get device for command queue");

      std::string svm_source( ( *GetKernelSources().e_SvmKernelFunctions)->GetSource( PRECISION, device.Extensions()));
      std::string reduction_source( ( *GetKernelSources().e_Linal)->GetSource( PRECISION, device.Extensions()));
      if( svm_source.empty() || reduction_source.empty())
      {
        return CL_INVALID_KERNEL_DEFINITION;
      }
      static const std::string s_actual_kernel_string( "PLACE_HOLDER_FOR_ACTUAL_KERNEL_CALL");
      svm_source.replace( svm_source.find( s_actual_kernel_string), s_actual_kernel_string.size(), m_Kernel->GetDeviceKernelFunctionCallString());
      source.push_back( std::make_pair( svm_source.c_str(), svm_source.length()));
      source.push_back( std::make_pair( reduction_source.c_str(), reduction_source.length()));

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

  } // namespace opencl
} // namespace bcl
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
#include "opencl/bcl_opencl_tools.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_operations.h"
#include "opencl/bcl_opencl_kernel_source_file.h"
#include "opencl/bcl_opencl_kernel_sources.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Tools::Tools()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Tools::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief access to the command line platform
    //! @return Platform reference to the platform selected on the command line
    const Platform &Tools::GetPlatform()
    {
      return m_Platform;
    }

    //! @brief access the default context with the commandline platform and devices
    //! @return Context reference to the context created from the command line platform and devices
    const Context &Tools::GetContext()
    {
      return m_Context;
    }

    //! @brief access to the command line devices
    //! @return storage::Vector< Device> vector of devices selected on the command line
    const storage::Vector< Device> &Tools::GetDevices()
    {
      return m_Devices;
    }

    //! @brief access to the command line queues
    //! @return storage::Vector< CommandQueue> vector of command queues selected on the command line
    const storage::Vector< CommandQueue> &Tools::GetCommandQueues()
    {
      return m_CommandQueues;
    }

    //! @brief ceck if any opencl commandqueue is available
    //! @return true if the is at least one command queue available, false if not e.g. not platform or no device on the platform was found
    bool Tools::HasCommandQueues()
    {
      return !m_CommandQueues.IsEmpty();
    }

    //! @brief convenience function to access the first command queue
    //! @return CommandQueue the first command queue
    const CommandQueue &Tools::GetFirstCommandQueue()
    {
      return m_CommandQueues.FirstElement();
    }

    //! @brief return a buffer program for that datatype
    //! @param DATA_TYPE cpp data type
    //! @return program associated with first command queues device
    const cl::Program &Tools::GetBufferProgram( const util::CPPDataTypes::Types &DATA_TYPE, const CommandQueue &QUEUE)
    {
      std::pair< util::CPPDataTypes::Types, util::SiPtr< const CommandQueue> >
        search_pair( DATA_TYPE, util::SiPtr< const CommandQueue>( QUEUE));
      std::map
      <
        std::pair< util::CPPDataTypes::Types, util::SiPtr< const CommandQueue> >,
        cl::Program
      >::const_iterator itr( m_BufferPrograms.find( search_pair));

      if( itr != m_BufferPrograms.end())
      {
        return itr->second;
      }

      // try to compile program
      cl_int error_number( CL_SUCCESS);
      itr =
        m_BufferPrograms.insert
        (
          std::make_pair
          (
            search_pair,
            KernelSources::Compile( GetKernelSources().e_Buffer, DATA_TYPE, QUEUE, std::string(), &error_number)
          )
        ).first;
      BCL_Assert( error_number == CL_SUCCESS, "cannot compile buffer program for: " + util::Format()( DATA_TYPE));

      return itr->second;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Initialize the platform, devices and the commandqueue from the command line flag
    //! This function should be called by all member functions that access queues, platform, or context
    void Tools::UpdateCurrentPlatformDevicesQueuesFromCommandLineFlag()
    {
      GetTools().UpdateCurrentPlatformDevicesQueuesFromCommandLineFlagImpl();
    }

    //! @brief Initialize the platform, devices and the commandqueue from the command line flag
    //! This function should be called by all member functions that access queues, platform, or context
    void Tools::UpdateCurrentPlatformDevicesQueuesFromCommandLineFlagImpl()
    {
      static storage::Vector< std::string> s_platform;

      // if no platform is available, then just return
      if( !Platform::GetPlatformFlag().IsDefined())
      {
        return;
      }

      // if opencl was disabled, then just return
      if( Platform::GetPlatformFlag()->GetFirstParameter()->GetValue() == "Disable")
      {
        return;
      }

      // first, check whether anything has changed
      storage::Vector< std::string> platform_info( Platform::GetPlatformFlag()->GetStringList());
      platform_info.Append( KernelSources::GetKernelsBinaryPathFlag()->GetStringList());
      platform_info.Append( KernelSources::GetKernelsSourcePathFlag()->GetStringList());
      if( platform_info == s_platform)
      {
        return;
      }

      s_platform = platform_info;

      BCL_MessageVrb( "Updating current platform device queues for OpenCL");

      m_BufferPrograms.clear();
      cl_int error( CL_SUCCESS);

      Platform::InitializeFromCommandLine( m_Platform, m_Devices, &error);
      if( error != CL_SUCCESS)
      {
        BCL_MessageCrt( "unable to initialize platform for devices with error: " + ErrorString( error));
        return;
      }

      for( size_t gpus( 0); gpus < m_Devices.GetSize(); ++gpus)
      {
        BCL_MessageDbg
        (
          "UpdateCurrentPlatformDevicesQueuesFromCommandLineFlag() gives this as device( " + util::Format()( gpus) + "): \n\n"
          + util::Format()( m_Devices( gpus).GetDescription())
        );
      }

      BCL_MessageVrb( "platform created: " + m_Platform.Name());

      // update the context
      m_Context = Context( m_Devices, &error);
      if( error != CL_SUCCESS)
      {
        BCL_MessageCrt( "unable to initialize context from devices with error: " + ErrorString( error));
        m_Platform = Platform();
        m_Devices.Reset();
        return;
      }

      // update s_platform, since the device may have been updated
      s_platform = Platform::GetPlatformFlag()->GetStringList();
      s_platform.Append( KernelSources::GetKernelsBinaryPathFlag()->GetStringList());
      s_platform.Append( KernelSources::GetKernelsSourcePathFlag()->GetStringList());

      BCL_MessageDbg( "Initialized opencl context");

      // reset the command queues
      m_CommandQueues.Reset();

      // iterate over devices and create a command queue each
      m_CommandQueues.AllocateMemory( m_Devices.GetSize());
      for
      (
        storage::Vector< Device>::const_iterator dev_itr( m_Devices.Begin()), dev_itr_end( m_Devices.End());
        dev_itr != dev_itr_end;
        ++dev_itr
      )
      {
        m_CommandQueues.PushBack( CommandQueue( m_Context, *dev_itr, &error));
        if( error != CL_SUCCESS)
        {
          BCL_MessageCrt( "unable to initialize command queue from context and device with error: " + ErrorString( error));
          m_CommandQueues.Reset();
          m_Context = Context();
          m_Devices.Reset();
          m_Platform = Platform();

          return;
        }
      }

      // notify listeners that command queue has changed
      if( Platform::GetPlatformFlag().IsDefined())
      {
        m_QueueUpdateSignal.Emit( *this);
      }
      if( GetPlatform().GetPlatformFlag()->GetFlag())
      {
        linal::GetOperationsNonConst< float>().SetDefaultOperationsType( linal::Operations< float>::Operator( "OpenCL"));
        if( GetFirstCommandQueue().GetDevice( NULL).Extensions( NULL).Contains( GetExtensions().e_khr_fp64))
        {
          linal::GetOperationsNonConst< double>().SetDefaultOperationsType( linal::Operations< double>::Operator( "OpenCL"));
        }
      }
    }

    //! list all platforms and devices
    //! @param MESSAGE_LEVEL the message level, to use
    void Tools::ListPlatformsWithDevices( const util::Message::MessageLevel &MESSAGE_LEVEL)
    {
      cl_int error( CL_SUCCESS);

      // get all platforms
      const storage::Vector< Platform> &platform_list( Platform::QueryPlatforms( &error));
      BCL_Assert( error == CL_SUCCESS, "unable to acquire opencl platforms: " + ErrorString( error));

      // iterate over all platforms and display info
      for( storage::Vector< Platform>::const_iterator itr( platform_list.Begin()), itr_end( platform_list.End()); itr != itr_end; ++itr)
      {
        BCL_Message( MESSAGE_LEVEL, "OpenCL platform:\n" + itr->Description());

        // list all devices for that platform
        const storage::Vector< Device> devices( itr->Devices( CL_DEVICE_TYPE_ALL, &error));

        if( error != CL_SUCCESS)
        {
          BCL_MessageCrt( "unable to query devices");
          continue;
        }
        // iterate over devices and display properties
        for( storage::Vector< Device>::const_iterator dev_itr( devices.Begin()), dev_itr_end( devices.End()); dev_itr != dev_itr_end; ++dev_itr)
        {
          BCL_Message( MESSAGE_LEVEL, "device:\n" + dev_itr->GetDescription());
        }
      }
    }

    //! Get and log the binary (PTX) from the OpenCL compiler for the requested program & device
    //! @param cpProgram                   OpenCL program
    //! @param cPtxFileName   optional PTX file name
    void Tools::LogPtx( const cl::Program &cpProgram, const std::string &cPtxFileName)
    {
      cl_int error_number; // Error code var

      std::vector< char *> binaries( cpProgram.getInfo< CL_PROGRAM_BINARIES>( &error_number));
      if( error_number != CL_SUCCESS)
      {
        BCL_MessageCrt( "could not acquire the binaries error: " + ErrorString( error_number));
        return;
      }

      // successfully acquired the binary information
      BCL_MessageDbg( "writing " + util::Format()( binaries.size()) + " binaries for opencl operations to file: " + cPtxFileName);

      // write binaries to file
      io::OFStream write;
      io::File::MustOpenOFStream( write, cPtxFileName);
      // iterate over all binaries and write
      for( std::vector< char *>::iterator itr( binaries.begin()), itr_end( binaries.end()); itr != itr_end; ++itr)
      {
        if( *itr != NULL)
        {
          write.write( *itr, std::strlen( *itr));
          delete[] ( *itr);
        }
      }
      io::File::CloseClearFStream( write);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Tools::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &Tools::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief assign error to a pointer
    //! @param ERROR_PTR pointer to error storage - can be NULL
    //! @param ACTUAL_ERROR the actual error to be assigned
    void Tools::AssignError( cl_int *ERROR_PTR, const cl_int ACTUAL_ERROR)
    {
      if( ERROR_PTR != NULL)
      {
        *ERROR_PTR = ACTUAL_ERROR;
      }
    }

    //! @brief calculate smallest multiple of block size larger equal than global_size needed
    //! @param GROUP_SIZE the number of threads in that block dimension
    //! @param GLOBAL_SIZE the total number of threads needed in that dimension
    //! @return the number of threads needed
    int Tools::RoundUp( const int GROUP_SIZE, const int GLOBAL_SIZE)
    {
      return ( ( ( GLOBAL_SIZE - 1) / GROUP_SIZE) + 1) * GROUP_SIZE;
    }

    //! @brief calculate the padding to add
    //! @param BLOCK_SIZE
    //! @param CURRENT_DIMENSION
    //! @return the amount of padding to add
    size_t Tools::CalcPadding( const size_t BLOCK_SIZE, const size_t CURRENT_DIMENSION)
    {
      return ( BLOCK_SIZE - ( CURRENT_DIMENSION % BLOCK_SIZE)) % BLOCK_SIZE;
    }

    //! @brief compile programs for given precision
    //! @param PRECISION float or double
    //! @param QUEUE the command queue
    //! @param PROGRAM the cl program
    //! @param KERNEL_FILENAME the filename of the kernel
    //! @return ERROR error that occured, CL_SUCCESS if no error
    cl_int Tools::CompilePrograms
    (
      const util::CPPDataTypes::Types &PRECISION,
      const CommandQueue &QUEUE,
      cl::Program &PROGRAM,
      const std::string &KERNEL_FILENAME
    )
    {
      cl_int error_number( CL_SUCCESS);

      // compile the program
      cl::Program::Sources source;
      KernelSourceFile source_file( KERNEL_FILENAME);

      // construct opencl device
      const Device device( QUEUE.GetDevice( &error_number));
      BCL_Assert( error_number == CL_SUCCESS, "cannot get device for command queue");

      BCL_MessageDbg( "device description: " + util::Format()( device.GetDescription()));

      // construct kernel source strings
      const std::string kernels_source( source_file.GetSource( PRECISION, device.Extensions()));

      // check if kernel strings are empty
      if( kernels_source.empty())
      {
        return CL_INVALID_KERNEL_DEFINITION;
      }

      // pushback strings to program sources vector
      source.push_back( std::make_pair( kernels_source.c_str(), kernels_source.length()));

      // create the program
      cl::Program &current_program( PROGRAM);
      current_program = cl::Program( QUEUE.GetContext(), source, &error_number);
      if( error_number != CL_SUCCESS)
      {
        BCL_MessageCrt( "create program error: " + Tools::ErrorString( error_number));
        return error_number;
      }

      // build the program
      error_number = current_program.build( std::vector< cl::Device>( 1, device), s_CLCompilerOptions);
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

    //! @brief common compiler options
    const char *Tools::s_CLCompilerOptions( "-cl-mad-enable -cl-fast-relaxed-math");

    // Helper function to get error string
    const std::string &Tools::ErrorString( cl_int error)
    {
      static const std::string s_error_string[] = {
        "CL_SUCCESS",                         //   0
        "CL_DEVICE_NOT_FOUND",                //  -1
        "CL_DEVICE_NOT_AVAILABLE",            //  -2
        "CL_COMPILER_NOT_AVAILABLE",          //  -3
        "CL_MEM_OBJECT_ALLOCATION_FAILURE",   //  -4
        "CL_OUT_OF_RESOURCES",                //  -5
        "CL_OUT_OF_HOST_MEMORY",              //  -6
        "CL_PROFILING_INFO_NOT_AVAILABLE",    //  -7
        "CL_MEM_COPY_OVERLAP",                //  -8
        "CL_IMAGE_FORMAT_MISMATCH",           //  -9
        "CL_IMAGE_FORMAT_NOT_SUPPORTED",      // -10
        "CL_BUILD_PROGRAM_FAILURE",           // -11
        "CL_MAP_FAILURE",                     // -12
        "",                                   // -13
        "",                                   // -14
        "",                                   // -15
        "",                                   // -16
        "",                                   // -17
        "",                                   // -18
        "",                                   // -19
        "",                                   // -20
        "",                                   // -21
        "",                                   // -22
        "",                                   // -23
        "",                                   // -24
        "",                                   // -25
        "",                                   // -26
        "",                                   // -27
        "",                                   // -28
        "",                                   // -29
        "CL_INVALID_VALUE",                   // -30
        "CL_INVALID_DEVICE_TYPE",             // -31
        "CL_INVALID_PLATFORM",                // -32
        "CL_INVALID_DEVICE",                  // -33
        "CL_INVALID_CONTEXT",                 // -34
        "CL_INVALID_QUEUE_PROPERTIES",        // -35
        "CL_INVALID_COMMAND_QUEUE",           // -36
        "CL_INVALID_HOST_PTR",                // -37
        "CL_INVALID_MEM_OBJECT",              // -38
        "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR", // -39
        "CL_INVALID_IMAGE_SIZE",              // -40
        "CL_INVALID_SAMPLER",                 // -41
        "CL_INVALID_BINARY",                  // -42
        "CL_INVALID_BUILD_OPTIONS",           // -43
        "CL_INVALID_PROGRAM",                 // -44
        "CL_INVALID_PROGRAM_EXECUTABLE",      // -45
        "CL_INVALID_KERNEL_NAME",             // -46
        "CL_INVALID_KERNEL_DEFINITION",       // -47
        "CL_INVALID_KERNEL",                  // -48
        "CL_INVALID_ARG_INDEX",               // -49
        "CL_INVALID_ARG_VALUE",               // -50
        "CL_INVALID_ARG_SIZE",                // -51
        "CL_INVALID_KERNEL_ARGS",             // -52
        "CL_INVALID_WORK_DIMENSION",          // -53
        "CL_INVALID_WORK_GROUP_SIZE",         // -54
        "CL_INVALID_WORK_ITEM_SIZE",          // -55
        "CL_INVALID_GLOBAL_OFFSET",           // -56
        "CL_INVALID_EVENT_WAIT_LIST",         // -57
        "CL_INVALID_EVENT",                   // -58
        "CL_INVALID_OPERATION",               // -59
        "CL_INVALID_GL_OBJECT",               // -60
        "CL_INVALID_BUFFER_SIZE",             // -61
        "CL_INVALID_MIP_LEVEL",               // -62
        "CL_INVALID_GLOBAL_WORK_SIZE"         // -63
      };

      static const int s_error_count( sizeof( s_error_string) / sizeof( s_error_string[ 0]));

      const int index( -error);

      static const std::string s_unknwon( "ERROR_OUT_OF_KNOWN_RANGE");

      // valid error
      if( index >= 0 && index < s_error_count)
      {
        return s_error_string[ index];
      }
      else
      {
        if( error == CL_PLATFORM_NOT_FOUND_KHR)
        {
          static const std::string s_CL_PLATFORM_NOT_FOUND_KHR_string( "CL_PLATFORM_NOT_FOUND_KHR");
          return s_CL_PLATFORM_NOT_FOUND_KHR_string;
        }
        BCL_MessageCrt( "given opencl error out of known range: " + util::Format()( error));
        return s_unknwon;
      }
    }

    //! @brief access to the tools singleton
    //! @return reference to the only instance of Tools
    Tools &GetTools()
    {
      static Tools s_tools;
      return s_tools;
    }

  } // namespace opencl
} // namespace bcl
