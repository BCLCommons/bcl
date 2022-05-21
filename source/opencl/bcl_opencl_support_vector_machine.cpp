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
