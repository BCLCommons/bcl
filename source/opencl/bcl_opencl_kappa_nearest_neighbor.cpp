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
