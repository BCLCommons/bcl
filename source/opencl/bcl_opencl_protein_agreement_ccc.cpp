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
