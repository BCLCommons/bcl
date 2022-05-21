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
