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
