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
