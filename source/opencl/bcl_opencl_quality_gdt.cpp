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
