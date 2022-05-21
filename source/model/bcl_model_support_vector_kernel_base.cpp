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
#include "model/bcl_model_support_vector_kernel_base.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_range.h"
#include "model/bcl_model_feature_data_set.h"
#include "sched/bcl_sched_scheduler_interface.h"
#include "sched/bcl_sched_tertiary_function_job_with_data.h"
#include "storage/bcl_storage_triplet.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //! @brief compute kernel values for every pair input_vector_i with each other vector in data set
    //! @param TRAINING_DATA featurevector data set without labels
    //! @param INPUT_VECTOR_I feature vector
    //! @return vector with all combinations towards INPUT_VECTOR_I
    linal::Vector< float> SupportVectorKernelBase::GetInputVectorIKernelMatrixMultiOutput
    (
      const FeatureDataSetInterface< float> &TRAINING_DATA,
      const size_t &INPUT_VECTOR_I
    ) const
    {
      // size of data set
      const size_t data_set_size( TRAINING_DATA.GetNumberFeatures());
      BCL_Assert( INPUT_VECTOR_I < data_set_size, "vector i is outside the data set");

      // create storage::Vector of size corresponding to entries in TRAINING_DATA initialized with 0.0
      linal::Vector< float> kernel_data( data_set_size, float( 0.0));
      // reference to index_vector_i in training data
      const FeatureReference< float> input_vector_i( TRAINING_DATA( INPUT_VECTOR_I));

      // define a group id for threads that will be created
      const size_t group_id( 0);

      // number of available cpus
      const size_t number_available_cpus( sched::GetNumberCPUs());

      const size_t interval( data_set_size / number_available_cpus);

      BCL_Assert( interval != 0, "Interval is zero! Less data points then processors!");

      // number of jobs
      size_t number_jobs( number_available_cpus);

      storage::Triplet
      <
        util::SiPtr< const FeatureDataSetInterface< float> >,
        const FeatureReference< float>,
        linal::Vector< float>
      > train_index_result
      (
        util::SiPtr< const FeatureDataSetInterface< float> >( TRAINING_DATA),
        input_vector_i,
        kernel_data
      );

      util::ShPtrList< sched::JobInterface> schedule;

      // ranges for jobs
      storage::Vector< math::Range< size_t> > ranges;
      ranges.AllocateMemory( number_available_cpus);

      // create ranges
      for
      (
        size_t job_number( 0), start( 0), between( interval);
        job_number < number_jobs;
        ++job_number, start = between + 1, between += interval
      )
      {
        // if last job the extend till dataset size
        if( job_number + 1 == number_jobs)
        {
          between = data_set_size - 1;
        }
        // push back result
        ranges.PushBack( math::Range< size_t>( start, between));
      }

      // for every interval create a job
      for
      (
        storage::Vector< math::Range< size_t> >::const_iterator itr_ranges( ranges.Begin()), itr_ranges_end( ranges.End());
        itr_ranges != itr_ranges_end;
        ++itr_ranges
      )
      {
        // create and pushback jobs into scheduler
        schedule.PushBack
        (
          util::ShPtr< sched::JobInterface>
          (
            new sched::TertiaryFunctionJobWithData
            <
              storage::Triplet
              <
                util::SiPtr< const FeatureDataSetInterface< float> >,
                const FeatureReference< float>,
                linal::Vector< float>
              >,
              const size_t,
              const size_t,
              void,
              SupportVectorKernelBase
            >
            (
              group_id,
              *this,
              &SupportVectorKernelBase::ComputePartialInputVectorIKernelMatrix,
              train_index_result,
              itr_ranges->GetMin(),
              itr_ranges->GetMax(),
              sched::JobInterface::e_READY,
              NULL
            )
          )
        );
      }

      // submit jobs
      for
      (
        util::ShPtrList< sched::JobInterface>::iterator itr_jobs( schedule.Begin()), itr_jobs_end( schedule.End());
        itr_jobs != itr_jobs_end;
        ++itr_jobs
      )
      {
        // submit job for descriptor combination
        sched::GetScheduler().RunJob( *itr_jobs);
      }

      // join all jobs ( need to be finished)
      for
      (
        util::ShPtrList< sched::JobInterface>::iterator itr_jobs( schedule.Begin()), itr_jobs_end( schedule.End());
        itr_jobs != itr_jobs_end;
        ++itr_jobs
      )
      {
        sched::GetScheduler().Join( *itr_jobs);
      }

      // initialize kernel data
      const linal::Vector< float> &kernel_values( train_index_result.Third());
      return kernel_values;
    }

        //! @brief compute kernel values for every pair input_vector_i with each other vector in data set
    //! @param TRAINING_DATA featurevector data set without labels
    //! @param INPUT_VECTOR_I feature vector
    //! @param SIGNS
    //! @param PROBLENGTH
    //! @return vector with all combinations towards INPUT_VECTOR_I
    linal::Vector< float> SupportVectorKernelBase::GetInputVectorIKernelMatrix
    (
      const FeatureDataSetInterface< float> &TRAINING_DATA,
      const size_t &INPUT_VECTOR_I,
      const storage::Vector< int> &SIGNS,
      const size_t &PROBLENGTH
    ) const
    {
      // size of data set
      const size_t data_set_size( TRAINING_DATA.GetNumberFeatures());

      // create storage::Vector of size corresponding to entries in TRAINING_DATA initialized with 0.0
      linal::Vector< float> kernel_data( data_set_size, float( 0.0));

      // calculate index value for input vector i
      size_t index_vector_i( 0);

      // adjust index for input vector i
      if( INPUT_VECTOR_I > ( data_set_size - 1))
      {
        index_vector_i = ( INPUT_VECTOR_I - PROBLENGTH / 2);
      }
      else
      {
        index_vector_i = INPUT_VECTOR_I;
      }

      // reference to index_vector_i in training data
      const FeatureReference< float> input_vector_i( TRAINING_DATA( index_vector_i));

      // define a group id for threads that will be created
      const size_t group_id( 0);

      // number of available cpus
      const size_t number_available_cpus( sched::GetNumberCPUs());

      const size_t interval( data_set_size / number_available_cpus);

      BCL_Assert( interval != 0, "Interval is zero! Less data points then processors!");

      // number of jobs
      size_t number_jobs( number_available_cpus);

      storage::Triplet
      <
        util::SiPtr< const FeatureDataSetInterface< float> >,
        const FeatureReference< float>,
        linal::Vector< float>
      > train_index_result
      (
        util::SiPtr< const FeatureDataSetInterface< float> >( TRAINING_DATA),
        input_vector_i,
        kernel_data
      );

      util::ShPtrList< sched::JobInterface> schedule;

      // ranges for jobs
      storage::Vector< math::Range< size_t> > ranges;
      ranges.AllocateMemory( number_available_cpus);

      // create ranges
      for
      (
        size_t job_number( 0), start( 0), between( interval);
        job_number < number_jobs;
        ++job_number, start = between + 1, between += interval
      )
      {
        // if last job the extend till dataset size
        if( job_number + 1 == number_jobs)
        {
          between = data_set_size - 1;
        }
        // push back result
        ranges.PushBack( math::Range< size_t>( start, between));
      }

      // for every interval create a job
      for
      (
        storage::Vector< math::Range< size_t> >::const_iterator itr_ranges( ranges.Begin()), itr_ranges_end( ranges.End());
        itr_ranges != itr_ranges_end;
        ++itr_ranges
      )
      {
        // create and pushback jobs into scheduler
        schedule.PushBack
        (
          util::ShPtr< sched::JobInterface>
          (
            new sched::TertiaryFunctionJobWithData
            <
              storage::Triplet
              <
                util::SiPtr< const FeatureDataSetInterface< float> >,
                const FeatureReference< float>,
                linal::Vector< float>
              >,
              const size_t,
              const size_t,
              void,
              SupportVectorKernelBase
            >
            (
              group_id,
              *this,
              &SupportVectorKernelBase::ComputePartialInputVectorIKernelMatrix,
              train_index_result,
              itr_ranges->GetMin(),
              itr_ranges->GetMax(),
              sched::JobInterface::e_READY,
              NULL
            )
          )
        );
      }

      // submit jobs
      for
      (
        util::ShPtrList< sched::JobInterface>::iterator itr_jobs( schedule.Begin()), itr_jobs_end( schedule.End());
        itr_jobs != itr_jobs_end;
        ++itr_jobs
      )
      {
        // submit job for descriptor combination
        sched::GetScheduler().RunJob( *itr_jobs);
      }

      // join all jobs ( need to be finished)
      for
      (
        util::ShPtrList< sched::JobInterface>::iterator itr_jobs( schedule.Begin()), itr_jobs_end( schedule.End());
        itr_jobs != itr_jobs_end;
        ++itr_jobs
      )
      {
        sched::GetScheduler().Join( *itr_jobs);
      }

      // initialize kernel data
      const linal::Vector< float> &kernel_values( train_index_result.Third());

      // initialize vector for kernel results
      linal::Vector< float> result_data( PROBLENGTH, 0.0);
      // iterator over all signs on training data
      storage::Vector< int>::const_iterator itr_begin_signs( SIGNS.Begin());

      size_t index( 0);
      const size_t problem_real_size( PROBLENGTH / 2);
      const int sign_vector_i( SIGNS( INPUT_VECTOR_I));

      for
      (
        linal::Vector< float>::iterator
          itr_begin_result_data( result_data.Begin()),
          itr_end_result_data( result_data.End());
        itr_begin_result_data != itr_end_result_data;
        ++itr_begin_result_data, ++itr_begin_signs, ++index
      )
      {
        // decide whether you reach outside of the
        if( index >= problem_real_size)
        {
          *itr_begin_result_data = sign_vector_i * ( *itr_begin_signs) * kernel_values( index - problem_real_size);
        }
        else
        {
          *itr_begin_result_data = sign_vector_i * ( *itr_begin_signs) * kernel_values( index);
        }
      }

      return result_data;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief compute kernel matrix for one particular feature with all feature values in training data
    //! @param TRAIN_INDEX_RESULT data structure which contains pointer to training data, current vector,
    //!        and result vector
    void SupportVectorKernelBase::ComputePartialInputVectorIKernelMatrix
    (
      storage::Triplet
      <
        util::SiPtr< const FeatureDataSetInterface< float> >,
        const FeatureReference< float>,
        linal::Vector< float>
      > &TRAIN_INDEX_RESULT,
      const size_t &START_INDEX,
      const size_t &END_INDEX
    ) const
    {
      // starting index of current feature vector in training dataset
      size_t index_vector_j( START_INDEX);

      // iterate over range of result vector
      for
      (
        float *itr_result( TRAIN_INDEX_RESULT.Third().Begin() + START_INDEX);
        index_vector_j <= END_INDEX;
        ++itr_result, ++index_vector_j
      )
      {
        *itr_result = operator()
        (
          TRAIN_INDEX_RESULT.Second(), TRAIN_INDEX_RESULT.First()->operator()( index_vector_j)
        );
      }
    }

  } // namespace model
} // namespace bcl
