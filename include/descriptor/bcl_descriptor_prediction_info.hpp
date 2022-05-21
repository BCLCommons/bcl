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

// include header of this class
#include "bcl_descriptor_prediction_info.h"
// includes from bcl - sorted alphabetically
#include "bcl_descriptor_combine.h"
#include "bcl_descriptor_iterator.h"
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_vector_reference.h"
#include "math/bcl_math_running_average.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_running_min_max.h"
#include "model/bcl_model_feature_data_set.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_call_stack.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Statistic as string
    //! @param STAT the statistic
    //! @return the string for the stat
    template< typename t_DataType>
    const std::string &PredictionInfo< t_DataType>::GetStatisticName( const Statistic &STAT)
    {
      static const std::string s_names[] =
      {
        "Min",
        "Max",
        "Mean",
        "StandardDeviation",
        GetStaticClassName< Statistic>()
      };
      return s_names[ size_t( STAT)];
    }

    //! @brief virtual copy constructor
    template< typename t_DataType>
    PredictionInfo< t_DataType> *PredictionInfo< t_DataType>::Clone() const
    {
      return new PredictionInfo( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    template< typename t_DataType>
    const std::string &PredictionInfo< t_DataType>::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    template< typename t_DataType>
    const std::string &PredictionInfo< t_DataType>::GetAlias() const
    {
      static const std::string s_name( "PredictionInfo");
      return s_name;
    }

    //! @brief get the type of this descriptor
    //! @return the type of this descriptor (should ignore dimension setting)
    template< typename t_DataType>
    Type PredictionInfo< t_DataType>::GetType() const
    {
      return m_Prediction.GetType();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    template< typename t_DataType>
    size_t PredictionInfo< t_DataType>::GetNormalSizeOfFeatures() const
    {
      return m_Prediction.GetModelsNumberOutputs() * ( m_Statistics.GetSize() + m_Measures.GetSize());
    }

    //! @brief function to return derived-class-held implementations to this interface
    //! This allows this base class to handle mundane details like calling SetDimension and SetObject on all internal
    //! implementations
    template< typename t_DataType>
    iterate::Generic< Base< t_DataType, float> > PredictionInfo< t_DataType>::GetInternalDescriptors()
    {
      return iterate::Generic< Base< t_DataType, float> >( &m_Prediction, &m_Prediction + 1);
    }

    //! @brief load the models: called the 1st time that recalculate is called
    template< typename t_DataType>
    void PredictionInfo< t_DataType>::LoadModels()
    {
      m_Prediction.LoadModels();

      // load roc curve if there are any ROC measures
      if( !m_Measures.IsEmpty())
      {
        // get all the roc curves
        storage::Vector< math::ROCCurve> roc_curves
        (
          m_Prediction.m_ModelStorage->GetMergedIndependentROCCurves()
        );
        // compute roc maps out of them
        m_RocCurves.Reset();
        m_RocCurves.AllocateMemory( roc_curves.GetSize());
        for
        (
          storage::Vector< math::ROCCurve>::const_iterator itr( roc_curves.Begin()), itr_end( roc_curves.End());
          itr != itr_end;
          ++itr
        )
        {
          // compute local roc curve
          m_RocCurves.PushBack( itr->ToMap());
        }

        if( m_Measures.Find( math::ContingencyMatrixMeasures::e_LocalPPV) < m_Measures.GetSize())
        {
          m_LocalPPVs.Reset();
          m_LocalPPVs.AllocateMemory( m_RocCurves.GetSize());
          for
          (
            storage::Vector< math::ROCCurve>::const_iterator itr( roc_curves.Begin()), itr_end( roc_curves.End());
            itr != itr_end;
            ++itr
          )
          {
            // compute local roc curve
            m_LocalPPVs.PushBack( itr->GetLocalPPVCurve());
          }
        }
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    template< typename t_DataType>
    io::Serializer PredictionInfo< t_DataType>::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Computes prediction-related information, including standard deviation, min, max, and ROC-curve derived "
        "statistics such as the PPV or local-PPV associated with a prediction. "
        "Output will be statistics first, in the order provided; followed by metrics, in the order they are provided"
      );
      parameters.AddInitializer
      (
        "predictor",
        "predictor to obtain derived info from",
        io::Serialization::GetAgent( &m_Prediction)
      );
      parameters.AddOptionalInitializer
      (
        "statistics",
        "statistics to compute for the predictions",
        io::Serialization::GetAgent( &m_Statistics)
      );
      parameters.AddOptionalInitializer
      (
        "metrics",
        "ROC-curve based metrics to compute for each model's output. "
        "e.g. LocalPPV will give you the likelihood that a prediction is a true positive",
        io::Serialization::GetAgent( &m_Measures)
      );

      return parameters;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief calculate the descriptors for a given sequence iterator position, ignoring cache
    //! @param ITR sequence iterator pointing to the tuple of sequence elements of interest; will be for the native dimension of the class
    //! @param STORAGE storage for the descriptor
    //! The inheriting class can assume that ITR is of native dimension and that STORAGE is zeroed and of the proper
    //! dimension
    template< typename t_DataType>
    void PredictionInfo< t_DataType>::RecalculateImpl
    (
      const Iterator< t_DataType> &ITR,
      linal::VectorReference< float> &STORAGE
    )
    {
      // load the models, if they are not already loaded
      if( m_Prediction.m_Models.IsEmpty())
      {
        LoadModels();
      }

      math::RunningAverageSD< linal::Vector< float> > ave_sd;
      math::RunningMinMax< linal::Vector< float> > min_max;
      storage::Vector< math::RunningAverage< linal::Vector< float> > > ave_roc_measures( m_Measures.GetSize());

      const size_t number_models( m_Prediction.m_NumberModels);
      size_t features_per_model
      (
        std::max( size_t( 1), m_Prediction.m_NumberOutputs)
        / number_models
      );

      linal::VectorConstReference< float> all_values_ref( m_Prediction( ITR));
      size_t storage_pos( 0);
      linal::Vector< float> measure_values( features_per_model, float( 0.0));
      for( size_t model_number( 0); model_number < number_models; ++model_number, storage_pos += features_per_model)
      {
        // get the feature reference
        linal::VectorConstReference< float> pred_this_model
        (
          all_values_ref.CreateSubVectorConstReference( features_per_model, storage_pos)
        );
        ave_sd += pred_this_model;
        min_max += pred_this_model;
        for( size_t measure_id( 0), n_measures( m_Measures.GetSize()); measure_id < n_measures; ++measure_id)
        {
          const math::ContingencyMatrixMeasures measure( m_Measures( measure_id));
          // compute the measure for each of the outputs
          if( measure.GetMeasure() == math::ContingencyMatrixMeasures::e_LocalPPV)
          {
            for( size_t output_id( 0); output_id < features_per_model; ++output_id)
            {
              measure_values( output_id) = util::IsDefined( pred_this_model( output_id)) ? m_LocalPPVs( output_id)( pred_this_model( output_id)) : 0;
            }
          }
          else
          {
            for( size_t output_id( 0); output_id < features_per_model; ++output_id)
            {
              // prediction for this model, for this output
              const double prediction( pred_this_model( output_id));

              // compute lower bound, which is first element that is >= the given value
              storage::Map< double, math::ContingencyMatrix>::const_iterator itr_map_gte
              (
                m_RocCurves( output_id).LowerBound( prediction)
              );

              // compute first element <= the given value
              storage::Map< double, math::ContingencyMatrix>::const_iterator itr_map_lte( itr_map_gte);

              // check whether we need to decrement
              if( itr_map_lte != m_RocCurves( output_id).Begin() && itr_map_lte->first > prediction)
              {
                --itr_map_lte;
              }

              // compute the measure at each point
              const double measure_lte( measure( itr_map_lte->second));
              const double measure_gte( measure( itr_map_gte->second));

              measure_values( output_id)
                = math::LinearFunction( itr_map_lte->first, measure_lte, itr_map_gte->first, measure_gte)( prediction);
            }
          }
          ave_roc_measures( measure_id) += measure_values;
        }
      }

      // store statistics first
      size_t storage_place( 0);
      for
      (
        typename storage::Vector< StatisticEnum>::const_iterator
          itr_statistic( m_Statistics.Begin()), itr_statistic_end( m_Statistics.End());
        itr_statistic != itr_statistic_end;
        ++itr_statistic, storage_place += features_per_model
      )
      {
        linal::VectorReference< float> storage_ref
        (
          STORAGE.CreateSubVectorReference( features_per_model, storage_place)
        );
        switch( itr_statistic->GetEnum())
        {
          case e_Mean: storage_ref.CopyValues( ave_sd.GetAverage()); break;
          case e_StandardDeviation: storage_ref.CopyValues( ave_sd.GetStandardDeviation()); break;
          case e_Min: storage_ref.CopyValues( min_max.GetMin()); break;
          case e_Max: storage_ref.CopyValues( min_max.GetMax()); break;
          default: break;
        }
      }
      for
      (
        size_t measure_id( 0), n_measures( m_Measures.GetSize());
        measure_id < n_measures;
        ++measure_id, storage_place += features_per_model
      )
      {
        linal::VectorReference< float> storage_ref
        (
          STORAGE.CreateSubVectorReference( features_per_model, storage_place)
        );
        storage_ref.CopyValues( ave_roc_measures( measure_id).GetAverage());
      }
    }

  } // namespace descriptor
} // namespace bcl
