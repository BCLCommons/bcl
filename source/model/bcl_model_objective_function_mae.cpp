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
#include "model/bcl_model_objective_function_mae.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_running_average.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionMae::s_Instance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionMae( e_None)
      )
    );
    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionMae::s_FractionExplainedInstance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionMae( e_FractionExplained)
      )
    );
    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionMae::s_MADInstance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionMae( e_MeanAbsoluteDeviation)
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ObjectiveFunctionMae::ObjectiveFunctionMae( const Normalization &NORM) :
      m_Normalization( NORM)
    {
    }

    //! copy constructor
    ObjectiveFunctionMae *ObjectiveFunctionMae::Clone() const
    {
      return new ObjectiveFunctionMae( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ObjectiveFunctionMae::GetAlias() const
    {
      static const std::string s_name( "MAE"), s_fe_mae( "MAE_FractionExplained"), s_mad( "MAE_NMAD");
      return m_Normalization == e_None ? s_name : m_Normalization == e_FractionExplained ? s_fe_mae : s_mad;
    }

    //! @brief set dataset; some objective functions may need to setup internal data structures based on this function
    //! @param DATA monitoring dataset results, non-scaled
    //! @param IDS ids; can be used by the objective function
    void ObjectiveFunctionMae::SetData
    (
      const FeatureDataSet< float> &DATA,
      const FeatureDataSet< char> &IDS
    )
    {
      const size_t n_results( DATA.GetFeatureSize());
      storage::Vector< float> averages( n_results);
      if( m_Normalization == e_MeanAbsoluteDeviation || m_Normalization == e_FractionExplained)
      {
        storage::Vector< math::RunningAverage< float> > average_result_vals( n_results);
        for( size_t row_id( 0), n_rows( DATA.GetNumberFeatures()); row_id < n_rows; ++row_id)
        {
          for( size_t res_id( 0); res_id < n_results; ++res_id)
          {
            if( util::IsDefined( DATA( row_id)( res_id)))
            {
              average_result_vals( res_id) += DATA( row_id)( res_id);
            }
          }
        }
        for( size_t res_id( 0); res_id < n_results; ++res_id)
        {
          averages( res_id) = average_result_vals( res_id).GetAverage();
          average_result_vals( res_id).Reset();
        }
        for( size_t row_id( 0), n_rows( DATA.GetNumberFeatures()); row_id < n_rows; ++row_id)
        {
          for( size_t res_id( 0); res_id < n_results; ++res_id)
          {
            if( util::IsDefined( DATA( row_id)( res_id)))
            {
              average_result_vals( res_id) += math::Absolute( DATA( row_id)( res_id) - averages( res_id));
            }
          }
        }
        m_AverageDeviations.Resize( n_results);
        for( size_t res_id( 0); res_id < n_results; ++res_id)
        {
          m_AverageDeviations( res_id) = average_result_vals( res_id).GetAverage();
        }
      }
    }

    //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return objective function value based on the given data of experimental and predicted values
    float ObjectiveFunctionMae::operator()
    (
      const FeatureDataSetInterface< float> &EXPERIMENTAL,
      const FeatureDataSetInterface< float> &PREDICTED
    ) const
    {
      BCL_Assert
      (
        EXPERIMENTAL.GetNumberFeatures() == PREDICTED.GetNumberFeatures(),
        "Experimental and predicted values do not have the same number of elements!"
      );

      BCL_Assert
      (
        EXPERIMENTAL.GetFeatureSize() == PREDICTED.GetFeatureSize(),
        "Experimental and predicted values do not have the same number of elements!"
      );

      // number of experimental values
      const size_t data_set_size( EXPERIMENTAL.GetNumberFeatures());

      // number of predicted result columns
      const size_t result_size( PREDICTED.GetFeatureSize());

      if( data_set_size == 0 || result_size == 0)
      {
        // no data available, no result column given
        return float( 0);
      }

      // runnung average rmsd value
      storage::Vector< math::RunningAverage< float> > avg_compare( result_size);

      // iterate over all rows in dataset
      for( size_t counter( 0); counter < data_set_size; ++counter)
      {
        // iterate over all rows in dataset
        for( size_t result_index( 0); result_index < result_size; ++result_index)
        {
          if( util::IsDefined( EXPERIMENTAL( counter)( result_index)))
          {
            // get difference between experimental value and predicted values
            avg_compare( result_index)
              += math::Absolute( EXPERIMENTAL( counter)( result_index) - PREDICTED( counter)( result_index));
          }
        }
      }
      math::RunningAverage< float> avg_sum;
      for( size_t result_index( 0); result_index < result_size; ++result_index)
      {
        if( m_Normalization == e_None)
        {
          avg_sum += avg_compare( result_index).GetAverage();
          BCL_MessageStd
          (
            "MAE Detail Result Index: " + util::Format()( result_index)
            + " MAE: " + util::Format()( avg_compare( result_index).GetAverage())
          );
        }
        else if( m_Normalization == e_MeanAbsoluteDeviation)
        {
          avg_sum += avg_compare( result_index).GetAverage() / m_AverageDeviations( result_index);
          BCL_MessageStd
          (
            "MAE Detail Result Index: " + util::Format()( result_index)
            + " MAE: " + util::Format()( avg_compare( result_index).GetAverage())
            + " MAD: " + util::Format()( m_AverageDeviations( result_index))
            + " MAE_NormByMAD: " +
            util::Format()
            (
              avg_compare( result_index).GetAverage() / m_AverageDeviations( result_index)
            )
          );
        }
        else if( m_Normalization == e_FractionExplained)
        {
          avg_sum += 1.0 - avg_compare( result_index).GetAverage() / m_AverageDeviations( result_index);
          BCL_MessageStd
          (
            "MAE Detail Result Index: " + util::Format()( result_index)
            + " MAE: " + util::Format()( avg_compare( result_index).GetAverage())
            + " MAD: " + util::Format()( m_AverageDeviations( result_index))
            + " MAD_FractionExplained: " +
            util::Format()
            (
              1.0 - avg_compare( result_index).GetAverage() / m_AverageDeviations( result_index)
            )
          );
        }
      }

      return avg_sum.GetAverage();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionMae::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        m_Normalization == e_None
        ? "Calculates the mean-absolute-deviation between predicted and actual results"
        : (
             m_Normalization == e_FractionExplained
             ? "Calculates the fraction of absolute deviation explained"
             : "Mean absolute error divided by mean absolute deviation of each result column in the dataset"
           )
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
