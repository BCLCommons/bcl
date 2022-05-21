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
#include "model/bcl_model_objective_function_rmsd.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_running_min_max.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionRmsd::s_RawInstance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionRmsd( ObjectiveFunctionRmsd::e_None)
      )
    );
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionRmsd::s_MinMaxInstance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionRmsd( ObjectiveFunctionRmsd::e_MinMax)
      )
    );
    const util::SiPtr< const util::ObjectInterface> ObjectiveFunctionRmsd::s_StdInstance
    (
      util::Enumerated< ObjectiveFunctionInterface>::AddInstance
      (
        new ObjectiveFunctionRmsd( ObjectiveFunctionRmsd::e_Std)
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor from type of normalization
    ObjectiveFunctionRmsd::ObjectiveFunctionRmsd( const Normalization &NORM) :
      m_Normalization( NORM)
    {
    }

    //! copy constructor
    ObjectiveFunctionRmsd *ObjectiveFunctionRmsd::Clone() const
    {
      return new ObjectiveFunctionRmsd( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &ObjectiveFunctionRmsd::GetAlias() const
    {
      static const std::string s_name( "RMSD"), s_normalized_name( "NRMSD"), s_normalized_std_name( "RMSD_NSTD");
      return m_Normalization == e_None ? s_name
             : m_Normalization == e_MinMax ? s_normalized_name
             : s_normalized_std_name;
    }

    //! @brief evaluate for a given datasets with experimental and predicted values the implemented objective function
    //! @param EXPERIMENTAL feature dataset with experimental values
    //! @param PREDICTED feature dataset with predicted values
    //! @return objective function value based on the given data of experimental and predicted values
    float ObjectiveFunctionRmsd::operator()
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
        "Experimental and predicted values do not have the same element size!" + util::Format()( EXPERIMENTAL.GetFeatureSize()) + " " + util::Format()( PREDICTED.GetFeatureSize())
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

      // handle simple normalization
      if( m_Normalization == e_None)
      {
        // runnung average rmsd value
        math::RunningAverage< float> avg_compare;

        // iterate over all rows in dataset
        for( size_t counter( 0); counter < data_set_size; ++counter)
        {
          // runnung rmsd value
          float compare( 0);

          // iterate over all rows in dataset
          for( size_t result_index( 0); result_index < result_size; ++result_index)
          {
            if( util::IsDefined( EXPERIMENTAL( counter)( result_index)))
            {
              // get difference between experimental value and predicted values
              float tmp( EXPERIMENTAL( counter)( result_index) - PREDICTED( counter)( result_index));
              compare += tmp * tmp;
            }
          }

          // add to avg rmsd accumulation
          avg_compare += compare;
        }

        return math::Sqrt( avg_compare.GetAverage() / float( result_size));
      }
      else if( m_Normalization == e_MinMax)
      {
        // track the average difference in each column
        storage::Vector< math::RunningAverage< float> > avg_difference_sq( result_size);

        // track the mins and maxes for each column
        storage::Vector< math::RunningMinMax< float> > min_max( result_size);

        // iterate over all rows in dataset
        for( size_t counter( 0); counter < data_set_size; ++counter)
        {
          // iterate over all rows in dataset
          for( size_t result_index( 0); result_index < result_size; ++result_index)
          {
            if( util::IsDefined( EXPERIMENTAL( counter)( result_index)))
            {
              // get difference between experimental value and predicted values
              float tmp( EXPERIMENTAL( counter)( result_index) - PREDICTED( counter)( result_index));
              avg_difference_sq( result_index) += tmp * tmp;
              min_max( result_index) += EXPERIMENTAL( counter)( result_index);
            }
          }
        }
        // get the average difference, normalized by max-min of each input column
        math::RunningAverage< float> ave_difference_norm;
        for( size_t result_index( 0); result_index < result_size; ++result_index)
        {
          // get this column's range
          const float range( min_max( result_index).GetMax() - min_max( result_index).GetMin());
          if( range > float( 0.0))
          {
            ave_difference_norm += avg_difference_sq( result_index) / math::Sqr( range);
          }
        }
        return math::Sqrt( ave_difference_norm.GetAverage());
      }
      else if( m_Normalization == e_Std)
      {
        // track the average difference in each column
        storage::Vector< math::RunningAverage< float> > avg_difference_sq( result_size);

        // track the mins and maxes for each column
        storage::Vector< math::RunningAverageSD< float> > variances( result_size);

        // iterate over all rows in dataset
        for( size_t counter( 0); counter < data_set_size; ++counter)
        {
          // iterate over all rows in dataset
          for( size_t result_index( 0); result_index < result_size; ++result_index)
          {
            if( util::IsDefined( EXPERIMENTAL( counter)( result_index)))
            {
              // get difference between experimental value and predicted values
              float tmp( EXPERIMENTAL( counter)( result_index) - PREDICTED( counter)( result_index));
              avg_difference_sq( result_index) += tmp * tmp;
              variances( result_index) += EXPERIMENTAL( counter)( result_index);
            }
          }
        }
        // get the average difference, normalized by max-min of each input column
        math::RunningAverage< float> ave_difference_norm;
        std::stringstream oss;
        for( size_t result_index( 0); result_index < result_size; ++result_index)
        {
          // get this column's range
          const float range( variances( result_index).GetVariance());
          if( range > float( 0.0))
          {
            ave_difference_norm += avg_difference_sq( result_index) / range;
            oss << "Result " << result_index << " RMSD_NSTD " << math::Sqrt( avg_difference_sq( result_index) / range)
                << " Variance explained " << 1.0 - avg_difference_sq( result_index) / range << '\n';
          }
        }
        BCL_MessageVrb( "RMSD detail: " + oss.str());
        return math::Sqrt( ave_difference_norm.GetAverage());
      }
      return 0.0;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ObjectiveFunctionRmsd::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Calculates the root-mean-squared-deviation between predicted and actual results"
        "normalized by # of outputs "
        + std::string
          (
            m_Normalization == e_None ? "" :
            m_Normalization == e_MinMax ? "and each result column's range (max-min) " :
            m_Normalization == e_Std ? "and each result column's standard deviation " : ""
          )
      );
      return parameters;
    }

  } // namespace model
} // namespace bcl
