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

#ifndef BCL_MODEL_DATA_SET_STATISTICS_H_
#define BCL_MODEL_DATA_SET_STATISTICS_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_feature_label_set.h"
#include "bcl_model_retrieve_data_set_base.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "linal/bcl_linal_vector.h"
#include "linal/bcl_linal_vector_operations.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_running_min_max.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataSetStatistics
    //! @brief interface to retrieve data sets independent of source
    //!
    //! @see @link example_model_data_set_statistics.cpp @endlink
    //! @author mendenjl
    //! @date Aug 17, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DataSetStatistics :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      math::RunningAverageSD< linal::Vector< float> > m_AveStdFeatures; //!< Averages and standard deviations of the features
      math::RunningMinMax< linal::Vector< float> >    m_MinMaxFeatures; //!< Min and max of the features
      math::RunningAverageSD< linal::Vector< float> > m_AveStdResults;  //!< Averages and standard deviations of the results
      math::RunningMinMax< linal::Vector< float> >    m_MinMaxResults;  //!< Min and max of the results

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new DataSetStatistics
      DataSetStatistics *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief Get the averages and stddevs of each feature
      //! @return the averages and stddevs of each feature
      const math::RunningAverageSD< linal::Vector< float> > &GetAveStdFeatures() const
      {
        return m_AveStdFeatures;
      }

      //! @brief Get the averages and stddevs of each result
      //! @return the averages and stddevs of each result
      const math::RunningAverageSD< linal::Vector< float> > &GetAveStdResults() const
      {
        return m_AveStdResults;
      }

      //! @brief Get the mins and maxs of each feature
      //! @return the mins and maxs of each feature
      const math::RunningMinMax< linal::Vector< float> > &GetMinMaxFeatures() const
      {
        return m_MinMaxFeatures;
      }

      //! @brief Get the mins and maxs of each result
      //! @return the mins and maxs of each result
      const math::RunningMinMax< linal::Vector< float> > &GetMinMaxResults() const
      {
        return m_MinMaxResults;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief compute the data set statistics from a dataset retriever
      //! @param DATASET some means of making a dataset retriever
      void AddDataSet( util::Implementation< RetrieveDataSetBase> DATASET);

      //! @brief compute the data set statistics from a dataset
      //! @param DATASET a dataset
      void AddDataSet( const descriptor::Dataset &DATASET);

      //! @brief compute the data set statistics as the difference between two other datasets
      //! @param DATASET_A, DATASET_B actual datasets
      void AddDifference
      (
        const descriptor::Dataset &DATASET_A,
        const descriptor::Dataset &DATASET_B
      );

      //! @brief compute the data set statistics as the difference between two other datasets
      //! @param DATASET_A, DATASET_B some means of making a dataset
      void AddDifference
      (
        util::Implementation< RetrieveDataSetBase> DATASET_A,
        util::Implementation< RetrieveDataSetBase> DATASET_B
      );

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class DataSetStatistics

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_DATA_SET_STATISTICS_H_
