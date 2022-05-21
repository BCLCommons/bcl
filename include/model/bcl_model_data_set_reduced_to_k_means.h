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

#ifndef BCL_MODEL_DATA_SET_REDUCED_TO_K_MEANS_H_
#define BCL_MODEL_DATA_SET_REDUCED_TO_K_MEANS_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_retrieve_data_set_base.h"
#include "bcl_model_training_schedule.h"
#include "linal/bcl_linal_matrix.h"
#include "linal/bcl_linal_vector.h"
#include "linal/bcl_linal_vector_const_reference.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataSetReducedToKMeans
    //! @brief reduces a data set k-clusters, chosen to minimize the standard deviation between the features and clusters
    //!
    //! @author loweew, mendenjl
    //! @see @link example_model_data_set_reduced_to_k_means.cpp @endlink
    //! @date 12/03/10
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DataSetReducedToKMeans :
      public RetrieveDataSetBase
    {
    private:

    //////////
    // data //
    //////////

      size_t m_NumberClusters; //!< # of clusters desired
      size_t m_MaxReclusterAttempts; //!< max # attempts to improve cluster
      bool   m_Autoscale; //!< true if the input data need to be scaled to 0-1 before clustering
      util::Implementation< RetrieveDataSetBase> m_Retriever; //!< Retrieval method for the data set to be clustered

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param NUMBER_CLUSTERS the desired number of clusters
      //! @param MAX_RECLUSTER_ATTEMPTS the maximum number of times to try improving the clusters
      //! @param AUTOSCALE true if the data need to be normalized before clustering
      DataSetReducedToKMeans
      (
        const size_t &NUMBER_CLUSTERS = 0,
        const size_t &MAX_RECLUSTER_ATTEMPTS = 20,
        const bool &AUTOSCALE = false
      );

      //! @brief Clone function
      //! @return pointer to new DataSetReducedToKMeans
      DataSetReducedToKMeans *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

      //! @brief return the maximum number of clusters that will be returned
      //! @return the maximum number of clusters that will be returned
      //! the actual number of clusters will be no larger than the size of the data set
      size_t GetNumberClusters() const;

      //! @brief Set the code / label for the feature (1st part) of the data set
      //! @param CODE the new code
      //! @return the code / label for the feature (1st part) of the data set
      void SelectFeatures( const util::ObjectDataLabel &CODE);

      //! @brief Set the code / label for the result (2nd part) of the data set
      //! @return the code / label for the result (2nd part) of the data set
      void SelectResults( const util::ObjectDataLabel &CODE);

      //! @brief Set which id columns to retrieve
      //! @param CODE the id column names to retrieve
      void SelectIds( const util::ObjectDataLabel &CODE);

      //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
      //! @return the code / label for the feature (1st part) of the data set with sizes of each property
      //! the feature code set
      FeatureLabelSet GetFeatureLabelsWithSizes() const;

      //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
      //! @return the code / label for the result (2nd part) of the data set with sizes of each property
      //! the feature code set
      FeatureLabelSet GetResultCodeWithSizes() const;

      //! @brief Get the code / label for the ids of the data set with sizes of each property
      //! @return the code / label for the ids of the data set with sizes of each property
      //! the feature code set
      FeatureLabelSet GetIdCodeWithSizes() const;

      //! @brief get whether dataset generation requires labels
      //! @return true if dataset generation requires labels
      bool RequiresFeatureLabels() const
      {
        return m_Retriever->RequiresFeatureLabels();
      }

      //! @brief get whether dataset generation requires result labels
      //! @return true if dataset generation requires result labels
      bool RequiresResultLabels() const
      {
        return m_Retriever->RequiresResultLabels();
      }

      //! @brief get the number of partitions requested by the user, along with the partition ids
      //! @return the number of partitions requested by the user, along with the partition ids
      storage::Pair< size_t, math::RangeSet< size_t> > GetNumberPartitionsAndIds() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief compute the square of the distance between two vectors, so long as it is no greater than LIMIT
      //! @param VEC_A, VEC_B the two vectors
      //! @param LIMIT the maximum distance of interest
      //! @return the distance, so long as it is less than LIMIT, otherwise, a value larger than LIMIT
      static float LimitedSquareDistance
      (
        const linal::VectorConstInterface< float> &VEC_A,
        const linal::VectorConstInterface< float> &VEC_B,
        const float &LIMIT
      );

    ////////////////
    // operations //
    ////////////////

      //! @brief generate dataset
      //! @return generated dataset
      util::ShPtr< descriptor::Dataset> GenerateDataSet();

    ///////////////
    // operators //
    ///////////////

      //! @brief reduce a data set
      //! @param DATA the data set
      //! @return the reduced data set
      util::ShPtr< descriptor::Dataset> operator ()
      (
        const util::ShPtr< descriptor::Dataset> &DATA,
        const TrainingSchedule &BALANCING = TrainingSchedule()
      ) const;

      //!@brief returns indices of objects in cluster
      storage::Vector< size_t> IdentifyClusterIndices( util::ShPtr< descriptor::Dataset> &REDUCED_MATRIX) const;

    protected:

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief refines the current set of clusters
      //! @param FEATURES the features to cluster
      //! @param FEATURE_CLUSTER_INDICES the index of the closest mean in the clusters last round, for each feature
      //! @param CLUSTERS the current vector of clusters, which will be updated
      //! @return true if the clusters changed
      bool RefineClusters
      (
        const storage::Vector< storage::VectorND< 2, linal::VectorConstReference< float> > > &FEATURES,
        storage::Vector< size_t> &FEATURE_CLUSTER_INDICES,
        storage::Vector< math::RunningAverage< linal::Vector< float> > > &CLUSTERS,
        const TrainingSchedule &TRAINING_SCHEDULE
      ) const;

      //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
      //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
      virtual size_t GetNominalSize() const;

    }; // class DataSetReducedToKMeans

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_DATA_SET_REDUCED_TO_K_MEANS_H_

