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

#ifndef BCL_MODEL_DATA_SET_REDUCED_TO_CLUSTER_CENTERS_H_
#define BCL_MODEL_DATA_SET_REDUCED_TO_CLUSTER_CENTERS_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_retrieve_data_set_base.h"
#include "math/bcl_math_range.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataSetReducedToClusterCenters
    //! @brief reduces a data set to a series of cluster centers that cover the same space
    //! There may be a smaller data set that covers the same space, which may be obtained using k-means
    //! This method should be faster, and can be used on its own or to provide the initial cluster points
    //! for the k-means algorithm
    //!
    //! @author mendenjl
    //! @see @link example_model_data_set_reduced_to_cluster_centers.cpp @endlink
    //! @date 12/01/10
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DataSetReducedToClusterCenters :
      public RetrieveDataSetBase
    {
    public:

      //! Which member of the cluster to select
      enum Selection
      {
        e_First,  //!< Fastest option, just picks the 1st observed member of the cluster
        e_Center, //!< Computes the centermost based on feature similarity
        e_Max,    //!< Chooses the maximum result; sorted by first element (then second element if 1st elements are equal, etc.)
        e_Min,    //!< Chooses the minimum result, sorted by first element (then second element if 1st elements are equal, etc.)
        s_NumberSelections
      };

      //! @brief get the name associated with a particular selection
      //! @param SELECTION the selection type of interest
      //! @return the name associated with that selection type
      static const std::string &GetSelectionName( const Selection &SELECTION);

      //! @brief Selection enum I/O helper
      typedef util::WrapperEnum< Selection, &GetSelectionName, s_NumberSelections> SelectionEnum;

    private:

    //////////
    // data //
    //////////

      double m_RMSD;      //!< rmsd's between feature vectors that are less than m_Rmsd are part of the same cluster
      math::Range< size_t> m_ClusterSizeRange; //!< the min and max # of clusters that are be acceptable
      size_t m_MaxSteps;  //!< the maximum # of attempts to change the rmsd to get the # of clusters into the range
      bool   m_Autoscale; //!< true if the input data need to be scaled to 0-1 before clustering
      util::Implementation< RetrieveDataSetBase> m_Retriever; //!< Retrieval method for the data set to be clustered
      SelectionEnum m_Selection; //!< True if the most similar member of the cluster should be located

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from parameters
      //! @param RMSD the rmsd to use in the clustering algorithm
      //! @param CLUSTER_RANGE the min and max # of clusters that would be acceptable
      //! @param MAX_STEPS_TO_REACH_CLUSTER_SIZE the maximum # of attempts to change the rmsd to get the # of clusters into the range
      //! @param AUTOSCALE true if the data need to be normalized before clustering
      DataSetReducedToClusterCenters
      (
        const double &RMSD = 0.0,
        const math::Range< size_t> &RANGE = math::Range< size_t>( 1, std::numeric_limits< size_t>::max()),
        const size_t &MAX_STEPS_TO_REACH_CLUSTER_SIZE = 10,
        const bool &AUTOSCALE = false,
        const Selection &SELECTION = e_First
      );

      //! @brief Clone function
      //! @return pointer to new DataSetReducedToClusterCenters
      DataSetReducedToClusterCenters *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the rmsd parameter
      //! @return the rmsd parameter
      double GetRMSDParameter() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

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

      //! @brief generate dataset
      //! @return generated dataset
      util::ShPtr< descriptor::Dataset> GenerateDataSet();

    ///////////////
    // operators //
    ///////////////

      //! @brief reduce a data set
      //! @param DATA the data set
      //! @return the reduced data set
      void operator ()( util::ShPtr< descriptor::Dataset> DATA) const;

      //! @brief reduce a data set
      //! @param DATA the data set
      //! @return the reduced data set
      util::ShPtr< descriptor::Dataset> operator ()( const descriptor::Dataset &DATA) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief get a set of features that cover the same space as FEATURES at resolution RMSD
      //! @param FEATURES the set of features to consider
      //! @param RMSD the rmsd; vectors closer than this rmsd will be considered representative of each other
      //! @param CLUSTER_IDS vector that will store cluster ids (1-offset), used for refinement if m_ClusterIDs was set
      //! @return indices of features that cover the same space as FEATURES at resolution RMSD
      storage::Vector< size_t> GetRepresentativeFeatures
      (
        const linal::MatrixConstInterface< float> &FEATURES,
        const double &RMSD,
        storage::Vector< size_t> &CLUSTER_IDS
      ) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
      //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
      virtual size_t GetNominalSize() const;

      //! @brief refine the id set to be the center-most element of each cluster
      //! @param CLUSTER_IDS vector that will store cluster ids (1-offset), used for refinement if m_ClusterIDs was set
      //! @param FEATURES set of features to consider
      //! @return the updated feature ids
      storage::Vector< size_t> GetClusterCenters
      (
        storage::Vector< size_t> &CLUSTER_IDS,
        const linal::MatrixConstInterface< float> &FEATURES
      ) const;

    }; // class DataSetReducedToClusterCenters

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_DATA_SET_REDUCED_TO_CLUSTER_CENTERS_H_

