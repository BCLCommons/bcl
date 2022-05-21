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

#ifndef BCL_MODEL_RETRIEVE_DATA_SET_BY_FEATURE_H_
#define BCL_MODEL_RETRIEVE_DATA_SET_BY_FEATURE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_retrieve_data_set_base.h"
#include "math/bcl_math_range_set.h"
#include "util/bcl_util_function_interface_serializable.h"
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
    //! @class RetrieveDataSetByFeature
    //! @brief reduces a data set to a series of cluster centers that cover the same space
    //! There may be a smaller data set that covers the same space, which may be obtained using k-means
    //! This method should be faster, and can be used on its own or to provide the initial cluster points
    //! for the k-means algorithm
    //!
    //! @author mendenjl
    //! @see @link example_model_retrieve_data_set_by_feature.cpp @endlink
    //! @date Jul 24, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RetrieveDataSetByFeature :
      public RetrieveDataSetBase,
      public util::FunctionInterfaceSerializable< util::ShPtr< descriptor::Dataset>, util::ShPtr< descriptor::Dataset> >
    {
    private:

    //////////
    // data //
    //////////

      math::RangeSet< float>                 m_FeatureRanges;    //!< ranges of results to accept
      std::string                            m_FeatureName;      //!< Name of the feature
      size_t                                 m_FeatureIndex;     //!< Index of the feature
      util::Implementation< RetrieveDataSetBase> m_Retriever;    //!< Retrieval method for the data set
      bool                                   m_RetrieveFeature;  //!< True if retrieving features (false for results)
      bool                                   m_FeatureInDataset; //!< True if the feature is one that is also in the dataset
      util::Implementation< RetrieveDataSetBase> m_RetrieverHardCopy; //!< Retriever hard copy, needed if the feature is not in the final dataset
      bool                                   m_Invert;           //!< True to invert the comparison; e.g. select all outside the range
      bool                                   m_JustFilter;       //!< True to just filter the results

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_FeaturesInstance;
      static const util::SiPtr< const util::ObjectInterface> s_ResultsInstance;
      static const util::SiPtr< const util::ObjectInterface> s_FeaturesFilterInstance;
      static const util::SiPtr< const util::ObjectInterface> s_ResultsFilterInstance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from parameters
      //! @param RANGES the ranges of features to consider
      //! @param RETRIEVER the dataset retriever to use
      RetrieveDataSetByFeature
      (
        const bool &RETRIEVE_FEATURES = true,
        const bool &JUST_FILTER = false,
        const size_t &FEATURE_INDEX = util::GetUndefined< size_t>(),
        const math::RangeSet< float> &RANGES = math::RangeSet< float>::GetCompleteRange(),
        const util::Implementation< RetrieveDataSetBase> &RETRIEVER = util::Implementation< RetrieveDataSetBase>(),
        const std::string &DESCRIPTOR = std::string()
      );

      //! @brief Clone function
      //! @return pointer to new RetrieveDataSetByFeature
      RetrieveDataSetByFeature *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the feature ranges of interest
      //! @return the feature ranges of interest
      const math::RangeSet< float> &GetFeatureRanges() const
      {
        return m_FeatureRanges;
      }

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

      //! @brief operator() can be used to filter a dataset
      //! @param DATASET the dataset to filter
      //! @return the filtered dataset
      util::ShPtr< descriptor::Dataset> operator()( const util::ShPtr< descriptor::Dataset> &DATASET) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief test whether this retriever can generate sub-ranges of datasets without loading the entire dataset
      //! @return true if this retriever can generate sub-ranges of datasets without loading the entire dataset
      bool SupportsEfficientSubsetLoading() const;

      //! @brief load a range of data from the dataset
      //! @param SUBSET the range of data to load
      //! @param FEATURES_STORAGE where to store features that are loaded, must be large enough to hold the subset without resizing
      //! @param RESULTS_STORAGE where to store the corresponding results, must be large enough to hold the subset without resizing
      //! @param START_FEATURE_NUMBER position to store the first feature in FEATURES_STORAGE
      //! @return # of features actually loaded
      //! Note: Implementations should overload this and SupportsEfficientSubsetLoading together
      size_t GenerateDataSubset
      (
        const math::Range< size_t> &SUBSET,
        linal::MatrixInterface< float> &FEATURES_STORAGE,
        linal::MatrixInterface< float> &RESULTS_STORAGE,
        linal::MatrixInterface< char> &IDS_STORAGE,
        const size_t &START_FEATURE_NUMBER
      );

      //! @brief get the feature indices in the specified range
      //! @param FEATURES the set of features to consider
      //! @return indices of features for which the desired feature is in the specified range set
      storage::Vector< size_t> GetFeaturesInRange
      (
        const linal::MatrixConstInterface< float> &FEATURES,
        const size_t &FEATURE_INDEX
      ) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
      //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
      size_t GetNominalSize() const;

      //! @brief update the m_FeatureIndex
      void UpdateFeatureIndex();

    }; // class RetrieveDataSetByFeature

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_RETRIEVE_DATA_SET_BY_FEATURE_H_

