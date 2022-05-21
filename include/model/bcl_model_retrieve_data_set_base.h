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

#ifndef BCL_MODEL_RETRIEVE_DATA_SET_BASE_H_
#define BCL_MODEL_RETRIEVE_DATA_SET_BASE_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "descriptor/bcl_descriptor.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_interface.h"
#include "descriptor/bcl_descriptor_dataset.h"
#include "linal/bcl_linal_matrix_const_reference.h"
#include "linal/bcl_linal_matrix_reference.h"
#include "math/bcl_math_range_set.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_object_interface.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RetrieveDataSetBase
    //! @brief interface to retrieve data sets independent of source
    //!
    //! @remarks example unnecessary
    //! @author mendenjl, butkiem1
    //! @date Jan 29, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RetrieveDataSetBase :
      public virtual util::SerializableInterface
    {

    private:

    //////////
    // data //
    //////////

      util::ObjectDataLabel m_FeatureCodeLabel; //!< label for the feature code
      util::ObjectDataLabel m_ResultCodeLabel;  //!< label for the result code
      util::ObjectDataLabel m_IdCodeLabel;      //!< label for the id code

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new RetrieveDataSetBase
      virtual RetrieveDataSetBase *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief Set the code / label for the feature (1st part) of the data set
      //! @param CODE the new code
      virtual void SelectFeatures( const util::ObjectDataLabel &CODE);

      //! @brief Set the code / label for the result (2nd part) of the data set
      //! @param CODE the new code
      virtual void SelectResults( const util::ObjectDataLabel &CODE);

      //! @brief Set the code / label for the ids (3rd part) of the data set
      //! @param CODE the new code
      virtual void SelectIds( const util::ObjectDataLabel &CODE);

      //! @brief Convenience function for apps that use implementations of this interface
      //! This function checks loads the feature labels and/or columns from the given filename flag
      //! and asserts if the filename flag was not given but was needed
      //! @param FLAG flag which can contain filename of the features and/or feature indices
      void SelectFeaturesGivenFilenameFlag( const command::FlagInterface &FLAG);

      //! @brief Convenience function for apps that use implementations of this interface
      //! This function checks loads the result labels and/or columns from the given filename flag
      //! and asserts if the filename flag was not given but was needed
      //! @param FLAG flag which can contain filename of the feature labels, if they are needed
      void SelectResultsGivenFilenameFlag( const command::FlagInterface &FLAG);

      //! @brief Convenience function for apps that use implementations of this interface
      //! This function checks loads the ids labels and/or columns from the given filename flag
      //! and asserts if the filename flag was not given but was needed
      //! @param FLAG which can contain filename of the id labels, if they are needed
      void SelectIdsGivenFilenameFlag( const command::FlagInterface &FLAG);

      //! @brief Get the code / label for the feature (1st part) of the data set
      //! @return the code / label for the feature (1st part) of the data set
      const util::ObjectDataLabel &GetFeatureCode() const;

      //! @brief Get the code / label for the result (2nd part) of the data set
      //! @return the code / label for the result (2nd part) of the data set
      const util::ObjectDataLabel &GetResultCode() const;

      //! @brief Get the code / label for the ids (3rd part) of the data set
      //! @return the code / label for the ids (3rd part) of the data set
      const util::ObjectDataLabel &GetIdCode() const;

      //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
      //! @return the code / label for the feature (1st part) of the data set with sizes of each property
      //! the feature code set
      virtual FeatureLabelSet GetFeatureLabelsWithSizes() const = 0;

      //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
      //! @return the code / label for the result (2nd part) of the data set with sizes of each property
      //! the feature code set
      virtual FeatureLabelSet GetResultCodeWithSizes() const = 0;

      //! @brief Get the code / label for the ids of the data set with sizes of each property
      //! @return the code / label for the ids of the data set with sizes of each property
      //! the feature code set
      virtual FeatureLabelSet GetIdCodeWithSizes() const = 0;

      //! @brief get the number of partitions requested by the user, along with the partition ids
      //! @return the number of partitions requested by the user, along with the partition ids
      virtual storage::Pair< size_t, math::RangeSet< size_t> > GetNumberPartitionsAndIds() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief generate dataset
      //! @return generated dataset
      virtual util::ShPtr< descriptor::Dataset> GenerateDataSet() = 0;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief determine the start and end index position given the size, the range id and the number of ranges
      //! @param RANGE_OF_IDS range of range ids
      //! @param NUMBER_RANGES number of total ranges
      //! @param TOTAL_SIZE_OF_DATASET total size of dataset
      //! @return range with start and end position
      static math::Range< size_t> GetStartEndPositionOfRange
      (
        const math::Range< size_t> RANGE_OF_IDS,
        const size_t NUMBER_RANGES,
        const size_t TOTAL_SIZE_OF_DATASET
      );

      //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
      //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
      virtual size_t GetNominalSize() const = 0;

      //! @brief get whether dataset generation requires labels
      //! @return true if dataset generation requires labels
      virtual bool RequiresFeatureLabels() const
      {
        return false;
      }

      //! @brief get whether dataset generation requires result labels
      //! @return true if dataset generation requires result labels
      virtual bool RequiresResultLabels() const
      {
        return false;
      }

      //! @brief test whether this retriever can generate sub-ranges of datasets without loading the entire dataset
      //! @return true if this retriever can generate sub-ranges of datasets without loading the entire dataset
      virtual bool SupportsEfficientSubsetLoading() const
      {
        // by default, dataset retriever do not support efficient subset loading
        return false;
      }

      //! @brief load a range of data from the dataset
      //! @param SUBSET the range of data to load
      //! @param FEATURES_STORAGE where to store features that are loaded
      //! @param RESULTS_STORAGE where to store the corresponding results
      //! @param IDS_STORAGE wehre to put the corresponding ids
      //! @param START_FEATURE_NUMBER position to store the first feature in FEATURES_STORAGE
      //! @return # of features actually loaded
      //! Note: Implementations should overload this and SupportsEfficientSubsetLoading together
      virtual size_t GenerateDataSubset
      (
        const math::Range< size_t> &SUBSET,
        linal::MatrixInterface< float> &FEATURES_STORAGE,
        linal::MatrixInterface< float> &RESULTS_STORAGE,
        linal::MatrixInterface< char> &IDS_STORAGE,
        const size_t &START_FEATURE_NUMBER
      );

    protected:

      //! @brief get the size of the dataset resulting from a selected range of chunks of a dataset with a given size
      //! @param RANGES_OF_IDS range of range ids
      //! @param NUMBER_RANGES number of total ranges
      //! @param TOTAL_SIZE_OF_DATASET total size of dataset
      //! @return the size of the dataset resulting from a selected range of chunks of a dataset with a given size
      static size_t GetSubsetSize
      (
        const math::RangeSet< size_t> RANGES_OF_IDS,
        const size_t NUMBER_RANGES,
        const size_t TOTAL_SIZE_OF_DATASET
      );

    }; // class RetrieveDataSetBase

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_RETRIEVE_DATA_SET_BASE_H_
