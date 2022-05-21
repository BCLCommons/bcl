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

#ifndef BCL_MODEL_RETRIEVE_DATASET_SUBSET_H_
#define BCL_MODEL_RETRIEVE_DATASET_SUBSET_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "io/bcl_io.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_data_set_select_columns.h"
#include "bcl_model_retrieve_data_set_base.h"
#include "math/bcl_math_range_set.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RetrieveDatasetSubset
    //! @brief retrieve data sets from a file
    //!
    //! @see @link example_model_retrieve_dataset_subset.cpp @endlink
    //! @author mendenjl
    //! @date Apr 14, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RetrieveDatasetSubset :
      public RetrieveDataSetBase
    {
    private:

    //////////
    // data //
    //////////

      std::string m_Filename;                //!< file from which the data set can be retrieved
      size_t m_NumberChunks;                 //!< split the file into this many parts, to retrieve the chunks specified
      math::RangeSet< size_t> m_ChunkRanges; //!< ranges of chunks to extract from the file
      DataSetSelectColumns    m_Features;    //!< Feature columns that will be selected
      DataSetSelectColumns    m_Results;     //!< Result columns that will be selected
      DataSetSelectColumns    m_Ids;         //!< Id columns that will be selected

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param FILENAME the file to retrieve the data set from
      //! @param NUMBER_CHUNKS # of chunks the file should be split into (conceptually), used to divide the file into disparate datasets
      //! @param CHUNKS chunks to load
      RetrieveDatasetSubset
      (
        const std::string &FILENAME = std::string(),
        const size_t &NUMBER_CHUNKS = 1,
        const math::RangeSet< size_t> &CHUNKS = math::RangeSet< size_t>( math::Range< size_t>( 0, 0))
      );

      //! @brief Clone function
      //! @return pointer to new RetrieveDatasetSubset
      RetrieveDatasetSubset *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

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

    ////////////////
    // operations //
    ////////////////

      //! @brief generate dataset
      //! @return generated dataset
      util::ShPtr< descriptor::Dataset> GenerateDataSet();

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief store the data from a dataset into a master dataset
      //! @param FILENAME file to store the dataset in
      //! @param DATA_SET_RETRIEVER object used to retrieve the dataset
      //! @param BLOCK_SIZE_MB maximum MB of data to generate per thread before writing to disk
      //! @param ALLOW_INCOMPLETE_RECORDS true to allow incomplete results to be stored; by default, all results
      //!        must be given for each record
      static void StoreMasterDataset
      (
        const std::string &FILENAME,
        RetrieveDataSetBase &DATA_SET_RETRIEVER,
        const double &BLOCK_SIZE_MB = double( 32.0),
        const bool &ALLOW_INCOMPLETE_RECORDS = false
      );

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      virtual io::Serializer GetSerializer() const;

      //! @brief get the size of the complete dataset without chunking
      //! @return the size of the complete dataset without chunking
      size_t GetTotalDatasetSize() const;

      //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
      //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
      size_t GetNominalSize() const;

      //! @brief test whether this retriever can generate sub-ranges of datasets without loading the entire dataset
      //! @return true if this retriever can generate sub-ranges of datasets without loading the entire dataset
      bool SupportsEfficientSubsetLoading() const
      {
        return true;
      }

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

      //! @brief get the number of partitions requested by the user, along with the partition ids
      //! @return the number of partitions requested by the user, along with the partition ids
      storage::Pair< size_t, math::RangeSet< size_t> > GetNumberPartitionsAndIds() const;

    private:

      //! @brief read a feature label set written to a binary file
      static FeatureLabelSet ReadFeatureLabelSetBinary( std::ifstream &READ);

      //! @brief read a feature label set written to a binary file
      static void WriteFeatureLabelSetBinary( std::ostream &WRITE, const FeatureLabelSet &LABELS);

      //! @brief read the feature and result labels and create DataSetSelectColumns from them
      //! @param INPUT the input stream to use; should not be open initially
      //!        after this function is called, it will be on the first line of the dataset
      //! @return the # of features in this file, which is determined here
      size_t GetTotalSizeAndParseLabels( std::ifstream &INPUT);

    }; // class RetrieveDatasetSubset

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_RETRIEVE_DATASET_SUBSET_H_
