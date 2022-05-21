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

#ifndef BCL_MODEL_RETRIEVE_DATA_SET_CHUNK_H_
#define BCL_MODEL_RETRIEVE_DATA_SET_CHUNK_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
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
    //! @class RetrieveDataSetChunk
    //! @brief retrieve one or multiple ranges of rows from a given dataset retriever.
    //! @detail this retriever allows for extracting a number of rows from a dataset given a range [0,max] eg. [0,1000].
    //!         It is assumed that there are as many chunks available as there are rows in the dataset
    //!         If the specified range is out of bound then the overhead from [max, range max] will be subtracted.
    //!
    //! @see @link example_model_retrieve_data_set_chunk.cpp @endlink
    //! @author mendenjl
    //! @date Oct 31, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RetrieveDataSetChunk :
      public RetrieveDataSetBase
    {
    private:

    //////////
    // data //
    //////////

      util::Implementation< RetrieveDataSetBase> m_Implementation;

      //! number of partitions / chunks that the dataset should be considered to have
      size_t m_NumberChunks;

      math::RangeSet< size_t> m_ChunkRanges; //!< ranges of chunks to extract from the file

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new RetrieveDataSetChunk
      RetrieveDataSetChunk *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

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

      //! @brief get the number of partitions requested by the user, along with the partition ids
      //! @return the number of partitions requested by the user, along with the partition ids
      storage::Pair< size_t, math::RangeSet< size_t> > GetNumberPartitionsAndIds() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief generate dataset
      //! @return generated dataset
      virtual util::ShPtr< descriptor::Dataset> GenerateDataSet();

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

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

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      virtual io::Serializer GetSerializer() const;

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

      //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
      //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
      virtual size_t GetNominalSize() const;

      //! @brief get the size of the complete dataset without chunking
      //! @return the size of the complete dataset without chunking
      size_t GetTotalDatasetSize() const;

      //! @brief get the sizes of the features and results in the dataset without loading the whole dataset
      //! @return the sizes of the features and results in the dataset without loading the whole dataset
      storage::VectorND< 2, size_t> GetFeatureResultSizes() const;

    }; // class RetrieveDataSetChunk

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_RETRIEVE_DATA_SET_CHUNK_H_
