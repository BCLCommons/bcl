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
#include "model/bcl_model_retrieve_data_set_chunk.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetChunk::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new RetrieveDataSetChunk())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorAtomFurthest which is a copy of this
    RetrieveDataSetChunk *RetrieveDataSetChunk::Clone() const
    {
      return new RetrieveDataSetChunk( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RetrieveDataSetChunk::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RetrieveDataSetChunk::GetAlias() const
    {
      static const std::string s_Name( "Chunks");
      return s_Name;
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    //! @return the code / label for the feature (1st part) of the data set
    void RetrieveDataSetChunk::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectFeatures( CODE);
      m_Implementation->SelectFeatures( CODE);
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @return the code / label for the result (2nd part) of the data set
    void RetrieveDataSetChunk::SelectResults( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectResults( CODE);
      m_Implementation->SelectResults( CODE);
    }

    //! @brief Set which id columns to retrieve
    //! @param CODE the id column names to retrieve
    void RetrieveDataSetChunk::SelectIds( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectIds( CODE);
      m_Implementation->SelectIds( CODE);
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetChunk::GetFeatureLabelsWithSizes() const
    {
      return m_Implementation->GetFeatureLabelsWithSizes();
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetChunk::GetResultCodeWithSizes() const
    {
      return m_Implementation->GetResultCodeWithSizes();
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetChunk::GetIdCodeWithSizes() const
    {
      return m_Implementation->GetIdCodeWithSizes();
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > RetrieveDataSetChunk::GetNumberPartitionsAndIds() const
    {
      return storage::Pair< size_t, math::RangeSet< size_t> >( m_NumberChunks, m_ChunkRanges);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate dataset
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset>
    RetrieveDataSetChunk::GenerateDataSet()
    {
      const size_t total_number_feature_results( m_Implementation->GetNominalSize());

      math::RangeSet< size_t> fixed_chunk_ranges( m_ChunkRanges);
      // if more chunks were specified than there are points in the dataset, reduce the number to the max allowed for
      // this dataset
      const size_t number_chunks( std::min( m_NumberChunks, total_number_feature_results));
      // if range [0,0] is given assume user wants to use entire data set
      if( fixed_chunk_ranges.IsEmpty())
      {
        fixed_chunk_ranges += math::Range< size_t>( 0, number_chunks);
      }

      // eliminate range above total_number_feature_results
      if( fixed_chunk_ranges.GetMax() >= number_chunks)
      {
        fixed_chunk_ranges -= math::Range< size_t>( number_chunks, fixed_chunk_ranges.GetMax());
      }

      // determine the # of features/results in the set of data that will be loaded
      const size_t number_feature_results
      (
        RetrieveDataSetBase::GetSubsetSize( fixed_chunk_ranges, number_chunks, total_number_feature_results)
      );

      // initialize a new data set
      util::ShPtr< descriptor::Dataset> complete
      (
        new descriptor::Dataset
        (
          number_feature_results,
          GetFeatureLabelsWithSizes(),
          GetResultCodeWithSizes(),
          GetIdCodeWithSizes()
        )
      );

      linal::MatrixReference< float> features( complete->GetFeaturesReference());
      linal::MatrixReference< float> results( complete->GetResultsReference());
      linal::MatrixReference< char> ids( complete->GetIdsReference());

      // set the vector up with the given ranges
      size_t features_so_far( 0);
      for
      (
        storage::Set< math::Range< size_t> >::const_iterator
          itr( fixed_chunk_ranges.GetRanges().Begin()), itr_end( fixed_chunk_ranges.GetRanges().End());
        itr != itr_end;
        ++itr
      )
      {
        math::Range< size_t> start_end_chunk
        (
          RetrieveDataSetBase::GetStartEndPositionOfRange( *itr, number_chunks, total_number_feature_results)
        );
        features_so_far +=
          m_Implementation->GenerateDataSubset
          (
            start_end_chunk.StandardizeRange(),
            features,
            results,
            ids,
            features_so_far
          );
      }

      complete->ShrinkRows( features_so_far);

      return complete;
    } // GenerateDataSet

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RetrieveDataSetChunk::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Retrieves chunks (approximately equal-sized partitions) of a specific dataset");
      parameters.AddInitializer
      (
        "number chunks",
        "number of equal sized chunks (approximately equal-sized partitions) to consider",
        io::Serialization::GetAgentWithMin( &m_NumberChunks, 1)
      );
      parameters.AddInitializer
      (
        "chunks",
        "ranges of chunks to load, e.g. chunks=\"[ 0, 5) (7,10)\"",
        io::Serialization::GetAgent( &m_ChunkRanges)
      );
      parameters.AddInitializer
      (
        "dataset",
        "dataset retriever to call to get the entire data set",
        io::Serialization::GetAgent( &m_Implementation)
      );

      return parameters;
    } // GetParameters

    //! @brief test whether this retriever can generate sub-ranges of datasets without loading the entire dataset
    //! @return true if this retriever can generate sub-ranges of datasets without loading the entire dataset
    bool RetrieveDataSetChunk::SupportsEfficientSubsetLoading() const
    {
      return m_Implementation->SupportsEfficientSubsetLoading();
    }

    //! @brief load a range of data from the dataset
    //! @param SUBSET the range of data to load
    //! @param FEATURES_STORAGE where to store features that are loaded, must be large enough to hold the subset without resizing
    //! @param RESULTS_STORAGE where to store the corresponding results, must be large enough to hold the subset without resizing
    //! @param START_FEATURE_NUMBER position to store the first feature in FEATURES_STORAGE
    //! @return # of features actually loaded
    //! Note: Implementations should overload this and SupportsEfficientSubsetLoading together
    size_t RetrieveDataSetChunk::GenerateDataSubset
    (
      const math::Range< size_t> &SUBSET,
      linal::MatrixInterface< float> &FEATURES_STORAGE,
      linal::MatrixInterface< float> &RESULTS_STORAGE,
      linal::MatrixInterface< char> &IDS_STORAGE,
      const size_t &START_FEATURE_NUMBER
    )
    {
      // GenerateDataset in these wrapper functions just translates the external view of the outside dataset retriever
      // into the internal view of the internally held dataset
      const size_t total_number_feature_results( m_Implementation->GetNominalSize());

      math::RangeSet< size_t> fixed_chunk_ranges( m_ChunkRanges);
      // if more chunks were specified than there are points in the dataset, reduce the number to the max allowed for
      // this dataset
      const size_t number_chunks( std::min( m_NumberChunks, total_number_feature_results));
      // if range [0,0] is given assume user wants to use entire data set
      if( fixed_chunk_ranges.IsEmpty())
      {
        fixed_chunk_ranges += math::Range< size_t>( 0, number_chunks);
      }

      // eliminate range above total_number_feature_results
      if( fixed_chunk_ranges.GetMax() >= number_chunks)
      {
        fixed_chunk_ranges -= math::Range< size_t>( number_chunks, fixed_chunk_ranges.GetMax());
      }

      // standardize the range, so the left side of it has a closed border and the right has an open border
      const math::Range< size_t> standard_subset( SUBSET.StandardizeRange());

      // ensure that the end position is within the storage matrix
      BCL_Assert
      (
        START_FEATURE_NUMBER + standard_subset.GetWidth() <= FEATURES_STORAGE.GetNumberRows(),
        "Storage was too small to contain entire subset"
      );

      // determine the # of features/results in the set of data that will be loaded
      const size_t number_feature_results
      (
        RetrieveDataSetBase::GetSubsetSize( fixed_chunk_ranges, number_chunks, total_number_feature_results)
      );
      math::Range< size_t> subset_fixed( standard_subset);
      if( standard_subset.GetMax() > number_feature_results)
      {
        subset_fixed.SetMax( number_feature_results);
      }

      // get the actual ranges that will be taken from the dataset
      math::RangeSet< size_t> ranges;

      for
      (
        storage::Set< math::Range< size_t> >::const_iterator
          itr_range( fixed_chunk_ranges.GetRanges().Begin()), itr_range_end( fixed_chunk_ranges.GetRanges().End());
        itr_range != itr_range_end;
        ++itr_range
      )
      {
        ranges += RetrieveDataSetBase::GetStartEndPositionOfRange( *itr_range, number_chunks, total_number_feature_results);
      }
      ranges = ranges.GetMappedSubset( subset_fixed);

      // keep track of how many rows have been loaded
      size_t features_so_far( START_FEATURE_NUMBER);

      // walk over the lines desired, loading only those columns and results that we want
      for
      (
        storage::Set< math::Range< size_t> >::const_iterator
          itr_range( ranges.GetRanges().Begin()), itr_range_end( ranges.GetRanges().End());
        itr_range != itr_range_end;
        ++itr_range
      )
      {
        features_so_far +=
          m_Implementation->GenerateDataSubset
          (
            *itr_range,
            FEATURES_STORAGE,
            RESULTS_STORAGE,
            IDS_STORAGE,
            features_so_far
          );
      }

      return features_so_far - START_FEATURE_NUMBER;
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool RetrieveDataSetChunk::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      return true;
    }

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t RetrieveDataSetChunk::GetNominalSize() const
    {
      const size_t total_number_feature_results( m_Implementation->GetNominalSize());
      // determine the # of features/results in the set of data that will be loaded
      const size_t number_feature_results
      (
        RetrieveDataSetBase::GetSubsetSize( m_ChunkRanges, m_NumberChunks, total_number_feature_results)
      );

      return number_feature_results;
    }

  } // namespace model
} // namespace bcl
