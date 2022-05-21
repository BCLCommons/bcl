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
#include "model/bcl_model_retrieve_data_set_rows.h"

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
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetRows::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new RetrieveDataSetRows())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorAtomFurthest which is a copy of this
    RetrieveDataSetRows *RetrieveDataSetRows::Clone() const
    {
      return new RetrieveDataSetRows( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &RetrieveDataSetRows::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RetrieveDataSetRows::GetAlias() const
    {
      static const std::string s_Name( "Rows");
      return s_Name;
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    //! @return the code / label for the feature (1st part) of the data set
    void RetrieveDataSetRows::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectFeatures( CODE);
      m_Implementation->SelectFeatures( CODE);
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @return the code / label for the result (2nd part) of the data set
    void RetrieveDataSetRows::SelectResults( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectResults( CODE);
      m_Implementation->SelectResults( CODE);
    }

    //! @brief Set which id columns to retrieve
    //! @param CODE the id column names to retrieve
    void RetrieveDataSetRows::SelectIds( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectIds( CODE);
      m_Implementation->SelectIds( CODE);
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetRows::GetFeatureLabelsWithSizes() const
    {
      return m_Implementation->GetFeatureLabelsWithSizes();
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetRows::GetResultCodeWithSizes() const
    {
      return m_Implementation->GetResultCodeWithSizes();
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetRows::GetIdCodeWithSizes() const
    {
      return m_Implementation->GetIdCodeWithSizes();
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > RetrieveDataSetRows::GetNumberPartitionsAndIds() const
    {
      return m_Implementation->GetNumberPartitionsAndIds();
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief generate dataset
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset>
    RetrieveDataSetRows::GenerateDataSet()
    {
      const size_t total_number_feature_results( m_Implementation->GetNominalSize());

      // determine the # of features/results in the set of data that will be loaded
      const size_t number_feature_results
      (
        RetrieveDataSetBase::GetSubsetSize( m_ChunkRanges, total_number_feature_results, total_number_feature_results)
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
          itr( m_ChunkRanges.GetRanges().Begin()), itr_end( m_ChunkRanges.GetRanges().End());
        itr != itr_end;
        ++itr
      )
      {
        features_so_far +=
          m_Implementation->GenerateDataSubset( itr->StandardizeRange(), features, results, ids, features_so_far);
      }

      complete->ShrinkRows( features_so_far);

      return complete;
    } // GenerateDataSet

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RetrieveDataSetRows::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription( "Retrieves row ranges of a specific dataset");
      parameters.AddInitializer
      (
        "rows",
        "ranges of rows to load, e.g. rows=\"[ 0, 5) (7,10)\"",
        io::Serialization::GetAgent( &m_ChunkRanges),
        "\"[0]\""
      );
      parameters.AddInitializer
      (
        "dataset",
        "dataset retriever to call to get the entire data set",
        io::Serialization::GetAgent( &m_Implementation)
      );

      return parameters;
    } // GetParameters

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool RetrieveDataSetRows::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      const size_t total_number_feature_results( m_Implementation->GetNominalSize());

      // if range [0,0] is given assume user wants to use entire data set
      if( m_ChunkRanges.IsEmpty())
      {
        m_ChunkRanges += math::Range< size_t>( 0, total_number_feature_results);
      }

      // eliminate range above total_number_feature_results
      if( m_ChunkRanges.GetMax() >= total_number_feature_results)
      {
        m_ChunkRanges -= math::Range< size_t>( total_number_feature_results, m_ChunkRanges.GetMax());
      }
      return true;
    }

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t RetrieveDataSetRows::GetNominalSize() const
    {
      const size_t total_number_feature_results( m_Implementation->GetNominalSize());
      // determine the # of features/results in the set of data that will be loaded
      const size_t number_feature_results
      (
        RetrieveDataSetBase::GetSubsetSize( m_ChunkRanges, total_number_feature_results, total_number_feature_results)
      );

      return number_feature_results;
    }

  } // namespace model
} // namespace bcl
