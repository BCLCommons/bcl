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

// include forward header of this class
#include "model/bcl_model_retrieve_data_set_by_feature.h"

// include other forward headers - sorted alphabetically

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
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetByFeature::s_FeaturesInstance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new RetrieveDataSetByFeature( true))
    );
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetByFeature::s_ResultsInstance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new RetrieveDataSetByFeature( false))
    );
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetByFeature::s_FeaturesFilterInstance
    (
      util::Enumerated< util::FunctionInterfaceSerializable< util::ShPtr< descriptor::Dataset>, util::ShPtr< descriptor::Dataset> > >::AddInstance
      (
        new RetrieveDataSetByFeature( true, true)
      )
    );
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetByFeature::s_ResultsFilterInstance
    (
      util::Enumerated< util::FunctionInterfaceSerializable< util::ShPtr< descriptor::Dataset>, util::ShPtr< descriptor::Dataset> > >::AddInstance
      (
        new RetrieveDataSetByFeature( false, true)
      )
    );

    //! @brief constructor from parameters
    //! @param RANGES the ranges of features to consider
    //! @param RETRIEVER the dataset retriever to use
    RetrieveDataSetByFeature::RetrieveDataSetByFeature
    (
      const bool &RETRIEVE_FEATURES,
      const bool &JUST_FILTER,
      const size_t &FEATURE_INDEX,
      const math::RangeSet< float> &RANGES,
      const util::Implementation< RetrieveDataSetBase> &RETRIEVER,
      const std::string &DESCRIPTOR
    ) :
      m_FeatureRanges( RANGES),
      m_FeatureName( DESCRIPTOR),
      m_FeatureIndex( FEATURE_INDEX),
      m_Retriever( RETRIEVER),
      m_RetrieveFeature( RETRIEVE_FEATURES),
      m_FeatureInDataset( false),
      m_Invert( false),
      m_JustFilter( JUST_FILTER)
    {
    }

    //! @brief Clone function
    //! @return pointer to new RetrieveDataSetByFeature
    RetrieveDataSetByFeature *RetrieveDataSetByFeature::Clone() const
    {
      return new RetrieveDataSetByFeature( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RetrieveDataSetByFeature::GetClassIdentifier() const
    {
      return GetStaticClassName< RetrieveDataSetByFeature>();
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RetrieveDataSetByFeature::GetAlias() const
    {
      static const std::string s_feature_name( "FeatureRange"), s_results_name( "ResultRange");
      return m_RetrieveFeature ? s_feature_name : s_results_name;
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    void RetrieveDataSetByFeature::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectFeatures( CODE);
      m_Retriever->SelectFeatures( CODE);
      if( m_RetrieveFeature)
      {
        m_FeatureIndex = util::GetUndefined< size_t>();
      }
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @return the code / label for the result (2nd part) of the data set
    void RetrieveDataSetByFeature::SelectResults( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectResults( CODE);
      m_Retriever->SelectResults( CODE);
      if( !m_RetrieveFeature)
      {
        m_FeatureIndex = util::GetUndefined< size_t>();
      }
    }

    //! @brief Set which id columns to retrieve
    //! @param CODE the id column names to retrieve
    void RetrieveDataSetByFeature::SelectIds( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectIds( CODE);
      m_Retriever->SelectIds( CODE);
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetByFeature::GetFeatureLabelsWithSizes() const
    {
      return m_Retriever->GetFeatureLabelsWithSizes();
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetByFeature::GetResultCodeWithSizes() const
    {
      return m_Retriever->GetResultCodeWithSizes();
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetByFeature::GetIdCodeWithSizes() const
    {
      return m_Retriever->GetIdCodeWithSizes();
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > RetrieveDataSetByFeature::GetNumberPartitionsAndIds() const
    {
      return m_Retriever->GetNumberPartitionsAndIds();
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief generate dataset, reduced to the desired # of cluster centers
    //! @return generated dataset
    util::ShPtr< descriptor::Dataset>
      RetrieveDataSetByFeature::GenerateDataSet()
    {
      // get the base data set
      util::ShPtr< descriptor::Dataset> dataset( m_Retriever->GenerateDataSet());

      // update the feature index if necessary
      UpdateFeatureIndex();

      storage::Vector< size_t> features_to_keep;
      if( m_FeatureInDataset)
      {
        // select the desired features
        features_to_keep =
          GetFeaturesInRange( m_RetrieveFeature ? dataset->GetFeaturesReference() : dataset->GetResultsReference(), m_FeatureIndex);
      }
      else
      {
        util::ShPtr< descriptor::Dataset> dataset_tmp( m_RetrieverHardCopy->GenerateDataSet());
        BCL_Assert
        (
          dataset_tmp->GetSize() == dataset->GetSize(),
          "Dataset filtering failed because dataset generated without the filter descriptor has a different size than "
          "the dataset generated with the descriptor"
        );
        features_to_keep =
          GetFeaturesInRange( m_RetrieveFeature ? dataset_tmp->GetFeaturesReference() : dataset_tmp->GetResultsReference(), m_FeatureIndex);
      }
      dataset->KeepRows( features_to_keep);

      // return
      return dataset;
    }

    //! @brief operator() can be used to filter a dataset
    //! @param DATASET the dataset to filter
    //! @return the filtered dataset
    util::ShPtr< descriptor::Dataset> RetrieveDataSetByFeature::operator()( const util::ShPtr< descriptor::Dataset> &DATASET) const
    {
      // update the feature index if necessary
      FeatureLabelSet feature_label_set
      (
        m_RetrieveFeature ? *DATASET->GetFeaturesPtr()->GetFeatureLabelSet() : *DATASET->GetResultsPtr()->GetFeatureLabelSet()
      );
      feature_label_set = feature_label_set.SplitFeatureLabelSet( true);
      size_t feature_position( feature_label_set.GetMemberLabels().Find( m_FeatureName));
      BCL_Assert
      (
        feature_position < feature_label_set.GetMemberLabels().GetSize(),
        "Filter descriptor must be present in the dataset for pass-through filtering to work (e.g. model::TrainMCCV)"
      );
      BCL_MessageStd( "Selecting from property index: " + util::Format()( feature_position));
      return
        DATASET->GetRows
        (
          GetFeaturesInRange
          (
            m_RetrieveFeature ? DATASET->GetFeaturesReference() : DATASET->GetResultsReference(),
            feature_position
          )
        );
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get the feature indices in the specified range
    //! @param FEATURES the set of features to consider
    //! @return indices of features for which the desired feature is in the specified range set
    storage::Vector< size_t> RetrieveDataSetByFeature::GetFeaturesInRange
    (
      const linal::MatrixConstInterface< float> &FEATURES,
      const size_t &FEATURE_INDEX
    ) const
    {
      storage::Vector< size_t> chosen_features;
      chosen_features.AllocateMemory( FEATURES.GetNumberRows());
      BCL_Assert( FEATURE_INDEX < FEATURES.GetNumberCols(), "Requested feature not in this dataset");
      for( size_t row( 0), n_rows( FEATURES.GetNumberRows()); row < n_rows; ++row)
      {
        if( m_FeatureRanges.IsWithin( FEATURES[ row][ FEATURE_INDEX]) != m_Invert)
        {
          chosen_features.PushBack( row);
        }
      }
      return chosen_features;
    }

    //! @brief test whether this retriever can generate sub-ranges of datasets without loading the entire dataset
    //! @return true if this retriever can generate sub-ranges of datasets without loading the entire dataset
    bool RetrieveDataSetByFeature::SupportsEfficientSubsetLoading() const
    {
      return m_Retriever->SupportsEfficientSubsetLoading();
    }

    //! @brief load a range of data from the dataset
    //! @param SUBSET the range of data to load
    //! @param FEATURES_STORAGE where to store features that are loaded, must be large enough to hold the subset without resizing
    //! @param RESULTS_STORAGE where to store the corresponding results, must be large enough to hold the subset without resizing
    //! @param START_FEATURE_NUMBER position to store the first feature in FEATURES_STORAGE
    //! @return # of features actually loaded
    //! Note: Implementations should overload this and SupportsEfficientSubsetLoading together
    size_t RetrieveDataSetByFeature::GenerateDataSubset
    (
      const math::Range< size_t> &SUBSET,
      linal::MatrixInterface< float> &FEATURES_STORAGE,
      linal::MatrixInterface< float> &RESULTS_STORAGE,
      linal::MatrixInterface< char> &IDS_STORAGE,
      const size_t &START_FEATURE_NUMBER
    )
    {
      // get the base data set
      const size_t n_rows( SUBSET.GetWidth() + 1);
      linal::Matrix< float> features( n_rows, FEATURES_STORAGE.GetNumberCols(), util::GetUndefined< float>());
      linal::Matrix< float> results( n_rows, RESULTS_STORAGE.GetNumberCols(), util::GetUndefined< float>());
      linal::Matrix< char> ids( n_rows, IDS_STORAGE.GetNumberCols(), ' ');
      size_t n_created( m_Retriever->GenerateDataSubset( SUBSET, features, results, ids, size_t( 0)));

      // update the feature index if necessary
      UpdateFeatureIndex();

      if( !n_created)
      {
        BCL_MessageStd( "No features kept in this block pre-filter");
        return 0;
      }
      features.ShrinkRows( n_created);
      results.ShrinkRows( n_created);
      ids.ShrinkRows( n_created);

      storage::Vector< size_t> features_to_keep;
      if( m_FeatureInDataset)
      {
        // select the desired features
        linal::MatrixConstReference< float> reference( m_RetrieveFeature ? features : results);
        features_to_keep = GetFeaturesInRange( reference, m_FeatureIndex);
      }
      else
      {
        linal::Matrix< float> filter( n_rows, 1, util::GetUndefined< float>()), other( filter);
        linal::Matrix< char> ids_tmp;
        size_t n_created_tmp( 0);
        if( m_RetrieveFeature)
        {
          n_created_tmp = m_RetrieverHardCopy->GenerateDataSubset( SUBSET, filter, other, ids_tmp, size_t( 0));
        }
        else
        {
          n_created_tmp = m_RetrieverHardCopy->GenerateDataSubset( SUBSET, other, filter, ids, size_t( 0));
        }
        BCL_Assert
        (
          n_created_tmp == n_created,
          "Dataset filtering failed because dataset generated without the filter descriptor has a different size than "
          "the dataset generated with the descriptor"
        );
        features_to_keep = GetFeaturesInRange( filter, m_FeatureIndex);
      }
      features.KeepRows( features_to_keep);
      results.KeepRows( features_to_keep);
      ids.KeepRows( features_to_keep);

      BCL_Assert
      (
        features_to_keep.GetSize() == features.GetNumberRows(),
        "All features in the given range should have been kept"
      );

      if( features_to_keep.IsEmpty())
      {
        BCL_MessageStd( "No features kept in this block, post-filter");
        return 0;
      }

      // copy the code vectors into the matrices
      std::copy( features.Begin(), features.End(), FEATURES_STORAGE[ START_FEATURE_NUMBER]);
      // copy the result
      std::copy( results.Begin(), results.End(), RESULTS_STORAGE[ START_FEATURE_NUMBER]);
      // copy the ids
      std::copy( ids.Begin(), ids.End(), IDS_STORAGE[ START_FEATURE_NUMBER]);

      // return
      return features_to_keep.GetSize();
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RetrieveDataSetByFeature::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription
      (
        m_RetrieveFeature
        ? "Select features (and corresponding results) within a particular range set"
        : "Select results (and corresponding features) within a particular range set"
      );
      if( !m_JustFilter)
      {
        member_data.AddInitializer
        (
          "dataset",
          "dataset retriever to call to get the entire data set",
          io::Serialization::GetAgent( &m_Retriever)
        );
      }
      member_data.AddInitializer
      (
        "range",
        "ranges of values to load, e.g. range=\"[ 0, 5) (7,10)\"",
        io::Serialization::GetAgent( &m_FeatureRanges),
        math::RangeSet< float>::GetCompleteRange().AsString()
      );
      member_data.AddInitializer
      (
        "invert",
        "true to load all values except those specified in range",
        io::Serialization::GetAgent( &m_Invert),
        "False"
      );

      member_data.AddInitializer
      (
        m_RetrieveFeature ? "feature" : "result",
        "Name of " + std::string( m_RetrieveFeature ? "feature" : "result") + " to consider",
        io::Serialization::GetAgent( &m_FeatureName)
      );

      return member_data;
    }

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t RetrieveDataSetByFeature::GetNominalSize() const
    {
      return m_Retriever->GetNominalSize();
    }

    //! @brief update the m_FeatureIndex
    void RetrieveDataSetByFeature::UpdateFeatureIndex()
    {
      // if the feature index is undefined, determine it from the feature string
      if( !util::IsDefined( m_FeatureIndex))
      {
        m_FeatureInDataset = false;
        storage::Vector< size_t> property_indices;
        FeatureLabelSet feature_label_set
        (
          m_RetrieveFeature ? m_Retriever->GetFeatureLabelsWithSizes() : m_Retriever->GetResultCodeWithSizes()
        );
        feature_label_set = feature_label_set.SplitFeatureLabelSet( true);
        size_t feature_position( feature_label_set.GetMemberLabels().Find( m_FeatureName));

        if( feature_position >= feature_label_set.GetMemberLabels().GetSize())
        {
          m_RetrieverHardCopy = m_Retriever.HardCopy();
          if( m_RetrieveFeature)
          {
            BCL_Assert
            (
              !m_Retriever->RequiresFeatureLabels(),
              "Dataset filtering requires that either the dataset have already been generated, or that the filter descriptor "
              "is included in the dataset"
            );
            m_RetrieverHardCopy->SelectFeatures( util::ObjectDataLabel( "", "Combine", m_FeatureName));
            m_RetrieverHardCopy->SelectResults
            (
              util::ObjectDataLabel
              (
                "",
                "Combine",
                m_RetrieverHardCopy->GetResultCodeWithSizes().SplitFeatureLabelSet().GetMemberLabels().FirstElement()
              )
            );
          }
          else
          {
            BCL_Assert
            (
              !m_Retriever->RequiresResultLabels(),
              "Dataset filtering requires that either the dataset have already been generated, or that the filter descriptor "
              "is included in the dataset"
            );
            m_RetrieverHardCopy->SelectResults( util::ObjectDataLabel( "", "Combine", m_FeatureName));
            m_RetrieverHardCopy->SelectFeatures
            (
              util::ObjectDataLabel
              (
                "",
                "Combine",
                m_RetrieverHardCopy->GetFeatureLabelsWithSizes().SplitFeatureLabelSet().GetMemberLabels().FirstElement()
              )
            );
          }
          m_RetrieverHardCopy->SelectIds( util::ObjectDataLabel());
          m_FeatureIndex = 0;
          m_FeatureInDataset = false;
          BCL_MessageStd
          (
            "Filter descriptor was not in dataset; finding it using "
            + m_RetrieverHardCopy->GetFeatureCode().ToString()
            + " features and " + m_RetrieverHardCopy->GetResultCode().ToString() + " results"
          );
        }
        else
        {
          BCL_MessageStd( "Selecting from property index: " + util::Format()( feature_position));
          m_FeatureIndex = feature_position;
          m_FeatureInDataset = true;
        }
      }
    }

  } // namespace model
} // namespace bcl
