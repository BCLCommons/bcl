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
#include "model/bcl_model_retrieve_data_set_by_id.h"

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
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetById::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new RetrieveDataSetById())
    );
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetById::s_FilterInstance
    (
      util::Enumerated< util::FunctionInterfaceSerializable< util::ShPtr< descriptor::Dataset>, util::ShPtr< descriptor::Dataset> > >::AddInstance
      (
        new RetrieveDataSetById( storage::Set< std::string>(), util::Implementation< RetrieveDataSetBase>(), true)
      )
    );

    //! @brief constructor from parameters
    //! @param IDS the ids to consider
    //! @param RETRIEVER the dataset retriever to use
    //! @param JUST_FILTER set this to have this object run in pure-filter mode, in which case the dataset will not be
    //!        requested from the command line
    RetrieveDataSetById::RetrieveDataSetById
    (
      const storage::Set< std::string> &IDS,
      const util::Implementation< RetrieveDataSetBase> &RETRIEVER,
      const bool &JUST_FILTER
    ) :
      m_IdName(),
      m_Ids( IDS),
      m_Retriever( RETRIEVER),
      m_Invert( false),
      m_JustFilter( JUST_FILTER)
    {
    }

    //! @brief Clone function
    //! @return pointer to new RetrieveDataSetById
    RetrieveDataSetById *RetrieveDataSetById::Clone() const
    {
      return new RetrieveDataSetById( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RetrieveDataSetById::GetClassIdentifier() const
    {
      return GetStaticClassName< RetrieveDataSetById>();
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RetrieveDataSetById::GetAlias() const
    {
      static const std::string s_name( "SelectIDs");
      return s_name;
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    //! @return the code / label for the feature (1st part) of the data set
    void RetrieveDataSetById::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectFeatures( CODE);
      m_Retriever->SelectFeatures( CODE);
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @return the code / label for the result (2nd part) of the data set
    void RetrieveDataSetById::SelectResults( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectResults( CODE);
      m_Retriever->SelectResults( CODE);
    }

    //! @brief Set which id columns to retrieve
    //! @param CODE the id column names to retrieve
    void RetrieveDataSetById::SelectIds( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectIds( CODE);
      m_Retriever->SelectIds( CODE);
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetById::GetFeatureLabelsWithSizes() const
    {
      return m_Retriever->GetFeatureLabelsWithSizes();
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetById::GetResultCodeWithSizes() const
    {
      return m_Retriever->GetResultCodeWithSizes();
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetById::GetIdCodeWithSizes() const
    {
      return m_Retriever->GetIdCodeWithSizes();
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > RetrieveDataSetById::GetNumberPartitionsAndIds() const
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
      RetrieveDataSetById::GenerateDataSet()
    {
      UpdateIDsIndex();

      // get the base data set
      util::ShPtr< descriptor::Dataset> dataset( m_Retriever->GenerateDataSet());

      if( m_IdIndices.IsEmpty())
      {
        for( size_t i( 0), n_ids( dataset->GetIdSize()); i < n_ids; ++i)
        {
          m_IdIndices.PushBack( i);
        }
      }

      // filter by the desired ids
      dataset->KeepRows( GetDesiredIDs( dataset->GetIdsReference(), m_IdIndices));

      // return
      return dataset;
    }

    //! @brief operator() can be used to filter a dataset
    //! @param DATASET the dataset to filter
    //! @return the filtered dataset
    util::ShPtr< descriptor::Dataset> RetrieveDataSetById::operator()( const util::ShPtr< descriptor::Dataset> &DATASET) const
    {
      storage::Vector< size_t> id_indices;
      const size_t expected_size( m_Ids.Begin()->size());
      if
      (
        m_IdName.empty()
        || !DATASET->GetIdsPtr()->GetFeatureLabelSet().IsDefined()
        || !DATASET->GetIdsPtr()->GetFeatureLabelSet()->GetSize()
        || DATASET->GetIdsPtr()->GetFeatureSize() == expected_size
      )
      {
        id_indices = storage::CreateIndexVector( expected_size);
      }
      else
      {
        id_indices = DATASET->GetIdsPtr()->GetFeatureLabelSet()->GetPropertyIndices( util::ObjectDataLabel( m_IdName));
      }
      return DATASET->GetRows( GetDesiredIDs( DATASET->GetIdsReference(), id_indices));
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get indices of the desired id rows
    //! @param IDS the matrix of ids
    //! @return indices of the desired ids
    storage::Vector< size_t> RetrieveDataSetById::GetDesiredIDs
    (
      const linal::MatrixConstInterface< char> &IDS,
      const storage::Vector< size_t> &ID_INDICES
    ) const
    {
      storage::Vector< size_t> chosen_results;
      chosen_results.AllocateMemory( IDS.GetNumberRows());
      BCL_Assert( ID_INDICES.LastElement() < IDS.GetNumberCols(), "Out of bounds id column!");
      const size_t n_id_cols( ID_INDICES.GetSize());
      std::string tmp( ID_INDICES.GetSize(), ' ');
      for( size_t row( 0), n_rows( IDS.GetNumberRows()); row < n_rows; ++row)
      {
        // copy the desired character into the temp string
        const char *row_ptr( IDS[ row]);
        for( size_t col( 0); col < n_id_cols; ++col)
        {
          tmp[ col] = row_ptr[ ID_INDICES( col)];
        }
        // check whether the id was selected by the user
        if( m_Ids.Contains( util::TrimString( tmp)) != m_Invert)
        {
          chosen_results.PushBack( row);
        }
      }
      return chosen_results;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RetrieveDataSetById::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription
      (
        "Select rows of the dataset with the given ids"
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
      member_data.AddOptionalInitializer
      (
        "id type",
        "type of id to select (e.g. AASeqID, AAPdbID).  "
        "If not given, the complete set of ids selected by -id_labels is used",
        io::Serialization::GetAgent( &m_IdName)
      );
      member_data.AddInitializer
      (
        "ids",
        "actual ids to select",
        io::Serialization::GetAgent( &m_Ids)
      );
      member_data.AddInitializer
      (
        "invert",
        "select only rows that lack this id",
        io::Serialization::GetAgent( &m_Invert),
        "False"
      );
      return member_data;
    }

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t RetrieveDataSetById::GetNominalSize() const
    {
      return m_Retriever->GetNominalSize();
    }

    //! @brief test whether this retriever can generate sub-ranges of datasets without loading the entire dataset
    //! @return true if this retriever can generate sub-ranges of datasets without loading the entire dataset
    bool RetrieveDataSetById::SupportsEfficientSubsetLoading() const
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
    size_t RetrieveDataSetById::GenerateDataSubset
    (
      const math::Range< size_t> &SUBSET,
      linal::MatrixInterface< float> &FEATURES_STORAGE,
      linal::MatrixInterface< float> &RESULTS_STORAGE,
      linal::MatrixInterface< char> &IDS_STORAGE,
      const size_t &START_FEATURE_NUMBER
    )
    {
      UpdateIDsIndex();
      if( m_IdIndices.IsEmpty())
      {
        for( size_t i( 0), n_ids( IDS_STORAGE.GetNumberCols()); i < n_ids; ++i)
        {
          m_IdIndices.PushBack( i);
        }
      }

      // get the base data set
      const size_t n_rows( SUBSET.GetWidth() + 1);
      linal::Matrix< float> features( n_rows, FEATURES_STORAGE.GetNumberCols(), util::GetUndefined< float>());
      linal::Matrix< float> results( n_rows, RESULTS_STORAGE.GetNumberCols(), util::GetUndefined< float>());
      linal::Matrix< char> ids( n_rows, IDS_STORAGE.GetNumberCols(), ' ');
      size_t n_created( m_Retriever->GenerateDataSubset( SUBSET, features, results, ids, size_t( 0)));

      features.ShrinkRows( n_created);
      results.ShrinkRows( n_created);
      ids.ShrinkRows( n_created);

      if( !n_created)
      {
        BCL_MessageStd( "No features kept in this block pre-filter");
        return 0;
      }

      // select the desired features
      linal::MatrixConstReference< char> reference( ids);
      storage::Vector< size_t> features_to_keep( GetDesiredIDs( reference, m_IdIndices));
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

    //! @brief update the m_IdIndices
    void RetrieveDataSetById::UpdateIDsIndex()
    {
      if( m_IdIndices.IsEmpty())
      {
        if( !m_IdName.empty())
        {
          // create a label from the given name
          util::ObjectDataLabel label( m_IdName);
          if( RetrieveDataSetBase::GetIdCode().IsEmpty())
          {
            this->SelectIds( label);
          }
          m_IdIndices = this->GetIdCodeWithSizes().GetPropertyIndices( label);
        }
      }
    }

    //! @brief Set the members of this object from the given LABEL
    //! @param LABEL the label containing members that should be read of this class
    //! @param ERROR_STREAM stream with which to write errors
    bool RetrieveDataSetById::ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM)
    {
      // trim all the input strings
      storage::Set< std::string> ids_l;
      for( auto itr( m_Ids.Begin()), itr_end( m_Ids.End()); itr != itr_end; ++itr)
      {
        ids_l.Insert( util::TrimString( *itr));
      }
      m_Ids.InternalData().swap( ids_l.InternalData());
      return true;
    }

  } // namespace model
} // namespace bcl
