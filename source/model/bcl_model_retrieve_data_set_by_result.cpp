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
#include "model/bcl_model_retrieve_data_set_by_result.h"

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
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetByResult::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new RetrieveDataSetByResult())
    );
    const util::SiPtr< const util::ObjectInterface> RetrieveDataSetByResult::s_FilterInstance
    (
      util::Enumerated< util::FunctionInterfaceSerializable< util::ShPtr< descriptor::Dataset>, util::ShPtr< descriptor::Dataset> > >::AddInstance
      (
        new RetrieveDataSetByResult
        (
          math::RangeSet< float>::GetCompleteRange(),
          util::Implementation< RetrieveDataSetBase>(),
          true
        )
      )
    );

    //! @brief constructor from parameters
    //! @param RANGES the ranges of results to consider
    //! @param RETRIEVER the dataset retriever to use
    //! @param FILTER_ONLY true if the instance will only be used to filter a dataset generated elsewhere in the code
    RetrieveDataSetByResult::RetrieveDataSetByResult
    (
      const math::RangeSet< float> &RANGES,
      const util::Implementation< RetrieveDataSetBase> &RETRIEVER,
      const bool &FILTER_ONLY
    ) :
      m_ResultRanges( RANGES),
      m_Retriever( RETRIEVER),
      m_JustFilter( FILTER_ONLY)
    {
    }

    //! @brief Clone function
    //! @return pointer to new RetrieveDataSetByResult
    RetrieveDataSetByResult *RetrieveDataSetByResult::Clone() const
    {
      return new RetrieveDataSetByResult( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &RetrieveDataSetByResult::GetClassIdentifier() const
    {
      return GetStaticClassName< RetrieveDataSetByResult>();
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &RetrieveDataSetByResult::GetAlias() const
    {
      static const std::string s_name( "ResultsInRange");
      return s_name;
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    //! @return the code / label for the feature (1st part) of the data set
    void RetrieveDataSetByResult::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectFeatures( CODE);
      m_Retriever->SelectFeatures( CODE);
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @return the code / label for the result (2nd part) of the data set
    void RetrieveDataSetByResult::SelectResults( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectResults( CODE);
      m_Retriever->SelectResults( CODE);
    }

    //! @brief Set which id columns to retrieve
    //! @param CODE the id column names to retrieve
    void RetrieveDataSetByResult::SelectIds( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectIds( CODE);
      m_Retriever->SelectIds( CODE);
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetByResult::GetFeatureLabelsWithSizes() const
    {
      return m_Retriever->GetFeatureLabelsWithSizes();
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetByResult::GetResultCodeWithSizes() const
    {
      return m_Retriever->GetResultCodeWithSizes();
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet RetrieveDataSetByResult::GetIdCodeWithSizes() const
    {
      return m_Retriever->GetIdCodeWithSizes();
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > RetrieveDataSetByResult::GetNumberPartitionsAndIds() const
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
      RetrieveDataSetByResult::GenerateDataSet()
    {
      // get the base data set and use the clustering algorithm on it
      util::ShPtr< descriptor::Dataset> dataset( m_Retriever->GenerateDataSet());

      // get the desired results
      const storage::Vector< size_t> desired_result_ids( GetResultsInRange( dataset->GetResultsReference()));

      // filter by those results
      dataset->KeepRows( desired_result_ids);

      // return
      return dataset;
    }

    //! @brief operator() can be used to filter a dataset
    //! @param DATASET the dataset to filter
    //! @return the filtered dataset
    util::ShPtr< descriptor::Dataset> RetrieveDataSetByResult::operator()( const util::ShPtr< descriptor::Dataset> &DATASET) const
    {
      return DATASET->GetRows( GetResultsInRange( DATASET->GetResultsReference()));
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief get a set of features that cover the same space as FEATURES at resolution RMSD
    //! @param FEATURES the set of features to consider
    //! @param RMSD the rmsd; vectors closer than this rmsd will be considered representative of each other
    //! @return indices of features that cover the same space as FEATURES at resolution RMSD
    storage::Vector< size_t> RetrieveDataSetByResult::GetResultsInRange( const linal::MatrixConstInterface< float> &RESULTS) const
    {
      storage::Vector< size_t> chosen_results;
      chosen_results.AllocateMemory( RESULTS.GetNumberRows());
      const size_t result_size( RESULTS.GetNumberCols());
      for( size_t row( 0), n_rows( RESULTS.GetNumberRows()); row < n_rows; ++row)
      {
        const float *row_ptr( RESULTS[ row]);
        for( size_t col( 0); col < result_size; ++col)
        {
          if( m_ResultRanges.IsWithin( row_ptr[ col]))
          {
            chosen_results.PushBack( row);
            break;
          }
        }
      }
      return chosen_results;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer RetrieveDataSetByResult::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription
      (
        "Select results (and corresponding features) within a particular range set"
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
        "ranges of result values to load, e.g. range=\"[ 0, 5)+(7,10)\"",
        io::Serialization::GetAgent( &m_ResultRanges),
        math::RangeSet< float>::GetCompleteRange().AsString()
      );

      return member_data;
    }

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t RetrieveDataSetByResult::GetNominalSize() const
    {
      return m_Retriever->GetNominalSize();
    }

  } // namespace model
} // namespace bcl
