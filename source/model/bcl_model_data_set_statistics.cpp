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
#include "model/bcl_model_data_set_statistics.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new DataSetStatistics
    DataSetStatistics *DataSetStatistics::Clone() const
    {
      return new DataSetStatistics( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DataSetStatistics::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief compute the data set statistics from a dataset retriever
    //! @param DATASET some means of making a dataset retriever
    void DataSetStatistics::AddDataSet( util::Implementation< RetrieveDataSetBase> DATASET)
    {
      util::ShPtr< descriptor::Dataset> dataset( DATASET->GenerateDataSet());
      this->AddDataSet( *dataset);
    }

    //! @brief compute the data set statistics from a dataset
    //! @param DATASET a dataset
    void DataSetStatistics::AddDataSet( const descriptor::Dataset &DATASET)
    {
      const size_t n_features( DATASET.GetSize());
      const size_t n_feature_cols( DATASET.GetFeatureSize());
      const size_t n_result_cols( DATASET.GetResultSize());

      for( size_t feature( 0); feature < n_features; ++feature)
      {
        linal::VectorConstReference< float> feature_ref( n_feature_cols, DATASET.GetFeaturesReference()[ feature]);
        linal::VectorConstReference< float> result_ref( n_result_cols, DATASET.GetResultsReference()[ feature]);
        m_AveStdFeatures += feature_ref;
        m_MinMaxFeatures += feature_ref;
        m_AveStdResults += result_ref;
        m_MinMaxResults += result_ref;
      }
    }

    //! @brief compute the data set statistics as the difference between two other datasets
    //! @param DATASET_A, DATASET_B actual datasets
    void DataSetStatistics::AddDifference
    (
      const descriptor::Dataset &DATASET_A,
      const descriptor::Dataset &DATASET_B
    )
    {
      descriptor::Dataset difference( DATASET_A);
      linal::MatrixReference< float> features( difference.GetFeaturesReference());
      linal::MatrixReference< float> results( difference.GetResultsReference());
      features -= DATASET_B.GetFeaturesReference();
      results -= DATASET_B.GetResultsReference();
      AddDataSet( difference);
    }

    //! @brief compute the data set statistics as the difference between two other datasets
    //! @param DATASET_A, DATASET_B some means of making a dataset
    void DataSetStatistics::AddDifference
    (
      util::Implementation< RetrieveDataSetBase> DATASET_A,
      util::Implementation< RetrieveDataSetBase> DATASET_B
    )
    {
      util::ShPtr< descriptor::Dataset> dataset_a( DATASET_A->GenerateDataSet());
      util::ShPtr< descriptor::Dataset> dataset_b( DATASET_B->GenerateDataSet());
      linal::MatrixReference< float> features( dataset_a->GetFeaturesReference());
      linal::MatrixReference< float> results( dataset_a->GetResultsReference());
      features -= dataset_b->GetFeaturesReference();
      results -= dataset_b->GetResultsReference();
      AddDataSet( *dataset_a);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DataSetStatistics::Read( std::istream &ISTREAM)
    {
      io::Serialize::Read( m_AveStdFeatures, ISTREAM);
      io::Serialize::Read( m_MinMaxFeatures, ISTREAM);
      io::Serialize::Read( m_AveStdResults, ISTREAM);
      io::Serialize::Read( m_MinMaxResults, ISTREAM);
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DataSetStatistics::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      io::Serialize::Write( m_AveStdFeatures, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinMaxFeatures, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_AveStdResults, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MinMaxResults, OSTREAM, INDENT);
      return OSTREAM;
    }

  } // namespace model
} // namespace bcl
