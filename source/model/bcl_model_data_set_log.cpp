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
#include "model/bcl_model_data_set_log.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "linal/bcl_linal_matrix_operations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // add the interface to the set of known implementations
    const util::SiPtr< const util::ObjectInterface> DataSetLog::s_Instance
    (
      util::Enumerated< RetrieveDataSetBase>::AddInstance( new DataSetLog())
    );

    //! @brief default constructor
    DataSetLog::DataSetLog() :
        m_Label()
    {
    }

    //! @brief constructor from the rmsd
    //! @param RMSD the rmsd to use in the clustering algorithm
    //! @param RANGE the min and max # of clusters that would be acceptable
    //! @param MAX_STEPS_TO_REACH_CLUSTER_SIZE the maximum # of attempts to change the rmsd to get the # of clusters into the range
    //! @param AUTOSCALE true if the data need to be normalized before clustering
    DataSetLog::DataSetLog
    (
      const std::string &ID_LABEL
    ) :
      m_Label( ID_LABEL)
    {
    }

    //! @brief Clone function
    //! @return pointer to new DataSetLog
    DataSetLog *DataSetLog::Clone() const
    {
      return new DataSetLog( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &DataSetLog::GetClassIdentifier() const
    {
      return GetStaticClassName< DataSetLog>();
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &DataSetLog::GetAlias() const
    {
      static const std::string s_Name( "Log");
      return s_Name;
    }

    //! @brief Set the code / label for the feature (1st part) of the data set
    //! @param CODE the new code
    //! @return the code / label for the feature (1st part) of the data set
    void DataSetLog::SelectFeatures( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectFeatures( CODE);
      m_Retriever->SelectFeatures( CODE);
    }

    //! @brief Set the code / label for the result (2nd part) of the data set
    //! @return the code / label for the result (2nd part) of the data set
    void DataSetLog::SelectResults( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectResults( CODE);
      m_Retriever->SelectResults( CODE);
    }

    //! @brief Set which id columns to retrieve
    //! @param CODE the id column names to retrieve
    void DataSetLog::SelectIds( const util::ObjectDataLabel &CODE)
    {
      RetrieveDataSetBase::SelectIds( CODE);
      m_Retriever->SelectIds( CODE);
    }

    //! @brief get the number of partitions requested by the user, along with the partition ids
    //! @return the number of partitions requested by the user, along with the partition ids
    storage::Pair< size_t, math::RangeSet< size_t> > DataSetLog::GetNumberPartitionsAndIds() const
    {
      return m_Retriever->GetNumberPartitionsAndIds();
    }

    //! @brief Get the code / label set for the feature (1st part) of the data set with sizes of each property
    //! @return the code / label for the feature (1st part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet DataSetLog::GetFeatureLabelsWithSizes() const
    {
      return m_Retriever->GetFeatureLabelsWithSizes();
    }

    //! @brief Get the code / label for the result (2nd part) of the data set with sizes of each property
    //! @return the code / label for the result (2nd part) of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet DataSetLog::GetResultCodeWithSizes() const
    {
      return m_Retriever->GetResultCodeWithSizes();
    }

    //! @brief Get the code / label for the ids of the data set with sizes of each property
    //! @return the code / label for the ids of the data set with sizes of each property
    //! the feature code set
    FeatureLabelSet DataSetLog::GetIdCodeWithSizes() const
    {
      return m_Retriever->GetIdCodeWithSizes();
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
      DataSetLog::GenerateDataSet()
    {
      // get the base data set and use the clustering algorithm on it
      util::ShPtr< descriptor::Dataset> changeable( m_Retriever->GenerateDataSet());
      this->operator()( changeable);
      return changeable;
    }

    //! @brief reduce a data set
    //! @param DATA the data set
    //! @return the reduced data set
    void DataSetLog::operator ()( util::ShPtr< descriptor::Dataset> DATA) const
    {
      linal::MatrixReference< float> feature_matrix( DATA->GetFeatures().GetMatrix());
      linal::MatrixReference< float> result_matrix( DATA->GetResults().GetMatrix());
      for( linal::MatrixReference< float>::iterator itr_m( feature_matrix.Begin()), itr_m_end( feature_matrix.End());
          itr_m != itr_m_end;
          ++itr_m
         )
      {
        *itr_m = ( *itr_m < 0 ? -1.0 : 1.0) * std::log( std::fabs( *itr_m));
      }

      if( m_Label.empty())
      {
        return;
      }

      linal::MatrixReference< char> id_mat( DATA->GetIds().GetMatrix());
      BCL_Assert
      (
        id_mat.GetNumberCols() >= m_Label.length(),
        "Specified label \"" + m_Label + "\" is too long for the given dataset (no fewer than "
          + util::Format()( id_mat.GetNumberCols()) + " characters allowed)"
      );

      std::string padded_label( m_Label);
      padded_label += std::string( id_mat.GetNumberCols() - m_Label.length(), ' ');
      for( size_t r( 0), end_r( id_mat.GetNumberRows()); r < end_r; ++r)
      {
        std::copy( padded_label.begin(), padded_label.end(), id_mat[ r]);
      }
    }

    //! @brief reduce a data set
    //! @param DATA the data set
    //! @return the reduced data set
    util::ShPtr< descriptor::Dataset> DataSetLog::operator()( const descriptor::Dataset &DATA) const
    {
      util::ShPtr< descriptor::Dataset> changeable( DATA.HardCopy());
      this->operator()( changeable);
      return changeable;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer DataSetLog::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription
      (
        "Multiplies the feature and result values of a dataset by a fixed value"
      );
      member_data.AddInitializer
      (
        "dataset",
        "dataset retriever to call to get the entire data set",
        io::Serialization::GetAgent( &m_Retriever)
      );

      member_data.AddOptionalInitializer
      (
        "id label",
        "the id label to add to multiplied data points",
        io::Serialization::GetAgent( &m_Label)
      );

      return member_data;
    }

    //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
    size_t DataSetLog::GetNominalSize() const
    {
      return m_Retriever->GetNominalSize();
    }

  } // namespace model
} // namespace bcl
