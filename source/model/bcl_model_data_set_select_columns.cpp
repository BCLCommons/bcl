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
#include "io/bcl_io_serialization.h"
#include "model/bcl_model_data_set_select_columns.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> DataSetSelectColumns::s_Instance
    (
      GetObjectInstances().AddInstance( new DataSetSelectColumns())
    );

    //! @brief default constructor
    //! @param FEATURE_SIZE the expected size of the features
    //! @param INDICES indices of the features to keep
    DataSetSelectColumns::DataSetSelectColumns
    (
      const size_t &FEATURE_SIZE,
      const storage::Vector< size_t> &INDICES
    ) :
      m_ColumnIndices( INDICES),
      m_FeatureSize( FEATURE_SIZE)
    {
    }

    //! @brief Clone function
    //! @return pointer to new DataSetSelectColumns
    DataSetSelectColumns *DataSetSelectColumns::Clone() const
    {
      return new DataSetSelectColumns( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &DataSetSelectColumns::GetClassIdentifier() const
    {
      return GetStaticClassName< DataSetSelectColumns>();
    }

    //! @brief return the size of vectors expected by operator()
    //! @return the size of vectors expected by operator()
    size_t DataSetSelectColumns::GetInputFeatureSize() const
    {
      return m_FeatureSize;
    }

    //! @brief return the size of vectors returned by operator()
    //! @return the size of vectors returned by operator()
    size_t DataSetSelectColumns::GetOutputFeatureSize() const
    {
      return m_ColumnIndices.GetSize();
    }

    //! @brief return the indices of the input features returned by operator()
    //! @return the indices of the input features returned by operator()
    const storage::Vector< size_t> &DataSetSelectColumns::GetColumnIndices() const
    {
      return m_ColumnIndices;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief reduce a feature to have only the selected features
    //! @param FEATURES a row of a data set
    //! @param STORAGE a pointer to where begin storing the selected values
    void DataSetSelectColumns::operator ()
    (
      const linal::VectorConstInterface< float> &FEATURES,
      float *STORAGE
    ) const
    {
      // make sure the feature has the right size
      BCL_Assert
      (
        FEATURES.GetSize() == m_FeatureSize || !m_FeatureSize,
        "Feature given to " + GetClassIdentifier() + " with incorrect size: " + util::Format()( FEATURES.GetSize())
        + " expected size: " + util::Format()( m_FeatureSize)
      );

      // get an iterator to the beginning of the old feature
      const float *old_feature( FEATURES.Begin());

      // copy the features over in the desired order
      for
      (
        size_t new_feature_index( 0), new_feature_size( GetOutputFeatureSize());
        new_feature_index < new_feature_size;
        ++new_feature_index
      )
      {
        STORAGE[ new_feature_index] = old_feature[ m_ColumnIndices( new_feature_index)];
      }
    }

    //! @brief reduce a feature to have only the selected features
    //! @param FEATURES a row of a data set
    //! @param STORAGE a pointer to where begin storing the selected values
    void DataSetSelectColumns::operator ()
    (
      const linal::VectorConstInterface< char> &FEATURES,
      char *STORAGE
    ) const
    {
      // make sure the feature has the right size
      BCL_Assert
      (
        FEATURES.GetSize() == m_FeatureSize || !m_FeatureSize,
        "Feature given to " + GetClassIdentifier() + " with incorrect size: " + util::Format()( FEATURES.GetSize())
        + " expected size: " + util::Format()( m_FeatureSize)
      );

      // get an iterator to the beginning of the old feature
      const char *old_feature( FEATURES.Begin());

      // copy the features over in the desired order
      for
      (
        size_t new_feature_index( 0), new_feature_size( GetOutputFeatureSize());
        new_feature_index < new_feature_size;
        ++new_feature_index
      )
      {
        STORAGE[ new_feature_index] = old_feature[ m_ColumnIndices( new_feature_index)];
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer DataSetSelectColumns::GetSerializer() const
    {
      io::Serializer member_data;
      member_data.SetClassDescription
      (
        "Selects columns given indices"
      );

      member_data.AddInitializer
      (
        "",
        "indices of columns to select",
        io::Serialization::GetAgent( &m_ColumnIndices)
      );

      return member_data;
    } // GetParameters

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief derived classes can change the indices that will be retained
    //! @param INDICES indices of the features to keep
    void DataSetSelectColumns::SetColumnIndices( const storage::Vector< size_t> &INDICES)
    {
      m_ColumnIndices = INDICES;
    }

  } // namespace model
} // namespace bcl
