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

#ifndef BCL_MODEL_DATA_SET_SELECT_COLUMNS_H
#define BCL_MODEL_DATA_SET_SELECT_COLUMNS_H

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DataSetSelectColumns
    //! @brief Forms a new dataset row from selected columns of another data set row
    //!
    //! @author mendenjl
    //! @see @link example_model_data_set_select_columns.cpp @endlink
    //! @date Dec 17, 2010
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DataSetSelectColumns :
      public util::SerializableInterface
    {
    private:

    //////////
    // data //
    //////////

      storage::Vector< size_t> m_ColumnIndices; //!< indices of the feature that should be kept
      size_t m_FeatureSize;                     //!< op() will exit with feature vectors not of this length

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param FEATURE_SIZE the expected size of the features
      //! @param INDICES indices of the features to keep
      DataSetSelectColumns
      (
        const size_t &FEATURE_SIZE = 0,
        const storage::Vector< size_t> &INDICES = storage::Vector< size_t>()
      );

      //! @brief Clone function
      //! @return pointer to new DataSetSelectColumns
      DataSetSelectColumns *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the size of vectors expected by operator()
      //! @return the size of vectors expected by operator()
      size_t GetInputFeatureSize() const;

      //! @brief return the size of vectors returned by operator()
      //! @return the size of vectors returned by operator()
      size_t GetOutputFeatureSize() const;

      //! @brief return the indices of the input features returned by operator()
      //! @return the indices of the input features returned by operator()
      const storage::Vector< size_t> &GetColumnIndices() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief reduce a feature to have only the selected features
      //! @param FEATURES a row of a data set
      //! @param STORAGE a pointer to where begin storing the selected values
      void operator ()
      (
        const linal::VectorConstInterface< float> &FEATURES,
        float *STORAGE
      ) const;

      //! @brief reduce a feature to have only the selected features
      //! @param FEATURES a row of a data set
      //! @param STORAGE a pointer to where begin storing the selected values
      void operator ()
      (
        const linal::VectorConstInterface< char> &FEATURES,
        char *STORAGE
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      virtual io::Serializer GetSerializer() const;

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief derived classes can change the indices that will be retained
      //! @param INDICES indices of the features to keep
      void SetColumnIndices( const storage::Vector< size_t> &INDICES);

    }; // class DataSetSelectColumns

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_DATA_SET_SELECT_COLUMNS_H

