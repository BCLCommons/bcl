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

#ifndef BCL_MODEL_RETRIEVE_DATA_SET_BY_ID_H_
#define BCL_MODEL_RETRIEVE_DATA_SET_BY_ID_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_model_retrieve_data_set_base.h"
#include "math/bcl_math_range_set.h"
#include "util/bcl_util_function_interface_serializable.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_object_data_label.h"
#include "util/bcl_util_si_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RetrieveDataSetById
    //! @brief Selects rows of a dataset based on a given id
    //!
    //! @author mendenjl
    //! @see @link example_model_retrieve_data_set_by_id.cpp @endlink
    //! @date Oct 04, 2013
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RetrieveDataSetById :
      public RetrieveDataSetBase,
      public util::FunctionInterfaceSerializable< util::ShPtr< descriptor::Dataset>, util::ShPtr< descriptor::Dataset> >
    {
    private:

    //////////
    // data //
    //////////

      std::string                                m_IdName;       //!< Name of the id
      storage::Set< std::string>                 m_Ids;          //!< ids of interest
      util::Implementation< RetrieveDataSetBase> m_Retriever;    //!< Retrieval method for the data set
      storage::Vector< size_t>                   m_IdIndices;    //!< Start index for the ids
      bool                                       m_Invert;       //!< Retrieve only ids that are not in m_Ids
      bool                                       m_JustFilter;   //!< True if this instance is only expected to filter datasets

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
      static const util::SiPtr< const util::ObjectInterface> s_FilterInstance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from parameters
      //! @param IDS the ids to consider
      //! @param RETRIEVER the dataset retriever to use
      //! @param JUST_FILTER set this to have this object run in pure-filter mode, in which case the dataset will not be
      //!        requested from the command line
      RetrieveDataSetById
      (
        const storage::Set< std::string> &IDS = storage::Set< std::string>(),
        const util::Implementation< RetrieveDataSetBase> &RETRIEVER = util::Implementation< RetrieveDataSetBase>(),
        const bool &JUST_FILTER = false
      );

      //! @brief Clone function
      //! @return pointer to new RetrieveDataSetById
      RetrieveDataSetById *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the ids of interest
      //! @return the ids of interest
      const storage::Set< std::string> &GetSelectedIds() const
      {
        return m_Ids;
      }

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

      //! @brief get whether dataset generation requires labels
      //! @return true if dataset generation requires labels
      bool RequiresFeatureLabels() const
      {
        return m_Retriever->RequiresFeatureLabels();
      }

      //! @brief get whether dataset generation requires result labels
      //! @return true if dataset generation requires result labels
      bool RequiresResultLabels() const
      {
        return m_Retriever->RequiresResultLabels();
      }

      //! @brief get the number of partitions requested by the user, along with the partition ids
      //! @return the number of partitions requested by the user, along with the partition ids
      storage::Pair< size_t, math::RangeSet< size_t> > GetNumberPartitionsAndIds() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief generate dataset
      //! @return generated dataset
      util::ShPtr< descriptor::Dataset> GenerateDataSet();

      //! @brief operator() can be used to filter a dataset
      //! @param DATASET the dataset to filter
      //! @return the filtered dataset
      util::ShPtr< descriptor::Dataset> operator()( const util::ShPtr< descriptor::Dataset> &DATASET) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

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

      //! @brief get indices of the desired id rows
      //! @param IDS the matrix of ids
      //! @return indices of the desired ids
      storage::Vector< size_t> GetDesiredIDs
      (
        const linal::MatrixConstInterface< char> &IDS,
        const storage::Vector< size_t> &ID_INDICES
      ) const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
      //! @return the nominal (e.g. best estimate without generating the entire dataset) size of the dataset
      size_t GetNominalSize() const;

      //! @brief update the m_IdIndices
      void UpdateIDsIndex();

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    }; // class RetrieveDataSetById

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_RETRIEVE_DATA_SET_BY_ID_H_

