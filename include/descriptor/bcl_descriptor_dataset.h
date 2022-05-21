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

#ifndef BCL_DESCRIPTOR_DATASET_H_
#define BCL_DESCRIPTOR_DATASET_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "model/bcl_model_feature_data_set.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Dataset
    //! @brief Mutable dataset containing features, results, identification, and rescale
    //!
    //! @see @link example_descriptor_dataset.cpp @endlink
    //! @author mendenjl
    //! @date Dec 13, 2012
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Dataset :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

      //! Features, used to make predictions
      util::ShPtr< model::FeatureDataSet< float> > m_Features;

      //! Results, may be predicted or experimental
      util::ShPtr< model::FeatureDataSet< float> > m_Results;

      //! Identification for each row in the feature dataset
      //! This too has to be a matrix to maintain fixed width
      util::ShPtr< model::FeatureDataSet< char> > m_Identification;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      Dataset();

      //! @brief constructor from sizes
      //! @param NR_EXAMPLES number of training examples the dataset will contain
      //! @param FEATURE_SIZE the size of each feature in the dataset
      //! @param RESULT_SIZE the size of each result in the dataset
      //! @param ID_SIZE the size of each id in the dataset
      Dataset
      (
        const size_t &NR_EXAMPLES,
        const size_t &FEATURE_SIZE,
        const size_t &RESULT_SIZE,
        const size_t &ID_SIZE
      );

      //! @brief constructor from labels
      //! @param NR_EXAMPLES number of training examples the dataset will contain
      //! @param FEATURE_LABELS labels for all features in the dataset
      //! @param RESULT_LABELS labels for each result in the dataset
      //! @param ID_LABELS labels for each id in the dataset
      Dataset
      (
        const size_t &NR_EXAMPLES,
        const model::FeatureLabelSet &FEATURE_LABELS,
        const model::FeatureLabelSet &RESULT_LABELS,
        const model::FeatureLabelSet &ID_LABELS
      );

      //! @brief constructor
      //! @param FEATURES the features
      //! @param RESULTS the results
      Dataset
      (
        const util::ShPtr< model::FeatureDataSet< float> > &FEATURES,
        const util::ShPtr< model::FeatureDataSet< float> > &RESULTS
      );

      //! @brief constructor
      //! @param FEATURES the features
      //! @param RESULTS the results
      Dataset
      (
        const linal::MatrixConstInterface< float> &FEATURES,
        const linal::MatrixConstInterface< float> &RESULTS
      );

      //! @brief constructor
      //! @param FEATURES the features
      //! @param RESULTS the results
      //! @param IDS the ids
      Dataset
      (
        const linal::MatrixConstInterface< float> &FEATURES,
        const linal::MatrixConstInterface< float> &RESULTS,
        const linal::MatrixConstInterface< char> &IDS
      );

      //! @brief constructor from vector of features and results
      //! @param DATA_SET the vector of features and results
      //! @note this constructor is DEPRECATED : use matrices in newly written code
      Dataset
      (
        const storage::Vector< storage::VectorND< 2, linal::Vector< float> > > &DATA_SET
      );

      //! @brief Clone function
      //! @return pointer to new Dataset
      Dataset *Clone() const;

      //! @brief Hard copy -> copies internal data as well
      //! @return pointer to new Dataset
      Dataset *HardCopy() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return number of examples in the data set
      //! @return returns number of examples in the Dataset
      //size_t GetSize() const;
      //! @brief return number of examples in the data set
      //! @return returns number of examples in the Dataset
      size_t GetSize() const;

      //! @brief return the size of each feature in the dataset
      //! @return return the size of each feature in the dataset
      size_t GetFeatureSize() const;

      //! @brief return the size of each id (# characters) in the dataset
      //! @return return the size of each id (# characters) in the dataset
      size_t GetIdSize() const;

      //! @brief return the size of each result in the dataset
      //! @return return the size of each result in the dataset
      size_t GetResultSize() const;

      //! @brief check whether features or results contain a number of zero rows
      //! @return true if features or results have zero rows in matrix
      bool IsEmpty() const;

      //! @brief returns the features pointer
      //! @return features pointer
      const util::ShPtr< model::FeatureDataSet< float> > &GetFeaturesPtr() const;

      //! @brief returns the results pointer
      //! @return the results pointer
      const util::ShPtr< model::FeatureDataSet< float> > &GetResultsPtr() const;

      //! @brief returns the ids pointer
      //! @return the ids pointer
      const util::ShPtr< model::FeatureDataSet< char> > &GetIdsPtr() const;

      //! @brief returns the features pointer
      //! @return features pointer
      model::FeatureDataSet< float> &GetFeatures();

      //! @brief returns the results pointer
      //! @return the results pointer
      model::FeatureDataSet< float> &GetResults();

      //! @brief returns the ids dataset
      //! @return the ids dataset
      model::FeatureDataSet< char> &GetIds();

      //! @brief access to the features matrix
      //! @return a features matrix reference
      linal::MatrixReference< float> GetFeaturesReference();

      //! @brief access to the features matrix
      //! @return a features matrix reference
      linal::MatrixReference< float> GetResultsReference();

      //! @brief access to the features matrix
      //! @return an id matrix reference
      linal::MatrixReference< char> GetIdsReference();

      //! @brief const access to the features matrix
      //! @return a features matrix reference
      linal::MatrixConstReference< float> GetFeaturesReference() const;

      //! @brief const access to the features matrix
      //! @return a features matrix reference
      linal::MatrixConstReference< float> GetResultsReference() const;

      //! @brief const access to the features matrix
      //! @return an id matrix reference
      linal::MatrixConstReference< char> GetIdsReference() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief remove all undefined examples
      //! @return the number of examples actually removed
      size_t RemoveUndefinedExamples();

      //! @brief append another data set onto this one
      //! @param DATASET the dataset to append
      void Append( const Dataset &DATASET);

      //! @brief add a specified number of extra rows to this dataset
      //! @param N_ROWS the number of rows to add
      void AddRows( const size_t &N_ROWS);

      //! @brief shrink rows (remove unused rows at the end of the dataset)
      //! @param NEW_SIZE the new size of the matrix
      void ShrinkRows( const size_t &NEW_SIZE);

      //! @brief keep only the specified rows
      //! @param KEEPERS indices of rows to keep
      void KeepRows( const storage::Vector< size_t> &KEEPERS);

      //! @brief get a new dataset consisting of only the specified rows
      //! @param ROWS the rows desired
      util::ShPtr< Dataset> GetRows( const storage::Vector< size_t> &KEEPERS) const;

      //! @brief Shuffle the rows randomly
      void Shuffle();

      //! @brief Shuffle only the results (y-scrambling)
      void YScramble();

      //! @brief add data for a specific row
      //! @param ROW the row to add data to
      //! @param FEATURE the feature to add
      //! @param RESULT the result to add
      //! @param ID the id information to add
      void AddData
      (
        const size_t &ROW,
        const linal::VectorConstInterface< float> &FEATURE,
        const linal::VectorConstInterface< float> &RESULT,
        const linal::VectorConstInterface< char> &ID
      );

      //! @brief set features
      //! @param FEATURES new FeatureDataSet
      void SetIds( const util::ShPtr< model::FeatureDataSet< char> > &IDS);

      //! @brief set features
      //! @param FEATURES new FeatureDataSet
      void SetFeatures( const util::ShPtr< model::FeatureDataSet< float> > &FEATURES);

      //! @brief set results
      //! @param RESULTS new FeatureDataSet
      void SetResults( const util::ShPtr< model::FeatureDataSet< float> > &RESULTS);

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class Dataset

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_DATASET_H_
