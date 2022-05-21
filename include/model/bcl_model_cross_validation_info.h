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

#ifndef BCL_MODEL_CROSS_VALIDATION_INFO_H_
#define BCL_MODEL_CROSS_VALIDATION_INFO_H_

// include the namespace header
#include "bcl_model.h"

// include other forward headers - sorted alphabetically
#include "descriptor/bcl_descriptor.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_model_objective_function_interface.h"
#include "util/bcl_util_serializable_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CrossValidationInfo
    //! @brief Storage class that stores the dataset labels used for monitoring and independently assessing a
    //!        model and the final objective function result and parity, and generically all meta-data associated with
    //!        training of the model, except the feature and result descriptors.
    //!
    //!
    //! @see @link example_model_cross_validation_info.cpp @endlink
    //! @author mendenjl
    //! @date May 16, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CrossValidationInfo :
      public util::SerializableInterface
    {
    public:

    //////////
    // enum //
    //////////

      // how the prediction file is formatted
      enum PredictionOutputFormat
      {
        e_ToBeDetermined,   //!< Could be any of the below
        e_Matrix,           //!< Output as a simple linal::Matrix< float>, with the first N-result columns being the
                            //!< experimental/target values, the second N columns being the predicted
        e_Delimited,        //!< Output as ID,experimental result 1,experimental result 2,...,predicted result 1
        s_NumberPredictionOutputFormats
      };

      //! @brief Name for each prediction file format
      //! @param FORMAT the format to convert to a string
      //! @return string representing FORMAT
      static const std::string &GetPredictionOutputFormatString( const PredictionOutputFormat &FORMAT);

      //! typedef for wrapper enum for this class; assists with I/O
      typedef util::WrapperEnum< PredictionOutputFormat, &GetPredictionOutputFormatString, s_NumberPredictionOutputFormats>
        PredictionOutputFormatEnum;

    private:

    //////////
    // data //
    //////////

      //! result value on the independent dataset
      float m_Result;
      //! get the directionality of result improvement
      opti::ImprovementTypeEnum m_ImprovementType;
      //! independent dataset retriever
      util::ObjectDataLabel m_Independent;
      //! monitoring dataset retriever
      util::ObjectDataLabel m_Monitoring;
      //! training dataset retriever
      util::ObjectDataLabel m_Training;
      //! Final objective function
      util::ObjectDataLabel m_Objective;
      //! Actual iterate used
      util::ObjectDataLabel m_Iterate;
      //! Where independent predictions are stored, empty if they were not stored anywhere
      std::string           m_IndependentPredictions;
      //! Number of output columns for this model
      size_t                m_NumberOutputs;
      //! ID-descriptors that are included in each of the predictions files
      util::ObjectDataLabel m_IdLabels;
      //! Characters in each ID label
      storage::Vector< size_t> m_CharsPerId;
      //! Output format
      PredictionOutputFormatEnum m_PredictionOutputFormat;
      //! filename where this cross-validation info is located
      std::string m_Filename;

      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      CrossValidationInfo();

      //! @brief constructor from filename
      explicit CrossValidationInfo( const std::string &FILENAME);

      //! @brief constructor from members
      CrossValidationInfo
      (
        const float &RESULT,
        const opti::ImprovementTypeEnum &IMPROVEMENT_TYPE,
        const util::ObjectDataLabel &INDEPENDENT,
        const util::ObjectDataLabel &MONITORING,
        const util::ObjectDataLabel &TRAINING,
        const util::ObjectDataLabel &OBJECTIVE,
        const util::ObjectDataLabel &ITERATE,
        const std::string &INDEPENDENT_PREDICTIONS,
        const FeatureLabelSet &ID_LABELS,
        const size_t &NUMBER_OUTPUTS
      );

      //! @brief Clone function
      //! @return pointer to new CrossValidationInfo
      CrossValidationInfo *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief get the result value
      //! @return the result
      float GetResult() const
      {
        return m_Result;
      }

      //! @brief get the improvement type
      //! @return the improvement type
      const opti::ImprovementTypeEnum &GetImprovementType() const
      {
        return m_ImprovementType;
      }

      //! @brief get the independent data chunk id
      //! @return the independent data chunk id
      const util::ObjectDataLabel &GetIndependentDatasetRetriever() const
      {
        return m_Independent;
      }

      //! @brief get the monitoring data chunk id
      //! @return the monitoring data chunk id
      const util::ObjectDataLabel &GetMonitoringDatasetRetriever() const
      {
        return m_Monitoring;
      }

      //! @brief get the training data chunk id
      //! @return the training data chunk id
      const util::ObjectDataLabel &GetTrainingDatasetRetriever() const
      {
        return m_Training;
      }

      //! @brief get the training data chunk id
      //! @return the training data chunk id
      const util::ObjectDataLabel &GetIterate() const
      {
        return m_Iterate;
      }

      //! @brief get the objective function label
      //! @return the objective function label
      const util::ObjectDataLabel &GetObjective() const
      {
        return m_Objective;
      }

      //! @brief get the filename for the independent predictions (blank if no independent predictions were made)
      //! @return the filename for the independent predictions
      const std::string &GetIndependentPredictionsFilename() const
      {
        return m_IndependentPredictions;
      }

      //! @brief set the filename for the independent predictions (blank if no independent predictions were made)
      //! @param FILENAME the independent predictions filename
      //! @return *this
      CrossValidationInfo &SetIndependentPredictionsFilename( const std::string &FILENAME)
      {
        m_IndependentPredictions = FILENAME;
        return *this;
      }

      //! @brief set the filename for the independent predictions (blank if no independent predictions were made)
      //! @param IND_DATASET the independent dataset label
      //! @return *this
      CrossValidationInfo &SetIndependentDataset( const util::ObjectDataLabel &IND_DATASET)
      {
        m_Independent = IND_DATASET;
        return *this;
      }

      //! @brief set the predictions format used by *this for reading and writing predictions
      //! @param FORMAT the desired format
      //! @return *this
      CrossValidationInfo &SetPredictionsFormat( const PredictionOutputFormat &FORMAT)
      {
        m_PredictionOutputFormat = FORMAT;
        return *this;
      }

      //! @brief get the number of output values for this model
      //! @return the number of output values for this model (0 if unknown)
      const size_t &GetNumberOutputs() const
      {
        return m_NumberOutputs;
      }

      //! @brief get the id feature label set
      //! @return the id feature label set
      FeatureLabelSet GetIDsFeatureLabelSet() const;

      //! @brief keep only the information in this cv-info that is held in common with OTHER
      //! @param OTHER the other cross-validation info object
      void KeepSharedInfo( const CrossValidationInfo &OTHER);

      //! @brief read the independent dataset predictions into a dataset object, whose features are the predicted values
      //!        results are the same result values, and ids are the given ids
      util::ShPtr< descriptor::Dataset> ReadIndependentPredictions();

      //! @brief read predictions from any file
      //! @param FILENAME the filename
      static util::ShPtr< descriptor::Dataset> ReadPredictions( const std::string &FILENAME);

      //! @brief write independent predictions using the current format, to the independent prediction file
      //! @param EXPERIMENTAL matrix with experimental data; non-const to allow descaling
      //! @param PREDICTIONS predicted results; non-const to allow descaling
      //! @param IDS ids matrix
      //! @param FILENAME file to write the predictions out to
      //! @param FORMAT the desired format to write the predictions in
      static void WritePredictions
      (
        FeatureDataSet< float> &EXPERIMENTAL,
        FeatureDataSet< float> &PREDICTIONS,
        const FeatureDataSetInterface< char> &IDS,
        const std::string &FILENAME,
        const PredictionOutputFormat &FORMAT = e_Delimited
      );

      //! @brief write independent predictions using the current format, to the independent prediction file
      //! @param EXPERIMENTAL matrix with experimental data; non-const to allow descaling
      //! @param PREDICTIONS predicted results; non-const to allow descaling
      //! @param IDS ids matrix
      //! @param FILENAME file to write the predictions out to
      //! @param FORMAT the desired format to write the predictions in
      static void WritePredictions
      (
        const linal::MatrixConstInterface< float> &EXPERIMENTAL,
        const linal::MatrixConstInterface< float> &PREDICTIONS,
        const linal::MatrixConstInterface< char> &IDS,
        const std::string &FILENAME,
        const PredictionOutputFormat &FORMAT = e_Delimited
      );

    private:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief hook that derived classes can override, called after TryRead successfully reads a serializer containing only initializer info (no data variables)
      //! @param SERIALIZER the serializer object with initialization information
      //! @param ERR_STREAM stream to write out errors to
      //! @return true, unless there were new errors
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &SERIALIZER, std::ostream &ERR_STREAM);

    }; // class CrossValidationInfo

  } // namespace model
} // namespace bcl

#endif // BCL_MODEL_CROSS_VALIDATION_INFO_H_
