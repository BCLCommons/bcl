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
#include "model/bcl_model_retrieve_interface.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_dataset.h"
#include "math/bcl_math_roc_curve.h"
#include "math/bcl_math_running_average.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace model
  {
    //! @brief Retrieves all descriptors, asserts if they are not all identical
    //! @return the unique descriptor
    util::ObjectDataLabel RetrieveInterface::RetrieveUniqueDescriptor() const
    {
      storage::List< util::ObjectDataLabel> all_descriptors( RetrieveEnsembleDescriptors());
      BCL_Assert( !all_descriptors.IsEmpty(), "No descriptors retrieved!");
      const size_t n_descriptor_sets( all_descriptors.GetSize());
      BCL_Assert
      (
        n_descriptor_sets == size_t( 1)
        ||
        std::count_if
        (
          ++all_descriptors.Begin(),
          all_descriptors.End(),
          std::bind2nd( std::equal_to< util::ObjectDataLabel>(), all_descriptors.FirstElement())
        )
        == int( n_descriptor_sets - 1),
        "Multiple descriptor sets found!"
      );
      return all_descriptors.FirstElement();
    }

    //! @brief get the cv info that is held in common across all CV in the dataset
    //! @return the commonly-held cross-validation info
    CrossValidationInfo RetrieveInterface::RetrieveCommonCVInfo() const
    {
      storage::List< CrossValidationInfo> all_cv_info( RetrieveEnsembleCVInfo());
      if( all_cv_info.IsEmpty())
      {
        return CrossValidationInfo();
      }

      storage::List< CrossValidationInfo>::const_iterator itr( all_cv_info.Begin()), itr_end( all_cv_info.End());
      CrossValidationInfo common_cv( *itr);
      while( ++itr != itr_end)
      {
        common_cv.KeepSharedInfo( *itr);
      }
      return common_cv;
    }

    //! @brief read the merged independent dataset predictions into a dataset object, whose features are the predicted values
    //!        results are the same result values, and ids are the given ids
    util::ShPtr< descriptor::Dataset> RetrieveInterface::ReadMergedIndependentPredictions()
    {
      // map from cross-validation ids to filenames
      storage::Map< std::string, storage::Vector< CrossValidationInfo> > cv_filenames;
      // filenames need to be retrieved from the model storage
      storage::List< CrossValidationInfo> cv_info( this->RetrieveEnsembleCVInfo());
      storage::Vector< size_t> counts( cv_info.GetSize(), size_t( 0));
      for
      (
        storage::List< CrossValidationInfo>::const_iterator itr( cv_info.Begin()), itr_end( cv_info.End());
        itr != itr_end;
        ++itr
      )
      {
        storage::Vector< CrossValidationInfo> &vec
        (
          cv_filenames[ itr->GetIndependentDatasetRetriever().ToString()]
        );
        // track the # of times each cv independent id is seen
        ++counts( vec.GetSize());
        vec.PushBack( *itr);
      }

      // test whether all cv's have the same # of files. This command adds the number of times the counts histogram was
      // 0, which should be counts.GetSize() - 1 if all the cvs have the same # of files.
      const size_t count
      (
        std::count( counts.Begin(), counts.End(), size_t( 0))
        + std::count( counts.Begin(), counts.End(), counts( 0))
      );
      if( count != counts.GetSize())
      {
        // Not equal number of counts. Either there really were different numbers of cross validations run per model, or
        // else the cross validations were run from different clusters, in which case the apparent absolute path may
        // differ
        storage::Map< std::string, storage::Vector< CrossValidationInfo> > cv_filenames_diff;
        counts.SetAllElements( 0);
        for
        (
          storage::List< CrossValidationInfo>::const_iterator itr( cv_info.Begin()), itr_end( cv_info.End());
          itr != itr_end;
          ++itr
        )
        {
          // use the difference function to compute the difference of the independent retriever from the training retriever
          // The reason this is needed is because the underlying training dataset files may have been computed from
          // identical files at different locations (e.g. using scratch on a cluster vs. an explicit path when done locally)
          storage::Vector< CrossValidationInfo> &vec
          (
            cv_filenames_diff
            [
              itr->GetIndependentDatasetRetriever().Difference( itr->GetTrainingDatasetRetriever()).ToString()
            ]
          );
          // track the # of times each cv independent id is seen
          ++counts( vec.GetSize());
          vec.PushBack( *itr);
        }
        const size_t count_new
        (
          std::count( counts.Begin(), counts.End(), size_t( 0))
          + std::count( counts.Begin(), counts.End(), counts( 0))
        );
        if( count_new == counts.GetSize())
        {
          BCL_MessageCrt
          (
            "Warning; differing number of files for different cross validation independent sets. Fixed by "
            "assuming that the cross validations were run from different locations"
          );
          cv_filenames = cv_filenames_diff;
        }
        BCL_MessageStd
        (
          "Found " + util::Format()( cv_info.GetSize()) + " models with "
          + util::Format()( cv_filenames.GetSize()) + " independent chunks"
        );
      }

      // final container of prediction matrices
      storage::Vector< linal::Matrix< float> > final_prediction_matrices;
      storage::Vector< linal::Matrix< char> > final_id_matrices;
      storage::Vector< linal::Matrix< float> > final_result_matrices;

      CrossValidationInfo common_cv_info;
      bool first_in_all_cv( true);
      // iterate over all input file names and construct prediction matrices for each entry
      for
      (
        storage::Map< std::string, storage::Vector< CrossValidationInfo> >::iterator
          itr_independent_to_filenames( cv_filenames.Begin()), itr_independent_to_filenames_end( cv_filenames.End());
        itr_independent_to_filenames != itr_independent_to_filenames_end;
        ++itr_independent_to_filenames
      )
      {
        // container for prediction matrices
        storage::Vector< linal::Matrix< float> > prediction_matrices;
        FeatureDataSet< char> ids_matrix;
        FeatureDataSet< float> results_matrix;

        math::RunningAverage< linal::Matrix< float> > avg_prediction_matrices;
        bool is_first( true);
        for
        (
          storage::Vector< CrossValidationInfo>::iterator
            itr_cv_info( itr_independent_to_filenames->second.Begin()),
            itr_cv_info_end( itr_independent_to_filenames->second.End());
          itr_cv_info != itr_cv_info_end;
          ++itr_cv_info
        )
        {
          util::ShPtr< descriptor::Dataset> sp_dataset( itr_cv_info->ReadIndependentPredictions());
          BCL_Assert( sp_dataset.IsDefined(), "Could not read dataset given by " + itr_cv_info->GetLabel().ToString());

          // retain only the common cv information
          if( first_in_all_cv)
          {
            common_cv_info = *itr_cv_info;
            first_in_all_cv = false;
          }
          else
          {
            common_cv_info.KeepSharedInfo( *itr_cv_info);
          }
          if( is_first)
          {
            ids_matrix = sp_dataset->GetIds();
            results_matrix = sp_dataset->GetResults();
            is_first = false;
          }
          else
          {
            BCL_Assert( ids_matrix.GetMatrix() == sp_dataset->GetIdsReference(), "Incorrect ids matrix");
            BCL_Assert( results_matrix.GetMatrix() == sp_dataset->GetResultsReference(), "Incorrect results matrix");
          }

          avg_prediction_matrices += sp_dataset->GetFeaturesReference();
        }
        // append result for this cross validation to the final_prediction_matrices vector
        if( avg_prediction_matrices.GetWeight())
        {
          final_prediction_matrices.PushBack( avg_prediction_matrices.GetAverage());
        }
        final_result_matrices.PushBack( results_matrix.GetMatrix());
        final_id_matrices.PushBack( ids_matrix.GetMatrix());
      }

      // combine all the matrices together
      linal::Matrix< float> combined_prediction_matrix, combined_result_matrix;
      linal::Matrix< char> combined_id_matrix;
      combined_prediction_matrix.Append( final_prediction_matrices);
      combined_result_matrix.Append( final_result_matrices);
      combined_id_matrix.Append( final_id_matrices);

      util::ShPtr< descriptor::Dataset> sp_dataset
      (
        new descriptor::Dataset
        (
          combined_prediction_matrix,
          combined_result_matrix,
          combined_id_matrix
        )
      );

      if( common_cv_info.GetIDsFeatureLabelSet().GetSize())
      {
        util::ShPtr< FeatureDataSet< char> > sp_id_ptr( sp_dataset->GetIdsPtr());
        sp_id_ptr->SetFeatureLabelSet( common_cv_info.GetIDsFeatureLabelSet());
      }
      return sp_dataset;
    }

    //! @brief read the merged independent dataset predictions into a dataset object, whose features are the predicted values
    //!        results are the same result values, and ids are the given ids
    storage::Vector< math::ROCCurve> RetrieveInterface::GetMergedIndependentROCCurves()
    {
      // get the merged independent predictions
      util::ShPtr< descriptor::Dataset> sp_dataset( ReadMergedIndependentPredictions());
      const linal::MatrixReference< float> exp_ref( sp_dataset->GetResultsReference());
      const linal::MatrixReference< float> predictions_ref( sp_dataset->GetFeaturesReference());

      // get the objective function
      util::Implementation< ObjectiveFunctionInterface> objective( this->RetrieveCommonCVInfo().GetObjective());
      BCL_Assert
      (
        objective->GetGoalType() == ObjectiveFunctionInterface::e_RankClassification,
        "GetMergedIndependentROCCurve requires a model trained against a rank-classification-based classifier"
      );

      // construct roc curves
      return ROCCurvesFromDataset( exp_ref, predictions_ref, objective->GetThreshold(), objective->GetRankingParity());
    }

    //! @brief transform an experimental/predicted matrix into a series of ROC curves
    //! @param EXP experimental results
    //! @param PRED predicted results
    //! @param THRESHOLD threshold for defining Positive/Negative class split
    //! @param POS_ABOVE_THRESH true if positives are those above the threshold
    storage::Vector< math::ROCCurve> RetrieveInterface::ROCCurvesFromDataset
    (
      const linal::MatrixConstInterface< float> &EXP,
      const linal::MatrixConstInterface< float> &PRED,
      const float &THRESHOLD,
      const bool &POS_ABOVE_THRESH
    )
    {
      const size_t n_results( EXP.GetNumberCols());
      const size_t n_predictions( EXP.GetNumberRows());
      // create ROC curves for each result
      storage::Vector< math::ROCCurve> roc_curves( n_results);
      BCL_Assert( util::IsDefined( THRESHOLD), "A threshold is necessary for ROC curves");
      BCL_Assert( n_results == PRED.GetNumberCols(), "ConstructROCCurves given matrices of different sizes");
      for( size_t result_number( 0); result_number < n_results; ++result_number)
      {
        // collect all results for the result_number-th result
        storage::List< storage::Pair< double, double> > exp_predicted;
        for( size_t pred_number( 0); pred_number < n_predictions; ++pred_number)
        {
          exp_predicted.PushBack
          (
            storage::Pair< double, double>
            (
              PRED( pred_number, result_number),
              EXP( pred_number, result_number)
            )
          );
        }
        roc_curves( result_number) = math::ROCCurve( exp_predicted, THRESHOLD, POS_ABOVE_THRESH);
      }
      return roc_curves;
    }

  } // namespace model
} // namespace bcl
