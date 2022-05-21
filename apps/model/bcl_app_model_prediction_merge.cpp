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
#include "bcl_app_model_prediction_merge.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "descriptor/bcl_descriptor_dataset.h"
#include "io/bcl_io_file.h"
#include "linal/bcl_linal_matrix_reference.h"
#include "math/bcl_math_piecewise_function.h"
#include "math/bcl_math_roc_curve.h"
#include "math/bcl_math_running_min_max.h"
#include "model/bcl_model_feature_data_reference.h"
#include "model/bcl_model_objective_function_interface.h"
#include "model/bcl_model_retrieve_interface.h"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {
  ///////////////////////////////////
  // construction and destruction //
  ///////////////////////////////////

    //! @brief Clone function
    //! @return pointer to new FoldProtein
    ModelPredictionMerge *ModelPredictionMerge::Clone() const
    {
      return new ModelPredictionMerge( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ModelPredictionMerge::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> ModelPredictionMerge::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // filenames for prediction matrices
      sp_cmd->AddFlag( m_InputFilenames);
      // flag for appending or merging input matrices
      sp_cmd->AddFlag( m_AppendMerge);
      // flag for using model storage to find input files and deciding which to merge
      sp_cmd->AddFlag( m_ModelStorage);
      // filename for output matrix
      sp_cmd->AddFlag( m_OutputFilename);
      // flag for computing median prediction
      sp_cmd->AddFlag( m_Median);
      // flag for computing jury-based prediction
      sp_cmd->AddFlag( m_Jury);
      // flag for computing min prediction
      sp_cmd->AddFlag( m_Min);
      // flag for computing max prediction
      sp_cmd->AddFlag( m_Max);
      // flag for computing max prediction
      sp_cmd->AddFlag( m_LocalPPV);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int ModelPredictionMerge::Main() const
    {
      int count_merge_type_flags( m_Min->GetFlag() + m_Max->GetFlag() + m_Jury->GetFlag() + m_Median->GetFlag());
      BCL_Assert( count_merge_type_flags <= 1, "At most, 1 of -min, -max, -jury, or -median may be given");

      int count_source_flags( int( m_InputFilenames->GetFlag()) + int( m_ModelStorage->GetFlag()));
      BCL_Assert
      (
        count_source_flags == 1,
        "Either a model storage (-input_model_storage) or input filenames (-input_filenames) must be given"
      );

      // map from cross-validation ids to filenames
      storage::Map< std::string, storage::Vector< model::CrossValidationInfo> > cv_filenames;
      if( m_InputFilenames->GetFlag())
      {
        if( m_AppendMerge->GetFlag())
        {
          // all files will be treated as if from a different cv
          storage::Vector< std::string> filenames( m_InputFilenames->GetStringList());
          for( size_t i( 0), n_files( filenames.GetSize()); i < n_files; ++i)
          {
            cv_filenames[ util::Format()( i)].PushBack
            (
              model::CrossValidationInfo().SetIndependentPredictionsFilename( filenames( i))
            );
          }
        }
        else
        {
          storage::Vector< std::string> filenames( m_InputFilenames->GetStringList());
          for( size_t i( 0), n_files( filenames.GetSize()); i < n_files; ++i)
          {
            cv_filenames[ ""].PushBack
            (
              model::CrossValidationInfo().SetIndependentPredictionsFilename( filenames( i))
            );
          }
        }
      }
      else
      {
        // filenames need to be retrieved from the model storage
        util::Implementation< model::RetrieveInterface> model_retriever( m_ModelStorage->GetFirstParameter()->GetValue());
        storage::List< model::CrossValidationInfo> cv_info( model_retriever->RetrieveEnsembleCVInfo());
        storage::Vector< size_t> counts( cv_info.GetSize(), size_t( 0));
        for
        (
          storage::List< model::CrossValidationInfo>::const_iterator itr( cv_info.Begin()), itr_end( cv_info.End());
          itr != itr_end;
          ++itr
        )
        {
          storage::Vector< model::CrossValidationInfo> &vec
          (
            cv_filenames[ itr->GetIndependentDatasetRetriever().ToString()]
          );
          // track the # of times each cv independent id is seen
          ++counts( vec.GetSize());
          vec.PushBack( *itr);
        }
        // test whether all cv's have the same # of files
        const size_t count
        (
          std::count( counts.Begin(), counts.End(), size_t( 0)) + std::count( counts.Begin(), counts.End(), counts( 0))
        );
        if( count != counts.GetSize())
        {
          storage::Map< std::string, storage::Vector< model::CrossValidationInfo> > cv_filenames_diff;
          counts.SetAllElements( 0);
          for
          (
            storage::List< model::CrossValidationInfo>::const_iterator itr( cv_info.Begin()), itr_end( cv_info.End());
            itr != itr_end;
            ++itr
          )
          {
            // use the difference function to compute the difference of the independent retriever from the training retriever
            // The reason this is needed is because the underlying training dataset files may have been computed from
            // identical files at different locations (e.g. using scratch on a cluster vs. an explicit path when done locally)
            storage::Vector< model::CrossValidationInfo> &vec
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
            std::count( counts.Begin(), counts.End(), size_t( 0)) + std::count( counts.Begin(), counts.End(), counts( 0))
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

        // set the jury objective function, if it was not already given
        if( ( m_Jury->GetFlag() || m_LocalPPV->GetFlag()) && !m_Jury->GetFirstParameter()->GetWasSetInCommandLine())
        {
          m_Jury->GetParameterList()( 0)->SetParameter
          (
            cv_info.FirstElement().GetObjective().ToString(),
            util::GetLogger()
          );
        }
      }

      // final container of prediction matrices
      storage::Vector< linal::Matrix< float> > final_prediction_matrices;
      storage::Vector< linal::Matrix< char> > final_id_matrices;
      storage::Vector< linal::Matrix< float> > final_result_matrices;

      model::CrossValidationInfo common_cv_info;
      bool first_in_all_cv( true);
      // iterate over all input file names and construct prediction matrices for each entry
      for
      (
        storage::Map< std::string, storage::Vector< model::CrossValidationInfo> >::iterator
          itr_independent_to_filenames( cv_filenames.Begin()), itr_independent_to_filenames_end( cv_filenames.End());
        itr_independent_to_filenames != itr_independent_to_filenames_end;
        ++itr_independent_to_filenames
      )
      {
        // container for prediction matrices
        storage::Vector< linal::Matrix< float> > prediction_matrices;
        model::FeatureDataSet< char> ids_matrix;
        model::FeatureDataSet< float> results_matrix;

        math::RunningAverage< linal::Matrix< float> > avg_prediction_matrices;
        math::RunningMinMax< linal::Matrix< float> > minmax_prediction_matrices;
        bool is_first( true);
        for
        (
          storage::Vector< model::CrossValidationInfo>::iterator
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
            BCL_Assert
            (
              linal::EqualIgnoringNans( results_matrix.GetMatrix(), sp_dataset->GetResultsReference()),
              "Incorrect results matrix from " + itr_cv_info->GetIndependentPredictionsFilename()
            );
          }

          linal::MatrixConstReference< float> prediction( sp_dataset->GetFeaturesReference());
          if( m_AppendMerge->GetFlag())
          {
            final_prediction_matrices.PushBack( prediction);
          }
          else if( m_Median->GetFlag() || m_Jury->GetFlag())
          {
            prediction_matrices.PushBack( prediction);
            BCL_Assert( prediction.GetNumberCols() == prediction_matrices( 0).GetNumberCols(), "Inconsistent # columns");
            BCL_Assert( prediction.GetNumberRows() == prediction_matrices( 0).GetNumberRows(), "Inconsistent # rows");
          }
          else if( m_Min->GetFlag() || m_Max->GetFlag())
          {
            minmax_prediction_matrices += prediction;
          }
          else
          {
            BCL_MessageDbg( "Merge ..");
            avg_prediction_matrices += prediction;
          }
        }
        // append result for this cross validation to the final_prediction_matrices vector
        if( m_Median->GetFlag())
        {
          final_prediction_matrices.PushBack( ComputeMedian( prediction_matrices));
        }
        else if( m_Jury->GetFlag())
        {
          final_prediction_matrices.PushBack( ComputeJury( prediction_matrices, results_matrix, ids_matrix));
        }
        else if( m_Min->GetFlag())
        {
          final_prediction_matrices.PushBack( minmax_prediction_matrices.GetMin());
        }
        else if( m_Max->GetFlag())
        {
          final_prediction_matrices.PushBack( minmax_prediction_matrices.GetMax());
        }
        else if( avg_prediction_matrices.GetWeight())
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

      // if doing local ppv, transform accordingly.
      if( m_LocalPPV->GetFlag())
      {
        util::Implementation< model::ObjectiveFunctionInterface> obj( m_Jury->GetParameterList()( 0)->GetValue());
        storage::Vector< math::ROCCurve> roc_curves
        (
          model::RetrieveInterface::ROCCurvesFromDataset
          (
            combined_result_matrix,
            combined_prediction_matrix,
            obj->GetThreshold(),
            obj->GetRankingParity()
          )
        );
        storage::Vector< math::PiecewiseFunction> localppvs( roc_curves.GetSize());
        for( size_t result( 0), n_results( localppvs.GetSize()); result < n_results; ++result)
        {
          localppvs( result) = roc_curves( result).GetLocalPPVCurve();
        }
        for( size_t row( 0), n_rows( combined_prediction_matrix.GetNumberRows()); row < n_rows; ++row)
        {
          for( size_t result( 0), n_results( localppvs.GetSize()); result < n_results; ++result)
          {
            combined_prediction_matrix( row, result) = localppvs( result)( combined_prediction_matrix( row, result));
          }
        }
      }
      if( common_cv_info.GetIDsFeatureLabelSet().GetSize())
      {
        io::OFStream output;
        io::File::MustOpenOFStream( output, m_OutputFilename->GetFirstParameter()->GetValue() + ".info");
        common_cv_info.SetIndependentPredictionsFilename( m_OutputFilename->GetFirstParameter()->GetValue());
        output << common_cv_info;
        io::File::CloseClearFStream( output);
      }

      // output file
      common_cv_info.WritePredictions
      (
        combined_result_matrix,
        combined_prediction_matrix,
        combined_id_matrix,
        m_OutputFilename->GetFirstParameter()->GetValue()
      );

      BCL_MessageStd
      (
        "Final prediction matrix written in "
        + m_OutputFilename->GetFirstParameter()->GetValue()
      );

      // end
      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ModelPredictionMerge::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ModelPredictionMerge::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief compute median for each result column, output into a matrix
    //! @param RESULTS vector of predictions, all corresponding to the same dataset
    //! @return median prediction values for each result column
    linal::Matrix< float> ModelPredictionMerge::ComputeMedian( const storage::Vector< linal::Matrix< float> > &RESULTS) const
    {
      // compute median value
      const size_t dataset_size( RESULTS( 0).GetNumberRows());
      const size_t n_results( RESULTS( 0).GetNumberCols());
      const size_t n_models( RESULTS.GetSize());
      linal::Matrix< float> combined_matrices( dataset_size, n_results);
      storage::Vector< float> models_values( n_models);
      const size_t median_high_value( ( n_models + 1) / 2);
      const size_t median_low_value( median_high_value - !( n_models & size_t( 1)));

      // determine the median
      for( size_t datum( 0); datum < dataset_size; ++datum)
      {
        for( size_t result( 0); result < n_results; ++result)
        {
          for( size_t model( 0); model < n_models; ++model)
          {
            models_values( model) = RESULTS( model)( datum, result);
          }
          models_values.Sort( std::less< float>());
          combined_matrices( datum, result)
            = 0.5 * ( models_values( median_high_value) + models_values( median_low_value));
        }
      }
      return combined_matrices;
    }

    //! @brief compute jury for each result column, output into a matrix
    //! @param PREDICTIONS vector of predictions, all corresponding to the same dataset
    //! @param RESULTS actual results
    //! @param RESULT_IDS ids for each result, needed by some objective functions
    //! @return jury prediction values for each result column
    linal::Matrix< float> ModelPredictionMerge::ComputeJury
    (
      const storage::Vector< linal::Matrix< float> > &PREDICTIONS,
      const model::FeatureDataSet< float> &RESULTS,
      const model::FeatureDataSet< char> &RESULT_IDS
    ) const
    {
      const size_t dataset_size( RESULTS.GetNumberFeatures());
      const size_t n_results( RESULTS.GetFeatureSize());
      const size_t n_models( PREDICTIONS.GetSize());
      util::Implementation< model::ObjectiveFunctionInterface> obj_function( m_Jury->GetFirstParameter()->GetValue());
      storage::Vector< linal::Matrix< char> > prediction_classifications( n_models);
      obj_function->SetData( RESULTS, RESULT_IDS);
      // classify each prediction
      for( size_t model( 0); model < n_models; ++model)
      {
        model::FeatureDataReference< float> reference( PREDICTIONS( model), RESULTS.GetFeatureLabelSet());
        prediction_classifications( model) = obj_function->GetFeaturePredictionClassifications( RESULTS, reference);
      }
      storage::Vector< size_t> counts( size_t( 5), size_t( 0));
      linal::Matrix< float> combined_matrices( dataset_size, n_results);
      for( size_t datum( 0); datum < dataset_size; ++datum)
      {
        for( size_t result( 0); result < n_results; ++result)
        {
          counts.SetAllElements( 0);
          for( size_t model( 0); model < n_models; ++model)
          {
            switch( prediction_classifications( model)( datum, result))
            {
              case '\0': ++counts( 0); break;
              case 'P': ++counts( 1); break;
              case 'N': ++counts( 2); break;
              case 'p': ++counts( 3); break;
              case 'n': ++counts( 4); break;
              default: break;
            }
          }
          // determine which prediction letters are most popular
          const size_t target_count( *std::max_element( counts.Begin(), counts.End()));
          std::string target_chars;
          if( counts( 0) == target_count)
          {
            target_chars += '\0';
          }
          if( counts( 1) == target_count)
          {
            target_chars += 'P';
          }
          if( counts( 2) == target_count)
          {
            target_chars += 'N';
          }
          if( counts( 3) == target_count)
          {
            target_chars += 'p';
          }
          if( counts( 4) == target_count)
          {
            target_chars += 'n';
          }
          math::RunningAverage< float> ave;
          for( size_t model( 0); model < n_models; ++model)
          {
            if( target_chars.find( prediction_classifications( model)( datum, result)) != std::string::npos)
            {
              ave += PREDICTIONS( model)( datum, result);
            }
          }
          combined_matrices( datum, result) = ave.GetAverage();
        }
      }
      return combined_matrices;
    }

    //! @brief standard constructor
    ModelPredictionMerge::ModelPredictionMerge() :
      m_InputFilenames
      (
        new command::FlagDynamic
        (
          "input",
          "input prediction matrices",
          command::Parameter // the template parameter
          (
            "input",
            "input prediction matrix"
          ),
          0,  // min # parameters for this flag
          50  // max # parameters for this flag
        )
      ),
      m_AppendMerge
      (
        new command::FlagStatic
        (
          "modus",
          "flag for appending (flag is set) or merging (flag is not set) prediction matrices"
        )
      ),
      m_ModelStorage
      (
        new command::FlagStatic
        (
          "input_model_storage",
          "Use a model storage to compute all the input filenames, which to merge, which to append, etc., rather than"
          "multiple calls to this app with m_InputFilenames",
          command::Parameter
          (
            "model_storage",
            "",
            command::ParameterCheckSerializable( util::Implementation< model::RetrieveInterface>()),
            ""
          )
        )
      ),
      m_OutputFilename
      (
        new command::FlagStatic
        (
          "output",
          "\toutput filename for final prediction matrix",
          command::Parameter
          (
            "output",
            "output filename for final prediction matrix"
          )
        )
      ),
      m_Median
      (
        new command::FlagStatic
        (
          "median",
          "flag computing median prediction prediction matrices"
        )
      ),
      m_Jury
      (
        new command::FlagStatic
        (
          "jury",
          "flag computing jury prediction based on an objective function",
          command::Parameter
          (
            "obj_function",
            "objective function to use for computation of positives vs negatives. The final objective function for model"
            " training will be used if this isn't specified and -input_model_storage was given",
            command::ParameterCheckSerializable( util::Implementation< model::ObjectiveFunctionInterface>()),
            ""
          )
        )
      ),
      m_Min
      (
        new command::FlagStatic
        (
          "min",
          "flag computing min prediction matrices"
        )
      ),
      m_Max
      (
        new command::FlagStatic
        (
          "max",
          "flag computing max prediction matrices"
        )
      ),
      m_LocalPPV
      (
        new command::FlagStatic
        (
          "local_ppv",
          "average local PPV values instead of raw model outputs",
          command::Parameter
          (
            "obj_function",
            "objective function to use for computation of positives vs negatives. The final objective function for model"
            " training will be used if this isn't specified and -input_model_storage was given",
            command::ParameterCheckSerializable( util::Implementation< model::ObjectiveFunctionInterface>()),
            ""
          )
        )
      )
    {
    }

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ModelPredictionMerge::GetReadMe() const
    {
      static const std::string readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL model:PredictionMerge, terms of use, appropriate citation, installation "
        "procedures, BCL model:PredictionMerge execution, technical support, and future research directions.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL model:PredictionMerge?\n"
        "BCL model:PredictionMerge is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Commons. BCL model:PredictionMerge reads in multiple prediction outputs"
        " of a trained machine learning models and merges or appends these into one prediction outputs. The functionality is beneficial "
        " for cross-validation use case scenarios. This application is used in conjunction with BCL model:Train and model:ComputeStatistics."
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL model:PredictionMerge.\n"
        "When using BCL model:PredictionMerge in a publication, please cite the following publication describing the application's "
        "development:\n"
        "\n"
        "Butkiewicz, M.; Lowe, E.W., Jr.; Mueller, R.; Mendenhall, J.L.; Teixeira, P.L.; Weaver, C.D.; Meiler, J.\n"
        "'Benchmarking Ligand-Based Virtual High-Throughput Screening with the PubChem Database'. Molecules 2013, 18, 735-756.\n"
        "\n"
        + DefaultSectionSeparator() +
        "V. INSTALLATION PROCEDURES.\n"
        + DefaultInstallationProcedure() +
        "\n"
        + DefaultSectionSeparator() +
        "VI. RUNNING BCL model:PredictionMerge.\n"
        "Running BCL model:PredictionMerge consists of the following steps.\n"
        "\n"
        "1) At a command prompt, navigate to the location of your BCL model:PredictionMerge executable program.\n"
        "\n"
        "2) Run BCL model:PredictionMerge on one or multiple <experimental predicted value file>s\n"
        "\n"
        "3) Obtain one merged <experimental predicted value file> which can be used as input for BCL model:ComputeStatistics.\n"
        "\n"
        "ADDITIONAL INFORMATION\n"
        "\n"
        "Information concerning syntax and flags  can be obtained by typing: <bcl.exe>  model:PredictionMerge -help\n"
        "\n"
        "For further explanation, examples of the flags, example command lines, input and output format information,\n"
        "and example output please see the documentation file.\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL model:PredictionMerge.\n"
        "BCL model:PredictionMerge is under ongoing further development.\n"
        "\n"
        + DefaultSectionSeparator() +
        "IX. EXAMPLE COMMANDLINE\n"
        "Example 1: Using csv or BCL-linal::Matrix style files\n"
        "<bcl.exe> model:PredictionMerge -input <experimental predicted values file> <experimental predicted values file> "
        "-output <new merged experimental predicted values file>\n"
        "\n"
        "The <experimental predicted values file> can be in simple csv format, e.g.:\n"
        "experimental_value1,predicted_value1\n"
        "experimental_value2,predicted_value2\n"
        "or for multiple outputs:\n"
        "experimental_value_1_1,experimental_value_1_2,predicted_value_1_1,predicted_value_1_2\n"
        "experimental_value_2_1,experimental_value_2_2,predicted_value_2_1,predicted_value_2_2\n"
        "Additionally, one fixed width ID column is allowed, so a file could look like:\n"
        "1ABC D F,0,1,0.3,0.9\n"
        "1ABC X F,0,1,0.25,0.8\n"
        "\nAlternatively, the BCL Matrix format can be used for some or all of the input files:\n"
        "\n"
        "bcl::linal::Matrix<float>\n"
        "number_rows number_columns\n"
        "experimental_value1 experimental_value2 predicted_value1 predicted_value2 ...\n"
        "experimental_value1 experimental_value2 predicted_value1 predicted_value2 ...\n"
        "...\n"
        "\n"
        "The <experimental predicted values file> of a model with only one output (eg. predicted biological activity)\n"
        "would have the following format:\n"
        "\n"
        "bcl::linal::Matrix<float>\n"
        "number_rows number_columns\n"
        "experimental_value1 predicted_value1\n"
        "experimental_value2 predicted_value2\n"
        "...\n"
        "Example 2: Integration with model:Train\n"
        "Using the commands in the example above requires that N+1 runs to model:PredictionMerge in an "
        "N X M-fold cross-validation (once for each independent dataset plus once to concatenate all the predictions.)\n"
        "Often the input files have been written using model:Train -print_independent_predictions and the "
        "models themselves have been written out (using -storage_model), then only a single command is needed: \n"
        "<bcl.exe> model:PredictionMerge -input_model_storage <model-storage> -output <merged file> -jury\n"
        " <model-storage> should be the same string as that used for -storage_model in model:Train, e.g. "
        "'File(directory=models/test_run,prefix=\"model\")'\n"
        " <merged_file> will be the output file, containing all the merged predictions\n"
        " Note that for this command-line, if the objective function is not specified after -jury, model:PredictionMerge "
        " automatically uses the same objective function used during model:Train\n"
        "\n"
        + DefaultSectionSeparator()
      );

      return readme;
    }

    //! @brief returns web text information
    //! @return text (html allowed but not required) that will be displayed on the website
    //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
    const std::string &ModelPredictionMerge::GetWebText() const
    {
      // create a static string to hold readme information
      static const std::string s_web_text
      (
        "BCL::ModelPredictionMerge: Merging of Machine Learning Predictions\n\n"
        "As the part of the application suite BCL::ChemInfo, BCL::ModelPredictionMerge allows for handling of "
        "prediction results generated by machine learning algorithms. It is an helper application that provides "
        "functionality to post-process prediction results from BCL::ModelTrain."
        "\n"
        "Necessary for NxM cross-validation calculations, this applications provides functionality to append and merge"
        "prediction matrices, given a list of matrix files, csv files, or a model storage\n"
        "!bcl_model_prediction_merge.png!\n"
        "\nFig. 1:\n"
        "Schematic of the functionality provided by BCL::ModelPreditionMerge.\n"
        ""
      );

      return s_web_text;
    }

    // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
    // but which should not be displayed, e.g. for help
    storage::Vector< std::string> ModelPredictionMerge::GetDeprecatedAppNames() const
    {
      return storage::Vector< std::string>( size_t( 1), "TrainModelPredictionMerge");
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string ModelPredictionMerge::GetDescription() const
    {
      return "Merge cross-validation/prediction matrices horizontally or simply append";
    }

    const ApplicationType ModelPredictionMerge::ModelPredictionMerge_Instance
    (
      GetAppGroups().AddAppToGroup( new ModelPredictionMerge(), GetAppGroups().e_Model)
    );

  } // namespace app
} // namespace bcl
