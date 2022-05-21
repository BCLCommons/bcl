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
#include "bcl_app_model_compute_statistics.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "command/bcl_command_parameter_check_serializable.h"
#include "descriptor/bcl_descriptor_dataset.h"
#include "descriptor/bcl_descriptor_example_string_sequence.h"
#include "descriptor/bcl_descriptor_iterator.h"
#include "io/bcl_io_directory_entry.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_contingency_matrix_measures.h"
#include "math/bcl_math_roc_curve.h"
#include "math/bcl_math_running_average.h"
#include "math/bcl_math_statistics.h"
#include "math/bcl_math_template_instantiations.h"
#include "model/bcl_model_cross_validation_info.h"
#include "model/bcl_model_objective_function_enrichment_average.h"
#include "model/bcl_model_objective_function_integral_precision_fraction_predicted.h"
#include "storage/bcl_storage_table.hpp"
#include "util/bcl_util_implementation.h"

namespace bcl
{
  namespace app
  {
    //! @brief Clone function
    //! @return pointer to new FoldProtein
    ModelComputeStatistics *ModelComputeStatistics::Clone() const
    {
      return new ModelComputeStatistics( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &ModelComputeStatistics::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> ModelComputeStatistics::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      sp_cmd->AddFlag( m_InputFilenamesFlag);
      sp_cmd->AddFlag( m_MaxFilesPerConsensusFlag);
      sp_cmd->AddFlag( m_TableNameFlag);
      sp_cmd->AddFlag( m_SortByFlag);
      sp_cmd->AddFlag( m_PlotXFlag);
      sp_cmd->AddFlag( m_PlotLogXFlag);
      sp_cmd->AddFlag( m_PlotYFlag);
      sp_cmd->AddFlag( m_NoPlotFlag);
      sp_cmd->AddFlag( m_TakeLog10Flag);
      sp_cmd->AddFlag( m_CorrelationFlag);
      sp_cmd->AddFlag( m_PotencyCutOffFlag);
      sp_cmd->AddFlag( m_ActivesBelowCutoffFlag);
      sp_cmd->AddFlag( m_OutputDirectory);
      sp_cmd->AddFlag( m_ObjectiveFunction);
      sp_cmd->AddFlag( m_EvaluateObjectiveFunctionFilename);
      sp_cmd->AddFlag( m_ImageFormatFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int ModelComputeStatistics::Main() const
    {
      if( !m_NoPlotFlag->GetFlag())
      {
        // test whether gnuplot is available
        m_GnuplotAvailable = system( "gnuplot -V") == 0;
      }
      else
      {
        m_GnuplotAvailable = false;
      }

      // if the objective function flag was given, see if other flags should be set too
      if( m_ObjectiveFunction->GetFlag())
      {
        // get the first objective function
        util::Implementation< model::ObjectiveFunctionInterface> first_obj
        (
          m_ObjectiveFunction->GetFirstParameter()->GetValue()
        );

        // set the potency cutoff if the objective function has one and the flag was not set
        if( !m_PotencyCutOffFlag->GetFlag() && util::IsDefined( first_obj->GetThreshold()))
        {
          m_PotencyCutOffFlag->GetParameterList().FirstElement()->SetParameter
          (
            util::Format()( first_obj->GetThreshold()), util::GetLogger()
          );
          m_PotencyCutOffFlag->SetFlag();
        }

        // set ranking parity
        if( !m_ActivesBelowCutoffFlag->GetFlag() && !first_obj->GetRankingParity())
        {
          m_ActivesBelowCutoffFlag->SetFlag();
        }

        if( !m_NoPlotFlag->GetFlag())
        {
          // check whether any plotting flags were given
          if( !m_CorrelationFlag->GetFlag() && !m_PlotXFlag->GetFlag() && !m_PlotLogXFlag->GetFlag() && !m_PlotYFlag->GetFlag())
          {
            // set -correlation if no other plotting flags were given and the objective function is regression
            if( first_obj->GetGoalType() == model::ObjectiveFunctionInterface::e_Regression)
            {
              m_CorrelationFlag->SetFlag();
            }
            else
            {
              // setup the plot-x function
              m_PlotXFlag->ReadFromList( storage::Vector< std::string>::Create( "FPR", "Cutoff"), util::GetLogger());
              m_PlotXFlag->SetFlag();
              m_PlotLogXFlag->ReadFromList( storage::Vector< std::string>::Create( "FPR"), util::GetLogger());
              m_PlotLogXFlag->SetFlag();
            }
          }
          // check whether an plot x was given, but no plot y was given
          if( m_PlotXFlag->GetFlag() && !m_PlotYFlag->GetFlag())
          {
            if( first_obj->GetGoalType() == model::ObjectiveFunctionInterface::e_Classification)
            {
              m_PlotYFlag->ReadFromList( storage::Vector< std::string>::Create( "Accuracy", "TPR", "TNR", "MCC"), util::GetLogger());
            }
            else if( first_obj->GetGoalType() == model::ObjectiveFunctionInterface::e_RankClassification)
            {
              m_PlotYFlag->ReadFromList
              (
                storage::Vector< std::string>::Create( "InformationGainRatio", "TPR", "PPV", "LocalPPV", "Ideal-PPV", "Ideal-PPV_FPRelative"),
                util::GetLogger()
              );
            }
          }
        }

        if( !m_OutputDirectory->GetFlag())
        {
          const std::string output_dir
          (
            io::File::SplitToPathAndFileName( m_InputFilenamesFlag->GetFirstParameter()->GetValue())( 0)
          );
          if( output_dir.size())
          {
            BCL_MessageStd( "Output path set to" + util::Format()( output_dir));
            m_OutputDirectory->GetParameterList()( 0)->SetParameter( output_dir, util::GetLogger());
            m_OutputDirectory->SetFlag();
          }
        }
      }

      // validation
      BCL_Assert
      (
        !m_PlotYFlag->GetFlag() || m_PlotXFlag->GetFlag() || m_PlotLogXFlag->GetFlag(),
        "-plot_y cannot be given without -plot_x or -plot_log_x!"
      );

      // ostreams
      io::IFStream input;

      // vector of predictions and experimental values
      storage::Vector< storage::Pair< std::string, linal::Matrix< float> > > data_sets;
      linal::Matrix< float> data;

      // file names
      const storage::Vector< std::string> files( m_InputFilenamesFlag->GetStringList());

      // load in predictions files
      size_t number_result_cols( 0), nr_points( 0);
      linal::Matrix< float> experimental, experimental_transposed;
      model::FeatureDataSet< char> ids;
      const storage::TableHeader table_header( GetTableHeader( m_CorrelationFlag->GetFlag()));
      storage::Table< float> table( table_header);

      if( m_ObjectiveFunction->GetFlag() && m_EvaluateObjectiveFunctionFilename->GetFlag())
      {
        // remove any existing file, since the file may be appended to multiple times
        io::DirectoryEntry entry( m_EvaluateObjectiveFunctionFilename->GetFirstParameter()->GetValue());
        if( entry.DoesExist())
        {
          entry.Remove();
        }
      }

      for
      (
        storage::Vector< std::string>::const_iterator itr_files( files.Begin()), itr_files_end( files.End());
        itr_files != itr_files_end;
        ++itr_files
      )
      {
        util::ShPtr< descriptor::Dataset> sp_dataset( model::CrossValidationInfo::ReadPredictions( *itr_files));
        BCL_Assert( sp_dataset.IsDefined(), "Undefined dataset at " + *itr_files);
        data_sets.PushBack();
        storage::Pair< std::string, linal::Matrix< float> > &filename_result( data_sets.LastElement());
        if( files.GetSize() == size_t( 1))
        {
          filename_result.First() = io::File::RemovePath( *itr_files);
        }
        else if( itr_files != files.Begin())
        {
          filename_result.First() = io::File::RemoveCommonComponents( *itr_files, files.FirstElement());
        }
        else
        {
          filename_result.First() = io::File::RemoveCommonComponents( *itr_files, files.LastElement());
        }
        filename_result.Second() = sp_dataset->GetFeaturesReference();
        const size_t new_number_result_cols( filename_result.Second().GetNumberCols());
        const size_t new_nr_points( filename_result.Second().GetNumberRows());
        if( number_result_cols && number_result_cols != new_number_result_cols)
        {
          BCL_Exit( "Prediction files contained different #s of result columns", -1);
        }
        number_result_cols = new_number_result_cols;
        if( nr_points && nr_points != new_nr_points)
        {
          BCL_MessageCrt
          (
            "Tried combining result files that have results in different #s of points! Files: "
            + *itr_files + ", " + files.FirstElement() + "; skipping"
          );
          data_sets.PopBack();
          continue;
        }
        filename_result.Second().Transpose();
        if( !nr_points)
        {
          // copy the first experimental columns
          experimental_transposed = experimental = sp_dataset->GetResultsReference();
          experimental.Transpose();
          nr_points = new_nr_points;
          ids = sp_dataset->GetIds();
          if( m_TakeLog10Flag->GetFlag())
          {
            for( size_t row( 0); row < nr_points; ++row)
            {
              for( size_t col( 0); col < number_result_cols; ++col)
              {
                experimental( row, col) = log10( experimental( row, col));
              }
            }
          }
        }
        else
        {
          if( sp_dataset->GetResultsReference() != experimental_transposed || sp_dataset->GetIdsReference() != ids.GetMatrix())
          {
            BCL_MessageCrt
            (
              "Tried combining result files that have results in different orders! Files: "
              + *itr_files + ", " + files.FirstElement() + "; skipping"
            );
            data_sets.PopBack();
            continue;
          }
        }
        linal::Matrix< float> &mat( data_sets.LastElement().Second());
        if( m_TakeLog10Flag->GetFlag())
        {
          for( size_t row( 0); row < nr_points; ++row)
          {
            for( size_t col( 0); col < number_result_cols; ++col)
            {
              mat( row, col) = log10( mat( row, col));
            }
          }
        }
      }

      const size_t max_files_at_a_time
      (
        std::min( m_MaxFilesPerConsensusFlag->GetFirstParameter()->GetNumericalValue< size_t>(), data_sets.GetSize())
      );
      model::FeatureDataSet< float> exp_fds;
      if( m_ObjectiveFunction->GetFlag())
      {
        m_ObjectiveFunctions.Reset();
        m_ObjectiveFunctions.AllocateMemory( m_ObjectiveFunction->GetSize());

        // construct FeatureDataSets for given experimental values
        exp_fds = model::FeatureDataSet< float>( experimental_transposed);
        for
        (
          util::ShPtrVector< command::ParameterInterface>::const_iterator
            itr( m_ObjectiveFunction->GetParameterList().Begin()), itr_end( m_ObjectiveFunction->GetParameterList().End());
          itr != itr_end;
          ++itr
        )
        {
          m_ObjectiveFunctions.PushBack( util::Implementation< model::ObjectiveFunctionInterface>( ( *itr)->GetValue()));
          m_ObjectiveFunctions.LastElement()->SetData( exp_fds, ids);
        }
      }

      if( max_files_at_a_time >= data_sets.GetSize())
      {
        for( size_t number( 1), final( 1 << data_sets.GetSize()); number < final; ++number)
        {
          math::RunningAverage< linal::Matrix< float> > mat_ave;
          std::string concat_filenames;
          size_t tmp( number);
          for( size_t count( 0); tmp != 0; ++count)
          {
            if( tmp & 1)
            {
              mat_ave += data_sets( count).Second();
              if( !concat_filenames.empty())
              {
                concat_filenames += '+';
              }
              concat_filenames += data_sets( count).First();
            }
            tmp >>= 1;
          }
          if( m_ObjectiveFunction->GetFlag())
          {
            EvaluateObjectiveFunctions( mat_ave.GetAverage().Transposed(), exp_fds, concat_filenames);
            util::GetLogger().LogStatus( util::Format()( 100.0 * float( number + 1) / float( final)) + " % complete");
          }
          if( !m_ObjectiveFunction->GetFlag() || m_PlotXFlag->GetFlag() || m_CorrelationFlag->GetFlag())
          {
            CreateDesiredPlots( table, mat_ave.GetAverage(), experimental, concat_filenames);
          }
        }
      }
      else
      {
        // use the descriptor iterator, which supports multi-dimensional iteration over unique elements
        descriptor::StringSequence sequence( std::string( data_sets.GetSize(), ' '));
        descriptor::Iterator< char> itr
        (
          descriptor::Type( max_files_at_a_time, true, descriptor::Type::e_Symmetric),
          sequence
        );
        BCL_MessageStd( "# of combinations: " + util::Format()( itr.GetSize()));
        for( ; itr.NotAtEnd(); ++itr)
        {
          math::RunningAverage< linal::Matrix< float> > mat_ave;
          std::string concat_filenames;
          size_t last_position( util::GetUndefined< size_t>());
          for( size_t i( 0); i < max_files_at_a_time; ++i)
          {
            size_t next_position( itr( i).GetPosition());
            if( next_position != last_position)
            {
              mat_ave += data_sets( next_position).Second();
              if( !concat_filenames.empty())
              {
                concat_filenames += '+';
              }
              concat_filenames += data_sets( next_position).First();
              last_position = next_position;
            }
          }
          if( m_ObjectiveFunction->GetFlag())
          {
            EvaluateObjectiveFunctions( mat_ave.GetAverage().Transposed(), exp_fds, concat_filenames);
            util::GetLogger().LogStatus( util::Format()( 100.0 * float( itr.GetPosition() + 1) / float( itr.GetSize())) + " % complete");
          }
          if( !m_ObjectiveFunction->GetFlag() || m_PlotXFlag->GetFlag() || m_CorrelationFlag->GetFlag())
          {
            CreateDesiredPlots( table, mat_ave.GetAverage(), experimental, concat_filenames);
          }
        }
      }

      if( m_ObjectiveFunction->GetFlag())
      {
        return 0;
      }

      std::string sort_by_this( m_SortByFlag->GetFirstParameter()->GetValue());
      table.SortByColumn( sort_by_this);

      std::string table_name
      (
        m_OutputDirectory->GetFirstParameter()->GetValue() + "/" + m_TableNameFlag->GetFirstParameter()->GetValue()
      );

      io::OFStream table_out;
      io::File::MustOpenOFStream( table_out, table_name);
      table.WriteFormatted( table_out);
      io::File::CloseClearFStream( table_out);

      // end
      return 0;
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ModelComputeStatistics::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return output stream which was written to
    std::ostream &ModelComputeStatistics::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

    //! @brief create the raw data output to create a file that can be parsed by gnuplot
    //! @param TABLE table of interest that will be filled up with results data
    //! @param PRED predicted data matrix
    //! @param EXP experimental data matrix
    //! @param NAMES base name to create appropriate output filenames
    void ModelComputeStatistics::CreateDesiredPlots
    (
      storage::Table< float> &TABLE,
      const linal::Matrix< float> &PRED,
      const linal::Matrix< float> &EXP,
      const std::string &NAMES
    ) const
    {
      if( m_NoPlotFlag->GetFlag())
      {
        return;
      }
      const size_t number_result_cols( PRED.GetNumberRows());
      if( number_result_cols == size_t( 1))
      {
        // handle the trivial case with one output
        CreateDesiredPlot( TABLE, PRED.GetRow( 0), EXP.GetRow( 0), NAMES);
      }
      else
      {
        // handle multiple outputs by appending to NAMES
        const std::string name_col( NAMES + "_col");
        for( size_t i( 0); i < number_result_cols; ++i)
        {
          CreateDesiredPlot( TABLE, PRED.GetRow( i), EXP.GetRow( i), name_col + util::Format()( i));
        }
      }
    }

    //! @brief create the raw data output to create a file that can be parsed by gnuplot
    //! @param TABLE table of interest that will be filled up with results data
    //! @param PRED predicted data
    //! @param EXP experimental data
    //! @param NAMES base name to create appropriate output filenames
    void ModelComputeStatistics::CreateDesiredPlot
    (
      storage::Table< float> &TABLE,
      const linal::Vector< float> &PRED,
      const linal::Vector< float> &EXP,
      const std::string &NAMES
    ) const
    {
      std::string file_name( m_OutputDirectory->GetFirstParameter()->GetValue() + "/" + NAMES);

      if( m_CorrelationFlag->GetFlag())
      {
        float correlation( math::Statistics::Correlation( EXP.Begin(), EXP.End(), PRED.Begin(), PRED.End()));
        float correlation_spearman( math::Statistics::CorrelationSpearman( EXP.Begin(), EXP.End(), PRED.Begin(), PRED.End()));
        float r2( math::Statistics::RSquared( EXP.Begin(), EXP.End(), PRED.Begin(), PRED.End()));
        float sum( 0);
        for( size_t count( 0); count < PRED.GetSize(); ++count)
        {
          sum += math::Sqr( EXP( count) - PRED( count));
        }

        float rmsd( math::Sqrt( sum / PRED.GetSize()));

        TABLE.InsertRow( file_name, storage::Vector< float>::Create( rmsd, correlation, correlation_spearman, r2));
        file_name.append( ".gnuplot");
        io::OFStream out;
        io::File::MustOpenOFStream( out, file_name);
        const std::string script_name( file_name);
        file_name.append( "_txt");
        WriteCorrelationHeader( out, file_name);
        WriteCorrelationLegend( out, rmsd, r2);
        io::File::CloseClearFStream( out);
        io::File::MustOpenOFStream( out, file_name);
        WriteFormattedCorrelationInput( out, EXP, PRED);
        io::File::CloseClearFStream( out);
        if( m_GnuplotAvailable)
        {
          if( system( ( "gnuplot " + script_name).c_str()))
          {
            BCL_MessageStd( "Failed to run gnuplot " + script_name);
          }
        }
      }
      else
      {
        // set up the potency cutoff to classify ROC curve
        const float potency_cutoff( m_PotencyCutOffFlag->GetFirstParameter()->GetNumericalValue< float>());

        // compute parity
        const bool parity( !m_ActivesBelowCutoffFlag->GetFlag());

        //convert to data type for ROC curve
        storage::List< storage::Pair< double, double> > dataset_pair;

        // loop over data
        for
        (
          linal::Vector< float>::const_iterator itr_pred( PRED.Begin()), itr_pred_end( PRED.End()),
          itr_exp( EXP.Begin());
          itr_pred != itr_pred_end;
          ++itr_pred, ++itr_exp
        )
        {
          // stores data according to predicted, experimental
          storage::Pair< double, double> dataset_list( *itr_pred, *itr_exp);
          dataset_pair.PushBack( dataset_list);
        }

        //BCL_MessageStd( "next roc!");

        // initialize ROC Curve
        math::ROCCurve roc_dataset( dataset_pair, potency_cutoff, !m_ActivesBelowCutoffFlag->GetFlag());

        const float enr_max( 0.05);
        const float enr_step_size( 0.005);

        model::FeatureDataSet< float> exp_fds( linal::Matrix< float>( EXP.GetSize(), 1, EXP.Begin()));
        model::FeatureDataSet< float> pred_fds( linal::Matrix< float>( PRED.GetSize(), 1, PRED.Begin()));

        model::ObjectiveFunctionEnrichmentAverage obj_avg_enr
        (
          potency_cutoff, enr_max, enr_step_size, !m_ActivesBelowCutoffFlag->GetFlag()
        );
        obj_avg_enr.SetData( exp_fds);

        model::ObjectiveFunctionIntegralPrecisionFractionPredicted obj_ppv_fpp_0_1pct_1pct
        (
          potency_cutoff,
          math::Range< double>( 0.001, 0.01),
          parity
        );
        obj_ppv_fpp_0_1pct_1pct.SetData( exp_fds);

        model::ObjectiveFunctionIntegralPrecisionFractionPredicted obj_ppv_fpp_0_1pct_2pct
        (
          potency_cutoff,
          math::Range< double>( 0.001, 0.02),
          parity
        );
        obj_ppv_fpp_0_1pct_2pct.SetData( exp_fds);

        model::ObjectiveFunctionIntegralPrecisionFractionPredicted obj_ppv_fpp_0_1pct_3pct
        (
          potency_cutoff,
          math::Range< double>( 0.001, 0.03),
          parity
        );
        obj_ppv_fpp_0_1pct_3pct.SetData( exp_fds);

        model::ObjectiveFunctionIntegralPrecisionFractionPredicted obj_ppv_fpp_0pct_100pct
        (
          potency_cutoff,
          math::Range< double>( 0, 1),
          parity
        );
        obj_ppv_fpp_0pct_100pct.SetData( exp_fds);

        const float avg_enr( obj_avg_enr( exp_fds, pred_fds));
        const float roc_integral( roc_dataset.Integral());
        const float integral_ppv_fpp_0_1pct_1pct( obj_ppv_fpp_0_1pct_1pct( exp_fds, pred_fds));
        const float integral_ppv_fpp_0_1pct_2pct( obj_ppv_fpp_0_1pct_2pct( exp_fds, pred_fds));
        const float integral_ppv_fpp_0_1pct_3pct( obj_ppv_fpp_0_1pct_3pct( exp_fds, pred_fds));
        const float integral_ppv_fpp_0pct_100pct( obj_ppv_fpp_0pct_100pct( exp_fds, pred_fds));

        float sum( 0);
        for( size_t count( 0); count < PRED.GetSize(); ++count)
        {
          sum += math::Sqr( EXP( count) - PRED( count));
        }

        float rmsd( math::Sqrt( sum / PRED.GetSize()));

        BCL_MessageStd
        (
          "name: " + NAMES +
          " RMSD: " + util::Format()( rmsd) +
          " Enr:" + util::Format()( avg_enr) +
          " AUC: " + util::Format()( roc_integral) +
          " Int_FppVsPPV(0.1%,1%): " + util::Format()( integral_ppv_fpp_0_1pct_1pct) +
          " Int_FppVsPPV(0.1%,2%): " + util::Format()( integral_ppv_fpp_0_1pct_2pct) +
          " Int_FppVsPPV(0.1%,3%): " + util::Format()( integral_ppv_fpp_0_1pct_3pct) +
          " Int_FppVsPPV(0%,100%): " + util::Format()( integral_ppv_fpp_0pct_100pct)
        );

        TABLE.InsertRow
        (
          NAMES,
          storage::Vector< float>::Create
          (
            rmsd,
            avg_enr,
            roc_integral,
            integral_ppv_fpp_0_1pct_1pct,
            integral_ppv_fpp_0_1pct_2pct,
            integral_ppv_fpp_0_1pct_3pct,
            integral_ppv_fpp_0pct_100pct
          )
        );

        //roc_dataset.WriteRatePlottingTable( script, util::Format());

        if( !m_PlotXFlag->GetFlag() && !m_PlotLogXFlag->GetFlag())
        {
          // open output stream
          io::OFStream script;
          file_name.append( ".gnuplot");
          io::File::MustOpenOFStream( script, file_name);
          const std::string script_name( file_name);
          file_name.append( "_txt");
          WriteROCHeader( script, file_name, avg_enr, roc_integral);
          io::File::CloseClearFStream( script);
          io::File::MustOpenOFStream( script, file_name);

          // print plotting table for ppv vs fraction positive predicted values
          roc_dataset.WriteRatePlottingTableGeneric
          (
            script,
            &math::ContingencyMatrix::GetFalsePositiveRate,
            &math::ContingencyMatrix::GetTruePositiveRate,
            math::Range< double>( 0.0, 1.0),
            util::Format(),
            std::string( "FALSE_POSITIVE_RATE\tTRUE_POSITIVE_RATE")
          );

          // insert newline in script
          script << std::endl;

          // print plotting table for ppv vs fraction positive predicted values
          roc_dataset.WriteRatePlottingTableGeneric
          (
            script,
            &math::ContingencyMatrix::GetFalsePositiveRate,
            &math::ContingencyMatrix::GetFractionPredictedPositives,
            math::Range< double>( 0.0, 1.0),
            util::Format(),
            std::string( "FPR\tTPR")
          );

          script << std::endl;

          // print plotting table for ppv vs fraction positive predicted values
          roc_dataset.WriteRatePlottingTableGeneric
          (
            script,
            &math::ContingencyMatrix::GetFractionPredictedPositives,
            &math::ContingencyMatrix::GetPositivePredictiveValue,
            math::Range< double>( 0.0, 1.0),
            util::Format(),
            std::string( "PPV vs FractionPositivePredicted")
          );

          // insert newline in script

          const size_t total_results( roc_dataset.GetNumberResults());
          const double num_actual_positives( roc_dataset.GetNumberActualPositives());
          const double step_size( 1.0 / double( total_results - 1));

          script << std::endl;
          script << "IDEAL PPV vs FractionPositivePredicted" << std::endl;

          // iterate over sorted counts of roc curve and plot ideal curve for PPV vs FractionPositivePredicted
          for
          (
            storage::Vector< math::ROCCurve::Point>::const_iterator
              itr( roc_dataset.GetSortedCounts().Begin()), itr_end( roc_dataset.GetSortedCounts().End());
            itr != itr_end;
            ++itr
          )
          {
            const size_t counts( itr->GetNumberPredictedPositives());
            const double progress( step_size * counts);
            script << progress << '\t' << std::min( 1.0, num_actual_positives / ( counts + 1)) << '\n';
          }
          io::File::CloseClearFStream( script);
          if( m_GnuplotAvailable)
          {
            if( system( ( "gnuplot " + script_name).c_str()))
            {
              BCL_MessageStd( "Failed to run gnuplot " + script_name);
            }
          }
        }
        else
        {
          const std::string script_name( file_name + ".gnuplot");
          const std::string data_file_name( file_name + ".data");

          storage::Vector< std::string> x_measures_from_flag( m_PlotXFlag->GetStringList());
          storage::Vector< std::string> log_x_measures_from_flag( m_PlotLogXFlag->GetStringList());
          if( !log_x_measures_from_flag.IsEmpty())
          {
            for
            (
              storage::Vector< std::string>::const_iterator
                itr( log_x_measures_from_flag.Begin()), itr_end( log_x_measures_from_flag.End());
              itr != itr_end;
              ++itr
            )
            {
              if( x_measures_from_flag.Find( *itr) >= x_measures_from_flag.GetSize())
              {
                x_measures_from_flag.PushBack( *itr);
              }
            }
          }
          storage::Vector< std::string> measures_from_flag( m_PlotYFlag->GetStringList());
          // ensure that the x-axis values are also stored in the data file
          if( !measures_from_flag.IsEmpty())
          {
            for
            (
              storage::Vector< std::string>::const_iterator
                itr( x_measures_from_flag.Begin()), itr_end( x_measures_from_flag.End());
              itr != itr_end;
              ++itr
            )
            {
              if( measures_from_flag.Find( *itr) >= measures_from_flag.GetSize())
              {
                measures_from_flag.PushBack( *itr);
              }
            }
          }

          storage::Vector< std::string> measures
          (
            roc_dataset.WriteRatePlottingTableComplete( data_file_name, measures_from_flag)
          );

          // open output stream
          io::OFStream script;
          io::File::MustOpenOFStream( script, script_name);
          WriteROCCompleteHeader( script, data_file_name, measures);
          io::File::CloseClearFStream( script);
          if( m_GnuplotAvailable)
          {
            if( system( ( "gnuplot " + script_name).c_str()))
            {
              BCL_MessageStd( "Failed to run gnuplot " + script_name);
            }
          }
        }
      }
    }

    //! @brief generate the table header for the final results table
    //! @param CORRELATION flag whether a correlation plot will be considered
    //! @return table header
    storage::TableHeader ModelComputeStatistics::GetTableHeader( bool CORRELATION) const
    {
      if( CORRELATION)
      {
        return storage::TableHeader
        (
          storage::Vector< std::string>::Create( "RMSD", "CorrelationPearson", "CorrelationSpearman", "R^2")
        );
      }
      else
      {
        return storage::TableHeader
        (
          storage::Vector< std::string>::Create
          (
            "RMSD", "Enrichment", "AUC", " Int FppVsPPV (0.1%,1%) ", " Int FppVsPPV (0.1%,2%) ",
            " Int FppVsPPV (0.1%,3%) ", " Int FppVsPPV (0%,100%) "
          )
        );
      }
    }

    //! @brief gets the string for the terminal type to use for gnuplot graphs
    //! @return a string with the gnuplot terminal type
    std::string ModelComputeStatistics::GetGnuplotTerminalType() const
    {
      std::string img_type( m_ImageFormatFlag->GetFirstParameter()->GetValue());
      if( img_type == "svg")
      {
        return std::string( "svg");
      }
      return std::string( "pngcairo");
    }

    //! @brief get the extension for the given image type
    //! @return the three-letter extension for the image output type
    std::string ModelComputeStatistics::GetImageFormatExtension() const
    {
      return m_ImageFormatFlag->GetFirstParameter()->GetValue();
    }

    //! @brief write the correlation plot header for a gnuplot file to output stream
    //! @param OSTREAM output stream of interest
    //! @param OUTPUT_FILENAME output file name OSTREAM will write to
    void ModelComputeStatistics::WriteCorrelationHeader
    (
      io::OFStream &OSTREAM,
      std::string &OUTPUT_FILENAME
    ) const
    {
      // write header to plot ROC curve with gnuplot
      OSTREAM << "set terminal " << GetGnuplotTerminalType() << " enhanced size 500,500 # transparent size 2160,2160\n";
      OSTREAM << "set output \"" << OUTPUT_FILENAME << "_cross_correlation." << GetImageFormatExtension() << "\"\n";
      OSTREAM << "set encoding iso\n";
      OSTREAM << "set title \"Cross Correlation\"\n";
      OSTREAM << "set xlabel \"Experimental\"\n";
      OSTREAM << "set autoscale fix\n";
      OSTREAM << "plot x with lines linetype -1 linewidth 1 linecolor rgb 'black', \\\n";
      OSTREAM << "  \"" << OUTPUT_FILENAME << "\" using 1:2  with points pointsize 0.5 linecolor rgb 'red'\n";
    }

    //! @brief write the correlation plot legend for a gnuplot file to output stream
    //! @param OSTREAM output stream of interest
    //! @param RMSD rmsd value that is written on the plot
    //! @param RSQAURED r^2 value that is written on the plot
    void ModelComputeStatistics::WriteCorrelationLegend
    (
      io::OFStream &OSTREAM,
      double RMSD,
      double RSQAURED
    ) const
    {
      // text for legend in gnuplot
      OSTREAM << "set title \"Cross correlation: RMSD " << RMSD << ' ' << "R^2 " << RSQAURED << "\"";
    }

    //! @brief write the correlation plot data to a filestream that is readable by gnuplot
    //! @param OSTREAM output stream of interest
    //! @param EXP experimental data
    //! @param PRED predicted data
    void ModelComputeStatistics::WriteFormattedCorrelationInput
    (
      io::OFStream &OSTREAM,
      const linal::Vector< float> &EXP,
      const linal::Vector< float> &PRED
    ) const
    {
      OSTREAM << "EXPERIMENTAL\tPREDICTED\t" << PRED.GetSize() << "\t" << PRED.GetSize() << "\n";

      for( size_t element( 0); element < PRED.GetSize(); ++element)
      {
        OSTREAM << EXP( element) << "\t" << PRED( element) << "\n";
      }
    }

    //! @brief write ROC curve plot header for a gnuplot file to output stream
    //! @param OSTREAM output stream of interest
    //! @param OUTPUT_FILENAME output file name OSTREAM will write to
    //! @param ENRICHMENT enrichment value that is written on the plot
    //! @param AUC area under the curve value that is written on the plot
    void ModelComputeStatistics::WriteROCHeader
    (
      io::OFStream &OSTREAM,
      const std::string &OUTPUT_FILENAME,
      const float ENRICHMENT,
      const float AUC
    ) const
    {
      // write header to plot ROC curve with gnuplot
      OSTREAM << "set terminal " << GetGnuplotTerminalType() << " enhanced size 1200,600 # transparent size 2160,640\n";
      OSTREAM << "set output \"" << OUTPUT_FILENAME << "." << GetImageFormatExtension() << "\"\n";
      OSTREAM << "set encoding iso\n";
      OSTREAM << "set key right bottom\n";
      OSTREAM << "set title \"roc curve\"\n";
      OSTREAM << "set grid\n";
      OSTREAM << "set nologscale\n";
      OSTREAM << "set xlabel \"false positive rate\"\n";
      OSTREAM << "set xrange [0:1]\n";
      OSTREAM << "set ylabel \"true positive rate\"\n";
      OSTREAM << "set yrange [0:1]\n";
      OSTREAM << "set nokey\n";
      OSTREAM << "set multiplot\n";
      OSTREAM << "set origin 0.0,0.0\n";
      OSTREAM << "set size 0.5,1\n";
      OSTREAM << "plot x with lines linetype 1 linewidth 1 linecolor rgb 'black'\n";
      OSTREAM << "set style line 1 linetype 1 linewidth 0 pointtype 1 linecolor rgb 'red'\n";
      OSTREAM << "set key\n";
      OSTREAM << "plot \\\n";
      // text for legend in gnuplot
      OSTREAM << "\"" << OUTPUT_FILENAME << "\" every :::0::0 using 1:2 with line linestyle 1 title \"";
      OSTREAM << "average enrichment " << ENRICHMENT << ' ' << " area under curve " << AUC << "\"\n";
      OSTREAM << "set logscale x 10\n";
      OSTREAM << "set xrange [0.0001:1]\n";
      OSTREAM << "set yrange [0.0001:1]\n";
      OSTREAM << "set xlabel \"\"\n";
      OSTREAM << "set ylabel \"\"\n";
      OSTREAM << "set grid xtics y2tics\n";
      OSTREAM << "set origin 0.23,0.15\n";
      OSTREAM << "set size 0.25,0.5\n";
      OSTREAM << "clear\n";
      OSTREAM << "set notitle\n";
      OSTREAM << "set nokey\n";
      OSTREAM << "set noytics\n";
      OSTREAM << "set y2tics\n";
      OSTREAM << "replot\n";
      OSTREAM << "plot x with lines linetype 1 linewidth 1 linecolor rgb 'black'\n";
      OSTREAM << "set origin 0.5,0.0\n";
      OSTREAM << "set size 0.5,1\n";

      OSTREAM << "set title \"Precision vs Fraction Positive Predicted Values\"\n";
      OSTREAM << "set xlabel \"Fraction Positive Predicted Values\"\n";
      OSTREAM << "set ylabel \"Precision\"\n";
      OSTREAM << "plot \\\n";

      // text for legend in gnuplot
      OSTREAM << "\"" << OUTPUT_FILENAME << "\" every :::1::1 using 1:2 with line linestyle 1 title \"\"\n";

      OSTREAM << "set style line 1 linetype 1 linewidth 0 pointtype 1 linecolor rgb 'black'\n";
      OSTREAM << "plot \\\n";
      // text for legend in gnuplot
      OSTREAM << "\"" << OUTPUT_FILENAME << "\" every :::2::2 using 1:2 with line linestyle 1 title \"\"\n";
      OSTREAM << "set nomultiplot\n";
    }

    //! @brief write ROC curve plot header for a gnuplot file with all contingency matrix measures to output stream
    //! @param OSTREAM output stream of interest
    //! @param OUTPUT_FILENAME output file name OSTREAM will write to
    //! @param MEASURES a list of all measures that were outuput, in the order they were output
    void ModelComputeStatistics::WriteROCCompleteHeader
    (
      io::OFStream &OSTREAM,
      const std::string &OUTPUT_FILENAME,
      const storage::Vector< std::string> &MEASURES
    ) const
    {
      // get all the x-measures out of the flag
      storage::Vector< std::string> x_axes( m_PlotXFlag->GetStringList());
      const size_t number_unscaled_x_axes( x_axes.GetSize());
      x_axes.Append( m_PlotLogXFlag->GetStringList());
      storage::Vector< std::string> y_axes( m_PlotYFlag->GetStringList());
      storage::Vector< size_t> is_y_axis( MEASURES.GetSize(), size_t( 0));
      for
      (
        storage::Vector< std::string>::const_iterator itr( y_axes.Begin()), itr_end( y_axes.End());
        itr != itr_end;
        ++itr
      )
      {
        is_y_axis( MEASURES.Find( *itr)) = size_t( 1);
      }
      const std::string ideal_ppv_hit_relative
      (
        math::ContingencyMatrixMeasures::GetMeasureInfo( math::ContingencyMatrixMeasures::e_IdealPPVRelativeToHitRate).First()
      );
      const std::string ideal_ppv_fpr_relative
      (
        math::ContingencyMatrixMeasures::GetMeasureInfo( math::ContingencyMatrixMeasures::e_IdealPPVRelativeToFPR).First()
      );
      size_t x_axes_id( 0);
      for
      (
        storage::Vector< std::string>::const_iterator itr( x_axes.Begin()), itr_end( x_axes.End());
        itr != itr_end;
        ++itr, ++x_axes_id
      )
      {
        // get the measure for the x-axis
        util::Implementation< util::FunctionInterfaceSerializable< math::ContingencyMatrix, double> > x_measure( *itr);

        bool is_fpr_sensitive( *itr != "Cutoff" && *itr != "HitRate");

        // find the measure
        const size_t x_measure_index( MEASURES.Find( x_measure.GetAlias()));

        BCL_Assert
        (
          x_measure_index < MEASURES.GetSize(),
          "Could not find measure: " + x_measure.GetAlias() + " in measures " + util::Format()( MEASURES)
        );

        const bool is_log_axis( x_axes_id >= number_unscaled_x_axes);

        //output the desired columns
        OSTREAM << "set terminal " << GetGnuplotTerminalType() << " enhanced size 600,600\n"
                << "set output \"" << OUTPUT_FILENAME << '.' << ( is_log_axis ? "log" : "") << *itr << "." << GetImageFormatExtension() << "\"\n"
                << "set origin 0.0,0.0\n";
        if( !is_log_axis)
        {
          OSTREAM << "set nologscale\n"
                  << "set grid\n"
                  << "set autoscale fix\n";
        }
        else
        {
          OSTREAM << "set xrange [0.0001:1]\n"
                  << "set logscale x\n"
                  << "set grid\n"
                  << "set autoscale fix\n";
        }
        OSTREAM << "set key left bottom\n"
                << "set title \"Contingency Matrix Measures\"\n"
                << "set xlabel \"" << x_measure.GetAlias() << "\"\n"
                << "set ylabel \"\"\n";
        OSTREAM << "plot ";
        size_t measure_index( 0);
        bool have_printed_first_measure( false);
        for
        (
          storage::Vector< std::string>::const_iterator itr( MEASURES.Begin()), itr_end( MEASURES.End());
          itr != itr_end;
          ++itr, ++measure_index
        )
        {
          if( !is_y_axis( measure_index))
          {
            continue;
          }
          // for FPR sensitive measures, the Ideal PPV is a different curve than for Hit Rate
          if( is_fpr_sensitive && *itr == ideal_ppv_hit_relative)
          {
            continue;
          }
          else if( !is_fpr_sensitive && *itr == ideal_ppv_fpr_relative)
          {
            continue;
          }
          if( have_printed_first_measure)
          {
            OSTREAM << ", \\\n  ";
          }
          std::string measure_name( MEASURES( measure_index));
          if( measure_name == ideal_ppv_hit_relative || measure_name == ideal_ppv_fpr_relative)
          {
            measure_name = "Ideal-PPV";
          }
          have_printed_first_measure = true;
          OSTREAM << "\"" << OUTPUT_FILENAME << "\" using " << ( x_measure_index + 1) << ':' << ( measure_index + 1)
                  << " with line title \"" << measure_name << "\" ls " << measure_index + 1;
        }
        OSTREAM << "\nreset\n";
      }
    }

    //! @brief evaluate objective function values for given obj function labels per commandline flag
    //! @param PRED predicted data
    //! @param EXP experimental data
    //! @param NAME base name for writing out
    //! @param IDS any ids associated with the predictions
    void ModelComputeStatistics::EvaluateObjectiveFunctions
    (
      const linal::Matrix< float> &PRED, const model::FeatureDataSet< float> &EXP, const std::string &NAME
    ) const
    {
      // construct FeatureDataSets for given experimental and predicted values
      const model::FeatureDataSet< float> pred_fds( PRED);

      // string to capture obj function values and names
      std::stringstream obj_values_output;

      // iterate over all objective function labels, instantiate obj function and evaluate
      const storage::Vector< std::string> obj_function_names( m_ObjectiveFunction->GetStringList());
      storage::Vector< std::string>::const_iterator itr_obj_name( obj_function_names.Begin());
      for
      (
        storage::Vector< util::Implementation< model::ObjectiveFunctionInterface> >::iterator
          itr( m_ObjectiveFunctions.Begin()), itr_end( m_ObjectiveFunctions.End());
        itr != itr_end;
        ++itr, ++itr_obj_name
      )
      {
        // write out the actual objective function value
        obj_values_output << ( *itr)->operator()( EXP, pred_fds)
                          // write name, improvement type, goal type, threshold, and name
                          << '\t' << *itr_obj_name
                          << '\t' << opti::GetImprovementTypeName( ( *itr)->GetImprovementType())
                          << '\t' << NAME << '\n';
      }

      // write obj function names and values to file otherwise to terminal
      if( m_EvaluateObjectiveFunctionFilename->GetFlag())
      {
        // remove any existing file, since the file may be appended to multiple times
        io::OFStream ofstream;
        io::File::MustOpenOFStream( ofstream, m_EvaluateObjectiveFunctionFilename->GetFirstParameter()->GetValue(), std::ios::app | std::ios::out);
        ofstream << obj_values_output.str() << '\n';
        io::File::CloseClearFStream( ofstream);
      }
      else
      {
        BCL_MessageStd( obj_values_output.str());
      }
    }

    //! @brief standard constructor
    ModelComputeStatistics::ModelComputeStatistics() :
      m_InputFilenamesFlag
      (
        new command::FlagDynamic
        (
          "input",
          "list of input filenames",
          command::Parameter
          (
            "filename",
            "file that contains experimental and predicted data in either csv format (prediction columns (usually 1), "
            "experimental columns) or BCL format linal::Matrix<float> (for backwards compatibility only)",
            command::ParameterCheckFileExistence()
          ),
          1
        )
      ),
      m_MaxFilesPerConsensusFlag
      (
        new command::FlagStatic
        (
          "consensus_size_limit",
          "maximum number of files to consider together for making a consensus model.",
          command::Parameter
          (
            "consensus_size_limit",
            "maximum number of files to consider together for making a consensus model.",
            command::ParameterCheckRanged< size_t>(),
            util::Format()( util::GetUndefined< size_t>())
          )
        )
      ),
      m_TableNameFlag
      (
        new command::FlagStatic
        (
          "table_name", "name of the table containing the quality measures",
          command::Parameter
          (
            "table_name", "name of the table containing the quality measures", "table.txt"
          )
        )
      ),
      m_SortByFlag
      (
        new command::FlagStatic
        (
          "sort_by", "the quality measure to sort the table by",
          command::Parameter
          (
            "sort_by", "the quality measure to sort the table by",
            command::ParameterCheckAllowed
            (
              storage::Vector< std::string>::Create
              (
                "RMSD", "Correlation", "CorrelationSpearman", "R^2", "AUC", "Enrichment"
              )
            ),
            "RMSD"
          )
        )
      ),
      m_PlotXFlag
      (
        new command::FlagDynamic
        (
          "plot_x",
          "choose which contingency matrix measure to plot on the x-axis",
          command::Parameter
          (
            "measure",
            "measure to plot on the x-axis",
            command::ParameterCheckSerializable
            (
              util::Implementation< util::FunctionInterfaceSerializable< math::ContingencyMatrix, double> >()
            ),
            "FPR"
          )
        )
      ),
      m_PlotLogXFlag
      (
        new command::FlagDynamic
        (
          "plot_log_x",
          "choose which contingency matrix measure to plot on the log x-axis",
          command::Parameter
          (
            "measure",
            "measure to plot on the log x-axis",
            command::ParameterCheckSerializable
            (
              util::Implementation< util::FunctionInterfaceSerializable< math::ContingencyMatrix, double> >()
            ),
            "FPR"
          )
        )
      ),
      m_PlotYFlag
      (
        new command::FlagDynamic
        (
          "plot_y",
          "choose which contingency matrix measure to plot on the y-axis; if none are given, uses all measures except "
          "the x-axis.  This flag has no effect unless -plot_x is given",
          command::Parameter
          (
            "measures",
            "measures to plot on the y-axis",
            command::ParameterCheckSerializable
            (
              util::Implementation< util::FunctionInterfaceSerializable< math::ContingencyMatrix, double> >()
            ),
            ""
          )
        )
      ),
      m_NoPlotFlag
      (
        new command::FlagStatic
        (
          "no_plot",
          "Disable plotting or writing out data files for each output; instead, just calculate objective functions"
        )
      ),
      m_TakeLog10Flag
      (
        new command::FlagStatic
        (
          "take_log10", "takes log base 10 of the input values",
          command::Parameter
          (
            "take_log10", "takes log base 10 of the input values", ""
          )
        )
      ),
      m_CorrelationFlag
      (
        new command::FlagStatic
        (
          "correlation", "indicates that correlation plots and not ROC curves are desired"
        )
      ),
      m_PotencyCutOffFlag
      (
        new command::FlagStatic
        (
          "potency_cutoff", "potency cutoff for plotting roc curves",
          command::Parameter
          (
            "potency_cutoff", "potency cutoff for plotting roc curves",
            command::ParameterCheckRanged< float>( -std::numeric_limits< float>::max(), std::numeric_limits< float>::max()), ""
          )
        )
      ),
      m_ActivesBelowCutoffFlag
      (
        new command::FlagStatic
        (
          "actives_below_cutoff", "indicates that those values predicted below cutoff are treated as actives"
        )
      ),
      m_OutputDirectory
      (
        new command::FlagStatic
        (
          "output_directory", "directory where the output will be re-directed; if not specified the current directory is chosen",
          command::Parameter
          (
            "output_directory", "directory where the output will be re-directed; if not specified the current directory is chosen",
            "./"
          )
        )
      ),
      m_ObjectiveFunction
      (
        new command::FlagDynamic
        (
          "obj_function",
          "list of objective function labels.  The first objective function given is used to set -plot_x, -plot_y, -correlation, "
          "-actives_below_cutoff, and -potency_cutoff, if they are otherwise not given over the command line",
          command::Parameter
          (
            "obj_function",
            "list of objective function labels",
            command::ParameterCheckSerializable( util::Implementation< model::ObjectiveFunctionInterface>())
          )
        )
      ),
      m_EvaluateObjectiveFunctionFilename
      (
        new command::FlagStatic
        (
          "filename_obj_function", "file where all specified objective function evaluation are stored",
          command::Parameter
          (
            "filename_obj_function", "file where all specified objective function evaluation are stored",
            "eval_obj_results.txt"
          )
        )
      ),
      m_ImageFormatFlag
      (
        new command::FlagStatic
        (
          "image_format",
          "specify the image type to display, if outputting graphs of quality measures",
          command::Parameter
          (
            "format",
            "the format to use",
            command::ParameterCheckAllowed( storage::Vector< std::string>::Create( "png", "svg")),
            "png"
          )
        )
      )
    {
    }

    const ApplicationType ModelComputeStatistics::ModelComputeStatistics_Instance
    (
      GetAppGroups().AddAppToGroup( new ModelComputeStatistics(), GetAppGroups().e_Model)
    );

    //! @brief returns readme information
    //! @return string containing information about application
    const std::string &ModelComputeStatistics::GetReadMe() const
    {
      static const std::string readme
      (
        DefaultSectionSeparator() +
        "I. OVERVIEW.\n"
        "This document provides a description of BCL model:ComputeStatistics, terms of use, appropriate citation, installation "
        "procedures, BCL model:ComputeStatistics execution, technical support, and future research directions.\n"
        "\n"
        + DefaultSectionSeparator() +
        "II. WHAT IS BCL model:ComputeStatistics?\n"
        "BCL model:ComputeStatistics is a C++ based application, created by Vanderbilt University's Meiler Laboratory, which is part "
        "of a larger library of applications called BCL::Commons. BCL model:ComputeStatistics reads in the prediction output "
        "of a trained machine learning model and evaluates an array of quality measures and objective functions to "
        "determine the prediction performance of the trained model. BCL model:ComputeStatistics outputs a table file with all"
        "calculated information and provides the user with specific ROC curve type and correlation plots for visualization "
        "purposes. In addition, BCL model:ComputeStatistics is able to read in multiple prediction result files to perform"
        " a consensus prediction for all prediction combinations."
        "\n"
        + DefaultSectionSeparator() +
        "III. TERMS OF USE.\n"
        + DefaultTermsOfUseString() +
        "\n"
        + DefaultSectionSeparator() +
        "IV. APPROPRIATE CITATIONS FOR USING BCL model:ComputeStatistics.\n"
        "When using BCL model:ComputeStatistics in a publication, please cite the following publication describing the application's "
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
        "VI. RUNNING BCL model:ComputeStatistics.\n"
        "Running BCL model:ComputeStatistics consists of the following steps.\n"
        "\n"
        "1) At a command prompt, navigate to the location of your BCL model:ComputeStatistics executable program.\n"
        "\n"
        "2) Run BCL model:ComputeStatistics on one or multiple <experimental predicted value file>s\n"
        "\n"
        "3) Obtain a gnuplot plot of your chosen quality measures or a correlation plot."
        "\n"
        "ADDITIONAL INFORMATION\n"
        "\n"
        "Information concerning syntax and flags  can be obtained by typing: <bcl.exe>  model:ComputeStatistics -help\n"
        "\n"
        "For further explanation, examples of the flags, example command lines, input and output format information,\n"
        "and example output please see the documentation file.\n"
        "\n"
        + DefaultSectionSeparator() +
        "VII. TECHNICAL SUPPORT.\n"
        + DefaultTechnicalSupportString() +
        "\n"
        + DefaultSectionSeparator() +
        "VIII. FUTURE DEVELOPMENT OF BCL model:ComputeStatistics.\n"
        "BCL model:ComputeStatistics is under ongoing further development.\n"
        "\n"
        + DefaultSectionSeparator() +
        "IX. EXAMPLE COMMANDLINE"
        "<bcl.exe> model:ComputeStatistics -input <experimental predicted values file> <experimental predicted values file> "
        "-potency_cutoff <your experimental actives/inactives cutoff value> "
        "-table_name results_table.txt -plot_x FPR -plot_y InformationGainRatio TPR PPV Ideal-PPV\n"
        "\n"
        "Get an overview about all flags and options when typing: <bcl.exe>  model:ComputeStatistics -help\n"
        "\n"
        "Alternatively, most of the flags can be omitted by providing the -objective_function flag, which causes "
        "only that objective to be evaluated, and potency cutoffs and plots will be made based on the type of objective "
        "function provided and its parameters\n"
        "The <experimental predicted values file> should be in csv format:\n"
        "predicted_value1,prediction_value2,...,experimental_value1,experimental_value2,...\n"
        "predicted_value1,prediction_value2,...,experimental_value1,experimental_value2,...\n"
        "...\n"
        "Where each row provides the predicted and experimental values for training point. "
        "optionally, the first column may be an ID column, which will be ignored except for certain "
        "ID-dependent objective functions, such as SegmentOverlap\n"
        "For backwards compatibility, a BCL Matrix format is also supported, which looks like:\n"
        "\n"
        "bcl::linal::Matrix<float>\n"
        "number_rows number_columns\n"
        "experimental_value1 predicted_value1 experimental_value2 predicted_value2 ...\n"
        "experimental_value1 predicted_value1 experimental_value2 predicted_value2 ...\n"
        "...\n"
        "\n"
        "An <experimental predicted values file> can be obtained through the BCL application model:Train. "
        "and model::PredictionMerge\n"
        "\n"
        "When providing multiple <experimental predicted values file>s with your -input flag, an averaged consenus prediction"
        " for all input combinations is calculated, and consensus model pairs or higher order combinations of models "
        " can be computed\n"
        "\n"
        "The resulting *.gnuplot file are processed by the following command to obtain an png file:\n"
        "gnuplot *.gnuplot\n"
        "if gnuplot is available and in your path.  If not, may run gnuplot manually to obtain the pngs\n"
        + DefaultSectionSeparator()
      );

      return readme;
    }

    //! @brief get a description for the app
    //! @return a brief (no more than 3 line) description for the application
    std::string ModelComputeStatistics::GetDescription() const
    {
      const std::string description
      (
        "Evaluate quality measures of qsar model preditions or just experimental/predicted values and present results"
        "in table and as gnuplot graphics"
      );

      return description;
    }

    // This function is used to retrieve prior names for applications that are still valid (for backwards compatibility)
    // but which should not be displayed, e.g. for help
    storage::Vector< std::string> ModelComputeStatistics::GetDeprecatedAppNames() const
    {
      return storage::Vector< std::string>( size_t( 1), "ComputeJuryStatistics");
    }

  } // namespace app
} // namespace bcl
