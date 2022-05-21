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

// includes from bcl - sorted alphabetically
#include "app/bcl_app_apps.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_binary_function_bind_second.h"
#include "math/bcl_math_contingency_matrix.h"
#include "math/bcl_math_contingency_matrix_measures.h"
#include "math/bcl_math_function_adapter.h"
#include "math/bcl_math_mutate_vector.h"
#include "math/bcl_math_roc_curve.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_statistics.h"
#include "math/bcl_math_template_instantiations.h"
#include "mc/bcl_mc_approximator.h"
#include "mc/bcl_mc_temperature_accepted.h"
#include "opti/bcl_opti_criterion_combine.h"
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "opti/bcl_opti_criterion_unimproved.h"
#include "sched/bcl_sched_binary_function_job_with_data.h"
#include "sched/bcl_sched_sum_function.h"
#include "score/bcl_score_consensus_enrichment.h"
#include "storage/bcl_storage_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FunctionCombine
    //! @brief combines multiple functions and returns their result in a vector
    //!
    //! @author woetzen
    //! @date 02/06/2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ArgumentType, typename t_ResultType>
    class FunctionCombine :
      public math::FunctionInterfaceSerializable< t_ArgumentType, storage::Vector< t_ResultType> >
    {
    private:

    //////////
    // data //
    //////////

      //! store all functions in a sharedpointerlist
      util::ShPtrList< math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> > m_Functions;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FunctionCombine()
      {
      }

      //! @brief Clone function
      //! @return pointer to new FunctionCombine< t_ArgumentType, t_ResultType>
      FunctionCombine< t_ArgumentType, t_ResultType> *Clone() const
      {
        return new FunctionCombine< t_ArgumentType, t_ResultType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get functions list
      //! @return functions list
      const util::ShPtrList< math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> > &GetFunctions() const
      {
        return m_Functions;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief insert additional function object
      void PushBack( const util::ShPtr< math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> > &SP_FUNCTION)
      {
        m_Functions.PushBack( SP_FUNCTION);
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator, that takes the argument and returns result
      //! @param ARGUMENT the argument the function is calculated for
      //! @return list of results, the function value for each function in the FunctionCombine object
      storage::Vector< t_ResultType> operator()( const t_ArgumentType &ARGUMENT) const
      {
        storage::Vector< t_ResultType> results;
        results.AllocateMemory( m_Functions.GetSize());

        // iterate over all functions and insert the result for each function into the results
        for
        (
          typename util::ShPtrList< math::FunctionInterfaceSerializable< t_ArgumentType, t_ResultType> >::const_iterator
            sp_func_itr( m_Functions.Begin()), sp_func_itr_end( m_Functions.End());
          sp_func_itr != sp_func_itr_end;
          ++sp_func_itr
        )
        {
          results.PushBack( ( *sp_func_itr)->operator()( ARGUMENT));
        }

        // end
        return results;
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        // read member
        io::Serialize::Read( m_Functions, ISTREAM);

        // end
        return ISTREAM;
      }

      //! write to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_Functions, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class FunctionCombine

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class NthBest
    //! @brief sorts a given vector using std::less and returns Nth smallest value.
    //!
    //! @author woetzen
    //! @date Jun 08, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ResultType>
    class NthBest :
      public math::FunctionInterfaceSerializable< storage::Vector< t_ResultType>, t_ResultType>
    {
    private:

    //////////
    // data //
    //////////

      //! the position to be used after sorting the results ascending
      size_t m_PositionN;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      NthBest() :
        m_PositionN( 0)
      {
      }

      //! @brief construct from position
      NthBest( const size_t POSITION_N) :
        m_PositionN( POSITION_N)
      {
      }

      //! @brief virtual copy constructor
      //! @return pointer to new copy of this
      NthBest< t_ResultType> *Clone() const
      {
        return new NthBest< t_ResultType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator, that takes the list of arguments and returns the result at PositionN
      //! @param ARGUMENT the argument the function is calculated for
      //! @return result at PositionN after sorting
      t_ResultType operator()( const storage::Vector< t_ResultType> &ARGUMENT) const
      {
        // check that the argument
        BCL_Assert
        (
          m_PositionN < ARGUMENT.GetSize(),
          "The requested position " + util::Format()( m_PositionN) + " is larger than size of argument " +
          util::Format()( ARGUMENT.GetSize())
        );

        // make a copy
        storage::Vector< t_ResultType> results( ARGUMENT.Begin(), ARGUMENT.End());

        // sort
        results.Sort( std::less< t_ResultType>());

        // end
        return results( m_PositionN);
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        // read member
        io::Serialize::Read( m_PositionN, ISTREAM);

        // end
        return ISTREAM;
      }

      //! write to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // write member
        io::Serialize::Write( m_PositionN, OSTREAM, INDENT);

        // end
        return OSTREAM;
      }

    }; // template class NthBest

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SumVector
    //! @brief sums all values in a vector
    //!
    //! @author woetzen
    //! @date Jun 08, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_ResultType>
    class SumVector :
      public math::FunctionInterfaceSerializable< storage::Vector< t_ResultType>, t_ResultType>
    {
    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SumVector()
      {
      }

      //! @brief virtual copy constructor
      //! @return pointer to new copy of this
      SumVector< t_ResultType> *Clone() const
      {
        return new SumVector< t_ResultType>( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator, that takes the list of arguments and returns the result at PositionN
      //! @param ARGUMENT the argument the function is calculated for
      //! @return result at PositionN after sorting
      t_ResultType operator()( const storage::Vector< t_ResultType> &ARGUMENT) const
      {
        // copy
        storage::Vector< t_ResultType> copy( ARGUMENT);
        std::transform( copy.Begin(), copy.End(), copy.Begin(), std::ptr_fun< double, double>( std::sqrt));

        return std::accumulate( copy.Begin(), copy.End(), t_ResultType( 0.0));
      }

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! read from std::istream
      std::istream &Read( std::istream &ISTREAM)
      {
        // end
        return ISTREAM;
      }

      //! write to std::ostream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        // end
        return OSTREAM;
      }

    }; // template class NthBest

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MinimizeScoreWeightSet
    //! @brief This is application is used for minimizing a given weightset using various criteria. It reads in a score
    //! table, picks the users specified columns to optimize weights for as well as a column which specifies the
    //! objective value. Depending on options given, it can optimize the same weightset using multiple different score
    //! tables as well as creating multiple cross-validation tables for each score table.
    //!
    //! @author woetzen, karakam
    //! @date Jul 29, 2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MinimizeScoreWeightSet :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! input score table file list
      util::ShPtr< command::FlagInterface> m_ScoreTableFileList;

      //! input a single table name
      util::ShPtr< command::FlagInterface> m_ScoreTableFilename;

      //!max number of rejected steps and max iterations
      util::ShPtr< command::FlagStatic>         m_McMaxIterationsUnimprovedFlag;
      util::ShPtr< command::ParameterInterface> m_McMaxIterationsParam;
      util::ShPtr< command::ParameterInterface> m_MCMaxStepsUnimprovedParam;

      //! flag for cross validation
      util::ShPtr< command::FlagStatic>         m_CrossValidationFlag;
      util::ShPtr< command::ParameterInterface> m_CriteriaNameParam;
      util::ShPtr< command::ParameterInterface> m_CriteriaCutoffParam;
      util::ShPtr< command::ParameterInterface> m_FractionParam;
      util::ShPtr< command::ParameterInterface> m_NumberCrossValidationParam;
      util::ShPtr< command::ParameterInterface> m_MinimalTableSizeParam;

      //! flag to pass file containing math::Vector with weight set to combine score
      util::ShPtr< command::FlagInterface> m_WeightSetStartFlag;

      //! flag to identify where the final optimized weightset should be written to
      util::ShPtr< command::FlagInterface> m_WeightSetWriteFlag;

      //! flag to specify the given ranges as percentages instead of absolute ranges
      util::ShPtr< command::FlagStatic> m_UsePercentageRangesFlag;

      //! flag that writes the roccurve for the given tables using the weightset
      util::ShPtr< command::FlagInterface> m_WriteRocFlag;

      //! flag to indicate if smaller values or larger values for the criteria are true
      util::ShPtr< command::FlagInterface> m_SortOrderFlag;

      //! number of repetitions for the minimization
      util::ShPtr< command::FlagInterface> m_NumberRepeatsFlag;

      //! Nth worst table to improve the enrichment upon
      util::ShPtr< command::FlagInterface> m_TargetTableIndexFlag;

      //! number of weights mutated in a step
      util::ShPtr< command::FlagInterface> m_NumberWeightsMutatedInOneStepFlag;

      //! keep the weights positive for each iteration
      util::ShPtr< command::FlagInterface> m_KeepWeightsPositiveFlag;

      //! score to keep constant during each iteration. If left blank, no scores will be kept constant
      //! for sure
      util::ShPtr< command::FlagInterface> m_ConstantFlag;

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
      MinimizeScoreWeightSet();

      //! @brief Clone function
      //! @return pointer to new MinimizeScoreWeightSet
      MinimizeScoreWeightSet *Clone() const
      {
        return new MinimizeScoreWeightSet( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      int Main() const;

      storage::List< storage::Table< double> > CreateCrossValidationTablesFromFile
      (
        const std::string &FILE_NAME,
        const storage::Table< double> &WEIGHT_RANGE_TABLE
      ) const;

      storage::List< storage::Table< double> > CreateCrossValidationTables
      (
        const storage::Table< double> &TRUE_ENTRIES,
        const storage::Table< double> &FALSE_ENTRIES
      ) const;

      //! @brief create objective function that combines the enrichments for multiple tables by balancing and summing
      //! @param TABLE_FILE_NAMES list of filenames to each of the tables to be considered
      //! @param CRITERIA_NAME    name of the column that contains the criteria for the enrichment
      //! @param WEIGHT_RANGE_TABLE table of starting weights and ranges - passed in to ensure that all tables have the
      //!        same header and we have weights for each of their columns
      //! @param CONSENSUS_MAP map that will store the ConsensusEnrichment for each cross validation for each table file
      //! @return sum enrichment function for balanced tables
      util::ShPtr< util::FunctionInterface< linal::Vector< double>, double> > EnrichmentObjective
      (
        const storage::Vector< std::string> &TABLE_FILE_NAMES,
        const std::string &CRITERIA_NAME,
        const storage::Table< double> &WEIGHT_RANGE_TABLE,
        storage::Map< std::string, util::ShPtr< FunctionCombine< linal::Vector< double>, double> > > &CONSENSUS_MAP,
        math::RunningAverageSD< storage::Vector< double> > &AVE_SD
      ) const;

      //! @brief create objective function that calculates the ROCIntergral for one table and a given criteria
      //! @param TABLE_FILE_NAME file name of the table to be used
      //! @param CRITERIA_NAME    name of the column that contains the criteria for the enrichment
      //! @param WEIGHT_RANGE_TABLE table of starting weights and ranges - passed in to ensure that all tables have the
      //!        same header and we have weights for each of their columns
      //! @param COMPARISON binary predicate that classifies true or false if above criteria threshold
      //! @return ROCIntergral function on the given table
      util::ShPtr< util::FunctionInterface< linal::Vector< double>, double> > ROCIntegralObjective
      (
        const std::string &TABLE_FILE_NAME,
        const std::string &CRITERIA_NAME,
        const storage::Table< double> &WEIGHT_RANGE_TABLE,
        const math::Comparisons< double>::Comparison &COMPARISON,
        math::RunningAverageSD< storage::Vector< double> > &AVE_SD
      ) const;

      //! @brief prints out the individual enrichments in a table format
      //! @param WEIGHTS vector of weights
      //! @param CONSENSUS_MAP map that will store the ConsensusEnrichment for each cross validation for each table file
      void PrintIndividualEnrichments
      (
        const linal::Vector< double> &WEIGHTS,
        const storage::Map< std::string, util::ShPtr< FunctionCombine< linal::Vector< double>, double> > > &CONSENSUS_MAP
      ) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    private:

      static const ApplicationType MinimizeScoreWeightSet_Instance;

    }; // MinimizeScoreWeightSet

    util::ShPtr< command::Command> MinimizeScoreWeightSet::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // input Score table file list
      sp_cmd->AddFlag( m_ScoreTableFileList);

      // input a single table
      sp_cmd->AddFlag( m_ScoreTableFilename);

      // flags for monte carlo minimization
      // adjust start and end condition
      sp_cmd->AddFlag( mc::TemperatureAccepted::GetFlagTemperature());

      // max number of rejected steps and max iterations
      sp_cmd->AddFlag( m_McMaxIterationsUnimprovedFlag);

      // cross validation
      sp_cmd->AddFlag( m_CrossValidationFlag);

      // flag for weight set defined in given file
      sp_cmd->AddFlag( m_WeightSetStartFlag);

      // flag to identify where the final optimized weightset should be written to
      sp_cmd->AddFlag( m_WeightSetWriteFlag);

      // flag to specify the given ranges as percentages instead of absolute ranges
      sp_cmd->AddFlag( m_UsePercentageRangesFlag);

      // flag for writing the actual roccurve
      sp_cmd->AddFlag( m_WriteRocFlag);

      // flag to determine the order
      sp_cmd->AddFlag( m_SortOrderFlag);

      // flag to determine number of repetitions of the minimization
      sp_cmd->AddFlag( m_NumberRepeatsFlag);

      // Nth worst table to improve the enrichment upon
      sp_cmd->AddFlag( m_TargetTableIndexFlag);

      // number of weights mutated in a step
      sp_cmd->AddFlag( m_NumberWeightsMutatedInOneStepFlag);

      // keep weigths positive
      sp_cmd->AddFlag( m_KeepWeightsPositiveFlag);

      // constant
      sp_cmd->AddFlag( m_ConstantFlag);

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief main function
    //! @return 0 on success of error code
    int MinimizeScoreWeightSet::Main() const
    {
      // check that the provided table index is valid, smaller than number cross validations
      BCL_Assert
      (
        m_TargetTableIndexFlag->GetFirstParameter()->GetNumericalValue< size_t>() <
        m_NumberCrossValidationParam->GetNumericalValue< size_t>(),
        "The target table index has to be smaller than number of cross validations"
      )

      // acquire comparison function
      const math::Comparisons< double>::Comparison comp_function( m_SortOrderFlag->GetFirstParameter()->GetValue());

      // name of criteria
      const std::string &criteria_name( m_CriteriaNameParam->GetValue());

      // setup ranges - parameters to be mutated
      storage::TableHeader temp_header;
      storage::Table< double> weightset_table_read( temp_header);

      // if a specified weightset file was given
      if( m_WeightSetStartFlag->GetFlag())
      {
        // read it in
        io::IFStream weight_file_stream;
        io::File::MustOpenIFStream
        (
          weight_file_stream,
          m_WeightSetStartFlag->GetFirstParameter()->GetValue()
        );

        // read starting weights and mutate ranges
        weightset_table_read.ReadFormatted( weight_file_stream);

        // make sure the table has a weights and a ranges rows
        BCL_Assert
        (
          weightset_table_read.GetSize() == 2 &&
          weightset_table_read.HasRow( "weights") &&
          weightset_table_read.HasRow( "ranges"),
          "The weightset table given needs to have rows named \"weights\" and \"ranges\""
        );

        // close stream
        io::File::CloseClearFStream( weight_file_stream);
      }
      // otherwise write a template
      else
      {
        BCL_MessageCrt
        (
          "no weight file was found - write default to file: " + m_WeightSetStartFlag->GetFirstParameter()->GetValue()
          + " " + "and add column names according to the score tables"
        );
        weightset_table_read.InsertRow( "weights");
        weightset_table_read.InsertRow( "ranges");
        io::OFStream write;
        io::File::MustOpenOFStream( write, m_WeightSetStartFlag->GetFirstParameter()->GetValue());
        weightset_table_read.WriteFormatted( write);
        io::File::CloseClearFStream( write);
        return 0;
      }

      // get the weights vector and range vector
      storage::Vector< double> weights_vector_read( weightset_table_read[ "weights"].GetData());
      storage::Vector< double> ranges_vector_read( weightset_table_read[ "ranges"].GetData());

      // initialize a vector of valid column names
      storage::Vector< std::string> valid_column_names;

      // iterate over column indices
      for( size_t index( 0), index_max( weightset_table_read.GetNumberColumns()); index < index_max; ++index)
      {
        // if the weights and the ranges are not 0 for this column
        if( weights_vector_read( index) != 0.0 || ranges_vector_read( index) != 0.0)
        {
          valid_column_names.PushBack( weightset_table_read.GetHeader()( index));
        }
      }

      // create a new weight set table that has only columns that have non-zero weights and/or non-zero ranges
      storage::Table< double> weightset_table( weightset_table_read.ExtractSubTableFromColumns( valid_column_names));

      // make sure the weight set table does not have the criteria
      BCL_Assert
      (
        !weightset_table.GetHeader().HasColumn( criteria_name),
        "The weightset can't have an assigned weight or range to the criteria to be used: " + criteria_name
      );

      // insert the criteria into weight set with 0 weight and 0 range
      storage::Table< double> transposed_table( weightset_table.GetTransposedTable());
      transposed_table.InsertRow( criteria_name, storage::Vector< double>::Create( 0.0, 0.0));
      weightset_table = transposed_table.GetTransposedTable();

      // create sum function object
      util::ShPtr< math::FunctionInterfaceSerializable< linal::Vector< double>, double> >
        sp_objective_function;

      // create consensus enrichment list
      storage::Map< std::string, util::ShPtr< FunctionCombine< linal::Vector< double>, double> > > consensus_map;

      math::RunningAverageSD< storage::Vector< double> > average_sd_scores;

      // if a score table was given
      if( m_ScoreTableFileList->GetFirstParameter()->GetWasSetInCommandLine())
      {
        // read all score tables and create one sum objective function
        io::IFStream score_table_file_list;
        io::File::MustOpenIFStream( score_table_file_list, m_ScoreTableFileList->GetFirstParameter()->GetValue());

        // list of filenames for scoring tables considered
        storage::Vector< std::string> table_file_names( util::StringListFromIStream( score_table_file_list));

        // close stream
        io::File::CloseClearFStream( score_table_file_list);

        // append single table to file name list
        if( m_ScoreTableFilename->GetFirstParameter()->GetWasSetInCommandLine())
        {
          BCL_MessageStd( "append single table passed on command line to the list of tables");
          table_file_names.PushBack( m_ScoreTableFilename->GetFirstParameter()->GetValue());
        }

        // initialize objective function
        sp_objective_function =
          EnrichmentObjective( table_file_names, criteria_name, weightset_table, consensus_map, average_sd_scores);
        // print the initial scores
        BCL_MessageStd( "The enrichment table before the minimization");

        PrintIndividualEnrichments( weightset_table[ "weights"].GetData(), consensus_map);
      }
      // mode for single table
      else if( m_ScoreTableFilename->GetFirstParameter()->GetWasSetInCommandLine())
      {
        sp_objective_function = ROCIntegralObjective
        (
          m_ScoreTableFilename->GetFirstParameter()->GetValue(),
          criteria_name,
          weightset_table,
          comp_function,
          average_sd_scores
        );
      }
      else
      {
        BCL_MessageCrt( "neither a table list nor a single table was given!");
        return 1;
      }
      if( 1)
      {
        linal::VectorConstReference< double> ref
        (
          average_sd_scores.GetStandardDeviation().GetSize(),
          &*average_sd_scores.GetStandardDeviation().Begin()
        );
        auto itr_ref( ref.Begin());
        for
        (
          auto itr( weightset_table[ "weights"].Begin()), itr_end( weightset_table[ "weights"].End());
          itr != itr_end;
          ++itr, ++itr_ref
        )
        {
          *itr = 1.0 / ( 0.000001 + *itr_ref);
        }
        weightset_table[ "weights"]( ref.GetSize() - 1) = 0.0;
      }

      // do not optimize if roc curves have been requested
      if( m_WriteRocFlag->GetFlag())
      {
        return 0;
      }

      // print out the score before the minimizations
      const double starting_enrich_roc( sp_objective_function->operator ()( weightset_table[ "weights"].GetData()));
      BCL_MessageStd
      (
        "Pre-minimization enrichment or ROC-integral: " + util::Format()( starting_enrich_roc)
      );

      // initialize mutate object
      util::ShPtr< math::MutateInterface< linal::Vector< double> > > sp_mutate
      (
        new math::MutateVector
        (
          linal::Vector< double>( weightset_table[ "ranges"].GetData()),
          m_NumberWeightsMutatedInOneStepFlag->GetFirstParameter()->GetNumericalValue< size_t>(),
          m_UsePercentageRangesFlag->GetFlag(),
          m_KeepWeightsPositiveFlag->GetFlag(),
          m_ConstantFlag->GetFlag()
          ? weightset_table.GetHeader()[ m_ConstantFlag->GetFirstParameter()->GetValue()]
          : util::GetUndefinedSize_t()
        )
      );

      // create a table to hold all weightsets with their enrichment sum
      storage::Vector< std::string> all_table_header( weightset_table.GetHeader());
      all_table_header.PushBack( "enrich_sum");

      storage::Table< double> all_table( all_table_header);
      {
        storage::Vector< double> starting_weights( weightset_table[ "weights"].GetData());
        starting_weights.PushBack( starting_enrich_roc);
        all_table.InsertRow( "weights_start", starting_weights);
      }

      // for the requested number of repeats
      for
      (
        size_t num_repeats( 0), max_number_repeats( m_NumberRepeatsFlag->GetFirstParameter()->GetNumericalValue< size_t>());
        num_repeats < max_number_repeats;
        ++num_repeats
      )
      {
        // create the temperature control
        util::ShPtr< mc::TemperatureInterface> sp_temperature
        (
          new mc::TemperatureAccepted
          (
            mc::TemperatureAccepted::GetParameterStartFraction()->GetNumericalValue< double>(),
            mc::TemperatureAccepted::GetParameterEndFraction()->GetNumericalValue< double>(),
            m_McMaxIterationsParam->GetNumericalValue< size_t>(),
            mc::TemperatureAccepted::GetParameterStartTemperature()->GetNumericalValue< double>(),
            mc::TemperatureAccepted::GetParameterUpdateInterval()->GetNumericalValue< size_t>()
          )
        );

        // create the metropolis
        util::ShPtr< mc::Metropolis< double> > sp_metropolis( new mc::Metropolis< double>( sp_temperature, true));

        // create the termination criterion
        opti::CriterionCombine< linal::Vector< double>, double> criterion_combine;

        // add termination criterion that depends on the total number of MC iterations
        criterion_combine.InsertCriteria
        (
          opti::CriterionNumberIterations< linal::Vector< double>, double>
          (
            m_McMaxIterationsParam->GetNumericalValue< size_t>()
          )
        );

        // add termination criterion that depends on the number of unimproved steps in a row
        criterion_combine.InsertCriteria
        (
          opti::CriterionUnimproved< linal::Vector< double>, double>
          (
            m_MCMaxStepsUnimprovedParam->GetNumericalValue< size_t>()
          )
        );

        // create the approximator
        mc::Approximator< linal::Vector< double>, double> approximator
        (
          *sp_objective_function,
          *sp_mutate,
          *sp_metropolis,
          criterion_combine,
          linal::Vector< double>( weightset_table[ "weights"].GetData())
        );

        // approximate
        BCL_MessageStd( "starting optimization");
        approximator.Approximate();

        // get the tracker
        const opti::Tracker< linal::Vector< double>, double> &tracker( approximator.GetTracker());

        // get the approximation result
        const util::ShPtr< storage::Pair< linal::Vector< double>, double> > sp_approximation_result( tracker.GetBest());

        // print total number of MC iterations
        BCL_MessageStd( "Total number of mc iterations: " + util::Format()( tracker.GetIteration()));

        // retrieve very best entry with its score
        BCL_MessageStd
        (
          "optimization number " + util::Format()( num_repeats) +
          " completed.\nThis is the enrichment or ROC integral achieved: " +
          util::Format()( sp_approximation_result->Second())
        );

        // final weight set
        storage::Table< double> final_weightset( weightset_table.GetHeader());
        final_weightset.InsertRow
        (
          "weights",
          storage::Vector< double>
          (
            sp_approximation_result->First().Begin(), sp_approximation_result->First().End()
          )
        );

        // if multiple tables were given
        if( m_ScoreTableFileList->GetFirstParameter()->GetWasSetInCommandLine())
        {
          // print out the enrichments
          BCL_MessageStd( "The enrichment table after the minimization");
          PrintIndividualEnrichments( sp_approximation_result->First(), consensus_map);
        }

        BCL_MessageStd( "this is the optimized weight set: ");
        storage::Vector< double> normalized_data( final_weightset[ "weights"].GetData());
        math::Statistics::SetToSum( normalized_data.Begin(), normalized_data.End(), double( 1.0));
        // insert a row with normalized weights
        final_weightset.InsertRow( "normalized_weights", normalized_data);
        storage::Vector< double> influence_weights( final_weightset[ "weights"].GetData());
        auto itr_id( influence_weights.Begin());
        for
        (
          auto itr( average_sd_scores.GetStandardDeviation().Begin()),
               itr_end( average_sd_scores.GetStandardDeviation().End());
          itr != itr_end;
          ++itr, ++itr_id
        )
        {
          *itr_id *= *itr;
        }
        math::Statistics::SetToSum( influence_weights.Begin(), influence_weights.End(), double( 1.0));
        final_weightset.InsertRow( "influence_weights", influence_weights);
        final_weightset.WriteFormatted( util::GetLogger());

        // if flag was given to write to file
        if( m_WeightSetWriteFlag->GetFlag())
        {
          // open output
          io::OFStream write_weightset;
          io::File::MustOpenOFStream
          (
            write_weightset,
            m_WeightSetWriteFlag->GetFirstParameter()->GetValue() + util::Format()( num_repeats) + ".weights"
          );
          final_weightset.WriteFormatted( write_weightset);
          io::File::CloseClearFStream( write_weightset);
        }

        storage::Vector< double> weights_and_enrich( final_weightset[ "weights"].GetData());
        weights_and_enrich.PushBack( sp_approximation_result->Second());
        all_table.InsertRow( "weights_" + util::Format()( num_repeats), weights_and_enrich);
      } // number repeats

      // print all weights and enrichment table
      BCL_MessageStd( "table of all repetitions, with their weights and enrichments:");
      all_table.WriteFormatted( util::GetLogger());

      // end
      return 0;
    }

    storage::List< storage::Table< double> > MinimizeScoreWeightSet::CreateCrossValidationTablesFromFile
    (
      const std::string &FILE_NAME,
      const storage::Table< double> &WEIGHT_RANGE_TABLE
    ) const
    {
      // open the file stream
      io::IFStream score_table_stream;
      io::File::MustOpenIFStream( score_table_stream, FILE_NAME);

      BCL_MessageStd( "reading score table: " + FILE_NAME);
      // read table from stream including rows with similar names
      storage::Table< double> current_table;
      current_table.ReadFormatted( score_table_stream, true);
      // close stream
      io::File::CloseClearFStream( score_table_stream);
      BCL_MessageStd( "finished reading score table: " + FILE_NAME);

      // check if they have the same headers
      const bool same_headers( current_table.GetHeader() == WEIGHT_RANGE_TABLE.GetHeader());

      // make sure it is same or at least has all the necessary columns
      BCL_Assert
      (
        same_headers || current_table.GetHeader().HasColumns( WEIGHT_RANGE_TABLE.GetHeader()),
        "The given table has to either have all the same headers or contain all the ones used in weight optimization"
        + util::Format()( current_table.GetHeader()) + "\nvs\n" + util::Format()( WEIGHT_RANGE_TABLE.GetHeader())
      );

      // if they were not the same then extract the required sub-table
      if( !same_headers)
      {
        current_table = current_table.ExtractSubTableFromColumns( WEIGHT_RANGE_TABLE.GetHeader());
      }

      const std::string criteria_name( m_CriteriaNameParam->GetValue());
      // acquire comparison function
      const double criteria_cutoff( m_CriteriaCutoffParam->GetNumericalValue< double>());
      const math::Comparisons< double>::Comparison comp_func( m_SortOrderFlag->GetFirstParameter()->GetValue());
      const util::ShPtr< math::FunctionInterfaceSerializable< double, bool> > sp_unary_classifier
      (
        math::Comparisons< double>::GetEnums().CreateUnaryPredicate( comp_func, criteria_cutoff).Clone()
      );

      // sort by criteria
      current_table.SortByColumn( criteria_name, **comp_func);

      // find first iterator that is larger than given criteria threshold
      const size_t criteria_index( current_table.GetHeader()[ criteria_name]);

      // if roc flag was given
      if( m_WriteRocFlag->GetFlag())
      {
        // initialize vector of weights
        const linal::Vector< double> weights( WEIGHT_RANGE_TABLE[ "weights"].GetData());

        // initialize list to store the classified results
        storage::List< storage::Pair< double, bool> > classified_results;

        // iterate over this table
        for
        (
          storage::Table< double>::const_iterator itr( current_table.Begin()), itr_end( current_table.End());
          itr != itr_end;
          ++itr
        )
        {
          // get the scores
          const linal::Vector< double> scores( itr->Second().GetData());
          // store the weighted scores with the associated classification
          classified_results.PushBack
          (
            storage::Pair< double, bool>
            (
              scores * weights,
              sp_unary_classifier->operator()( itr->Second()( criteria_index))
            )
          );
        }

        // sor the results
        classified_results.Sort
        (
          storage::PairBinaryPredicateFirst< double, bool>
          (
            util::BinaryFunctionSTLWrapper< std::less< double> >()
          )
        );

        // create a roc curve from the sorted classified results
        const math::ROCCurve roc_curve( classified_results);

        // output the roc curve to a file
        io::OFStream write_roc;
        io::File::MustOpenOFStream( write_roc, FILE_NAME + m_WriteRocFlag->GetFirstParameter()->GetValue());
        roc_curve.GetThinnedRocCurveTotal( 10000).WriteRatePlottingTable( write_roc);
        io::File::CloseClearFStream( write_roc);

        // end
        return storage::List< storage::Table< double> >();
      }

      // initialize split iterators
      storage::Table< double>::const_iterator row_itr_split( current_table.Begin()), row_itr_end( current_table.End());
      for
      (
        ;
        row_itr_split != row_itr_end &&
        sp_unary_classifier->operator()( row_itr_split->Second()( criteria_index));
        ++row_itr_split
      );

      // if all entries satisfied the cutoff (larger than cutoff with Less used) then skip this table
      if( row_itr_split == row_itr_end)
      {
        BCL_MessageCrt
        (
          "could not find any rows with criteria larger than given criteria cutoff: " + util::Format()( criteria_cutoff)
        );

        // end
        return storage::List< storage::Table< double> >();
      }

      // if there were no entries satisfying the cutoff (smaller than cutoff with Less used) then skip this table
      if( row_itr_split == current_table.Begin())
      {
        BCL_MessageCrt
        (
          "could only find structures with criteria smaller than given criteria cutoff: " + util::Format()( criteria_cutoff)
        );

        // end
        return storage::List< storage::Table< double> >();
      }

      BCL_MessageStd( "finished reading score table: " + FILE_NAME);

      storage::Table< double> true_entries( current_table.Begin(), row_itr_split);
      storage::Table< double> false_entries( row_itr_split, row_itr_end);
      storage::Vector< storage::Pair< std::string, storage::Row< double> > >
        true_shuffle( true_entries.Begin(), true_entries.End());
      true_shuffle.Shuffle();
      storage::List< storage::Pair< std::string, storage::Row< double> > > tsl
      (
        true_shuffle.Begin(), true_shuffle.End()
      );
      true_entries = storage::Table< double>( tsl.Begin(), tsl.End());
      storage::Vector< storage::Pair< std::string, storage::Row< double> > >
        false_shuffle( false_entries.Begin(), false_entries.End());
      false_shuffle.Shuffle();
      storage::List< storage::Pair< std::string, storage::Row< double> > > fsl
      (
        false_shuffle.Begin(), false_shuffle.End()
      );
      false_entries = storage::Table< double>( fsl.Begin(), fsl.End());
      if( true_entries.GetSize() < size_t( 10) || false_entries.GetSize() < size_t( 10))
      {
        return storage::List< storage::Table< double> >();
      }

      return CreateCrossValidationTables( true_entries, false_entries);
    }

    storage::List< storage::Table< double> > MinimizeScoreWeightSet::CreateCrossValidationTables
    (
      const storage::Table< double> &TRUE_ENTRIES,
      const storage::Table< double> &FALSE_ENTRIES
    ) const
    {
      // initialize table to store the true entries and false entries
      const size_t nr_true_entries( TRUE_ENTRIES.GetSize());
      const size_t nr_false_entries( FALSE_ENTRIES.GetSize());

      BCL_MessageStd
      (
        "found " + util::Format()( nr_true_entries) + " true and " + util::Format()( nr_false_entries) +
        " false entries in that table and will try to balance to table target ratio: " + m_FractionParam->GetValue()
      );

      // store the target fraction
      const double target_fraction( m_FractionParam->GetNumericalValue< double>());

      // Initialize the list to store cross validation lists
      storage::List< storage::Table< double> > cross_validation_list;

      BCL_MessageStd( "creating cross validation lists");

      // store number of cross validations
      const size_t nr_cross_validations( m_NumberCrossValidationParam->GetNumericalValue< size_t>());

      // if too many false entries
      if( double( nr_true_entries) / double( nr_true_entries + nr_false_entries) < target_fraction)
      {
        const size_t target_nr_false_entries( size_t( double( nr_true_entries) / target_fraction) - nr_true_entries);

        // are there enough false entries?
        if( target_nr_false_entries < nr_cross_validations)
        {
          return cross_validation_list;
        }

        // will the table contain enough entries
        if( target_nr_false_entries + nr_true_entries < m_MinimalTableSizeParam->GetNumericalValue< size_t>())
        {
          return cross_validation_list;
        }

        // initialize false entries
        cross_validation_list = FALSE_ENTRIES.CreateDifferentSubTables( nr_cross_validations, target_nr_false_entries);

        // append true to each false entries in cross validation list
        for
        (
          storage::List< storage::Table< double> >::iterator
            itr( cross_validation_list.Begin()), itr_end( cross_validation_list.End());
          itr != itr_end; ++itr
        )
        {
          itr->Append( TRUE_ENTRIES);
        }
      }
      // if too many true entries
      else
      {
        // the number of true entries we should have after removing the excess ones
        size_t target_nr_true_entries( size_t( double( nr_false_entries) / ( 1.0 - target_fraction)) - nr_false_entries);

        // are their enough true entries to make all tables?
        if( target_nr_true_entries < nr_cross_validations)
        {
          return cross_validation_list;
        }

        // make sure that target_nr_true_entries is less than nr_true_entries - (number_cross_validation - 1)
        // otherwise it wouldn't be able to make the specified number of distinct subtables
        if( nr_true_entries - target_nr_true_entries < nr_cross_validations - 1)
        {
          // update number true entries to have a smaller size
          target_nr_true_entries = nr_true_entries - nr_cross_validations - 1;
        }

        // will the table contain enough entries
        if( target_nr_true_entries + nr_false_entries < m_MinimalTableSizeParam->GetNumericalValue< size_t>())
        {
          return cross_validation_list;
        }

        // initialize the true entries
        cross_validation_list = TRUE_ENTRIES.CreateDifferentSubTables( nr_cross_validations, target_nr_true_entries);

        // append false to each true entries in cross validation list
        for
        (
          storage::List< storage::Table< double> >::iterator
            itr( cross_validation_list.Begin()), itr_end( cross_validation_list.End());
          itr != itr_end; ++itr
        )
        {
          itr->Append( FALSE_ENTRIES);
        }
      }

      return cross_validation_list;
    }

    //! @brief create objective function that combines the enrichments for multiple tables by balancing and summing
    //! @param TABLE_FILE_NAMES list of filenames to each of the tables to be considered
    //! @param CRITERIA_NAME    name of the column that contains the criteria for the enrichment
    //! @param WEIGHT_RANGE_TABLE table of starting weights and ranges - passed in to ensure that all tables have the
    //!        same header and we have weights for each of their columns
    //! @param CONSENSUS_MAP map that will store the ConsensusEnrichment for each cross validation for each table file
    //! @return sum enrichment function for balanced tables
    util::ShPtr< util::FunctionInterface< linal::Vector< double>, double> > MinimizeScoreWeightSet::EnrichmentObjective
    (
      const storage::Vector< std::string> &TABLE_FILE_NAMES,
      const std::string &CRITERIA_NAME,
      const storage::Table< double> &WEIGHT_RANGE_TABLE,
      storage::Map< std::string, util::ShPtr< FunctionCombine< linal::Vector< double>, double> > > &CONSENSUS_MAP,
      math::RunningAverageSD< storage::Vector< double> > &AVE_SD
    ) const
    {
      BCL_MessageStd
      (
        "reading all tables from given listfile: " + m_ScoreTableFileList->GetFirstParameter()->GetValue()
      );

      storage::List
      <
        storage::Triplet
        <
          util::ShPtr< sched::JobInterface>,
          std::string,
          storage::List< storage::Table< double> >
        >
      > jobs_results;

      // iterate over individual table file names
      for
      (
        storage::Vector< std::string>::const_iterator
          filename_itr( TABLE_FILE_NAMES.Begin()), filename_itr_end( TABLE_FILE_NAMES.End());
        filename_itr != filename_itr_end;
        ++filename_itr
      )
      {
        // initialize table to store the true entries and false entries
        jobs_results.PushBack
        (
          storage::Triplet
          <
            util::ShPtr< sched::JobInterface>,
            std::string,
            storage::List< storage::Table< double> >
          >()
        );

        // reference on current job
        storage::Triplet
        <
          util::ShPtr< sched::JobInterface>,
          std::string,
          storage::List< storage::Table< double> >
        > &current_job( jobs_results.LastElement());

        current_job.Second() = *filename_itr;
        current_job.First() = util::ShPtr< sched::JobInterface>
        (
          new sched::BinaryFunctionJobWithData< const std::string, const storage::Table< double>, storage::List< storage::Table< double> >, MinimizeScoreWeightSet>
          (
            0,
            *this,
            &MinimizeScoreWeightSet::CreateCrossValidationTablesFromFile,
            *filename_itr,
            WEIGHT_RANGE_TABLE,
            sched::JobInterface::e_READY,
            &current_job.Third()
          )
        );

        // submit job
        sched::GetScheduler().SubmitJob( current_job.First());
      }

      // create sum function object
      util::ShPtr< sched::SumFunction< linal::Vector< double>, double> >
        objective_function( new sched::SumFunction< linal::Vector< double>, double>());

      const double criteria_cutoff( m_CriteriaCutoffParam->GetNumericalValue< double>());
      const math::Comparisons< double>::Comparison comp_func( m_SortOrderFlag->GetFirstParameter()->GetValue());
      const util::ShPtr< math::FunctionInterfaceSerializable< double, bool> > sp_unary_classifier
      (
        math::Comparisons< double>::GetEnums().CreateUnaryPredicate( comp_func, criteria_cutoff).Clone()
      );

      // iterate over jobs and join
      for
      (
        storage::List
        <
          storage::Triplet
          <
            util::ShPtr< sched::JobInterface>,
            std::string,
            storage::List< storage::Table< double> >
          >
        >::iterator job_itr( jobs_results.Begin()), job_itr_end( jobs_results.End());
        job_itr != job_itr_end;
        ++job_itr
      )
      {
        sched::GetScheduler().Join( job_itr->First());

        // Initialize the list to store cross validation lists
        const storage::List< storage::Table< double> > &cross_validation_list( job_itr->Third());
        if( cross_validation_list.IsEmpty())
        {
          BCL_MessageCrt( "unable to create crossvalidation tables from: " + job_itr->Second());
          continue;
        }

        // calculate prefactor
        //const double prefactor( 1.0 / double( cross_validation_list.GetSize()));

        // initialize the FunctionCombine
        util::ShPtr< FunctionCombine< linal::Vector< double>, double> >
          sp_combine( new FunctionCombine< linal::Vector< double>, double>());

        // iterate over all tables in cross_validation_list
        for
        (
          storage::List< storage::Table< double> >::const_iterator
            itr( cross_validation_list.Begin()), itr_end( cross_validation_list.End());
          itr != itr_end; ++itr
        )
        {
          // create the ShPtr to consensus function
          util::ShPtr< score::ConsensusEnrichment> sp_consensus
          (
            new score::ConsensusEnrichment( *itr, CRITERIA_NAME, score::ConsensusEnrichment::e_Enrichment, sp_unary_classifier)
          );
          math::RunningAverageSD< storage::Vector< double> > this_table;
          for( auto itr_row( itr->Begin()), itr_row_end( itr->End()); itr_row != itr_row_end; ++itr_row)
          {
            this_table += itr_row->Second().InternalData();
          }
          AVE_SD += this_table.GetStandardDeviation();

          // create single objective function
          //*objective_function -= prefactor * score::ConsensusEnrichment( *itr, CRITERIA_NAME, criteria_cutoff, **comp_function);
          sp_combine->PushBack( sp_consensus);
        }
        // insert ShPtr to the consensus map given to this function
        CONSENSUS_MAP[ job_itr->Second()] = sp_combine;

        if( m_TargetTableIndexFlag->GetFlag())
        {
          math::FunctionAdapter< linal::Vector< double>, storage::Vector< double>, double> best_nth_from_combined
          (
            sp_combine,
            util::ShPtr< math::FunctionInterfaceSerializable< storage::Vector< double>, double> >
            (
              new NthBest< double>( m_TargetTableIndexFlag->GetFirstParameter()->GetNumericalValue< size_t>())
            )
          );
          objective_function->NewOperand( best_nth_from_combined, -1.0);
        }
        else
        {
          math::FunctionAdapter< linal::Vector< double>, storage::Vector< double>, double> sum_combined
          (
            sp_combine,
            util::ShPtr< math::FunctionInterfaceSerializable< storage::Vector< double>, double> >
            (
              new SumVector< double>()
            )
          );
          objective_function->NewOperand( sum_combined, -1.0);
        }
        BCL_MessageStd( "finished creating objective functions");
      }

      BCL_MessageStd
      (
        "finished reading all tables from given list file: " + m_ScoreTableFileList->GetFirstParameter()->GetValue()
      );

      // end
      return objective_function;
    }

    //! @brief create objective function that calculates the ROCIntergral for one table and a given criteria
    //! @param TABLE_FILE_NAME file name of the table to be used
    //! @param CRITERIA_NAME    name of the column that contains the criteria for the enrichment
    //! @param WEIGHT_RANGE_TABLE table of starting weights and ranges - passed in to ensure that all tables have the
    //!        same header and we have weights for each of their columns
    //! @param COMPARISON binary predicate that classifies true or false if above criteria threshold
    //! @return ROCIntergral function on the given table
    util::ShPtr< util::FunctionInterface< linal::Vector< double>, double> > MinimizeScoreWeightSet::ROCIntegralObjective
    (
      const std::string &TABLE_FILE_NAME,
      const std::string &CRITERIA_NAME,
      const storage::Table< double> &WEIGHT_RANGE_TABLE,
      const math::Comparisons< double>::Comparison &COMPARISON,
      math::RunningAverageSD< storage::Vector< double> > &AVE_SD
    ) const
    {
      const double criteria_cutoff( m_CriteriaCutoffParam->GetNumericalValue< double>());
      const math::Comparisons< double>::Comparison comp_func( m_SortOrderFlag->GetFirstParameter()->GetValue());
      const util::ShPtr< math::FunctionInterfaceSerializable< double, bool> > sp_unary_classifier
      (
        math::Comparisons< double>::GetEnums().CreateUnaryPredicate( comp_func, criteria_cutoff).Clone()
      );

      // read table
      io::IFStream score_table_stream;
      io::File::MustOpenIFStream( score_table_stream, TABLE_FILE_NAME);

      BCL_MessageStd( "reading score table: " + TABLE_FILE_NAME);
      // read table from stream including rows with similar names
      storage::Table< double> current_table;
      current_table.ReadFormatted( score_table_stream, true);
      // close stream
      io::File::CloseClearFStream( score_table_stream);

      // check if they have the same headers
      const bool same_headers( current_table.GetHeader() == WEIGHT_RANGE_TABLE.GetHeader());

      // make sure it is same or at least has all the necessary columns
      BCL_Assert
      (
        same_headers || current_table.GetHeader().HasColumns( WEIGHT_RANGE_TABLE.GetHeader()),
        "The given table has to either have all the same headers or contain all the ones used in weight optimization"
        + util::Format()( current_table.GetHeader()) + "\nvs\n" + util::Format()( WEIGHT_RANGE_TABLE.GetHeader())
      );

      // if they were not the same then extract the required sub-table
      if( !same_headers)
      {
        current_table = current_table.ExtractSubTableFromColumns( WEIGHT_RANGE_TABLE.GetHeader());
      }

      for( auto itr( current_table.Begin()), itr_end( current_table.End()); itr != itr_end; ++itr)
      {
        AVE_SD += itr->Second().InternalData();
      }

      // sort by criteria
      current_table.SortByColumn( CRITERIA_NAME, **COMPARISON);

      if( m_WriteRocFlag->GetFlag())
      {
        // index of the criteria column
        const size_t criteria_index( WEIGHT_RANGE_TABLE.GetHeader()[ CRITERIA_NAME]);

        // get the weights as a vector
        const linal::Vector< double> weights( WEIGHT_RANGE_TABLE[ "weights"].GetData());

        // initialize list for classified results
        storage::List< storage::Pair< double, bool> > classified_results;

        // iterate over the score table
        for
        (
          storage::Table< double>::const_iterator itr( current_table.Begin()), itr_end( current_table.End());
          itr != itr_end;
          ++itr
        )
        {
          // get the vector of scores
          const linal::Vector< double> scores( itr->Second().GetData());

          // create a new classification for this row and insert it into classified_results
          classified_results.PushBack
          (
            storage::Pair< double, bool>
            (
              scores * weights,
              sp_unary_classifier->operator()( itr->Second()( criteria_index))
            )
          );
        }

        // sort the results by the predicate
        classified_results.Sort
        (
          storage::PairBinaryPredicateFirst< double, bool>
          (
            util::BinaryFunctionSTLWrapper< std::less< double> >()
          )
        );

        // create a ROC curve
        const math::ROCCurve roc_curve( classified_results);

        // first, get the last element of the curve, which is needed for computing the
        const math::ROCCurve::Point &roc_end_point( roc_curve.GetSortedCounts().LastElement());

        std::pair< storage::Vector< math::ROCCurve::Point>::const_iterator, double>
          best_cutoff_and_metric_value_acc
          (
            roc_curve.GetMaxima( math::ContingencyMatrixMeasures( math::ContingencyMatrixMeasures::e_Accuracy))
          );
        std::pair< storage::Vector< math::ROCCurve::Point>::const_iterator, double>
          best_cutoff_and_metric_value_mcc
          (
            roc_curve.GetMaxima( math::ContingencyMatrixMeasures( math::ContingencyMatrixMeasures::e_MatthewsCorrelationCoefficient))
          );
        const math::ContingencyMatrix best_contingency_matrix_acc
        (
          best_cutoff_and_metric_value_acc.first->GetContingencyMatrix( roc_end_point)
        );
        const math::ContingencyMatrix best_contingency_matrix_mcc
        (
          best_cutoff_and_metric_value_mcc.first->GetContingencyMatrix( roc_end_point)
        );
        const float best_metric_value_acc( best_cutoff_and_metric_value_acc.second);
        const float best_metric_value_mcc( best_cutoff_and_metric_value_mcc.second);
        const float best_cutoff_acc( best_cutoff_and_metric_value_acc.first->GetCutoff());
        const float best_cutoff_mcc( best_cutoff_and_metric_value_mcc.first->GetCutoff());

        std::ostringstream output_info;
        output_info << "Best Accuracy Contingency Matrix TP/TN/FP/FN/Accuracy: "
                    << best_contingency_matrix_acc.GetNumberTruePositives() << ' '
                    << best_contingency_matrix_acc.GetNumberTrueNegatives() << ' '
                    << best_contingency_matrix_acc.GetNumberFalsePositives() << ' '
                    << best_contingency_matrix_acc.GetNumberFalseNegatives() << ' '
                    << best_metric_value_acc
                    << " @cutoff " << best_cutoff_acc << '\n'
                    << "Best MCC Contingency Matrix TP/TN/FP/FN/MCC: "
                    << best_contingency_matrix_mcc.GetNumberTruePositives() << ' '
                    << best_contingency_matrix_mcc.GetNumberTrueNegatives() << ' '
                    << best_contingency_matrix_mcc.GetNumberFalsePositives() << ' '
                    << best_contingency_matrix_mcc.GetNumberFalseNegatives() << ' '
                    << best_metric_value_mcc
                    << " @cutoff " << best_cutoff_mcc << '\n';
        BCL_MessageStd( "roc integral: " + util::Format()( roc_curve.Integral()));
        BCL_MessageStd( output_info.str());

        // write the ROC curve
        io::OFStream write_roc;
        io::File::MustOpenOFStream( write_roc, TABLE_FILE_NAME + m_WriteRocFlag->GetFirstParameter()->GetValue());
        roc_curve.GetThinnedRocCurveTotal( 10000).WriteRatePlottingTable( write_roc);
        io::File::CloseClearFStream( write_roc);

        // end
        return util::ShPtr< math::FunctionInterfaceSerializable< linal::Vector< double>, double> >();
      }

      // create objective function
      util::ShPtr< math::SumFunctionMixin< math::FunctionInterfaceSerializable< linal::Vector< double>, double> > > objective_function
      (
        new math::SumFunctionMixin< math::FunctionInterfaceSerializable< linal::Vector< double>, double> >()
      );

      // initialize the objective function with the negative of the enrichment
      *objective_function -=
        score::ConsensusEnrichment
        (
          current_table, CRITERIA_NAME, score::ConsensusEnrichment::e_ROCIntegeral, sp_unary_classifier
        );

      // end
      return objective_function;
    }

    //! @brief prints out the individual enrichments in a table format
    //! @param WEIGHTS vector of weights
    //! @param CONSENSUS_MAP map that will store the ConsensusEnrichment for each cross validation for each table file
    void MinimizeScoreWeightSet::PrintIndividualEnrichments
    (
      const linal::Vector< double> &WEIGHTS,
      const storage::Map< std::string, util::ShPtr< FunctionCombine< linal::Vector< double>, double> > > &CONSENSUS_MAP
    ) const
    {
      // initialize table header strings
      storage::Vector< std::string> table_header_strings;
      for( size_t i = 0, i_max( m_NumberCrossValidationParam->GetNumericalValue< size_t>()); i < i_max; ++i)
      {
        table_header_strings.PushBack( "table_" + util::Format()( i));
      }
      // initialize table
      storage::Table< double> results_table( table_header_strings);

      // iterate over all table files in the map
      for
      (
        storage::Map< std::string, util::ShPtr< FunctionCombine< linal::Vector< double>, double> > >::const_iterator
          map_itr( CONSENSUS_MAP.Begin()), map_itr_end( CONSENSUS_MAP.End());
        map_itr != map_itr_end; ++map_itr
      )
      {
        // store the table file name
        const std::string table_file_name( map_itr->first);
        // store the results
        linal::Vector< double> results( map_itr->second->operator()( WEIGHTS));
        // sort the results
        std::sort( results.Begin(), results.End());
        // insert value into table
        results_table.InsertRow( table_file_name, storage::Vector< double>( results.Begin(), results.End()));
      }

      // print out the table
      results_table.WriteFormatted( util::GetLogger(), util::Format().W( 10).R().FFP( 3));
    }

    //! @brief default constructor
    MinimizeScoreWeightSet::MinimizeScoreWeightSet() :
      m_ScoreTableFileList
      (
        new command::FlagStatic
        (
          "list",
          "\tlist of score tables",
          command::Parameter
          (
            "scoring_tables_list_file", "file containing list of score table files",
            command::ParameterCheckFileExistence(),
            ""
          )
        )
      ),
      m_ScoreTableFilename
      (
        new command::FlagStatic
        (
          "table",
          "\tsingle table to be considered",
          command::Parameter
          (
            "scoring_table_filename",
            "file containing a score table",
            command::ParameterCheckFileExistence(),
            ""
          )
        )
      ),
      m_McMaxIterationsUnimprovedFlag
      (
        new command::FlagStatic
            (
              "mc_tot_unimproved",
              "\tmodify the number of iterations and max steps without improvement"
            )
      ),
      m_McMaxIterationsParam
      (
        new command::Parameter( "total_iterations", "\ttotal number of iterations for minimization", "100")
      ),
      m_MCMaxStepsUnimprovedParam
      (
        new command::Parameter( "steps_without_improvement", "\tnumber of steps without improvement before terminating", "50")
      ),
      m_CrossValidationFlag( new command::FlagStatic( "enrichment", "calculate enrichment for scores")),
      m_CriteriaNameParam( new command::Parameter( "criteria", "column name for criteria", "criteria")),
      m_CriteriaCutoffParam( new command::Parameter( "cutoff", "cutoff for classification", "0.5")),
      m_FractionParam( new command::Parameter( "fraction", "define fraction of true vs. flase", command::ParameterCheckRanged< double>( 0.0, 1.0), "0.5")),
      m_NumberCrossValidationParam
      (
        new command::Parameter
        (
          "number_cross_validation",
          "to achieve the fraction either true or false entries get removed systematically. To vary that and later calculate a standard deviation, the enrichment calculation will be repeated n times",
          command::ParameterCheckRanged< size_t>( 1, 100), "10"
        )
      ),
      m_MinimalTableSizeParam
      (
        new command::Parameter
        (
          "minimal_table_size",
          "the minimal number of entries per table, so that the minimzation does not specialize too much",
          command::ParameterCheckRanged< size_t>( 1, 10000), "10"
        )
      ),
      m_WeightSetStartFlag
      (
        new command::FlagStatic
        (
          "weight_set",
          "file containing bcl::storage::Table<double> with two header identical to given tables and one col with \"weights\" and a second with \"ranges\"",
          command::Parameter( "weight_set_file_name", "filename of file containing weight set and mutate ranges", "weights.table")
        )
      ),
      m_WeightSetWriteFlag
      (
        new command::FlagStatic
        (
          "weight_set_write",
          "file to which the optimized weightset will be written to",
          command::Parameter
          (
            "weight_set_write_prefix",
            "prefix for filename to which the optimized weightset will be written to, repeat number and .weights extensions will be added",
            "minimized_"
          )
        )
      ),
      m_UsePercentageRangesFlag
      (
        new command::FlagStatic
        (
          "use_percentage_ranges", "flag of using the given ranges in the weightset file as percentages rather than absolute mutate ranges"
        )
      ),
      m_WriteRocFlag
      (
        new command::FlagStatic
        (
          "write_roc",
          "uses the given weight set to calculate the roc curve and prints them to a .roc file for each of the given tables in the list file",
          command::Parameter( "roc_file_extension", "add the extension to the output roc file", ".roc")
        )
      ),
      m_SortOrderFlag
      (
        new command::FlagStatic
        (
          "sort_order",
          "defines the order the criteria is sorted - less means that small values come first and considered true",
          command::Parameter
          (
            "order",
            "a selection of order functions",
            command::ParameterCheckEnumerate< math::Comparisons< double> >(),
            math::Comparisons< double>::GetEnums().e_Less.GetName()
          )
        )
      ),
      m_NumberRepeatsFlag
      (
        new command::FlagStatic
        (
          "number_repeats",
          "\tdefines the number of repetitions the minimizations should be run",
          command::Parameter
          (
            "number_repeats",
            "the number of repetitions the minimizations should be run",
            command::ParameterCheckRanged< size_t>( 1, 100),
            "1"
          )
        )
      ),
      m_TargetTableIndexFlag
      (
        new command::FlagStatic
        (
          "target_table_index",
          "\tindex of the target table, (nth worst) to optimize for enrichment at each round of minimization",
          command::Parameter
          (
            "target_table_index_param",
            "\tindex of the target table, (nth worst) to optimize for enrichment at each round of minimization, should be less than total number of tables",
            command::ParameterCheckRanged< size_t>( 0, 99),
            "0"
          )
        )
      ),
      m_NumberWeightsMutatedInOneStepFlag
      (
        new command::FlagStatic
        (
          "number_weights_mutated",
          "\tnumber of weights altered in one iteration",
          command::Parameter
          (
            "number_weights_mutated_param",
            "\tnumber of weights the are modified in by the given ranges in the weight table",
            command::ParameterCheckRanged< size_t>( 1, 99),
            "2"
          )
        )
      ),
      m_KeepWeightsPositiveFlag
      (
        new command::FlagStatic
        (
          "keep_positive",
          "keep the weight set positive during minimization"
        )
      ),
      m_ConstantFlag
      (
        new command::FlagStatic
        (
          "constant",
          "name of the score to keep constant during minimization; to prevent scoreweight blowup",
          command::Parameter( "constant_score", "a score name", "")
        )
      )
    {
      // attach parameters to flags
      m_McMaxIterationsUnimprovedFlag->PushBack( m_McMaxIterationsParam);
      m_McMaxIterationsUnimprovedFlag->PushBack( m_MCMaxStepsUnimprovedParam);
      m_CrossValidationFlag->PushBack( m_CriteriaNameParam);
      m_CrossValidationFlag->PushBack( m_CriteriaCutoffParam);
      m_CrossValidationFlag->PushBack( m_FractionParam);
      m_CrossValidationFlag->PushBack( m_NumberCrossValidationParam);
      m_CrossValidationFlag->PushBack( m_MinimalTableSizeParam);
    }

    const ApplicationType MinimizeScoreWeightSet::MinimizeScoreWeightSet_Instance
    (
      GetAppGroups().AddAppToGroup( new MinimizeScoreWeightSet(), GetAppGroups().e_Bcl)
    );

  } // namespace app
} // namespace bcl
