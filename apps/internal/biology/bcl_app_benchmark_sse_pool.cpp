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
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_factory_mc.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_file_existence.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "io/bcl_io_file.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "sspred/bcl_sspred_method_handler.h"
#include "sspred/bcl_sspred_sse_factory_threshold.h"
#include "storage/bcl_storage_template_instantiations.h"

namespace bcl
{
  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class BenchmarkSSEPool
    //! @brief This application benchmarks SSEPool generation by trying all different combinations specified
    //! @details This class benchmarks different ways of SSEPool generation. It tries each combination of following
    //! parameters specified in the command line
    //! - secondary structure prediction methods (also tries internal combinations)
    //! - use consensus prediction or not
    //! - chop long SSEs or not
    //! - minimum helix length from pool( all combinations specified by given range and number of intervals)
    //! - minimum strand length from pool( all combinations specified by given range and number of intervals)
    //! - minimum helix prediction threshold( all combinations specified by given range and number of intervals)
    //! - minimum strand prediction threshold( all combinations specified by given range and number of intervals)
    //! For each combination, the program iterates over given list of pdbs and creates a statistics table that
    //! have a variety of measures that measure how accurate the predicted pool is with respect to native SSE
    //! definitions from the native structure. It outputs these only if m_OutputAllTables flag is set.
    //! Then it generates a tag that describes the values of the variables in the combination and adds the averaged
    //! values for these over 54 proteins, to a table as a row with the tag as its name. The averaging also has some
    //! values which are weighted averaged. At the end of the benchmark, this table that contains average and weighted
    //! accuracy measures for all combinations to a specified file.
    //! The values in tags are separated by plus signs, so they can easily be converted to columns in Excel. The invidual
    //! methods within that combination are separated by underscore signs.
    //!
    //! @author karakam
    //! @date 03/22/2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BenchmarkSSEPool :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! prefix where to secondary structure predictions should be found
      util::ShPtr< command::FlagStatic> m_PdbListFlag;
      util::ShPtr< command::Parameter> m_PdbPathParam;
      util::ShPtr< command::Parameter> m_Pdb5LetterCodeWithChainsFileParam;
      util::ShPtr< command::Parameter> m_PdbPathHierarchyParam;

      //! output path flag
      util::ShPtr< command::FlagInterface> m_OutputPrefixFlag;

      //! ss methods to be used in pool generation
      util::ShPtr< command::FlagInterface> m_SsMethodsFlag;

      //! ranges for sse thresholds
      util::ShPtr< command::FlagStatic> m_SSEThresholdRangesFlag;
      util::ShPtr< command::ParameterInterface> m_SSEThresholdMinParam;
      util::ShPtr< command::ParameterInterface> m_SSEThresholdMaxParam;
      util::ShPtr< command::ParameterInterface> m_SSEThresholdNumberIntervalsParam;

      //! ranges for helix pool lengths
      util::ShPtr< command::FlagStatic> m_PoolLengthHelixFlag;
      util::ShPtr< command::ParameterInterface> m_PoolLengthHelixMinParam;
      util::ShPtr< command::ParameterInterface> m_PoolLengthHelixMaxParam;

      //! ranges for strand pool lengths
      util::ShPtr< command::FlagStatic> m_PoolLengthStrandFlag;
      util::ShPtr< command::ParameterInterface> m_PoolLengthStrandMinParam;
      util::ShPtr< command::ParameterInterface> m_PoolLengthStrandMaxParam;

      //! flag whether all the pools should be output
      util::ShPtr< command::FlagInterface> m_OutputAllTablesFlag;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      BenchmarkSSEPool();

    public:

      //! @brief Clone function
      //! @return pointer to new BenchmarkSSEPool
      BenchmarkSSEPool *Clone() const
      {
        return new BenchmarkSSEPool( *this);
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

      //! @brief function to calculate the sum row for a single combination from given table
      //! @param TABLE table of interest
      //! @return vector of summed values that
      storage::Vector< double> CalculateSumRow( const storage::Table< double> &TABLE) const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const
      {
        // initialize ShPtr to a new command object
        util::ShPtr< command::Command> sp_cmd( new command::Command());

        // prefix where to secondary structure predictions should be found
        sp_cmd->AddFlag( m_PdbListFlag);

        // output path flag
        sp_cmd->AddFlag( m_OutputPrefixFlag);
        // ss methods to be used in pool generation
        sp_cmd->AddFlag( m_SsMethodsFlag);

        // ranges for sse thresholds
        sp_cmd->AddFlag( m_SSEThresholdRangesFlag);

        // ranges for helix pool lengths
        sp_cmd->AddFlag( m_PoolLengthHelixFlag);

        // ranges for strand pool lengths
        sp_cmd->AddFlag( m_PoolLengthStrandFlag);

        // flag whether all the tables should be output
        sp_cmd->AddFlag( m_OutputAllTablesFlag);

        // adjust minimal sse lengths
        pdb::Factory::GetFlagMinSSESize()->GetParameterList()( biol::GetSSTypes().HELIX)->SetDefaultParameter( "7");
        pdb::Factory::GetFlagMinSSESize()->GetParameterList()( biol::GetSSTypes().STRAND)->SetDefaultParameter( "4");
        pdb::Factory::GetFlagMinSSESize()->GetParameterList()( biol::GetSSTypes().COIL)->SetDefaultParameter( "999");
        sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());

        // add default bcl parameters
        command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

        // return assembled Command object
        return sp_cmd;
      }

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

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      static const ApplicationType BenchmarkSSEPool_Instance;

    }; // class BenchmarkSSEPool

    //! @brief function to calculate the sum row for a single combination from given table
    //! @param TABLE table of interest
    //! @return vector of summed values that
    storage::Vector< double> BenchmarkSSEPool::CalculateSumRow( const storage::Table< double> &TABLE) const
    {
      // construct vector
      storage::Vector< double> sum_vector( TABLE.GetNumberColumns(), 0.0);

      // create reference on the header
      const storage::TableHeader &header( TABLE.GetHeader());

      // iterate over rows
      for
      (
        storage::List< storage::Pair< std::string, storage::Row< double> > >::const_iterator
          row_itr( TABLE.Begin()), row_itr_end( TABLE.End());
        row_itr != row_itr_end; ++row_itr
      )
      {
        // create reference on the data for this row
        const storage::Vector< double> &this_row( row_itr->Second().GetData());

        // sum up values
        std::transform
        (
          sum_vector.Begin(),
          sum_vector.End(),
          this_row.Begin(),
          sum_vector.Begin(),
          std::plus< double>()
        );
      }

      // store the counts
      // remember indices for things needed for weighted average calculation
      const double nr_identified_sses( sum_vector( header[ "nr_identified_sses"]));
      const double nr_identified_helix( sum_vector( header[ "nr_identified_helix"]));
      const double nr_identified_strand( sum_vector( header[ "nr_identified_strand"]));
      const double nr_overlaps( sum_vector( header[ "nr_overlaps"]));
      const double nr_helix_overlaps( sum_vector( header[ "nr_helix_overlaps"]));
      const double nr_strand_overlaps( sum_vector( header[ "nr_strand_overlaps"]));
      const double nr_rows( TABLE.GetNumberRows());

      // divide the counts by number rows for values upto the column where sum start
      const size_t index_sums_start( header[ "best_sum_overlap"]);
      std::transform
      (
        sum_vector.Begin(),
        sum_vector.Begin() + index_sums_start,
        sum_vector.Begin(),
        std::bind2nd( std::divides< double>(), nr_rows)
      );

      // create itr
      storage::Vector< double>::iterator itr( sum_vector.Begin() + index_sums_start);

      // for the sum of measures normalize by corresponding summed counts instead of number rows
      // sums for best matches
      *itr = ( nr_identified_sses   == 0.0 ? 0.0 : *itr / nr_identified_sses);   ++itr; // best_sum_overlap
      *itr = ( nr_identified_helix  == 0.0 ? 0.0 : *itr / nr_identified_helix);  ++itr; // best_sum_helix_overlap
      *itr = ( nr_identified_strand == 0.0 ? 0.0 : *itr / nr_identified_strand); ++itr; // best_sum_strand_overlap
      *itr = ( nr_identified_sses   == 0.0 ? 0.0 : *itr / nr_identified_sses);   ++itr; // best_sum_shift
      *itr = ( nr_identified_helix  == 0.0 ? 0.0 : *itr / nr_identified_helix);  ++itr; // best_sum_helix_shift
      *itr = ( nr_identified_strand == 0.0 ? 0.0 : *itr / nr_identified_strand); ++itr; // best_sum_strand_shift
      *itr = ( nr_identified_sses   == 0.0 ? 0.0 : *itr / nr_identified_sses);   ++itr; // best_sum_length_dev
      *itr = ( nr_identified_helix  == 0.0 ? 0.0 : *itr / nr_identified_helix);  ++itr; // best_sum_helix_length_dev
      *itr = ( nr_identified_strand == 0.0 ? 0.0 : *itr / nr_identified_strand); ++itr; // best_sum_strand_length_dev

      // sums for all overlaps
      *itr = ( nr_overlaps        == 0.0 ? 0.0 : *itr / nr_overlaps);        ++itr; // all_sum_overlap
      *itr = ( nr_helix_overlaps  == 0.0 ? 0.0 : *itr / nr_helix_overlaps);  ++itr; // all_sum_helix_overlap
      *itr = ( nr_strand_overlaps == 0.0 ? 0.0 : *itr / nr_strand_overlaps); ++itr; // all_sum_strand_overlap
      *itr = ( nr_overlaps        == 0.0 ? 0.0 : *itr / nr_overlaps);        ++itr; // all_sum_overlap
      *itr = ( nr_helix_overlaps  == 0.0 ? 0.0 : *itr / nr_helix_overlaps);  ++itr; // all_sum_helix_overlap
      *itr = ( nr_strand_overlaps == 0.0 ? 0.0 : *itr / nr_strand_overlaps); ++itr; // all_sum_strand_overlap
      *itr = ( nr_overlaps        == 0.0 ? 0.0 : *itr / nr_overlaps);        ++itr; // all_sum_overlap
      *itr = ( nr_helix_overlaps  == 0.0 ? 0.0 : *itr / nr_helix_overlaps);  ++itr; // all_sum_helix_overlap
      *itr = ( nr_strand_overlaps == 0.0 ? 0.0 : *itr / nr_strand_overlaps); ++itr; // all_sum_strand_overlap

      // end
      return sum_vector;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int BenchmarkSSEPool::Main() const
    {
      // initialize write and read stream objects
      io::OFStream write;
      io::IFStream read;

      // get the list of methods
      const storage::Set< sspred::Method> ss_methods( m_SsMethodsFlag->GetObjectSet< sspred::Method>());

      // store data prefix
      const std::string prefix( m_PdbPathParam->GetValue() + PATH_SEPARATOR);
      const bool pdb_hierarchy( m_PdbPathHierarchyParam->GetNumericalValue< bool>());

      // construct maps to store the sequences and the pdbs
      storage::Map< std::string, assemble::ProteinModel> models_map;

      {
        // read the pdb list
        const std::string pdb_list_filename( m_Pdb5LetterCodeWithChainsFileParam->GetValue());
        io::File::MustOpenIFStream( read, pdb_list_filename);
        const storage::Vector< std::string> pdb_list( util::StringListFromIStream( read));
        io::File::CloseClearFStream( read);
        BCL_MessageCrt( util::Format()( pdb_list.GetSize()) + " pdbs found");

        // construct a set to ensure they are not redundant
        const storage::Set< std::string> pdb_list_set( pdb_list.Begin(), pdb_list.End());
        BCL_Assert
        (
          pdb_list.GetSize() == pdb_list_set.GetSize(),
          "The given pdb ids are not unique, " + util::Format()( pdb_list.GetSize() - pdb_list_set.GetSize()) +
          " of them are redundant!"
        );

        // iterate over all the pdbs given in the list
        for
        (
          storage::Vector< std::string>::const_iterator pdb_itr( pdb_list.Begin()), pdb_itr_end( pdb_list.End());
          pdb_itr != pdb_itr_end; ++pdb_itr
        )
        {
          // read the fasta and construct an AASequence assume the first 4 characters need to be in lower case
          std::string pdb_id( *pdb_itr);
          std::transform( pdb_id.begin(), pdb_id.end() - 1, pdb_id.begin(), tolower);

          // read the native structure
          const std::string this_prefix( prefix + ( pdb_hierarchy ? ( pdb_id.substr( 1, 2) + PATH_SEPARATOR) : ""));
          const std::string pdb_filename( this_prefix + pdb_id + ".pdb");
          pdb::Factory pdb_factory( biol::GetAAClasses().e_AABackBone);

          io::File::MustOpenIFStream( read, pdb_filename);
          pdb::Handler pdb_handler( read);
          io::File::CloseClearFStream( read);

          // read in the sequences from pdb file
          models_map[ pdb_id] = pdb_factory.ProteinModelFromPDB( pdb_handler);
          BCL_MessageCrt( "Reading the ss methods");
          BCL_Assert
          (
            sspred::MethodHandler::ReadPredictionsForProteinModel( ss_methods, models_map[ pdb_id], pdb_id, this_prefix),
            "reading sspredictions have failed!"
          );
        }
      }

    /////////////////////////////////
    // combinations initialization //
    /////////////////////////////////

      // now construct the combinations for each parameter
      const storage::Vector< size_t> chop_sses_combine( storage::Vector< size_t>::Create( 0, 1));
      storage::Vector< double> helix_thresholds;
      storage::Vector< double> strand_thresholds;
      storage::Vector< size_t> helix_pool_lengths;
      storage::Vector< size_t> strand_pool_lengths;

      // get the min max and nr_intervals
      const double threshold_min( m_SSEThresholdMinParam->GetNumericalValue< double>());
      const double threshold_max( m_SSEThresholdMaxParam->GetNumericalValue< double>());
      const double threshold_nr_intervals( m_SSEThresholdNumberIntervalsParam->GetNumericalValue< size_t>());

      // make sure max is larger than min
      BCL_Assert( threshold_min <= threshold_max, "The max threshold should be larger than/equal to min_threshold");
      // if they are the same
      if( threshold_min == threshold_max)
      {
        BCL_MessageCrt( "The min max threshold is the same, therefore using a single value");
        helix_thresholds.PushBack( threshold_min);
        strand_thresholds.PushBack( threshold_min);
      }
      else
      {
        // divide the threshold diff by the given number
        const double interval_size( ( threshold_max - threshold_min) / double( threshold_nr_intervals));

        double this_threshold( threshold_min);
        // iterate
        for( size_t i( 0); i <= threshold_nr_intervals; ++i)
        {
          helix_thresholds.PushBack( this_threshold);
          strand_thresholds.PushBack( this_threshold);
          this_threshold += interval_size;
        }
      }

      // min helix lengths for the pool
      const size_t pool_helix_range_min( m_PoolLengthHelixMinParam->GetNumericalValue< size_t>());
      const size_t pool_helix_range_max( m_PoolLengthHelixMaxParam->GetNumericalValue< size_t>());
      BCL_Assert
      (
        pool_helix_range_min <= pool_helix_range_max,
        "The min length should be smaller than max length for helices pool lengths"
      );
      // iterate
      for( size_t i( pool_helix_range_min); i <= pool_helix_range_max; ++i)
      {
        helix_pool_lengths.PushBack( i);
      }

      // min strand lengths for the pool
      const size_t pool_strand_range_min( m_PoolLengthStrandMinParam->GetNumericalValue< size_t>());
      const size_t pool_strand_range_max( m_PoolLengthStrandMaxParam->GetNumericalValue< size_t>());
      BCL_Assert
      (
        pool_strand_range_min <= pool_strand_range_max,
        "The min length should be smaller than max length for helices pool lengths"
      );
      // iterate
      for( size_t i( pool_strand_range_min); i <= pool_strand_range_max; ++i)
      {
        strand_pool_lengths.PushBack( i);
      }

      // calculate total number of combinations
      const size_t nr_combinations
      (
        ss_methods.GetSize() *
        chop_sses_combine.GetSize() *
        helix_thresholds.GetSize() *
        strand_thresholds.GetSize() *
        helix_pool_lengths.GetSize() *
        strand_pool_lengths.GetSize()
      );

      BCL_MessageCrt
      (
        " starting benchmark on " + util::Format()( nr_combinations) +
        " combinations x " + util::Format()( models_map.GetSize()) + " proteins"
      );
      BCL_MessageCrt( "The number of combinations for each variable");
      BCL_MessageCrt
      (
        " ss_methods:\t" + util::Format()( ss_methods.GetSize())
      );
      BCL_MessageCrt
      (
        " chop_sse:\t" + util::Format()( chop_sses_combine.GetSize())
      );
      BCL_MessageCrt
      (
        " helix_tresholds :\t" + util::Format()( helix_thresholds.GetSize())
      );
      BCL_MessageCrt
      (
        " strand_tresholds:\t" + util::Format()( strand_thresholds.GetSize())
      );
      BCL_MessageCrt
      (
        " helix_lengths:\t" + util::Format()( helix_pool_lengths.GetSize())
      );
      BCL_MessageCrt
      (
        " strand_lengths:\t " + util::Format()( strand_pool_lengths.GetSize()));

      // construct a table header to be used
      util::ShPtr< storage::TableHeader> sp_table_header
      (
        new storage::TableHeader( assemble::SSEPool::GetStatisticsTableHeaders())
      );

      // start counter for combinations
      size_t combination_ctr( 0);

      // construct a table to hold the summed statistics from each combination
      storage::Table< double> stats_table( sp_table_header);
      storage::Table< double> stats_table_mc( sp_table_header);

      // set output prefix
      const std::string output_prefix( m_OutputPrefixFlag->GetFirstParameter()->GetValue());

    //////////////////
    // ssefactorymc //
    //////////////////

      // methods
      for
      (
        storage::Set< sspred::Method>::const_iterator ssitr( ss_methods.Begin()), ssitr_end( ss_methods.End());
        ssitr != ssitr_end;
        ++ssitr
      )
      {
        // construct factory
        const double confidence_threshold( 0.5);
        const assemble::SSEFactoryMC ss_factory( *ssitr, confidence_threshold);

        BCL_MessageVrb( "processing method: " + ssitr->GetName() + " for factory mc");

        storage::Map< std::string, assemble::SSEPool> pdb_pools;

        // iterate over all the pdbs
        for
        (
          storage::Map< std::string, assemble::ProteinModel>::const_iterator
            pdb_itr( models_map.Begin()), pdb_itr_end( models_map.End());
          pdb_itr != pdb_itr_end; ++pdb_itr
        )
        {
          BCL_MessageDbg( "processing pdb: " + pdb_itr->first);

          // create sse pool and get the statistics
          pdb_pools[ pdb_itr->first] = ss_factory( *pdb_itr->second.GetChains().FirstElement()->GetSequence());
        } // pdb itr

        // helix pool length
        for
        (
          storage::Vector< size_t>::const_iterator
            helix_min_itr( helix_pool_lengths.Begin()), helix_min_itr_end( helix_pool_lengths.End());
            helix_min_itr != helix_min_itr_end; ++helix_min_itr
        )
        {
          // strand pool length
          for
          (
            storage::Vector< size_t>::const_iterator
              strand_min_itr( strand_pool_lengths.Begin()), strand_min_itr_end( strand_pool_lengths.End());
            strand_min_itr != strand_min_itr_end; ++strand_min_itr
          )
          {
            storage::Map< biol::SSType, size_t> min_sse_size;
            min_sse_size[ biol::GetSSTypes().HELIX] = *helix_min_itr;
            min_sse_size[ biol::GetSSTypes().STRAND] = *strand_min_itr;

            // tag
            const std::string this_tag
            (
              "ssefactory_mc_" + ssitr->GetName() + "+" +
              util::Format()( *helix_min_itr) + "+" +
              util::Format()( *strand_min_itr)
            );

            // construct a table to hold the statistics for this combination
            storage::Table< double> this_combination_table( sp_table_header);

            // iterate over pools
            for
            (
              storage::Map< std::string, assemble::SSEPool>::const_iterator
                pool_itr( pdb_pools.Begin()), pool_itr_end( pdb_pools.End());
              pool_itr != pool_itr_end;
              ++pool_itr
            )
            {
              assemble::SSEPool current_pool( pool_itr->second);
              current_pool.Prune( min_sse_size);

              // get the statistics and store the row
              const storage::Row< double> this_row
              (
                current_pool.CalculateStatistics( models_map[ pool_itr->first]).InternalData().FirstElement().Second()
              );

              // get the data and add to the table for this combination
              this_combination_table.InsertRow( pool_itr->first, this_row.GetData(), true);
            }

            // calculate sum statistics and insert to final table
            stats_table_mc.InsertRow( this_tag, CalculateSumRow( this_combination_table), true);

            // if output requested
            if( m_OutputAllTablesFlag->GetFlag())
            {
              const std::string this_table_filename( output_prefix + this_tag + ".table");
              io::File::MustOpenOFStream( write, this_table_filename);
              this_combination_table.WriteFormatted( write);
              io::File::CloseClearFStream( write);
            }
          }
        }
      }

      {
        // output the stats table
        const std::string stats_table_filename( output_prefix + "factory_mc_stats.table");
        BCL_MessageCrt( "Outputting the statistics table to " + stats_table_filename);
        io::File::MustOpenOFStream( write, stats_table_filename);
        stats_table_mc.WriteFormatted( write);
        io::File::CloseClearFStream( write);
      }

    ///////////////////////////////////
    // ssefactorysspred combinations //
    ///////////////////////////////////

      // methods
      for
      (
        storage::Set< sspred::Method>::const_iterator
          methods_itr( ss_methods.Begin()), methods_itr_end( ss_methods.End());
        methods_itr != methods_itr_end; ++methods_itr
      )
      {
        // construct methods tag
        std::string methods_tag( methods_itr->GetName());

        // chop
        for
        (
          storage::Vector< size_t>::const_iterator chop_itr( chop_sses_combine.Begin()), chop_itr_end( chop_sses_combine.End());
          chop_itr != chop_itr_end; ++chop_itr
        )
        {
          // create reference
          const bool chop_sses( *chop_itr);

          // helix threshold
          for
          (
            storage::Vector< double>::const_iterator
              helix_thr_itr( helix_thresholds.Begin()), helix_thr_itr_end( helix_thresholds.End());
            helix_thr_itr != helix_thr_itr_end; ++helix_thr_itr
          )
          {
            // strand threshold
            for
            (
              storage::Vector< double>::const_iterator
                strand_thr_itr( strand_thresholds.Begin()), strand_thr_itr_end( strand_thresholds.End());
              strand_thr_itr != strand_thr_itr_end; ++strand_thr_itr
            )
            {
              // thresholds map
              storage::Map< biol::SSType, double> thresholds_map;
              thresholds_map[ biol::GetSSTypes().HELIX] = *helix_thr_itr;
              thresholds_map[ biol::GetSSTypes().STRAND] = *strand_thr_itr;

              // helix pool length
              for
              (
                storage::Vector< size_t>::const_iterator
                  helix_min_itr( helix_pool_lengths.Begin()), helix_min_itr_end( helix_pool_lengths.End());
                  helix_min_itr != helix_min_itr_end; ++helix_min_itr
              )
              {
                // strand pool length
                for
                (
                  storage::Vector< size_t>::const_iterator
                    strand_min_itr( strand_pool_lengths.Begin()), strand_min_itr_end( strand_pool_lengths.End());
                  strand_min_itr != strand_min_itr_end; ++strand_min_itr
                )
                {

                  // increment combination counter
                  ++combination_ctr;

                  // if a multiple of 100 print out
                  if( combination_ctr % 100 == 0)
                  {
                    BCL_MessageStd
                    (
                      "Combination #" + util::Format()( combination_ctr) + "/" + util::Format()( nr_combinations)
                    );
                  }

                  // construct tag for this combination
                  const std::string this_tag
                  (
                    "ssefactory_sspred_" + methods_tag + "+" +
                    util::Format()( size_t( chop_sses)) + "+" +
                    util::Format()( *helix_thr_itr) + "+" +
                    util::Format()( *strand_thr_itr) + "+" +
                    util::Format()( *helix_min_itr) + "+" +
                    util::Format()( *strand_min_itr)
                  );

                  BCL_MessageStd( "tag: " + this_tag);

                  // create min pool lengths vector
                  storage::Map< biol::SSType, size_t> min_pool_lengths;
                  min_pool_lengths[ biol::GetSSTypes().HELIX] = ( *helix_min_itr);
                  min_pool_lengths[ biol::GetSSTypes().STRAND] = ( *strand_min_itr);
                  min_pool_lengths[ biol::GetSSTypes().COIL] = 0;

                  // construct SSEFactoryThreshold
                  const sspred::SSEFactoryThreshold ss_factory
                  (
                    *methods_itr, thresholds_map, false
                  );

                  // construct a table to hold the statistics for this combination
                  storage::Table< double> this_combination_table( sp_table_header);

                  // iterate over all the pdbs
                  for
                  (
                    storage::Map< std::string, assemble::ProteinModel>::const_iterator
                      pdb_itr( models_map.Begin()), pdb_itr_end( models_map.End());
                    pdb_itr != pdb_itr_end; ++pdb_itr
                  )
                  {
                    // create sse pool and get the statistics
                    assemble::SSEPool this_pool( ss_factory( *pdb_itr->second.GetChains().FirstElement()->GetSequence()));

                    // prune
                    this_pool.Prune( min_pool_lengths);

                    // chop
                    if( chop_sses)
                    {
                      this_pool.ChopSSEs( min_pool_lengths);
                    }

                    // get the statistics and store the row
                    storage::Row< double> this_row
                    (
                      this_pool.CalculateStatistics( pdb_itr->second).InternalData().FirstElement().Second()
                    );

                    // get the data and add to the table for this combination
                    this_combination_table.InsertRow( pdb_itr->first, this_row.GetData(), true);
                  } // pdb itr

                  // calculate sum statistics and insert to final table
                  stats_table.InsertRow( this_tag, CalculateSumRow( this_combination_table), true);

                  // if output requested
                  if( m_OutputAllTablesFlag->GetFlag())
                  {
                    const std::string this_table_filename( output_prefix + this_tag + ".table");
                    io::File::MustOpenOFStream( write, this_table_filename);
                    this_combination_table.WriteFormatted( write);
                    io::File::CloseClearFStream( write);
                  }

                } // strand min length
              } // helix min length
            } // strand threshold
          } // helix threshold
        } // chop
      } // methods

      BCL_MessageCrt( "All combinations have been sampled");

      // output the stats table
      const std::string stats_table_filename( output_prefix + "factory_sspred_stats.table");
      BCL_MessageCrt( "Outputting the statistics table to " + stats_table_filename);
      io::File::MustOpenOFStream( write, stats_table_filename);
      stats_table.WriteFormatted( write);
      io::File::CloseClearFStream( write);

      return 0;
    }

    //! default constructor
    BenchmarkSSEPool::BenchmarkSSEPool() :
      m_PdbListFlag
      (
        new command::FlagStatic
        (
          "pdblist",
          "\tfrom a given path and a list of pdb 5 letter codes, statistics are derived"
        )
      ),
      m_PdbPathParam
      (
        new command::Parameter
        (
          "pdb_path",
          "path of pdb files, where for each 5 letter code, there are corresponding secondary structure prediction files present"
        )
      ),
      m_Pdb5LetterCodeWithChainsFileParam
      (
        new command::Parameter( "list", "a list file of pdb 5 letter codes", command::ParameterCheckFileExistence())
      ),
      m_PdbPathHierarchyParam
      (
        new command::Parameter
        (
          "hierarchy",
          "boolean to indicate whether a pdb hierarchy is used in input paths so for pdbtag 1abc.pdb it looks at {path}/ab/1abc.pdb",
          command::ParameterCheckRanged< size_t>( 0, 1), "0"
        )
      ),
      m_OutputPrefixFlag
      (
        new command::FlagStatic
        (
          "output_prefix",
          "\tflag for providing an output prefix for files to be written such as pool",
          command::Parameter
          (
            "output_prefix_param",
            "\tfprefix for output files to be written such as pool",
            ""
          )
        )
      ),
      m_SsMethodsFlag
      (
        new command::FlagDynamic
        (
          "ssmethods", "\tone or more ssmethods to be used in generation of the pool",
          command::Parameter
          (
            "ssmethod", "any ssmethod from the list", command::ParameterCheckEnumerate< sspred::Methods>()
          ),
          1,
          sspred::GetMethods().GetEnumCount()
        )
      ),
      m_SSEThresholdRangesFlag
      (
        new command::FlagStatic
        (
          "sse_threshold", "\tthreshold to be used to identify whether a residue in SSE or not, between 0 and 1"
        )
      ),
      m_SSEThresholdMinParam
      (
        new command::Parameter
        (
          "sse_threshold_min", "\tmin threshold", command::ParameterCheckRanged< double>( 0.0, 1.0), "0.2"
        )
      ),
      m_SSEThresholdMaxParam
      (
        new command::Parameter
        (
          "sse_threshold_max", "\tmax threshold", command::ParameterCheckRanged< double>( 0.0, 1.0), "0.6"
        )
      ),
      m_SSEThresholdNumberIntervalsParam
      (
        new command::Parameter
        (
          "sse_threshold_nr_intervals", "\tnumber of intervals to be used",
          command::ParameterCheckRanged< size_t>( 0, 10), "5"
        )
      ),
      m_PoolLengthHelixFlag
      (
        new command::FlagStatic( "pool_length_helix", "min and max pool lengths for helices to be tested")
      ),
      m_PoolLengthHelixMinParam
      (
        new command::Parameter
        (
          "pool_min_length_helix", "minimum helix pool length",
          command::ParameterCheckRanged< size_t>( 2, 15), "7"
        )
      ),
      m_PoolLengthHelixMaxParam
      (
        new command::Parameter
        (
          "pool_max_length_helix", "maximum helix pool length",
          command::ParameterCheckRanged< size_t>( 2, 15), "7"
        )
      ),
      m_PoolLengthStrandFlag
      (
        new command::FlagStatic( "pool_length_strand", "min and max pool lengths for strands to be tested")
      ),
      m_PoolLengthStrandMinParam
      (
        new command::Parameter
        (
          "pool_min_length_strand", "minimum strand pool length",
          command::ParameterCheckRanged< size_t>( 2, 15), "7"
        )
      ),
      m_PoolLengthStrandMaxParam
      (
        new command::Parameter
        (
          "pool_max_length_strand", "maximum strand pool length",
          command::ParameterCheckRanged< size_t>( 2, 15), "7"
        )
      ),
      m_OutputAllTablesFlag
      (
        new command::FlagStatic
        (
          "output_all_tables",
          "flag to determine whether all tables should be output"
        )
      )
    {
      m_PdbListFlag->PushBack( m_PdbPathParam);
      m_PdbListFlag->PushBack( m_Pdb5LetterCodeWithChainsFileParam);
      m_PdbListFlag->PushBack( m_PdbPathHierarchyParam);
      m_SSEThresholdRangesFlag->PushBack( m_SSEThresholdMinParam);
      m_SSEThresholdRangesFlag->PushBack( m_SSEThresholdMaxParam);
      m_SSEThresholdRangesFlag->PushBack( m_SSEThresholdNumberIntervalsParam);
      m_PoolLengthHelixFlag->PushBack( m_PoolLengthHelixMinParam);
      m_PoolLengthHelixFlag->PushBack( m_PoolLengthHelixMaxParam);
      m_PoolLengthStrandFlag->PushBack( m_PoolLengthStrandMinParam);
      m_PoolLengthStrandFlag->PushBack( m_PoolLengthStrandMaxParam);
    }

    const ApplicationType BenchmarkSSEPool::BenchmarkSSEPool_Instance
    (
      GetAppGroups().AddAppToGroup( new BenchmarkSSEPool(), GetAppGroups().e_InternalBiol)
    );

  } // namespace app
} // namespace bcl
