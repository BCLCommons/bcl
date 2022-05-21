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
#include "assemble/bcl_assemble_protein_storage_file.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_running_min_max.h"
#include "pdb/bcl_pdb_handler.h"
#include "pdb/bcl_pdb_printer_score.h"
#include "storage/bcl_storage_table.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FoldAnalysis
    //! @brief Application for analyzing protein models generated during a folding run
    //! @details Reads protein models and calculates requested statistics
    //!
    //! @author weinerbe
    //! @date Dec 14, 2011
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API FoldAnalysis :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      util::ShPtr< command::FlagInterface> m_TableFromFileFlag;       //!< read score table from file flag
      util::ShPtr< command::FlagInterface> m_TableFromPDBListFlag;    //!< score table from pdb list flag
      util::ShPtr< command::FlagInterface> m_TableFromPrefixListFlag; //!< score table from prefix string list flag
      util::ShPtr< command::FlagInterface> m_TableFromStagesFlag;     //!< score table only from models from some stages

      util::ShPtr< command::FlagInterface> m_SortFlag;                //!< sorting flag
      util::ShPtr< command::FlagInterface> m_SortCriteriaFlag;        //!< sorting method param

      util::ShPtr< command::FlagStatic>    m_FilterFlag;              //!< filtering flag

      util::ShPtr< command::FlagInterface> m_OutputTableFlag;         //!< filename to write score table to
      util::ShPtr< command::FlagInterface> m_OutputPdbListFlag;       //!< filename to write the (filtered) pdb list to

      util::ShPtr< command::FlagStatic>         m_OutputValueFlag;    //!< output a particular value (column, row)
      util::ShPtr< command::ParameterInterface> m_ColumnParam;        //!< column to be evaluated param
      util::ShPtr< command::ParameterInterface> m_RowParam;           //!< row param
      util::ShPtr< command::ParameterInterface> m_FilenameParam;      //!< file to write information to

      util::ShPtr< command::FlagStatic> m_FlagUseUnweightedScores;    //!< flag to use unweighted scores

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      FoldAnalysis();

      //! @brief Clone function
      //! @return pointer to new FoldAnalysis
      FoldAnalysis *Clone() const;

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the Command object
      util::ShPtr< command::Command> InitializeCommand() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to std::ostream
      //! @param OSTREAM outputstream to write to
      //! @param INDENT number of indentations
      //! @return outputstream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief create a score table
      //! @return score table
      storage::Table< double> GetScoreTable() const;

      //! @param LINES pdb lines to analyze
      //! @param ROW_NAME row name
      //! @param TABLE table to write to
      //! @param USE_UNWEIGHTED_SCORES unweighted scores will be gotten
      static void AddTableRow
      (
        const util::ShPtrList< pdb::Line> &LINES,
        const std::string &ROW_NAME,
        storage::Table< double> &TABLE,
        const bool USE_UNWEIGHTED_SCORES
      );

      //! @brief sort the passed score table according to criteria passed on command line
      //! @param SCORE_TABLE the table to sort
      void SortScoreTable( storage::Table< double> &SCORE_TABLE) const;

      //! @brief filter the passed score table according to criteria passed on command line
      //! @param SCORE_TABLE the table to filter
      void FilterScoreTable( storage::Table< double> &SCORE_TABLE) const;

      //! @brief extract and print information about a particular row and column as passed on command line
      //! @param SCORE_TABLE the table to extract from
      void PrintColumnRowInformation( storage::Table< double> &SCORE_TABLE) const;

      static const ApplicationType FoldAnalysis_Instance; //!< application instance

    }; // class FoldAnalysis

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    FoldAnalysis::FoldAnalysis() :
      m_TableFromFileFlag
      (
        new command::FlagStatic
        (
          "table_from_file",
          "file to read score table from",
          command::Parameter( "file", "score table file", "")
        )
      ),
      m_TableFromPDBListFlag
      (
        new command::FlagStatic
        (
          "table_from_pdb_list",
          "file containing list of PDB files to create a score table",
          command::Parameter( "file", "pdb list file containing list of PDB files, 1 per line", "")
        )
      ),
      m_TableFromPrefixListFlag
      (
        new command::FlagDynamic
        (
          "table_from_prefix",
          "all prefixes that contain all strings and not contain any !strings will be considered",
          command::Parameter( "string_list", "list of strings contained and !strings not contained in desired prefixes")
        )
      ),
      m_TableFromStagesFlag
      (
        new command::FlagDynamic
        (
          "table_from_stages",
          "specify stages to use (not together with table from file and table from pdb list)",
          command::Parameter( "stage", "stage to use")
        )
      ),
      m_SortFlag
      (
        new command::FlagDynamic
        (
          "sort",
          "sort score table with by column name(s), if multiple given, the sum will be considered",
          command::Parameter( "sort_column", "name of column to be sorted")
        )
      ),
      m_SortCriteriaFlag
      (
        new command::FlagStatic
        (
          "sort_criteria",
          "sorting criteria",
          command::Parameter
          (
            "sort_criteria",
            "sorting criteria",
            command::ParameterCheckEnumerate< math::Comparisons< double> >(),
            math::Comparisons< double>::GetEnums().e_Less.GetName()
          )
        )
      ),
      m_FilterFlag
      (
        new command::FlagStatic
        (
          "filter",
          "filter score table to remove incomplete models",
          command::Parameter( "range", "minimum completeness is average minus range", "0.2")
        )
      ),
      m_OutputTableFlag
      (
        new command::FlagStatic
        (
          "output_table",
          "filename for table output",
          command::Parameter( "output", "filename for table output", "score.tbl")
        )
      ),
      m_OutputPdbListFlag
      (
        new command::FlagStatic
        (
          "output_pdb_list",
          "filename for (filtered) pdb list output",
          command::Parameter( "output", "filename for list output", "pdbs-analysed.ls")
        )
      ),
      m_OutputValueFlag( new command::FlagStatic( "output_value", "value (column, row) to evaluate")),
      m_ColumnParam
      (
        new command::Parameter
        (
          "column",
          "name of column to be evaluated",
          quality::GetMeasures().e_RMSD.GetName()
        )
      ),
      m_RowParam
      (
        new command::Parameter
        (
          "row",
          "row to be evaluated (i.e. 10 will print the value at the 10th row and statistics up to the 10th row)",
          "1"
        )
      ),
      m_FilenameParam
      (
        new command::Parameter
        (
          "filename",
          "name of file that will contain output data",
          "fold_analysis.out"
        )
      ),
      m_FlagUseUnweightedScores
      (
        new command::FlagStatic
        (
          "use_unweighted_scores",
          "If this flag is set, the unweighted scores will be gotten from the pdbs"
        )
      )
    {
      // add parameters to flags
      m_OutputValueFlag->PushBack( m_ColumnParam);
      m_OutputValueFlag->PushBack( m_RowParam);
      m_OutputValueFlag->PushBack( m_FilenameParam);
    }

    //! @brief Clone function
    //! @return pointer to new FoldAnalysis
    FoldAnalysis *FoldAnalysis::Clone() const
    {
      return new FoldAnalysis( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &FoldAnalysis::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the Command object
    util::ShPtr< command::Command> FoldAnalysis::InitializeCommand() const
    {
      // create command
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add member flags
      sp_cmd->AddFlag( m_TableFromFileFlag);
      sp_cmd->AddFlag( m_TableFromPDBListFlag);
      sp_cmd->AddFlag( m_TableFromPrefixListFlag);
      sp_cmd->AddFlag( m_TableFromStagesFlag);

      sp_cmd->AddFlag( m_SortFlag);
      sp_cmd->AddFlag( m_SortCriteriaFlag);
      sp_cmd->AddFlag( m_FilterFlag);

      sp_cmd->AddFlag( m_OutputTableFlag);
      sp_cmd->AddFlag( m_OutputPdbListFlag);
      sp_cmd->AddFlag( m_OutputValueFlag);

      // add storage flag
      sp_cmd->AddFlag( assemble::ProteinStorageFile::GetDefaultStorageFlag());

      // add flag to specify getting unweighted scores
      sp_cmd->AddFlag( m_FlagUseUnweightedScores);

      // add default bcl flags
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled command
      return sp_cmd;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief the Main function
    //! @return error code - 0 for success
    int FoldAnalysis::Main() const
    {
      // get the score table from whatever the user specified
      storage::Table< double> table( GetScoreTable());

      // sort the score table
      SortScoreTable( table);

      // filter the score table
      FilterScoreTable( table);

      // if the table is to be written, write it
      if( m_OutputTableFlag->GetFlag())
      {
        io::OFStream write;
        io::File::MustOpenOFStream( write, m_OutputTableFlag->GetFirstParameter()->GetValue());
        table.WriteFormatted( write, pdb::PrinterScore::GetTableColumnFormat());
        io::File::CloseClearFStream( write);
      }

      // if the (filtered, sorted) pdb list is to be written, write it
      if( m_OutputPdbListFlag->GetFlag())
      {
        io::OFStream write;
        io::File::MustOpenOFStream( write, m_OutputPdbListFlag->GetFirstParameter()->GetValue());
        storage::List< std::string> pdbs( table.GetRowNames());
        for( storage::List< std::string>::const_iterator itr( pdbs.Begin()), itr_end( pdbs.End()); itr != itr_end; ++itr)
        {
          write << *itr << "\n";
        }
        io::File::CloseClearFStream( write);
      }

      // print information about a particular column and row if requested
      PrintColumnRowInformation( table);

      // end
      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &FoldAnalysis::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &FoldAnalysis::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief create a score table
    //! @return score table
    storage::Table< double> FoldAnalysis::GetScoreTable() const
    {
      // make sure the user set reasonable arguments
      const size_t how_many_flags
      (
        m_TableFromFileFlag->GetFlag() + m_TableFromPDBListFlag->GetFlag() + m_TableFromPrefixListFlag->GetFlag()
      );
      BCL_Assert
      (
        how_many_flags < 2,
        "At most one of score table, pdb list, prefix and prefix strings flags must be supplied!"
      );

      // give a warning as stages flag is ignored when table_from_file or table_from_pdb_list flags are set
      if( ( m_TableFromFileFlag->GetFlag() || m_TableFromPDBListFlag->GetFlag()) && m_TableFromStagesFlag->GetFlag())
      {
        BCL_MessageStd
        (
          "Stages flag is ignored when score table is read from a file or created from a pdb list!"
        );
      }

      storage::Table< double> table; // initialize table

      if( m_TableFromFileFlag->GetFlag()) // if the user passed a filename for an existing score table, read it
      {
        const std::string &score_filename( m_TableFromFileFlag->GetFirstParameter()->GetValue());
        BCL_MessageStd( "reading score table from file " + util::Format()( score_filename));

        io::IFStream read;
        io::File::MustOpenIFStream( read, score_filename);
        table.ReadFormatted( read);
        io::File::CloseClearFStream( read);
      }
      else if( m_TableFromPDBListFlag->GetFlag()) // otherwise create a new table, e.g. from a pdb list
      {
        const std::string &list_filename( m_TableFromPDBListFlag->GetFirstParameter()->GetValue());
        BCL_MessageStd( "reading pdb from list file " + util::Format()( list_filename));

        // iterate over the pdbs listed in the passed file
        io::IFStream read;
        io::File::MustOpenIFStream( read, list_filename);
        while( !read.eof())
        {
          // read the filename
          std::string pdb_filename;
          read >> pdb_filename;
          util::TrimString( pdb_filename);

          // if it is not valid, ignore
          if( pdb_filename.empty())
          {
            continue;
          }

          // open the file
          io::IFStream read_pdb;
          io::File::MustOpenIFStream( read_pdb, pdb_filename);

          // construct a handler
          const pdb::Handler handler( read_pdb);
          io::File::CloseClearFStream( read_pdb);

          // add the row from the remark lines
          AddTableRow( handler.GetLines( pdb::GetLineTypes().REMARK), pdb_filename, table, m_FlagUseUnweightedScores->GetFlag());
        }
        io::File::CloseClearFStream( read);
      }
      else // or create a new table from all files or db-entries that match the prefix or prefix-contains flags
      // or all that are found if none of the two flags are given
      {
        BCL_MessageStd( "reading matching pdb content from storage");

        // create the storage object
        const util::ShPtr< assemble::ProteinStorageFile>
          sp_storage( assemble::ProteinStorageFile::GetDefaultStorage());

        // initialize prefixes
        storage::Set< std::string> prefixes;

        // get all prefixes stored
        const storage::Set< std::string> all_prefixes( sp_storage->GetAllSources());

        if( m_TableFromPrefixListFlag->GetFlag()) // if strings are passed
        {
          // add parameter values
          const storage::Set< std::string> prefix_strings( m_TableFromPrefixListFlag->GetObjectSet< std::string>());

          // iterate over all the prefixes in the storage
          for
          (
            storage::Set< std::string>::const_iterator prefix_itr( all_prefixes.Begin()),
              prefix_itr_end( all_prefixes.End());
            prefix_itr != prefix_itr_end; ++prefix_itr
          )
          {
            bool matches_strings( true); // initialize bool for storing whether all strings are contained in the prefix

            // iterate over the prefix strings
            for
            (
              storage::Set< std::string>::const_iterator string_itr( prefix_strings.Begin()),
                string_itr_end( prefix_strings.End());
              string_itr != string_itr_end; ++string_itr
            )
            {
              if( string_itr->length() > 0 && ( *string_itr)[ 0] == '!') // if the prefix starts with '!'
              // it may be necessary to escape the '!' on the command line
              {
                if( prefix_itr->find( string_itr->substr( 1)) != std::string::npos) // and prefix is found
                {
                  matches_strings = false; // set bool to false
                  break;
                }
              }
              else // prefix does NOT start with '!'
              {
                if( prefix_itr->find( *string_itr) == std::string::npos) // and prefix does NOT contain the string
                {
                  matches_strings = false; // set bool to false
                  break;
                }
              }
            } // for

            // add the prefix
            if( matches_strings)
            {
              prefixes.Insert( *prefix_itr);
            }
          } // for
        } // if prefix-flag
        else // no prefix flag, so use all prefixes
        {
          prefixes = all_prefixes;
        }

        BCL_Assert( prefixes.GetSize() != 0, "No protein models found in the storage match the given prefix(es)!");

        // iterate over the prefixes
        for
        (
          storage::Set< std::string>::const_iterator prefix_itr( prefixes.Begin()), prefix_itr_end( prefixes.End());
          prefix_itr != prefix_itr_end; ++prefix_itr
        )
        {
          // get all keys using the given prefix
          const storage::Vector< std::string> all_keys( sp_storage->GetAllKeys( *prefix_itr));

          // iterate over the keys
          for
          (
            storage::Vector< std::string>::const_iterator key_itr( all_keys.Begin()), key_itr_end( all_keys.End());
            key_itr != key_itr_end; ++key_itr
          )
          {
            // check if this key has a stage
            const bool has_stage( key_itr->size() >= 6 && key_itr->at( key_itr->length() - 6) == '_');

            // if specific stages are requested and it has a stage
            if( m_TableFromStagesFlag->GetFlag() && has_stage)
            {
              // iterate over the passed stages
              for
              (
                util::ShPtrVector< command::ParameterInterface>::const_iterator
                  stage_itr( m_TableFromStagesFlag->GetParameterList().Begin()),
                  stage_itr_end( m_TableFromStagesFlag->GetParameterList().End());
                stage_itr != stage_itr_end; ++stage_itr
              )
              {
                const size_t stage_number
                (
                  util::ConvertStringToNumericalValue< size_t>( key_itr->substr( key_itr->length() - 5, 5))
                );
                if( stage_number == ( *stage_itr)->GetNumericalValue< size_t>())
                {
                  AddTableRow( sp_storage->RetrieveRemarkLines( *prefix_itr, *key_itr), *prefix_itr + *key_itr, table, m_FlagUseUnweightedScores->GetFlag());
                }
              }
            }
            else if( !m_TableFromStagesFlag->GetFlag() && !has_stage)
            {
              AddTableRow( sp_storage->RetrieveRemarkLines( *prefix_itr, *key_itr), *prefix_itr + *key_itr, table, m_FlagUseUnweightedScores->GetFlag());
            }
          } // for keys
        } // for prefixes
      } // else

      // make sure table exists i.e. is not empty
      BCL_Assert( table.GetNumberRows() != 0, "Score table is empty!");

      // end
      return table;
    }

    //! @param LINES pdb lines to analyze
    //! @param ROW_NAME row name
    //! @param TABLE table to write to
    //! @param USE_UNWEIGHTED_SCORES unweighted scores will be gotten
    void FoldAnalysis::AddTableRow
    (
      const util::ShPtrList< pdb::Line> &LINES,
      const std::string &ROW_NAME,
      storage::Table< double> &TABLE,
      const bool USE_UNWEIGHTED_SCORES
    )
    {
      BCL_MessageStd( "adding pdb file " + util::Format()( ROW_NAME));

      // initialize vectors
      storage::Vector< std::string> column_names;
      storage::Vector< double> values;

      static const size_t weighted_score_column_index( 3);
      static const size_t unweighted_score_column_index( 2);
      const size_t score_column_index
      (
        USE_UNWEIGHTED_SCORES ? unweighted_score_column_index : weighted_score_column_index
      );

      // iterate over the remark lines
      for
      (
        util::ShPtrList< pdb::Line>::const_iterator line_itr( LINES.Begin()),
          line_itr_end( LINES.End());
        line_itr != line_itr_end; ++line_itr
      )
      {
        // if this is a score table line
        if
        (
          util::ConvertStringToNumericalValue< size_t>( ( *line_itr)->GetString( pdb::GetEntryTypes().REMARK_Number)) ==
            pdb::PrinterScore::s_RemarkNumber
        )
        {
          // get the remark string and split it
          const storage::Vector< std::string> split_line
          (
            util::SplitString( ( *line_itr)->GetString( pdb::GetEntryTypes().REMARK_String))
          );

          // if this is a valid line
          if( split_line.GetSize() == 4 && util::IsNumerical( split_line( score_column_index)))
          {
            // pushback the name and value
            column_names.PushBack( split_line( 0));
            values.PushBack( util::ConvertStringToNumericalValue< double>( split_line( score_column_index)));
          }
        }
      }

      // if a header has not been stored
      if( TABLE.GetSize() == 0)
      {
        TABLE = storage::Table< double>( storage::TableHeader( column_names));
      }

      // insert the row into the table
      TABLE.InsertRow( ROW_NAME, values);
    }

    //! @brief sort the passed score table according to criteria passed on command line
    //! @param SCORE_TABLE the table to sort
    void FoldAnalysis::SortScoreTable( storage::Table< double> &SCORE_TABLE) const
    {
      if( !m_SortFlag->GetFlag()) // if no sorting is requested, do not change the score table at all
      {
        return;
      }

      // if sorting is requested, but no column is passed, just sort by row name (protein model name) as default
      if( m_SortFlag->GetParameterList().GetSize() == 0)
      {
        SCORE_TABLE.SortByRowName( std::less_equal< std::string>());
        return;
      }

      if( m_SortFlag->GetParameterList().GetSize() == 1) // if only one parameter is passed
      {
        // sort the table by the given column name
        SCORE_TABLE.SortByColumn
        (
          m_SortFlag->GetFirstParameter()->GetValue(),
          **math::Comparisons< double>::Comparison( m_SortCriteriaFlag->GetFirstParameter()->GetValue())
        );
        return;
      }

      // multiple columns need to be summed
      storage::TableHeader new_header( SCORE_TABLE.GetHeader()); // copy the header

      // add a new column
      const std::string sum_column_name( "sum_for_sort");
      new_header.PushBack( sum_column_name);

      // create a new table
      storage::Table< double> new_table( new_header);

      // initialize vector of indices
      storage::Vector< size_t> indices;

      // iterate over the passed column names
      for
      (
        util::ShPtrVector< command::ParameterInterface>::const_iterator
          param_itr( m_SortFlag->GetParameterList().Begin()),
          param_itr_end( m_SortFlag->GetParameterList().End());
        param_itr != param_itr_end; ++param_itr
      )
      {
        // check if the table has this column
        if( new_header.HasColumn( ( *param_itr)->GetValue()))
        {
          indices.PushBack( new_header[ ( *param_itr)->GetValue()]);
        }
        else
        {
          BCL_MessageStd
          (
            "No column found for " + ( *param_itr)->GetValue() + ", so skipping it for sorting"
          );
        }
      }

      // iterate over the original table
      for
      (
        storage::Table< double>::const_iterator table_itr( SCORE_TABLE.Begin()), table_itr_end( SCORE_TABLE.End());
        table_itr != table_itr_end; ++table_itr
      )
      {
        // get the data from the row
        storage::Vector< double> row_data( table_itr->Second().GetData());

        // initialize sum
        double sum( 0);

        // iterate over the indices
        for
        (
          storage::Vector< size_t>::const_iterator index_itr( indices.Begin()), index_itr_end( indices.End());
          index_itr != index_itr_end; ++index_itr
        )
        {
          sum += row_data( *index_itr);
        }

        // add the sum to the vector
        row_data.PushBack( sum);

        // insert a new row
        new_table.InsertRow( table_itr->First(), row_data);
      }

      // set the table
      SCORE_TABLE = new_table;

      // now sort the table by the new sum column
      SCORE_TABLE.SortByColumn
      (
        sum_column_name,
        **math::Comparisons< double>::Comparison( m_SortCriteriaFlag->GetFirstParameter()->GetValue())
      );
    }

    //! @brief filter the passed score table according to criteria passed on command line
    //! @param SCORE_TABLE the table to filter
    void FoldAnalysis::FilterScoreTable( storage::Table< double> &SCORE_TABLE) const
    {
      if( !m_FilterFlag->GetFlag()) // do not filter if not requested
      {
        return;
      }

      // try to find a completeness column
      static const std::string completeness_column_header( "completeness");
      static const std::string completeness_estimate_column_header( "completeness_estimate");

      // if completeness header isn't found, return
      if( !SCORE_TABLE.GetHeader().HasColumn( completeness_column_header) && !SCORE_TABLE.GetHeader().HasColumn( completeness_estimate_column_header))
      {
        BCL_MessageStd( "Cannot filter, score table is missing a completeness column!");
        return;
      }

      // determine which column header to use
      const std::string column_header
      (
        SCORE_TABLE.GetHeader().HasColumn( completeness_column_header) ?
          completeness_column_header : completeness_estimate_column_header
      );

      // filter
      BCL_MessageStd( "Filter for complete models.");

      math::RunningAverageSD< double> stats_mean; // to calculate the mean
      math::RunningMinMax< double> stats_min_max; // to calculate the min and max

      // collect data from table
      for
      (
        storage::Table< double>::const_iterator row_itr( SCORE_TABLE.Begin()), row_itr_end( SCORE_TABLE.End());
        row_itr != row_itr_end; ++row_itr
      )
      {
        const double &current_completeness( row_itr->Second()[ column_header]);
        stats_mean += current_completeness;
        stats_min_max += current_completeness;
      }

      BCL_MessageStd( "Maximal completeness: " + util::Format()( stats_min_max.GetMax()));
      BCL_MessageStd( "Minimal completeness: " + util::Format()( stats_min_max.GetMin()));
      BCL_MessageStd( "Average completeness: " + util::Format()( stats_mean.GetAverage()));

      double min_completeness_threshold
      (
        stats_mean.GetAverage() - m_FilterFlag->GetFirstParameter()->GetNumericalValue< double>()
      );

      BCL_MessageStd
      (
        "Allowed minimal completeness for filtering: " + util::Format()( min_completeness_threshold)
      );

      storage::Table< double> new_table( SCORE_TABLE.GetHeader()); // create a new table

      // iterate again over the original table and save models passing criteria into new table
      for
      (
        storage::Table< double>::const_iterator row_itr( SCORE_TABLE.Begin()), row_itr_end( SCORE_TABLE.End());
        row_itr != row_itr_end; ++row_itr
      )
      {
        const double &current_completeness( row_itr->Second()[ column_header]);

        if( current_completeness >= min_completeness_threshold)
        {
          new_table.InsertRow( row_itr->First(), row_itr->Second().GetData());
        }
      }

      // print how many models passed the criteria
      BCL_MessageStd
      (
        util::Format()( new_table.GetSize()) + "/" + util::Format()( SCORE_TABLE.GetSize())
          + " models passed the filter criteria."
      );

      SCORE_TABLE = new_table; // set the table
    }

    //! @brief extract and print information about a particular row and column as passed on command line
    //! @param SCORE_TABLE the table to extract from
    void FoldAnalysis::PrintColumnRowInformation( storage::Table< double> &SCORE_TABLE) const
    {
      if( !m_OutputValueFlag->GetFlag()) // do not print anything if flag was not set
      {
        return;
      }

      // print the requested info
      const size_t index_end( m_RowParam->GetNumericalValue< size_t>());
      size_t index( 0);
      math::RunningAverageSD< double> stats_mean_sd;
      math::RunningMinMax< double> stats_min_max;
      double final_value( 0);
      for
      (
        storage::Table< double>::const_iterator row_itr( SCORE_TABLE.Begin()), row_itr_end( SCORE_TABLE.End());
        row_itr != row_itr_end && index != index_end; ++row_itr, ++index
      )
      {
        stats_mean_sd += row_itr->Second()[ m_ColumnParam->GetValue()];
        stats_min_max += row_itr->Second()[ m_ColumnParam->GetValue()];
        if( index + 1 == index_end)
        {
          final_value = row_itr->Second()[ m_ColumnParam->GetValue()];
        }
      }
      BCL_MessageCrt
      (
        m_ColumnParam->GetValue() + " at row number " + util::Format()( index_end) + " is " +
          util::Format()( final_value) + " and mean of 1.." + util::Format()( index_end) + " is " +
          util::Format()( stats_mean_sd.GetAverage()) + " +/- " + util::Format()( stats_mean_sd.GetStandardDeviation())
      );

      io::OFStream write;
      io::File::MustOpenOFStream( write, m_FilenameParam->GetValue());
      write << m_ColumnParam->GetValue() << "\tPosition:\t" << index_end << "\tMax:\t" << SCORE_TABLE.GetSize()
            << "\tValue:\t" << final_value << "\tMean:\t" << stats_mean_sd.GetAverage()
            << "\tS.D.:\t" << stats_mean_sd.GetStandardDeviation() << "\tMin:\t" << stats_min_max.GetMin()
            << "\tMax:\t" << stats_min_max.GetMax() << '\n';
      io::File::CloseClearFStream( write);
    }

    const ApplicationType FoldAnalysis::FoldAnalysis_Instance( GetAppGroups().AddAppToGroup( new FoldAnalysis(), GetAppGroups().e_Protein));

  } // namespace app
} // namespace bcl
