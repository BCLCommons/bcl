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
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_static.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_gnuplot.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_statistics.h"
#include "storage/bcl_storage_table.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class DistributionPlotFromTables
    //! @brief TODO: add a brief comment
    //! @details TODO: add an detailed description to this class
    //!
    //! @see @link example_app_DistributionPlotFromTables.cpp @endlink
    //! @author alexanns, weinerbe
    //! @date Mar 7, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API DistributionPlotFromTables :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      util::ShPtr< command::FlagStatic> m_FlagWriteHistogramsToFile;
      util::ShPtr< command::ParameterInterface> m_ParamWriteHistogramsToFile;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      DistributionPlotFromTables();

      //! @brief Clone function
      //! @return pointer to new DistributionPlotFromTables
      DistributionPlotFromTables *Clone() const;

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

    ///////////////
    // operators //
    ///////////////

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

      static const ApplicationType DistributionPlotFromTables_Instance;

    }; // class DistributionPlotFromTables

    //! @brief default constructor
    DistributionPlotFromTables::DistributionPlotFromTables() :
      m_FlagWriteHistogramsToFile
      (
        new command::FlagStatic
        (
          "histogram_output_filename",
          "The name of the file which will contain the histograms making up the plot."
        )
      ),
      m_ParamWriteHistogramsToFile
      (
        new command::Parameter
        (
          "histogram_output_filename",
          "string which is the name of the file to contain the histograms",
          "histograms.txt"
        )
      )
    {
      m_FlagWriteHistogramsToFile->PushBack( m_ParamWriteHistogramsToFile);
    }

    //! @brief Clone function
    //! @return pointer to new FoldAnalysis
      DistributionPlotFromTables *DistributionPlotFromTables::Clone() const
    {
      return new DistributionPlotFromTables( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &DistributionPlotFromTables::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the Command object
    util::ShPtr< command::Command> DistributionPlotFromTables::InitializeCommand() const
    {
      // create command
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add member flags
      sp_cmd->AddFlag( math::Gnuplot::GetFlagTableInputFilenames());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagTableColumnsX());
      sp_cmd->AddFlag( math::Histogram::GetFlagMin());
      sp_cmd->AddFlag( math::Histogram::GetFlagBinSize());
      sp_cmd->AddFlag( math::Histogram::GetFlagNumberOfBins());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagGnuplotOutputFilename());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagXPixels());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagYPixels());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagTitle());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagSetKey());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagXLabel());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagYMax());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagSeriesNames());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagTicIntervalX());
      sp_cmd->AddFlag( m_FlagWriteHistogramsToFile);

      // add default bcl flags
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled command
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int DistributionPlotFromTables::Main() const
    {
      // get all the table filenames
      storage::Vector< std::string> table_filenames( math::Gnuplot::GetFlagTableInputFilenames()->GetStringList());

      // to hold all the tables
      storage::Vector< storage::Table< double> > tables;

      // iterate over the table filenames and create the tables
      for
      (
        storage::Vector< std::string>::const_iterator
          filename_itr( table_filenames.Begin()), filename_itr_end( table_filenames.End());
        filename_itr != filename_itr_end;
        ++filename_itr
      )
      {
        // tmp table
        storage::Table< double> table;

        // read it in
        io::IFStream read;
        io::File::MustOpenIFStream( read, *filename_itr);
        table.ReadFormatted( read);
        io::File::CloseClearFStream( read);

        // add the table to the vector of tables
        tables.PushBack( table);
      }

      // get the list of desired columns
      const storage::Vector< std::string> table_columns( math::Gnuplot::GetFlagTableColumnsX()->GetStringList());

      // tables with desired columns
      storage::Vector< storage::Table< double> > tables_columns;

      // data for individual lines/distributions in the plot
      storage::Vector< math::Histogram> plot_series
      (
        tables.GetSize(),
        math::Histogram
        (
          math::Histogram::GetFlagMin()->GetFirstParameter()->GetNumericalValue< double>(),
          math::Histogram::GetFlagBinSize()->GetFirstParameter()->GetNumericalValue< double>(),
          math::Histogram::GetFlagNumberOfBins()->GetFirstParameter()->GetNumericalValue< double>()
        )
      );

      // get the columns from the tables
      for
      (
        storage::Vector< storage::Table< double> >::const_iterator
          table_itr( tables.Begin()), table_itr_end( tables.End());
          table_itr != table_itr_end;
        ++table_itr
      )
      {
        // current table
        const storage::Table< double> &table( *table_itr);

        // get current histogram corresponding to this table
        math::Histogram &current_histogram( plot_series( table_itr - tables.Begin()));

        // create new subtable out of columns
        const storage::Table< double> new_table( table.ExtractSubTableFromColumns( storage::TableHeader( table_columns)));

        // iterate over the sub table to sum the columns
        for
        (
          storage::Table< double>::const_iterator
            new_table_itr( new_table.Begin()), new_table_itr_end( new_table.End());
          new_table_itr != new_table_itr_end;
          ++new_table_itr
        )
        {
          // add column sum to the current histogram
          current_histogram.PushBack
          (
            math::Statistics::Sum( new_table_itr->Second().GetData().Begin(), new_table_itr->Second().GetData().End())
          );
        }

        current_histogram.Normalize();
      }

      // get the gnuplot output filename
      const std::string gnuplot_filename( math::Gnuplot::GetFlagGnuplotOutputFilename()->GetFirstParameter()->GetValue());
      // open file
      io::OFStream write;
      io::File::MustOpenOFStream( write, gnuplot_filename);

      write << "# BCL generated linear plot from histogram\n";

      write << "set terminal png transparent enhanced size ";
      write << math::Gnuplot::GetFlagXPixels()->GetFirstParameter()->GetNumericalValue< size_t>();
      write << ",";
      write << math::Gnuplot::GetFlagYPixels()->GetFirstParameter()->GetNumericalValue< size_t>();
      write << "\n";

      write << "set output \"" << gnuplot_filename << ".png\"\n";
      write << "set encoding iso\n";
      write << "set view map\n";
      if( math::Gnuplot::GetFlagTitle()->GetFlag())
      {
        write << "set title \"" << math::Gnuplot::GetFlagTitle()->GetFirstParameter()->GetValue() << "\"\n";
      }
      if( !math::Gnuplot::GetFlagSetKey()->GetFlag())
      {
        write << "unset key\n\n";
      }

      write << "set xlabel \"" << math::Gnuplot::GetFlagXLabel()->GetFirstParameter()->GetValue() << "\"\n";
      write << "set xrange [" << plot_series.FirstElement().GetBoundaries().First() << ":"
            << plot_series.FirstElement().GetBoundaries().Second() << "]\n";
      write << "set ylabel \"Fraction\"\n";
      if( math::Gnuplot::GetFlagTicIntervalX()->GetFlag())
      {
        write << "set xtics "
              << math::Gnuplot::GetFlagTicIntervalX()->GetFirstParameter()->GetNumericalValue< double>()
              << '\n';
      }
      // true if the max value of y - axis is not user specified
      if( !math::Gnuplot::GetFlagYMax()->GetFlag())
      {
        write << "set autoscale y\n";
      }
      else //< user specified max y value
      {
        write << "set yrange [0 : "
              << math::Gnuplot::GetFlagYMax()->GetFirstParameter()->GetNumericalValue< double>()<< "]\n";
      }

      write << "plot ";

//      // get the numbers specifying what the plot lines should look like
//      storage::Vector< int> line_types( m_FlagLineTypes->GetNumericalList< int>());

      // get all the series names
      storage::Vector< std::string> series_names( math::Gnuplot::GetFlagSeriesNames()->GetStringList());

      // iterate over the data histograms to write the plot command
      for
      (
        storage::Vector< std::string>::const_iterator
          series_names_itr( series_names.Begin()), series_names_itr_end( series_names.End());
        series_names_itr != series_names_itr_end;
        ++series_names_itr
      )
      {
        // dont' print comma if at first series
        if( series_names_itr != series_names.Begin())
        {
          write << ", ";
        }

        // print command to plot current series
        write << " '-' with linespoints lw 2 title \"" << *series_names_itr << "\"";
      }

      write << '\n';

      io::OFStream histogram_file;
      if( m_FlagWriteHistogramsToFile->GetFlag())
      {
        io::File::MustOpenOFStream( histogram_file, m_ParamWriteHistogramsToFile->GetValue());
      }

      // iterate through the histograms to write out the data from the histograms
      for
      (
        storage::Vector< math::Histogram>::const_iterator
          histogram_itr( plot_series.Begin()), histogram_itr_end( plot_series.End());
        histogram_itr != histogram_itr_end;
        ++histogram_itr
      )
      {
        math::Histogram::WriteGnuPlotLinePlotFormatted( write, histogram_itr->GetBinning(), histogram_itr->GetHistogram());

        if( m_FlagWriteHistogramsToFile->GetFlag())
        {
          histogram_file << m_ParamWriteHistogramsToFile->GetValue()
                         << " " << series_names( histogram_itr - plot_series.Begin()) << " ";
          histogram_itr->WriteHorizontally( histogram_file, util::Format().W( 8).FFP( 5), util::Format().W( 8).FFP( 5));
        }
      }

      // end
      return 0;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &DistributionPlotFromTables::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &DistributionPlotFromTables::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

    const ApplicationType DistributionPlotFromTables::DistributionPlotFromTables_Instance( GetAppGroups().AddAppToGroup( new DistributionPlotFromTables(), GetAppGroups().e_Bcl));

  } // namespace app
} // namespace bcl
