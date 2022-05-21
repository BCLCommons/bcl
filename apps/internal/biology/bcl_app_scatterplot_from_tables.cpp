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
#include "bcl_app_scatterplot_from_tables.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "io/bcl_io_file.h"
#include "math/bcl_math_gnuplot.h"
#include "math/bcl_math_statistics.h"
#include "storage/bcl_storage_table.h"
#include "storage/bcl_storage_vector_nd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace app
  {

    const ApplicationType ScatterplotFromTables::ScatterplotFromTables_Instance( GetAppGroups().AddAppToGroup( new ScatterplotFromTables(), GetAppGroups().e_Bcl));

    //! @brief default constructor
    ScatterplotFromTables::ScatterplotFromTables()
    {
    }

    //! @brief Clone function
    //! @return pointer to new FoldAnalysis
    ScatterplotFromTables *ScatterplotFromTables::Clone() const
    {
      return new ScatterplotFromTables( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ScatterplotFromTables::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the Command object
    util::ShPtr< command::Command> ScatterplotFromTables::InitializeCommand() const
    {
      // create command
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      // add member flags
      sp_cmd->AddFlag( math::Gnuplot::GetFlagTableInputFilenames());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagTableColumnsX());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagTableColumnsY());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagGnuplotOutputFilename());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagXPixels());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagYPixels());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagSetKey());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagXLabel());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagYLabel());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagXMin());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagXMax());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagYMin());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagYMax());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagTicIntervalX());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagTicIntervalY());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagSeriesNames());
      sp_cmd->AddFlag( math::Gnuplot::GetFlagTitle());

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
    int ScatterplotFromTables::Main() const
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
      const storage::Vector< std::string> table_columns_x( math::Gnuplot::GetFlagTableColumnsX()->GetStringList());
      const storage::Vector< std::string> table_columns_y( math::Gnuplot::GetFlagTableColumnsY()->GetStringList());

      // data for individual x y values for every table
      storage::Vector< storage::Vector< storage::VectorND< 2, double> > > plot_series( tables.GetSize());

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
        storage::Vector< storage::VectorND< 2, double> > &current_x_y( plot_series( table_itr - tables.Begin()));

        // create new subtable out of columns
        const storage::Table< double> new_table_x( table.ExtractSubTableFromColumns( storage::TableHeader( table_columns_x)));
        const storage::Table< double> new_table_y( table.ExtractSubTableFromColumns( storage::TableHeader( table_columns_y)));

        // iterate over the sub table to sum the columns
        for
        (
          storage::Table< double>::const_iterator
            new_table_x_itr( new_table_x.Begin()), new_table_itr_x_end( new_table_x.End()),
            new_table_y_itr( new_table_y.Begin()), new_table_itr_y_end( new_table_y.End());
            new_table_x_itr != new_table_itr_x_end && new_table_y_itr != new_table_itr_y_end;
          ++new_table_x_itr, ++new_table_y_itr
        )
        {
          storage::VectorND< 2, double> x_y
          (
            math::Statistics::Sum( new_table_x_itr->Second().GetData().Begin(), new_table_x_itr->Second().GetData().End()),
            math::Statistics::Sum( new_table_y_itr->Second().GetData().Begin(), new_table_y_itr->Second().GetData().End())
          );

          current_x_y.PushBack( x_y);
        }
      }

      // get the gnuplot output filename
      const std::string gnuplot_filename( math::Gnuplot::GetFlagGnuplotOutputFilename()->GetFirstParameter()->GetValue());
      // open file
      io::OFStream write;
      io::File::MustOpenOFStream( write, gnuplot_filename);

      write << "# BCL generated linear plot from histogram\n";

      write << "set terminal png transparent enhanced size "
            << math::Gnuplot::GetFlagXPixels()->GetFirstParameter()->GetNumericalValue< size_t>() << ","
            << math::Gnuplot::GetFlagYPixels()->GetFirstParameter()->GetNumericalValue< size_t>() << "\n";
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
      write << "set ylabel \"" << math::Gnuplot::GetFlagYLabel()->GetFirstParameter()->GetValue() << "\"\n";

      if( math::Gnuplot::GetFlagXMin()->GetFlag() && math::Gnuplot::GetFlagXMax()->GetFlag())
      {
        write << "set xrange [" << math::Gnuplot::GetFlagXMin()->GetFirstParameter()->GetNumericalValue< double>()
              << ":"
              << math::Gnuplot::GetFlagXMax()->GetFirstParameter()->GetNumericalValue< double>() << "]\n";
      }
      if( math::Gnuplot::GetFlagYMin()->GetFlag() && math::Gnuplot::GetFlagYMax()->GetFlag())
      {
        write << "set yrange [" << math::Gnuplot::GetFlagYMin()->GetFirstParameter()->GetNumericalValue< double>()
              << ":"
              << math::Gnuplot::GetFlagYMax()->GetFirstParameter()->GetNumericalValue< double>() << "]\n";
      }

      if( math::Gnuplot::GetFlagTicIntervalX()->GetFlag())
      {
        write << "set xtics "
              << math::Gnuplot::GetFlagTicIntervalX()->GetFirstParameter()->GetNumericalValue< double>()
              << '\n';
      }

      if( math::Gnuplot::GetFlagTicIntervalY()->GetFlag())
      {
        write << "set ytics "
              << math::Gnuplot::GetFlagTicIntervalY()->GetFirstParameter()->GetNumericalValue< double>()
              << '\n';
      }

      write << "set style data points\n";

      write << "plot ";

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
        write << " '-' title \"" << *series_names_itr << "\"";
      }

      write << '\n';

      for
      (
        storage::Vector< storage::Vector< storage::VectorND< 2, double> > >::const_iterator
          series_itr( plot_series.Begin()), series_itr_end( plot_series.End());
          series_itr != series_itr_end;
          ++series_itr
      )
      {
        // get current series x y values
        const storage::Vector< storage::VectorND< 2, double> > &current_x_y_values( *series_itr);

        // iterate over current values to write them to the script
        for
        (
          storage::Vector< storage::VectorND< 2, double> >::const_iterator
            values_itr( current_x_y_values.Begin()), values_itr_end( current_x_y_values.End());
            values_itr != values_itr_end;
            ++values_itr
        )
        {
          write << values_itr->First() << " " << values_itr->Second() << '\n';
        }
        write << "e\n";
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
    std::istream &ScatterplotFromTables::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &ScatterplotFromTables::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  } // namespace app
} // namespace bcl
