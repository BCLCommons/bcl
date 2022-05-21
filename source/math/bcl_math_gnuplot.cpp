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
#include "math/bcl_math_gnuplot.h"

// includes from bcl - sorted alphabetically
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter.h"
#include "coord/bcl_coord_axes.h"
#include "linal/bcl_linal_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    //! @brief flag for setting number of pixels in x direction
    //! @return flag for setting number of pixels in x direction
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagXPixels()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_pixels_x",
          "The number of pixels in x direction.",
          command::Parameter( "number_pixels", "size_t which is the number of pixels", "900")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting number of pixels in y direction
    //! @return flag for setting number of pixels in y direction
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagYPixels()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_pixels_y",
          "The number of pixels in y direction.",
          command::Parameter( "number_pixels", "size_t which is the number of pixels", "600")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the title of the plot
    //! @return flag for setting the title of the plot
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagTitle()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_title",
          "The title of the gnuplot.",
          command::Parameter( "title", "string which is the title of the plot", "YourPlotTitleHere")
        )
      );

      return s_flag;
    }

    //! @brief flag for determining whether the key is shown or not in the plot
    //! @return flag for determining whether the key is shown or not in the plot
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagSetKey()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_set_key",
          "Determines whether the key is shown or not in the plot."
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the label on the x axis
    //! @return flag for setting the label on the x axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagXLabel()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_x_label",
          "The label of the x axis.",
          command::Parameter( "label", "string which is the label of the x axis", "x-axis")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the label on the y axis
    //! @return flag for setting the label on the y axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagYLabel()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_y_label",
          "The label of the y axis.",
          command::Parameter( "label", "string which is the label of the y axis", "y-axis")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the maximum value on the x axis
    //! @return flag for setting the maximum value on the x axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagXMax()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_max_x",
          "The maximum value shown on the x axis.",
          command::Parameter( "max_value", "double which is the maximum value shown on the x axis", "25.0")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the maximum value on the y axis
    //! @return flag for setting the maximum value on the y axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagYMax()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_max_y",
          "The maximum value shown on the y axis.",
          command::Parameter( "max_value", "double which is the maximum value shown on the y axis", "25.0")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the minimum value on the x axis
    //! @return flag for setting the minimum value on the x axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagXMin()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_min_x",
          "The minimum value shown on the x axis.",
          command::Parameter( "min_value", "double which is the minimum value shown on the x axis", "0.0")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting the minimum value on the y axis
    //! @return flag for setting the minimum value on the y axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagYMin()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_min_y",
          "The minimum value shown on the y axis.",
          command::Parameter( "min_value", "double which is the minimum value shown on the y axis", "0.0")
        )
      );

      return s_flag;
    }

    //! @brief flag for setting series names to be plotted
    //! @return flag for setting series names to be plotted
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagSeriesNames()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "gnuplot_series_names",
          "The series names to be plotted. Zero or more can be passed",
          command::Parameter( "series_name", "string which is a series names to be plotted", "series_a")
        )
      );

      return s_flag;
    }

    //! @brief flag for specifying filenames of tables to be included in distribution plot
    //! @return flag for specifying filenames of tables to be included in distribution plot
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagTableInputFilenames()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "gnuplot_input_table_filenames",
          "The filenames of tables to be included in distribution plot. Zero or more can be passed",
          command::Parameter( "table_filenames", "string which filename of a table to be plotted", "table_a")
        )
      );

      return s_flag;
    }

    //! @brief flag for specifying columns of the tables to be used for the plot by column name for x axis
    //! @return flag for specifying columns of the tables to be used for the plot by column name for x axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagTableColumnsX()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "gnuplot_table_columns_x",
          "The columns of the tables to be used for the plot by column name for x axis. Zero or more can be passed."
          " If more than one is passed they will be summed to create one value per row (datapoint).",
          command::Parameter( "column_name", "string which is the column names for plotting", "RMSD100")
        )
      );

      return s_flag;
    }

    //! @brief flag for specifying columns of the tables to be used for the plot by column name for y axis
    //! @return flag for specifying columns of the tables to be used for the plot by column name for y axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagTableColumnsY()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagDynamic
        (
          "gnuplot_table_columns_y",
          "The columns of the tables to be used for the plot by column name for y axis. Zero or more can be passed."
          " If more than one is passed they will be summed to create one value per row (datapoint).",
          command::Parameter( "column_name", "string which is the column names for plotting", "RMSD100")
        )
      );

      return s_flag;
    }

    //! @brief flag for specifying the name of the outputted gnuplot file
    //! @return flag for specifying the name of the outputted gnuplot file
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagGnuplotOutputFilename()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_output_filename",
          "The name of the outputted gnuplot file.",
          command::Parameter
          (
            "output_filename", "string which is the name of the outputted gnuplot file", "output.gnuplot"
          )
        )
      );

      return s_flag;
    }

    //! @brief flag for specifying the interval of tics on the x axis
    //! @return flag for specifying the interval of tics on the x axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagTicIntervalX()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_tics_x",
          "The interval of tics on the x axis.",
          command::Parameter( "tick_interval", "size_t which is the interval of tics on the x axis", "1")
        )
      );

      return s_flag;
    }

    //! @brief flag for specifying the interval of tics on the y axis
    //! @return flag for specifying the interval of tics on the y axis
    const util::ShPtr< command::FlagInterface> &Gnuplot::GetFlagTicIntervalY()
    {
      static const util::ShPtr< command::FlagInterface> s_flag
      (
        new command::FlagStatic
        (
          "gnuplot_tics_y",
          "The interval of tics on the y axis.",
          command::Parameter( "tick_interval", "size_t which is the interval of tics on the y axis", "1")
        )
      );

      return s_flag;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Gnuplot::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief string for axis
    //! @param AXIS the axis
    //! @return the string for the axis
    const std::string &Gnuplot::GetAxisString( const coord::Axis &AXIS)
    {
      static const std::string s_axis_string[] = { "x", "y", "cb"};
      return s_axis_string[ AXIS.GetIndex()];
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief tics from a vector of binning
    //! @param BINNING vector of bin centers
    //! @param BINS_PER_TIC number of bins that should be associated with one tic
    //! @param FORMAT format converting double to string
    //! @return vector of tic labels
    storage::Vector< std::string> Gnuplot::TicsFromBinning
    (
      const linal::Vector< double> &BINNING,
      const size_t BINS_PER_TIC,
      const util::Format &FORMAT
    )
    {
      storage::Vector< std::string> tics;
      for( const double *ptr( BINNING.Begin()), *ptr_end( BINNING.End()); ptr < ptr_end; ptr += BINS_PER_TIC)
      {
        tics.PushBack( FORMAT( *ptr));
      }

      // end
      return tics;
    }

  } // namespace math
} // namespace bcl
