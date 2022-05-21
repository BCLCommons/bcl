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

#ifndef BCL_MATH_GNUPLOT_H_
#define BCL_MATH_GNUPLOT_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically
#include "command/bcl_command.fwd.hh"
#include "coord/bcl_coord.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Gnuplot
    //! @brief Interface class for plots in gnuplot
    //! @details Allows setting the basic options available to many types of plots in gnuplot
    //!
    //! @remarks example unnecessary
    //! @author alexanns
    //! @date Jan 5, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API Gnuplot :
      public util::ObjectInterface
    {

    private:

    //////////
    // data //
    //////////

    public:

      //! @brief flag for setting number of pixels in x direction
      //! @return flag for setting number of pixels in x direction
      static const util::ShPtr< command::FlagInterface> &GetFlagXPixels();

      //! @brief flag for setting number of pixels in y direction
      //! @return flag for setting number of pixels in y direction
      static const util::ShPtr< command::FlagInterface> &GetFlagYPixels();

      //! @brief flag for setting the title of the plot
      //! @return flag for setting the title of the plot
      static const util::ShPtr< command::FlagInterface> &GetFlagTitle();

      //! @brief flag for determining whether the key is shown or not in the plot
      //! @return flag for determining whether the key is shown or not in the plot
      static const util::ShPtr< command::FlagInterface> &GetFlagSetKey();

      //! @brief flag for setting the label on the x axis
      //! @return flag for setting the label on the x axis
      static const util::ShPtr< command::FlagInterface> &GetFlagXLabel();

      //! @brief flag for setting the label on the y axis
      //! @return flag for setting the label on the y axis
      static const util::ShPtr< command::FlagInterface> &GetFlagYLabel();

      //! @brief flag for setting the maximum value on the x axis
      //! @return flag for setting the maximum value on the x axis
      static const util::ShPtr< command::FlagInterface> &GetFlagXMax();

      //! @brief flag for setting the maximum value on the y axis
      //! @return flag for setting the maximum value on the y axis
      static const util::ShPtr< command::FlagInterface> &GetFlagYMax();

      //! @brief flag for setting the minimum value on the x axis
      //! @return flag for setting the minimum value on the x axis
      static const util::ShPtr< command::FlagInterface> &GetFlagXMin();

      //! @brief flag for setting the minimum value on the y axis
      //! @return flag for setting the minimum value on the y axis
      static const util::ShPtr< command::FlagInterface> &GetFlagYMin();

      //! @brief flag for setting series names to be plotted
      //! @return flag for setting series names to be plotted
      static const util::ShPtr< command::FlagInterface> &GetFlagSeriesNames();

      //! @brief flag for specifying filenames of tables to be included in distribution plot
      //! @return flag for specifying filenames of tables to be included in distribution plot
      static const util::ShPtr< command::FlagInterface> &GetFlagTableInputFilenames();

      //! @brief flag for specifying columns of the tables to be used for the plot by column name for x axis
      //! @return flag for specifying columns of the tables to be used for the plot by column name for x axis
      static const util::ShPtr< command::FlagInterface> &GetFlagTableColumnsX();

      //! @brief flag for specifying columns of the tables to be used for the plot by column name for y axis
      //! @return flag for specifying columns of the tables to be used for the plot by column name for y axis
      static const util::ShPtr< command::FlagInterface> &GetFlagTableColumnsY();

      //! @brief flag for specifying the name of the outputted gnuplot file
      //! @return flag for specifying the name of the outputted gnuplot file
      static const util::ShPtr< command::FlagInterface> &GetFlagGnuplotOutputFilename();

      //! @brief flag for specifying the interval of tics on the x axis
      //! @return flag for specifying the interval of tics on the x axis
      static const util::ShPtr< command::FlagInterface> &GetFlagTicIntervalX();

      //! @brief flag for specifying the interval of tics on the y axis
      //! @return flag for specifying the interval of tics on the y axis
      static const util::ShPtr< command::FlagInterface> &GetFlagTicIntervalY();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief string for axis
      //! @param AXIS the axis
      //! @return the string for the axis
      static const std::string &GetAxisString( const coord::Axis &AXIS);

      //! @brief set tics for x axis - the two outside tics are placed at either the outer boundary or on the two outer bins
      //! @param TICS_X vector of strings
      //! @param CENTER_TICS center the tics (false, between two bins)
      //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
      //! @return true if the number of tics required matches the number of bins, so that tics are symmetric
      virtual bool SetTicsX( const storage::Vector< std::string> &TICS_X, const bool CENTER_TICS, const size_t NTH_BIN) = 0;

      //! @brief set tics for y axis - the two outside tics are placed at either the outer boundary or on the two outer bins
      //! @param TICS_Y vector of strings
      //! @param CENTER_TICS center the tics (false, to the top and bottom of the bin)
      //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
      //! @return true if the number of tics required matches the number of bins, so that tics are symmetric
      virtual bool SetTicsY( const storage::Vector< std::string> &TICS_Y, const bool CENTER_TICS, const size_t NTH_BIN) = 0;

      //! @brief set the pixel size and ratio
      //! @param PIXEL_X number of pixel of x
      //! @param PIXEL_Y number of pixel of y
      //! @param RATIO ratio for plot y / x
      virtual void SetPixelAndRatio( const size_t PIXEL_X, const size_t PIXEL_Y, const double RATIO) = 0;

      //! @brief set the rotation of the xtics
      //! @param DEGREE_ROTATION degree of rotation - 0 non
      virtual void SetRotationXTics( const double DEGREE_ROTATION) = 0;

      //! @brief set the filename of the generated plot
      //! @param FILENAME filename
      virtual void SetFilename( const std::string &FILENAME) = 0;

      //! @brief set the font
      //! @param FONT string of font name
      //! @param FONT_SIZE size of font in pixel
      virtual void SetFont( const std::string &FONT, const size_t FONT_SIZE) = 0;

      //! @brief if set, then when WriteScript is called, the WritePreHeader string will be written
      virtual void SetWritePreHeader( const bool SET) = 0;

      //! @brief if set, then when WriteScript is called, the WriteHeader string will be written
      virtual void SetWriteHeader( const bool SET) = 0;

      //! @brief if set, then when WriteScript is called, the WriteBoxes string will be written
      virtual void SetWriteBoxes( const bool SET) = 0;

      //! @brief if set, then when WriteScript is called, the WriteData string will be written
      virtual void SetWriteData( const bool SET) = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief tics from a vector of binning
      //! @param BINNING vector of bin centers
      //! @param BINS_PER_TIC number of bin that should be associated with one bin
      //! @param FORMAT format converting double to string
      //! @return vector of tic labels
      static storage::Vector< std::string> TicsFromBinning
      (
        const linal::Vector< double> &BINNING,
        const size_t BINS_PER_TIC,
        const util::Format &FORMAT
      );

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write unchangeable information necessary for the plot
      //! @param OSTREAM stream to be written to
      //! @param TITLE title for the heat map
      virtual std::ostream &WritePreHeader( std::ostream &OSTREAM) const = 0;

      //! @brief write a header
      //! @param OSTREAM stream to be written to
      virtual std::ostream &WriteHeader( std::ostream &OSTREAM) const = 0;

      //! @brief write the whole information for an axis
      //! @param OSTREAM the stream to write to
      //! @param AXIS the axis that is written
      //! @param LABEL the axis label - if empty, it will be disabled
      //! @param MIN the minimal value for this axis
      //! @param MAX the max value on this axis
      //! @param TICS tics for the axis
      //! @param CENTER_TICS center the tics (false, to the top and bottom of the bin)
      //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
      virtual std::ostream &WriteAxisInfo
      (
        std::ostream &OSTREAM,
        const coord::Axis &AXIS,
        const std::string &LABEL,
        const double MIN,
        const double MAX,
        const storage::Vector< std::string> &TICS,
        const bool CENTER_TICS,
        const size_t NTH_BIN
      ) const = 0;

      //! @brief write the data
      //! @param OSTREAM the stream to write to
      //! @return the stream written to
      virtual std::ostream &WriteData( std::ostream &OSTREAM) const = 0;

      //! @brief write the gnuplot script
      //! @param OSTREAM the stream to write to
      virtual void WriteScript( std::ostream &OSTREAM) const = 0;

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

    private:

    }; // class Gnuplot

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_GNUPLOT_H_ 
