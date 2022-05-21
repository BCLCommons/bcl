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

#ifndef BCL_MATH_GNUPLOT_MULTIPLOT_H_
#define BCL_MATH_GNUPLOT_MULTIPLOT_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_math_gnuplot.h"
#include "util/bcl_util_sh_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GnuplotMultiplot
    //! @brief Allows putting multiple plots into a single output picture
    //! @details Plots can arranged in columns and rows as desired. Uses list of gnuplots to make the individual plots.
    //!
    //! @see @link example_math_gnuplot_multiplot.cpp @endlink
    //! @author alexanns
    //! @date Jan 5, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GnuplotMultiplot :
      public Gnuplot
    {

    private:

    //////////
    // data //
    //////////

      //! the list of plots making up this multiplot
      util::ShPtrList< Gnuplot> m_Plots;

      //! the number of rows the plots will be divided into
      size_t m_Rows;

      //! the number of columns the plots will be divided into
      size_t m_Cols;

      std::string m_Filename; //!< filename of png that will be generated

      std::string m_Font;     //!< font to be used
      size_t      m_FontSize; //!< size of font

      size_t m_PixelX; //!< size of produced plot in pixel x
      size_t m_PixelY; //!< size of produced plot in pixel y
      double m_Ratio;  //!< ratio of y and x axis length

      //! booleans for what to write
      bool m_WritePreHeader;
      bool m_WritePalettes;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      GnuplotMultiplot();

      //! @brief Clone function
      //! @return pointer to new GnuplotMultiplot
      GnuplotMultiplot *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief set tics for x axis - the two outside tics are placed at either the outer boundary or on the two outer bins
      //! @param TICS_X vector of strings
      //! @param CENTER_TICS center the tics (false, between two bins)
      //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
      //! @return true if the number of tics required matches the number of bins, so that tics are symmetric
      bool SetTicsX( const storage::Vector< std::string> &TICS_X, const bool CENTER_TICS, const size_t NTH_BIN);

      //! @brief set tics for y axis - the two outside tics are placed at either the outer boundary or on the two outer bins
      //! @param TICS_Y vector of strings
      //! @param CENTER_TICS center the tics (false, to the top and bottom of the bin)
      //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
      //! @return true if the number of tics required matches the number of bins, so that tics are symmetric
      bool SetTicsY( const storage::Vector< std::string> &TICS_Y, const bool CENTER_TICS, const size_t NTH_BIN);

      //! @brief set the pixel size and ratio
      //! @param PIXEL_X number of pixel of x
      //! @param PIXEL_Y number of pixel of y
      //! @param RATIO ratio for plot y / x
      void SetPixelAndRatio( const size_t PIXEL_X, const size_t PIXEL_Y, const double RATIO);

      //! @brief set the rotation of the xtics
      //! @param DEGREE_ROTATION degree of rotation - 0 non
      void SetRotationXTics( const double DEGREE_ROTATION);

      //! @brief set the filename of the generated plot
      //! @param FILENAME filename
      void SetFilename( const std::string &FILENAME);

      //! @brief set the font
      //! @param FONT string of font name
      //! @param FONT_SIZE size of font in pixel
      void SetFont( const std::string &FONT, const size_t FONT_SIZE);

      //! @brief if set, then when WriteScript is called, the WritePreHeader string will be written
      void SetWritePreHeader( const bool SET);

      //! @brief if set, then when WriteScript is called, the WriteHeader string will be written
      void SetWriteHeader( const bool SET);

      //! @brief if set, then when WriteScript is called, the WriteBoxes string will be written
      void SetWriteBoxes( const bool SET);

      //! @brief if set, then when WriteScript is called, the WriteData string will be written
      void SetWriteData( const bool SET);

      //! @brief set the row and column layout of the multiplot
      void SetRowsCols( size_t ROWS, size_t COLS);

    ////////////////
    // operations //
    ////////////////

      //! @brief add a plot to the multiplot
      //! @param PLOT the plot which will be added to the multiplot
      void Insert( const util::ShPtr< Gnuplot> &PLOT);

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write unchangeable information necessary for the plot
      //! @param OSTREAM stream to be written to
      //! @param TITLE title for the heat map
      std::ostream &WritePreHeader( std::ostream &OSTREAM) const;

      //! @brief write a header
      //! @param OSTREAM stream to be written to
      std::ostream &WriteHeader( std::ostream &OSTREAM) const;

      //! @brief write the whole information for an axis
      //! @param OSTREAM the stream to write to
      //! @param AXIS the axis that is written
      //! @param LABEL the axis label - if empty, it will be disabled
      //! @param MIN the minimal value for this axis
      //! @param MAX the max value on this axis
      //! @param TICS tics for the axis
      //! @param CENTER_TICS center the tics (false, to the top and bottom of the bin)
      //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
      std::ostream &WriteAxisInfo
      (
        std::ostream &OSTREAM,
        const coord::Axis &AXIS,
        const std::string &LABEL,
        const double MIN,
        const double MAX,
        const storage::Vector< std::string> &TICS,
        const bool CENTER_TICS,
        const size_t NTH_BIN
      ) const;

      //! @brief write the data
      //! @param OSTREAM the stream to write to
      //! @return the stream written to
      std::ostream &WriteData( std::ostream &OSTREAM) const;

      //! @brief write the gnuplot script
      //! @param OSTREAM the stream to write to
      void WriteScript( std::ostream &OSTREAM) const;

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

    }; // class GnuplotMultiplot

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_GNUPLOT_MULTIPLOT_H_ 
