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
#include "math/bcl_math_gnuplot_multiplot.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> GnuplotMultiplot::s_Instance
    (
      GetObjectInstances().AddInstance( new GnuplotMultiplot())
    );
  
  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    GnuplotMultiplot::GnuplotMultiplot() :
      m_Plots(),
      m_Rows( 1),
      m_Cols( 1),
      m_Filename( "heatmap"),
      m_Font( "Arial"),
      m_FontSize( 12),
      m_PixelX( 600),
      m_PixelY( 400),
      m_Ratio( util::GetUndefinedDouble()),
      m_WritePreHeader( true),
      m_WritePalettes( true)
    {
    }

    //! @brief Clone function
    //! @return pointer to new GnuplotMultiplot
    GnuplotMultiplot *GnuplotMultiplot::Clone() const
    {
      return new GnuplotMultiplot( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GnuplotMultiplot::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief set tics for x axis - the two outside tics are placed at either the outer boundary or on the two outer bins
    //! @param TICS_X vector of strings
    //! @param CENTER_TICS center the tics (false, between two bins)
    //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
    //! @return true if the number of tics required matches the number of bins, so that tics are symmetric
    bool GnuplotMultiplot::SetTicsX
    (
      const storage::Vector< std::string> &TICS_X, const bool CENTER_TICS, const size_t NTH_BIN
    )
    {
      return true;
    }

    //! @brief set tics for y axis - the two outside tics are placed at either the outer boundary or on the two outer bins
    //! @param TICS_Y vector of strings
    //! @param CENTER_TICS center the tics (false, to the top and bottom of the bin)
    //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
    //! @return true if the number of tics required matches the number of bins, so that tics are symmetric
    bool GnuplotMultiplot::SetTicsY
    (
      const storage::Vector< std::string> &TICS_Y, const bool CENTER_TICS, const size_t NTH_BIN
    )
    {
      return true;
    }

    //! @brief set the pixel size and ratio
    //! @param PIXEL_X number of pixel of x
    //! @param PIXEL_Y number of pixel of y
    //! @param RATIO ratio for plot y / x
    void GnuplotMultiplot::SetPixelAndRatio( const size_t PIXEL_X, const size_t PIXEL_Y, const double RATIO)
    {
      m_PixelX = PIXEL_X;
      m_PixelY = PIXEL_Y;
      m_Ratio  = RATIO;
    }

    //! @brief set the rotation of the xtics
    //! @param DEGREE_ROTATION degree of rotation - 0 non
    void GnuplotMultiplot::SetRotationXTics( const double DEGREE_ROTATION)
    {
      return;
    }

    //! @brief set the filename of the generated plot
    //! @param FILENAME filename
    void GnuplotMultiplot::SetFilename( const std::string &FILENAME)
    {
      m_Filename = FILENAME;
    }

    //! @brief set the font
    //! @param FONT string of font name
    //! @param FONT_SIZE size of font in pixel
    void GnuplotMultiplot::SetFont( const std::string &FONT, const size_t FONT_SIZE)
    {
      m_Font = FONT;
      m_FontSize = FONT_SIZE;
    }

    //! @brief if set, then when WriteScript is called, the WritePreHeader string will be written
    void GnuplotMultiplot::SetWritePreHeader( const bool SET)
    {
      m_WritePreHeader = SET;
    }

    //! @brief if set, then when WriteScript is called, the WriteHeader string will be written
    void GnuplotMultiplot::SetWriteHeader( const bool SET)
    {
      return;
    }

    //! @brief if set, then when WriteScript is called, the WriteBoxes string will be written
    void GnuplotMultiplot::SetWriteBoxes( const bool SET)
    {
      return;
    }

    //! @brief if set, then when WriteScript is called, the WriteData string will be written
    void GnuplotMultiplot::SetWriteData( const bool SET)
    {
      return;
    }

    //! @brief set the row and column layout of the multiplot
    void GnuplotMultiplot::SetRowsCols( size_t ROWS, size_t COLS)
    {
      m_Rows = ROWS;
      m_Cols = COLS;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief add a plot to the multiplot
    //! @param PLOT the plot which will be added to the multiplot
    void GnuplotMultiplot::Insert( const util::ShPtr< Gnuplot> &PLOT)
    {
      m_Plots.PushBack( PLOT);
    }

  ///////////////
  // operators //
  ///////////////

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write unchangeable information necessary for the plot
    //! @param OSTREAM stream to be written to
    //! @param TITLE title for the heat map
    std::ostream &GnuplotMultiplot::WritePreHeader( std::ostream &OSTREAM) const
    {
      // write comments
      OSTREAM << "# BCL generated heatmap\n";
      OSTREAM << "set terminal png transparent ";
      if( !m_Font.empty())
      {
        OSTREAM << "font \"" << m_Font << "\" " << m_FontSize;
      }
      OSTREAM << " size " << m_PixelX << ',' << m_PixelY << '\n';
      if( util::IsDefined( m_Ratio))
      {
        OSTREAM << "set size ratio " << m_Ratio << '\n';
      }
      OSTREAM << "set output \"" << m_Filename << ".png\"\n";
      OSTREAM << "set encoding iso_8859_1\n";
      OSTREAM << "set multiplot layout " << m_Rows << "," << m_Cols << "\n";

      return OSTREAM;
    }

    //! @brief write a header
    //! @param OSTREAM stream to be written to
    std::ostream &GnuplotMultiplot::WriteHeader( std::ostream &OSTREAM) const
    {
      return OSTREAM;
    }

    //! @brief write the whole information for an axis
    //! @param OSTREAM the stream to write to
    //! @param AXIS the axis that is written
    //! @param LABEL the axis label - if empty, it will be disabled
    //! @param MIN the minimal value for this axis
    //! @param MAX the max value on this axis
    //! @param TICS tics for the axis
    //! @param CENTER_TICS center the tics (false, to the top and bottom of the bin)
    //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
    std::ostream &GnuplotMultiplot::WriteAxisInfo
    (
      std::ostream &OSTREAM,
      const coord::Axis &AXIS,
      const std::string &LABEL,
      const double MIN,
      const double MAX,
      const storage::Vector< std::string> &TICS,
      const bool CENTER_TICS,
      const size_t NTH_BIN
    ) const
    {
      return OSTREAM;
    }

    //! @brief write the data
    //! @param OSTREAM the stream to write to
    //! @return the stream written to
    std::ostream &GnuplotMultiplot::WriteData( std::ostream &OSTREAM) const
    {
      return OSTREAM;
    }

    //! @brief write the gnuplot script
    //! @param OSTREAM the stream to write to
    void GnuplotMultiplot::WriteScript( std::ostream &OSTREAM) const
    {
      if( m_WritePreHeader)
      {
        WritePreHeader( OSTREAM);
      }

      // iterate through the plots to print them all out
      for
      (
        util::ShPtrList< Gnuplot>::const_iterator plot_itr( m_Plots.Begin()), plot_itr_end( m_Plots.End());
        plot_itr != plot_itr_end; ++plot_itr
      )
      {
        ( *plot_itr)->WriteScript( OSTREAM);
      }
    }

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &GnuplotMultiplot::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &GnuplotMultiplot::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

  } // namespace math
} // namespace bcl
