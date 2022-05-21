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

#ifndef BCL_MATH_GNUPLOT_HEATMAP_H_
#define BCL_MATH_GNUPLOT_HEATMAP_H_

// include the namespace header
#include "bcl_math.h"

// include other forward headers - sorted alphabetically
#include "coord/bcl_coord.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_math_gnuplot.h"
#include "linal/bcl_linal_matrix_operations.h"
#include "storage/bcl_storage_vector.h"
#include "storage/bcl_storage_vector_nd.h"
#include "util/bcl_util_object_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class GnuplotHeatmap
    //! @brief this class creates a script file to produce a graphical heatmap from a written gnuplot script
    //!
    //! @see @link example_math_gnuplot_heatmap.cpp @endlink
    //! @author woetzen
    //! @date Aug 26, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API GnuplotHeatmap :
      public Gnuplot
    {

    public:

      //! @enum PaletteTypes
      //! @brief enumerates different color palettes
      enum PaletteTypes
      {
        e_Default,
        e_GreeRedViolet,
        e_GreyScale,
        e_BlackBlueVioletYellowWhite,
        e_BlueGreenYellowRed,
        e_WhiteBlueGreenRed,
        e_BlueGreenRedWhite,
        e_None,
        s_NumberPaletteTypes
      };

      //! @brief GetEnumDescriptor provides the name of ENUM
      //! @param ENUM - the enum for which a name is desired
      static const std::string &GetPaletteTypeString( const PaletteTypes &ENUM);

      //! @brief Gnuplot Palette String provides the string to set the palette
      //! @param ENUM - the enum for which the palette is desired
      static const std::string &GetPaletteString( const PaletteTypes &ENUM);

    private:

    //////////
    // data //
    //////////

      linal::Matrix< double> m_Data; //!< data in heatmap; x - cols; y - rows
      double                 m_MinZ; //!< lowest value in z
      double                 m_MaxZ; //!< highest value in z

      std::string m_Title;  //!< title of heatmap
      std::string m_LabelX; //!< label y
      std::string m_LabelY; //!< label x
      std::string m_LabelZ; //!< label z

      storage::Vector< std::string> m_TicsX; //!< tics for x-axis
      storage::Vector< std::string> m_TicsY; //!< tics for y-axis
      bool                          m_TicsXCenter; //!< tics for x axis in center of bin
      bool                          m_TicsYCenter; //!< tics for x axis in center of bin
      size_t                        m_BinsPerTicX; //!< number of bins per tic x
      size_t                        m_BinsPerTicY; //!< number of bins per tic y

      size_t m_PixelX; //!< size of produced plot in pixel x
      size_t m_PixelY; //!< size of produced plot in pixel y
      double m_Ratio;  //!< ratio of y and x axis length
      double m_RotateXTics; //!< degree by which x tics are rotated

      std::string m_Filename; //!< filename of png that will be generated

      std::string m_Font;     //!< font to be used
      size_t      m_FontSize; //!< size of font

      //! [(x,y),(x,y)] coordinates for each box - one for the lower left corner, one for the upper right, respectively
      storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > > m_BoxCoords;
      //! strings specifying format of the boxes (linewidth, etc.)
      storage::Vector< std::string> m_BoxCoordsSpecifications;

      //! bool indicating if the color box should be displayed or not - the box that shows the color scale
      bool m_ShowColorBox;

      //! booleans for what to write
      bool m_WritePreHeader;
      bool m_WriteHeader;
      bool m_WriteBoxes;
      bool m_WriteData;

      PaletteTypes m_Palette; //!< palette that is desired
      bool m_ShowXTics;      //!< true if the x tics should be shown in the plot
      bool m_NoMirrorTics;   //!< true if the tics should not be shown on the opposite side of the plot
      double m_TopMargin;    //!< ratio to which this plot extends to the top    edge of the png
      double m_BottomMargin; //!< ratio to which this plot extends to the bottom edge of the png
      double m_RightMargin;  //!< ratio to which this plot extends to the right  edge of the png
      double m_LeftMargin;   //!< ratio to which this plot extends to the left   edge of the png

    public:

      //! symbol that represents Angstrom
      static const std::string s_AngstromSymbolGnuplot;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      GnuplotHeatmap();

      //! @brief virtual copy constructor
      GnuplotHeatmap *Clone() const
      {
        return new GnuplotHeatmap( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName< GnuplotHeatmap>();
      }

      //! @brief set from a matrix
      //! @param MATRIC the matrix of values to use
      void SetFromMatrix( const linal::Matrix< double> &MATRIX);

      //! @brief set heatmap from histogram 1D
      //! @param HISTOGRAM1D 1d histogram
      //! @param VERTICAL true - vertical, false horizontal
      //! @param CENTER_TICS center the tics (false, between two bins)
      void SetFromHistogram( const Histogram &HISTOGRAM, const bool VERTICAL, const bool CENTER_TICS);

      //! @brief set heatmap from histograms
      //! @param HISTOGRAMS siptr vector of histograms
      //! @param VERTICAL true - vertical, false horizontal
      //! @param CENTER_TICS center the tics (false, between two bins)
      void SetFromHistograms( const util::SiPtrVector< const Histogram> &HISTOGRAMS, const bool VERTICAL, const bool CENTER_TICS);

      //! @brief set heatmap from unary function
      //! @param FUNCTION the function to use
      //! @param NUMBER_FUNCTION_VALUES number of function values to be calculated
      //! @param START first argument
      //! @param DELTA difference between two consecutive arguments
      //! @param VERTICAL true - vertical, false - horizontal
      //! @param CENTER_TICS center the tics (false, between two bins)
      //! @param NORMALIZE_EACH_FUNCTION if true, sum of the function values will be normalized to one
      //! @param TICS_PER_DELTA number of delta-units per x-axis tic
      void SetFromFunction
      (
        const FunctionInterfaceSerializable< double, double> &FUNCTION,
        const size_t NUMBER_FUNCTION_VALUES,
        const double START,
        const double DELTA,
        const bool VERTICAL,
        const bool CENTER_TICS,
        const bool NORMALIZE_EACH_FUNCTION,
        const size_t TICS_PER_DELTA = 1
      );

      //! @brief set heatmap from cubic spline
      //! @param CUBIC_SPLINE the spline to use
      //! @param VERTICAL true - vertical, false horizontal
      //! @param CENTER_TICS center the tics (false, between two bins)
      void SetFromCubicSpline
      (
        const CubicSpline &CUBIC_SPLINE,
        const bool VERTICAL,
        const bool CENTER_TICS
      );

      //! @brief set heatmap from cubic spline
      //! @param CUBIC_SPLINE the spline to use
      //! @param VERTICAL true - vertical, false horizontal
      //! @param CENTER_TICS center the tics (false, between two bins)
      void SetFromCubicSpline
      (
        const CubicSplineDamped &CUBIC_SPLINE,
        const bool VERTICAL,
        const bool CENTER_TICS
      );

      //! @brief set heatmap from unary functions
      //! @param FUNCTIONS the functions to use
      //! @param NUMBER_FUNCTION_ARGUMENTS number of function values to be calculated
      //! @param START first argument
      //! @param DELTA difference between two consecutive arguments
      //! @param VERTICAL true - vertical, false horizontal
      //! @param CENTER_TICS center the tics (false, between two bins)
      //! @param NORMALIZE_EACH_FUNCTION if true, sum of the values for each function will be normalized to one
      //! @param TICS_PER_DELTA number of delta-units per x-axis tic
      void SetFromFunctions
      (
        const util::SiPtrVector< const FunctionInterfaceSerializable< double, double> > &FUNCTIONS,
        const size_t NUMBER_FUNCTION_VALUES,
        const double START,
        const double DELTA,
        const bool VERTICAL,
        const bool CENTER_TICS,
        const bool NORMALIZE_EACH_FUNCTION,
        const size_t TICS_PER_DELTA = 1
      );

      //! @brief set heatmap from cubic splines
      //! @param CUBIC_SPLINES the splines to use
      //! @param VERTICAL true - vertical, false horizontal
      //! @param CENTER_TICS center the tics (false, between two bins)
      void SetFromCubicSplines
      (
        const util::SiPtrVector< const CubicSpline> &CUBIC_SPLINES,
        const bool VERTICAL,
        const bool CENTER_TICS
      );

      //! @brief set heatmap from cubic splines
      //! @param CUBIC_SPLINES the splines to use
      //! @param VERTICAL true - vertical, false horizontal
      //! @param CENTER_TICS center the tics (false, between two bins)
      void SetFromCubicSplines
      (
        const util::SiPtrVector< const CubicSplineDamped> &CUBIC_SPLINES,
        const bool VERTICAL,
        const bool CENTER_TICS
      );

      //! @brief set heatmap from histogram 2D
      //! @param HISTOGRAM2D 2d histogram
      //! @param CENTER_TICS_X center the tics x axis (false, between two bins)
      //! @param CENTER_TICS_Y center the tics y axis (false, between two bins)
      void SetFromHistogram
      (
        const Histogram2D &HISTOGRAM2D,
        const bool CENTER_TICS_X,
        const bool CENTER_TICS_Y
      );

      //! @brief set heatmap from binary function
      //! @param BINARY_FUNCTION the function to use
      //! @param NUMBER_FUNCTION_VALUES_X number of function values to be calculated in x
      //! @param START_X first argument x
      //! @param DELTA_X difference between two consecutive arguments x
      //! @param NUMBER_FUNCTION_VALUES_Y number of function values to be calculated in y
      //! @param START_Y first argument y
      //! @param DELTA_Y difference between two consecutive arguments y
      //! @param CENTER_TICS_X center the tics x axis (false, between two bins)
      //! @param CENTER_TICS_Y center the tics y axis (false, between two bins)
      void SetFromBinaryFunction
      (
        const BinaryFunctionInterface< double, double, double> &BINARY_FUNCTION,
        const size_t NUMBER_FUNCTION_VALUES_X,
        const double START_X,
        const double DELTA_X,
        const size_t NUMBER_FUNCTION_VALUES_Y,
        const double START_Y,
        const double DELTA_Y,
        const bool CENTER_TICS_X,
        const bool CENTER_TICS_Y
      );

      //! @brief set heatmap from bicubic spline
      //! @param BICUBIC_SPLINE the bicubic spline to use
      //! @param CENTER_TICS_X center the tics x axis (false, between two bins)
      //! @param CENTER_TICS_Y center the tics y axis (false, between two bins)
      void SetFromBicubicSpline
      (
        const BicubicSpline &BICUBIC_SPLINE,
        const bool CENTER_TICS_X,
        const bool CENTER_TICS_Y
      );

      //! @brief set title and labels
      //! @details leave any of those empty to not label the axis
      //! @param TITLE title for the map
      //! @param LABEL_X label for the x axis
      //! @param LABEL_Y label for the y axis
      //! @param LABEL_Z label for the z axis
      void SetTitleAndLabel
      (
        const std::string &TITLE,
        const std::string &LABEL_X,
        const std::string &LABEL_Y,
        const std::string &LABEL_Z
      );

      //! @brief set the min and max of z
      //! @param MIN_Z minimum on z - if undefined, pick the lowest automatically
      //! @param MAX_Z maximum on z - if undefined, pick the highest automatically
      void SetMinMaxZ( const double MIN_Z, const double MAX_Z);

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
      void SetRotationXTics( const double DEGREE_ROTATION)
      {
        m_RotateXTics = DEGREE_ROTATION;
      }

      //! @brief get the data
      //! @return matrix with data
      const linal::Matrix< double> &GetData() const
      {
        return m_Data;
      }

      //! @brief set the filename of the generated plot
      //! @param FILENAME filename
      void SetFilename( const std::string &FILENAME)
      {
        m_Filename = FILENAME;
      }

      //! @brief set the font
      //! @param FONT string of font name
      //! @param FONT_SIZE size of font in pixel
      void SetFont( const std::string &FONT, const size_t FONT_SIZE);

      //! @brief sets the coordinates for boxes
      //! @param BOX_COORDS the lower left, upper right coordinates for each box [(x,y),(x,y)]
      //! @param X_RANGE the range in the x direction of the plot - the total size of the x direction
      //! @param Y_RANGE the range in the y direction of the plot -  the total size of the y direction
      //! @param X_MIN the minimum value in the x direction of the heat map
      //! @param Y_MIN the minimum value in the y direction of the heat map
      void SetBoxes
      (
        const storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > > &BOX_COORDS,
        const double X_RANGE, const double Y_RANGE, const double X_MIN, const double Y_MIN
      );

      //! @brief sets the specifications for boxes
      //! @param SPECIFICATIONS strings which specify the format of the boxes
      void SetBoxesSpecifications( const storage::Vector< std::string> &SPECIFICATIONS);

      //! @brief bool indicating if the color box should be displayed or not - the box that shows the color scale
      //! @param SET_COLOR_BOX bool true if color box should be set - false otherwise
      void SetShowColorBox( const bool SET_COLOR_BOX)
      {
        m_ShowColorBox = SET_COLOR_BOX;
      }

      //! @brief if set, then when WriteScript is called, the WritePreHeader string will be written. True by default
      void SetWritePreHeader( const bool SET);

      //! @brief if set, then when WriteScript is called, the WriteHeader    string will be written. True by default
      void SetWriteHeader( const bool SET);

      //! @brief if set, then when WriteScript is called, the WritePalettes  string will be written. True by default
      void SetWritePalettes( const bool SET);

      //! @brief if set, then when WriteScript is called, the WriteBoxes     string will be written. True by default
      void SetWriteBoxes( const bool SET);

      //! @brief if set, then when WriteScript is called, the WriteData      string will be written. True by default
      void SetWriteData( const bool SET);

      //! @brief set, the palette to use
      //! @param PALETTE the color palette used by the splot
      void SetPalette( const PaletteTypes PALETTE);

      //! @brief if set, then the x tics will be shown
      void SetShowXTics( const bool SET);

      //! @brief if set, then the the tics won't be mirrored on the opposite side of the axis
      void SetNoMirrorTics( const bool SET);

      //! @brief set the ratio to which the actual plot area extends to the edge of the png
      //! @param TOP    ratio to which this plot extends to the top    edge of the png
      //! @param BOTTOM ratio to which this plot extends to the bottom edge of the png
      //! @param RIGHT  ratio to which this plot extends to the right  edge of the png
      //! @param LEFT   ratio to which this plot extends to the left   edge of the png
      void SetMargins( const double TOP, const double BOTTOM, const double RIGHT, const double LEFT);

    ////////////////
    // operations //
    ////////////////

      //! @brief multiply the heatmap
      //! @details can be used to get similar scaling between different heat maps
      //! @param FACTOR
      void Muliply( const double FACTOR);

    //////////////////////
    // input and output //
    //////////////////////

      //! @brief write unchangeable information necessary for the plot
      //! @param OSTREAM stream to be written to
      //! @param TITLE title for the heat map
      std::ostream &WritePreHeader( std::ostream &OSTREAM) const;

      //! @brief write a choice of palettes
      //! @param OSTREAM stream to be written to
      std::ostream &WritePalettes( std::ostream &OSTREAM) const;

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

      //! @brief write any boxes to the script for display
      //! @param OSTREAM the stream the boxes will be written to
      std::ostream &WriteBoxes( std::ostream &OSTREAM) const;

      //! @brief write the data
      //! @param OSTREAM the stream to write to
      //! @return the stream written to
      std::ostream &WriteData( std::ostream &OSTREAM) const;

      //! @brief write the gnuplot script
      //! @param OSTREAM the stream to write to
      void WriteScript( std::ostream &OSTREAM) const;

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT - number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    }; // class GnuplotHeatmap

  } // namespace math
} // namespace bcl

#endif // BCL_MATH_GNUPLOT_HEATMAP_H_
