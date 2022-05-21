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
#include "math/bcl_math_gnuplot_heatmap.h"

// includes from bcl - sorted alphabetically
#include "coord/bcl_coord_axes.h"
#include "math/bcl_math_bicubic_spline.h"
#include "math/bcl_math_cubic_spline.h"
#include "math/bcl_math_cubic_spline_damped.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_histogram_2d.h"
#include "math/bcl_math_limits.h"
#include "math/bcl_math_running_average.h"
#include "util/bcl_util_logger_interface.h"
#include "util/bcl_util_si_ptr_vector.h"
// external includes - sorted alphabetically

namespace bcl
{
  namespace math
  {

    //! @brief GetEnumDescriptor provides the name of ENUM
    //! @param ENUM - the enum for which a name is desired
    const std::string &GnuplotHeatmap::GetPaletteTypeString( const PaletteTypes &ENUM)
    {
      static const std::string s_descriptors[] =
      {
        "Default",
        "GreeRedViolet",
        "GreyScale",
        "BlackBlueVioletYellowWhite",
        "BlueGreenYellowRed",
        "WhiteBlueGreenRed",
        "BlueGreenRedWhite",
        "None",
        GetStaticClassName< PaletteTypes>()
      };

      return s_descriptors[ ENUM];
    }

    //! @brief Gnuplot Palette String provides the string to set the palette
    //! @param ENUM - the enum for which the palette is desired
    const std::string &GnuplotHeatmap::GetPaletteString( const PaletteTypes &ENUM)
    {
      static const std::string s_palette_strings[] =
      {
        "rgbformulae 22, 13, -31",
        "rgbformulae 3, 11, 6 # green-red-violet",
        "grey negative",
        "rgbformulae 30,31,32 # color printable on gray (black-blue-violet-yellow-white)",
        "rgbformulae 33,13,10 # rainbow (blue-green-yellow-red)",
        "defined (0 1 1 1, 0.00001 0 0 1, 1 0 1 0, 2 1 0 0)",
        "defined (0 0 0 1, 1 0 1 0, 2 1 0 0, 2.00001 1 1 1)",
        "#None",
        "#no palette set"
      };

      return s_palette_strings[ ENUM];
    }

  //////////
  // data //
  //////////

    //! symbol that represents Angstrom
    const std::string GnuplotHeatmap::s_AngstromSymbolGnuplot = "\\305";

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    GnuplotHeatmap::GnuplotHeatmap() :
      m_Data( 0, 0, 0.0),
      m_MinZ( util::GetUndefined< double>()),
      m_MaxZ( util::GetUndefined< double>()),
      m_Title(  ""),
      m_LabelX( ""),
      m_LabelY( ""),
      m_LabelZ( ""),
      m_TicsX(    ),
      m_TicsY(    ),
      m_TicsXCenter( true),
      m_TicsYCenter( true),
      m_BinsPerTicX( 1),
      m_BinsPerTicY( 1),
      m_PixelX( 1080),
      m_PixelY(  800),
      m_Ratio( util::GetUndefined< double>()),
      m_RotateXTics( 0.0),
      m_Filename( "heatmap"),
      m_Font( "Arial"),
      m_FontSize( 12),
      m_BoxCoords(),
      m_BoxCoordsSpecifications(),
      m_ShowColorBox( true),
      m_WritePreHeader( true),
      m_WriteHeader( true),
      m_WriteBoxes( true),
      m_WriteData( true),
      m_Palette( e_Default),
      m_ShowXTics( true),
      m_NoMirrorTics( false),
      m_TopMargin   ( util::GetUndefinedDouble()),
      m_BottomMargin( util::GetUndefinedDouble()),
      m_RightMargin ( util::GetUndefinedDouble()),
      m_LeftMargin  ( util::GetUndefinedDouble())
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief set from a matrix
    //! @param MATRIC the matrix of values to use
    void GnuplotHeatmap::SetFromMatrix( const linal::Matrix< double> &MATRIX)
    {
      m_Data = MATRIX;
    }

    //! @brief set heatmap from histogram 1D
    //! @param HISTOGRAM1D 1d histogram
    //! @param VERTICAL true - vertical, false horizontal
    //! @param CENTER_TICS center the tics (false, between two bins)
    void GnuplotHeatmap::SetFromHistogram( const Histogram &HISTOGRAM, const bool VERTICAL, const bool CENTER_TICS)
    {
      SetFromHistograms( util::SiPtrVector< const Histogram>::Create( HISTOGRAM, HISTOGRAM), VERTICAL, CENTER_TICS);
    }

    //! @brief set heatmap from histograms
    //! @param HISTOGRAMS siptr vector of histograms
    //! @param VERTICAL true - vertical, false horizontal
    //! @param CENTER_TICS center the tics (false, between two bins)
    void GnuplotHeatmap::SetFromHistograms( const util::SiPtrVector< const Histogram> &HISTOGRAMS, const bool VERTICAL, const bool CENTER_TICS)
    {
      if( HISTOGRAMS.IsEmpty())
      {
        return;
      }

      const Histogram &template_hist( *HISTOGRAMS.FirstElement());
      const linal::Vector< double> binning
      (
        linal::FillVector< double>
        (
          template_hist.GetNumberOfBins() + size_t( !CENTER_TICS),
          template_hist.GetBoundaries().First() + ( CENTER_TICS ? 0.5 * template_hist.GetBinSize() : 0.0),
          template_hist.GetBinSize()
        )
      );
      const storage::Vector< std::string> tics( Gnuplot::TicsFromBinning( binning, 1, util::Format().W( 4)));

      if( VERTICAL)
      {
        m_Data = linal::Matrix< double>( template_hist.GetNumberOfBins(), HISTOGRAMS.GetSize(), 0.0);
        SetTicsY( tics, CENTER_TICS, 1);
      }
      else
      {
        m_Data = linal::Matrix< double>( HISTOGRAMS.GetSize(), template_hist.GetNumberOfBins(), 0.0);
        SetTicsX( tics, CENTER_TICS, 1);
      }

      size_t count( 0);
      // iterate over all histograms
      for
      (
        util::SiPtrVector< const Histogram>::const_iterator itr( HISTOGRAMS.Begin()), itr_end( HISTOGRAMS.End());
        itr != itr_end;
        ++itr, ++count
      )
      {
        VERTICAL ? m_Data.ReplaceCol( count, ( *itr)->GetHistogram()) : m_Data.ReplaceRow( count, ( *itr)->GetHistogram());
      }
    }

    //! @brief set heatmap from unary function
    //! @param FUNCTION the function to use
    //! @param NUMBER_FUNCTION_VALUES number of function values to be calculated
    //! @param START first argument
    //! @param DELTA difference between two consecutive arguments
    //! @param VERTICAL true - vertical, false - horizontal
    //! @param CENTER_TICS center the tics (false, between two bins)
    //! @param NORMALIZE if true, sum of the function values will be normalized to one
    //! @param TICS_PER_DELTA number of delta-units per x-axis tic
    void GnuplotHeatmap::SetFromFunction
    (
      const FunctionInterfaceSerializable< double, double> &FUNCTION,
      const size_t NUMBER_FUNCTION_VALUES,
      const double START,
      const double DELTA,
      const bool VERTICAL,
      const bool CENTER_TICS,
      const bool NORMALIZE,
      const size_t TICS_PER_DELTA
    )
    {
      SetFromFunctions
      (
        util::SiPtrVector< const FunctionInterfaceSerializable< double, double> >::Create( FUNCTION, FUNCTION),
        NUMBER_FUNCTION_VALUES, START, DELTA, VERTICAL, CENTER_TICS, NORMALIZE, TICS_PER_DELTA
      );
    }

    //! @brief set heatmap from cubic spline
    //! @param CUBIC_SPLINE the spline to use
    //! @param VERTICAL true - vertical, false horizontal
    //! @param CENTER_TICS center the tics (false, between two bins)
    void GnuplotHeatmap::SetFromCubicSpline
    (
      const CubicSpline &CUBIC_SPLINE,
      const bool VERTICAL,
      const bool CENTER_TICS
    )
    {
      SetFromFunction
      (
        CUBIC_SPLINE,
        CUBIC_SPLINE.GetValues().GetSize(),
        CUBIC_SPLINE.GetStart(),
        CUBIC_SPLINE.GetDelta(),
        VERTICAL,
        CENTER_TICS,
        false
      );
    }

    //! @brief set heatmap from cubic spline
    //! @param CUBIC_SPLINE the spline to use
    //! @param VERTICAL true - vertical, false horizontal
    //! @param CENTER_TICS center the tics (false, between two bins)
    void GnuplotHeatmap::SetFromCubicSpline
    (
      const CubicSplineDamped &CUBIC_SPLINE,
      const bool VERTICAL,
      const bool CENTER_TICS
    )
    {
      SetFromFunction
      (
        CUBIC_SPLINE,
        CUBIC_SPLINE.GetXValues().GetSize() * 4,
        CUBIC_SPLINE.GetXValues().First(),
        CUBIC_SPLINE.GetDelta() / 4.0,
        VERTICAL,
        CENTER_TICS,
        false,
        4
      );
    }

    //! @brief set heatmap from unary functions
    //! @param FUNCTIONS the functions to use
    //! @param NUMBER_FUNCTION_ARGUMENTS number of function values to be calculated
    //! @param START first argument
    //! @param DELTA difference between two consecutive arguments
    //! @param VERTICAL true - vertical, false horizontal
    //! @param CENTER_TICS center the tics (false, between two bins)
    //! @param NORMALIZE_EACH_FUNCTION if true, sum of the values for each function will be normalized to one
    //! @param TICS_PER_DELTA number of delta-units per x-axis tic
    void GnuplotHeatmap::SetFromFunctions
    (
      const util::SiPtrVector< const FunctionInterfaceSerializable< double, double> > &FUNCTIONS,
      const size_t NUMBER_FUNCTION_VALUES,
      const double START,
      const double DELTA,
      const bool VERTICAL,
      const bool CENTER_TICS,
      const bool NORMALIZE_EACH_FUNCTION,
      const size_t TICS_PER_DELTA
    )
    {
      if( FUNCTIONS.IsEmpty())
      {
        return;
      }

      const linal::Vector< double> binning
      (
        linal::FillVector< double>
        (
          NUMBER_FUNCTION_VALUES + size_t( !CENTER_TICS), START - ( CENTER_TICS ? 0.0 : 0.5 * DELTA), DELTA
        )
      );
      const storage::Vector< std::string> tics( TicsFromBinning( binning, TICS_PER_DELTA, util::Format()));

      if( VERTICAL)
      {
        m_Data = linal::Matrix< double>( NUMBER_FUNCTION_VALUES, FUNCTIONS.GetSize(), 0.0);
        SetTicsY( tics, CENTER_TICS, TICS_PER_DELTA);
      }
      else
      {
        m_Data = linal::Matrix< double>( FUNCTIONS.GetSize(), NUMBER_FUNCTION_VALUES, 0.0);
        SetTicsX( tics, CENTER_TICS, TICS_PER_DELTA);
      }

      size_t count( 0);
      // iterate over all histograms
      for
      (
        util::SiPtrVector< const FunctionInterfaceSerializable< double, double> >::const_iterator
          itr( FUNCTIONS.Begin()), itr_end( FUNCTIONS.End());
        itr != itr_end;
        ++itr, ++count
      )
      {
        const FunctionInterfaceSerializable< double, double> &function( **itr);
        // create the row from function values in the spline
        linal::Vector< double> val( NUMBER_FUNCTION_VALUES);
        double *val_ptr( val.Begin());
        for( size_t bin( 0); bin < NUMBER_FUNCTION_VALUES; ++bin, ++val_ptr)
        {
          *val_ptr = function( START + bin * DELTA);
        }
        if( NORMALIZE_EACH_FUNCTION)
        {
          val.SetToSum( 1.0);
        }
        VERTICAL ? m_Data.ReplaceCol( count, val) : m_Data.ReplaceRow( count, val);
      }
    }

    //! @brief set heatmap from cubic splines
    //! @param CUBIC_SPLINES the splines to use
    //! @param VERTICAL true - vertical, false horizontal
    //! @param CENTER_TICS center the tics (false, between two bins)
    void GnuplotHeatmap::SetFromCubicSplines
    (
      const util::SiPtrVector< const CubicSpline> &CUBIC_SPLINES,
      const bool VERTICAL,
      const bool CENTER_TICS
    )
    {
      // check if splines are defined
      if( CUBIC_SPLINES.IsEmpty() && !CUBIC_SPLINES.IsDefined())
      {
        return;
      }

      const CubicSpline &template_spline( *CUBIC_SPLINES.FirstElement());
      SetFromFunctions
      (
        util::SiPtrVector< const FunctionInterfaceSerializable< double, double> >( CUBIC_SPLINES),
        template_spline.GetValues().GetSize() * 4,
        template_spline.GetStart(),
        template_spline.GetDelta() / 4.0,
        VERTICAL,
        CENTER_TICS,
        false,
        4
      );
    }

    //! @brief set heatmap from cubic splines
    //! @param CUBIC_SPLINES the splines to use
    //! @param VERTICAL true - vertical, false horizontal
    //! @param CENTER_TICS center the tics (false, between two bins)
    void GnuplotHeatmap::SetFromCubicSplines
    (
      const util::SiPtrVector< const CubicSplineDamped> &CUBIC_SPLINES,
      const bool VERTICAL,
      const bool CENTER_TICS
    )
    {
      // check if splines are defined
      if( CUBIC_SPLINES.IsEmpty() && !CUBIC_SPLINES.IsDefined())
      {
        return;
      }

      RunningAverage< double> min_val;
      RunningAverage< double> max_size;
      for( auto itr( CUBIC_SPLINES.Begin()), itr_end( CUBIC_SPLINES.End()); itr != itr_end; ++itr)
      {
        max_size +=( *itr)->GetXValues().GetSize();
        min_val += ( *itr)->GetXValues().First();
      }
      const CubicSplineDamped &template_spline( *CUBIC_SPLINES.FirstElement());
      SetFromFunctions
      (
        util::SiPtrVector< const FunctionInterfaceSerializable< double, double> >( CUBIC_SPLINES),
        size_t( max_size.GetAverage()) * 4,
        min_val.GetAverage(),
        template_spline.GetDelta() / 4,
        VERTICAL,
        CENTER_TICS,
        false,
        4
      );
    }

    //! @brief set heatmap from histogram 2D
    //! @param HISTOGRAM2D 2d histogram
    //! @param CENTER_TICS_X center the tics x axis (false, between two bins)
    //! @param CENTER_TICS_Y center the tics y axis (false, between two bins)
    void GnuplotHeatmap::SetFromHistogram
    (
      const Histogram2D &HISTOGRAM2D,
      const bool CENTER_TICS_X,
      const bool CENTER_TICS_Y
    )
    {
      m_Data = HISTOGRAM2D.GetHistogram();

      // binning tics
      const linal::Vector< double> binning_x
      (
        linal::FillVector< double>
        (
          HISTOGRAM2D.GetNumberOfBinsX() + size_t( !CENTER_TICS_X),
          HISTOGRAM2D.GetBoundariesX().First() + ( CENTER_TICS_X ? 0.5 * HISTOGRAM2D.GetBinSizeXY().First() : 0.0),
          HISTOGRAM2D.GetBinSizeXY().First()
        )
      );
      const linal::Vector< double> binning_y
      (
        linal::FillVector< double>
        (
          HISTOGRAM2D.GetNumberOfBinsY() + size_t( !CENTER_TICS_Y),
          HISTOGRAM2D.GetBoundariesY().First() + ( CENTER_TICS_Y ? 0.5 * HISTOGRAM2D.GetBinSizeXY().Second() : 0.0),
          HISTOGRAM2D.GetBinSizeXY().Second()
        )
      );
      const storage::Vector< std::string> xtics( TicsFromBinning( binning_x, 1, util::Format().W( 4)));
      const storage::Vector< std::string> ytics( TicsFromBinning( binning_y, 1, util::Format().W( 4)));
      SetTicsX( xtics, CENTER_TICS_X, 1);
      SetTicsY( ytics, CENTER_TICS_Y, 1);
    }

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
    void GnuplotHeatmap::SetFromBinaryFunction
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
    )
    {
      // data
      m_Data = linal::Matrix< double>( NUMBER_FUNCTION_VALUES_Y, NUMBER_FUNCTION_VALUES_X, double( 0.0));

      // binning
      const linal::Vector< double> binning_x
      (
        linal::FillVector< double>
        (
          NUMBER_FUNCTION_VALUES_X, START_X - ( CENTER_TICS_X ? 0.0 : 0.5 * DELTA_X), DELTA_X
        )
      );
      const linal::Vector< double> binning_y
      (
        linal::FillVector< double>
        (
          NUMBER_FUNCTION_VALUES_Y, START_Y - ( CENTER_TICS_Y ? 0.0 : 0.5 * DELTA_Y), DELTA_Y
        )
      );

      for( size_t count_x( 0); count_x < NUMBER_FUNCTION_VALUES_Y; ++count_x)
      {
        for( size_t count_y( 0); count_y < NUMBER_FUNCTION_VALUES_X; ++count_y)
        {
          const double val_Y( START_Y + count_x * DELTA_Y);
          const double val_X( START_X + count_y * DELTA_X);

          m_Data( count_x, count_y) = BINARY_FUNCTION( val_X, val_Y);
        }
      }

      SetTicsX( TicsFromBinning( binning_x, 1, util::Format().FFP( 1)), CENTER_TICS_X, 1);
      SetTicsY( TicsFromBinning( binning_y, 1, util::Format().FFP( 2)), CENTER_TICS_Y, 1);
    }

    //! @brief set heatmap from bicubic spline
    //! @param BICUBIC_SPLINE the bicubic spline to use
    //! @param CENTER_TICS_X center the tics x axis (false, between two bins)
    //! @param CENTER_TICS_Y center the tics y axis (false, between two bins)
    void GnuplotHeatmap::SetFromBicubicSpline
    (
      const BicubicSpline &BICUBIC_SPLINE,
      const bool CENTER_TICS_X,
      const bool CENTER_TICS_Y
    )
    {
      SetFromBinaryFunction
      (
        BICUBIC_SPLINE,
        BICUBIC_SPLINE.GetValues().GetNumberRows() * 4.0, BICUBIC_SPLINE.GetStartX(), BICUBIC_SPLINE.GetDeltaX() / 4.0,
        BICUBIC_SPLINE.GetValues().GetNumberCols() * 4.0, BICUBIC_SPLINE.GetStartY(), BICUBIC_SPLINE.GetDeltaY() / 4.0,
        CENTER_TICS_X,
        CENTER_TICS_Y
      );
    }

    //! @brief set title and labels
    //! @param TITLE title for the map
    //! @param LABEL_X label for the x axis
    //! @param LABEL_Y label for the y axis
    //! @param LABEL_Z label for the z axis
    void GnuplotHeatmap::SetTitleAndLabel
    (
      const std::string &TITLE,
      const std::string &LABEL_X,
      const std::string &LABEL_Y,
      const std::string &LABEL_Z
    )
    {
      m_Title = TITLE;
      m_LabelX = LABEL_X;
      m_LabelY = LABEL_Y;
      m_LabelZ = LABEL_Z;
    }

    //! @brief set the min and max of z
    //! @param MIN_Z minimum on z - if undefined, pick the lowest automatically
    //! @param MAX_Z maximum on z - if undefined, pick the highest automatically
    void GnuplotHeatmap::SetMinMaxZ( const double MIN_Z, const double MAX_Z)
    {
      m_MinZ = MIN_Z;
      m_MaxZ = MAX_Z;
    }

    //! @brief set tics for x axis - the two outside tics are placed at either the outer boundary or on the two outer bins
    //! @param TICS_X vector of strings
    //! @param CENTER_TICS center the tics (false, between two bins)
    //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
    //! @return true if the number of tics required matches the number of bins, so that tics are symmetric
    bool GnuplotHeatmap::SetTicsX( const storage::Vector< std::string> &TICS_X, const bool CENTER_TICS, const size_t NTH_BIN)
    {
      // empty data or tics
      if( m_Data.IsEmpty() || TICS_X.IsEmpty() || NTH_BIN == 0)
      {
        BCL_MessageCrt( "m_Data.IsEmpty() || TICS_X.IsEmpty() || NTH_BIN == 0");
        return false;
      }

      if( CENTER_TICS)
      {
        if( m_Data.GetNumberCols() - 1 != ( TICS_X.GetSize() - 1) * NTH_BIN)
        {
          BCL_MessageCrt
          (
            "m_Data.GetNumberCols() - 1 != ( TICS_X.GetSize() - 1) * NTH_BIN but center tics"
          );
          BCL_MessageCrt
          (
            "m_Data.GetNumberCols() - 1 " + util::Format()( m_Data.GetNumberCols())
          );
          BCL_MessageCrt
          (
            "( TICS_X.GetSize() - 1) " + util::Format()( ( TICS_X.GetSize() - 1))
          );
          BCL_MessageCrt
          (
            "NTH_BIN " + util::Format()( NTH_BIN)
          );
          return false;
        }
      }
      else
      {
        if( m_Data.GetNumberCols() != ( TICS_X.GetSize() - 1) * NTH_BIN)
        {
          BCL_MessageCrt
          (
            "m_Data.GetNumberCols() != ( TICS_X.GetSize() - 1) * NTH_BIN but not center tics"
          );
          BCL_MessageCrt
          (
            "m_Data.GetNumberCols() " + util::Format()( m_Data.GetNumberCols())
          );
          BCL_MessageCrt
          (
            "( TICS_X.GetSize() - 1) " + util::Format()( ( TICS_X.GetSize() - 1))
          );
          BCL_MessageCrt
          (
            "NTH_BIN " + util::Format()( NTH_BIN)
          );
          return false;
        }
      }

      m_TicsX       = TICS_X;
      m_TicsXCenter = CENTER_TICS;
      m_BinsPerTicX = NTH_BIN;

      return true;
    }

    //! @brief set tics for y axis - the two outside tics are placed at either the outer boundary or on the two outer bins
    //! @param TICS_Y vector of strings
    //! @param CENTER_TICS center the tics (false, to the top and bottom of the bin)
    //! @param NTH_BIN every nth bin a tick is added - same number of bins to left and right are required
    //! @return true if the number of tics required matches the number of bins, so that tics are symmetric
    bool GnuplotHeatmap::SetTicsY( const storage::Vector< std::string> &TICS_Y, const bool CENTER_TICS, const size_t NTH_BIN)
    {
      // empty data or tics
      if( m_Data.IsEmpty() || TICS_Y.IsEmpty() || NTH_BIN == 0)
      {
        return false;
      }

      if( CENTER_TICS)
      {
        if( m_Data.GetNumberRows() - 1 != ( TICS_Y.GetSize() - 1) * NTH_BIN)
        {
          return false;
        }
      }
      else
      {
        if( m_Data.GetNumberRows() != ( TICS_Y.GetSize() - 1) * NTH_BIN)
        {
          return false;
        }
      }

      m_TicsY       = TICS_Y;
      m_TicsYCenter = CENTER_TICS;
      m_BinsPerTicY = NTH_BIN;

      return true;
    }

    //! @brief set the pixel size and ratio
    //! @param PIXEL_X number of pixel of x
    //! @param PIXEL_Y number of pixel of y
    //! @param RATIO ratio for plot y / x - if undefined, no roatio will be set
    void GnuplotHeatmap::SetPixelAndRatio( const size_t PIXEL_X, const size_t PIXEL_Y, const double RATIO)
    {
      m_PixelX = PIXEL_X;
      m_PixelY = PIXEL_Y;
      m_Ratio  = RATIO;
    }

    //! @brief set the font
    //! @param FONT string of font name
    //! @param FONT_SIZE size of font in pixel
    void GnuplotHeatmap::SetFont( const std::string &FONT, const size_t FONT_SIZE)
    {
      m_Font = FONT;
      m_FontSize = FONT_SIZE;
    }

    //! @brief sets the coordinates for boxes
    //! @param BOX_COORDS the lower left, upper right coordinates for each box [(x,y),(x,y)]
    //! @param X_RANGE the range in the x direction of the plot - the total size of the x direction
    //! @param Y_RANGE the range in the y direction of the plot -  the total size of the y direction
    //! @param X_MIN the minimum value in the x direction of the heat map
    //! @param Y_MIN the minimum value in the y direction of the heat map
    void GnuplotHeatmap::SetBoxes
    (
      const storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > > &BOX_COORDS,
      const double X_RANGE, const double Y_RANGE, const double X_MIN, const double Y_MIN
    )
    {
      BCL_Assert( !m_Data.IsEmpty(), "set the gnuplot before setting boxes");

      storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > > scaled_coords;

      const double x_range( m_Data.GetNumberCols());
      const double y_range( m_Data.GetNumberRows());
      static const double s_x_start( -0.5);
      static const double s_y_start( -0.5);

      // iterate over the boxes to scale the coordinates into the graph
      for
      (
        storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > >::const_iterator
          box_itr( BOX_COORDS.Begin()), box_itr_end( BOX_COORDS.End());
        box_itr != box_itr_end;
        ++box_itr
      )
      {
        const double lower_left_x( box_itr->First().First());
        const double lower_left_y( box_itr->First().Second());
        const double scaled_lower_left_x( ( lower_left_x - X_MIN) / X_RANGE * x_range + s_x_start);
        const double scaled_lower_left_y( ( lower_left_y - Y_MIN) / Y_RANGE * y_range + s_y_start);
        const storage::VectorND< 2, double> scaled_lower_left( scaled_lower_left_x, scaled_lower_left_y);

        const double upper_right_x( box_itr->Second().First());
        const double upper_right_y( box_itr->Second().Second());
        const double scaled_upper_right_x( ( upper_right_x - X_MIN) / X_RANGE * x_range + s_x_start);
        const double scaled_upper_right_y( ( upper_right_y - Y_MIN) / Y_RANGE * y_range + s_y_start);
        BCL_MessageStd( "upper_right_y " + util::Format()( upper_right_y));
        const storage::VectorND< 2, double> scaled_upper_right( scaled_upper_right_x, scaled_upper_right_y);

        scaled_coords.PushBack( storage::VectorND< 2, storage::VectorND< 2, double> >( scaled_lower_left, scaled_upper_right));
      }

      m_BoxCoords = scaled_coords;
    }

    //! @brief sets the specifications for boxes
    //! @param SPECIFICATIONS strings which specify the format of the boxes
    void GnuplotHeatmap::SetBoxesSpecifications( const storage::Vector< std::string> &SPECIFICATIONS)
    {
      m_BoxCoordsSpecifications = SPECIFICATIONS;
    }

    //! @brief if set, then when WriteScript is called, the WritePreHeader string will be written. True by default
    void GnuplotHeatmap::SetWritePreHeader( const bool SET)
    {
      m_WritePreHeader = SET;
    }

    //! @brief if set, then when WriteScript is called, the WriteHeader    string will be written. True by default
    void GnuplotHeatmap::SetWriteHeader( const bool SET)
    {
      m_WriteHeader = SET;
    }

    //! @brief if set, then when WriteScript is called, the WriteBoxes     string will be written. True by default
    void GnuplotHeatmap::SetWriteBoxes( const bool SET)
    {
      m_WriteBoxes = SET;
    }

    //! @brief if set, then when WriteScript is called, the WriteData      string will be written. True by default
    void GnuplotHeatmap::SetWriteData( const bool SET)
    {
      m_WriteData = SET;
    }

    //! @brief set, the palette to use
    //! @param PALETTE the color palette used by the splot
    void GnuplotHeatmap::SetPalette( const PaletteTypes PALETTE)
    {
      m_Palette = PALETTE;
    }

    //! @brief if set, then the x tics will be shown
    void GnuplotHeatmap::SetShowXTics( const bool SET)
    {
      m_ShowXTics = SET;
    }

    //! @brief if set, then the the tics won't be mirrored on the opposite side of the axis
    void GnuplotHeatmap::SetNoMirrorTics( const bool SET)
    {
      m_NoMirrorTics = SET;
    }

    //! @brief set the ratio to which the actual plot area extends to the edge of the png
    //! @param TOP    ratio to which this plot extends to the top    edge of the png
    //! @param BOTTOM ratio to which this plot extends to the bottom edge of the png
    //! @param RIGHT  ratio to which this plot extends to the right  edge of the png
    //! @param LEFT   ratio to which this plot extends to the left   edge of the png
    void GnuplotHeatmap::SetMargins( const double TOP, const double BOTTOM, const double RIGHT, const double LEFT)
    {
      m_TopMargin    = TOP;
      m_BottomMargin = BOTTOM;
      m_RightMargin  = RIGHT;
      m_LeftMargin   = LEFT;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief multiply the heatmap
    //! @details can be used to get similar scaling between different heat maps
    //! @param FACTOR
    void GnuplotHeatmap::Muliply( const double FACTOR)
    {
      m_Data *= FACTOR;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write unchangeable information necessary for the plot
    //! @param OSTREAM stream to be written to
    //! @param TITLE title for the heat map
    std::ostream &GnuplotHeatmap::WritePreHeader( std::ostream &OSTREAM) const
    {
      // write comments
      OSTREAM << "# BCL generated heatmap\n";
      OSTREAM << "set terminal png enhanced transparent ";
      if( !m_Font.empty())
      {
        OSTREAM << "font \"" << m_Font << ',' << m_FontSize << "\" ";
      }
      OSTREAM << "size " << m_PixelX << ',' << m_PixelY << '\n';
      OSTREAM << "set output \"" << m_Filename << ".png\"\n";
      OSTREAM << "set encoding iso_8859_1\n";

      return OSTREAM;
    }

    //! @brief write a choice of palettes
    //! @param OSTREAM stream to be written to
    std::ostream &GnuplotHeatmap::WritePalettes( std::ostream &OSTREAM) const
    {
      OSTREAM << "set palette " << GetPaletteString( m_Palette) << '\n';

      // end
      return OSTREAM;
    }

    //! @brief write a header
    //! @param OSTREAM stream to be written to
    //! @param TITLE title for the heat map
    std::ostream &GnuplotHeatmap::WriteHeader( std::ostream &OSTREAM) const
    {
      OSTREAM << "set view map\n";
      if( !m_Title.empty())
      {
        OSTREAM << "set title \"" << m_Title << "\"\n";
      }
      else
      {
        OSTREAM << "unset title\n";
      }
      OSTREAM << "unset key\n\n";

      if( util::IsDefined( m_Ratio))
      {
        OSTREAM << "set size ratio " << m_Ratio << '\n';
      }

      // axes
      WriteAxisInfo( OSTREAM, coord::GetAxes().e_X, m_LabelX,   -0.5, m_Data.GetNumberCols() - 0.5, m_TicsX, m_TicsXCenter, m_BinsPerTicX);
      WriteAxisInfo( OSTREAM, coord::GetAxes().e_Y, m_LabelY,   -0.5, m_Data.GetNumberRows() - 0.5, m_TicsY, m_TicsYCenter, m_BinsPerTicY);
      WriteAxisInfo( OSTREAM, coord::GetAxes().e_Z, m_LabelZ, m_MinZ,                       m_MaxZ, storage::Vector< std::string>(), false, 0);

      // true if the color box should not be shown
      if( !m_ShowColorBox)
      {
        OSTREAM << "unset colorbox\n";
      }

      if( !m_ShowXTics)
      {
        OSTREAM << "set noxtics\n";
      }
      if( !m_NoMirrorTics)
      {
        OSTREAM << "\n";
      }

      // true if the user set the margins
      if
      (
        util::IsDefined( m_TopMargin)   && util::IsDefined( m_BottomMargin) &&
        util::IsDefined( m_RightMargin) && util::IsDefined( m_LeftMargin)
      )
      {
        // set the margins
        OSTREAM << "set lmargin screen " << m_LeftMargin   << "\n";
        OSTREAM << "set rmargin screen " << m_RightMargin  << "\n";
        OSTREAM << "set bmargin screen " << m_BottomMargin << "\n";
        OSTREAM << "set tmargin screen " << m_TopMargin    << "\n";
      }

//      OSTREAM << "set colorbox horizontal user origin 0.1,0.9 size 0.8,0.01\n";

      // end
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
    std::ostream &GnuplotHeatmap::WriteAxisInfo
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
      // string for that axis
      const std::string &axis_string( Gnuplot::GetAxisString( AXIS));

      // label
      if( !LABEL.empty())
      {
        OSTREAM << "set " << axis_string << "label \"" << LABEL << "\"";
        // Move axis label further away if tick labels are rotated
        OSTREAM << ( ( AXIS == coord::GetAxes().e_X && m_RotateXTics != 0) ? " offset 0,-1\n" : "\n");
      }

      // range
      OSTREAM << "set " << axis_string << "range [ " << ( util::IsDefined( MIN) ? util::Format()( MIN) : "*") << " : " << ( util::IsDefined( MAX) ? util::Format()( MAX) : "*") << " ]\n";

      // tics
      if( TICS.IsEmpty())
      {
        if( AXIS != coord::GetAxes().e_Z)
        {
          OSTREAM << "set no" << axis_string << "tics\n";
        }
        else
        {
          OSTREAM << "#set cbtics 1\n";
          OSTREAM << "#set format cb \"%3.1f\"\n";
        }
      }
      else
      {
        OSTREAM << "set " << axis_string << "tics ";

        if( m_NoMirrorTics)
        {
          OSTREAM << "nomirror ";
        }

        if( AXIS == coord::GetAxes().e_X)
        {
          OSTREAM << "rotate by " << m_RotateXTics << ' ';
        }

        OSTREAM << "(";
        double current_bin( 0.0);
        if( !CENTER_TICS)
        {
          current_bin -= 0.5;
        }

        // prevent endless loop
        if( NTH_BIN == 0)
        {
          return OSTREAM;
        }

        for
        (
          storage::Vector< std::string>::const_iterator itr( TICS.Begin()), itr_end( TICS.End());
          itr != itr_end;
          current_bin += NTH_BIN
        )
        {
//          OSTREAM << "\"" << *itr << "     \" " << current_bin;
          OSTREAM << "\" " << *itr << " \" " << current_bin;
          ++itr;
          if( itr == itr_end)
          {
            break;
          }
          OSTREAM << ", ";
        }
        OSTREAM << ")\n";
      }

      // end
      return OSTREAM;
    }

    //! @brief write any boxes to the script for display
    //! @param OSTREAM the stream the boxes will be written to
    std::ostream &GnuplotHeatmap::WriteBoxes( std::ostream &OSTREAM) const
    {
      static const std::string s_box_color( "black");
      static const std::string s_fill_style( "empty");
      static const size_t s_border_color( 0);
      static const double s_border_linewidth( 0.9);
      static const std::string s_border_command
      (
        "border " + util::Format()( s_border_color) + " linewidth " + util::Format()( s_border_linewidth)
      );
      static const std::string s_default_specifications
      (
        " front fillcolor rgb \"" + s_box_color + "\" fs " + s_fill_style + " " + s_border_command
      );

      size_t counter( 0);
      // iterate through the boxes
      for
      (
        storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > >::const_iterator
          box_itr( m_BoxCoords.Begin()), box_itr_end( m_BoxCoords.End());
        box_itr != box_itr_end;
        ++box_itr, ++counter
      )
      {
        const double lower_left_x( box_itr->First().First());
        const double lower_left_y( box_itr->First().Second());
        const double upper_right_x( box_itr->Second().First());
        const double upper_right_y( box_itr->Second().Second());

        std::string specifications( s_default_specifications);
        if( m_BoxCoordsSpecifications.GetSize() > counter)
        {
          specifications = m_BoxCoordsSpecifications( counter);
        }

        OSTREAM << "set object rect from " << lower_left_x << "," << lower_left_y << " to "
          << upper_right_x << "," << upper_right_y
          << " " << specifications << "\n";
      }

      return OSTREAM;
    }

    //! @brief write the data
    //! @param OSTREAM the stream to write to
    //! @return the stream written to
    std::ostream &GnuplotHeatmap::WriteData( std::ostream &OSTREAM) const
    {
      OSTREAM << "splot '-' using 1:2:3 with image\n";
      OSTREAM << "# number x values " << m_Data.GetNumberCols() << '\n';
      OSTREAM << "# number y values " << m_Data.GetNumberRows() << '\n';

      // ptr to first element
      const double *ptr( m_Data.Begin());

      // iterate over rows
      for( size_t row( 0); row != m_Data.GetNumberRows(); ++row)
      {
        for( size_t col( 0); col != m_Data.GetNumberCols(); ++col, ++ptr)
        {
          OSTREAM << col << '\t' << row << " \t" << *ptr << '\n';
        }
        OSTREAM << '\n';
      }

      // indicate end of data for heatmap
      OSTREAM << "e\n";

      // end
      return OSTREAM;
    }

    //! @brief write the gnuplot script
    //! @param OSTREAM the stream to write to
    void GnuplotHeatmap::WriteScript( std::ostream &OSTREAM) const
    {
      if( m_WritePreHeader)
      {
        WritePreHeader( OSTREAM);
      }
      if( m_WriteHeader)
      {
        WriteHeader(    OSTREAM);
      }
      if( m_Palette != e_None)
      {
        WritePalettes( OSTREAM);
      }
      if( m_WriteBoxes)
      {
        WriteBoxes(     OSTREAM);
      }
      if( m_WriteData)
      {
        WriteData(      OSTREAM);
      }
    }

  } // namespace math
} // namespace bcl
