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
#include "restraint/bcl_restraint_analyze_coordinate_distance_distribution.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_gnuplot_multiplot.h"
#include "math/bcl_math_histogram.h"
#include "math/bcl_math_running_average_sd.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeCoordinateDistanceDistribution::s_Instance
    (
      util::Enumerated< assemble::AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeCoordinateDistanceDistribution())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeCoordinateDistanceDistribution::AnalyzeCoordinateDistanceDistribution() :
      m_OutFilePostFix( ".CoordinateDistanceDistribution"),
      m_CoordinateA( assemble::LocatorAtom( 'A', 1, biol::GetAtomTypes().CA)),
      m_CoordinateB( assemble::LocatorAtom( 'A', 1, biol::GetAtomTypes().CA)),
      m_ComparisonFunction(),
      m_HistogramMinimum( 0),
      m_HistogramBinSize( 2),
      m_HistogramNumberOfBins( 30),
      m_Title( "Title"),
      m_PixelX   ( 600),
      m_PixelY   ( 400),
      m_Font     ( "Arial"),
      m_FontSize ( 12),
      m_GreyScale( false),
      m_PlotRatio( util::GetUndefinedDouble()),
      m_MeanStdDevOutFile( "mean_stddev.txt"),
      m_BinsPerTic( 1),
      m_CenterTics( false),
      m_MinZ( util::GetUndefinedDouble()),
      m_MaxZ( util::GetUndefinedDouble()),
      m_Normalize( true)
      {
      }

    //! @brief Clone function
    //! @return pointer to new AnalyzeCoordinateDistanceDistribution
    AnalyzeCoordinateDistanceDistribution *AnalyzeCoordinateDistanceDistribution::Clone() const
    {
      return new AnalyzeCoordinateDistanceDistribution( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeCoordinateDistanceDistribution::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeCoordinateDistanceDistribution::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeCoordinateDistanceDistribution::GetAlias() const
    {
      static const std::string s_Name( "CoordinateDistanceDistribution");
      return s_Name;
    }

  ////////////////
  // operations //
  ////////////////

  ///////////////
  // operators //
  ///////////////

    //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
    //! @param ENSEMBLE the ensemble that will be analyzed
    //! @return string which has the analyzed information about the ensemble
    std::string AnalyzeCoordinateDistanceDistribution::operator()( const assemble::ProteinEnsemble &ENSEMBLE) const
    {
      // intialize
      math::GnuplotMultiplot multiplot;
      multiplot.SetRowsCols( 2, 1);
      multiplot.SetFont( m_Font, m_FontSize);
      multiplot.SetPixelAndRatio( m_PixelX, m_PixelY, util::GetUndefinedDouble());

      // ensemble histogram plot
      multiplot.Insert( GetHeatMap( GetDistanceHistogram( ENSEMBLE)));

      // comparison function histogram plot
      multiplot.Insert( GetReferenceHeatMap());

      // write heat map to string stream
      std::stringstream stream;
      multiplot.WriteScript( stream);
      std::string analysis( stream.str());

      // return analysis string of the gnuplot script
      return analysis;
    }

  //////////////////////
  // input and output //
  //////////////////////

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeCoordinateDistanceDistribution::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Creates heat map of showing frequency with which a distance between two coordinates is observed. The heat "
        "map of frequency with which a distance is observed for two coordinates located in a protein ensemble can "
        "be compared with a second heat map. The second heat map is given by a math function. This allows the "
        "distribution of distances coming from the ensemble to be compared with a second distribution which is "
        "perhaps the expected distribution. An example use of this is to see the distribution of distances coming "
        "from two residues in the ensemble whose distance has been measured by EPR. This ensemble distribution of "
        "distances can then be compared to the EPR distribution described by a gaussian function defined by the "
        "mean and standard deviation of the distance probability distribution measured by EPR."
      );

      parameters.AddInitializer
      (
        "filename_postfix",
        "the postfix that will be appended to the output file in order to identify this analysis",
        io::Serialization::GetAgent( &m_OutFilePostFix),
        ".ChiAnglePairDistribution"
      );

      parameters.AddInitializer
      (
        "coord_a",
        "the first coordinate that will be used in the distance",
        io::Serialization::GetAgent( &m_CoordinateA)
      );

      parameters.AddInitializer
      (
        "coord_b",
        "the second coordinate that will be used in the distance",
        io::Serialization::GetAgent( &m_CoordinateB)
      );

      parameters.AddInitializer
      (
        "function",
        "the function that the coordinate heat map will be compared against",
        io::Serialization::GetAgent( &m_ComparisonFunction)
      );

      parameters.AddInitializer
      (
        "histogram_minimum",
        "the minimal value representing the left boundary of the score histogram",
        io::Serialization::GetAgent( &m_HistogramMinimum),
        "-1.0"
      );

      parameters.AddInitializer
      (
        "histogram_binsize",
        "the width of one bin of the score histograms",
        io::Serialization::GetAgent( &m_HistogramBinSize),
        "0.1"
      );

      parameters.AddInitializer
      (
        "histogram_num_bins",
        "the number of bins in the score histograms",
        io::Serialization::GetAgent( &m_HistogramNumberOfBins),
        "10"
      );

      parameters.AddInitializer
      (
        "title",
        "the title that will label the resulting heat map",
        io::Serialization::GetAgent( &m_Title),
        "Title"
      );

      parameters.AddInitializer
      (
        "pixel_x",
        "The size of the plot png in the x direction",
        io::Serialization::GetAgent( &m_PixelX),
        "600"
      );

      parameters.AddInitializer
      (
        "pixel_y",
        "The size of the plot png in the y direction",
        io::Serialization::GetAgent( &m_PixelY),
        "400"
      );

      parameters.AddInitializer
      (
        "font",
        "The font to be used for text in the plot",
        io::Serialization::GetAgent( &m_Font),
        "Arial"
      );

      parameters.AddInitializer
      (
        "font_size",
        "The size of the font for the plot",
        io::Serialization::GetAgent( &m_FontSize),
        "12"
      );

      parameters.AddInitializer
      (
        "grey_scale",
        "boolean true the color palette for the gradient should be grey scale. 1=true;0=false",
        io::Serialization::GetAgent( &m_GreyScale),
        "0"
      );

      parameters.AddInitializer
      (
        "plot_ratio",
        "How much of the png area the plot will cover top to bottom and left to write."
        " (this does not include any labels or tics). So if you say 0.5, plot will go from 0 to 0.5 from bottom"
        " to top and left to right",
        io::Serialization::GetAgent( &m_PlotRatio),
        util::Format()( util::GetUndefinedDouble())
      );

      parameters.AddInitializer
      (
        "mean_stddev_out_file",
        "Full path and filename for holding the mean and standard deviation of the distribution",
        io::Serialization::GetAgent( &m_MeanStdDevOutFile),
        "mean_stddev.txt"
      );

      parameters.AddInitializer
      (
        "bins_per_tic",
        "how many bins per tic in the plot. 1 means every bin is labeled, 2 means every other, etc.",
        io::Serialization::GetAgent( &m_BinsPerTic),
        "1"
      );

      parameters.AddInitializer
      (
        "center_tics",
        "boolean if true, tics will be centered on the bins, if false, will be on edges of bins. 1=true;0=false",
        io::Serialization::GetAgent( &m_CenterTics),
        "0"
      );

      parameters.AddInitializer
      (
        "min_z",
        "the minimum value used in pymol for the z axis (i.e. color). Values below this will just get the minimum value color.",
        io::Serialization::GetAgent( &m_MinZ),
        util::Format()( util::GetUndefinedDouble())
      );
      parameters.AddInitializer
      (
        "max_z",
        "the maximum value used in pymol for the z axis (i.e. color). Values above this will just get the maximum value color.",
        io::Serialization::GetAgent( &m_MaxZ),
        util::Format()( util::GetUndefinedDouble())
      );

      parameters.AddInitializer
      (
        "normalize_histogram",
        "boolean true indicates that the distance distribution histogram will be normalized. 1=true;0=false",
        io::Serialization::GetAgent( &m_Normalize),
        "1"
      );

      return parameters;
    }

    //! @brief gives the histogram of the distance distribution of the protein ensemble
    //! @param ENSEMBLE the ensemble for which the distance distribution will be calculated
    //! @return histogram that has the distance distribution calculated from the provided protein ensemble
    math::Histogram AnalyzeCoordinateDistanceDistribution::GetDistanceHistogram
    (
      const assemble::ProteinEnsemble &ENSEMBLE
    ) const
    {
      // initialize the histogram distance distribution
      math::Histogram distance_distribution( m_HistogramMinimum, m_HistogramBinSize, m_HistogramNumberOfBins);

      // to keep track of the mean and standard deviation of the distances in the ensemble
      math::RunningAverageSD< double> mean_stddev;

      // iterate through the ensemble to get the distribution of distances
      for
      (
        assemble::ProteinEnsemble::const_iterator
          ensemble_itr( ENSEMBLE.Begin()), ensemble_itr_end( ENSEMBLE.End());
        ensemble_itr != ensemble_itr_end;
        ++ensemble_itr
      )
      {
        // locate the two coordinates and make sure they are defined
        const linal::Vector3D coord_a( m_CoordinateA->Locate( **ensemble_itr));
        BCL_Assert( coord_a.IsDefined(), "coord_a could not be located with " + util::Format()( *m_CoordinateA));
        const linal::Vector3D coord_b( m_CoordinateB->Locate( **ensemble_itr));
        BCL_Assert( coord_b.IsDefined(), "coord_b could not be located with " + util::Format()( *m_CoordinateB));

        // calculate the distance between the two coordinates and make sure it is defined
        const double distance( linal::Distance( coord_a, coord_b));
        BCL_Assert( util::IsDefined( distance), "distance is not defined");

        // add the distance into the distance distribution and the statistics object
        distance_distribution.PushBack( distance);
        mean_stddev += distance;
      }

      // write the mean and standard deviation of the observed distances to the desired output file
      io::OFStream write;
      io::File::MustOpenOFStream( write, m_MeanStdDevOutFile);
      write << m_Title << "\tmean " << mean_stddev.GetAverage() << "\tstddev " << mean_stddev.GetStandardDeviation();
      io::File::CloseClearFStream( write);

      // true if desired to normalize
      if( m_Normalize)
      {
        // normalize distribution
        distance_distribution.Normalize();
      }

      // return the histogram of the distances
      return distance_distribution;
    }

    //! @brief gives the heat map corresponding to the provided histogram made from the protein ensemble
    //! @param HISTGRAM the histogram for which the heat map will be made
    //! @return heat map representing the distance distribution of the protein ensemble
    util::ShPtr< math::GnuplotHeatmap>
    AnalyzeCoordinateDistanceDistribution::GetHeatMap( const math::Histogram &HISTOGRAM) const
    {
      // make new heatmap
      util::ShPtr< math::GnuplotHeatmap> gnuplot( new math::GnuplotHeatmap());

      // set the heatmap from the provided histogram
      gnuplot->SetFromHistogram( HISTOGRAM, false, false);

      // set options on the heat map
      gnuplot->SetShowColorBox( false);
      gnuplot->SetTitleAndLabel( "", "", "", "");
      gnuplot->SetMargins( 0.5, 0.5 - m_PlotRatio, 0.99, 0.01);
      gnuplot->SetWritePreHeader( false);
      gnuplot->SetWriteHeader( true);
      gnuplot->SetPalette( math::GnuplotHeatmap::e_GreyScale);
      gnuplot->SetRotationXTics( 90);
      gnuplot->SetNoMirrorTics( true);

      // set the tics of the heat map according to the preferences of the user
      {
        const bool center_tic_x( m_CenterTics);
        const linal::Vector< double> binning_x
        (
          linal::FillVector< double>
          (
            HISTOGRAM.GetNumberOfBins() + size_t( !center_tic_x),
            HISTOGRAM.GetBoundaries().First() + ( center_tic_x ? 0.5 * HISTOGRAM.GetBinSize() : 0.0),
            HISTOGRAM.GetBinSize()
          )
        );
        BCL_Assert
        (
          gnuplot->SetTicsX
          (
            math::GnuplotHeatmap::TicsFromBinning( binning_x, m_BinsPerTic, util::Format().W( 4)),
            center_tic_x,
            m_BinsPerTic
          ), "unable to set tics x"
        );
      }

      // set the min max z values
      gnuplot->SetMinMaxZ( m_MinZ, m_MaxZ);

      // return the heat map
      return gnuplot;
    }

    //! @brief gives the heat map corresponding to the reference function
    //! @return shptr to a heat map that represents the reference function
    util::ShPtr< math::GnuplotHeatmap> AnalyzeCoordinateDistanceDistribution::GetReferenceHeatMap() const
    {
      // make new heatmap
      util::ShPtr< math::GnuplotHeatmap> gnuplot( new math::GnuplotHeatmap());

      // set the heat map from the reference function
      gnuplot->SetFromFunction
      (
        *m_ComparisonFunction, m_HistogramNumberOfBins,
        m_HistogramMinimum + ( m_HistogramBinSize / 2.0), m_HistogramBinSize, false, false, m_Normalize
      );

      // set options on the heat map
      gnuplot->SetMargins( 0.5 + m_PlotRatio, 0.5, 0.99, 0.01);
      gnuplot->SetNoMirrorTics( true);
      gnuplot->SetWritePreHeader( false);
      gnuplot->SetWriteHeader( true);
      gnuplot->SetTitleAndLabel( m_Title, "", "", "");
      gnuplot->SetPalette( math::GnuplotHeatmap::e_GreyScale);
      gnuplot->SetShowXTics( false);
      // set the min max z values
      gnuplot->SetMinMaxZ( m_MinZ, m_MaxZ);

      // return the heat map
      return gnuplot;
    }
  } // namespace restraint
} // namespace bcl
