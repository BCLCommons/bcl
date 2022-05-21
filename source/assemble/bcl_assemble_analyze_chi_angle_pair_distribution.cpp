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
#include "assemble/bcl_assemble_analyze_chi_angle_pair_distribution.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "biol/bcl_biol_rotamer.h"
#include "io/bcl_io_file.h"
#include "io/bcl_io_serialization.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_histogram_2d.h"
#include "util/bcl_util_string_functions.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> AnalyzeChiAnglePairDistribution::s_Instance
    (
      util::Enumerated< AnalyzeProteinEnsembleInterface>::AddInstance( new AnalyzeChiAnglePairDistribution())
    );

    //! @brief provides the min and max values contained in the heat map
    //! @return min and max values contained in the heat map
    const storage::VectorND< 2, double> &AnalyzeChiAnglePairDistribution::GetMinMax()
    {
      static const storage::VectorND< 2, double> s_min_max( -180.0, 180.0);
      return s_min_max;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    AnalyzeChiAnglePairDistribution::AnalyzeChiAnglePairDistribution() :
      m_OutFilePostFix( ".ChiAnglePairDistribution"),
      m_HistogramBinSize( 60),
      m_LocatorAA( CollectorAAType( storage::Set< biol::AAType>::Create( biol::GetAATypes().ALA, biol::GetAATypes().ARG))),
      m_ChiA( biol::ChiAngle::e_One),
      m_ChiB( biol::ChiAngle::e_Two),
      m_SetTitleAndLabel( true),
      m_Title( "Title"),
      m_ShowColorBox( true),
      m_FontSize( 10),
      m_ResolutionX( 640),
      m_ResolutionY( 640),
      m_ReferenceChiFilename( "native_chi_pairs.ls"),
      m_Font     ( "Arial"),
      m_GreyScale( false),
      m_PlotRatio( util::GetUndefinedDouble()),
      m_BinsPerTic( 1),
      m_CenterTics( false),
      m_MinZ( util::GetUndefinedDouble()),
      m_MaxZ( util::GetUndefinedDouble()),
      m_Normalize( true)
    {
    }

    //! @brief Clone function
    //! @return pointer to new AnalyzeChiAnglePairDistribution
    AnalyzeChiAnglePairDistribution *AnalyzeChiAnglePairDistribution::Clone() const
    {
      return new AnalyzeChiAnglePairDistribution( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &AnalyzeChiAnglePairDistribution::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gives the string to append to the the end of a filename to identify this analysis
    //! @return gives the string to append to the the end of a filename to identify this analysis
    const std::string &AnalyzeChiAnglePairDistribution::GetOutFilePostfix() const
    {
      return m_OutFilePostFix;
    }

    //! @brief returns the name used for this class in an object data label
    //! @return the name used for this class in an object data label
    const std::string &AnalyzeChiAnglePairDistribution::GetAlias() const
    {
      static const std::string s_Name( "ChiAnglePairDistribution");
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
    std::string AnalyzeChiAnglePairDistribution::operator()( const ProteinEnsemble &ENSEMBLE) const
    {
      // calculate normalized chi angle statistics
      const math::Histogram2D chi_stats
      (
        CalculateChiAngleStatistics( ENSEMBLE, m_HistogramBinSize, m_LocatorAA, m_ChiA, m_ChiB)
      );

      // get the gnuplot script object
      const math::GnuplotHeatmap heatmap
      (
        GetHeatMap
        (
          chi_stats, m_ChiA, m_ChiB, m_SetTitleAndLabel, m_Title, m_ShowColorBox,
          m_FontSize, m_ResolutionX, m_ResolutionY, m_ReferenceChiFilename
        )
      );

      // write heat map to string stream
      std::stringstream stream;
      heatmap.WriteScript( stream);

      // return the string
      return stream.str();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &AnalyzeChiAnglePairDistribution::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_OutFilePostFix, ISTREAM);
      io::Serialize::Read( m_HistogramBinSize, ISTREAM);
      io::Serialize::Read( m_LocatorAA, ISTREAM);
      io::Serialize::Read( m_ChiA, ISTREAM);
      io::Serialize::Read( m_ChiB, ISTREAM);

      io::Serialize::Read( m_SetTitleAndLabel, ISTREAM);
      io::Serialize::Read( m_Title, ISTREAM);
      io::Serialize::Read( m_ShowColorBox, ISTREAM);
      io::Serialize::Read( m_FontSize, ISTREAM);
      io::Serialize::Read( m_ResolutionX, ISTREAM);
      io::Serialize::Read( m_ResolutionY, ISTREAM);
      io::Serialize::Read( m_ReferenceChiFilename, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &AnalyzeChiAnglePairDistribution::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_OutFilePostFix, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_HistogramBinSize, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_LocatorAA, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_ChiA, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_ChiB, OSTREAM, INDENT) << "\n";

      io::Serialize::Write( m_SetTitleAndLabel, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_Title, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_ShowColorBox, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_FontSize, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_ResolutionX, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_ResolutionY, OSTREAM, INDENT) << "\n";
      io::Serialize::Write( m_ReferenceChiFilename, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AnalyzeChiAnglePairDistribution::GetSerializer() const
    {
      io::Serializer parameters;
      parameters.SetClassDescription
      (
        "Outputs gnuplot script to create 2D heat map showing frequency with which chi angles are observed."
        "Uses a collector interface to find the residues of interest from the protein ensemble. The chi angles"
        "are calculated for two desired chi angles for the residues of interest. The frequency with which"
        "a pair of chi angles is observed is shown in the heat map. The possibility to highlight up to two"
        "different sets of desired chi angle combinations is available (e.g. a) combinations available in a"
        "rotamer library and b) experimentally observed combinations). These are highlighted via boxes of"
        "specified dimensions."
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
        "histogram_binsize",
        "the width of one bin of the score histograms",
        io::Serialization::GetAgent( &m_HistogramBinSize),
        "0.1"
      );

      parameters.AddInitializer
      (
        "collector_type",
        "the type of collector that should be used to get residues of interest",
        io::Serialization::GetAgent( &m_LocatorAA)
      );

      parameters.AddInitializer
      (
        "chi_angle_a",
        "the first of two chi angles that statistics should be done over",
        io::Serialization::GetAgent( &m_ChiA),
        "e_One"
      );

      parameters.AddInitializer
      (
        "chi_angle_b",
        "the second of two chi angles that statistics should be done over",
        io::Serialization::GetAgent( &m_ChiB),
        "e_Two"
      );

      parameters.AddInitializer
      (
        "set_title_and_label",
        "1 if the title and axis labels should be set on the resulting heat map - 0 otherwise",
        io::Serialization::GetAgent( &m_SetTitleAndLabel),
        "1"
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
        "show_color_box",
        "1 if the color box (i.e. coloring scale) should be shown - 0 otherwise",
        io::Serialization::GetAgent( &m_ShowColorBox),
        "1"
      );

      parameters.AddInitializer
      (
        "font_size",
        "the size text in the heat map should have",
        io::Serialization::GetAgent( &m_FontSize),
        "10"
      );

      parameters.AddInitializer
      (
        "x_pixels",
        "the number of pixels in the x direction",
        io::Serialization::GetAgent( &m_ResolutionX),
        "640"
      );

      parameters.AddInitializer
      (
        "y_pixels",
        "the number of pixels in the y direction",
        io::Serialization::GetAgent( &m_ResolutionY),
        "640"
      );

      parameters.AddInitializer
      (
        "reference_chi_filename",
        "The name of the file that contains reference chi angle pairs and optional gnuplot formatted specifications"
        "for the box that will indicate the location of the chi angle pair. The format of the input file should be"
        "\n<chi_a> <chi_b> <box_side_length> <specifications>\nAn example is\n"
        "-69.01 -56.55 60 fs empty border -1 linewidth 0.5 front\nor\n"
        "-76.44 -55.52 30 fs empty border  0 linewidth 0.5 front\n"
        "this file is always attempted to be read, but if it cannot be opened, the heatmap is made anyways.",
        io::Serialization::GetAgent( &m_ReferenceChiFilename),
        "reference_chi_angle_pairs.ls"
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

    //! @brief gives distribution for observed combinations of two chi angles for a desired residue in an ensemble
    //! @param ENSEMBLE the ensemble the distribution will be calculated over
    //! @param BIN_SIZE the size of bins that should be used for the histogram distribution
    //! @param LOCATOR the method for locating the residue of interest
    //! @param CHI_A the first chi angle of interest
    //! @param CHI_B the second chi angle of interest
    //! @return histogram 2D which has the distribution of chia and chib angles observed in the ensemble
    math::Histogram2D AnalyzeChiAnglePairDistribution::CalculateChiAngleStatistics
    (
      const ProteinEnsemble &ENSEMBLE, const double BIN_SIZE,
      const find::CollectorInterface
      <
        util::SiPtrList< const biol::AABase>, util::SiPtrVector< const biol::AABase>
       > &LOCATOR,
       const biol::ChiAngle::ChiEnum CHI_A, const biol::ChiAngle::ChiEnum CHI_B
    ) const
    {
      // to hold the distribution
      const size_t num_bins( ( GetMinMax().Second() - GetMinMax().First()) / BIN_SIZE);
      math::Histogram2D distribution
      (
        storage::VectorND< 2, double>( GetMinMax().First(), GetMinMax().First()),
        storage::VectorND< 2, double>( BIN_SIZE, BIN_SIZE), storage::VectorND< 2, size_t>( num_bins, num_bins)
      );

      // iterate through the ensemble
      for
      (
        ProteinEnsemble::const_iterator ensemble_itr( ENSEMBLE.Begin()), ensemble_itr_end( ENSEMBLE.End());
        ensemble_itr != ensemble_itr_end;
        ++ensemble_itr
      )
      {
        // locate the residue of interest
        const util::SiPtrList< const biol::AABase> resi( LOCATOR.Collect( ( *ensemble_itr)->GetAminoAcids()));

        // iterate over the located residues
        for
        (
          util::SiPtrList< const biol::AABase>::const_iterator
            resi_itr( resi.Begin()), resi_itr_end( resi.End());
          resi_itr != resi_itr_end;
          ++resi_itr
        )
        {
          // true if the resi could not be located
          if( !resi_itr->IsDefined())
          {
            BCL_MessageDbg( "resi " + util::Format()( LOCATOR) + " could not be located");

            // go to next model in ensemble
            continue;
          }

          // the dihedral angles
          const biol::Rotamer dihedral_angles( ( *resi_itr)->CalculateSideChainDihedralAngles());

          // the chi angles of interest
          const double chi_angle_a( dihedral_angles.GetAngle( CHI_A, math::Angle::e_Degree)),
            chi_angle_b( dihedral_angles.GetAngle( CHI_B, math::Angle::e_Degree));

          // true if the desired chi angles of the residue could not be calculated
          if( !util::IsDefined( chi_angle_a) || !util::IsDefined( chi_angle_b))
          {
            BCL_MessageDbg
            (
              "can not get dihedral angles chi " + CHI_A.GetString()
              + " and chi " + CHI_B.GetString() +
              " of residue " + ( *resi_itr)->GetIdentification() + " because only the \n" +
              util::Format()( dihedral_angles) + "\ndihedral angles could be calculated"
            );

            // go to next in ensemble
            continue;
          }

          distribution.PushBack( storage::VectorND< 2, double>( chi_angle_a, chi_angle_b));
        }

      }

      // true if desired to normalize
      if( m_Normalize)
      {
        // normalize distribution
        distribution.Normalize();
      }

      return distribution;
    }

    //! @brief creats a heat map based on the distribution of chi angles
    //! @param DISTRIBUTION the distribution of chi angles that will be turned into a heat map
    //! @param CHI_A the first chi angle of interest
    //! @param CHI_B the second chi angle of interest
    //! @param SET_TITLE_AND_LABEL true if the title and axis labels should be set on the resulting heat map or not
    //! @param TITLE the title that will label the resulting heat map
    //! @param SHOW_COLOR_BOX bool TRUE if the color box should be displayed or not - box that shows the color scale
    //! @param FONT_SIZE the size text in the heat map should have
    //! @param RESOLUTION_X the pixel height and width of the heat map graphic
    //! @param RESOLUTION_Y the pixel width and width of the heat map graphic
    //! @param CHI_FILENAME the name of the file that contains chi angle pairs of interest
    //! @return GnuplotHeatmap showing the distribution of chi angles
    math::GnuplotHeatmap AnalyzeChiAnglePairDistribution::GetHeatMap
    (
      const math::Histogram2D &DISTRIBUTION, const biol::ChiAngle::ChiEnum CHI_A,
      const biol::ChiAngle::ChiEnum CHI_B, const bool SET_TITLE_AND_LABEL, const std::string &TITLE,
      const bool SHOW_COLOR_BOX, const double FONT_SIZE, const size_t RESOLUTION_X, const size_t RESOLUTION_Y,
      const std::string &CHI_FILENAME
    ) const
    {
      math::GnuplotHeatmap heatmap;
      const bool center_tic_x( m_CenterTics);
      const bool center_tic_y( m_CenterTics);
      BCL_MessageStd( "center tic x is " + util::Format()( center_tic_x));
      BCL_MessageStd( "center tic y is " + util::Format()( center_tic_y));
      heatmap.SetFromHistogram( DISTRIBUTION, center_tic_x, center_tic_y);

      if( SET_TITLE_AND_LABEL)
      {
        heatmap.SetTitleAndLabel
        (
          TITLE,
          "X_" + util::Format()( CHI_A + 1) + " (\\260)",
          "X_" + util::Format()( CHI_B + 1) + " (\\260)",
          "Frequency (fraction of counts)"
        );
      }
      heatmap.SetRotationXTics( 90);
      heatmap.SetFilename( "AnalyzeChiAnglePairDistribution.gnuplot");
      heatmap.SetFont( m_Font, FONT_SIZE);
      heatmap.SetShowColorBox( SHOW_COLOR_BOX);
      heatmap.SetPixelAndRatio( RESOLUTION_X, RESOLUTION_Y, 1.0);
      heatmap.SetPalette( math::GnuplotHeatmap::e_GreyScale);
      heatmap.SetNoMirrorTics( true);

      if( util::IsDefined( m_PlotRatio) && m_PlotRatio != double( 0))
      {
        heatmap.SetMargins( m_PlotRatio, 0, m_PlotRatio, 0);
      }

      // binning tics
      BCL_MessageCrt( "manually setting tics");
      const linal::Vector< double> binning_x
      (
        linal::FillVector< double>
        (
          DISTRIBUTION.GetNumberOfBinsX() + size_t( !center_tic_x),
          DISTRIBUTION.GetBoundariesX().First() + ( center_tic_x ? 0.5 * DISTRIBUTION.GetBinSizeXY().First() : 0.0),
          DISTRIBUTION.GetBinSizeXY().First()
        )
      );
      BCL_Assert
      (
        heatmap.SetTicsX
        (
          math::GnuplotHeatmap::TicsFromBinning( binning_x, m_BinsPerTic, util::Format().W( 4)),
          center_tic_x,
          m_BinsPerTic
        ), "unable to set tics x"
      );
      const linal::Vector< double> binning_y
      (
        linal::FillVector< double>
        (
          DISTRIBUTION.GetNumberOfBinsY() + size_t( !center_tic_y),
          DISTRIBUTION.GetBoundariesY().First() + ( center_tic_y ? 0.5 * DISTRIBUTION.GetBinSizeXY().Second() : 0.0),
          DISTRIBUTION.GetBinSizeXY().Second()
        )
      );
      BCL_Assert
      (
        heatmap.SetTicsY
        (
          math::GnuplotHeatmap::TicsFromBinning( binning_y, m_BinsPerTic, util::Format().W( 4)),
          center_tic_y,
          m_BinsPerTic
        ), "unable to set tics y"
      );

      // set the min max z values
      heatmap.SetMinMaxZ( m_MinZ, m_MaxZ);

      const double x_min  ( DISTRIBUTION.GetBoundariesX().First());
      const double x_range( DISTRIBUTION.GetBoundariesX().Second() - DISTRIBUTION.GetBoundariesX().First());
      const double y_min  ( DISTRIBUTION.GetBoundariesY().First());
      const double y_range( DISTRIBUTION.GetBoundariesY().Second() - DISTRIBUTION.GetBoundariesY().First());
      AddBoxesToHeatMap( heatmap, CHI_FILENAME, x_range, y_range, x_min, y_min);

      return heatmap;
    }

    //! @brief provides the text necessary to add any desired boxes to the heat map
    //! @param HEAT_MAP the heat map the boxes will be added to
    //! @param CHI_FILENAME the name of the file that contains chi angle pairs of interest
    //! @param X_RANGE the range in the x direction of the plot - the total size of the x direction
    //! @param Y_RANGE the range in the y direction of the plot -  the total size of the y direction
    //! @param X_MIN the minimum value in the x direction of the heat map
    //! @param Y_MIN the minimum value in the y direction of the heat map
    void
    AnalyzeChiAnglePairDistribution::AddBoxesToHeatMap
    (
      math::GnuplotHeatmap &HEAT_MAP, const std::string &CHI_FILENAME,
      const double X_RANGE, const double Y_RANGE, const double X_MIN, const double Y_MIN
    )
    {
      // try to open the chi containing file
      io::IFStream read;
      const bool opened( io::File::TryOpenIFStream( read, CHI_FILENAME));
      if( !opened)
      {
        return;
      }

      // read in teh chi info
      storage::Vector< storage::Vector< std::string> > chi_box_info( util::SplittedStringLineListFromIStream( read));

      // to hold the box coordinates and any specifications of boxes
      storage::Vector< storage::VectorND< 2, storage::VectorND< 2, double> > > box_coords;
      storage::Vector< std::string> box_specifications;

      // iterate through the split lines to get the chi angle and box information
      for
      (
        storage::Vector< storage::Vector< std::string> >::const_iterator
          line_itr( chi_box_info.Begin()), line_itr_end( chi_box_info.End());
        line_itr != line_itr_end;
        ++line_itr
      )
      {
        const storage::Vector< std::string> &current_line( *line_itr);

        BCL_Assert
        (
          current_line.GetSize() > 2, "need to have at least three columns in chi file. But current line has " +
          util::Format()( current_line) + " See parameter reference_chi_filename for information."
        );

        const double chi_a    ( util::ConvertStringToNumericalValue< double>( current_line( 0)));
        const double chi_b    ( util::ConvertStringToNumericalValue< double>( current_line( 1)));
        const double side_size( util::ConvertStringToNumericalValue< double>( current_line( 2)));
        const double half_side_size( side_size / 2.0);

        // get box corner coordinates
        const double lower_x( chi_a - half_side_size);
        const double lower_y( chi_b - half_side_size);
        const double upper_x( chi_a + half_side_size);
        const double upper_y( chi_b + half_side_size);

        // get the specifications
        std::string current_specs;
        for
        (
          storage::Vector< std::string>::const_iterator
            line_itr_b( line_itr->Begin() + 3), line_itr_end_b( line_itr->End());
          line_itr_b != line_itr_end_b;
          ++line_itr_b
        )
        {
          current_specs += " " + *line_itr_b;
        }

        // add box
        {
          storage::VectorND< 2, storage::VectorND< 2, double> > box
          (
            storage::VectorND< 2, double>( lower_x, lower_y), storage::VectorND< 2, double>( upper_x, upper_y)
          );
          box_coords.PushBack( box);
          box_specifications.PushBack( current_specs);
        }

        // get parts of box that need to wrap around to other side of plot
        // x direction - wrap from left to right
        if( math::Absolute( lower_x) > math::Absolute( GetMinMax().First()))
        {
          const double remainder( math::Absolute( lower_x) - math::Absolute( GetMinMax().First()));
          const double remainder_lower_x( GetMinMax().Second() - remainder);
          const double remainder_upper_x( GetMinMax().Second());
          storage::VectorND< 2, storage::VectorND< 2, double> > box
          (
            storage::VectorND< 2, double>( remainder_lower_x, lower_y),
            storage::VectorND< 2, double>( remainder_upper_x, upper_y)
          );
          box_coords.PushBack( box);
          box_specifications.PushBack( current_specs);
        }
        // x direction - wrap from right to left
        if( math::Absolute( upper_x) > math::Absolute( GetMinMax().Second()))
        {
          const double remainder( math::Absolute( upper_x) - math::Absolute( GetMinMax().Second()));
          const double remainder_lower_x( GetMinMax().First());
          const double remainder_upper_x( GetMinMax().First() + remainder);
          storage::VectorND< 2, storage::VectorND< 2, double> > box
          (
            storage::VectorND< 2, double>( remainder_lower_x, lower_y),
            storage::VectorND< 2, double>( remainder_upper_x, upper_y)
          );
          box_coords.PushBack( box);
          box_specifications.PushBack( current_specs);
        }
        // y direction - wrap from bottom to top
        if( math::Absolute( lower_y) > math::Absolute( GetMinMax().First()))
        {
          const double remainder( math::Absolute( lower_y) - math::Absolute( GetMinMax().First()));
          const double remainder_lower_y( GetMinMax().Second() - remainder);
          const double remainder_upper_y( GetMinMax().Second());
          storage::VectorND< 2, storage::VectorND< 2, double> > box
          (
            storage::VectorND< 2, double>( lower_x, remainder_lower_y),
            storage::VectorND< 2, double>( upper_x, remainder_upper_y)
          );
          box_coords.PushBack( box);
          box_specifications.PushBack( current_specs);
        }

        // y direction - wrap from top to bottom
        if( math::Absolute( upper_y) > math::Absolute( GetMinMax().Second()))
        {
          const double remainder( math::Absolute( upper_y) - math::Absolute( GetMinMax().Second()));
          const double remainder_lower_y( GetMinMax().First());
          const double remainder_upper_y( GetMinMax().First() + remainder);
          storage::VectorND< 2, storage::VectorND< 2, double> > box
          (
            storage::VectorND< 2, double>( lower_x, remainder_lower_y),
            storage::VectorND< 2, double>( upper_x, remainder_upper_y)
          );
          box_coords.PushBack( box);
          box_specifications.PushBack( current_specs);
        }

        // lower left corner - wrap to upper right corner
        if
        (
          math::Absolute( lower_x) > math::Absolute( GetMinMax().First())
          && math::Absolute( lower_y) > math::Absolute( GetMinMax().First())
        )
        {
          const double remainder_x( math::Absolute( lower_x) - math::Absolute( GetMinMax().First()));
          const double remainder_lower_x( GetMinMax().Second() - remainder_x);
          const double remainder_upper_x( GetMinMax().Second());
          const double remainder_y( math::Absolute( lower_y) - math::Absolute( GetMinMax().First()));
          const double remainder_lower_y( GetMinMax().Second() - remainder_y);
          const double remainder_upper_y( GetMinMax().Second());
          storage::VectorND< 2, storage::VectorND< 2, double> > box
          (
            storage::VectorND< 2, double>( remainder_lower_x, remainder_lower_y),
            storage::VectorND< 2, double>( remainder_upper_x, remainder_upper_y)
          );
          box_coords.PushBack( box);
          box_specifications.PushBack( current_specs);
        }

        // upper left corner - wrap to lower right corner
        if
        (
          math::Absolute( lower_x) > math::Absolute( GetMinMax().First())
          && math::Absolute( upper_y) > math::Absolute( GetMinMax().Second())
        )
        {
          const double remainder_x( math::Absolute( lower_x) - math::Absolute( GetMinMax().First()));
          const double remainder_lower_x( GetMinMax().Second() - remainder_x);
          const double remainder_upper_x( GetMinMax().Second());
          const double remainder_y( math::Absolute( upper_y) - math::Absolute( GetMinMax().Second()));
          const double remainder_lower_y( GetMinMax().First());
          const double remainder_upper_y( GetMinMax().First() + remainder_y);
          storage::VectorND< 2, storage::VectorND< 2, double> > box
          (
            storage::VectorND< 2, double>( remainder_lower_x, remainder_lower_y),
            storage::VectorND< 2, double>( remainder_upper_x, remainder_upper_y)
          );
          box_coords.PushBack( box);
          box_specifications.PushBack( current_specs);
        }

        // upper right corner - wrap to lower left corner
        if
        (
          math::Absolute( upper_x) >  math::Absolute( GetMinMax().Second()) &&
          math::Absolute( upper_y) >  math::Absolute( GetMinMax().Second())
        )
        {
          const double remainder_x( math::Absolute( upper_x) - math::Absolute( GetMinMax().Second()));
          const double remainder_lower_x( GetMinMax().First());
          const double remainder_upper_x( GetMinMax().First() + remainder_x);

          const double remainder_y( math::Absolute( upper_y) - math::Absolute( GetMinMax().Second()));
          const double remainder_lower_y( GetMinMax().First());
          const double remainder_upper_y( GetMinMax().First() + remainder_y);
          storage::VectorND< 2, storage::VectorND< 2, double> > box
          (
            storage::VectorND< 2, double>( remainder_lower_x, remainder_lower_y),
            storage::VectorND< 2, double>( remainder_upper_x, remainder_upper_y)
          );
          box_coords.PushBack( box);
          box_specifications.PushBack( current_specs);
        }

        // lower right corner - wrap to upper left corner
        if
        (
          math::Absolute( upper_x) >  math::Absolute( GetMinMax().Second()) &&
          math::Absolute( lower_y) >  math::Absolute( GetMinMax().First())
        )
        {
          const double remainder_x( math::Absolute( upper_x) - math::Absolute( GetMinMax().Second()));
          const double remainder_lower_x( GetMinMax().First());
          const double remainder_upper_x( GetMinMax().First() + remainder_x);
          const double remainder_y( math::Absolute( lower_y) - math::Absolute( GetMinMax().First()));
          const double remainder_lower_y( GetMinMax().Second() - remainder_y);
          const double remainder_upper_y( GetMinMax().Second());
          storage::VectorND< 2, storage::VectorND< 2, double> > box
          (
            storage::VectorND< 2, double>( remainder_lower_x, remainder_lower_y),
            storage::VectorND< 2, double>( remainder_upper_x, remainder_upper_y)
          );
          box_coords.PushBack( box);
          box_specifications.PushBack( current_specs);
        }
      }

      HEAT_MAP.SetBoxes( box_coords, X_RANGE, Y_RANGE, X_MIN, Y_MIN);
      HEAT_MAP.SetBoxesSpecifications( box_specifications);
    }

  } // namespace assemble
} // namespace bcl
