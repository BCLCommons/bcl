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

#ifndef BCL_RESTRAINT_ANALYZE_COORDINATE_DISTANCE_DISTRIBUTION_H_
#define BCL_RESTRAINT_ANALYZE_COORDINATE_DISTANCE_DISTRIBUTION_H_

// include the namespace header
#include "bcl_restraint.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "assemble/bcl_assemble_locator_atom.h"
#include "find/bcl_find_locator_interface.h"
#include "io/bcl_io_directory_entry.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AnalyzeCoordinateDistanceDistribution
    //! @brief Creates heat map of showing frequency with which a distance between two coordinates is observed.
    //! @details The heat map of frequency with which a distance is observed for two coordinates located in a protein
    //!          ensemble can be compared with a second heat map. The second heat map is given by a math function.
    //!          This allows the distribution of distances coming from the ensemble to be compared with a second
    //!          distribution which is perhaps the expected distribution. An example use of this is to see the
    //!          distribution of distances coming from two residues in the ensemble whose distance has been measured
    //!          by EPR. This ensemble distribution of distances can then be compared to the EPR distribution described
    //!          by a gaussian function defined by the mean and standard deviation of the distance probability
    //!          distribution measured by EPR.
    //!
    //! @see @link example_restraint_AnalyzeCoordinateDistanceDistribution.cpp @endlink
    //! @author alexanns
    //! @date Sep 1, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AnalyzeCoordinateDistanceDistribution :
      public assemble::AnalyzeProteinEnsembleInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the postfix that will be appended to the filename that will hold this analysis
      std::string m_OutFilePostFix;

      //! the first coordinate that will be used in the distance
      util::Implementation< find::LocatorCoordinatesInterface< assemble::ProteinModel> > m_CoordinateA;

      //! the second coordinate that will be used in the distance
      util::Implementation< find::LocatorCoordinatesInterface< assemble::ProteinModel> > m_CoordinateB;

      //! the function that will be plotted to compare against the coordinate distance distribution
      util::Implementation< math::FunctionInterfaceSerializable< double, double> > m_ComparisonFunction;

      //! the minimal value representing the left boundary of the score histogram
      double m_HistogramMinimum;

      //! the width of one bin of the score histograms
      double m_HistogramBinSize;

      //! the number of bins in the score histograms
      size_t m_HistogramNumberOfBins;

      std::string m_Title;    //!< the title the plot will have
      size_t      m_PixelX;   //!< size of produced plot in pixel x
      size_t      m_PixelY;   //!< size of produced plot in pixel y
      std::string m_Font;     //!< font to be used
      size_t      m_FontSize; //!< size of font

      //! true if the palette desired is grey scale
      bool m_GreyScale;

      //! how much of the png area each plot will cover top to bottom (this does not include any labels or tics)
      double m_PlotRatio;

      //! filename that will hold the mean and standard deviation of the distribution
      std::string m_MeanStdDevOutFile;

      size_t m_BinsPerTic; //!< how many bins per tic in the plot
      bool   m_CenterTics; //!< if true, tics will be centered on the bins, if false, will be on edges of bins
      double m_MinZ;       //!< lowest value in z
      double m_MaxZ;       //!< highest value in z
      bool   m_Normalize;  //!< if true, indicates that the distance distribution histogram will be normalized

    public:

      //! single instance of that class
      static const util::SiPtr< const ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AnalyzeCoordinateDistanceDistribution();

      //! @brief Clone function
      //! @return pointer to new AnalyzeCoordinateDistanceDistribution
      AnalyzeCoordinateDistanceDistribution *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives the string to append to the the end of a filename to identify this analysis
      //! @return gives the string to append to the the end of a filename to identify this analysis
      const std::string &GetOutFilePostfix() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
      //! @param ENSEMBLE the ensemble that will be analyzed
      //! @return string which has the analyzed information about the ensemble
      std::string operator()( const assemble::ProteinEnsemble &ENSEMBLE) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief gives the histogram of the distance distribution of the protein ensemble
      //! @param ENSEMBLE the ensemble for which the distance distribution will be calculated
      //! @return histogram that has the distance distribution calculated from the provided protein ensemble
      math::Histogram GetDistanceHistogram( const assemble::ProteinEnsemble &ENSEMBLE) const;

      //! @brief gives the heat map corresponding to the provided histogram made from the protein ensemble
      //! @param HISTGRAM the histogram for which the heat map will be made
      //! @return heat map representing the distance distribution of the protein ensemble
      util::ShPtr< math::GnuplotHeatmap> GetHeatMap( const math::Histogram &HISTOGRAM) const;

      //! @brief gives the heat map corresponding to the reference function
      //! @return shptr to a heat map that represents the reference function
      util::ShPtr< math::GnuplotHeatmap> GetReferenceHeatMap() const;

    }; // class AnalyzeCoordinateDistanceDistribution

  } // namespace restraint
} // namespace bcl

#endif // BCL_RESTRAINT_ANALYZE_COORDINATE_DISTANCE_DISTRIBUTION_H_
