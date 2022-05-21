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

#ifndef BCL_ASSEMBLE_ANALYZE_CHI_ANGLE_PAIR_DISTRIBUTION_H_
#define BCL_ASSEMBLE_ANALYZE_CHI_ANGLE_PAIR_DISTRIBUTION_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_analyze_protein_ensemble_interface.h"
#include "bcl_assemble_collector_aa_type.h"
#include "biol/bcl_biol_chi_angle.h"
#include "find/bcl_find_collector_interface.h"
#include "find/bcl_find_locator_interface.h"
#include "util/bcl_util_sh_ptr.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AnalyzeChiAnglePairDistribution
    //! @brief outputs gnuplot script to create 2D heat map showing frequency with which chi angles are observed
    //! @details Uses a collector interface to find the residues of interest from the protein ensemble. The chi angles
    //!          are calculated for two desired chi angles for the residues of interest. The frequency with which
    //!          a pair of chi angles is observed is shown in the heat map. The possibility to highlight up to two
    //!          different sets of desired chi angle combinations is available (e.g. a) combinations available in a
    //!          rotamer library and b) experimentally observed combinations). These are highlighted via boxes of
    //!          specified dimensions.
    //!
    //! @see @link example_assemble_analyze_chi_angle_pair_distribution.cpp @endlink
    //! @author alexanns
    //! @date Aug 12, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AnalyzeChiAnglePairDistribution :
      public AnalyzeProteinEnsembleInterface
    {

    private:

    //////////
    // data //
    //////////

      //! the postfix that will be appended to the filename that will hold this analysis
      std::string m_OutFilePostFix;

      //! the width of one bin of the chi angle histograms in degrees
      double m_HistogramBinSize;

      //! locator for identifying the residue of interest whose chi angles will be analyzed
      CollectorAAType m_LocatorAA;

      //! first chi angle of interest
      biol::ChiAngle::ChiEnum m_ChiA;

      //! second chi angle of interest
      biol::ChiAngle::ChiEnum m_ChiB;

      //! true if the title and axis labels should be set on the resulting heat map or not
      bool m_SetTitleAndLabel;

      //! the title that will label the resulting heat map
      std::string m_Title;

      //! bool indicating if the color box should be displayed or not - the box that shows the color scale
      bool m_ShowColorBox;

      double m_FontSize; //!< the size text in the heat map should have

      //! the pixel height and width of the heat map graphic
      size_t m_ResolutionX;
      size_t m_ResolutionY;

      //! the name of the file that contains chi angle pairs of interest
      std::string m_ReferenceChiFilename;

      std::string m_Font; //!< font to be used
      bool m_GreyScale;   //!< true if the palette desired is grey scale

      //! how much of the png area the plot will cover top to bottom and left to right
      //! (this does not include any labels or tics)
      double m_PlotRatio;

      size_t m_BinsPerTic; //!< how many bins per tic in the plot
      bool   m_CenterTics; //!< if true, tics will be centered on the bins, if false, will be on edges of bins
      double m_MinZ;       //!< lowest value in z
      double m_MaxZ;       //!< highest value in z
      bool   m_Normalize;  //!< if true, indicates that the distance distribution histogram will be normalized

    public:

      //! single instance of that class
      static const util::SiPtr< const ObjectInterface> s_Instance;

      //! @brief provides the min and max values contained in the heat map
      //! @return min and max values contained in the heat map
      static const storage::VectorND< 2, double> &GetMinMax();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AnalyzeChiAnglePairDistribution();

      //! @brief Clone function
      //! @return pointer to new AnalyzeChiAnglePairDistribution
      AnalyzeChiAnglePairDistribution *Clone() const;

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

      //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
      //! @param ENSEMBLE the ensemble that will be analyzed
      //! @return string which has the analyzed information about the ensemble
      std::string operator()( const ProteinEnsemble &ENSEMBLE) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

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

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief gives distribution for observed combinations of two chi angles for a desired residue in an ensemble
      //! @param ENSEMBLE the ensemble the distribution will be calculated over
      //! @param BIN_SIZE the size of bins that should be used for the histogram distribution
      //! @param LOCATOR the method for locating the residue of interest
      //! @param CHI_A the first chi angle of interest
      //! @param CHI_B the second chi angle of interest
      //! @return histogram 2D which has the distribution of chia and chib angles observed in the ensemble
      math::Histogram2D CalculateChiAngleStatistics
      (
        const ProteinEnsemble &ENSEMBLE, const double BIN_SIZE,
        const find::CollectorInterface
        <
          util::SiPtrList< const biol::AABase>, util::SiPtrVector< const biol::AABase>
         > &LOCATOR,
         const biol::ChiAngle::ChiEnum CHI_A, const biol::ChiAngle::ChiEnum CHI_B
      ) const;

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
      math::GnuplotHeatmap GetHeatMap
      (
        const math::Histogram2D &DISTRIBUTION, const biol::ChiAngle::ChiEnum CHI_A,
        const biol::ChiAngle::ChiEnum CHI_B, const bool SET_TITLE_AND_LABEL, const std::string &TITLE,
        const bool SHOW_COLOR_BOX, const double FONT_SIZE, const size_t RESOLUTION_X, const size_t RESOLUTION_Y,
        const std::string &CHI_FILENAME
      ) const;

      //! @brief provides the text necessary to add any desired boxes to the heat map
      //! @param HEAT_MAP the heat map the boxes will be added to
      //! @param CHI_FILENAME the name of the file that contains chi angle pairs of interest
      //! @param X_RANGE the range in the x direction of the plot - the total size of the x direction
      //! @param Y_RANGE the range in the y direction of the plot -  the total size of the y direction
      //! @param X_MIN the minimum value in the x direction of the heat map
      //! @param Y_MIN the minimum value in the y direction of the heat map
      static void AddBoxesToHeatMap
      (
        math::GnuplotHeatmap &HEAT_MAP, const std::string &CHI_FILENAME,
        const double X_RANGE, const double Y_RANGE, const double X_MIN, const double Y_MIN
      );

    }; // class AnalyzeChiAnglePairDistribution

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_ANALYZE_CHI_ANGLE_PAIR_DISTRIBUTION_H_
