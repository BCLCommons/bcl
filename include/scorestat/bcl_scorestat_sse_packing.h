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

#ifndef BCL_SCORESTAT_SSE_PACKING_H_
#define BCL_SCORESTAT_SSE_PACKING_H_

// include the namespace header
#include "bcl_scorestat.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_analyze_protein_ensemble_interface.h"
#include "assemble/bcl_assemble_protein_model.h"
#include "math/bcl_math_histogram.h"
#include "storage/bcl_storage_table.h"
#include "storage/bcl_storage_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace scorestat
  {
    /////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPacking
    //! @brief extracts sse packing statistics from protein models
    //!
    //! @see @link example_sse_packing_statistics.cpp @endlink
    //! @author lib14
    //! @date Nov 25, 2014
    //////////////////////////////////////////////////////////////////////

    class BCL_API SSEPacking :
      public assemble::AnalyzeProteinEnsembleInterface
    {

    public:
        //! output options
        enum OutputOption
        {
          e_Table,
          e_Histogram,
          s_NumberOutputOptions
        };

        //! @brief OutputOption as string
        //! @param OUTPUT_OPTION the OutputOption
        //! @return the string for the OutputOption
        static const std::string &GetOutputOptionName( const OutputOption &OUTPUT_OPTION);

        //! @brief Output filename as string
        //! @param OutputOption the desired Output Type
        //! @return the string for the output file extension
        static const std::string &GetOutputFileName( const OutputOption &OUTPUT_OPTION);

        //! @brief OutputOptionEnum enum I/O helper
        typedef util::WrapperEnum< OutputOption, &GetOutputOptionName, s_NumberOutputOptions> OutputOptionEnum;

    private:

    //////////
    // data //
    //////////

        //! output options
        OutputOptionEnum m_OutputOption;

        //! sse distance bin size
        double m_SSEDistanceBinSize;

        //! number of bins for sse twist angle
        size_t m_SSEAngleNumberBins;

        //! strand distance bin size
        double m_StrandDistanceBinSize;

        //! strand twist angle bin size
        size_t m_StrandAngleNumberBins;

        //! fragment minimal interface length
        double m_FragmentMinInterfaceLength;

        //! chain ids
        std::string m_ChainIds;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SSEPacking();

      //! @brief virtual copy constructor
      SSEPacking *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as constant reference to std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gives the string to append to the the end of a filename to identify this analysis
      //! @return gives the string to append to the the end of a filename to identify this analysis
      const std::string &GetOutFilePostfix() const;

      //! @brief returns the sse distance bin size for the 2D histogram
      //! @return the sse distance bin size for the 2D histogram
      const double &GetSSEDistanceBinSize() const;

      //! @brief returns the sse twist angle bin size for the 2D histogram
      //! @return the sse twist angle bin size for the 2D histogram
      const size_t &GetSSEAngleNumberBins() const;

      //! @brief returns the strand distance bin size for the 2D histogram
      //! @return the strand distance bin size for the 2D histogram
      const double &GetStrandDistanceBinSize() const;

      //! @brief returns the strand twist angle bin size for the 2D histogram
      //! @return the strand twist angle bin size for the 2D histogram
      const size_t &GetStrandAngleNumberBins() const;

      //! @brief returns the fragment minimum interface length
      //! @return the fragment minimum interface length
      const double &GetFragmentMinInterfaceLength() const;

      //! @brief returns chain ids
      //! @return chain ids
      const std::string &GetChainIds() const;

      //! @brief returns the name used for this class in an object data label
      //! @return the name used for this class in an object data label
      const std::string &GetAlias() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief takes an ensemble and does an analysis resulting in string containing formatted information
      //! @param ENSEMBLE the protein ensemble that will be analyzed
      //! @return string which has the analyzed information about the ensemble
      std::string operator()( const assemble::ProteinEnsemble &ENSEMBLE) const;

    protected:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief collects twist angle and shortest distance for sse packing types
      //! @param SSE_PACK the sse packing type
      //!        SSEPACKING_ANGLE_DISTANCE a vector of 2D histograms that store twist angle and shortest distance
      //!        of sse packing types
      //!        STRAND_STRAND_ANGLE_DISTANCE a 2D histogram that stores the twist angle and shortest distance of
      //!        strand_strand packing
      //! @return a vector of 2D histograms that store twist angle and shortest distance of sse packing types
      storage::Vector< math::Histogram2D> &SSEPackingAngleDistance
      (
        const assemble::SSEGeometryPacking &SSE_PACK,
        storage::Vector< math::Histogram2D> &SSEPACKING_ANGLE_DISTANCE,
        math::Histogram2D &STRAND_STRAND_ANGLE_DISTANCE
      ) const;

    }; // end class SSEPackingStatistics
  } // namespace scorestat
} // namespace bcl

#endif // BCL_SCORESTAT_SSE_PACKING_H_
