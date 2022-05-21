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

#ifndef BCL_SCORE_SSE_PAIR_PACKING_H_
#define BCL_SCORE_SSE_PAIR_PACKING_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_sse_pack_interface.h"
#include "contact/bcl_contact_types.h"
#include "math/bcl_math_bicubic_spline.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPairPacking
    //! @brief This is a Function derived template class for scoring secondary structure element pairs
    //!
    //! @see @link example_score_sse_pair_packing.cpp @endlink
    //! @author woetzen, karakam
    //! @date 03.05.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEPairPacking :
      public SSEPackInterface
    {

    private:

    //////////
    // data //
    //////////

      //! scheme to be used in outputting schemes
      std::string m_Scheme;

      //! path to file where the statistics and in consequence the energy potentials are read from
      std::string m_HistogramFileName;

      //! map that contains distance ranges for valid contact types
      storage::Map< contact::Type, math::Range< double> > m_DistanceRangeMap;

      //! minimal interface length for packing to get a full weight
      double m_MinimalInterfaceLength;

      //! ShPtrVector of energy functions to be used
      storage::Map< contact::Type, util::ShPtr< math::BicubicSpline> > m_EnergyFunctions;

    public:

    //////////
    // data //
    //////////

      //! @brief returns default file where the statistics and in consequence the energy potentials are read from
      //! @return default file where the statistics and in consequence the energy potentials are read from
      static const std::string &GetDefaultHistogramFilename();

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a specified histogram file
      //! @param SCHEME scheme to be used
      //! @param HISTOGRAM_FILENAME filename of the histogram to be used
      //! @param DISTANCE_RANGE_MAP map of distance ranges for valid contact types
      SSEPairPacking
      (
        const std::string &SCHEME = GetDefaultScheme(),
        const std::string &HISTOGRAM_FILENAME = GetDefaultHistogramFilename(),
        const storage::Map< contact::Type, math::Range< double> > &DISTANCE_RANGE_MAP = contact::GetTypes().GetValidDistanceRanges()
      );

      //! @brief virtual copy constructor
      SSEPairPacking *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns scheme being used
      //! @return scheme being used
      const std::string &GetScheme() const
      {
        return m_Scheme;
      }

      //! @brief returns filename of the histogram being used
      //! @return filename of the histogram being used
      const std::string &GetHistogramFilename() const
      {
        return m_HistogramFileName;
      }

      //! @brief return distance range map
      //! @return distance range map
      const storage::Map< contact::Type, math::Range< double> > &GetDistanceRangeMap() const
      {
        return m_DistanceRangeMap;
      }

      //! @brief returns the minimal interface length used for calculating packing
      //! @return the minimal interface length used for calculating packing
      const double GetMinimalInterfaceLength() const
      {
        return m_MinimalInterfaceLength;
      }

      //! @brief get the energy functions
      //! @return map for each contact type an associated BicubicSpline with angle as x and distance as y
      const storage::Map< contact::Type, util::ShPtr< math::BicubicSpline> > &GetEnergyFunctions() const
      {
        return m_EnergyFunctions;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns whether the given pair of SSEs are valid
      //! @param SSE_A first SSE of interest
      //! @param SSE_B second SSE of interest
      //! @return whether the given pair of SSEs are valid
      bool AreValidSSEs( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const
      {
        return true;
      }

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that score packing of a pair of SSEs of interest
      //! @param SSE_PACK SSEGeometryPacking of pair of SSEs of interest
      //! @return potential of interaction
      double operator()( const assemble::SSEGeometryPacking &SSE_PACK) const;

    //////////////////////
    // input and output //
    //////////////////////
    protected:

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    public:

      //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
      //! @param SSE_PACK SSEGeometryPacking of pair of SSEs of interest
      //! @param OSTREAM std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &WriteDetailedSchemeAndValues
      (
        const assemble::SSEGeometryPacking &SSE_PACK,
        std::ostream &OSTREAM
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief read energy distribution for scoring pairs of sses
      void ReadEnergyVector();

    }; // class SSEPairPacking

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_SSE_PAIR_PACKING_H_
