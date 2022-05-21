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

#ifndef BCL_SCORE_STRAND_PAIRING_H_
#define BCL_SCORE_STRAND_PAIRING_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_sse_pack_interface.h"
#include "contact/bcl_contact_types.h"
#include "math/bcl_math_bicubic_spline.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StrandPairing
    //! @brief This is a Function derived class for scoring pairing of strands within one sheet
    //!
    //! @see @link example_score_strand_pairing.cpp @endlink
    //! @author woetzen, karakam
    //! @date 03.05.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StrandPairing :
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

      //! distance range
      math::Range< double> m_DistanceRange;

      //! minimal interface length for packing to get a full weight
      double m_MinimalInterfaceLength;

      //! ShPtr to energy function to be used
      util::ShPtr< math::BicubicSpline> m_EnergyFunction;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;
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

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a specified histogram file
      //! @param SCHEME scheme to be used
      //! @param HISTOGRAM_FILENAME filename of the histogram to be used
      //! @param DISTANCE_RANGE distance range
      StrandPairing
      (
        const std::string &SCHEME = GetDefaultScheme(),
        const std::string &HISTOGRAM_FILENAME = GetDefaultHistogramFilename(),
        const math::Range< double> &DISTANCE_RANGE = contact::GetTypes().STRAND_STRAND->GetDistanceRange()
      );

      //! @brief virtual copy constructor
      //! @return pointer to a new StrandPairing object that is copied from this one
      StrandPairing *Clone() const;

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

      //! @brief returns the minimal interface length used for calculating packing
      //! @return the minimal interface length used for calculating packing
      const double GetMinimalInterfaceLength() const
      {
        return m_MinimalInterfaceLength;
      }

      //! @brief returns distance range
      //! @return distance range
      const math::Range< double> &GetDistanceRange() const
      {
        return m_DistanceRange;
      }

      //! @brief return reference to the energy function
      //! @return const ref to a BicubicSpline where x is the twist angle and y the distance
      const math::BicubicSpline &GetEnergyFunction() const
      {
        return *m_EnergyFunction;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief returns whether the given pair of SSEs are valid
      //! @param SSE_A first SSE of interest
      //! @param SSE_B second SSE of interest
      //! @return whether the given pair of SSEs are valid
      bool AreValidSSEs( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const;

      //! @brief operator that calculates strand pair packing potential
      //! @param SSE_PACK SSEGeometryPacking of pair of SSEs of interest
      //! @return strand pairing potential
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
      //! @param OSTREAM the std::ostream to be written to
      //! @return std::ostream which was written to
      std::ostream &
      WriteDetailedSchemeAndValues
      (
        const assemble::SSEGeometryPacking &SSE_PACK,
        std::ostream &OSTREAM
      ) const;

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief read energy distribution for scoring pairs of  of strands within one sheet
      void ReadEnergyVector();

      //! @brief Set the members of this property from the given LABEL
      //! @param LABEL the label to parse
      //! @param ERR_STREAM stream to write out errors to
      bool ReadInitializerSuccessHook
      (
        const util::ObjectDataLabel &LABEL,
        std::ostream &ERR_STREAM
      );

    }; //class StrandPairing
  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_STRAND_PAIRING_H_
