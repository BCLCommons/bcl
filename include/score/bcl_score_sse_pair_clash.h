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

#ifndef BCL_SCORE_SSE_PAIR_CLASH_H_
#define BCL_SCORE_SSE_PAIR_CLASH_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_sse_pack_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPairClash
    //! @brief SSPackInterface derived class to identify and score clashing SSEs
    //! @details Depending on the distance in the SSEGeometryPacking object passed to the operator, it will return 0
    //! if the distance is larger than the minimal SSEDistance for that contact type.
    //!
    //! @see @link example_score_sse_pair_clash.cpp @endlink
    //! @author woetzen, karakam
    //! @date 26.05.2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API SSEPairClash :
      public SSEPackInterface
    {

    private:

    //////////
    // data //
    //////////

      //! single instance of this class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

      //! scheme to be used in outputting schemes
      std::string m_Scheme;

      //! minimal interface length for packing to get a full weight
      double m_MinimalInterfaceLength;

      //! width of the sigmoidal function to be used
      double m_SigmoidWidth;

    public:

    //////////
    // data //
    //////////

      //! @brief returns default scheme
      //! @return default scheme
      static const std::string &GetDefaultScheme();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a scheme
      //! @param SIGMOID_WIDTH width of the sigmoid repulsive term, before it reaches 1.0
      //! @param SCHEME scheme to be used
      SSEPairClash
      (
        const double SIGMOID_WIDTH = 1.0,
        const std::string &SCHEME = GetDefaultScheme()
      );

      //! @brief virtual copy constructor
      SSEPairClash *Clone() const;

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

      //! @brief returns the minimal interface length used for calculating packing
      //! @return the minimal interface length used for calculating packing
      const double GetMinimalInterfaceLength() const
      {
        return m_MinimalInterfaceLength;
      }

      //! @brief gets the width of the sigmoidal function
      //! @return the width of the sigmoidal function
      const double GetSigmoidWidth() const
      {
        return m_SigmoidWidth;
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

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

      //! @brief repulsive term from distance and shortest observed distance
      //! @param DISTANCE the actual distance between sses of interest
      //! @param SHORTEST_OBSERVED_DISTANCE shortest distance observed
      double CalculateRepulsiveTerm
      (
        const double DISTANCE,
        const double SHORTEST_OBSERVED_DISTANCE
      ) const;

    }; //class SSEPairClash

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_SSE_PAIR_CLASH_H_
