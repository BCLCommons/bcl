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
#include "score/bcl_score_sse_pair_angle_distance.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor - initialize with memory usage
    SSEPairAngleDistance::SSEPairAngleDistance() :
      m_SSEPackScoringFunction()
    {
    }

    //! @brief constructor from a SSEPackInterface
    //! @param SSE_PACK_INTERFACE reference to a SSEPackInterface derived class
    SSEPairAngleDistance::SSEPairAngleDistance( const SSEPackInterface &SSE_PACK_INTERFACE) :
      m_SSEPackScoringFunction( SSE_PACK_INTERFACE.Clone())
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new SSEPairAngleDistance copied from this one
    SSEPairAngleDistance *SSEPairAngleDistance::Clone() const
    {
      return new SSEPairAngleDistance( *this);
    }

    //! @brief destructor
    SSEPairAngleDistance::~SSEPairAngleDistance()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEPairAngleDistance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &SSEPairAngleDistance::GetScheme() const
    {
      return m_SSEPackScoringFunction->GetScheme();
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &SSEPairAngleDistance::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSEPackScoringFunction, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &SSEPairAngleDistance::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSEPackScoringFunction, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief write the Scheme and the function value for the ARGUMENT to the STREAM
    //! @param SSE_A first SSE of interest
    //! @param SSE_B second SSE of interest
    //! @param OSTREAM the std::ostream to be written to
    //! @return std::ostream which was written to
    std::ostream &
    SSEPairAngleDistance::WriteDetailedSchemeAndValues
    (
      const assemble::SSE &SSE_A,
      const assemble::SSE &SSE_B,
      std::ostream &OSTREAM
    ) const
    {

      // calculate the SSEPack
      const assemble::SSEGeometryPacking sse_pack
      (
        SSE_A, SSE_B, m_SSEPackScoringFunction->GetMinimalInterfaceLength()
      );

      // print SSE information
      // write sstype, Angel to membrane plane first and last aminoacid and value to the STREAM
      OSTREAM << SSE_A.GetIdentification() << '\t'
              << SSE_B.GetIdentification() << '\n';

      // now write the scheme for the scoring function scored
      return m_SSEPackScoringFunction->WriteDetailedSchemeAndValues( sse_pack, OSTREAM);
    }

  } // namespace score
} // namespace bcl
