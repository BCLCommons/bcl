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
#include "score/bcl_score_sse_pair_gap.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEPairGap::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEPairGap( 0.0))
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &SSEPairGap::GetDefaultScheme()
    {
      static const std::string s_scheme( "sse_pair_gap");
      return s_scheme;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a scheme
    //! @param NON_COIL_PENALTY additional penalty to add, if the ss type is not COIL
    //! @param SCHEME Scheme to be used
    SSEPairGap::SSEPairGap
    (
      const double NON_COIL_PENALTY,
      const std::string &SCHEME
    ) :
      m_Scheme( SCHEME),
      m_NonCoilPenalty( NON_COIL_PENALTY)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSEPairGap
    SSEPairGap *SSEPairGap::Clone() const
    {
      return new SSEPairGap( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEPairGap::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &SSEPairGap::GetScheme() const
    {
      return m_Scheme;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &SSEPairGap::GetAlias() const
    {
      static const std::string s_name( "SSEPairGap");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SSEPairGap::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Returns the number of residues between two SSEs as penalty.");
      serializer.AddInitializer
      (
        "non coil penalty",
        "penalty if both ends are not coils",
        io::Serialization::GetAgent( &m_NonCoilPenalty)
      );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief gives the score for two sses based on how close they are to being connected, if they are adjecent
    //! @param SSE_A the first sse whose closeness to being peptide bonded to the second sse will be scored
    //! @param SSE_B the second sse whose closeness to being peptide bonded to the first sse will be scored
    //! @return the score for two sses based on how close SSE_A and SSE_B are to being connected by a peptide bond
    double SSEPairGap::operator()( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const
    {
      size_t sequence_separation( biol::SequenceSeparation( *SSE_A.GetLastAA(), *SSE_B.GetFirstAA()));

      // no gap, score is zero
      if( !util::IsDefined( sequence_separation) || sequence_separation == 0)
      {
        return 0.0;
      }

      // score is proportinal to sequence separation
      double score( sequence_separation);

      // gap between structured sses is penalized more
      score += SSE_A.GetType()->IsStructured() * m_NonCoilPenalty;
      score += SSE_B.GetType()->IsStructured() * m_NonCoilPenalty;

      BCL_MessageDbg
      (
        "gap score between structured SSEs: |" + SSE_A.GetIdentification() + '|' + SSE_B.GetIdentification() + '|' +
        util::Format()( score)
      );

      return double( score);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEPairGap::Read( std::istream &ISTREAM)
    {
      // read members

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEPairGap::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members

      // return the stream
      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
