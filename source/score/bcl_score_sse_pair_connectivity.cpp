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
#include "score/bcl_score_sse_pair_connectivity.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "biol/bcl_biol_atom.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEPairConnectivity::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEPairConnectivity())
    );

    //! @brief returns default scheme
    //! @return default scheme
    const std::string &SSEPairConnectivity::GetDefaultScheme()
    {
      // static string
      static const std::string s_default_scheme( "sse_connectivity");

      // end
      return s_default_scheme;
    }

    //! @brief GetExtendedResidueLength gives the length a completely extended residue can cover
    //! @return the length a completely extended residue can cover
    const double &SSEPairConnectivity::GetExtendedResidueLength()
    {
      // the length of an extended residue
      static const double s_extended_residue_length( biol::GetSSTypes().STRAND->GetRiseInZPerResidue());

      return s_extended_residue_length;
    }

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a scheme
    //! @param SCHEME Scheme to be used
    SSEPairConnectivity::SSEPairConnectivity
    (
      const std::string &SCHEME
    ) :
      m_Scheme( SCHEME)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSEPairConnectivity
    SSEPairConnectivity *SSEPairConnectivity::Clone() const
    {
      return new SSEPairConnectivity( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEPairConnectivity::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &SSEPairConnectivity::GetScheme() const
    {
      return m_Scheme;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief gives the score for two sses based on how close they are to being connected, if they are adjecent
    //! @param SSE_A the first sse whose closeness to being peptide bonded to the second sse will be scored
    //! @param SSE_B the second sse whose closeness to being peptide bonded to the first sse will be scored
    //! @return the score for two sses based on how close SSE_A and SSE_B are to being connected by a peptide bond
    double SSEPairConnectivity::operator()( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const
    {
      // message which sses are currently being used to calculate connectivity score
      BCL_MessageDbg
      (
        "CalculateConnectivityScore for sse " + SSE_A.GetIdentification() + " and " +
        SSE_B.GetIdentification()
      );

      const size_t sequence_separation( biol::SequenceSeparation( *SSE_A.GetLastAA(), *SSE_B.GetFirstAA()));

      if( !util::IsDefined( sequence_separation))
      {
        return 0.0;
      }

      // calculate how far away the C atom of the last residue of SSE_A is from the N atom of the first residue of SSE_B
      const double c_n_distance
      (
        linal::Distance
        (
          SSE_A.GetLastAA()->GetAtom( biol::GetAtomTypes().C).GetCoordinates(),
          SSE_B.GetFirstAA()->GetAtom( biol::GetAtomTypes().N).GetCoordinates()
        )
      );

      // check if distance is defined
      if( !util::IsDefined( c_n_distance))
      {
        return 0.0;
      }

      // peptide bond length
      static const double s_peptide_bond_length( biol::GetAtomTypes().C->GetBondLength( biol::GetAtomTypes().N) * 1.02);

      // difference of actual distance to maximal distance that could be bridged
      const double safety( c_n_distance - ( sequence_separation * GetExtendedResidueLength() + s_peptide_bond_length));

      // if there is safety margin, return no penalty
      if( safety < 0)
      {
        return 0.0;
      }
      // saftey is crossed; return square as penalty
      else
      {
        return math::Sqr( safety);
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEPairConnectivity::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Scheme, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEPairConnectivity::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Scheme, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief write the Scheme and the Function value for the ARGUMENT to the STREAM
    //! @param SSE_A the first sse whose closeness to being peptide bonded to the second sse will be scored
    //! @param SSE_B the second sse whose closeness to being peptide bonded to the first sse will be scored
    //! @return std::ostream which was written to
    std::ostream &
    SSEPairConnectivity::WriteDetailedSchemeAndValues
    (
      const assemble::SSE &SSE_A, const assemble::SSE &SSE_B,
      std::ostream &OSTREAM
    ) const
    {
      // write the detailed values and score
      OSTREAM << SSE_A.GetIdentification() << '\t' << SSE_B.GetIdentification() << '\t'
              << operator()( SSE_A, SSE_B) << '\n';

      return OSTREAM;
    }

  } // namespace score
} // namespace bcl
