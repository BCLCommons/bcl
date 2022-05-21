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
#include "assemble/bcl_assemble_sse_geometry_compare.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "assemble/bcl_assemble_sse_geometry_phi_psi.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSEGeometryWithinSizeTolerance::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEGeometryWithinSizeTolerance())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief construct from helix and strand length tolerances
    //! @param HELIX_TOLERANCE helix tolerance
    //! @param STRAND_TOLERANCE strand tolerance
    SSEGeometryWithinSizeTolerance::SSEGeometryWithinSizeTolerance
    (
      const size_t HELIX_TOLERANCE,
      const size_t STRAND_TOLERANCE
    ) :
      m_SizeTolerances()
    {
      m_SizeTolerances[ biol::GetSSTypes().HELIX] = HELIX_TOLERANCE;
      m_SizeTolerances[ biol::GetSSTypes().STRAND] = STRAND_TOLERANCE;
    }

    //! @brief Clone function
    //! @return pointer to new SSEGeometryWithinSizeTolerance
    SSEGeometryWithinSizeTolerance *SSEGeometryWithinSizeTolerance::Clone() const
    {
      return new SSEGeometryWithinSizeTolerance( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSEGeometryWithinSizeTolerance::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief returns whether two SSE geometry interface objects are within a fragment length of each other
    //! @param SSE_TO_COMPARE SSE to be compared
    //! @param GEOMTRY_PHI_PSI SSEGeometryPhiPsi to be compared
    //! @return bool whether two SSE geometry interface objects are within a fragment length of each other
    bool SSEGeometryWithinSizeTolerance::operator()
    (
      const SSE &SSE_TO_COMPARE,
      const SSEGeometryPhiPsi &GEOMTRY_PHI_PSI
    ) const
    {
      // if the types don't match or are not helix or strand
      if( SSE_TO_COMPARE.GetType() != GEOMTRY_PHI_PSI.GetType() || !SSE_TO_COMPARE.GetType()->IsStructured())
      {
        return false;
      }

      // return whether the difference it is less than or equal to the tolerance
      return math::Absolute( int( SSE_TO_COMPARE.GetSize()) - int( GEOMTRY_PHI_PSI.GetPhiPsi()->GetAngles().GetSize()))
        <= int( m_SizeTolerances.GetValue( SSE_TO_COMPARE.GetType()));
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEGeometryWithinSizeTolerance::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SizeTolerances, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEGeometryWithinSizeTolerance::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SizeTolerances, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
