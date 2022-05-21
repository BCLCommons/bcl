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
#include "assemble/bcl_assemble_sse_pair_template.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SSEPairTemplate::s_Instance
    (
      GetObjectInstances().AddInstance( new SSEPairTemplate())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSEPairTemplate::SSEPairTemplate() :
      m_FirstSSEGeometry(),
      m_SecondSSEGeometry(),
      m_LoopLength( util::GetUndefined< size_t>()),
      m_Packing()
    {
    }

    //! @brief construct from two SSE geometries
    //! @param SP_FIRST_SSE_GEOMETRY first SSE geometry
    //! @param SP_SECOND_SSE_GEOMETRY second SSE geometry
    //! @param LOOP_LENGTH length of the loop connecting the two geometries
    SSEPairTemplate::SSEPairTemplate
    (
      const util::ShPtr< SSEGeometryInterface> &SP_FIRST_SSE_GEOMETRY,
      const util::ShPtr< SSEGeometryInterface> &SP_SECOND_SSE_GEOMETRY,
      const size_t LOOP_LENGTH
    ) :
      m_FirstSSEGeometry( SP_FIRST_SSE_GEOMETRY),
      m_SecondSSEGeometry( SP_SECOND_SSE_GEOMETRY),
      m_LoopLength( LOOP_LENGTH),
      m_Packing( *SP_FIRST_SSE_GEOMETRY, *SP_SECOND_SSE_GEOMETRY)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSEPairTemplate
    SSEPairTemplate *SSEPairTemplate::Clone() const
    {
      return new SSEPairTemplate( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &SSEPairTemplate::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSEPairTemplate::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_FirstSSEGeometry, ISTREAM);
      io::Serialize::Read( m_SecondSSEGeometry, ISTREAM);
      io::Serialize::Read( m_LoopLength, ISTREAM);
      io::Serialize::Read( m_Packing, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM outputstream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSEPairTemplate::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_FirstSSEGeometry, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SecondSSEGeometry, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_LoopLength, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Packing, OSTREAM, INDENT) << '\n';

      // end
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
