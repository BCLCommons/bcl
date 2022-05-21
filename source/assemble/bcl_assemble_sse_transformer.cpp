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
#include "assemble/bcl_assemble_sse_transformer.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> SSETransformer::s_Instance
    (
      GetObjectInstances().AddInstance( new SSETransformer())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    SSETransformer::SSETransformer() :
      m_Sequence()
    {
    }

    //! @brief construct from a transformation matrix and chain IDs
    //! @param SP_SEQUENCE sequence used to generate new SSE
    //! @param SP_TRANSFORMATION the transformation to apply to sse
    SSETransformer::SSETransformer
    (
      const util::ShPtr< biol::AASequence> &SP_SEQUENCE,
      const util::ShPtr< math::TransformationMatrix3D> &SP_TRANSFORMATION
    ) :
      m_Sequence( SP_SEQUENCE),
      m_Transformation( SP_TRANSFORMATION)
    {
    }

    //! @brief Clone function
    //! @return pointer to new SSETransformer
    SSETransformer *SSETransformer::Clone() const
    {
      return new SSETransformer( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &SSETransformer::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief returns SSE after transformation
    //! @param INITIAL_SSE SSE to be transformed
    //! @return SSE after transformation
    util::ShPtr< SSE> SSETransformer::operator ()( const SSE &INITIAL_SSE) const
    {
      // create a new SSE
      util::ShPtr< SSE> sp_sse( INITIAL_SSE.Clone());

      // connect the aa data and set to the correct chain id
      m_Sequence->ConnectAADataByPdbID( *sp_sse, false);

      // transform
      sp_sse->Transform( *m_Transformation);

      // end
      return sp_sse;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &SSETransformer::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_Sequence, ISTREAM);
      io::Serialize::Read( m_Transformation, ISTREAM);

      // return the stream
      return ISTREAM;
    }

    //! @brief write to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT number of indentations
    //! @return outputstream which was written to
    std::ostream &SSETransformer::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_Sequence, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Transformation, OSTREAM, INDENT);

      // return the stream
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
