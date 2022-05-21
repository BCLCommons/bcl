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
#include "assemble/bcl_assemble_sse_compare_type.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> SSECompareType::s_Instance
    (
      util::Enumerated< util::BinaryFunctionInterfaceSerializable< SSE, SSE, bool> >::AddInstance( new SSECompareType())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief virtual copy constructor
    SSECompareType *SSECompareType::Clone() const
    {
      return new SSECompareType( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &SSECompareType::GetAlias() const
    {
      static const std::string s_alias( "SSECompareType");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer SSECompareType::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Compares the type of SSEs");

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that tests whether the ss type of two sses is the same
    //! @param SSE_A first sse to be tested
    //! @param SSE_B second sse to be tested
    //! @return bool whether the ss type is the same
    bool SSECompareType::operator()( const SSE &SSE_A, const SSE &SSE_B) const
    {
      // return whether the ss type of the two sses is the same
      return SSE_A.GetType() == SSE_B.GetType();
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace assemble
} // namespace bcl
