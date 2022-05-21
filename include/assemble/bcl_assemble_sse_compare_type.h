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

#ifndef BCL_ASSEMBLE_SSE_COMPARE_TYPE_H_
#define BCL_ASSEMBLE_SSE_COMPARE_TYPE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse.h"
#include "util/bcl_util_binary_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSECompareType
    //! @brief to decide whether two sses are of the same type
    //! @details This class' operator checks whether two sses are of the same type.
    //!
    //! @see @link example_assemble_sse_compare_type.cpp @endlink
    //! @author linders
    //! @date Aug 17, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSECompareType :
      public util::BinaryFunctionInterfaceSerializable< SSE, SSE, bool>
    {

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new SSEsMatchType
      SSECompareType *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator that tests whether the ss type of two sses is the same
      //! @param SSE_A first sse to be tested
      //! @param SSE_B second sse to be tested
      //! @return bool whether the ss type is the same
      bool operator()( const SSE &SSE_A, const SSE &SSE_B) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class SSECompareType

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_COMPARE_TYPE_H_
