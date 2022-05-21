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

#ifndef BCL_ASSEMBLE_LOCATOR_SSE_RANDOM_H_
#define BCL_ASSEMBLE_LOCATOR_SSE_RANDOM_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_locator_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorSSERandom
    //! @brief locator implementation that randomly selects an SSE from a given domain.
    //! @details It collects all sses from a domain. Generate a random index on the collection and returns a SiPtr to this
    //! SSE. The SiPtr will be undefined, when the domain does not contain SSEs.
    //!
    //! @see @link example_assemble_locator_sse_random.cpp @endlink
    //! @author karakam
    //! @date Mar 2, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorSSERandom :
      public find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface>
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

      //! default constructor
      LocatorSSERandom()
      {
      }

      //! clone constructor
      LocatorSSERandom *Clone() const
      {
        return new LocatorSSERandom( *this);
      }

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Locate returns a random SSE from the domain argument
      //! @param SSE_DOMAIN domain which the LocatorSSE refers to
      //! @return returns SiPtr to the SSE denoted by the LocatorSSE
      util::SiPtr< const SSE> Locate( const DomainInterface &SSE_DOMAIN) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class LocatorSSERandom

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_LOCATOR_SSE_RANDOM_H_
