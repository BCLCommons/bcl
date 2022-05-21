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

#ifndef BCL_ASSEMBLE_LOCATOR_SSES_RANDOM_H_
#define BCL_ASSEMBLE_LOCATOR_SSES_RANDOM_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_collector_interface.h"
#include "find/bcl_find_locator_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorSSEsRandom
    //! @brief selects a given number of random SSEs from collected SSEs from a given domain.
    //! @details It collects all sses from a domain using the specified collector and pick the specified number of
    //! random SSEs and returns them
    //! If the domain does not have enough number of SSEs it will return empty
    //!
    //! @see @link example_assemble_locator_sses_random.cpp @endlink
    //! @author karakam
    //! @date Mar 2, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API LocatorSSEsRandom :
      public find::LocatorInterface< util::SiPtrList< const SSE>, DomainInterface>
    {

    private:

    //////////
    // data //
    //////////

      //! number of SSEs to be located
      size_t m_NumberSSEs;

      //! ShPtr to SSE collector to be used
      util::ShPtr< find::CollectorInterface< util::SiPtrList< const SSE>, DomainInterface> > m_SSECollector;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      LocatorSSEsRandom();

      //! @brief constructor from number of SSEs to collect and a ShPtr to SSE Collector
      //! @param NUMBER_SSES Number of SSEs to be located
      //! @param SP_COLLECTOR_SSE SSE ShPtr to SSE Collector to be used
      LocatorSSEsRandom
      (
        const size_t NUMBER_SSES,
        const util::ShPtr< find::CollectorInterface< util::SiPtrList< const SSE>, DomainInterface> > &SP_COLLECTOR_SSE
      );

      //! @brief clone constructor
      LocatorSSEsRandom *Clone() const;

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Locate returns a random SSE from the domain argument
      //! @param SSE_DOMAIN domain which the LocatorSSE refers to
      //! @return returns SiPtrList to the SSEs denoted by the collector
      util::SiPtrList< const SSE> Locate( const DomainInterface &SSE_DOMAIN) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief read from std::ostream
      //! @param OSTREAM input stream
      //! @return ostream which was read from
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class LocatorSSEsRandom

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_LOCATOR_SSES_RANDOM_H_
