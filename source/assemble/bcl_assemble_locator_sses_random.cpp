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
#include "assemble/bcl_assemble_locator_sses_random.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_sse.h"
#include "assemble/bcl_assemble_domain_interface.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> LocatorSSEsRandom::s_Instance
    (
      GetObjectInstances().AddInstance( new LocatorSSEsRandom())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    LocatorSSEsRandom::LocatorSSEsRandom() :
      m_NumberSSEs( 2),
      m_SSECollector( new CollectorSSE())
    {
    }

    //! @brief constructor from number of SSEs to collect and a ShPtr to SSE Collector
    //! @param NUMBER_SSES Number of SSEs to be located
    //! @param SP_COLLECTOR_SSE SSE ShPtr to SSE Collector to be used
    LocatorSSEsRandom::LocatorSSEsRandom
    (
      const size_t NUMBER_SSES,
      const util::ShPtr< find::CollectorInterface< util::SiPtrList< const SSE>, DomainInterface> > &SP_COLLECTOR_SSE
    ) :
      m_NumberSSEs( NUMBER_SSES),
      m_SSECollector( SP_COLLECTOR_SSE)
    {
    }

    //! @brief clone constructor
    LocatorSSEsRandom *LocatorSSEsRandom::Clone() const
    {
      return new LocatorSSEsRandom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorSSEsRandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Locate returns a random SSE from the domain argument
    //! @param SSE_DOMAIN domain which the LocatorSSE refers to
    //! @return returns SiPtrList to the SSEs denoted by the collector
    util::SiPtrList< const SSE> LocatorSSEsRandom::Locate( const DomainInterface &SSE_DOMAIN) const
    {
      // collect the SSEs
      util::SiPtrList< const SSE> collected_sses( m_SSECollector->Collect( SSE_DOMAIN));

      // initialize list to return
      util::SiPtrList< const SSE> sse_list;

      // if there are not enough collected sses
      if( collected_sses.GetSize() < m_NumberSSEs)
      {
        // return empty list
        return sse_list;
      }

      // iterate until the requested number of SSEs is reached
      for( size_t sse_ctr( 0); sse_ctr < m_NumberSSEs; ++sse_ctr)
      {
        // get random iterator on the list
        util::SiPtrList< const SSE>::iterator random_itr
        (
          random::GetGlobalRandom().Iterator( collected_sses.Begin(), collected_sses.End(), collected_sses.GetSize())
        );

        // double check that the random iterator is not equal to end
        BCL_Assert
        (
          random_itr != collected_sses.End(),
          "The random iterator returns the end which probably indicates an error with the list"
        );

        // insert it into sse_list
        sse_list.PushBack( *random_itr);

        // remove from collected_sses list
        collected_sses.Remove( random_itr);
      }

      // return the list
      return sse_list;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &LocatorSSEsRandom::Read( std::istream &ISTREAM)
    {
      // read members
      io::Serialize::Read( m_SSECollector, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read from std::ostream
    //! @param OSTREAM input stream
    //! @param INDENT indentation
    //! @return ostream which was read from
    std::ostream &LocatorSSEsRandom::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_SSECollector, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
