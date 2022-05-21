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

#ifndef BCL_ASSEMBLE_LOCATOR_SSE_FURTHEST_H_
#define BCL_ASSEMBLE_LOCATOR_SSE_FURTHEST_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_domain_interface.h"
#include "bcl_assemble_pick_sse_furthest_euclidean_center.h"
#include "find/bcl_find_collector_criteria_interface.h"
#include "find/bcl_find_locator_criteria_interface.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class LocatorSSEFurthest
    //! @brief locates the SSE of a domain which is most distant from the center of the domain
    //! @details This class find the SSE for the given protein model that is farthest from the center of the domain
    //!
    //! @see @link example_assemble_locator_sse_furthest.cpp @endlink
    //! @author alexanns, woetzen
    //! @date 03/18/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    template< typename t_CriteriaType>
    class LocatorSSEFurthest :
      public find::LocatorCriteriaInterface< util::SiPtr< const SSE>, DomainInterface, t_CriteriaType>
    {

    //////////
    // data //
    //////////

      //! collector to be used for locating sse furthest
      util::Implementation< find::CollectorCriteriaInterface< util::SiPtrList< const SSE>, DomainInterface, t_CriteriaType> > m_Collector;

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
      LocatorSSEFurthest() :
        m_Collector()
      {
      }

      //! @brief constructor from a collector
      //! @param COLLECTOR CollectorCriteriaInterface to be used
      LocatorSSEFurthest
      (
        const find::CollectorCriteriaInterface< util::SiPtrList< const SSE>, DomainInterface, t_CriteriaType> &COLLECTOR
      ) :
        m_Collector( COLLECTOR)
      {
      }

      //! @brief constructor from a collector
      //! @param SP_COLLECTOR CollectorCriteriaInterface to be used
      LocatorSSEFurthest
      (
        const util::ShPtr< find::CollectorCriteriaInterface< util::SiPtrList< const SSE>, DomainInterface, t_CriteriaType> > &SP_COLLECTOR
      ) :
        m_Collector( *SP_COLLECTOR)
      {
      }

      //! @brief Clone is the virtual Clone constructor
      //! @return a pointer to new LocatorSSEFurthest which is a copy of this
      LocatorSSEFurthest< t_CriteriaType> *Clone() const
      {
        return new LocatorSSEFurthest< t_CriteriaType>( *this);
      }

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const
      {
        static const std::string s_name( "LocatorSSEFurthest");
        return s_name;
      }

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        io::Serializer serializer;
        serializer.SetClassDescription( "Locates the SSE of a domain that is most distant from the center of the domain.");
        serializer.AddInitializer
        (
          "collector",
          "collector for finding the most distant SSE",
          io::Serialization::GetAgent( &m_Collector)
        );

        return serializer;
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief Locate returns the SSE which is furthest from the center of the domain argument
      //! @param SSE_DOMAIN domain which the LocatorSSE refers to
      //! @param CRITERIA
      //! @return returns SiPtr to the SSE which is furthest from the center of SSE_DOMAIN; empty if no SSEs in domain
      util::SiPtr< const SSE> Locate( const DomainInterface &SSE_DOMAIN, const t_CriteriaType &CRITERIA) const
      {
        // use the pick sse furthest to pick and return
        return PickSSEFurthestEuclideanCenter< t_CriteriaType>().Pick
        (
          m_Collector->Collect( SSE_DOMAIN, CRITERIA),
          CRITERIA
        );
      }

    //////////////////////
    // input and output //
    //////////////////////

    }; // class LocatorSSEFurthest

    // instantiate s_Instance
    template< typename t_CriteriaType>
    const util::SiPtr< const util::ObjectInterface> LocatorSSEFurthest< t_CriteriaType>::s_Instance
    (
      util::Enumerated< find::LocatorCriteriaInterface< util::SiPtr< const SSE>, DomainInterface, t_CriteriaType> >::AddInstance
      (
        new LocatorSSEFurthest< t_CriteriaType>()
      )
    );

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_LOCATOR_SSE_FURTHEST_H_
