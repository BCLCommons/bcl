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

#ifndef BCL_ASSEMBLE_COLLECTOR_ALL_POSSIBLE_DOMAINS_H_
#define BCL_ASSEMBLE_COLLECTOR_ALL_POSSIBLE_DOMAINS_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_domain_interface.h"
#include "find/bcl_find_collector_interface.h"
#include "math/bcl_math_function_interface_serializable.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorAllPossibleDomains
    //! @brief collects domains from a domain
    //! @details takes a domain and returns all the possible domains of the requested size
    //!
    //! @see @link example_assemble_collector_all_possible_domains.cpp @endlink
    //! @author weinerbe
    //! @date Sep 16, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CollectorAllPossibleDomains :
      public find::CollectorInterface< util::ShPtrVector< Domain>, DomainInterface>
    {

    private:

    //////////
    // data //
    //////////

      //! size of domains will be randomly chosen using this range
      math::Range< size_t> m_DomainSize;

      //! bool whether the remaining SSEs in the domain should also be a domain
      bool m_ForceDomainConnectivity;

      //! packing criteria - SSEs must be packed with at least one other SSE in the domain to be considered
      util::Implementation< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> > m_PackingCriteria;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief construct from the domain size and a packing criteria
      //! @param DOMAIN_SIZE domain size range
      //! @param CONNECTED_DOMAIN bool whether the SSEs in the domain should remain connected
      //! @param PACKING_CRITERIA packing criteria to be used
      CollectorAllPossibleDomains
      (
        const math::Range< size_t> DOMAIN_SIZE = math::Range< size_t>( 0, std::numeric_limits< size_t>::max()),
        const bool CONNECTED_DOMAIN = false,
        const util::ShPtr< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> > &PACKING_CRITERIA =
          GetDefaultPackingCriteria()
      );

      //! @brief Clone function
      //! @return pointer to new CollectorAllPossibleDomains
      CollectorAllPossibleDomains *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return packing criteria function
      //! @return packing criteria function
      static const util::ShPtr< math::FunctionInterfaceSerializable< SSEGeometryPacking, bool> > &GetDefaultPackingCriteria();

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Collect returns all domains in the domain argument
      //! @param SSE_DOMAIN domain from which Domains will be collected
      //! @return returns ShPtrVector of domains found in the given domain
      util::ShPtrVector< Domain> Collect( const DomainInterface &SSE_DOMAIN) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class CollectorAllPossibleDomains

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_COLLECTOR_ALL_POSSIBLE_DOMAINS_H_ 
