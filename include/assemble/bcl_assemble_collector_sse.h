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

#ifndef BCL_ASSEMBLE_COLLECTOR_SSE_H_
#define BCL_ASSEMBLE_COLLECTOR_SSE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_ss_types.h"
#include "find/bcl_find_collector_interface.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorSSE
    //! @brief collects all SSEs of a domain
    //! @details Collects and returns all SSEs of given SSTypes for the provided domain
    //!
    //! @see @link example_assemble_collector_sse.cpp @endlink
    //! @author karakam, alexanns
    //! @date 04/15/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API CollectorSSE :
      public find::CollectorInterface< util::SiPtrList< const SSE>, DomainInterface>
    {

    //////////
    // data //
    //////////

      storage::Set< biol::SSType> m_SSTypes; //!< Set of SSTypes that are going to be collected

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
      CollectorSSE();

      //! @brief constructor from a provided sstype
      //! @param SS_TYPES Set of SSTypes to be collected
      CollectorSSE( const storage::Set< biol::SSType> &SS_TYPES);

      //! @brief constructor from a provided sstype
      //! @param SS_TYPE to be collected
      CollectorSSE( const biol::SSType &SS_TYPE);

      //! @brief Clone is the virtual Clone constructor
      //! @return a pointer to new LocatorSSEFurthest which is a copy of this
      CollectorSSE *Clone() const;

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief returns SSTypes Set
      //! @return SSTypes Set
      const storage::Set< biol::SSType> &GetSSTypes() const;

      //! @brief set SSTypes Set
      //! @param SS_TYPES set of ss types
      void SetSSTypes( const storage::Set< biol::SSType> &SS_TYPES);

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Collect returns all SSEs in the domain argument
      //! @param SSE_DOMAIN domain from which SSEs will be collected
      //! @return returns SiPtrList of the collected SSEs objects
      util::SiPtrList< const SSE> Collect( const DomainInterface &SSE_DOMAIN) const;

      //! Collect returns all SSEs in the SiPtrVector that match the type
      //! @param SSES all SSEs
      //! @return returns SiPtrList of the collected SSEs objects
      util::SiPtrList< const SSE> Collect( const util::SiPtrVector< const SSE> &SSES) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class CollectorSSE

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_COLLECTOR_SSE_H_
