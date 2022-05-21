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

#ifndef BCL_ASSEMBLE_COLLECTOR_SSE_UNPAIRED_H_
#define BCL_ASSEMBLE_COLLECTOR_SSE_UNPAIRED_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_domain_interface.h"
#include "contact/bcl_contact_types.h"
#include "find/bcl_find_collector_interface.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorSSEUnpaired
    //! @brief collects all the unpaired SSEs of a domain.
    //! @details collects all SSE pairs of the specified contact type that are farther than the given maximum distance
    //!
    //! @see @link example_assemble_collector_sse_unpaired.cpp @endlink
    //! @author alexanns, woetzen
    //! @date 03/27/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CollectorSSEUnpaired :
      public find::CollectorInterface< util::SiPtrList< const SSE>, DomainInterface>
    {

    //////////
    // data //
    //////////

    private:

      contact::Type m_ContactType; //!< the contact type of interest; an SSE not having it is desired
      double        m_MaxDistance; //!< the maximum distance between SSEs for them to still be considered paired

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
      CollectorSSEUnpaired();

      //! @brief constructor taking a contact type
      //! @param CONTACT_TYPE is the contact type for which an SSE not having it is desired
      //! @param MAX_DISTANCE is the maximum distance
      CollectorSSEUnpaired( const contact::Type &CONTACT_TYPE, const double MAX_DISTANCE);

      //! @brief Clone is the virtual Clone constructor
      //! @return a pointer to new LocatorSSEFurthest which is a copy of this
      CollectorSSEUnpaired *Clone() const;

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief GetContactType gives the type of SSE contact of interest
      //! @return returns "m_ContactType" which is the contact type that the locator searches for SSEs not having
      const contact::Type &GetContactType() const;

      //! @brief SetContactType changes the type of SSE contact of interest
      //! @param CONTACT_TYPE is the new type of contact that will be searched for SSEs not having
      void SetContactType( const contact::Type &CONTACT_TYPE);

      //! @brief GetMaxDistance gives the maximum distance between SSEs for them to still be considered paired
      //! @return returns "m_MaxDistance" the maximum distance between SSEs for them to still be considered paired
      double GetMaxDistance() const;

      //! @brief SetMaxDistance changes the type of SSE contact of interest
      //! @param MAX_DISTANCE is the new maximum distance between SSEs for them to still be considered paired
      void SetMaxDistance( const double MAX_DISTANCE);

      //! @brief get a short name for this class
      //! @return a short name for this class
      const std::string &GetAlias() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief Collect returns the SSEs which are unpaired in the domain argument
      //! @param SSE_DOMAIN domain from which unpaired SSEs will be collected
      //! @return returns SiPtrList of the collected unpaired SSE objects
      util::SiPtrList< const SSE> Collect( const DomainInterface &SSE_DOMAIN) const;

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

    }; // class CollectorSSEUnpaired

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_COLLECTOR_SSE_UNPAIRED_H_
