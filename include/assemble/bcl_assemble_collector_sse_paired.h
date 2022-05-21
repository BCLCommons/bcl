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

#ifndef BCL_ASSEMBLE_COLLECTOR_SSE_PAIRED_H_
#define BCL_ASSEMBLE_COLLECTOR_SSE_PAIRED_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_geometry_packing_pickers.h"
#include "contact/bcl_contact_types.h"
#include "find/bcl_find_collector_interface.h"
#include "storage/bcl_storage_set.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorSSEPaired
    //! @brief collects all the paired SSEs of a domain.
    //! @details Collects all pairs of SSE that satisfy the given contact types, the maximum distance and orthogonal
    //! connection boolean
    //!
    //! @see @link example_assemble_collector_sse_paired.cpp @endlink
    //! @author karakam, woetzen
    //! @date 05/22/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CollectorSSEPaired :
      public find::CollectorInterface< storage::List< storage::VectorND< 2, util::SiPtr< const SSE> > >, DomainInterface>
    {

    //////////
    // data //
    //////////

    private:

      //! packer to be used
      SSEGeometryPackingPicker m_PackingPicker;

      //! set of contact types that will be collected
      storage::Set< contact::Type> m_ContactTypes;

      //! the maximum distance between SSEs for them to still be considered paired
      double m_MaxDistance;

      //! boolean value to indicate whether to include only orthogonal contacts
      bool m_ConsiderOrthogonalConnection;

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
      CollectorSSEPaired();

      //! @brief constructor taking a contact type
      //! @param PACKING_PICKER SSEGeometryPackingPicker to be used
      //! @param CONTACT_TYPES is the set of contact type for which sse pairs will be collected
      //! @param MAX_DISTANCE maximum distance to define to sses as paired
      //! @param ORTHOGONAL_CONNECTION boolean whether to include non-orthogonal contacts
      CollectorSSEPaired
      (
        const SSEGeometryPackingPicker &PACKING_PICKER,
        const storage::Set< contact::Type> &CONTACT_TYPES,
        const double MAX_DISTANCE,
        const bool ORTHOGONAL_CONNECTION
      );

      //! @brief Clone is the virtual Clone constructor
      //! @return a pointer to new CollectorSSEPaired which is a copy of this
      CollectorSSEPaired *Clone() const;

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief GetContactType gives the types of SSE contact of interest
      //! @return the types of SSE contact of interest
      const storage::Set< contact::Type> &GetContactTypes() const;

      //! @brief SetContactType changes the type of SSE contact of interest
      //! @param CONTACT_TYPES is the new type of contact that will be searched for SSEs not having
      void SetContactTypes( const storage::Set< contact::Type> &CONTACT_TYPES);

      //! @brief GetMaxDistance gives the maximum distance between SSEs for them to still be considered paired
      //! @return returns "m_MaxDistance" the maximum distance between SSEs for them to still be considered paired
      double GetMaxDistance() const;

      //! @brief SetMaxDistance changes the type of SSE contact of interest
      //! @param MAX_DISTANCE is the new maximum distance between SSEs for them to still be considered paired
      void SetMaxDistance( const double MAX_DISTANCE);

      //! @brief returns boolean to whether include non-orthogonal-contacts
      //! @return boolean to whether include non-orthogonal-contacts
      bool GetOrthogonalConnection() const;

      //! @brief sets boolean to whether include non-orthogonal-contacts
      //! @param ORTHOGONAL_CONNECTION boolean to whether include non-orthogonal-contacts
      void SetOrthogonalConnection( const bool ORTHOGONAL_CONNECTION);

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns the SSEs which are paired in the domain argument
      //! @param SSE_DOMAIN domain from which paired SSEs will be collected
      //! @return returns List of pair of paired sses
      storage::List< storage::VectorND< 2, util::SiPtr< const SSE> > >
      Collect( const DomainInterface &SSE_DOMAIN) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class CollectorSSEPaired

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_COLLECTOR_SSE_PAIRED_H_
