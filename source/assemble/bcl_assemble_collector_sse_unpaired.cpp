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
#include "assemble/bcl_assemble_collector_sse_unpaired.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "io/bcl_io_serialization.h"
#include "restraint/bcl_restraint_group.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> CollectorSSEUnpaired::s_Instance
    (
      GetObjectInstances().AddInstance( new CollectorSSEUnpaired())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CollectorSSEUnpaired::CollectorSSEUnpaired()
    {
    }

    //! @brief constructor taking a contact type
    //! @param CONTACT_TYPE is the contact type for which an SSE not having it is desired
    CollectorSSEUnpaired::CollectorSSEUnpaired
    (
      const contact::Type &CONTACT_TYPE, const double MAX_DISTANCE
    ) :
      m_ContactType( CONTACT_TYPE),
      m_MaxDistance( MAX_DISTANCE)
    {
    }

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorSSEFurthest which is a copy of this
    CollectorSSEUnpaired *CollectorSSEUnpaired::Clone() const
    {
      return new CollectorSSEUnpaired( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CollectorSSEUnpaired::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetContactType gives the type of SSE contact of interest
    //! @return returns "m_ContactType" which is the contact type that the locator searches for SSEs not having
    const contact::Type &CollectorSSEUnpaired::GetContactType() const
    {
      // return "m_ContactType"
      return m_ContactType;
    }

    //! @brief SetContactType changes the type of SSE contact of interest
    //! @param CONTACT_TYPE is the new type of contact that will be searched for SSEs not having
    void CollectorSSEUnpaired::SetContactType( const contact::Type &CONTACT_TYPE)
    {
      // set "m_ContactType" to "CONTACT_TYPE"
      m_ContactType = CONTACT_TYPE;
    }

    //! @brief GetMaxDistance gives the maximum distance between SSEs for them to still be considered paired
    //! @return returns "m_MaxDistance" the maximum distance between SSEs for them to still be considered paired
    double CollectorSSEUnpaired::GetMaxDistance() const
    {
      // return "m_MaxDistance"
      return m_MaxDistance;
    }

    //! @brief SetMaxDistance changes the type of SSE contact of interest
    //! @param MAX_DISTANCE is the new maximum distance between SSEs for them to still be considered paired
    void CollectorSSEUnpaired::SetMaxDistance( const double MAX_DISTANCE)
    {
      // set "m_MaxDistance" to "MAX_DISTANCE"
      m_MaxDistance = MAX_DISTANCE;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &CollectorSSEUnpaired::GetAlias() const
    {
      static const std::string s_name( "CollectorSSEUnpaired");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorSSEUnpaired::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Collects all unpaired SSEs in a domain.");
      serializer.AddInitializer
      (
        "contact type",
        "SSEs not having this contact are collected",
        io::Serialization::GetAgent( &m_ContactType)
      );
      serializer.AddInitializer
      (
        "maximum distance",
        "maximum distance between SSEs to still be considered paired",
        io::Serialization::GetAgent( &m_MaxDistance)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Collect returns the SSEs which are unpaired in the domain argument
    //! @param SSE_DOMAIN domain from which unpaired SSEs will be collected
    //! @return returns SiPtrList of the collected unpaired SSE objects
    util::SiPtrList< const SSE> CollectorSSEUnpaired::Collect( const DomainInterface &SSE_DOMAIN) const
    {
      // create Vector3D "protein_center" and initialize with the coordinates of the center of "PROTEIN_MODEL"
      linal::Vector3D protein_center( SSE_DOMAIN.GetCenter());

      // create SiPtrVector "sse_vector" and initialize to all the SSEs that could be involved in "m_ContactType"
      util::SiPtrVector< const SSE> sse_vector( SSE_DOMAIN.GetSSEs( m_ContactType->GetSSTypes()));

      // create SiPtrList "unpaired_sses"
      util::SiPtrList< const SSE> unpaired_sses;

      // iterate through "sse_vector"
      for
      (
        storage::Vector< util::SiPtr< const SSE> >::const_iterator
          sse_itr_a( sse_vector.Begin()), sse_itr_end( sse_vector.End());
        sse_itr_a != sse_itr_end;
        ++sse_itr_a
      )
      {
        // boolean to indicate whether current sse is paired or not
        bool is_paired( false);

        // iterate through "sse_vector"
        for
        (
          storage::Vector< util::SiPtr< const SSE> >::const_iterator sse_itr_b( sse_vector.Begin());
          sse_itr_b != sse_itr_end;
          ++sse_itr_b
        )
        {
          if( sse_itr_a != sse_itr_b)
          {
            // check that the distance between SSEs denoted by "sse_itr_a" and "sse_itr_b" is less than or equal to
            // "m_MaxDistance"
            if( linal::Distance( ( *sse_itr_a)->GetCenter(), ( *sse_itr_b)->GetCenter()) <= m_MaxDistance)
            {
              // create SSEGeometryPacking "current_packing"; initialize with SSEs denoted by  "sse_itr_a" and "sse_itr_b"
              SSEGeometryPacking current_packing( **sse_itr_a, **sse_itr_b);

              // check if the SSEs denoted by "sse_itr" and "sse_itr_b" are paired together with type "m_ContactType"
              if( current_packing.GetContactType() == m_ContactType)
              {
                // SSE of "sse_itr" has correct pairing so set "is_paired" to true
                is_paired = true;
                // SSE of "sse_itr" has correct pairing so no need to check it with other SSEs
                break;
              }
            }
          }
        }
        // check if SSE denoted by "sse_itr_a" was not found to be correctly paired with any other SSE
        if( !is_paired)
        {
          // SSE denoted by "sse_itr" does not have "m_ContactType" pairing with any SSEs so add to "unpaired_sses"
          unpaired_sses.PushBack( **sse_itr_a);
        }
      }

      // return "unpaired_sses"
      return unpaired_sses;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &CollectorSSEUnpaired::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_ContactType, ISTREAM);
      io::Serialize::Read( m_MaxDistance, ISTREAM);

      // return
      return ISTREAM;
    }

    //! @brief read from std::ostream
    //! @param OSTREAM input stream
    //! @param INDENT indentation
    //! @return ostream which was read from
    std::ostream &CollectorSSEUnpaired::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ContactType, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MaxDistance, OSTREAM, INDENT);

      // return
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
