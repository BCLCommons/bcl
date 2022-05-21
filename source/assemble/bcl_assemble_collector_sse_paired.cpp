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
#include "assemble/bcl_assemble_collector_sse_paired.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "assemble/bcl_assemble_sse_geometry_packing.h"
#include "io/bcl_io_serialization.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> CollectorSSEPaired::s_Instance
    (
      util::Enumerated< find::CollectorInterface< storage::List< storage::VectorND< 2, util::SiPtr< const SSE> > >, DomainInterface> >::AddInstance( new CollectorSSEPaired())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CollectorSSEPaired::CollectorSSEPaired() :
      m_PackingPicker( GetSSEGeometryPackingPickers().e_Undefined),
      m_ContactTypes(),
      m_MaxDistance( 0.0),
      m_ConsiderOrthogonalConnection( false)
    {
    }

    //! @brief constructor taking a contact type
    //! @param PACKING_PICKER SSEGeometryPackingPicker to be used
    //! @param CONTACT_TYPES is the set of contact type for which sse pairs will be collected
    //! @param MAX_DISTANCE maximum distance to define to sses as paired
    //! @param ORTHOGONAL_CONNECTION boolean whether to include non-orthogonal contacts
    CollectorSSEPaired::CollectorSSEPaired
    (
      const SSEGeometryPackingPicker &PACKING_PICKER,
      const storage::Set< contact::Type> &CONTACT_TYPES,
      const double MAX_DISTANCE,
      const bool ORTHOGONAL_CONNECTION
    ) :
      m_PackingPicker( PACKING_PICKER),
      m_ContactTypes( CONTACT_TYPES),
      m_MaxDistance( MAX_DISTANCE),
      m_ConsiderOrthogonalConnection( ORTHOGONAL_CONNECTION)
    {
    }

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new CollectorSSEPaired which is a copy of this
    CollectorSSEPaired *CollectorSSEPaired::Clone() const
    {
      return new CollectorSSEPaired( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CollectorSSEPaired::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetContactType gives the types of SSE contact of interest
    //! @return the types of SSE contact of interest
    const storage::Set< contact::Type> &CollectorSSEPaired::GetContactTypes() const
    {
      // return "m_ContactTypes"
      return m_ContactTypes;
    }

    //! @brief SetContactType changes the type of SSE contact of interest
    //! @param CONTACT_TYPES is the new type of contact that will be searched for SSEs not having
    void CollectorSSEPaired::SetContactTypes( const storage::Set< contact::Type> &CONTACT_TYPES)
    {
      // set "m_ContactType" to "CONTACT_TYPE"
      m_ContactTypes = CONTACT_TYPES;
    }

    //! @brief GetMaxDistance gives the maximum distance between SSEs for them to still be considered paired
    //! @return returns "m_MaxDistance" the maximum distance between SSEs for them to still be considered paired
    double CollectorSSEPaired::GetMaxDistance() const
    {
      // return "m_MaxDistance"
      return m_MaxDistance;
    }

    //! @brief SetMaxDistance changes the type of SSE contact of interest
    //! @param MAX_DISTANCE is the new maximum distance between SSEs for them to still be considered paired
    void CollectorSSEPaired::SetMaxDistance( const double MAX_DISTANCE)
    {
      // set "m_MaxDistance" to "MAX_DISTANCE"
      m_MaxDistance = MAX_DISTANCE;
    }

    //! @brief returns boolean to whether include non-orthogonal-contacts
    //! @return boolean to whether include non-orthogonal-contacts
    bool CollectorSSEPaired::GetOrthogonalConnection() const
    {
      return m_ConsiderOrthogonalConnection;
    }

    //! @brief sets boolean to whether include non-orthogonal-contacts
    //! @param ORTHOGONAL_CONNECTION boolean to whether include non-orthogonal-contacts
    void CollectorSSEPaired::SetOrthogonalConnection( const bool ORTHOGONAL_CONNECTION)
    {
      m_ConsiderOrthogonalConnection = ORTHOGONAL_CONNECTION;
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &CollectorSSEPaired::GetAlias() const
    {
      static const std::string s_alias( "CollectorSSEPaired");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorSSEPaired::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Collects all paired SSEs.");
      serializer.AddInitializer
      (
        "packing picker",
        "picker for protein geometries",
        io::Serialization::GetAgent( &m_PackingPicker)
      );
      serializer.AddInitializer
      (
        "contact types",
        "types of contacts to consider",
        io::Serialization::GetAgent( &m_ContactTypes)
      );
      serializer.AddInitializer
      (
        "max distance",
        "maximum distance between SSEs to be considered paired",
        io::Serialization::GetAgent( &m_MaxDistance)
      );
      serializer.AddInitializer
      (
        "orthogonal",
        "include orthogonal connections",
        io::Serialization::GetAgent( &m_ConsiderOrthogonalConnection)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief returns the SSEs which are paired in the domain argument
    //! @param SSE_DOMAIN domain from which paired SSEs will be collected
    //! @return returns List of pair of paired sses
    storage::List< storage::VectorND< 2, util::SiPtr< const SSE> > >
    CollectorSSEPaired::Collect( const DomainInterface &SSE_DOMAIN) const
    {
      // get vector of sses from the protein model
      const util::SiPtrVector< const SSE> sse_vector( SSE_DOMAIN.GetSSEs());

      // initialize paired_sses container
      storage::List< storage::VectorND< 2, util::SiPtr< const SSE> > > paired_sses;

      // iterate through "sse_vector"
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr_a( sse_vector.Begin()), sse_itr_end( sse_vector.End());
        sse_itr_a != sse_itr_end;
        ++sse_itr_a
      )
      {
        // iterate through "sse_vector"
        for
        (
          util::SiPtrVector< const SSE>::const_iterator sse_itr_b( sse_itr_a + 1);
          sse_itr_b != sse_itr_end;
          ++sse_itr_b
        )
        {
          // create SSEGeometryPacking "current_packing" from sse pair
          const SSEGeometryPacking this_pack
          (
            ( *m_PackingPicker)->operator ()( **sse_itr_a, **sse_itr_b)
          );

          // if this packing has the correct type, distance below the threshold and is orthogonal if required
          if
          (
            m_ContactTypes.Find( this_pack.GetContactType()) != m_ContactTypes.End() &&
            this_pack.GetShortestConnection().GetLength() <= m_MaxDistance &&
            ( !m_ConsiderOrthogonalConnection || this_pack.GetOrthogonalConnection())
          )
          {
            // insert into paired_sses
            paired_sses.PushBack( storage::VectorND< 2, util::SiPtr< const SSE> >( **sse_itr_a, **sse_itr_b));
          }
        }
      }

      // return
      return paired_sses;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace assemble
} // namespace bcl
