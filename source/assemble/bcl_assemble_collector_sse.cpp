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
#include "assemble/bcl_assemble_collector_sse.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"
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
    const util::SiPtr< const util::ObjectInterface> CollectorSSE::s_Instance
    (
      util::Enumerated< find::CollectorInterface< util::SiPtrList< const SSE>, DomainInterface> >::AddInstance
      (
        new CollectorSSE()
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    CollectorSSE::CollectorSSE() :
      m_SSTypes( biol::GetSSTypes().HELIX.GetIterator(), biol::GetSSTypes().COIL.GetIterator() + 1)
    {
    }

    //! @brief constructor from a provided sstype
    //! @param SS_TYPES Set of SSTypes to be collected
    CollectorSSE::CollectorSSE( const storage::Set< biol::SSType> &SS_TYPES) :
      m_SSTypes( SS_TYPES)
    {
    }

    //! @brief constructor from a provided sstype
    //! @param SS_TYPE to be collected
    CollectorSSE::CollectorSSE( const biol::SSType &SS_TYPE) :
      m_SSTypes( SS_TYPE)
    {
    }

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorSSEFurthest which is a copy of this
    CollectorSSE *CollectorSSE::Clone() const
    {
      return new CollectorSSE( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CollectorSSE::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief returns SSTypes Set
    //! @return SSTypes Set
    const storage::Set< biol::SSType> &CollectorSSE::GetSSTypes() const
    {
      return m_SSTypes;
    }

    //! @brief set SSTypes Set
    //! @param SS_TYPES set of ss types
    void CollectorSSE::SetSSTypes( const storage::Set< biol::SSType> &SS_TYPES)
    {
      m_SSTypes = SS_TYPES;
    }

    //! @brief get a short name for this class
    //! @return a short name for this class
    const std::string &CollectorSSE::GetAlias() const
    {
      static const std::string s_name( "CollectorSSE");
      return s_name;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorSSE::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Collects all SSEs in a domain.");
      serializer.AddInitializer
      (
        "sse types",
        "SSE types to collect",
        io::Serialization::GetAgent( &m_SSTypes)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Collect returns all SSEs in the domain argument
    //! @param SSE_DOMAIN domain from which SSEs will be collected
    //! @return returns SiPtrList of the collected SSEs objects
    util::SiPtrList< const SSE> CollectorSSE::Collect( const DomainInterface &SSE_DOMAIN) const
    {
      // get SSEs from the PROTEIN_MODEL( get all types of SSEs if m_SSTypes is empty)
      util::SiPtrVector< const SSE> sses
      (
        SSE_DOMAIN.GetSSEs( m_SSTypes)
      );

      // construct SiPtrList on the go and return it
      return util::SiPtrList< const SSE>( sses.Begin(), sses.End());
    }

    //! Collect returns all SSEs in the SiPtrVector that match the type
    //! @param SSES all SSEs
    //! @return returns SiPtrList of the collected SSEs objects
    util::SiPtrList< const SSE> CollectorSSE::Collect( const util::SiPtrVector< const SSE> &SSES) const
    {
      // initialize list to return
      util::SiPtrList< const SSE> sse_list;

      // iterate through the SSEs
      for
      (
        util::SiPtrVector< const SSE>::const_iterator sse_itr( SSES.Begin()), sse_itr_end( SSES.End());
        sse_itr != sse_itr_end; ++sse_itr
      )
      {
        // if the SSE is the right type
        if( m_SSTypes.Contains( ( *sse_itr)->GetType()))
        {
          // add the sse to the list
          sse_list.PushBack( *sse_itr);
        }
      }

      // end
      return sse_list;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace assemble
} // namespace bcl
