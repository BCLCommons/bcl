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
#include "assemble/bcl_assemble_collector_sse_size.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_domain_interface.h"
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

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> CollectorSSESize::s_Instance
    (
      util::Enumerated< find::CollectorInterface< util::SiPtrList< const SSE>, DomainInterface> >::AddInstance
      (
        new CollectorSSESize( math::Range< size_t>( 0, util::GetUndefined< size_t>()))
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from sse size range, applies to all SSTypes
    //! @param SIZE_RANGE size range for sse sequence length
    CollectorSSESize::CollectorSSESize( const math::Range< size_t> &SIZE_RANGE) :
      m_SizeRangeMap()
    {
      m_SizeRangeMap[ biol::GetSSTypes().HELIX]  = SIZE_RANGE;
      m_SizeRangeMap[ biol::GetSSTypes().STRAND] = SIZE_RANGE;
      m_SizeRangeMap[ biol::GetSSTypes().COIL]   = SIZE_RANGE;
    }

    //! @brief map of SSTypes and corresponding size ranges
    //! @param SIZE_RANGE_MAP size range for sse sequence length
    CollectorSSESize::CollectorSSESize( const storage::Map< biol::SSType, math::Range< size_t> > &SIZE_RANGE_MAP) :
      m_SizeRangeMap( SIZE_RANGE_MAP)
    {
    }

    //! @brief Clone function
    //! @return pointer to new CollectorSSESize
    CollectorSSESize *CollectorSSESize::Clone() const
    {
      return new CollectorSSESize( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &CollectorSSESize::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &CollectorSSESize::GetAlias() const
    {
      static const std::string s_alias( "CollectorSSESize");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer CollectorSSESize::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Collects SSEs with specified sizes.");
      serializer.AddInitializer
      (
        "size map",
        "size ranges for each SSE type",
        io::Serialization::GetAgent( &m_SizeRangeMap)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief collect all sses from domain that fulfill the desired size range
    //! @param SSE_DOMAIN domain to collect from
    //! @return sses is a SiPtrList that have a sequence length according to range
    util::SiPtrList< const SSE> CollectorSSESize::Collect( const DomainInterface &SSE_DOMAIN) const
    {
      // initialize collected SSEs list
      util::SiPtrList< const SSE> collected_sses;

      // get all SSEs from the domain
      util::SiPtrVector< const SSE> domain_sses( SSE_DOMAIN.GetSSEs());

      // iterate over domain sses
      for
      (
        util::SiPtrVector< const SSE>::const_iterator itr( domain_sses.Begin()), itr_end( domain_sses.End());
        itr != itr_end;
        ++itr
      )
      {
        // check to find if there is range for this SSE
        storage::Map< biol::SSType, math::Range< size_t> >::const_iterator range_itr
        (
          m_SizeRangeMap.Find( ( *itr)->GetType())
        );

        // if range not found continue
        if( range_itr == m_SizeRangeMap.End())
        {
          continue;
        }

        // if sequence length is with range
        if( range_itr->second.IsWithin( ( *itr)->GetSize()))
        {
          collected_sses.PushBack( *itr);
        }
      }

      // end
      return collected_sses;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace assemble
} // namespace bcl
