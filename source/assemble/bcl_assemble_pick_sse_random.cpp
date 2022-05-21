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
#include "assemble/bcl_assemble_pick_sse_random.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
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
    const util::SiPtr< const util::ObjectInterface> PickSSERandom::s_Instance
    (
      util::Enumerated< find::PickInterface< util::SiPtr< const SSE>, util::SiPtrList< const SSE> > >::AddInstance( new PickSSERandom())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PickSSERandom::PickSSERandom() :
      m_SSTypes()
    {
    }

    //! @brief constructor from a set of SSTypes
    //! @param SS_TYPES Set of SSTypes
    PickSSERandom::PickSSERandom( const storage::Set< biol::SSType> &SS_TYPES) :
      m_SSTypes( SS_TYPES)
    {
    }

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorSSEFurthest which is a copy of this
    PickSSERandom *PickSSERandom::Clone() const
    {
      return new PickSSERandom( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PickSSERandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &PickSSERandom::GetAlias() const
    {
      static const std::string s_alias( "PickSSERandom");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PickSSERandom::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Randomly pick an SSE of a given type.");
      serializer.AddInitializer
      (
        "sse_types",
        "which types of SSEs to pick",
        io::Serialization::GetAgent( &m_SSTypes)
      );

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Picks the SSE object which is furthest away from "m_ReferencePoint"
    //! @param SSE_LIST is the SiPtrList which provides the pool of assemble::SSE to pick from
    //! @return returns SiPtr to the assemble::SSE object which is furthest away from "m_ReferencePoint"
    util::SiPtr< const SSE>
    PickSSERandom::Pick
    (
      const util::SiPtrList< const SSE> &SSE_LIST
    ) const
    {
      // if not found, return empty SiPtr< SSE>
      if( SSE_LIST.IsEmpty())
      {
        return util::SiPtr< const SSE>();
      }

      // make a copy of the list
      util::SiPtrList< const SSE> sse_list;

      // if SSTypes is not empty meaning a specific set of SSTypes was requested
      if( !m_SSTypes.IsEmpty())
      {
        // iterate over the sse_list and collect SSEs of the correct types.
        for
        (
          util::SiPtrList< const SSE>::const_iterator sse_itr( SSE_LIST.Begin()), sse_itr_end( SSE_LIST.End());
          sse_itr != sse_itr_end; ++sse_itr
        )
        {
          // if of correct type
          if( m_SSTypes.Find( ( *sse_itr)->GetType()) != m_SSTypes.End())
          {
            // insert it
            sse_list.PushBack( *sse_itr);
          }
        }
      }
      else
      {
        sse_list = SSE_LIST;
      }

      // if not found, return empty SiPtr< SSE>
      if( sse_list.IsEmpty())
      {
        return util::SiPtr< const SSE>();
      }

      // get random iterator on the list
      util::SiPtrList< const SSE>::const_iterator random_iterator
      (
        random::GetGlobalRandom().Iterator
        (
          sse_list.Begin(),
          sse_list.End(),
          sse_list.GetSize()
        )
      );

      // double check that the random iterator is not equal to end
      BCL_Assert
      (
        random_iterator != sse_list.End(),
        "The random iterator returns the end which probably indicates an error with the list"
      );

      // return
      return *random_iterator;
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace assemble
} // namespace bcl
