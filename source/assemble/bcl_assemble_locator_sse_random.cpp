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
#include "assemble/bcl_assemble_locator_sse_random.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_collector_sse.h"
#include "assemble/bcl_assemble_pick_sse_random.h"
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
    const util::SiPtr< const util::ObjectInterface> LocatorSSERandom::s_Instance
    (
      util::Enumerated< find::LocatorInterface< util::SiPtr< const SSE>, DomainInterface> >::AddInstance( new LocatorSSERandom())
    );

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &LocatorSSERandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &LocatorSSERandom::GetAlias() const
    {
      static const std::string s_alias( "LocatorSSERandom");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer LocatorSSERandom::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Locates a random SSE in a given domain.");

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Locate returns a random SSE from the domain argument
    //! @param SSE_DOMAIN domain which the LocatorSSE refers to
    //! @return returns SiPtr to the SSE denoted by the LocatorSSE
    util::SiPtr< const SSE> LocatorSSERandom::Locate( const DomainInterface &SSE_DOMAIN) const
    {
      // initialize static pick and static collectorsse
      static const PickSSERandom s_picker;
      static const CollectorSSE s_collector;

      // return a random sse from all sses in the protein model
      return s_picker.Pick( s_collector.Collect( SSE_DOMAIN));
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace assemble
} // namespace bcl
