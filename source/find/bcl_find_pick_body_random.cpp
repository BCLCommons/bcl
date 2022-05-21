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
#include "find/bcl_find_pick_body_random.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_geometry_interface.h"
#include "io/bcl_io_serialization.h"
#include "util/bcl_util_enumerated.h"
#include "util/bcl_util_sh_ptr_vector.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> PickBodyRandom::s_Instance
    (
      util::Enumerated< PickInterface< util::ShPtr< assemble::SSEGeometryInterface>, util::ShPtrVector< assemble::SSEGeometryInterface> > >::AddInstance( new PickBodyRandom())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PickBodyRandom::PickBodyRandom()
    {
    }

    //! virtual copy constructor
    PickBodyRandom *PickBodyRandom::Clone() const
    {
      return new PickBodyRandom( *this);
    }

    //! virtual destructor
    PickBodyRandom::~PickBodyRandom()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PickBodyRandom::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief get the name of the object when used in a dynamic context
    //! @return the name of the object when used in a dynamic context
    const std::string &PickBodyRandom::GetAlias() const
    {
      static const std::string s_alias( "PickBodyRandom");
      return s_alias;
    }

    //! @brief returns parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer PickBodyRandom::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription( "Picks a random SSE geometry.");

      return serializer;
    }

  ////////////////
  // operations //
  ////////////////

    //! Pick the util::ShPtr< assemble::SSEGeometryInterface> object from restraint::Body based on an SSE
    //! @param BODIES is the object which will provide the pool of util::ShPtr< assemble::SSEGeometryInterface> to pick from
    //! @return return a util::ShPtr< assemble::SSEGeometryInterface> which can accommodate "SSE_CRITERIA"
    util::ShPtr< assemble::SSEGeometryInterface>
    PickBodyRandom::Pick( const util::ShPtrVector< assemble::SSEGeometryInterface> &BODIES) const
    {
      BCL_MessageVrb( "PickBodyRandom::Pick");
      size_t bodies_size( BODIES.GetSize());

      if( !bodies_size)
      {
        BCL_MessageVrb( "PickBodyRandom::Pick returning empy ShPtr since no bodies\n");
        return util::ShPtr< assemble::SSEGeometryInterface>();
      }
      else
      {
        BCL_MessageVrb( "PickBodyRandom::Pick get a body out of BODIES");
        const size_t random_size_t( random::GetGlobalRandom().Random< size_t>( bodies_size - 1));
        BCL_MessageVrb( "PickBodyRandom::Pick random size_t is " + util::Format()( random_size_t) + " and bodies size is " + util::Format()( bodies_size));
        util::ShPtr< assemble::SSEGeometryInterface> picked_body( BODIES( random_size_t));
        BCL_MessageVrb( "PickBodyRandom::Pick picked_body created");
        BCL_MessageVrb( "PickBodyRandom::Pickoutputting body extents: \n");
        BCL_MessageVrb
        (
          "PickBodyRandom::Pick body extents: \n "
          + util::Format()( picked_body->GetExtents())
        );
        BCL_MessageVrb( "PickBodyRandom::Pick end");
        // return an element accessed by a random index within the range of [0, restraint_size - 1]
        return picked_body;
      }
    }

  //////////////////////
  // input and output //
  //////////////////////

  } // namespace find
} // namespace bcl
