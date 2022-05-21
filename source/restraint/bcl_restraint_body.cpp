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
#include "restraint/bcl_restraint_body.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"
#include "restraint/bcl_restraint_assignment.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace restraint
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> Body::s_Instance
    (
      GetObjectInstances().AddInstance( new Body())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    Body::Body() :
      m_Bodies(),
      m_DetermineOccupancy()
    {
    }

    //! @brief construct from bodies
    //! @param BODIES the ShPtrVector of bodies which will serve as restraints
    //! @param DETERMINE_OCCUPANCY is a BinaryFunctionInterface which is used to determine if a restraint body is occupied by a SSE
    Body::Body
    (
      const util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > &BODIES,
      const util::ShPtr< util::BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool> >
        &DETERMINE_OCCUPANCY
    ) :
      m_Bodies( BODIES),
      m_DetermineOccupancy( DETERMINE_OCCUPANCY)
    {
    }

    //! @brief virtual copy constructor
    Body *Body::Clone() const
    {
      return new Body( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &Body::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief GetBody returns a const reference to "m_Bodies"
    //! @return returns a const reference to "m_Bodies"
    const util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > &Body::GetBody() const
    {
      return m_Bodies;
    }

    //! @brief SetBody "m_Bodies" to a new ShPtr of bodies
    //! @param BODIES is the ShPtrVector of bodies which "m_Bodies" will be changed to
    void Body::SetBody( const util::ShPtr< util::ShPtrVector< assemble::SSEGeometryInterface> > &BODIES)
    {
      m_Bodies = BODIES;
    }

  ///////////////
  // operation //
  ///////////////

    //! @brief GetUnoccupied determines which of the bodies in "m_Bodies" is not occupied
    //! @param SSES a ShPtrVector of sses which could occupy "m_Bodies"
    //! @return returns a ShPtrVector of bodies which are not occupied by any of "BODIES"
    util::ShPtrVector< assemble::SSEGeometryInterface>
    Body::GetUnoccupied( const util::SiPtrVector< const assemble::SSE> &SSES) const
    {
      BCL_Assert( m_DetermineOccupancy.IsDefined(), "restraint::Body::GetUnoccupied m_DetermineOccupancy is not defined");

      // create ShPtrVector "unoccupied" to hold bodies of "m_Bodies" that are not occupied by "BODIES"
      util::ShPtrVector< assemble::SSEGeometryInterface> unoccupied;

      // iterate through bodies
      for
      (
        util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator
          itr_body( m_Bodies->Begin()), itr_body_end( m_Bodies->End());
        itr_body != itr_body_end;
        ++itr_body
      )
      {
        bool is_occupied( false);

        // iterate over all SSEs and check if it occupies current body
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator
            itr_sse( SSES.Begin()), itr_sse_end( SSES.End());
          itr_sse != itr_sse_end;
          ++itr_sse
        )
        {
          // true if the body denoted by "itr_sse" falls within the body denoted by "itr_body"
          if( m_DetermineOccupancy->operator()( **itr_body, **itr_sse))
          {
            // set "is_occupied" to true since the body denoted by "itr_body" is occupied
            is_occupied = true;

            // leave this for loop
            break;
          }
        }

        // true if the body denoted by "itr_body" is unoccupied
        if( !is_occupied)
        {
          // add the unoccupied body denoted by "itr_body" to "unoccupied"
          unoccupied.PushBack( *itr_body);
        }
      }

      BCL_MessageDbg
      (
        "restraint::Body::GetUnoccupied unoccupied size " + util::Format()( unoccupied.GetSize())
      );

      // return "unoccupied"
      return unoccupied;
    }

    //! @brief get body that is occupied by given SSE
    //! @param SSE the sse to be considered to identify the body that is occupied by it
    //! @return SHPtr to occupied body - will be undefined if there is non
    util::ShPtr< assemble::SSEGeometryInterface> Body::GetOccupiedBody( const assemble::SSE &SSE) const
    {
      // iterate through bodies to find which one is occupied by the SSE
      for
      (
        util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator
          itr_body( m_Bodies->Begin()), itr_body_end( m_Bodies->End());
        itr_body != itr_body_end;
        ++itr_body
      )
      {
        // true if the body denoted by "itr_body" is occupied by "BODY"
        if( m_DetermineOccupancy->operator()( **itr_body, SSE))
        {
          return *itr_body;
        }
      }

      return util::ShPtr< assemble::SSEGeometryInterface>();
    }

    //! @brief GenerateAssignment creates the assignment of "m_Bodies" with other assemble::SSEs
    //! @param SSES SiPtrVector of SSE which will be assigned with "m_Bodies"
    //! @return returns an Assignment which assigns "m_Bodies" with the appropriate members of "SSES"
    const SSEAssignment
    Body::GenerateAssignment( const util::SiPtrVector< const assemble::SSE> &SSES) const
    {
      // create GroupCollection "group_collection" to hold elements of BODIES which occupy "m_Bodies"
      GroupCollection< size_t, assemble::SSE> occupied_group_collection;
      BCL_Assert( m_DetermineOccupancy.IsDefined(), "restraint::Body::GenerateAssignment m_DetermineOccupancy is not defined");

      // iterate through "m_Bodies" for instance the density rods
      for
      (
        util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator itr_body( m_Bodies->Begin()), itr_body_end( m_Bodies->End());
        // util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator itr_body( m_Bodies.Begin()), itr_body_end( m_Bodies.End());
        itr_body != itr_body_end;
        ++itr_body
      )
      {
        // iterate through SSES
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator
            itr_sse( SSES.Begin()), itr_sse_end( SSES.End());
          itr_sse != itr_sse_end;
          ++itr_sse
        )
        {
          const bool occupants( m_DetermineOccupancy->operator()( **itr_body, **itr_sse));
          BCL_MessageDbg
          (
            "origin of restraint bodies (e.g. density rod): " + util::Format()( ( *itr_body)->GetCenter()) +
            ", origin of argument (e.g. protein model SSE): " + util::Format()( ( *itr_sse)->GetCenter()) +
            ", occupied?: " + util::Format()( occupants)
          );
          if( occupants)
          {
            BCL_MessageDbg
            (
              "start point of density rod " + util::Format()( ( *itr_body)->GetCenter() - ( *itr_body)->GetAxis( coord::GetAxes().e_Z) * ( *itr_body)->GetExtent( coord::GetAxes().e_Z)) +
              " end point of density rod " + util::Format()( ( *itr_body)->GetCenter() + ( *itr_body)->GetAxis( coord::GetAxes().e_Z) * ( *itr_body)->GetExtent( coord::GetAxes().e_Z)) +
              " radius of density rod " + util::Format()( ( *itr_body)->GetExtent( coord::GetAxes().e_X))
            );
            BCL_MessageDbg
            (
              "start point of sse " + util::Format()( ( *itr_sse)->GetCenter() - ( *itr_sse)->GetAxis( coord::GetAxes().e_Z) * ( *itr_sse)->GetExtent( coord::GetAxes().e_Z)) +
              " end point of sse " + util::Format()( ( *itr_sse)->GetCenter() + ( *itr_sse)->GetAxis( coord::GetAxes().e_Z) * ( *itr_sse)->GetExtent( coord::GetAxes().e_Z)) +
              " radius of density rod " + util::Format()( ( *itr_sse)->GetExtent( coord::GetAxes().e_X))
            );
          }
          BCL_MessageDbg( "occupancy gotten\n ");

          // true if the body denoted by "itr_sse" occupies the body denoted by "itr_body" where occupancy is determined by "m_DetermineOccupancy"
          // for example if the origins of the body are the same
          if( occupants)
          {
            //BCL_Assert( math::EqualWithinTolerance( ( *itr_body)->GetOrigin(), ( *itr_sse)->GetOrigin()), "origin of bodies is not the same");
            BCL_MessageDbg( "getting index\n ");

            const size_t index( itr_body - m_Bodies->Begin());
            BCL_MessageDbg( "index is " + util::Format()( index));

            // insert the element of "BODIES" into "occupied_group_collection"
            occupied_group_collection.Insert( index, Group< assemble::SSE>( 1, *itr_sse));

            // leave this for loop since the occupant of the body denoted by "itr_body" has been found
            break;
          }
        }
      }

//        BCL_Assert
//        (
//          occupied_group_collection.TotalDepth() == SSES.GetSize(),
//            "\n \nevery sse in the protein model should be within a restraint body, therefore the number of bodies in SSES "
//            "(from the protein models) should equal the number of occupied restraint bodies \n but the number of occupied restraint bodies is \n"
//          + util::Format()( occupied_group_collection.TotalDepth())
//          + "\n and the number of sses in the protein model is "
//          + util::Format()( SSES.GetSize())
//          + "\n\n the occupied_group_collection (corresponding to occupied density rods) is "
//          + util::Format()( occupied_group_collection)
//          + "\n while the protein model contains the following SSEs "
//          + util::Format()( SSES)
//        );

      // return Assignment constructed with "m_Bodies" and "occupied_group_collection"
      return SSEAssignment
      (
        m_Bodies, occupied_group_collection
      );
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read distance from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &Body::Read( std::istream &ISTREAM)
    {
      // read member
      io::Serialize::Read( m_Bodies, ISTREAM);
      io::Serialize::Read( m_DetermineOccupancy, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief write distance to std::ostream
    //! @param OSTREAM output stream to write to
    //! @return output stream which was written to
    std::ostream &Body::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Bodies, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_DetermineOccupancy, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace restraint
} // namespace bcl

