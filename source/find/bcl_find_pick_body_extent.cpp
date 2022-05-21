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
#include "find/bcl_find_pick_body_extent.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_geometry_interface.h"
#include "find/bcl_find_pick_body_random.h"
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
    const util::SiPtr< const util::ObjectInterface> PickBodyExtent::s_Instance
    (
      GetObjectInstances().AddInstance( new PickBodyExtent())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    PickBodyExtent::PickBodyExtent() :
      m_ExtentUpperBoundTolerance( linal::Vector3D()),
      m_ExtentLowerBoundTolerance( linal::Vector3D())
    {
    }

    //! @brief construct from two linal::Vector3D which describe the upper and lower bound extent tolerance
    PickBodyExtent::PickBodyExtent( const linal::Vector3D &UPPER_TOLERANCE, const linal::Vector3D &LOWER_TOLERANCE) :
      m_ExtentUpperBoundTolerance( UPPER_TOLERANCE),
      m_ExtentLowerBoundTolerance( LOWER_TOLERANCE)
    {
    }

    //! @brief virtual copy constructor
    PickBodyExtent *PickBodyExtent::Clone() const
    {
      return new PickBodyExtent( *this);
    }

    //! virtual destructor
    PickBodyExtent::~PickBodyExtent()
    {
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PickBodyExtent::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! Pick the util::ShPtr< assemble::SSEGeometryInterface> object from restraint::Body based on an SSE
    //! @param BODIES is the object which will provide the pool of util::ShPtr< assemble::SSEGeometryInterface> to pick from
    //! @param CRITERIA the assemble::SSEGeometryInterface which provides the extent criteria by which a assemble::SSEGeometryInterface is chosen
    //! @return return a util::ShPtr< assemble::SSEGeometryInterface> which can accommodate "CRITERIA"
    util::ShPtr< assemble::SSEGeometryInterface>
    PickBodyExtent::Pick( const util::ShPtrVector< assemble::SSEGeometryInterface> &BODIES, const assemble::SSEGeometryInterface &CRITERIA) const
    {
      // create linal::Vector3D "criteria_extents" and initialize with the extents of "CRITERIA"
      linal::Vector3D criteria_extents( CRITERIA.GetExtents());

      BCL_MessageDbg( "criteria_extents " + util::Format()( criteria_extents));

      // create ShPtrVector "eligible_bodies" which will hold the bodies of "BODIES" that match "CRITERIA"
      // within the tolerances of "m_ExtentUpperBoundTolerance" and "m_ExtentLowerBoundTolerance"
      util::ShPtrVector< assemble::SSEGeometryInterface> eligible_bodies;

      // iterate through the bodies of "RESTRAINT"
      for
      (
        util::ShPtrVector< assemble::SSEGeometryInterface>::const_iterator
          itr( BODIES.Begin()), itr_end( BODIES.End());
        itr != itr_end;
        ++itr
      )
      {
        // true if extents of assemble::SSEGeometryInterface denoted by "itr" agree within tolerance with extents of "SSE_CRITERIA"
        if
        (
          ExtentsWithinTolerance
          (
            linal::Vector3D
            (
              ( *itr)->GetExtent( coord::GetAxes().e_X),
              ( *itr)->GetExtent( coord::GetAxes().e_Y),
              ( *itr)->GetExtent( coord::GetAxes().e_Z)
            ), criteria_extents
          )
        )
        {
          // add the body to "eligible_bodies"
          eligible_bodies.PushBack( *itr);
        }
      }

      // inform the user about the number of bodies that agree within tolerance with extents of "SSE_CRITERIA"
      BCL_MessageDbg( "BODIES.GetSize() " + util::Format()( BODIES.GetSize()));
      BCL_MessageDbg
      (
        "eligible_bodies.GetSize() " + util::Format()( eligible_bodies.GetSize())
      );

      //! if there are no bodes which can accommodate "SSE"
      if( eligible_bodies.IsEmpty())
      {
        // return empty ShPtr
        return util::ShPtr< assemble::SSEGeometryInterface>();
      }
      //! exactly one body can accommodate "SSE"
      else if( eligible_bodies.GetSize() == 1)
      {
        return eligible_bodies( 0);
      }
      // more than one body can accommodate "SSE"
      else
      {
        return PickBodyRandom().Pick( eligible_bodies);
      }
    }

    //! @brief ExtentsWithinTolerance checks to see if the extent of one Vector3D can accommodate the extent of a
    //!        second Vector3D with the tolerances added to the first Vector3D
    //! @param BODY_EXTENTS first Vector3D which is checked to see if it can accommodate the second Vector3D
    //! @param CRITERIA_EXTENTS second Vector3D which is checked to see if it can fit in RESTRAINT_EXTENTS
    //! @return returns a bool if CRITERIA_EXTENTS fits into BODY_EXTENTS given the member variable tolerances
    bool PickBodyExtent::ExtentsWithinTolerance
    (
      const linal::Vector3D &BODY_EXTENTS, const linal::Vector3D &CRITERIA_EXTENTS
    ) const
    {
      // true if "BODY_EXTENTS" match "CRITERIA_EXTENTS" within the allowed tolerances of
      // "m_ExtentUpperBoundTolerance" and "m_ExtentLowerBoundTolerance"
      if
      (
        CRITERIA_EXTENTS.X() + m_ExtentUpperBoundTolerance.X() >= BODY_EXTENTS.X()
        && CRITERIA_EXTENTS.X() - m_ExtentLowerBoundTolerance.X() <= BODY_EXTENTS.X()
        && CRITERIA_EXTENTS.Y() + m_ExtentUpperBoundTolerance.Y() >= BODY_EXTENTS.Y()
        && CRITERIA_EXTENTS.Y() - m_ExtentLowerBoundTolerance.Y() <= BODY_EXTENTS.Y()
        && CRITERIA_EXTENTS.Z() + m_ExtentUpperBoundTolerance.Z() >= BODY_EXTENTS.Z()
        && CRITERIA_EXTENTS.Z() - m_ExtentLowerBoundTolerance.Z() <= BODY_EXTENTS.Z()
      )
      {
        return true;
      }
      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PickBodyExtent::Read( std::istream &ISTREAM)
    {
      // read in members
      io::Serialize::Read( m_ExtentUpperBoundTolerance, ISTREAM);
      io::Serialize::Read( m_ExtentLowerBoundTolerance, ISTREAM);

      // end
      return ISTREAM;
    }

    //! @brief read from std::ostream
    //! @param OSTREAM input stream
    //! @return ostream which was read from
    std::ostream &PickBodyExtent::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write members
      io::Serialize::Write( m_ExtentUpperBoundTolerance, OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_ExtentLowerBoundTolerance, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

  } // namespace find
} // namespace bcl
