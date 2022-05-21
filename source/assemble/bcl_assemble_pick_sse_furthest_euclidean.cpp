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
#include "assemble/bcl_assemble_pick_sse_furthest_euclidean.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> PickSSEFurthestEuclidean::s_Instance
    (
      GetObjectInstances().AddInstance( new PickSSEFurthestEuclidean())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief Clone is the virtual Clone constructor
    //! @return a pointer to new LocatorSSEFurthest which is a copy of this
    PickSSEFurthestEuclidean *PickSSEFurthestEuclidean::Clone() const
    {
      return new PickSSEFurthestEuclidean( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &PickSSEFurthestEuclidean::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief Picks the assemble::SSE object which is furthest away from "m_ReferencePoint"
    //! @param SSE_LIST is the SiPtrList which provides the pool of assemble::SSE to pick from
    //! @param REFERENCE_POINT is the point from which the furthest SSE will be determined
    //! @return returns SiPtr to the assemble::SSE object which is furthest away from "m_ReferencePoint"
    util::SiPtr< const SSE>
    PickSSEFurthestEuclidean::Pick
    (
      const util::SiPtrList< const SSE> &SSE_LIST,
      const linal::Vector3D &REFERENCE_POINT
    ) const
    {
      // create double "max_distance" for holding the max distance between "m_ReferencePoint" and a assemble::SSE
      double max_distance( 0); //< initialize to zero

      // create SiPtr< const assemble::SSE> "furthest_moveable"
      util::SiPtr< const SSE> furthest_moveable;

      // iterate through "GROUP"
      for
      (
        util::SiPtrList< const SSE>::const_iterator itr( SSE_LIST.Begin()), itr_end( SSE_LIST.End());
        itr != itr_end;
        ++itr
      )
      {
        // create double "current_distance" initialize to the distance between "m_ReferencePoint" and the center of
        // the assemble::SSE currently denoted by "itr"
        double current_distance = linal::Distance( REFERENCE_POINT, ( *itr)->GetCenter());

        // check if "current_distance" is greater than "max_distance"
        if( current_distance > max_distance)
        {
          // set "furthest_sse" to the assemble::SSE currently denoted by "itr"
          furthest_moveable = ( *itr);

          // set "max_distance" to "current_distance"
          max_distance = current_distance;
        }
      }

      // return "furthest_moveable"
      return furthest_moveable;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &PickSSEFurthestEuclidean::Read( std::istream &ISTREAM)
    {
      // end
      return ISTREAM;
    }

    //! @brief read from std::ostream
    //! @param OSTREAM input stream
    //! @param INDENT indentation
    //! @return ostream which was read from
    std::ostream &PickSSEFurthestEuclidean::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // end
      return OSTREAM;
    }

  } // namespace assemble
} // namespace bcl
