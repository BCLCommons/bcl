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
#include "restraint/bcl_restraint_contains_body_origin.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace util
  {

    template class BinaryFunctionInterface< assemble::SSEGeometryInterface, assemble::SSE, bool>;

  } // namespace util

  namespace restraint
  {
  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> ContainsBodyOrigin::s_Instance
    (
      GetObjectInstances().AddInstance( new ContainsBodyOrigin())
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    ContainsBodyOrigin::ContainsBodyOrigin()
    {
    }

    //! @brief virtual copy constructor
    ContainsBodyOrigin *ContainsBodyOrigin::Clone() const
    {
      return new ContainsBodyOrigin( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ContainsBodyOrigin::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() determines if the BODY is occupied by any of the geometries of the SSE
    //! @details this is done by checking if any of the fragments in the body will contain the center of any of the
    //!          fragments in the sse
    //! @param BODY representing the density rod
    //! @param SSE represents the secondary structure element which contains geometries/fragments
    //! @return returns true if BODY contains any of the SSE's geometries
    bool ContainsBodyOrigin::operator()( const assemble::SSEGeometryInterface &BODY, const assemble::SSE &SSE) const
    {
      // create local copy of fragments of body (this is necessary because GetGeometries() returns only local copy of fragments)
      const util::SiPtrVector< const assemble::SSEGeometryInterface> body_geometries( BODY.GetSSEGeometries());

      // create local copy of fragments of sse
      const util::SiPtrVector< const assemble::SSEGeometryInterface> sse_geometries( SSE.GetSSEGeometries());

      // iterate over the fragments of BODY
      for
      (
        util::SiPtrVector< const assemble::SSEGeometryInterface>::const_iterator
          body_itr( body_geometries.Begin()), body_itr_end( body_geometries.End());
        body_itr != body_itr_end; ++body_itr
      )
      {
        // iterate over sse fragments
        for
        (
          util::SiPtrVector< const assemble::SSEGeometryInterface>::const_iterator
            sse_itr( sse_geometries.Begin()), sse_itr_end( sse_geometries.End());
          sse_itr != sse_itr_end; ++sse_itr
        )
        {
          // check whether origin of sse fragment is within the body_fragment
          if( ContainsBodyOriginCheck( **body_itr, ( *sse_itr)->GetCenter()))
          {
            return true;
          }
        }
      }

      // non of the sse fragments was found to be within one of the body fragments
      return false;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief read distance from std::istream
    //! @param ISTREAM input stream
    //! @return istream which was read from
    std::istream &ContainsBodyOrigin::Read( std::istream &ISTREAM)
    {
      return ISTREAM;
    }

    //! @brief write distance to std::ostream
    //! @param OSTREAM output stream to write to
    //! @param INDENT indentation
    //! @return output stream which was written to
    std::ostream &ContainsBodyOrigin::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      return OSTREAM;
    }

  //////////////////////
  // helper functions //
  //////////////////////

    //! @brief ContainsBodyOriginCheck determines if the coord::GeometryInterface envelops the origin
    //! @param BODY coord::GeometryInterface which is the body involved
    //! @param ORIGIN origin for which is checked whether it is within body
    //! @return returns true if Body envelops the origin
    bool ContainsBodyOrigin::ContainsBodyOriginCheck( const assemble::SSEGeometryInterface &BODY, const linal::Vector3D &ORIGIN) const
    {
      // calculate the distance of the origin from the z-axis of the body
      const double dist_origin_from_zaxis
      (
        linal::CalculateDistancePointFromLine( ORIGIN, BODY.GetCenter(), BODY.GetAxis( coord::GetAxes().e_Z))
      );

      // calculate the footpoint of the origin on the z-axis of the body
      const linal::Vector3D footpoint_origin_on_zaxis
      (
        linal::CalculateFootpoint( ORIGIN, BODY.GetCenter(), BODY.GetAxis( coord::GetAxes().e_Z))
      );

      // return true if distance of origin to footpoint is less than the sse radius and if distance of footpoint
      // to center of body is less than the z extent
      return
      (
        dist_origin_from_zaxis <= BODY.GetExtent( coord::GetAxes().e_X) &&
        linal::Distance( footpoint_origin_on_zaxis, BODY.GetCenter()) <= BODY.GetExtent( coord::GetAxes().e_Z)
      );
    }

  } // namespace restraint
} // namespace bcl
