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

#ifndef BCL_FIND_PICK_BODY_EXTENT_H_
#define BCL_FIND_PICK_BODY_EXTENT_H_

// include the namespace header
#include "bcl_find.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_find_pick_criteria_interface.h"
#include "linal/bcl_linal_vector_3d.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace find
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PickBodyExtent
    //! @brief to get assemble::SSEGeometryInterface out of a restraint::Body which fulfills an extent requirement
    //! @details For specific example, you want to place an SSE into a assemble::SSEGeometryInterface which has approximately the same extent
    //! in z-axis, so you want to pick the assemble::SSEGeometryInterface out of restraint::Body which has an extent matching the extent
    //! of the SSE within some tolerance.
    //!
    //! @see @link example_find_pick_body_extent.cpp @endlink
    //! @author alexanns
    //! @date 06/30/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PickBodyExtent :
      public PickCriteriaInterface< util::ShPtr< assemble::SSEGeometryInterface>, util::ShPtrVector< assemble::SSEGeometryInterface>, assemble::SSEGeometryInterface>
    {

    /////////////
    // friends //
    /////////////

    private:

    //////////
    // data //
    //////////

      //! Vector3D "m_ExtentUpperBoundTolerance" holds the tolerance for the upper bound of the extents
      linal::Vector3D m_ExtentUpperBoundTolerance;

      //! Vector3D "m_ExtentLowerBoundTolerance" holds the tolerance for the lower bound of the extents
      linal::Vector3D m_ExtentLowerBoundTolerance;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PickBodyExtent();

      //! @brief construct from two linal::Vector3D which describe the upper and lower bound extent tolerance
      PickBodyExtent( const linal::Vector3D &UPPER_TOLERANCE, const linal::Vector3D &LOWER_TOLERANCE);

      //! @brief virtual copy constructor
      PickBodyExtent *Clone() const;

      //! virtual destructor
      ~PickBodyExtent();

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! Pick the util::ShPtr< assemble::SSEGeometryInterface> object from restraint::Body based on an SSE
      //! @param BODIES is the object which will provide the pool of util::ShPtr< assemble::SSEGeometryInterface> to pick from
      //! @param CRITERIA the assemble::SSEGeometryInterface which provides the extent criteria by which a assemble::SSEGeometryInterface is chosen
      //! @return return a util::ShPtr< assemble::SSEGeometryInterface> which can accomodate "CRITERIA"
      util::ShPtr< assemble::SSEGeometryInterface>
      Pick( const util::ShPtrVector< assemble::SSEGeometryInterface> &BODIES, const assemble::SSEGeometryInterface &CRITERIA) const;

      //! @brief ExtentsWithinTolerance checks to see if the extent of one Vector3D can accomodate the extent of a
      //!        second Vector3D with the tolerances added to the first Vector3D
      //! @param BODY_EXTENTS first Vector3D which is checked to see if it can accomodate the second Vector3D
      //! @param CRITERIA_EXTENTS second Vector3D which is checked to see if it can fit in RESTRAINT_EXTENTS
      //! @return returns a bool if CRITERIA_EXTENTS fits into BODY_EXTENTS given the member variable tolerances
      bool ExtentsWithinTolerance( const linal::Vector3D &BODY_EXTENTS, const linal::Vector3D &CRITERIA_EXTENTS) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM);

      //! @brief read from std::ostream
      //! @param OSTREAM input stream
      //! @param INDENT indentation
      //! @return ostream which was read from
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // template class PickBodyExtent

  } // namespace find
} // namespace bcl

#endif //BCL_FIND_PICK_BODY_EXTENT_H_
