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

#ifndef BCL_ASSEMBLE_PICK_SSE_FURTHEST_EUCLIDEAN_H_
#define BCL_ASSEMBLE_PICK_SSE_FURTHEST_EUCLIDEAN_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "linal/bcl_linal.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_pick_criteria_interface.h"
#include "util/bcl_util_si_ptr_list.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class PickSSEFurthestEuclidean
    //! @brief used for picking the single furthest item from a point of reference.
    //! @brief this class picks from a given list of SSEs, the one that is furthest from the given coordinate
    //!
    //! @see @link example_assemble_pick_sse_furthest_euclidean.cpp @endlink
    //! @author alexanns
    //! @date 03/27/2008
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API PickSSEFurthestEuclidean :
      public find::PickCriteriaInterface< util::SiPtr< const SSE>, util::SiPtrList< const SSE>, linal::Vector3D>
    {

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      PickSSEFurthestEuclidean()
      {
      }

      //! virtual copy constructor
      PickSSEFurthestEuclidean *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    ////////////////
    // operations //
    ////////////////

      //! Picks the assemble::SSE object which is furthest away from "m_ReferencePoint"
      //! @param SSE_LIST is the SiPtrList which provides the pool of assemble::SSE to pick from
      //! @param REFERENCE_POINT is the point from which the furthest SSE will be determined
      //! @return returns SiPtr to the assemble::SSE object which is furthest away from "m_ReferencePoint"
      util::SiPtr< const SSE>
      Pick( const util::SiPtrList< const SSE> &SSE_LIST, const linal::Vector3D &REFERENCE_POINT) const;

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
      //! @return ostream which was read from
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

    }; // class PickSSEFurthestEuclidean

  } // namespace assemble
} // namespace bcl

#endif //BCL_ASSEMBLE_PICK_SSE_FURTHEST_EUCLIDEAN_H_
