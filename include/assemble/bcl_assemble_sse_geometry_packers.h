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

#ifndef BCL_ASSEMBLE_SSE_GEOMETRY_PACKERS_H_
#define BCL_ASSEMBLE_SSE_GEOMETRY_PACKERS_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_assemble_sse_geometry_packer_interface.h"
#include "bcl_assemble_sse_geometry_packing.h"
#include "math/bcl_math_binary_function_interface_serializable.h"
#include "storage/bcl_storage_list.h"
#include "util/bcl_util_enumerate.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEGeometryPackers
    //! @brief class enumerates different ways to calculate all packings between sub-geometries for a given pair of geometries
    //! @details This enumerator, enumerates different SSEGeometryPackerInterface derived classes with a return type
    //! of list of list of SSEGeometryPackings. For a given pair of GEOMETRY_A (composed of M sub-geometries) and
    //! GEOMETRY_B (composed of N sub-geometries), the return type is a list with M sublists with each having N packings.
    //! Each enumerator is repeated with a NoCache option, which makes sure even if the caching flag is provided
    //! that instance won't use caching.
    //!
    //! @see @link example_assemble_sse_geometry_packers.cpp @endlink
    //! @author karakam
    //! @date Jun 21, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEGeometryPackers :
      public util::Enumerate
             <
               util::ShPtr< SSEGeometryPackerInterface< storage::Vector< storage::List< SSEGeometryPacking> > > >,
               SSEGeometryPackers
             >
    {

      friend class
      util::Enumerate
      <
        util::ShPtr< SSEGeometryPackerInterface< storage::Vector< storage::List< SSEGeometryPacking> > > >,
        SSEGeometryPackers
      >;

    public:

    //////////
    // data //
    //////////

      const SSEGeometryPacker e_Default;          //!< default
      const SSEGeometryPacker e_Default_NoCache;  //!< default no cache
      const SSEGeometryPacker e_SSEClash;         //!< SSE clash
      const SSEGeometryPacker e_SSEClash_NoCache; //!< SSE clash no cache

    private:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief private constructor for constructing all SSEGeometryPackers
      SSEGeometryPackers();

    public:

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

    }; // class SSEGeometryFragmentsPacker

    //! @brief access static instance of SSEGeometryPackers enums
    //! @return static instance of SSEGeometryPackers enums
    BCL_API
    const SSEGeometryPackers &GetSSEGeometryPackers();

  } // namespace assemble

  namespace util
  {

  /////////////////////////////
  // explicit instantiations //
  /////////////////////////////

    BCL_EXPIMP_TEMPLATE template class BCL_API Enumerate< ShPtr< assemble::SSEGeometryPackerInterface< storage::Vector< storage::List< assemble::SSEGeometryPacking> > > >, assemble::SSEGeometryPackers>;

  } // namespace util
} // namespace bcl

#endif // BCL_ASSEMBLE_SSE_GEOMETRY_PACKERS_H_ 
