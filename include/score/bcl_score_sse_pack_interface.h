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

#ifndef BCL_SCORE_SSE_PACK_INTERFACE_H_
#define BCL_SCORE_SSE_PACK_INTERFACE_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "storage/bcl_storage.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_function_interface_serializable.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SSEPackInterface
    //! @brief is an interface class to be derived from in all score classes related to packing
    //! @details This class provides a interface for classes like SSEPairPacking and StrandPairing which score
    //! SSEGeometryPacking objects. In addition to operator() it also has a pure virtual function AreValidSSEs()
    //! that allows upper-level classes to check whether two SSEs are suitable for the specific score without without
    //! having to call the operator.
    //!
    //! @remarks example unnecessary
    //! @author karakam
    //! @date Apr 20, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SSEPackInterface :
      public math::FunctionInterfaceSerializable< assemble::SSEGeometryPacking, double>
    {

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to a new SSEPackInterface
      virtual SSEPackInterface *Clone() const = 0;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns the minimal interface length used for calculating packing
      //! @return the minimal interface length used for calculating packing
      virtual const double GetMinimalInterfaceLength() const = 0;

    ////////////////
    // operations //
    ////////////////

      //! @brief returns whether the given pair of SSEs are valid
      //! @param SSE_A first SSE of interest
      //! @param SSE_B first SSE of interest
      //! @return whether the given pair of SSEs are valid
      virtual bool AreValidSSEs( const assemble::SSE &SSE_A, const assemble::SSE &SSE_B) const = 0;

    }; // class SSEPackInterface

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_SSE_PACK_INTERFACE_H_
