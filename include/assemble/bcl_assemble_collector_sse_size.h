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

#ifndef BCL_ASSEMBLE_COLLECTOR_SSE_SIZE_H_
#define BCL_ASSEMBLE_COLLECTOR_SSE_SIZE_H_

// include the namespace header
#include "bcl_assemble.h"

// include other forward headers - sorted alphabetically
#include "biol/bcl_biol.fwd.hh"
#include "util/bcl_util.fwd.hh"

// includes from bcl - sorted alphabetically
#include "find/bcl_find_collector_interface.h"
#include "math/bcl_math_range.h"
#include "storage/bcl_storage_map.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class CollectorSSESize
    //! @brief a collector, that collects all sses in a Domain that have a size within a given Range
    //! @details A range size_t is supplied to the constructor, and the operator of the CollectorInterface implementation
    //!          returns a list of SiPtr to the sses that have that size
    //!
    //! @see @link example_assemble_collector_sse_size.cpp @endlink
    //! @author woetzen
    //! @date Jun 21, 2011
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API CollectorSSESize :
      public find::CollectorInterface< util::SiPtrList< const SSE>, DomainInterface>
    {

    private:

    //////////
    // data //
    //////////

      //! size range for each SSType in which sse is collected from domain
      storage::Map< biol::SSType, math::Range< size_t> > m_SizeRangeMap;

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from sse size range, applies to all SSTypes
      //! @param SIZE_RANGE size range for sse sequence length
      CollectorSSESize( const math::Range< size_t> &SIZE_RANGE);

      //! @brief map of SSTypes and corresponding size ranges
      //! @param SIZE_RANGE_MAP size range for sse sequence length
      CollectorSSESize( const storage::Map< biol::SSType, math::Range< size_t> > &SIZE_RANGE_MAP);

      //! @brief Clone function
      //! @return pointer to new CollectorSSESize
      CollectorSSESize *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the size range
      //! @return reference to the size range
      const storage::Map< biol::SSType, math::Range< size_t> > &GetRanges() const
      {
        return m_SizeRangeMap;
      }

      //! @brief get the name of the object when used in a dynamic context
      //! @return the name of the object when used in a dynamic context
      const std::string &GetAlias() const;

      //! @brief returns parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief collect all sses from domain that fulfill the desired size range
      //! @param SSE_DOMAIN domain to collect from
      //! @return sses is a SiPtrList that have a sequence length according to range
      util::SiPtrList< const SSE> Collect( const DomainInterface &SSE_DOMAIN) const;

    //////////////////////
    // input and output //
    //////////////////////

    }; // class CollectorSSESize

  } // namespace assemble
} // namespace bcl

#endif // BCL_ASSEMBLE_COLLECTOR_SSE_SIZE_H_ 
