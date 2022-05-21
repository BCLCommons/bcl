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
#include "contact/bcl_contact.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace contact
  {

    //! @brief identifier for the name space
    //! @return the name of the namespace
    const std::string &GetNamespaceIdentifier()
    {
      static const std::string *s_namespace_name( new std::string( util::ExtractNamespaceIdentifier( __PRETTY_FUNCTION__)));
      return *s_namespace_name;
    }

    //! @brief return default sequence separation range for contacts
    //! @return default sequence separation range contacts
    const math::Range< size_t> &GetDefaultSequenceSeparationRange()
    {
      // initialize static range
      static const math::Range< size_t> s_range( 6, std::numeric_limits< size_t>::max());

      // end
      return s_range;
    }

    //! @brief return sequence separation range for short-range contacts
    //! @return sequence separation range for short-range contacts
    const math::Range< size_t> &GetDefaultSequenceSeparationShortRange()
    {
      // initialize static range
      static const math::Range< size_t> s_range_short( 6, 11);

      // end
      return s_range_short;
    }

    //! @brief return sequence separation range for mid-range contacts
    //! @return sequence separation range for mid-range contacts
    const math::Range< size_t> &GetDefaultSequenceSeparationMidRange()
    {
      // initialize static range
      static const math::Range< size_t> s_range_mid( 12, 23);

      // end
      return s_range_mid;
    }

    //! @brief return sequence separation range for long-range contacts
    //! @return sequence separation range for long-range contacts
    const math::Range< size_t> &GetDefaultSequenceSeparationLongRange()
    {
      // initialize static range
      static const math::Range< size_t> s_range_long( 24, std::numeric_limits< size_t>::max());

      // end
      return s_range_long;
    }

  } // namespace contact
} // namespace bcl
