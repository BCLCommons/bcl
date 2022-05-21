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

#ifndef BCL_DESCRIPTOR_CACHE_PREFERENCE_H_
#define BCL_DESCRIPTOR_CACHE_PREFERENCE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "storage/bcl_storage_map.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_object_data_label.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @file bcl_descriptor_cache_preference.h
    //! @brief enum used by descriptors to declare whether they should be cached
    //!
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date Dec 11, 2012
    //!
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ///////////
  // enums //
  ///////////

    enum CachePreference
    {
      e_PreferCache,    //!< Try to retrieve a cached descriptor, recalculate only if necessary
      e_NeverCache,     //!< implies ignore cache, but also means that any derived descriptor of this should not be cached
                        //!< This should be used for descriptors that depend on a mutable property of the sequence in question
                        //!< Presently, for molecules, this should be used for an descriptor that depends on rotation/
                        //!< translation of the molecule. For protein sequences, it need not be used since we do not
                        //!< use descriptors and modify the protein in the same application
                        //!< Protein mutation classes copy the protein before doing the mutation as performance is not a
                        //!< concern there, so this cache preference doesn't need to be used for them
      e_IgnoreCache     //!< Default: Always recalculate the descriptor. Used for fast-to-compute descriptors except
                        //!< where e_NeverCache is appropriate
    };

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_CACHE_PREFERENCE_H_

