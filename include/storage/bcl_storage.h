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

#ifndef BCL_STORAGE_H_
#define BCL_STORAGE_H_

// include the namespace forward header
#include "bcl_storage.fwd.hh"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @file bcl_storage.h
  //! @brief TODO: document
  //! @details TODO: document
  //!
  //! @see @link example_storage.cpp @endlink
  //! @author woetzen
  //! @date Jul 22, 2008
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  namespace storage
  {
    //! @brief identifier for the name space
    //! @return the name of the namespace
    BCL_API
    const std::string &GetNamespaceIdentifier();

    //! @brief advances the provided iterator for given times while checking it does not reach end
    //! @param ITERATOR SequenceContainer iterator that will be advanced
    //! @param END_ITERATOR end iterator of the SequenceContainer that ITERATOR belongs to
    //! @param N number of times the iterator has to be increased
    template< typename t_IteratorType>
    inline t_IteratorType &AdvanceIterator
    (
      t_IteratorType &ITERATOR,
      const t_IteratorType &END_ITERATOR,
      const size_t N
    )
    {
      // while counter is smaller than N and the END is not reached
      for( size_t counter( 0); counter < N && ITERATOR != END_ITERATOR; ++counter)
      {
        ++ITERATOR;
      }

      // end
      return ITERATOR;
    }

  } // namespace storage
} // namespace bcl

#endif // BCL_STORAGE_H_
