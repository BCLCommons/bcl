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
#include "storage/bcl_storage_vector.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  namespace storage
  {
    //! @brief create an index vector consisting of numbers [0, N - 1), represnting indices
    //! @param N the upper end of the index vector range
    //! @return vector containing 0,...,N-1
    Vector< size_t> CreateIndexVector( const size_t &N)
    {
      Vector< size_t> indices;
      indices.AllocateMemory( N);
      for( size_t i( 0); i < N; ++i)
      {
        indices.PushBack( i);
      }
      return indices;
    }

    //! @brief create an index vector consisting of consecutive integers in range [MIN, MAX), representing indices
    //! @param MIN the lower end of the range
    //! @param MAX the upper end of the range (exclusive)
    //! @return vector containing MIN,...,MAX-1
    Vector< size_t> CreateIndexVector( const size_t &MIN, const size_t &MAX)
    {
      BCL_Assert( MIN <= MAX, "CreateIndexVector called with MAX < MIN");
      Vector< size_t> indices;
      indices.AllocateMemory( MAX - MIN);
      for( size_t i( MIN); i < MAX; ++i)
      {
        indices.PushBack( i);
      }
      return indices;
    }

    template class BCL_API Vector< std::string>;
    template class BCL_API Vector< double>;
    template class BCL_API Vector< float>;
    template class BCL_API Vector< size_t>;
  } // namespace storage
} // namespace bcl
