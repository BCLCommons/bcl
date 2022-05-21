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
#include "assemble/bcl_assemble_sse_factory_interface.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse_pool.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace assemble
  {

  ///////////////
  // operators //
  ///////////////

    //! @brief operator that returns a set of SSEs for the given vector of AASequences
    //! @brief SEQUENCES Vector of AASequences from which the SSEPool is going to be built
    //! @return SSEPool built from provided SEQUENCES
    SSEPool
    SSEFactoryInterface::operator()( const util::SiPtrVector< biol::AASequence> &SEQUENCES) const
    {
      // initialize SSEPool to be returned
      SSEPool sse_pool;

      // iterate over the provided sequences
      for
      (
        util::SiPtrVector< biol::AASequence>::const_iterator seq_itr( SEQUENCES.Begin()),
          seq_itr_end( SEQUENCES.End());
        seq_itr != seq_itr_end;
        ++seq_itr
      )
      {
        // create an ssepool from the sequence behind seq_itr
        const SSEPool this_sse_pool( operator()( **seq_itr));

        // append the pool for this sequence to sse_pool
        sse_pool.InsertElements( this_sse_pool.Begin(), this_sse_pool.End());
      }

      // end
      return sse_pool;
    }

  } // namespace assemble
} // namespace bcl
