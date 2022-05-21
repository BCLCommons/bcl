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

#ifndef BCL_DESCRIPTOR_SEQUENCE_SSE_COUNT_H_
#define BCL_DESCRIPTOR_SEQUENCE_SSE_COUNT_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "sspred/bcl_sspred.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_sequence.h"
#include "bcl_descriptor_type.h"
#include "biol/bcl_biol_aa_sequence.h"
#include "iterate/bcl_iterate_generic.h"
#include "storage/bcl_storage_vector.h"
#include "util/bcl_util_implementation.h"
#include "util/bcl_util_wrapper_enum.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class SequenceSSECount
    //! @brief SequenceSSECount gives the size of the in which this SSE resides, index of the SSE the given AA is in, the
    //! position (0-1.0) from N to C within the SSE, and the normalized distance from the center of the SSE
    //!
    //! @see @link example_descriptor_sequence_sse_count.cpp @endlink
    //! @author mendenjl
    //! @date Feb 13, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API SequenceSSECount :
      public BaseSequence< biol::AABase, float>
    {
    private:

    //////////
    // data //
    //////////

      //! min sse sizes for consideration
      size_t m_MinHelixSize;
      size_t m_MinStrandSize;
      size_t m_MinCoilSize;

      //! method to use to determine SSE info
      sspred::Method m_SSMethod;

      //! if true, only count SSEs in that are in the membrane
      bool m_OnlyMembraneCounts;

      //! if true, only count SSEs that are in the solution
      bool m_OnlySolutionCounts;

      //! method to use for defining the environment; used only if m_OnlyMembraneCounts or m_OnlySolutionCounts is true
      sspred::Method m_EnvironmentMethod;

      //! single instance of the three letter class
      static const util::SiPtr< const util::ObjectInterface> s_SSECountsInstance;

      //! single instance of the three letter class
      static const util::SiPtr< const util::ObjectInterface> s_MembraneSSECountsInstance;

      //! single instance of the three letter class
      static const util::SiPtr< const util::ObjectInterface> s_SolutionSSECountsInstance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      SequenceSSECount( const bool &ONLY_MEMBRANE_COUNTS, const bool &ONLY_SOLUTION_COUNTS);

      //! @brief virtual copy constructor
      SequenceSSECount *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

    //////////////////////
    // helper functions //
    //////////////////////

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

    private:

      //! @brief calculate the descriptors
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        linal::VectorReference< float> &STORAGE
      );

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const;

      //! @brief get the minimum sse size for a particular ss-type
      //! @param TYPE the type of sse to get the minimum sse size for
      size_t GetMinSSESize( const biol::SSType &TYPE) const;

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;
    }; // class SequenceSSECount

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_SEQUENCE_SSE_COUNT_H_
