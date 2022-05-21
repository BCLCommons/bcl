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

#ifndef BCL_DESCRIPTOR_AA_SSE_INFO_H_
#define BCL_DESCRIPTOR_AA_SSE_INFO_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
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
    //! @class AASSEInfo
    //! @brief AASSEInfo gives the size of the in which this SSE resides, index of the SSE the given AA is in, the
    //! position (0-1.0) from N to C within the SSE, and the normalized distance from the center of the SSE
    //!
    //! @see @link example_descriptor_aa_sse_info.cpp @endlink
    //! @author teixeipl, mendenjl
    //! @date Mar 09, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AASSEInfo :
      public BaseElement< biol::AABase, float>
    {
    public:

    //////////
    // data //
    //////////

      //! custom SSE info type
      enum SSEInfoType
      {
        e_AASSESize,    //!< size of the SSE in which this AA is located
        e_AASSEID,   //!< ID of the SSE in which this AA is located
        e_AAPositionInSSE, //!< Position within the SSE of this AA
        e_AADistanceFromSSECenter, //!< Distance from the center of the SSE in which this AA is located
        s_NumberSSEInfoTypes
      };

      //! @brief SSEInfoType as string
      //! @param SSE_INFO_TYPE the message level
      //! @return the SSEInfoType as string
      static const std::string &GetSSEInfoTypeString( const SSEInfoType &SSE_INFO_TYPE);

      //! SSEInfoTypeEnum simplifies the usage of the SSEInfoType enum of this class
      typedef util::WrapperEnum< SSEInfoType, &GetSSEInfoTypeString, s_NumberSSEInfoTypes> SSEInfoTypeEnum;

    private:

      //! Vector containing the mapping from all seq IDs to indices within the SSE info vectors
      storage::Map< int, size_t> m_SeqIDToIndexMap;

      //! Vector containing the sizes for each amino acid's SSE
      storage::Vector< size_t> m_SSESizeVec;

      //! Vector containing the index of the SSE the AA is in
      storage::Vector< size_t> m_SSEIDVec;

      //! Vector containing the indexed position of each amino acid within its SSE starting from 0
      storage::Vector< size_t> m_PositionIndexVec;

      //! sse type info desired
      SSEInfoTypeEnum m_SSEInfoTypeEnum;

      //! sspred method
      sspred::Method m_Method;

      //! single instance of the three letter class
      static const util::SiPtr< const util::ObjectInterface> s_AASSESizeInstance;

      //! single instance of the three letter class
      static const util::SiPtr< const util::ObjectInterface> s_AASSEIDInstance;

      //! single instance of the three letter class
      static const util::SiPtr< const util::ObjectInterface> s_AAPositionInSSEInstance;

      //! single instance of the distance from SSE center class
      static const util::SiPtr< const util::ObjectInterface> s_AADistanceFromSSECenterInstance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AASSEInfo( const SSEInfoTypeEnum &SSE_INFO_TYPE_ENUM);

      //! @brief virtual copy constructor
      AASSEInfo *Clone() const;

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
      //! @param ELEMENT_A the element of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const biol::AABase> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const;

      //! @brief function to load files; should only be called the first time Calculate is called with a new sequence
      //! since the results are often in the cache
      void LoadFiles();

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook();

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;
    }; // class AASSEInfo

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_SSE_INFO_H_
