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

#ifndef BCL_DESCRIPTOR_AA_DISTANCE_FROM_AA_TYPE_H_
#define BCL_DESCRIPTOR_AA_DISTANCE_FROM_AA_TYPE_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes
#include "model/bcl_model.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "bcl_descriptor_type.h"
#include "bcl_descriptor_window_alignment_type.h"
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
    //! @class AADistanceFromAAType
    //! @brief AADistanceFromAAType gives the normalized position of an amino acid within the membrane the values are from
    //! (0-1.0) starting at inner membrane going to outer membrane. Based on Octopus TM predictions i = 0 o = 1.0
    //!
    //! @see @link example_descriptor_aa_distance_from_aa_type.cpp @endlink
    //! @author mendenjl
    //! @date May 06, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AADistanceFromAAType :
      public BaseElement< biol::AABase, float>
    {
    public:

    //////////
    // data //
    //////////

    private:

      //! Set of SeqIDs for the current protein with the specified AA type
      storage::Set< int> m_SeqIDsOfType;

      //! bool of whether the type has been located in the protein yet
      bool m_HaveLocated;

      //! AA type of interest
      biol::AAType m_Type;

      //! whether to use BLAST probability (True) or true type (False)
      bool m_UseBlastProbability;

      //! max distance to consider
      size_t m_MaxDistance;

      //! Alignment type (whether to look both directions or only left or right)
      WindowAlignmentEnum m_Direction;

      //! Maximum # of different AAs to find, e.g. if set to 1, only find the first AA with that type, for 2 find the
      //! first two of that type, etc.
      size_t m_MaxToFind;

      //! instances of the class
      static const util::SiPtr< const util::ObjectInterface> s_BlastInstance;
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AADistanceFromAAType( const bool &USE_BLAST_PROBABILITY);

      //! @brief virtual copy constructor
      AADistanceFromAAType *Clone() const;

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

      //! @brief function to find the aas of the given type in the specified sequence
      void FindAAs();

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      virtual void SetObjectHook();

    protected:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;
    }; // class AADistanceFromAAType

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_DISTANCE_FROM_AA_TYPE_H_
