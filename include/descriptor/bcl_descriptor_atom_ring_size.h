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

#ifndef BCL_DESCRIPTOR_ATOM_RING_SIZE_H_
#define BCL_DESCRIPTOR_ATOM_RING_SIZE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AtomRingSize
    //! @brief calculates the # of rings with a given conjugation from the molecule
    //!
    //! @see @link example_descriptor_atom_ring_size.cpp @endlink
    //! @author mendenjl
    //! @date Apr 17, 2014
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AtomRingSize :
      public BaseElement< chemistry::AtomConformationalInterface, float>
    {

    private:

    //////////
    // data //
    //////////

      storage::Vector< size_t> m_AtomRingSizes; //!< Ring sizes in the current molecule
      bool m_PreferLargest; //!< If multiple rings are present for an atom, return the largest (false -> smallest)

    public:

      enum
      {
        s_ChainAtomsMaxRingSize = 2,    //!< Value returned if returning the largest ring and the atom is not in a ring
        s_ChainAtomsMinRingSize = 10000 //!< Value returned if returning the smallest ring and the atom is not in a ring
      };

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a bool
      //! @param LARGEST whether to prefer the larger of the rings that the atom is connected to,
      //!        if it is connected to multiple
      explicit AtomRingSize( const bool &LARGEST);

      //! @brief Clone function
      //! @return pointer to new AtomRingSize
      AtomRingSize *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief get name of the current class
      //! @return name of the class
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature siz  e under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return 1;
      }

    ////////////////
    // operations //
    ////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const chemistry::AtomConformationalInterface> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

      //! @brief hook that derived classes can override to add behavior after every time SetObject is called
      void SetObjectHook();

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

      //! @brief get the cache preference under the normal dimension setting (e.g. GetType().GetDimension())
      //! @return the cache preference, assuming this feature has its normal dimension setting
      CachePreference GetNormalCachePreference() const
      {
        return e_PreferCache;
      }

    }; // class AtomRingSize

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_ATOM_RING_SIZE_H_
