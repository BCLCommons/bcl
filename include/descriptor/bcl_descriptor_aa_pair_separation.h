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

#ifndef BCL_DESCRIPTOR_AA_PAIR_SEPARATION_H_
#define BCL_DESCRIPTOR_AA_PAIR_SEPARATION_H_

// include the namespace header
#include "bcl_descriptor.h"

// other forward includes

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_pair.h"
#include "bcl_descriptor_type.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAPairSeparation
    //! @brief AAPairSeparations the output of separation descriptors into one vector based on sequence ID
    //!
    //! @see @link example_descriptor_aa_pair_separation.cpp @endlink
    //! @author teixeipl, mendenjl
    //! @date Feb 04, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAPairSeparation :
      public BasePair< biol::AABase, float>
    {
    private:

    //////////
    // data //
    //////////

    public:

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief virtual copy constructor
      AAPairSeparation *Clone() const;

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

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    private:

      //! @brief calculate the descriptors
      //! @param ELEMENT_A, ELEMENT_B: the element pair of interest
      //! @param STORAGE storage for the descriptor
      virtual void Calculate
      (
        const iterate::Generic< const biol::AABase> &ELEMENT_A,
        const iterate::Generic< const biol::AABase> &ELEMENT_B,
        linal::VectorReference< float> &STORAGE
      );

    }; // class AAPairSeparation

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_PAIR_SEPARATION_H_
