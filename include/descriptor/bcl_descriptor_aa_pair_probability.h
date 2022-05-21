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

#ifndef BCL_DESCRIPTOR_AA_PAIR_PROBABILITY_H_
#define BCL_DESCRIPTOR_AA_PAIR_PROBABILITY_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically
#include "score/bcl_score.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_pair.h"
#include "biol/bcl_biol_aa_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAPairProbability
    //! @brief generates the descriptors using the given matrix by accessing it by AAType pair.
    //! @details This class can be used with a distance probability matrix, as well as BLOSUM, PHAM, PHAT matrices
    //! as well as any other matrix.
    //!
    //! @see @link example_descriptor_aa_pair_probability.cpp @endlink
    //! @author karakam, mendenjl
    //! @date Apr 24, 2009
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API AAPairProbability :
      public BasePair< biol::AABase, float>
    {

    private:

    //////////
    // data //
    //////////

      //! AAAssignment score
      score::AAAssignment m_Score;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new AAPairProbability
      AAPairProbability *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

    ////////////////
    // operations //
    ////////////////

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
      void Calculate
      (
        const iterate::Generic< const biol::AABase> &ELEMENT_A,
        const iterate::Generic< const biol::AABase> &ELEMENT_B,
        linal::VectorReference< float> &STORAGE
      );

    }; // class AAPairProbability

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_PAIR_PROBABILITY_H_
