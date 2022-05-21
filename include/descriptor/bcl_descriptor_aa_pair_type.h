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

#ifndef BCL_DESCRIPTOR_AA_PAIR_TYPE_H_
#define BCL_DESCRIPTOR_AA_PAIR_TYPE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"
#include "bcl_descriptor_window_alignment_type.h"
#include "biol/bcl_biol_aa_base.h"
#include "iterate/bcl_iterate_generic.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AAPairType
    //! @brief Generates the BlastProfile of the given amino acid
    //! @details This class returns the 20-value blast profile for a given amino acid
    //!
    //! @see @link example_descriptor_aa_pair_type.cpp @endlink
    //! @author mendenjl
    //! @date Oct 25, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API AAPairType :
      public BaseElement< biol::AABase, float>
    {
    private:

      //! Stencil of residues (does not include the central residue)
      //! e.g. 1,4 means that probabilities will be calculated
      //! P(In Helix | Res at -4 Res at 0) * P(In Helix | Res at -1 Res at 0) * P(In Helix | Res at 4 Res at 0) * P(In Helix | Res at 1 Res at 0)
      linal::Vector< size_t> m_Stencil;

      //! Weights for each component of the stencil (defaults to 1)
      linal::Vector< float> m_StencilWeights;

      //! Alignment of stencil (e.g. left, right, or center)
      WindowAlignmentEnum m_Alignment;

      //! true if the AB is the same as BA
      bool m_Symmetric;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_AsymmetricInstance;
      static const util::SiPtr< const util::ObjectInterface> s_SymmetricInstance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      AAPairType( const bool &SYMMETRIC = false);

      //! @brief Clone function
      //! @return pointer to new AAPairType
      AAPairType *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

    protected:

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      virtual void Calculate
      (
        const iterate::Generic< const biol::AABase> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      );

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class AAPairType

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_PAIR_TYPE_H_
