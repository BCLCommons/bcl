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

#ifndef BCL_DESCRIPTOR_AA_LETTER_CODE_H_
#define BCL_DESCRIPTOR_AA_LETTER_CODE_H_

// include the namespace header
#include "bcl_descriptor.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_descriptor_base_element.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class AALetterCode
    //! @brief A class that returns the correct letter code as a descriptor for identification (one or three letters)
    //!
    //! @remarks
    //! @see @link example_descriptor_aa_letter_code.cpp @endlink
    //! @author teixeipl, mendenjl
    //! @date Feb 23, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API AALetterCode :
      public BaseElement< biol::AABase, char>
    {
    public:

    //////////
    // data //
    //////////

      //! single instance of the one letter class
      static const util::SiPtr< const util::ObjectInterface> s_OneInstance;

      //! single instance of the three letter class
      static const util::SiPtr< const util::ObjectInterface> s_ThreeInstance;

    private:
      //! bool that stores which type of the letter code will be used, one or three
      bool m_UseOneLetter;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a list of properties
      //! @param USE_ONE_LETTER bool indicating whether one or three letter code will be used
      AALetterCode( const bool &USE_ONE_LETTER);

      //! @brief Clone function
      //! @return pointer to new BaseElement
      virtual AALetterCode *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const;

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const;

    /////////////////
    // data access //
    /////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      virtual void Calculate
      (
        const iterate::Generic< const biol::AABase> &ELEMENT,
        linal::VectorReference< char> &STORAGE
      );

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    }; // class AALetterCode

  } // namespace descriptor
} // namespace bcl

#endif // BCL_DESCRIPTOR_AA_LETTER_CODE_H_
