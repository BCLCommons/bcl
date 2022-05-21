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
#include "descriptor/bcl_descriptor_aa_letter_code.h"

// includes from bcl - sorted alphabetically
#include "biol/bcl_biol_aa_base.h"
#include "util/bcl_util_enumerated.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace descriptor
  {

  //////////
  // data //
  //////////

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AALetterCode::s_OneInstance
    (
      util::Enumerated< Base< biol::AABase, char> >::AddInstance
      (
        new AALetterCode( true)
      )
    );

    // instantiate s_Instance
    const util::SiPtr< const util::ObjectInterface> AALetterCode::s_ThreeInstance
    (
      util::Enumerated< Base< biol::AABase, char> >::AddInstance
      (
        new AALetterCode( false)
      )
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief constructor from a list of properties
    //! @param USE_ONE_LETTER bool indicating whether one or three letter code will be used
    AALetterCode::AALetterCode
    (
      const bool &USE_ONE_LETTER
    ) :
      m_UseOneLetter( USE_ONE_LETTER)
    {
    }

    //! @brief Clone function
    //! @return pointer to new BaseElement
    AALetterCode *AALetterCode::Clone() const
    {
      return new AALetterCode( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name of the object behind a pointer or the current object
    //! @return the class name
    const std::string &AALetterCode::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief return the data label
    //! @return data label as string
    const std::string &AALetterCode::GetAlias() const
    {
      static const std::string s_one_name( "AAOneLetterCode");
      static const std::string s_three_name( "AAThreeLetterCode");
      return m_UseOneLetter ? s_one_name : s_three_name;
    }

    //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
    //! @return the feature size, assuming this feature has its normal dimension setting
    size_t AALetterCode::GetNormalSizeOfFeatures() const
    {
      return m_UseOneLetter ? 1 : 3;
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief calculate the descriptors
    //! @param ELEMENT the element of the sequence of interest
    //! @param STORAGE storage for the descriptor
    void AALetterCode::Calculate
    (
      const iterate::Generic< const biol::AABase> &ELEMENT,
      linal::VectorReference< char> &STORAGE
    )
    {
      if( m_UseOneLetter)
      {
        STORAGE( 0) = ELEMENT->GetData()->GetType()->GetOneLetterCode();
      }
      else
      {
        const std::string &three_letter_code( ELEMENT->GetData()->GetType()->GetThreeLetterCode());
        std::copy( three_letter_code.begin(), three_letter_code.end(), STORAGE.Begin());
      }
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer AALetterCode::GetSerializer() const
    {
      io::Serializer parameters;
      const std::string one_params( "Returns the amino acid's one letter code");
      const std::string three_params( "Returns the amino acid's three letter code");
      parameters.SetClassDescription( m_UseOneLetter ? one_params : three_params);

      return parameters;
    }

  } // namespace descriptor
} // namespace bcl
