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

// include example header
#include "example.h"
// include the header of the class which this example is for
#include "descriptor/bcl_descriptor_window_pair_square.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_window_pair_square.cpp
  //!
  //! @author mendenjl
  //! @date Feb 04, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorWindowPairSquare :
    public ExampleInterface
  {
  public:

    ExampleDescriptorWindowPairSquare *Clone() const
    {
      return new ExampleDescriptorWindowPairSquare( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class Characters
    //! @brief A class that returns the pair of characters that are given
    //!
    //! @author mendenjl
    //! @date Feb 04, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class Characters :
      public descriptor::BasePair< char, char>
    {
    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief Clone function
      //! @return pointer to new BaseElement
      virtual Characters *Clone() const
      {
        return new Characters( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

      //! @brief return the data label
      //! @return data label as string
      const std::string &GetAlias() const
      {
        static const std::string s_name( "Characters");
        return s_name;
      }

      //! @brief return the type of symmetry this descriptor has
      //! @return e.g. symmetric if this descriptor returns the same value for A-B as B-A, Asymmetric otherwise
      descriptor::Type::Symmetry GetSymmetry() const
      {
        return descriptor::Type::e_Asymmetric;
      }

      //! @brief return whether this descriptor is valid if repeated elements are given
      //! @return true if this descriptor is valid if repeated elements are given
      //! This will be the case if the descriptor may have a legitimate value for A-A
      bool ConsiderRepeatedElements() const
      {
        return true;
      }

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return 2;
      }

    /////////////////
    // data access //
    /////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      //! @return true, if the recalculation was completed
      void Calculate
      (
        const iterate::Generic< const char> &ELEMENT_A,
        const iterate::Generic< const char> &ELEMENT_B,
        linal::VectorReference< char> &STORAGE
      )
      {
        STORAGE( 0) = *ELEMENT_A;
        STORAGE( 1) = *ELEMENT_B;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        return io::Serializer();
      }

    }; // class Character

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // retrieve a window of 2 characters around a string.
      util::Implementation< descriptor::Base< char, char> > window_char( "WindowPairSquare(size=1,Characters)");

      BCL_ExampleAssert( window_char.IsDefined(), true);

      BCL_ExampleCheck( window_char->GetSizeOfFeatures(), 18);

      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( window_char, "abc"),
        "       aaab babb  ;       aaabacbabbbc;       abac bbbc   ; "
        " aaab babb cacb   ; aaabacbabbbccacbcc; abac bbbc cbcc    ; "
        " babb cacb        ; babbbccacbcc      ; bbbc cbcc         ; "
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorWindowPairSquare

  const ExampleClass::EnumType ExampleDescriptorWindowPairSquare::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorWindowPairSquare())
  );

  const util::SiPtr< const util::ObjectInterface> ExampleDescriptorWindowPairSquare::Characters::s_Instance
  (
    util::Enumerated< descriptor::Base< char, char> >::AddInstance( new ExampleDescriptorWindowPairSquare::Characters())
  );
} // namespace bcl
