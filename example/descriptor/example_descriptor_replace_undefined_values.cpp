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
#include "descriptor/bcl_descriptor_replace_undefined_values.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_base_element.h"
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_replace_undefined_values.cpp
  //!
  //! @author mendenjl
  //! @date Feb 06, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorReplaceUndefinedValues :
    public ExampleInterface
  {
  public:

    ExampleDescriptorReplaceUndefinedValues *Clone() const
    {
      return new ExampleDescriptorReplaceUndefinedValues( *this);
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
    //! @class AlphabeticNumber
    //! @brief returns a character, converted to its alphabetic number (A=1,Z=26); undefined for all non-alpha characters
    //!
    //! @author mendenjl
    //! @date Feb 06, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class AlphabeticNumber :
      public descriptor::BaseElement< char, float>
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
      AlphabeticNumber *Clone() const
      {
        return new AlphabeticNumber( *this);
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
        static const std::string s_name( "AlphabeticNumber");
        return s_name;
      }

      //! @brief get the feature size under the normal dimension setting (e.g. GetNormalDimension())
      //! @return the feature size, assuming this feature has its normal dimension setting
      size_t GetNormalSizeOfFeatures() const
      {
        return 1;
      }

    /////////////////
    // data access //
    /////////////////

    protected:

      //! @brief calculate the descriptors
      //! @param ELEMENT the element of the sequence of interest
      //! @param STORAGE storage for the descriptor
      void Calculate
      (
        const iterate::Generic< const char> &ELEMENT,
        linal::VectorReference< float> &STORAGE
      )
      {
        if( isalpha( *ELEMENT))
        {
          STORAGE( 0) = float( toupper( *ELEMENT) - int( 'A') + 1);
        }
        else
        {
          STORAGE( 0) = util::GetUndefined< float>();
        }
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        return io::Serializer();
      }

    }; // class AlphabeticNumber

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // create a descriptor that will return undefined values for some numbers
      util::Implementation< descriptor::Base< char, float> > alphabetic_number( "AlphabeticNumber");

      // create a descriptor that always returns 0
      util::Implementation< descriptor::Base< char, float> > zero( "Constant(0)");

      // create a hybrid descriptor that replaces undefined values with a different descriptor
      descriptor::ReplaceUndefinedValues< char> alphabetic_number_else_zero( alphabetic_number, zero);

      // test that the descriptors created are defined
      BCL_ExampleAssert( zero.IsDefined() && alphabetic_number.IsDefined(), true);

      // test size of features; should be identical to the size of the core descriptor
      BCL_ExampleAssert( alphabetic_number_else_zero.GetSizeOfFeatures(), alphabetic_number->GetSizeOfFeatures());

      // test that the descriptor returns the expected values
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( alphabetic_number, "a.b", 1),
        "1.0 ; nan ; 2.0 ; "
      );

      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations
        (
          util::Implementation< descriptor::Base< char, float> >( alphabetic_number_else_zero),
          "a.b  ",
          1
        ),
        "1.0 ; 0.0 ; 2.0 ; 0.0 ; 0.0 ; "
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorReplaceUndefinedValues

  const ExampleClass::EnumType ExampleDescriptorReplaceUndefinedValues::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorReplaceUndefinedValues())
  );

  const util::SiPtr< const util::ObjectInterface> ExampleDescriptorReplaceUndefinedValues::AlphabeticNumber::s_Instance
  (
    util::Enumerated< descriptor::Base< char, float> >::AddInstance( new ExampleDescriptorReplaceUndefinedValues::AlphabeticNumber())
  );
} // namespace bcl
