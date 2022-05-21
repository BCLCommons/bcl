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
#include "descriptor/bcl_descriptor_window.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"
#include "linal/bcl_linal_vector_reference.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_window.cpp
  //!
  //! @author mendenjl
  //! @date Jan 31, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorWindow :
    public ExampleInterface
  {
  public:

    ExampleDescriptorWindow *Clone() const
    {
      return new ExampleDescriptorWindow( *this);
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
    //! @class Character
    //! @brief A class that returns a character, given a character
    //!
    //! @author mendenjl
    //! @date Jan 31, 2013
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class Character :
      public descriptor::BaseElement< char, char>
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
      virtual Character *Clone() const
      {
        return new Character( *this);
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
        static const std::string s_name( "Character");
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
        linal::VectorReference< char> &STORAGE
      )
      {
        STORAGE( 0) = *ELEMENT;
      }

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const
      {
        return io::Serializer();
      }

    }; // class AAIdentifier

    int Run() const
    {
    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // retrieve a window of 2 characters around a string.  Use the reflecting instance so that we don't have to worry
      // about undefined characters
      util::Implementation< descriptor::Base< char, char> > window_char( "ReflectingWindow(size=2,Character)");

      BCL_ExampleAssert( window_char.IsDefined(), true);

      BCL_ExampleCheck( window_char->GetSizeOfFeatures(), 5);

      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( window_char, "abcdefg"),
        "abcbc; babcd; cbade; dcbef; edcfg; fedgf; gfefe; "
      );

      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( window_char, "abc"),
        "abcbc; babcb; cbaba; "
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorWindow

  const ExampleClass::EnumType ExampleDescriptorWindow::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorWindow())
  );

  const util::SiPtr< const util::ObjectInterface> ExampleDescriptorWindow::Character::s_Instance
  (
    util::Enumerated< descriptor::Base< char, char> >::AddInstance( new ExampleDescriptorWindow::Character())
  );
} // namespace bcl
