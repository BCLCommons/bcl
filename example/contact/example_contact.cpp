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
#include "contact/bcl_contact.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_contact.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleContact :
    public ExampleInterface
  {
  public:

    ExampleContact *Clone() const
    {
      return new ExampleContact( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

      // check if get namespace identifier is working
      BCL_MessageStd
      (
        "this is the bcl::contact namespace identifier: " + contact::GetNamespaceIdentifier()
      );
      BCL_Example_Check
      (
        contact::GetNamespaceIdentifier() == "bcl::contact",
        "namespace identifier for contact is incorrect!"
      );

      // initialize ranges
      math::Range< size_t> range_default( 6, std::numeric_limits< size_t>::max());
      math::Range< size_t> range_short( 6, 11);
      math::Range< size_t> range_mid( 12, 23);
      math::Range< size_t> range_long( 24, std::numeric_limits< size_t>::max());

      // test GetDefaultSequenceSeparationRange()
      BCL_ExampleCheck( contact::GetDefaultSequenceSeparationRange(), range_default);

      // test GetDefaultSequenceSeparationShortRange()
      BCL_ExampleCheck( contact::GetDefaultSequenceSeparationShortRange(), range_short);

      // test GetDefaultSequenceSeparationMidRange()
      BCL_ExampleCheck( contact::GetDefaultSequenceSeparationMidRange(), range_mid);

      // test GetDefaultSequenceSeparationLongRange()
      BCL_ExampleCheck( contact::GetDefaultSequenceSeparationLongRange(), range_long);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleContact

  const ExampleClass::EnumType ExampleContact::s_Instance
  (
    GetExamples().AddEnum( ExampleContact())
  );

} // namespace bcl

