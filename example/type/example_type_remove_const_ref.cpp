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
#include "type/bcl_type_remove_const_ref.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_type_remove_const_ref.cpp
  //!
  //! @author mendenjl
  //! @date Nov 29, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleTypeRemoveConstRef :
    public ExampleInterface
  {
  public:

    ExampleTypeRemoveConstRef *Clone() const
    { return new ExampleTypeRemoveConstRef( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
       // check remove reference for a few different types
      BCL_ExampleCheck( GetStaticClassName< type::RemoveReference< const char &>::Type>(), "const-char");
      BCL_ExampleCheck( GetStaticClassName< type::RemoveReference< char>::Type>(), "char");

      // check remove const for a few different types
      BCL_ExampleCheck( GetStaticClassName< type::RemoveConst< const char &>::Type>(), "const-char&");
      BCL_ExampleCheck( GetStaticClassName< type::RemoveConst< const char>::Type>(), "char");
      BCL_ExampleCheck( GetStaticClassName< type::RemoveConst< const char *>::Type>(), "const-char*");

      // check remove const reference for a few different types
      BCL_ExampleCheck( GetStaticClassName< type::RemoveConstRef< const char &>::Type>(), "char");
      BCL_ExampleCheck( GetStaticClassName< type::RemoveConstRef< const char>::Type>(), "char");
      BCL_ExampleCheck( GetStaticClassName< type::RemoveConstRef< const char *const &>::Type>(), "const-char*");
      BCL_ExampleCheck( GetStaticClassName< type::RemoveConstRef< const char *&>::Type>(), "const-char*");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleTypeRemoveConstRef

  const ExampleClass::EnumType ExampleTypeRemoveConstRef::s_Instance
  (
    GetExamples().AddEnum( ExampleTypeRemoveConstRef())
  );

} // namespace bcl
