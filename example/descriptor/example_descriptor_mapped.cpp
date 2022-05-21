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
#include "descriptor/bcl_descriptor_mapped.h"

// includes from bcl - sorted alphabetically
#include "descriptor/bcl_descriptor_example_string_sequence.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_descriptor_mapped.cpp
  //!
  //! @author mendenjl
  //! @date Mar 14, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleDescriptorMapped :
    public ExampleInterface
  {
  public:

    ExampleDescriptorMapped *Clone() const
    {
      return new ExampleDescriptorMapped( *this);
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

    int Run() const
    {
      // construct an implementation that should map values of the property Mol_ID onto TestValues
      const util::Implementation< descriptor::Base< char, float> >
      mapper_float
      (
        "Mapped("
        "file=" + AddExampleInputPathToFilename( e_Descriptor, "input_map.object") + ","
        "key=Character,"
        "delimiter=;"
        ")"
      );

      // ensure that each is defined
      BCL_ExampleAssert( mapper_float.IsDefined(), true);

      // check the descriptor string
      BCL_ExampleCheck( mapper_float.GetLabel().GetValue(), "Mapped");

      // get the value that each id should map to in each molecule
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( mapper_float, "acnb", 1),
        "3.0 14.0 18.0 ; 5.0 6.0 12.0 ; nan nan nan ; 4.6 5.9 10.0 ; "
      );

      // now try a character mapper on the same file
      const util::Implementation< descriptor::Base< char, char> >
      mapper_char
      (
        "Mapped("
        "file=" + AddExampleInputPathToFilename( e_Descriptor, "input_map.object") + ","
        "key=Character,"
        "delimiter=;"
        ")"
      );
      // ensure that each is defined
      BCL_ExampleAssert( mapper_char.IsDefined(), true);

      // get the value that each id should map to in each molecule
      BCL_ExampleCheck
      (
        descriptor::StringSequence::WriteIterations( mapper_char, "acnb"),
        " 3.0, 14.0,   18.0;     5.0 6.0 12.0  ;                   ;  4.6 , 5.9,   10.0; "
      );

      return 0;
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleDescriptorMapped

  const ExampleClass::EnumType ExampleDescriptorMapped::s_Instance
  (
    GetExamples().AddEnum( ExampleDescriptorMapped())
  );

} // namespace bcl
