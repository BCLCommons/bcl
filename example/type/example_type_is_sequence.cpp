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
#include "type/bcl_type_is_sequence.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_type_is_sequence.cpp
  //!
  //! @author mendenjl
  //! @date Nov 29, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleTypeIsSequence :
    public ExampleInterface
  {
  public:

    ExampleTypeIsSequence *Clone() const
    { return new ExampleTypeIsSequence( *this);}

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
      // test is sequence on something with a const_iterator
      BCL_ExampleCheck( type::IsSequence< storage::Vector< float> >::value, 1);

      // test is sequence on something without a const iterator
      BCL_ExampleCheck( type::IsSequence< float>::value, 0);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleTypeIsSequence

  const ExampleClass::EnumType ExampleTypeIsSequence::s_Instance
  (
    GetExamples().AddEnum( ExampleTypeIsSequence())
  );

} // namespace bcl
