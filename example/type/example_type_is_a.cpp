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
#include "type/bcl_type_is_a.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_type_is_a.cpp
  //!
  //! @author mendenjl
  //! @date Nov 29, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleTypeIsA :
    public ExampleInterface
  {
  public:

    ExampleTypeIsA *Clone() const
    { return new ExampleTypeIsA( *this);}

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
      // test that doubles can be statically cast to doubles
      BCL_ExampleCheck( ( type::IsA< double, double>::value), 1);

      // test that doubles can be statically cast to floats
      BCL_ExampleCheck( ( type::IsA< double, float>::value), 0);

      // test that object interface cannot be statically cast to
      BCL_ExampleCheck( ( type::IsA< util::ObjectInterface, storage::Vector< double> >::value), 0);
      BCL_ExampleCheck( ( type::IsA< storage::Vector< double>, util::ObjectInterface>::value), 1);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleTypeIsA

  const ExampleClass::EnumType ExampleTypeIsA::s_Instance
  (
    GetExamples().AddEnum( ExampleTypeIsA())
  );

} // namespace bcl
