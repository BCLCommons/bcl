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
#include "type/bcl_type_is_derived_from.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_type_is_derived_from.cpp
  //!
  //! @author mendenjl
  //! @date Nov 29, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleTypeIsDerivedFrom :
    public ExampleInterface
  {
  public:

    ExampleTypeIsDerivedFrom *Clone() const
    { return new ExampleTypeIsDerivedFrom( *this);}

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
      // check that int's and longs are really different types
      BCL_ExampleCheck( ( type::IsDerivedFrom< storage::Vector< double>, util::ObjectInterface>::value), 1);
      BCL_ExampleCheck( ( type::IsDerivedFrom< util::ObjectInterface, storage::Vector< double> >::value), 0);
      BCL_ExampleCheck( ( type::IsDerivedFrom< storage::Vector< double>, ExampleTypeIsDerivedFrom>::value), 0);
      BCL_ExampleCheck( ( type::IsDerivedFrom< const double &, double &>::value), 0);
      BCL_ExampleCheck( ( type::IsDerivedFrom< double, double>::value), 0);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleTypeIsDerivedFrom

  const ExampleClass::EnumType ExampleTypeIsDerivedFrom::s_Instance
  (
    GetExamples().AddEnum( ExampleTypeIsDerivedFrom())
  );

} // namespace bcl
