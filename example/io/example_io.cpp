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
#include "io/bcl_io.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_io.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleIo :
    public ExampleInterface
  {
  public:

    ExampleIo *Clone() const
    {
      return new ExampleIo( *this);
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

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

      // check if get namespace identifier is working
      BCL_MessageStd
      (
        "this is the bcl::io namespace identifier: " + io::GetNamespaceIdentifier()
      );
      BCL_Example_Check
      (
        io::GetNamespaceIdentifier() == "bcl::io",
        "namespace identifier for io is incorrect!"
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleIo

  const ExampleClass::EnumType ExampleIo::s_Instance
  (
    GetExamples().AddEnum( ExampleIo())
  );

} // namespace bcl
