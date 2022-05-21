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
#include "command/bcl_command.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_command.cpp
  //! @details The example class for the  command line
  //! note: this also checks the functionality of other objects associated with the command line
  //!       (command flags, parameters, parameter checks), but the functionality of these objects
  //!       is more extensively checked in their own example classes
  //!
  //! @author woetzen, durhamea
  //! @date 11/02/2007, 11/06/2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleCommand :
    public ExampleInterface
  {
  public:

    ExampleCommand *Clone() const
    { return new ExampleCommand( *this);}

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! the function called when the examples are run
    //! @brief creates a correct and incorrect command line and calls Main_program with both of them
    //! @brief tests the read/write functions
    //! @see Main_program
    int Run() const
    {
      // check if get namespace identifier is working
      BCL_MessageStd
      (
        "this is the bcl::command namespace identifier: " + command::GetNamespaceIdentifier()
      );
      BCL_ExampleCheck( command::GetNamespaceIdentifier(), "bcl::command");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

    /////////////////
    // data access //
    /////////////////

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
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleCommand

  const ExampleClass::EnumType ExampleCommand::s_Instance( GetExamples().AddEnum( ExampleCommand()));

} // namespace bcl
