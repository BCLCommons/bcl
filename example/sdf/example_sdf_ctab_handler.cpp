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
#include "sdf/bcl_sdf_ctab_handler.h"

// includes from bcl - sorted alphabetically
#include "io/bcl_io_file.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_sdf_ctab_handler.cpp
  //!
  //! @author geanesar
  //! @date
  //! @remarks status incomplete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleSdfCTabHandler :
    public ExampleInterface
  {
  public:

    ExampleSdfCTabHandler *Clone() const
    {
      return new ExampleSdfCTabHandler( *this);
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
      return 0;
      /*

      io::IFStream read;
      BCL_ExampleMustOpenInputFile( read, AddExampleInputPathToFilename( e_Chemistry, "hexane.sdf"));

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // check default constructor
      sdf::CTabHandler handler_default;

      // check constructor with parameters
      sdf::CTabHandler handler( read);

      // check clone
      util::ShPtr< sdf::CTabHandler> sp_handler( handler.Clone());

      // close input file stream
      io::File::CloseClearFStream( read);

    /////////////////
    // data access //
    /////////////////

      // test GetMdlDescriptionLine
      BCL_ExampleCheck( handler.GetDescription(), "Untitled Document-1\n  ChemDraw06180812522D\n");

      // test GetMdlAtomLines
      BCL_ExampleCheck( handler.GetAtomInfo().GetSize(), size_t( 6));
      // test GetMdlBondLines
      BCL_ExampleCheck( handler.GetBondInfo().GetSize(), size_t( 5));

      // test IsValid
      BCL_ExampleCheck( handler.IsValid(), true);
      BCL_ExampleCheck( handler_default.IsValid(), false);

      // test GetMiscPropertiesLines
      BCL_ExampleCheck( handler_default.GetMiscProperties().IsEmpty(), true);
      BCL_ExampleCheck( handler.GetMiscProperties().GetSize(), size_t( 4));

      // check is terminal line with terminal line from sdf file written on windows, then linux
      BCL_ExampleCheck( sdf::CTabHandler::IsTerminalLine( "$$$$\r"), true);
      BCL_ExampleCheck( sdf::CTabHandler::IsTerminalLine( "$$$$"),   true);

      // try is terminal line with a line that is not terminal
      BCL_ExampleCheck( sdf::CTabHandler::IsTerminalLine( "$$$"),    false);

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
      */
    } // Run

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleSdfCTabHandler

  const ExampleClass::EnumType ExampleSdfCTabHandler::s_Instance
  (
    GetExamples().AddEnum( ExampleSdfCTabHandler())
  );

} // namespace bcl
