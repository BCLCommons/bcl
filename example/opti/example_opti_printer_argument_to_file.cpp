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
#include "opti/bcl_opti_printer_argument_to_file.h"

// includes from bcl - sorted alphabetically
#include "util/bcl_util_wrapper_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_opti_printer_argument_to_file.cpp
  //!
  //! @author butkiem1
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleOptiPrinterArgumentToFile :
    public ExampleInterface
  {
  public:

    ExampleOptiPrinterArgumentToFile *Clone() const
    {
      return new ExampleOptiPrinterArgumentToFile( *this);
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

      // check default printer
      opti::PrinterArgumentToFile< int, int> printer_default;

      // check printer with parameters
      opti::PrinterArgumentToFile< int, int> printer
      (
        AddExampleOutputPathToFilename( opti::GetNamespaceIdentifier(), "opti_printer_to_file.txt"),
        size_t( 1)
      );

    /////////////////
    // data access //
    /////////////////

      BCL_ExampleCheck
      (
        printer.GetClassIdentifier(),
        ( GetStaticClassName< opti::PrinterArgumentToFile< int, int> >())
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      BCL_ExampleCheck
      (
        printer.GetFilenamePrefix() + io::File::GetExtensionDelimiter() + "1" + printer.GetFilenameExtension(),
        AddExampleOutputPathToFilename( opti::GetNamespaceIdentifier(), "opti_printer_to_file.1.txt")
      );

    //////////////////////
    // input and output //
    //////////////////////

      // open output file
      WriteBCLObject( printer);

      // read bcl object
      opti::PrinterArgumentToFile< int, int> printer_read;
      ReadBCLObject( printer_read);

      // check read in printer
      BCL_ExampleCheck
      (
        printer_read.GetFilenamePrefix() + io::File::GetExtensionDelimiter() + "1" + printer_read.GetFilenameExtension(),
        AddExampleOutputPathToFilename( opti::GetNamespaceIdentifier(), "opti_printer_to_file.1.txt")
      );

      // test model result pair
      opti::Tracker< int, int> tracker;
      tracker.Track( util::CloneToShPtr( storage::Pair< int, int>( 2, 6)));

      // write out model into file
      printer.Print( tracker);

      // empty argument
      size_t argument_read;

      // open written out input file
      io::IFStream input_stream;
      BCL_ExampleMustOpenInputFile
      (
        input_stream,
        AddExampleOutputPathToFilename( opti::GetNamespaceIdentifier(), "opti_printer_to_file.1.txt")
      );

      // read model to file stream
      input_stream >> argument_read;
      io::File::CloseClearFStream( input_stream);

      // check read in argument
      BCL_ExampleCheck( argument_read, size_t( 2));

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleOptiPrinterArgumentToFile

  const ExampleClass::EnumType ExampleOptiPrinterArgumentToFile::s_Instance
  (
    GetExamples().AddEnum( ExampleOptiPrinterArgumentToFile())
  );

} // namespace bcl
