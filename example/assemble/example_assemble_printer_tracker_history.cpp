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
#include "assemble/bcl_assemble_printer_tracker_history.h"

// includes from bcl - sorted alphabetically
#include "opti/bcl_opti_criterion_number_iterations.h"
#include "util/bcl_util_wrapper_base.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_printer_tracker_history.cpp
  //! @brief this example tests the implementation of the assmble::PrinterTrackerHistory which prints the apllied
  //! mutates as well as their outcome and energy scores
  //!
  //! @author fischea
  //! @date Sep 22, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssemblePrinterTrackerHistory :
     public ExampleInterface
  {

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

  public:

    //! @brief Clone function
    //! @return pointer to a new ExampleAssemblePrinterTrackerHistory
    ExampleAssemblePrinterTrackerHistory *Clone() const
    {
      return new ExampleAssemblePrinterTrackerHistory( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns the class name
    //! @return the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

  ////////////////
  // operations //
  ////////////////

    //! @brief run routine
    //! this is performing the execution of the example
    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      assemble::PrinterTrackerHistory printer;

    /////////////////
    // data access //
    /////////////////

      // check if the correct class name is returned
      BCL_ExampleCheck
      (
        printer.GetClassIdentifier(),
        ( GetStaticClassName< assemble::PrinterTrackerHistory>())
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // class ExampleAssemblePrinterTrackerHistory

  //! single instance of that class
  const ExampleClass::EnumType ExampleAssemblePrinterTrackerHistory::s_Instance
  (
     GetExamples().AddEnum( ExampleAssemblePrinterTrackerHistory())
  );

} // namespace bcl
