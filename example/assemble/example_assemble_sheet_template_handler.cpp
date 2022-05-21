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
#include "assemble/bcl_assemble_sheet_template_handler.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_fold_template.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_sheet_template_handler.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleSheetTemplateHandler :
    public ExampleInterface
  {
  public:

    ExampleAssembleSheetTemplateHandler *Clone() const
    {
      return new ExampleAssembleSheetTemplateHandler( *this);
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
      // initialize a set of excluded pdbs
      storage::Set< std::string> exclude_set_2CRT;
      exclude_set_2CRT.Insert( "2CRT");

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      assemble::SheetTemplateHandler handler;

      // clone constructor
      util::ShPtr< assemble::SheetTemplateHandler> sp_handler( handler.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< assemble::SheetTemplateHandler>(), handler.GetClassIdentifier());

    ////////////////
    // operations //
    ////////////////

      // initialize number strands
      const size_t nr_strands( 3);

      BCL_MessageStd( "testing get random template excluding 2CRT");
      const assemble::FoldTemplate &template_a( handler.GetRandomTemplate( nr_strands));
      BCL_ExampleCheck( template_a.GetStrandGeometries().GetSize(), nr_strands);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleSheetTemplateHandler

  const ExampleClass::EnumType ExampleAssembleSheetTemplateHandler::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleSheetTemplateHandler())
  );

} // namespace bcl
