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
#include "assemble/bcl_assemble_fold_template_handler.h"

// includes from bcl - sorted alphabetically
#include "example_proteins.h"
#include "assemble/bcl_assemble_fold_template.h"
#include "assemble/bcl_assemble_protein_model.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_assemble_fold_template_handler.cpp
  //!
  //! @author weinerbe
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAssembleFoldTemplateHandler :
    public ExampleInterface
  {
  public:

    ExampleAssembleFoldTemplateHandler *Clone() const
    {
      return new ExampleAssembleFoldTemplateHandler( *this);
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
      storage::Map< biol::SSType, size_t> sse_min_sizes;
      sse_min_sizes[ biol::GetSSTypes().HELIX] = 5;
      sse_min_sizes[ biol::GetSSTypes().STRAND] = 3;
      sse_min_sizes[ biol::GetSSTypes().COIL] = 999;
      const assemble::ProteinModel protein_model
      (
        Proteins::GetModel
        (
          AddExampleInputPathToFilename( e_Biology, "1ubi.pdb"), biol::GetAAClasses().e_AABackBone, sse_min_sizes
        )
      );
      const util::SiPtrVector< const assemble::SSE> sses( protein_model.GetSSEs());

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      BCL_MessageStd( "testing default constructor");
      const assemble::FoldTemplateHandler def_construct;

      // test clone function
      BCL_MessageStd( "testing clone function");
      const util::ShPtr< assemble::FoldTemplateHandler> clone_construct( def_construct.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( GetStaticClassName< assemble::FoldTemplateHandler>(), def_construct.GetClassIdentifier());

    ////////////////
    // operations //
    ////////////////

      // test GetRandomTemplate
      BCL_MessageStd( "testing GetRandomTemplate function");
      const assemble::FoldTemplate &random_template( assemble::FoldTemplateHandler::GetRandomTemplate( 1, 4));
      BCL_ExampleIndirectCheck
      (
        random_template.GetHelicalGeometries().GetSize() == 1 && random_template.GetStrandGeometries().GetSize() == 4,
        true,
        "GetRandomTemplate produced a template with " + util::Format()( random_template.GetHelicalGeometries().GetSize()) +
        " helices and " + util::Format()( random_template.GetStrandGeometries().GetSize()) +
        " strands instead of 1 and 4, respectively"
      );

      // test GetRandomTemplate
      BCL_MessageStd( "testing GetRandomTemplate function with SSEs");
      const assemble::FoldTemplate &random_template_sses( assemble::FoldTemplateHandler::GetRandomTemplate( sses));
      BCL_ExampleIndirectCheck
      (
        random_template_sses.GetHelicalGeometries().GetSize() == 1 &&
        random_template_sses.GetStrandGeometries().GetSize() == 5,
        true,
        "GetRandomTemplate using SSEs produced a template with " +
        util::Format()( random_template_sses.GetHelicalGeometries().GetSize()) +
        " helices and " + util::Format()( random_template_sses.GetStrandGeometries().GetSize()) +
        " strands instead of 1 and 4, respectively"
      );

      // test GetRandomSubTemplate
      BCL_MessageStd( "testing GetRandomSubTemplate function");
      const assemble::FoldTemplate random_sub_template( assemble::FoldTemplateHandler::GetRandomSubTemplate( sses));
      BCL_ExampleIndirectCheck
      (
        random_sub_template.GetHelicalGeometries().GetSize() == 1 &&
        random_sub_template.GetStrandGeometries().GetSize() == 5,
        true,
        "GetRandomSubTemplate produced a template with " +
        util::Format()( random_sub_template.GetHelicalGeometries().GetSize()) +
        " helices and " + util::Format()( random_sub_template.GetStrandGeometries().GetSize()) +
        " strands instead of 1 and 4, respectively"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleAssembleFoldTemplateHandler

  const ExampleClass::EnumType ExampleAssembleFoldTemplateHandler::s_Instance
  (
    GetExamples().AddEnum( ExampleAssembleFoldTemplateHandler())
  );

} // namespace bcl
