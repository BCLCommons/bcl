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
#include "align/bcl_align_handler_classes.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_align_handler_classes.cpp
  //!
  //! @author heinzes1
  //! @date Jan 25, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignHandlerClasses :
    public ExampleInterface
  {
  public:

      ExampleAlignHandlerClasses *Clone() const
    {
      return new ExampleAlignHandlerClasses( *this);
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
      align::GetHandlerClasses< biol::AABase>();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct Handler from each class type
      BCL_MessageStd( "1: Constructing HandlerClass object from e_BLC");
      align::HandlerClasses< biol::AABase>::HandlerClass handlerclass_blc( align::GetHandlerClasses< biol::AABase>().e_BLC);
      BCL_ExampleIndirectCheck
      (
        handlerclass_blc,
        align::GetHandlerClasses< biol::AABase>().e_BLC,
        "1: Constructing Handler did not work for " + util::Format()( align::GetHandlerClasses< biol::AABase>().e_BLC->GetClassIdentifier())
      );

      BCL_MessageStd( "2: Constructing HandlerClass object from e_PIR");
      align::HandlerClasses< biol::AABase>::HandlerClass handlerclass_pir( align::GetHandlerClasses< biol::AABase>().e_PIR);
      BCL_ExampleIndirectCheck
      (
        handlerclass_pir,
        align::GetHandlerClasses< biol::AABase>().e_PIR,
        "2: Constructing Handler did not work for " + util::Format()( align::GetHandlerClasses< biol::AABase>().e_PIR->GetClassIdentifier())
      );

      // construct undefined HandlerClass
      BCL_MessageStd( "3: Constructing an undefined HandlerClass");
      align::HandlerClasses< biol::AABase>::HandlerClass handlerclass_undefined( util::GetUndefined< align::HandlerClasses< biol::AABase>::HandlerClass>());
      BCL_ExampleIndirectCheck
      (
        !handlerclass_undefined.IsDefined(),
        true,
        "3: Undefined Handler was not constructed correctly"
      );

      // use copy constructor
      BCL_MessageStd( "4: Copy constructor");
      align::HandlerClasses< biol::AABase>::HandlerClass handlerclass_pir_copy( handlerclass_pir);
      BCL_ExampleIndirectCheck( handlerclass_pir_copy, handlerclass_pir, "4: Copy constructed object incorrect");

    /////////////////
    // data access //
    /////////////////

      // test GetClassIdentifier()
      BCL_MessageStd( "5: GetClassIdentifier()");
      const std::string class_identifier( "bcl::align::HandlerClasses<bcl::biol::AABase>");
      std::string returned_class_identifier( align::GetHandlerClasses< biol::AABase>().GetClassIdentifier());
      BCL_ExampleIndirectCheck
      (
        returned_class_identifier,
        class_identifier,
        "5: returned_class_identifier==" + util::Format()( returned_class_identifier) + "!=" + util::Format()( class_identifier) + "==class_identifier"
      );

      // test GetEnumCount()
      BCL_MessageStd( "6: GetEnumCount()");
      const size_t enum_count( 4);
      size_t returned_enum_count( align::GetHandlerClasses< biol::AABase>().GetEnumCount());
      BCL_ExampleIndirectCheck
      (
        returned_enum_count,
        enum_count,
        "6: returned_enum_count==" + util::Format()( returned_enum_count) + "!=" + util::Format()( enum_count) + "==enum_count"
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

  }; //end ExampleAlignHandlerClasses

  const ExampleClass::EnumType ExampleAlignHandlerClasses::s_Instance
  (
    GetExamples().AddEnum( ExampleAlignHandlerClasses())
  );
  
} // namespace bcl

