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
#include "contact/bcl_contact_types.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_sse.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_contact_types.cpp
  //!
  //! @author karakam
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleContactTypes :
    public ExampleInterface
  {
  public:

    ExampleContactTypes *Clone() const
    {
      return new ExampleContactTypes( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! returns class name
    //! the class name as const ref std::string
    const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct AAClass from each class type
      BCL_MessageStd
      (
        "Constructing contact::Type"
      );

      contact::Type type_a( contact::GetTypes().HELIX_HELIX);

      // construct undefined AAClass
      BCL_ExampleCheck( util::GetUndefined< contact::Type>().IsDefined(), false);

      // use copy constructor
      BCL_ExampleAssert( contact::Type( type_a), type_a);

    /////////////////
    // data access //
    /////////////////

      // display total number of contact::Types
      BCL_MessageStd
      (
        "The total number of contact::Types is " + util::Format()( contact::GetTypes().GetEnumCount())
      );

      // test total number of contact::Types
      BCL_ExampleCheck( contact::GetTypes().GetEnumCount(), 10);

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

      // construct multiple AATypes
      contact::Type type_helix_helix( contact::GetTypes().HELIX_HELIX);
      contact::Type type_helix_sheet( contact::GetTypes().HELIX_SHEET);
      contact::Type type_sheet_helix( contact::GetTypes().SHEET_HELIX);
      contact::Type type_helix_strand( contact::GetTypes().HELIX_STRAND);
      contact::Type type_strand_helix( contact::GetTypes().STRAND_HELIX);
      contact::Type type_undef_helix_strand( contact::GetTypes().UNDEFINED_HELIX_STRAND);
      contact::Type type_undef_strand_helix( contact::GetTypes().UNDEFINED_STRAND_HELIX);

      // test reverse function
      BCL_ExampleCheck( contact::GetTypes().Reverse( contact::GetTypes().HELIX_HELIX), contact::GetTypes().HELIX_HELIX);
      BCL_ExampleCheck( contact::GetTypes().Reverse( contact::GetTypes().HELIX_SHEET), contact::GetTypes().SHEET_HELIX);
      BCL_ExampleCheck( contact::GetTypes().Reverse( contact::GetTypes().SHEET_HELIX), contact::GetTypes().HELIX_SHEET);

      BCL_ExampleCheck
      (
        contact::GetTypes().Reverse( contact::GetTypes().STRAND_HELIX),
        contact::GetTypes().HELIX_STRAND
      );
      BCL_ExampleCheck
      (
        contact::GetTypes().Reverse( contact::GetTypes().UNDEFINED_STRAND_HELIX),
        contact::GetTypes().UNDEFINED_HELIX_STRAND
      );

      // construct one helix one strand SSE
      const assemble::SSE sse_helix( biol::GetSSTypes().HELIX);
      const assemble::SSE sse_strand( biol::GetSSTypes().STRAND);
      const util::SiPtr< const assemble::SSE> sse_si_ptr_helix( sse_helix);
      const util::SiPtr< const assemble::SSE> sse_si_ptr_strand( sse_strand);

      // construct SSType pair
      storage::Pair< biol::SSType, biol::SSType> helix_strand_ss_type_pair
      (
        biol::GetSSTypes().HELIX, biol::GetSSTypes().STRAND
      );

      // construct SSE SiPtr VectorND
      const storage::VectorND< 2, util::SiPtr< const assemble::SSE> > helix_strand_vector_nd
      (
        sse_si_ptr_helix, sse_si_ptr_strand
      );

      // test TypeFromSSTypes that takes two sses
      BCL_ExampleCheck( contact::GetTypes().TypeFromSSTypes( sse_helix, sse_strand), contact::GetTypes().HELIX_SHEET);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleContactTypes

  const ExampleClass::EnumType ExampleContactTypes::s_Instance
  (
    GetExamples().AddEnum( ExampleContactTypes())
  );

} // namespace bcl

