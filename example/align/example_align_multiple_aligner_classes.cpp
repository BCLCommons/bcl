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
#include "align/bcl_align_multiple_aligner_classes.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_align_multiple_aligner_classes.cpp
  //!
  //! @author heinzes1
  //! @date Jan 25, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleAlignMultipleAlignerClasses :
    public ExampleInterface
  {
  public:

      ExampleAlignMultipleAlignerClasses *Clone() const
    {
      return new ExampleAlignMultipleAlignerClasses( *this);
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
      align::GetMultipleAlignerClasses< biol::AABase>();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // construct Aligner from each class type
      BCL_MessageStd( "1: Constructing MultipleAlignerClass object from e_AlignerMerge");
      align::MultipleAlignerClasses< biol::AABase>::MultipleAlignerClass
        alignerclass_merge( align::GetMultipleAlignerClasses< biol::AABase>().e_AlignerMerge);
      BCL_ExampleIndirectCheck
      (
        alignerclass_merge,
        align::GetMultipleAlignerClasses< biol::AABase>().e_AlignerMerge,
        "1: Constructing MultipleAligner did not work for "
          + util::Format()( align::GetMultipleAlignerClasses< biol::AABase>().e_AlignerMerge->GetClassIdentifier())
      );

      BCL_MessageStd( "2: Constructing MultipleAlignerClass object from e_AlignerProgressive");
      align::MultipleAlignerClasses< biol::AABase>::MultipleAlignerClass alignerclass_progressive
      (
        align::GetMultipleAlignerClasses< biol::AABase>().e_AlignerProgressiveDynamicProgramming
      );
      BCL_ExampleIndirectCheck
      (
        alignerclass_progressive,
        align::GetMultipleAlignerClasses< biol::AABase>().e_AlignerProgressiveDynamicProgramming,
        "2: Constructing MultipleAligner did not work for "
          + util::Format()( align::GetMultipleAlignerClasses< biol::AABase>().e_AlignerProgressiveDynamicProgramming->GetClassIdentifier())
      );

      // construct undefined MultipleAlignerClass
      BCL_MessageStd( "3: Constructing an undefined MultipleAlignerClass");
      align::MultipleAlignerClasses< biol::AABase>::MultipleAlignerClass
        alignerclass_undefined( util::GetUndefined< align::MultipleAlignerClasses< biol::AABase>::MultipleAlignerClass>());
      BCL_ExampleIndirectCheck
      (
        !alignerclass_undefined.IsDefined(),
        true,
        "3: Undefined MultipleAligner was not constructed correctly"
      );

      // use copy constructor
      BCL_MessageStd( "4: Copy constructor");
      align::MultipleAlignerClasses< biol::AABase>::MultipleAlignerClass alignerclass_merge_copy( alignerclass_merge);
      BCL_ExampleIndirectCheck( alignerclass_merge_copy, alignerclass_merge, "4: Copy constructed object incorrect");

    /////////////////
    // data access //
    /////////////////

      // test GetClassIdentifier()
      BCL_MessageStd( "5: GetClassIdentifier()");
      const std::string class_identifier( "bcl::align::MultipleAlignerClasses<bcl::biol::AABase>");
      std::string returned_class_identifier( align::GetMultipleAlignerClasses< biol::AABase>().GetClassIdentifier());
      BCL_ExampleIndirectCheck
      (
        returned_class_identifier,
        class_identifier,
        "5: returned_class_identifier==" + util::Format()( returned_class_identifier) + "!="
          + util::Format()( class_identifier) + "==class_identifier"
      );

      // test GetEnumCount()
      BCL_MessageStd( "6: GetEnumCount()");
      const size_t enum_count( 3);
      size_t returned_enum_count( align::GetMultipleAlignerClasses< biol::AABase>().GetEnumCount());
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

  }; //end ExampleAlignMultipleAlignerClasses

  const ExampleClass::EnumType ExampleAlignMultipleAlignerClasses::s_Instance
  (
    GetExamples().AddEnum( ExampleAlignMultipleAlignerClasses())
  );
  
} // namespace bcl
