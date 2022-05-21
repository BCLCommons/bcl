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
#include "util/bcl_util_object_data_label_tokenizer.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_util_object_data_label_tokenizer.cpp
  //!
  //! @author mendenjl
  //! @date Sep 28, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilObjectDataLabelTokenizer :
    public ExampleInterface
  {
  public:

    ExampleUtilObjectDataLabelTokenizer *Clone() const
    {
      return new ExampleUtilObjectDataLabelTokenizer( *this);
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

      // try constructors from empty strings, ensure that the types are correct
      BCL_ExampleCheck( util::ObjectDataLabelTokenizer( "").GetNextTokenType(), util::ObjectDataLabelTokenizer::e_End);
      BCL_ExampleCheck
      (
        util::ObjectDataLabelTokenizer( "  ").GetLastTokenType(),
        util::ObjectDataLabelTokenizer::e_Start
      );
      BCL_ExampleCheck
      (
        util::ObjectDataLabelTokenizer( "  ").GetNextTokenType(),
        util::ObjectDataLabelTokenizer::e_End
      );
      BCL_ExampleCheck( util::ObjectDataLabelTokenizer( "  ").Pop(), std::string());
      BCL_ExampleCheck( util::ObjectDataLabelTokenizer( "  ").GetScopeDepth(), 0);

      // try constructor from a single value
      util::ObjectDataLabelTokenizer value( " 0.1");
      BCL_ExampleCheck( util::ObjectDataLabelTokenizer( " 0.1").Pop(), "0.1");
      BCL_ExampleCheck
      (
        util::ObjectDataLabelTokenizer( " 0.1").GetNextTokenType(),
        util::ObjectDataLabelTokenizer::e_Scalar
      );
      BCL_ExampleCheck( util::ObjectDataLabelTokenizer( " 0.1").Pop(), "0.1");
      value.Pop();
      BCL_ExampleIndirectCheck
      (
        value.GetLastTokenType(),
        util::ObjectDataLabelTokenizer::e_Scalar,
        "last type is updated"
      );
      BCL_ExampleIndirectCheck
      (
        value.GetNextTokenType(),
        util::ObjectDataLabelTokenizer::e_End,
        "next type is updated"
      );

      util::ObjectDataLabelTokenizer call_square( "square(            0.1     )");
      BCL_ExampleCheck( call_square.Pop(), "square");
      BCL_ExampleCheck( call_square.Pop(), "(");
      BCL_ExampleIndirectCheck
      (
        call_square.GetNextTokenType(),
        util::ObjectDataLabelTokenizer::e_Scalar,
        "next type inside scope"
      );
      BCL_ExampleIndirectCheck( call_square.GetScopeDepth(), 1, "scope depth inside scope");
      BCL_ExampleIndirectCheck
      (
        call_square.GetLastTokenType(),
        util::ObjectDataLabelTokenizer::e_ScopeOpen,
        "last type inside scope"
      );
      BCL_ExampleCheck( call_square.Pop(), "0.1");
      BCL_ExampleCheck( call_square.Pop(), ")");
      BCL_ExampleIndirectCheck( call_square.GetScopeDepth(), 0, "scope depth after scope");
      BCL_ExampleIndirectCheck( call_square.Pop(), "", "Calling pop after all arguments have been popped");

      // test a harder example that includes aruments
      util::ObjectDataLabelTokenizer nested_equal( "nest(nest(equal( apples, oranges) , tag =  value(x)))");
      BCL_ExampleCheck( nested_equal.Pop(), "nest");
      BCL_ExampleCheck( nested_equal.Pop(), "(");
      BCL_ExampleCheck( nested_equal.Pop(), "nest");
      BCL_ExampleCheck( nested_equal.Pop(), "(");
      BCL_ExampleCheck( nested_equal.GetScopeDepth(), 2);
      BCL_ExampleCheck( nested_equal.Pop(), "equal");
      BCL_ExampleCheck( nested_equal.Pop(), "(");
      BCL_ExampleCheck( nested_equal.GetScopeDepth(), 3);
      BCL_ExampleCheck( nested_equal.Pop(), "apples");
      BCL_ExampleCheck( nested_equal.Pop(), ",");
      BCL_ExampleCheck( nested_equal.GetLastTokenType(), util::ObjectDataLabelTokenizer::e_ArgDelimiter);
      BCL_ExampleCheck( nested_equal.Pop(), "oranges");
      BCL_ExampleCheck( nested_equal.Pop(), ")");
      BCL_ExampleCheck( nested_equal.Pop(), ",");
      BCL_ExampleCheck( nested_equal.Pop(), "tag");
      BCL_ExampleCheck( nested_equal.Pop(), "=");
      BCL_ExampleCheck( nested_equal.GetScopeDepth(), 2);
      BCL_ExampleCheck( nested_equal.GetLastTokenType(), util::ObjectDataLabelTokenizer::e_TagDelimiter);

      // go until the end
      while( nested_equal.GetNextTokenType() != util::ObjectDataLabelTokenizer::e_End)
      {
        nested_equal.Pop();
      }

      // ensure that all scopes are closed
      BCL_ExampleCheck( nested_equal.GetScopeDepth(), 0);

      // ensure that quoted strings are maintained as scalars
      util::ObjectDataLabelTokenizer name_with_delimiter( "\"C=O\"=\"BondType=Aromatic\"(Single)");

      BCL_ExampleIndirectCheck( name_with_delimiter.Pop(), "C=O", "Pop on string beginning with \"C=O\"=...");
      BCL_ExampleIndirectCheck( name_with_delimiter.Pop(), "=", "Pop on string beginning with \"C=O\"=...");
      BCL_ExampleIndirectCheck( name_with_delimiter.Pop(), "BondType=Aromatic", "Pop on string continuing with \"BondType=Aromatic\"");
      BCL_ExampleIndirectCheck( name_with_delimiter.Pop(), "(", "Pop after quoted scalar");

      // try out the validation
      BCL_ExampleCheck( util::ObjectDataLabelTokenizer::Validate( "", util::GetLogger()), true);
      BCL_ExampleCheck( util::ObjectDataLabelTokenizer::Validate( "a=", util::GetLogger()), true);
      BCL_ExampleCheck( util::ObjectDataLabelTokenizer::Validate( "a=(x)", util::GetLogger()), true);
      BCL_ExampleCheck( util::ObjectDataLabelTokenizer::Validate( "\"a\"  = x ( x + y = x)", util::GetLogger()), true);
      BCL_ExampleCheck( util::ObjectDataLabelTokenizer::Validate( "\"a\"  = \"a\" = error", util::GetLogger()), false);
      BCL_ExampleCheck( util::ObjectDataLabelTokenizer::Validate( "(", util::GetLogger()), false);
      BCL_ExampleCheck( util::ObjectDataLabelTokenizer::Validate( ")", util::GetLogger()), false);
      BCL_ExampleCheck( util::ObjectDataLabelTokenizer::Validate( "()x", util::GetLogger()), false);
      BCL_ExampleCheck( util::ObjectDataLabelTokenizer::Validate( "a,b", util::GetLogger()), false);
      BCL_ExampleCheck( util::ObjectDataLabelTokenizer::Validate( "(a(),b)", util::GetLogger()), true);

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

  }; //end ExampleUtilObjectDataLabelTokenizer

  const ExampleClass::EnumType ExampleUtilObjectDataLabelTokenizer::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilObjectDataLabelTokenizer())
  );

} // namespace bcl
