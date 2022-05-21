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
#include "math/bcl_math_combination.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_combination.cpp
  //!
  //! @author weinerbe
  //! @date Sep 17, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathCombination :
    public ExampleInterface
  {
  public:

    ExampleMathCombination *Clone() const
    {
      return new ExampleMathCombination( *this);
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

      // test default constructor
      math::Combination< size_t> def_construct;
      BCL_ExampleIndirectCheck( def_construct.GetCombinationSize(), util::GetUndefinedSize_t(), "default constructor");

      // test constructor from a set of size_t's
      storage::Set< size_t> test_set;
      for( size_t i( 0); i != 10; ++i)
      {
        test_set.Insert( i);
      }
      const size_t combination_size( 4);
      math::Combination< size_t> set_construct( test_set, combination_size);

      // test clone function
      util::ShPtr< math::Combination< size_t> > clone_construct( set_construct.Clone());
      BCL_ExampleIndirectCheck( set_construct.GetData().GetSize(), clone_construct->GetData().GetSize(), "Clone");

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      BCL_ExampleCheck( clone_construct->GetClassIdentifier(), GetStaticClassName< math::Combination< size_t> >());

      // check GetData
      BCL_ExampleCheck( set_construct.GetData().GetSize(), test_set.GetSize());

      // check GetCombinationSize
      BCL_ExampleCheck( set_construct.GetCombinationSize(), combination_size);

    ////////////////
    // operations //
    ////////////////

      // test GetNumberOfCombinations
      const size_t nr_combinations( 210);
      BCL_ExampleCheck( set_construct.GetNumberOfCombinations(), nr_combinations);

      // test GetAllCombinations
      BCL_ExampleIndirectCheck
      (
        set_construct.GetAllCombinations().GetSize(),
        nr_combinations,
        "GetIndicesOfAllCombinations"
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write the object
      WriteBCLObject( set_construct);

      // read the object back in
      math::Combination< size_t> combination_read;
      ReadBCLObject( combination_read);

      BCL_ExampleIndirectCheck
      (
        combination_read.GetData().GetSize() == test_set.GetSize() &&
        combination_read.GetCombinationSize() == combination_size,
        true,
        "read and write"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathCombination

  const ExampleClass::EnumType ExampleMathCombination::s_Instance
  (
    GetExamples().AddEnum( ExampleMathCombination())
  );
  
} // namespace bcl
