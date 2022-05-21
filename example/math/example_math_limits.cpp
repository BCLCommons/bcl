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
#include "math/bcl_math_limits.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_limits.cpp
  //!
  //! @author mendenjl
  //! @date Oct 18, 2012
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathLimits :
    public ExampleInterface
  {
  public:

    ExampleMathLimits *Clone() const
    { return new ExampleMathLimits( *this);}

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
    //////////////
    // Limits //
    //////////////

      BCL_ExampleCheck( math::GetHighestBoundedValue< double>(), std::numeric_limits< double>::max());
      BCL_ExampleCheck( math::GetLowestBoundedValue< double>(), -std::numeric_limits< double>::max());
      BCL_ExampleCheck( math::GetHighestBoundedValue< int>(), std::numeric_limits< int>::max());
      BCL_ExampleCheck( math::GetLowestBoundedValue< int>() , std::numeric_limits< int>::min());
      BCL_ExampleCheck( math::GetLowestBoundedValue< int>() , math::GetLowestUnboundedValue< int>());
      BCL_ExampleCheck( math::GetHighestBoundedValue< size_t>(), std::numeric_limits< size_t>::max());
      BCL_ExampleCheck( math::GetLowestBoundedValue< size_t>(), 0);

      BCL_ExampleCheck( math::GetHighestUnboundedValue< double>() > math::GetHighestBoundedValue< double>(), true);
      BCL_ExampleCheck( math::GetHighestUnboundedValue< float>() > math::GetHighestBoundedValue< float>(), true);
      BCL_ExampleCheck( math::GetLowestBoundedValue< double>() < std::numeric_limits< double>::min(), true);
      BCL_ExampleCheck( math::GetLowestUnboundedValue< double>() < math::GetLowestBoundedValue< double>(), true);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathLimits

  const ExampleClass::EnumType ExampleMathLimits::s_Instance
  (
    GetExamples().AddEnum( ExampleMathLimits())
  );

} // namespace bcl
