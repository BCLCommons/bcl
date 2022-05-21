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
#include "util/bcl_util_format.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //
  //! @example example_util_format.cpp
  //! @details this examples demonstrates how to use the format object and what it can do
  //!
  //! @author woetzen, alexanns
  //! @date Dec 1, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilFormat :
    public ExampleInterface
  {
  public:

    ExampleUtilFormat *Clone() const
    { return new ExampleUtilFormat( *this);}

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

    ////////////////
    // operations //
    ////////////////

    ///////////////
    // operators //
    ///////////////

      //write a number with with = 7, and precision (significant digits after '.') to 3
      BCL_ExampleCheck( util::Format().W( 7).FFP( 3)( 5678.9243), "5678.924");

      //format a float using scientific representation with width of 12 and precesion 3 (using "e-11" or similar)

      // on mingw and possibly others, there may be a fixed width of the exponent
      // leading to util::Format().SFP( 3)( float( 0.00000000000056)) giving back 5.6e-013 rather than 5.6e-13
      // Both are fine though

      BCL_ExampleIndirectCheck
      (
        util::Format().SFP( 3)( float( 0.00000000000056)) == "5.600e-13"
        || util::Format().SFP( 3)( float( 0.00000000000056)) == "5.600e-013",
        true,
        "formatting, actual result = util::Format().SFP( 3)( float( 0.00000000000056))"
      );

      //format float using fixed number formatting
      BCL_ExampleCheck
      (
        util::Format().W( 7).FFP( 3)( float( 0.0000000056)),
        "  0.000"
      );
      //write string with defined width and left alignment
      BCL_ExampleCheck( util::Format().L().W(8)( "funny"), "funny   ");
      //write bool in string with width 8 and right alignment
      BCL_ExampleCheck( util::Format().R().W(8)( false), "       0");
      //write bool in string with width 8 and left alignment
      BCL_ExampleCheck( util::Format().L().W(8)( true), "1       ");
      //format float with additional SingleSpace
      BCL_ExampleCheck( util::Format().S().W( 7).FFP( 3)( float( 7.3)), "   7.300");
      BCL_ExampleCheck( util::Format().S().W( 7).FFP( 3)( 7.3), "   7.300");
      BCL_ExampleCheck( util::Format().S().W( 9)( 12), "        12");

      BCL_ExampleCheck( util::Format().W( 9).L()( 'C'), "C        ");
      BCL_ExampleCheck( util::Format().W( 9).R()( 'C'), "        C");

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleUtilFormat

  const ExampleClass::EnumType ExampleUtilFormat::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilFormat())
  );

} // namespace bcl
