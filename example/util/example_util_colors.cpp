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
#include "util/bcl_util_colors.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_util_colors.cpp
  //!
  //! @author karakam
  //! @date Nov 7, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleUtilColors :
    public ExampleInterface
  {
  public:

    ExampleUtilColors *Clone() const
    {
      return new ExampleUtilColors( *this);
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

      // construct colors from each class type
      BCL_MessageStd( "Constructing Color objects from all known Colors");

      // construct all colors
      util::Color color_black(     util::GetColors().e_Black);
      util::Color color_blue(      util::GetColors().e_Blue);
      util::Color color_cyan(      util::GetColors().e_Cyan);
      util::Color color_green(     util::GetColors().e_Green);
      util::Color color_magenta(   util::GetColors().e_Magenta);
      util::Color color_red(       util::GetColors().e_Red);
      util::Color color_white(     util::GetColors().e_White);
      util::Color color_yellow(    util::GetColors().e_Yellow);

      // check that the constructors were correct
      BCL_ExampleCheck( color_black,     util::GetColors().e_Black);
      BCL_ExampleCheck( color_blue,      util::GetColors().e_Blue);
      BCL_ExampleCheck( color_cyan,      util::GetColors().e_Cyan);
      BCL_ExampleCheck( color_green,     util::GetColors().e_Green);
      BCL_ExampleCheck( color_magenta,   util::GetColors().e_Magenta);
      BCL_ExampleCheck( color_red,       util::GetColors().e_Red);
      BCL_ExampleCheck( color_white,     util::GetColors().e_White);
      BCL_ExampleCheck( color_yellow,    util::GetColors().e_Yellow);

      // construct undefined Color
      BCL_MessageStd( "Constructing an undefined color");
      util::Color color_undefined( util::GetUndefined< util::Color>());
      BCL_ExampleCheck( color_undefined.IsDefined(), false);

      // use copy constructor
      BCL_MessageStd( "Calling copy constructor");
      util::Color copy_green( color_green);
      BCL_ExampleCheck( copy_green, color_green);

    /////////////////
    // data access //
    /////////////////

      // example check
      BCL_ExampleCheck( util::GetColors().GetClassIdentifier(), "bcl::util::Colors");

      // display total number of Measures
      BCL_MessageStd
      (
        "The total number of colors is " + util::Format()( util::GetColors().GetEnumCount())
      );

      // test total number of colors
      BCL_ExampleCheck( util::GetColors().GetEnumCount(), 8);

    ////////////////
    // operations //
    ////////////////

      // initialize rainbow colors vector
      const storage::Vector< util::Color> expected_rainbow
      (
        storage::Vector< util::Color>::Create
        (
          util::GetColors().e_Blue,
          util::GetColors().e_Green,
          util::GetColors().e_Yellow,
          util::GetColors().e_Red
        )
      );

      // test GetRainbo() function
      BCL_ExampleCheck( util::GetColors().GetRainbow(), expected_rainbow);

      // write enums to file
      WriteBCLObject( util::GetColors());

      return 0;
      }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilColors

  const ExampleClass::EnumType ExampleUtilColors::s_Instance
  (
    GetExamples().AddEnum( ExampleUtilColors())
  );
  
} // namespace bcl
