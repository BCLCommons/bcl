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
#include "math/bcl_math_angle.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_angle.cpp
  //! @brief The math angle class converts between degrees and radians and shows the unit name
  //!
  //! @author woetzen, putnamdk
  //! @date Aug 4, 2011
  //! @remarks status complete
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathAngle :
    public ExampleInterface
  {
  public:

    ExampleMathAngle *Clone() const
    {
      return new ExampleMathAngle( *this);
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

      // no constructor

    /////////////////
    // data access //
    /////////////////

      BCL_MessageStd
      (
        "this char can be used to add units for degree: "
        + std::string( 1, math::Angle::s_DegreeChar)
      );

      BCL_MessageStd
      (
        "this string can be inserted into gnuplot titles, if a degree char is desired: "
        + math::Angle::s_DegreeSymbolGnuplot
      );

    ///////////////
    // operators //
    ///////////////

    ////////////////
    // operations //
    ////////////////

      std::string radian( math::Angle::GetUnitName( math::Angle::e_Radian));
      BCL_ExampleCheck( radian, "radian");

      std::string degree( math::Angle::GetUnitName( math::Angle::e_Degree));
      BCL_ExampleCheck( degree, "degree");

      std::string undefined( math::Angle::GetUnitName( math::Angle::e_Undefined));
      BCL_ExampleCheck( undefined, "undefined");

      std::string angle( math::Angle::GetUnitName( math::Angle::s_NumberUnitTypes));
      BCL_ExampleCheck( angle, "Angle");

      // convert degree into radians
      BCL_ExampleCheckWithinTolerance( math::g_Pi, math::Angle::Radian( 180.0), 0.001);

      // convert radians to degree
      BCL_ExampleCheckWithinTolerance( 180.0, math::Angle::Degree( math::g_Pi), 0.001);

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathAngle

  const ExampleClass::EnumType ExampleMathAngle::s_Instance
  (
    GetExamples().AddEnum( ExampleMathAngle())
  );

} // namespace bcl
