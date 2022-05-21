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
#include "math/bcl_math_z_score.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_running_average_sd.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_z_score.cpp
  //!
  //! @author bitterd
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathZScore :
    public ExampleInterface
  {
  public:

    ExampleMathZScore *Clone() const
    {
      return new ExampleMathZScore( *this);
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
      math::RunningAverageSD< double> mean_sd;
      mean_sd += 3;
      mean_sd += 4;
      mean_sd += 5;
      // Mean should now be 4 and Sigma should be 0.8165

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //construct with default constructor
      math::ZScore zscore_default;

      //construct with mean and standard deviation
      math::ZScore zscore( mean_sd.GetAverage(), mean_sd.GetStandardDeviation());

      BCL_Example_Check
      (
        math::EqualWithinTolerance( zscore.GetMean(), mean_sd.GetAverage())
        && math::EqualWithinTolerance( zscore.GetSigma(), mean_sd.GetStandardDeviation()),
        "construction from ave/sd was not successful: " + util::Format()( zscore)
      );

      // clone
      util::ShPtr< math::ZScore> zscore_cloned( zscore.Clone());

      // check clone
      BCL_Example_Check
      (
        zscore_cloned->GetClassIdentifier() == GetStaticClassName< math::ZScore>() &&
        math::EqualWithinTolerance( zscore_cloned->GetMean(), zscore.GetMean()) &&
        math::EqualWithinTolerance( zscore_cloned->GetSigma(), zscore.GetSigma()),
        "cloned object is not pointing to the correct object"
      );

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      BCL_Example_Check
      (
        GetStaticClassName< math::ZScore>() == "bcl::math::ZScore"
        && zscore_cloned->GetClassIdentifier() == GetStaticClassName< math::ZScore>(),
        "wrong class name"
      );

    ///////////////
    // operators //
    ///////////////

      // calculate z(x) with mean 4 and sigma 1
      const double arg1( 3);
      const double expected_z1( -1.224745);
      const double zscore1( zscore( arg1));

      BCL_MessageStd
      (
        "Value of z-score with mean 4 sigma " + util::Format()( mean_sd.GetStandardDeviation()) + " at x="
        + util::Format()( arg1) + ": "
        + util::Format()( zscore1)
      );

      const double arg2( -7);
      const double expected_z2( -13.47219);
      const double zscore2( zscore( arg2));

      const double arg3( 5);
      const double expected_z3( 1.224745);
      const double zscore3( zscore( arg3));

      const double arg4( 3.8);
      const double expected_z4( -0.2449488);
      const double zscore4( zscore( arg4));

      BCL_MessageStd
      (
        "Value of z-score with mean 4 and sigma " + util::Format()( mean_sd.GetStandardDeviation()) + " at x="
        + util::Format()( arg2) + " and x="
        + util::Format()( arg3)
        + util::Format()( zscore2) + "; "+util::Format()( zscore3)
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance( zscore1, expected_z1)
        && math::EqualWithinTolerance( zscore2, expected_z2)
        && math::EqualWithinTolerance( zscore3, expected_z3)
        && math::EqualWithinTolerance( zscore4, expected_z4),
        "Z-Score with mean 4 and sigma 1 does not match at x=" + util::Format()( arg1) + ": "
        + util::Format()( expected_z1) +"!= "+ util::Format()( zscore1)
        + "or at x=" + util::Format()(arg2) + ": "
        + util::Format()( expected_z2) +"!= "+ util::Format()( zscore2)
        + "or at x=" + util::Format()(arg3) + ": "
        + util::Format()( expected_z3) +"!= "+ util::Format()( zscore3)
        + "or at x=" + util::Format()(arg4) + ": "
        + util::Format()( expected_z4) +"!= "+ util::Format()( zscore4)
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // write to File
      WriteBCLObject( zscore);

      // read from file
      math::ZScore zscore_read;
      ReadBCLObject( zscore_read);

      BCL_ExampleCheck( TestBCLObjectIOForSymmetry( zscore, zscore_read), true);

    //////////////////////
    // helper functions //
    //////////////////////

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleZScore

  const ExampleClass::EnumType ExampleMathZScore::s_Instance
  (
    GetExamples().AddEnum( ExampleMathZScore())
  );

} // namespace bcl
