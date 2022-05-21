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
#include "math/bcl_math_gaussian_function.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_gaussian_function.cpp
  //!
  //! @author woetzen
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathGaussianFunction :
    public ExampleInterface
  {
  public:

    ExampleMathGaussianFunction *Clone() const
    { return new ExampleMathGaussianFunction( *this);}

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
      static const double mean_1( 1);
      static const double sigma_2( 2);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // default constructor
      math::GaussianFunction gaussian_default;

      // construct from mean and sigma (standard deviation)
      math::GaussianFunction gaussian_1_2( mean_1, sigma_2);

      BCL_Example_Check
      (
        gaussian_1_2.GetMean() == mean_1
        && gaussian_1_2.GetSigma() == sigma_2,
        "construction from mean 1 and sigma 2 was not successful: " + util::Format()( gaussian_1_2)
      );

      // copy constructor
      math::GaussianFunction gaussian_copied( gaussian_1_2);

      BCL_Example_Check
      (
        gaussian_copied.GetMean() == mean_1
        && gaussian_copied.GetSigma() == sigma_2,
        "copied gaussian did not contain mean 1 and sigma 2: " + util::Format()( gaussian_copied)
      );

      // clone
      util::ShPtr< math::FunctionInterfaceSerializable< double, double> > gaussian_cloned( gaussian_1_2.Clone());

    /////////////////
    // data access //
    /////////////////

      // class identifiers
      BCL_MessageStd( "class name: " + gaussian_cloned->GetClassIdentifier());
      BCL_Example_Check
      (
        GetStaticClassName< math::GaussianFunction>() == "bcl::math::GaussianFunction"
        && gaussian_cloned->GetClassIdentifier() == GetStaticClassName< math::GaussianFunction>(),
        "incorrect class name"
      );

    ///////////////
    // operators //
    ///////////////

      //! calculate f(x)
      const double y_1( gaussian_1_2( mean_1));
      BCL_MessageStd
      (
        "function value of a gaussian with mean 1 and sigma 2 at x = 1: " + util::Format()( y_1)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( y_1, mean_1),
        "incorrect result for x = 1: " + util::Format()( y_1) + " != " + util::Format()( mean_1)
      );

      const double y_m5( gaussian_1_2( mean_1 - 3 * sigma_2));
      const double y_7( gaussian_1_2( mean_1 + 3 * sigma_2));
      const double y_3sigma_expected( 0.011109);
      BCL_MessageStd
      (
        "function value of a gaussian with mean 1 and sigma 2 at x = -5 and x = 7: "
        + util::Format()( y_m5) + "; " + util::Format()( y_7)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( y_m5, y_3sigma_expected)
        && math::EqualWithinTolerance( y_7, y_3sigma_expected),
        "incorrect result for x = -5 and x = 7: " + util::Format()( y_m5) +
        " != " + util::Format()( y_7) + " != " + util::Format()( y_3sigma_expected)
      );

    ////////////////
    // operations //
    ////////////////

      // get the range for n sigma
      const math::Range< double> range_mean1_3sigma( gaussian_1_2.GetRange( 3));
      BCL_Example_Check
      (
        range_mean1_3sigma.GetMin() == mean_1 - 3 * sigma_2
        && range_mean1_3sigma.GetMax() == mean_1 + 3 * sigma_2,
        "incorrect range for three sigma: " + util::Format()( range_mean1_3sigma)
      );

      // get the normalization factor to a get a distribution with the integral 1
      const double norm_factor_2( gaussian_1_2.GetNormalizationParameter());
      static const double expected_factor( double( 1) / ( sigma_2 * math::Sqrt( 2 * math::g_Pi)));
      BCL_Example_Check
      (
        expected_factor == norm_factor_2,
        "incorrect normalization parameter: " + util::Format()( norm_factor_2) + " != " + util::Format()( expected_factor)
      );

    //////////////////////
    // input and output //
    //////////////////////

      // write to file
      WriteBCLObject( gaussian_1_2);
      ReadBCLObject( gaussian_default);

      BCL_Example_Check
      (
        gaussian_default.GetMean() == mean_1
        && gaussian_default.GetSigma() == sigma_2,
        "read gaussian is different from written gaussian function: " + util::Format()( gaussian_default)
      );

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; // end ExampleMathGaussianFunction

  const ExampleClass::EnumType ExampleMathGaussianFunction::s_Instance
  (
    GetExamples().AddEnum( ExampleMathGaussianFunction())
  );

} // namespace bcl
