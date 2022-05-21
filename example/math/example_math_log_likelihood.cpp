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
#include "math/bcl_math_log_likelihood.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_running_average_sd.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_log_likelihood.cpp
  //!
  //! @author bitterd
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathLogLikelihood :
    public ExampleInterface
  {
  public:

      ExampleMathLogLikelihood *Clone() const
    {
      return new ExampleMathLogLikelihood( *this);
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
      // create a fictional mean and standard deviation for the scoring
      math::RunningAverageSD< double> mean_sd;
      mean_sd += 3;
      mean_sd += 4;
      mean_sd += 5;
      // Mean should now be 4 and Sigma should be 0.8165

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      math::LogLikelihood loglikelihood_default;

      //construct with mean and standard deviation
      math::LogLikelihood loglikelihood( mean_sd.GetAverage(), mean_sd.GetStandardDeviation());

      BCL_MessageStd
      (
        "Statistic Mean calculated by RunningAverageSD with values 3,4 and 5: "
        + util::Format()( mean_sd.GetAverage())
      );

      BCL_MessageStd
      (
        "Standard Deviation calculated by RunningAverageSD with values 3,4 and 5: "
        + util::Format()( mean_sd.GetStandardDeviation())
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance( loglikelihood.GetMean(), mean_sd.GetAverage())
        && math::EqualWithinTolerance( loglikelihood.GetSigma(), mean_sd.GetStandardDeviation()),
        "construction from ave/sd was not successful: " + util::Format()( loglikelihood)
       );

      // copy constructor
      math::LogLikelihood loglikelihood_copy( loglikelihood);

      BCL_Example_Check
      (
        math::EqualWithinTolerance( loglikelihood.GetMean(), mean_sd.GetAverage())
        && math::EqualWithinTolerance( loglikelihood.GetSigma(), mean_sd.GetStandardDeviation()),
        "Copied LogLikelihood did not contain RunningAverageSD: " + util::Format()( loglikelihood_copy)
      );

      // clone
      util::ShPtr< math::LogLikelihood> loglikelihood_cloned( loglikelihood.Clone());

      // check clone
      BCL_Example_Check
      (
        loglikelihood_cloned->GetClassIdentifier() == GetStaticClassName< math::LogLikelihood>() &&
        math::EqualWithinTolerance( loglikelihood_cloned->GetMean(), loglikelihood.GetMean()) &&
        math::EqualWithinTolerance( loglikelihood_cloned->GetSigma(), loglikelihood.GetSigma()),
        "cloned object is not pointing to the correct object"
      );

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      BCL_Example_Check
      (
        GetStaticClassName< math::LogLikelihood>() == "bcl::math::LogLikelihood"
        && loglikelihood_cloned->GetClassIdentifier() == GetStaticClassName< math::LogLikelihood>(),
        "wrong class name"
      );

    ///////////////
    // operators //
    ///////////////

      // calculate z(x) with mean 4 and sigma 1
      const double arg1( 3);
      const double expected_l1( -0.0425237);
      const double loglikelihood1( loglikelihood( arg1));

      BCL_MessageStd
      (
        "Value of z-score with mean 4 sigma " + util::Format()( mean_sd.GetStandardDeviation()) + " at x="
        + util::Format()( arg1) + ": "
        + util::Format()( loglikelihood1)
      );

      const double arg2( 2.0);
      const double expected_l2( -0.000266038);
      const double loglikelihood2( loglikelihood( arg2));

      const double arg3( 5.0);
      const double expected_l3( -3.1788799715);
      const double loglikelihood3( loglikelihood( arg3));

      const double arg4( 3.8);
      const double expected_l4( -0.453370324);
      const double loglikelihood4( loglikelihood( arg4));

      BCL_MessageStd
      (
        "Value of z-score with mean 4 and sigma " + util::Format()( mean_sd.GetStandardDeviation()) + " at x="
        + util::Format()( arg2) + " and x="
        + util::Format()( arg3)
        + util::Format()( loglikelihood2) + "; "+util::Format()( loglikelihood3)
      );

      BCL_Example_Check
      (
        math::EqualWithinTolerance( loglikelihood1, expected_l1, 0.05)
        && math::EqualWithinTolerance( loglikelihood2, expected_l2, 0.05)
        && math::EqualWithinTolerance( loglikelihood3, expected_l3, 0.05)
        && math::EqualWithinTolerance( loglikelihood4, expected_l4, 0.05),
        "Z-Score with mean 4 and sigma 1 does not match at x=" + util::Format()( arg1) + ": "
        + util::Format()( expected_l1) +"!= "+ util::Format()( loglikelihood1)
        + "or at x=" + util::Format()(arg2) + ": "
        + util::Format()( expected_l2) +"!= "+ util::Format()( loglikelihood2)
        + "or at x=" + util::Format()(arg3) + ": "
        + util::Format()( expected_l3) +"!= "+ util::Format()( loglikelihood3)
        + "or at x=" + util::Format()(arg4) + ": "
        + util::Format()( expected_l4) +"!= "+ util::Format()( loglikelihood4)
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      math::LogLikelihood loglikelihood_read;
      BCL_ExampleCheck( TestBCLObjectIOForSymmetry( loglikelihood, loglikelihood_read), true);

    //////////////////////
    // helper functions //
    //////////////////////

      // end
      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleLogLikelihood

  const ExampleClass::EnumType ExampleMathLogLikelihood::s_Instance
  (
    GetExamples().AddEnum( ExampleMathLogLikelihood())
  );

} // namespace bcl
