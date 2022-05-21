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
#include "score/bcl_score_residual_dipolar_coupling_histogram.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_residual_dipolar_coupling_histogram.cpp
  //!
  //! @author weinerbe
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreResidualDipolarCouplingHistogram :
    public ExampleInterface
  {
  public:

    ExampleScoreResidualDipolarCouplingHistogram *Clone() const
    {
      return new ExampleScoreResidualDipolarCouplingHistogram( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    virtual const std::string &GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    int Run() const
    {
      // create RDCs
      storage::Vector< double> exp_values;
      exp_values.PushBack(  8.2700);
      exp_values.PushBack(  10.4900);
      exp_values.PushBack(  9.8700);
      exp_values.PushBack(  9.1500);
      exp_values.PushBack(  3.7000);
      exp_values.PushBack(  6.4600);
      exp_values.PushBack(  -7.5300);
      exp_values.PushBack(  14.3400);
      exp_values.PushBack(  12.490);
      exp_values.PushBack(  1.6500);
      exp_values.PushBack(  -5.1500);
      exp_values.PushBack(  2.7800);
      exp_values.PushBack(  22.9600);
      exp_values.PushBack(  -14.4300);
      exp_values.PushBack(  5.4300);
      exp_values.PushBack(  21.5500);
      exp_values.PushBack(  -6.9800);
      exp_values.PushBack(  -5.2500);
      exp_values.PushBack(  3.8900);
      exp_values.PushBack(  12.4500);
      exp_values.PushBack(  -4.8900);
      exp_values.PushBack(  2.6500);
      exp_values.PushBack(  12.560);
      exp_values.PushBack(  1.2300);
      exp_values.PushBack(  -5.6300);
      exp_values.PushBack(  8.9500);
      exp_values.PushBack(  -13.2200);
      exp_values.PushBack(  -14.4500);

      storage::Vector< double> calc_values;
      calc_values.PushBack( 8.6605);
      calc_values.PushBack( 10.3457);
      calc_values.PushBack( 9.3999);
      calc_values.PushBack( 9.1155);
      calc_values.PushBack( 4.3006);
      calc_values.PushBack( 6.2136);
      calc_values.PushBack( -7.7609);
      calc_values.PushBack( 22.5605);
      calc_values.PushBack( 2.0457);
      calc_values.PushBack( 9.9999);
      calc_values.PushBack( -21.2155);
      calc_values.PushBack( 15.3006);
      calc_values.PushBack( 5.4236);
      calc_values.PushBack( -7.1002);
      calc_values.PushBack( 9.2515);
      calc_values.PushBack( 5.3695);
      calc_values.PushBack( 8.2346);
      calc_values.PushBack( 12.5454);
      calc_values.PushBack( 12.6545);
      calc_values.PushBack( 9.8524);
      calc_values.PushBack( -6.2532);
      calc_values.PushBack( 20.3133);
      calc_values.PushBack( 1.6262);
      calc_values.PushBack( -14.5632);
      calc_values.PushBack( -15.9582);
      calc_values.PushBack( -18.2666);
      calc_values.PushBack( 4.2589);
      calc_values.PushBack( -15.2222);

      nmr::RDCContainer rdc_container( exp_values, calc_values);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      score::ResidualDipolarCouplingHistogram default_constr;
      const double calculated_histogram( default_constr( rdc_container));

      // create const double "expected_histogram_bcl" and initialize with the value calculated by bcl
      const double expected_histogram_bcl( -2.45007);

       BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_histogram_bcl, calculated_histogram),
        "expected bcl histogram is " + util::Format()( expected_histogram_bcl) + " but calculated histogram is "
        + util::Format()( calculated_histogram) + " deviation is " +
        util::Format()( calculated_histogram - expected_histogram_bcl)
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreResidualDipolarCouplingHistogram

  const ExampleClass::EnumType ExampleScoreResidualDipolarCouplingHistogram::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreResidualDipolarCouplingHistogram())
  );
  
} // namespace bcl
