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
#include "score/bcl_score_residual_dipolar_coupling_q_value.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_score_residual_dipolar_coupling_q_value.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleScoreResidualDipolarCouplingQValue :
    public ExampleInterface
  {
  public:

    ExampleScoreResidualDipolarCouplingQValue *Clone() const
    {
      return new ExampleScoreResidualDipolarCouplingQValue( *this);
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
      const storage::Vector< double> calc_values
      (
        storage::Vector< double>::Create( 8.27, 10.49, 9.87, 9.15, 3.7, 6.46, -7.53)
      );
      const storage::Vector< double> exp_values
      (
        storage::Vector< double>::Create( 8.6605, 10.3457, 9.3999, 9.1155, 4.3006, 6.2136, -7.7609)
      );
      nmr::RDCContainer rdc_container( calc_values, exp_values);

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      // test default constructor
      score::ResidualDipolarCouplingQValue default_constr;
      const double calculated_q_value( default_constr( rdc_container));

      // create const double "expected_q_value" and initialize with the value given by dipocoup
      const double expected_q_value_dipocoup( -0.957047);
      // create const double "expected_q_value_bcl" and initialize with the value calculated by bcl
      const double expected_q_value_bcl( -0.957047);
      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance( expected_q_value_dipocoup, calculated_q_value, 0.03),
        "expected dipocoup q value is " + util::Format()( expected_q_value_dipocoup) + " but calculated q value is "
        + util::Format()( calculated_q_value) + " deviation is " +
        util::Format()( calculated_q_value - expected_q_value_dipocoup)
      );
      BCL_Example_Check
      (
        math::EqualWithinTolerance( expected_q_value_bcl, calculated_q_value),
        "expected bcl q value is " + util::Format()( expected_q_value_bcl) + " but calculated q value is "
        + util::Format()( calculated_q_value) + " deviation is " +
        util::Format()( calculated_q_value - expected_q_value_bcl)
      );
      BCL_MessageStd
      (
        "expected q dipocoup value is " + util::Format()( expected_q_value_dipocoup) +
        " and expected bcl calculated value is " + util::Format()( expected_q_value_bcl) +
        " calculated q value is " + util::Format()( calculated_q_value)
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleScoreResidualDipolarCouplingQValue

  const ExampleClass::EnumType ExampleScoreResidualDipolarCouplingQValue::s_Instance
  (
    GetExamples().AddEnum( ExampleScoreResidualDipolarCouplingQValue())
  );

} // namespace bcl
