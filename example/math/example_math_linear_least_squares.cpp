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
#include "math/bcl_math_linear_least_squares.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_template_instantiations.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_linear_least_squares.cpp
  //!
  //! @author mendenjl
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathLinearLeastSquares :
    public ExampleInterface
  {
  public:

    ExampleMathLinearLeastSquares *Clone() const
    {
      return new ExampleMathLinearLeastSquares( *this);
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

      // default constructor
      math::LinearLeastSquares lss_default;

      // copy constructor
      math::LinearLeastSquares lss_copy;

      // clone
      util::ShPtr< math::LinearLeastSquares> lss_clone( lss_copy.Clone());

    /////////////////
    // data access //
    /////////////////

      // class identifiers
      BCL_ExampleCheck( lss_clone->GetClassIdentifier(), GetStaticClassName< math::LinearLeastSquares>());

    ///////////////
    // operators //
    ///////////////

      double x_values[ 8] = { 1.0, 0.0, 1.0, 1.0, 1.0, 2.0, 1.0, 3.0};
      linal::Matrix< double> independent_data( 4, 2, x_values);
      double y_values[ 4] = { 3.0, 5.0, 7.0, 9.0};
      linal::Vector< double> dependent_data( 4, y_values);

      linal::Vector< double> parameters
      (
        lss_default.SolutionAndChiSquared( independent_data, dependent_data).First()
      );
      BCL_ExampleCheckWithinAbsTolerance
      (
        lss_default.SolutionAndChiSquared( independent_data, dependent_data).First(),
        linal::MakeVector< double>( 3, 2),
        0.0001
      );
      BCL_ExampleCheckWithinAbsTolerance
      (
        lss_default.SolutionAndChiSquared( independent_data, dependent_data).Second(),
        0.0,
        0.0001
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathLinearLeastSquares

  const ExampleClass::EnumType ExampleMathLinearLeastSquares::s_Instance
  (
    GetExamples().AddEnum( ExampleMathLinearLeastSquares())
  );

} // namespace bcl
