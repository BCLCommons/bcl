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
#include "math/bcl_math_linear_function.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_range.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_linear_function.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathLinearFunction :
    public ExampleInterface
  {
  public:

    ExampleMathLinearFunction *Clone() const
    {
      return new ExampleMathLinearFunction( *this);
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

      // test default constructor
      const math::LinearFunction default_constr;
      BCL_Example_Check
      (
        !util::IsDefined( default_constr.GetSlope()),
        "Slope should be undefined but is " + util::Format()( default_constr.GetSlope())
      );
      BCL_Example_Check
      (
        !util::IsDefined( default_constr.GetOrdinateIntercept()),
        "OrdinateIntercept should be undefined but is " + util::Format()( default_constr.GetOrdinateIntercept())
      );

      // test constructor taking a slope and a y-intercept
      const double slope( 2.0);
      const double ordinate_intercept( 1.0);
      const math::LinearFunction param_constr( slope, ordinate_intercept);
      BCL_Example_Check
      (
        param_constr.GetSlope() == slope,
        "Slope should be " + util::Format()( slope) + " but is " + util::Format()( param_constr.GetSlope())
      );
      BCL_Example_Check
      (
        param_constr.GetOrdinateIntercept() == ordinate_intercept,
        "OrdinateIntercept should be " + util::Format()( ordinate_intercept) + " but is " +
        util::Format()( param_constr.GetOrdinateIntercept())
      );

      // test LinearRegression calculation constructor
      const double x_values[ 3] = { 1.0, 2.0, 3.0};
      const double y_values[ 3] = { 3.0, 5.0, 7.0};
      const storage::List< double> x_data( 3, x_values);
      const storage::List< double> y_data( 3, y_values);
      const math::LinearFunction regression_constr( x_data.Begin(), x_data.End(), y_data.Begin(), y_data.End(), 3);
      BCL_Example_Check
      (
        regression_constr.GetSlope() == slope,
        "GetSlope function should have returned " + util::Format()( slope) + " but instead returned "
        + util::Format()( regression_constr.GetSlope())
      );
      BCL_Example_Check
      (
        regression_constr.GetOrdinateIntercept() == ordinate_intercept,
        "GetOrdinateIntercept function should have returned " + util::Format()( ordinate_intercept) +
        " but instead returned " + util::Format()( regression_constr.GetOrdinateIntercept())
      );
      // do another test where linear fit is not perfect
      const double x_values_b[ 5] = { 1.0, 2.0, 3.0, 4.0, 5.0};
      const double y_values_b[ 5] = { 4.0, 5.0, 6.0, 9.0, 11.0};
      const storage::List< double> x_data_b( 5, x_values_b);
      const storage::List< double> y_data_b( 5, y_values_b);
      const math::LinearFunction regression_constr_b
      (
        x_data_b.Begin(), x_data_b.End(), y_data_b.Begin(), y_data_b.End(), 5
      );
      const double correct_slope( 1.8);
      const double correct_ordinate_intercept( 1.6);
      BCL_ExampleIndirectCheck( regression_constr_b.GetSlope(), correct_slope, "Constructor from iterators");
      BCL_ExampleIndirectCheck( regression_constr_b.GetOrdinateIntercept(), correct_ordinate_intercept, "Constructor from iterators");

      // check copy constructor
      math::LinearFunction copy_constr( param_constr);
      BCL_Example_Check
      (
        copy_constr.GetSlope() == slope,
        "Slope should be " + util::Format()( slope) + " but is " + util::Format()( copy_constr.GetSlope())
      );
      BCL_Example_Check
      (
        copy_constr.GetOrdinateIntercept() == ordinate_intercept,
        "OrdinateIntercept should be " + util::Format()( ordinate_intercept) + " but is " +
        util::Format()( copy_constr.GetOrdinateIntercept())
      );

      // check clone constructor
      util::ShPtr< math::LinearFunction> clone_constr( param_constr.Clone());
      BCL_Example_Check
      (
        clone_constr->GetSlope() == slope,
        "Slope should be " + util::Format()( slope) + " but is " + util::Format()( clone_constr->GetSlope())
      );
      BCL_Example_Check
      (
        clone_constr->GetOrdinateIntercept() == ordinate_intercept,
        "OrdinateIntercept should be " + util::Format()( ordinate_intercept) + " but is " +
        util::Format()( clone_constr->GetOrdinateIntercept())
      );

      // try constructor from ranges
      math::Range< double> from_range( 0.1, 2.1);
      math::Range< double> to_range( 1.1, 17.1);
      math::LinearFunction linear_from_ranges( from_range, to_range);
      BCL_ExampleIndirectCheck( linear_from_ranges( 0.1), 1.1, "Constructor from ranges");
      BCL_ExampleIndirectCheck( linear_from_ranges( 2.1), 17.1, "Constructor from ranges");

      // try constructor from points that the line should pass through
      math::LinearFunction linear_from_points( 0.5, 1.2, 1.5, 6.4);
      BCL_ExampleIndirectCheck( linear_from_points( 0.5), 1.2, "Constructor from points");
      BCL_ExampleIndirectCheck( linear_from_points( 1.5), 6.4, "Constructor from points");

    /////////////////
    // data access //
    /////////////////

      // check GetStaticClassName
      const std::string correct_static_class_name( "bcl::math::LinearFunction");
      BCL_Example_Check
      (
        GetStaticClassName< math::LinearFunction>() == correct_static_class_name,
        "GetStaticClassName gives " + GetStaticClassName< math::LinearFunction>() + " but should give " +
        correct_static_class_name
      );

      // test GetSlope function
      BCL_Example_Check
      (
        param_constr.GetSlope() == slope,
        "GetSlope function should have returned " + util::Format()( slope) + " but instead returned "
        + util::Format()( param_constr.GetSlope())
      );

      // test GetOrdinateIntercept function
      BCL_Example_Check
      (
        param_constr.GetOrdinateIntercept() == ordinate_intercept,
        "GetOrdinateIntercept function should have returned " + util::Format()( ordinate_intercept) +
        " but instead returned " + util::Format()( param_constr.GetOrdinateIntercept())
      );

    ///////////////
    // operators //
    ///////////////

      // check operator()
      const double y_value( param_constr( 6));
      const double correct_y_value( 13);
      BCL_Example_Check
      (
        y_value == correct_y_value,
        "operator() should have returned " + util::Format()( correct_y_value) +
        " but instead returned " + util::Format()( y_value)
      );

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      WriteBCLObject( copy_constr);

      // read the file back
      math::LinearFunction read_linear_function;
      ReadBCLObject( read_linear_function);

      // make sure that the written LinearFunction and the read-in LinearFunction are the same
      BCL_Example_Check
      (
        copy_constr.GetSlope() == slope,
        "Slope should be " + util::Format()( slope) + " but is " + util::Format()( copy_constr.GetSlope())
      );
      BCL_Example_Check
      (
        copy_constr.GetOrdinateIntercept() == ordinate_intercept,
        "OrdinateIntercept should be " + util::Format()( ordinate_intercept) + " but is " +
        util::Format()( copy_constr.GetOrdinateIntercept())
      );

    //////////////////////
    // helper functions //
    //////////////////////

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleUtilRingElement

  const ExampleClass::EnumType ExampleMathLinearFunction::s_Instance
  (
    GetExamples().AddEnum( ExampleMathLinearFunction())
  );

} // namespace bcl
