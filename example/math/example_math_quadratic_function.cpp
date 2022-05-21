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
#include "math/bcl_math_quadratic_function.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_quadratic_function.cpp
  //!
  //! @author akinlr
  //! @date Nov 6, 2010 
  //! @remarks status complete
  //! @remarks reviewed by nobody on 
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathQuadraticFunction :
    public ExampleInterface
  {
  public:

    ExampleMathQuadraticFunction *Clone() const
    {
      return new ExampleMathQuadraticFunction( *this);
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
      const math::QuadraticFunction default_constr;
      BCL_Example_Check
      (
        !util::IsDefined( default_constr.GetX()),
        "Xcoord should be undefined but is " + util::Format()( default_constr.GetX())
      );
      BCL_Example_Check
      (
        !util::IsDefined( default_constr.GetY()),
        "Ycoord should be undefined but is " + util::Format()( default_constr.GetY())
      );
      BCL_Example_Check
      (
        !util::IsDefined( default_constr.GetA()),
        "A should be undefined but is " + util::Format()( default_constr.GetA())
      );

      //test constructor taking A, B, and C components from the general form
      const double avalue( 2.0);
      BCL_ExampleCheck( math::QuadraticFunction( 2.0, 3.0, 4.0).GetX(), -0.75);
      BCL_ExampleCheck( math::QuadraticFunction( 2.0, 3.0, 4.0).GetY(), 2.875);
      BCL_ExampleCheck( math::QuadraticFunction( avalue, 3.0, 4.0).GetA(), avalue);

      // test constructor using the x and y coordinates and the a variable from standard form
      storage::VectorND< 2, double> xy_coords( 2.0, 1.0);
      const double a_variable( 3.0);
      const math::QuadraticFunction param_constr( xy_coords, a_variable);
      BCL_ExampleCheck( math::QuadraticFunction( xy_coords, a_variable).GetX(), xy_coords.First());
      BCL_ExampleCheck( math::QuadraticFunction( xy_coords, a_variable).GetY(), xy_coords.Second());
      BCL_ExampleCheck( math::QuadraticFunction( xy_coords, a_variable).GetA(), a_variable);

      // check clone constructor
      util::ShPtr< math::QuadraticFunction> clone_constr( param_constr.Clone());
      BCL_Example_Check
      (
        clone_constr->GetX() == xy_coords.First(),
        "Xcoord should be " + util::Format()( xy_coords.First()) + " but is " + util::Format()( clone_constr->GetX())
      );
      BCL_Example_Check
      (
        clone_constr->GetY() == xy_coords.Second(),
        "Ycoord should be " + util::Format()( xy_coords.Second()) + " but is " +
        util::Format()( clone_constr->GetY())
      );
      BCL_Example_Check
      (
        clone_constr->GetA() == a_variable,
        "A should be " + util::Format()( a_variable) + " but is " +
        util::Format()( clone_constr->GetA())
      );

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      const std::string correct_static_class_name( "bcl::math::QuadraticFunction");
      BCL_ExampleCheck( GetStaticClassName< math::QuadraticFunction>(), clone_constr->GetClassIdentifier());

      // test GetX function
      BCL_Example_Check
      (
        param_constr.GetX() == xy_coords.First(),
        "GetX function should have returned " + util::Format()( xy_coords.First()) + " but instead returned "
        + util::Format()( param_constr.GetX())
      );

      // test GetY function
      BCL_Example_Check
      (
        param_constr.GetY() == xy_coords.Second(),
        "GetY function should have returned " + util::Format()( xy_coords.Second()) +
        " but instead returned " + util::Format()( param_constr.GetY())
      );

      // test GetA function
      BCL_Example_Check
      (
        param_constr.GetA() == a_variable,
        "GetA function should have returned " + util::Format()( a_variable) +
        " but instead returned " + util::Format()( param_constr.GetA())
      );

    ///////////////
    // operators //
    ///////////////

      // check operator()
      const double y_value( param_constr( 6));
      const double correct_y_value( 49);
      BCL_ExampleCheck( y_value, correct_y_value);

    ////////////////
    // operations //
    ////////////////

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      WriteBCLObject( param_constr);

      // read the file back
      math::QuadraticFunction read_quadratic_function;
      ReadBCLObject( read_quadratic_function);

      // make sure that the written LinearFunction and the read-in QuadraticFunction are the same
      BCL_Example_Check
      (
        read_quadratic_function.GetX() == xy_coords.First(),
        "Xcoord should be " + util::Format()( xy_coords.First()) +
        " but is " + util::Format()( read_quadratic_function.GetX())
      );
      BCL_Example_Check
      (
        read_quadratic_function.GetY() == xy_coords.Second(),
        "Ycoord should be " + util::Format()( xy_coords.Second()) +
        " but is " + util::Format()( read_quadratic_function.GetY())
      );
      BCL_Example_Check
      (
        read_quadratic_function.GetA() == a_variable,
        "A should be " + util::Format()( a_variable) + " but is " +
        util::Format()( read_quadratic_function.GetA())
      );

    //////////////////////
    // helper functions //
    //////////////////////

      //test CompletingTheSquare
      const double a_value( 2.0);
      const double b_value( 6.0);
      const double c_value( 13.0);
      storage::VectorND< 3, double> complete( math::QuadraticFunction::CompletingTheSquare( a_value, b_value, c_value));
      const double expected_av( 2.0);
      const double expected_x( -1.5);
      const double expected_y( 8.5);

      BCL_Example_Check
      (
        expected_av == complete.Third(),
        "First should be " + util::Format()( expected_av) + " but is " + util::Format()( complete.Third())
      );
      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance( expected_x, complete.First()),
        "Second should be " + util::Format()( expected_x) + " but is " + util::Format()( complete.First())
      );
      BCL_Example_Check
      (
        math::EqualWithinAbsoluteTolerance( expected_y, complete.Second()),
        "Third should be " + util::Format()( expected_y) + " but is " + util::Format()( complete.Second())
      );

      // test with the a_var as 0 to make sure all the values returned are undefined
      const double a_value_undef( 0.0);
      storage::VectorND< 3, double> gen_undef
      (
        math::QuadraticFunction::CompletingTheSquare( a_value_undef, b_value, c_value)
      );
      BCL_ExampleIndirectCheck( util::IsDefined( gen_undef.First()), false, "Completing the Square");
      BCL_ExampleIndirectCheck( util::IsDefined( gen_undef.Second()), false, "Completing the Square");
      BCL_ExampleIndirectCheck( util::IsDefined( gen_undef.Third()), false, "Completing the Square");

      //test Distribution (i.e. converting from standard form to general form)
      {
        const double a_var( 2.0);
        const double x_var( 3.0);
        const double y_var( 4.0);
        storage::VectorND< 3, double> gen( math::QuadraticFunction::Distribution( x_var, y_var, a_var));
        const double expected_a( 2.0);
        const double expected_b( -12.0);
        const double expected_c( 22.0);
        BCL_Example_Check
        (
          expected_a == gen.First(),
          "First should be " + util::Format()( expected_a) + " but is " + util::Format()( gen.First())
        );
        BCL_Example_Check
        (
          math::EqualWithinAbsoluteTolerance( expected_b, gen.Second()),
          "Second should be " + util::Format()( expected_b) + " but is " + util::Format()( gen.Second())
        );
        BCL_Example_Check
        (
          math::EqualWithinAbsoluteTolerance( expected_c, gen.Third()),
          "Third should be " + util::Format()( expected_c) + " but is " + util::Format()( gen.Third())
        );

      }

      // test GetDerivative
      {
        const util::ShPtr< math::FunctionInterfaceSerializable< double, double> > derivative( param_constr.GetDerivative());
        const double expected_ans( 6.0);
        BCL_Example_Check
        (
          expected_ans == derivative->operator()( 3.0),
          "GetDerivative should return" + util::Format()( expected_ans) + " but is " +
          util::Format()( derivative->operator()( 3.0))
        );
      }

      // test GetRoot
      BCL_ExampleCheckWithinTolerance( math::QuadraticFunction( 1.0, -3.0, 2.0).GetRoot().First(), 2.0, 0.001);
      BCL_ExampleCheckWithinTolerance( math::QuadraticFunction( 1.0, -3.0, 2.0).GetRoot().Second(), 1.0, 0.001);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end Exampleexample_name

  const ExampleClass::EnumType ExampleMathQuadraticFunction::s_Instance
  (
    GetExamples().AddEnum( ExampleMathQuadraticFunction())
  );
  
} // namespace bcl
