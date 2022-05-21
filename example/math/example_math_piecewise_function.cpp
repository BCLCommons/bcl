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
#include "math/bcl_math_piecewise_function.h"

// includes from bcl - sorted alphabetically
#include "math/bcl_math_quadratic_function.h"

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_piecewise_function.cpp
  //!
  //! @author alexanns
  //! @date Nov 6, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathPiecewiseFunction :
    public ExampleInterface
  {
  public:

    ExampleMathPiecewiseFunction *Clone() const
    {
      return new ExampleMathPiecewiseFunction( *this);
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
      const math::PiecewiseFunction default_constr;
      BCL_Example_Check
      (
        default_constr.GetRangesFunctions().IsEmpty(),
        "RangesFunctions should be empty but is " + util::Format()( default_constr.GetRangesFunctions())
      );

      // test constructor taking the list of Pairs of ranges and functions and passing it to the piecewise function
      // create the list to be used in the example

      // make a few Ranges
      const math::Range< double> range_a( math::RangeBorders::e_LeftOpen, 3.0, 4.0, math::RangeBorders::e_RightOpen),
        range_b( math::RangeBorders::e_LeftClosed, 4.0, 6.0, math::RangeBorders::e_RightOpen),
        range_c( math::RangeBorders::e_LeftClosed, 6.0, 10.0, math::RangeBorders::e_RightOpen);

      // make a few Functions
      const util::ShPtr< math::FunctionInterfaceSerializable< double, double> > function_a
      (
        new math::QuadraticFunction( 1.0, 2.0, 3.0)
      ),
        function_b( new math::LinearFunction( 2.0, 3.0)),
        function_c( new math::QuadraticFunction( 1.0, 2.0, 4.0));

      // store the ranges and functions in a list named "param_constr"
      storage::List
      <
        storage::Pair< math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> > >
      > param_constr;
      param_constr.PushBack
      (
        storage::Pair< math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> > >( range_a, function_a)
      ),
      param_constr.PushBack
      (
        storage::Pair< math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> > >( range_b, function_b)
      ),
      param_constr.PushBack
      (
        storage::Pair< math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> > >( range_c, function_c)
      );

      // now create PiecewiseFunction "constr" and initialize with "param_constr"
      const math::PiecewiseFunction constr( param_constr);

      // test piecewise function to ensure that it can get the list
      BCL_Example_Check
      (
        constr.GetRangesFunctions() == param_constr,
        "RangesFunctions should be " + util::Format()( param_constr) + " but is " + util::Format()
        (
          constr.GetRangesFunctions()
        )
      );

      // check copy constructor
      math::PiecewiseFunction copy_constr( constr);
      BCL_Example_Check
      (
        copy_constr.GetRangesFunctions() == param_constr,
        "RangesFunctions should be " + util::Format()( param_constr) + " but is " + util::Format()
        (
          copy_constr.GetRangesFunctions()
        )
      );

      // check clone constructor
      util::ShPtr< math::PiecewiseFunction> clone_constr( constr.Clone());
      BCL_Example_Check
      (
        clone_constr->GetRangesFunctions() == param_constr,
        "RangesFunctions should be " + util::Format()( param_constr) + " but is " + util::Format()
        (
          clone_constr->GetRangesFunctions()
        )
      );

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      const std::string correct_static_class_name( "bcl::math::PiecewiseFunction");
      BCL_Example_Check
      (
        GetStaticClassName< math::PiecewiseFunction>() == clone_constr->GetClassIdentifier() &&
        GetStaticClassName< math::PiecewiseFunction>() == correct_static_class_name,
        "GetClassIdentifier gives " + clone_constr->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

      // test GetRangesFunctions function
      BCL_Example_Check
      (
        constr.GetRangesFunctions() == param_constr,
        "GetXcoord function should have returned " + util::Format()( param_constr) + " but instead returned "
        + util::Format()( constr.GetRangesFunctions())
      );

    ///////////////
    // operators //
    ///////////////

      // check operator()
      const double y_value( constr( 5.2));
      const double correct_y_value( 13.4);
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
      math::PiecewiseFunction read_piecewise_function;
      ReadBCLObject( read_piecewise_function);

      // make sure that the written PiecewiseFunction and the read-in PiecewiseFunction are the same
      BCL_Example_Check
      (
        read_piecewise_function( 5.2) == correct_y_value,
        "RangesFunctions should be " + util::Format()( correct_y_value) + " but is " + util::Format()
        (
          read_piecewise_function( 5.2)
        )
      );

    //////////////////////
    // helper functions //
    //////////////////////

      // test RangesFunctionsValid function
      // check to make sure that the piecewise function stops if ranges overlap
      // create a list with overlapping ranges
      const math::Range< double> overlap_a( math::RangeBorders::e_LeftOpen, 3.0, 4.0, math::RangeBorders::e_RightOpen),
        overlap_b( math::RangeBorders::e_LeftClosed, 3.5, 6.0, math::RangeBorders::e_RightOpen),
        overlap_c( math::RangeBorders::e_LeftClosed, 6.0, 10.0, math::RangeBorders::e_RightOpen);
      storage::List
      <
        storage::Pair< math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> > >
      > param_constr_overlap;
      param_constr_overlap.PushBack
      (
        storage::Pair< math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> > >( overlap_a, function_a)
      ),
      param_constr_overlap.PushBack
      (
        storage::Pair< math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> > >( overlap_b, function_b)
      ),
      param_constr_overlap.PushBack
      (
        storage::Pair< math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> > >( overlap_c, function_c)
      );
      const bool is_valid( math::PiecewiseFunction::RangesFunctionsValid( param_constr_overlap));

      // perform the check that determines if the ranges overlap
      BCL_Example_Check
      (
        !is_valid,
        "RangesFunctionsValid should be false because of overlapping ranges "
      );

      // check to make sure that the piecewise function stops if the ShPtr points to nothing
      // create a ShPtr that points to nothing
      const util::ShPtr< math::FunctionInterfaceSerializable< double, double> > function_nothing_a;
      storage::List
      <
        storage::Pair< math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> > >
      > param_constr_nothing;

      param_constr_nothing.PushBack
      (
        storage::Pair< math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> > >
        (
          range_a, function_nothing_a
        )
      ),
      param_constr_nothing.PushBack
      (
        storage::Pair< math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> > >( range_b, function_b)
      ),
      param_constr_nothing.PushBack
      (
        storage::Pair< math::Range< double>, util::ShPtr< math::FunctionInterfaceSerializable< double, double> > >( range_c, function_c)
      );
      const bool constr_nothing( math::PiecewiseFunction::RangesFunctionsValid( param_constr_nothing));

      // perform the check that determines if the functions point to something
      BCL_Example_Check
      (
        !constr_nothing,
        "RangesFunctionsValid should be false because of nonexistent function"
      );

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

   }; //end example_math_piecewise_function

  const ExampleClass::EnumType ExampleMathPiecewiseFunction::s_Instance
  (
    GetExamples().AddEnum( ExampleMathPiecewiseFunction())
  );

} // namespace bcl

