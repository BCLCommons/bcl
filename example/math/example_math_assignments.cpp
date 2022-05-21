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
#include "math/bcl_math_assign.h"
#include "math/bcl_math_assignment_by_comparison.h"
#include "math/bcl_math_assignments.h"
#include "math/bcl_math_divide_equals.h"
#include "math/bcl_math_minus_equals.h"
#include "math/bcl_math_mod_equals.h"
#include "math/bcl_math_plus_equals.h"
#include "math/bcl_math_power_equals.h"
#include "math/bcl_math_times_equals.h"

// includes from bcl - sorted alphabetically
#include "linal/bcl_linal_vector_operations.h"

namespace bcl
{
  template< class t_DataType>
  t_DataType &DoAssignmentOp
  (
    t_DataType &LHS,
    const math::AssignmentOperationInterface< t_DataType> &OPERATION,
    const t_DataType &RHS
  )
  {
    OPERATION( LHS, RHS);
    return LHS;
  }

  template< class t_DataType>
  t_DataType &DoAssignmentOpViaEnum
  (
    t_DataType &LHS,
    const std::string &OPERATION,
    const t_DataType &RHS
  )
  {
    BCL_MessageStd( "OPERATION: " + OPERATION);
    util::Implementation< math::AssignmentOperationInterface< t_DataType> >( OPERATION)->operator()( LHS, RHS);
    return LHS;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_assignments.cpp
  //!
  //! @author mendenjl
  //! @date August 01, 2010
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathAssignments :
    public ExampleInterface
  {
  public:

    ExampleMathAssignments *Clone() const
    {
      return new ExampleMathAssignments( *this);
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

      math::Assign< float>       seq_float;
      math::PlusEquals< float>   peq_float;
      math::MinusEquals< float>  meq_float;
      math::TimesEquals< float>  teq_float;
      math::DivideEquals< float> deq_float;
      math::PowerEquals< float>  poweq_float;
      math::ModEquals< float>    modeq_float;
      math::MinusEquals< linal::Vector< float> > peq_vect_float;

    /////////////////
    // data access //
    /////////////////

      // test class identifier
      BCL_ExampleCheck( math::Assign< float>().GetClassIdentifier(), "bcl::math::Assign<float>");
      BCL_ExampleCheck( math::PlusEquals< float>().GetClassIdentifier(), "bcl::math::PlusEquals<float>");
      BCL_ExampleCheck( math::MinusEquals< float>().GetClassIdentifier(), "bcl::math::MinusEquals<float>");
      BCL_ExampleCheck( math::TimesEquals< float>().GetClassIdentifier(), "bcl::math::TimesEquals<float>");
      BCL_ExampleCheck( math::DivideEquals< float>().GetClassIdentifier(), "bcl::math::DivideEquals<float>");
      BCL_ExampleCheck( math::PowerEquals< float>().GetClassIdentifier(), "bcl::math::PowerEquals<float>");
      BCL_ExampleCheck
      (
        math::PlusEquals< std::string>().GetClassIdentifier(),
        "bcl::math::PlusEquals<std::string>"
      );

    ///////////////
    // operators //
    ///////////////

      float new_value( 0.0);

      BCL_ExampleCheck
      (
        DoAssignmentOpViaEnum( new_value = 1.5, "+", 10.0f),
        11.5
      );

      // test plus equals.
      BCL_ExampleCheck
      (
        DoAssignmentOpViaEnum( new_value = 11.5, "-", 10.0f),
        1.5
      );

      // test times equals.
      BCL_ExampleCheck
      (
        DoAssignmentOpViaEnum( new_value = 1.5, "*", 2.0f),
        2.0 * 1.5
      );

      BCL_ExampleCheck
      (
        DoAssignmentOpViaEnum( new_value = 1.5, "/", 2.0f),
        1.5 / 2.0
      );

      BCL_ExampleCheck
      (
        DoAssignmentOpViaEnum( new_value = 1.5, "^", 3.0f),
        1.5 * 1.5 * 1.5
      );

      // test assignment.
      BCL_ExampleCheck
      (
        DoAssignmentOp( new_value = 76.9, math::Assign< float>(), 10.0f),
        10.0
      );

      // test power equals.
      BCL_ExampleCheck
      (
        DoAssignmentOpViaEnum( new_value = 3.0, "^", 2.0f),
        9.0
      );

      // test mod equals
      BCL_ExampleCheck
      (
        DoAssignmentOpViaEnum( new_value = 100.0, "%", 8.0f),
        fmod( 100.0, 8.0)
      );

      math::Assign< std::string>     seq_string;
      math::PlusEquals< std::string> peq_string;

      const std::string bcl_str( "BCL");
      const std::string bcl_long_str( " - BioChemical Library");

      std::string new_str;

      // test minus equals
      BCL_Example_Check
      (
        DoAssignmentOp( new_str, seq_string, bcl_str) == bcl_str,
        "Assignment operation performing string = \"BCL\"failed"
      );

      new_str = bcl_str;

      BCL_Example_Check
      (
        DoAssignmentOp( new_str, peq_string, bcl_long_str) == bcl_str + bcl_long_str,
        "Assignment operation performing \"BCL\" += \" - BioChemical Library\" failed"
      );

      BCL_ExampleCheck
      (
        DoAssignmentOp( new_value = 1.5f, math::AssignmentByComparison< float, std::greater_equal>(), 10.0f),
        0.0f
      );

      BCL_ExampleCheck
      (
        DoAssignmentOp( new_value = 10.0f, math::AssignmentByComparison< float, std::greater_equal>(), 10.0f),
        1.0f
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
    } // Run

    static const ExampleClass::EnumType s_Instance;
  }; //end ExampleMathAssignments

  const ExampleClass::EnumType ExampleMathAssignments::s_Instance
  (
    GetExamples().AddEnum( ExampleMathAssignments())
  );
} // namespace bcl

