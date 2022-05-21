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
#include "math/bcl_math_assignment_unary_standard.h"

// includes from bcl - sorted alphabetically

namespace bcl
{
  template< class t_DataType>
  t_DataType &DoAssignmentOp
  (
    t_DataType &VALUE,
    const math::AssignmentUnaryInterface &OPERATION
  )
  {
    OPERATION( VALUE);
    return VALUE;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_assignment_unary_standard.cpp
  //!
  //! @author mendenjl
  //! @date Mar 11, 2013
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathAssignmentUnaryStandard :
    public ExampleInterface
  {
  public:

    ExampleMathAssignmentUnaryStandard *Clone() const
    {
      return new ExampleMathAssignmentUnaryStandard( *this);
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

      const math::AssignmentUnaryStandard cos_op( &cos);
      const math::AssignmentUnaryStandard log_op( &log);

    /////////////////
    // data access //
    /////////////////

      // test class identifier
      BCL_ExampleCheck( math::AssignmentUnaryStandard( &cos).GetAlias(), "Cos");
      BCL_ExampleCheck( math::AssignmentUnaryStandard( &log).GetAlias(), "Ln");

    ///////////////
    // operators //
    ///////////////

      double new_value( 0.0);

      BCL_ExampleCheck
      (
        DoAssignmentOp( new_value = 1.5, math::AssignmentUnaryStandard( &cos)),
        cos( 1.5)
      );

      // test plus equals.
      BCL_ExampleCheck
      (
        DoAssignmentOp( new_value = 1.5, math::AssignmentUnaryStandard( &log)),
        log( 1.5)
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
  }; //end ExampleMathAssignmentUnaryStandard

  const ExampleClass::EnumType ExampleMathAssignmentUnaryStandard::s_Instance
  (
    GetExamples().AddEnum( ExampleMathAssignmentUnaryStandard())
  );
} // namespace bcl

