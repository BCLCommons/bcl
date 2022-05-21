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
#include "math/bcl_math_trigonometric_transition.h"

// includes from bcl - sorted alphabetically

// external includes - sorted alphabetically

namespace bcl
{
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //!
  //! @example example_math_trigonometric_transition.cpp
  //!
  //! @author akinlr
  //! @date Jul 4, 2011
  //! @remarks status complete
  //! @remarks reviewed by nobody on
  //!
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  class ExampleMathTrigonometricTransition :
    public ExampleInterface
  {
  public:

      ExampleMathTrigonometricTransition *Clone() const
      {
        return new ExampleMathTrigonometricTransition( *this);
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
      const math::TrigonometricTransition default_constr;

      //test constructor taking A, B, and C components from the general form
      const double x0_value( -12.25);
      const double x1_value( 2.75);
      const double y0_value( 0.0);
      const double y1_value( -1.0);
      const math::TrigonometricTransition param_constr_a( x0_value, x1_value, y0_value, y1_value);

      // check copy constructor
      math::TrigonometricTransition copy_constr( param_constr_a);

      // check clone constructor
      util::ShPtr< math::TrigonometricTransition> clone_constr( param_constr_a.Clone());

    /////////////////
    // data access //
    /////////////////

      // check GetClassIdentifier
      const std::string correct_static_class_name( "bcl::math::TrigonometricTransition");
      BCL_Example_Check
      (
        GetStaticClassName< math::TrigonometricTransition>() == clone_constr->GetClassIdentifier() &&
        GetStaticClassName< math::TrigonometricTransition>() == correct_static_class_name,
        "GetClassIdentifier gives " + clone_constr->GetClassIdentifier() + " but should give " +
        correct_static_class_name
      );

    ///////////////
    // operators //
    ///////////////

      // check operator()
      const double y_value( param_constr_a( 1));
      const double correct_y_value( -0.9667902132);
      BCL_Example_Check
      (
        math::EqualWithinTolerance( y_value, correct_y_value),
        "operator() should have returned " + util::Format()( correct_y_value) +
        " but instead returned " + util::Format()( y_value)
      );

    ////////////////
    // operations //
    ////////////////

      // if message level is debug
      if( util::GetMessenger().IsSmallerEqualCurrentMessageLevel( util::Message::e_Debug))
      {
        // iterate over distances
        for( double i( 10); i >= -15; i -= 0.1)
        {
          if( i >= x0_value && i <= x1_value)
          // write out operator answer
          BCL_MessageDbg
          (
            util::Format()( i) + "\t" + util::Format()( param_constr_a( i))
          );
        }
      }

    //////////////////////
    // input and output //
    //////////////////////

      // test Write
      WriteBCLObject( copy_constr);

      // read the file back
      math::TrigonometricTransition read_trigonometric_transition;
      ReadBCLObject( read_trigonometric_transition);

      return 0;
    }

    static const ExampleClass::EnumType s_Instance;

  }; //end ExampleMathTrigonometricTransition

  const ExampleClass::EnumType ExampleMathTrigonometricTransition::s_Instance
  (
    GetExamples().AddEnum( ExampleMathTrigonometricTransition())
  );
} // namespace bcl
